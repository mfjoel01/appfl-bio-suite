"""Keep a published map's endpoint liveness current, without republishing the map.

``watch export --probe`` asks Globus Compute once and bakes the answer into
``network.map.json``. That is a photograph: an hour later the page still says what was
true when it was taken. A hosted map that claims a site is answering has to keep asking.

So liveness is a second, small file beside the map -- ``endpoints.status.json`` -- that a
loop rewrites every minute, and that the viewer polls. The map itself (institutions,
coordinates, stages, embedded results: megabytes) is exported once and left alone.

WHY THE PROBE LIST IS ITS OWN FILE
----------------------------------
The published site is drawn from ``website/federation.yaml``, whose endpoint UUIDs are
placeholders on purpose; the real ones live in ``local/federation.yaml``, next to partner
contacts, and ``scripts/check_site_preview.py`` exists to keep that file away from the
site build. A live map needs the real UUIDs and nothing else from that file. The targets
file is exactly that reduction -- public site id to endpoint UUIDs -- so the host that
serves the page never holds a contact, a path or a service account, and the page it
serves never holds a UUID: the status file says how many of a site's endpoints answered,
not which.

WHAT THE STATE MEANS
--------------------
``get_endpoint_status`` reports that the endpoint's daemon is connected to the Globus
Compute service. It is not the round trip ``endpoint smoke`` makes -- see
:mod:`appfl_bio_suite.core.endpoint` for everything that can be wrong behind an
``online`` -- and the viewer labels it accordingly. It is the check that is cheap enough
to repeat every minute.

WHO MAY ASK
-----------
Not everyone. A partner's multi-user endpoint answered a status request from an identity
that does not own it; a single-user endpoint owned by a service account refused one with
403. So the prober's identity decides what it can see, and an endpoint it may not ask
about is ``unknown`` on the map -- which is the true statement -- until the prober is
given an identity that endpoint accepts.

``unknown`` is not ``offline``, for the same reason and one more: if the prober's own
credential stops working, every request fails, and reporting that as a federation-wide
outage would be a false statement about ten institutions caused by one expired secret.
"""

from __future__ import annotations

import json
import logging
import os
import re
import time
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import yaml

from appfl_bio_suite.core.config import Federation

__all__ = [
    "LiveStatusError",
    "STATUS_FILENAME",
    "DEFAULT_INTERVAL",
    "load_targets",
    "targets_from_federation",
    "probe_targets",
    "write_status",
    "run",
]

logger = logging.getLogger(__name__)

# The viewer looks for this name beside network.map.json when the map's metadata says a
# live status is being published. Fixed, so a deployment has one fewer thing to agree on.
STATUS_FILENAME = "endpoints.status.json"

DEFAULT_INTERVAL = 60.0

# One request per endpoint per cycle against a shared public API. Faster than this buys
# nothing a person can see and is not ours to spend.
MIN_INTERVAL = 15.0

_UUID = re.compile(r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$")

# federation.yaml marks an endpoint that does not exist yet with an obviously synthetic
# UUID (dddddddd-0000-0000-0000-000000000002). Probing one returns 404, which would paint
# a site that is still onboarding as down.
_PLACEHOLDER = re.compile(r"^[0-9a-f]{8}-0000-0000-0000-")


class LiveStatusError(ValueError):
    """The probe list is missing or malformed."""


def is_placeholder(uuid: str) -> bool:
    return _PLACEHOLDER.match(uuid.lower()) is not None


def load_targets(path: Path | str) -> dict[str, list[str]]:
    """Read a targets file: ``{"sites": {"<site id>": ["<endpoint uuid>", ...]}}``.

    YAML or JSON. Site ids are the ones the *published* map uses, which need not be the
    ones in ``local/federation.yaml``.
    """
    path = Path(path)
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise LiveStatusError(f"cannot read probe targets {path}: {exc}") from exc
    sites = raw.get("sites") if isinstance(raw, dict) else None
    if not isinstance(sites, dict) or not sites:
        raise LiveStatusError(f"{path}: expected a non-empty `sites:` mapping")

    targets: dict[str, list[str]] = {}
    for site, uuids in sites.items():
        if isinstance(uuids, str):
            uuids = [uuids]
        if not isinstance(uuids, list) or not uuids:
            raise LiveStatusError(f"{path}: site {site!r} needs a list of endpoint UUIDs")
        cleaned = []
        for uuid in uuids:
            if not isinstance(uuid, str) or not _UUID.match(uuid.lower()):
                raise LiveStatusError(f"{path}: site {site!r}: not an endpoint UUID: {uuid!r}")
            cleaned.append(uuid.lower())
        targets[str(site)] = list(dict.fromkeys(cleaned))
    return targets


def targets_from_federation(
    federation: Federation, rename: dict[str, str] | None = None
) -> tuple[dict[str, list[str]], list[str]]:
    """The probe list a federation implies, and notes on what was left out.

    One entry per site, holding each distinct endpoint that site runs an enabled
    experiment on. ``rename`` maps a site id in this federation to the id the published
    map uses for the same institution.
    """
    from appfl_bio_suite.core.watch import participations

    rename = rename or {}
    targets: dict[str, list[str]] = {}
    notes: list[str] = []
    for site_id, parts in participations(federation).items():
        public_id = rename.get(site_id, site_id)
        for part in parts:
            label = f"{site_id} / {part.experiment}"
            if part.endpoint_uuid is None:
                notes.append(f"{label}: no endpoint registered yet")
            elif is_placeholder(part.endpoint_uuid):
                notes.append(f"{label}: placeholder UUID {part.endpoint_uuid}, not probed")
            else:
                uuids = targets.setdefault(public_id, [])
                if part.endpoint_uuid.lower() not in uuids:
                    uuids.append(part.endpoint_uuid.lower())
    return targets, notes


def _site_state(online: int, errors: int, total: int) -> str:
    if errors == total:
        return "unknown"
    if online == total:
        return "online"
    # Some answering and some not is neither up nor down, and is exactly the state worth
    # seeing at a glance -- the same reasoning as `idle` in `watch build --probe`.
    return "degraded" if online else "offline"


def probe_targets(
    targets: dict[str, list[str]],
    *,
    client: Any = None,
    interval: float = DEFAULT_INTERVAL,
    now: datetime | None = None,
) -> dict[str, Any]:
    """Ask Globus Compute about every endpoint once. Never raises for a failed request."""
    from appfl_bio_suite.core.endpoint import get_status

    sites: dict[str, dict[str, Any]] = {}
    for site, uuids in targets.items():
        results = [get_status(uuid, label=site, client=client) for uuid in uuids]
        online = sum(r.online for r in results)
        errors = sum(r.error is not None for r in results)
        for result in results:
            if result.error:
                # To the log, never to the status file: the message quotes the UUID.
                logger.warning("%s: %s", site, result.error)
        sites[site] = {
            "state": _site_state(online, errors, len(results)),
            "endpoints": len(results),
            "online": online,
        }
    return {
        "schema": 1,
        "checked_at": (now or datetime.now(UTC)).isoformat(timespec="seconds"),
        "interval_sec": interval,
        "sites": sites,
    }


def _unknown(targets: dict[str, list[str]], interval: float) -> dict[str, Any]:
    return {
        "schema": 1,
        "checked_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "interval_sec": interval,
        "sites": {
            site: {"state": "unknown", "endpoints": len(uuids), "online": 0}
            for site, uuids in targets.items()
        },
    }


def write_status(status: dict[str, Any], path: Path | str) -> Path:
    """Replace the status file atomically.

    A web server is reading this file while it is rewritten. Writing in place would hand
    some request a truncated document every so often; a rename never does.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    scratch = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    scratch.write_text(json.dumps(status, indent=1) + "\n", encoding="utf-8")
    # Readable by the web server whatever this process's umask is.
    scratch.chmod(0o644)
    os.replace(scratch, path)
    return path


def run(
    targets: dict[str, list[str]],
    out: Path | str,
    *,
    interval: float = DEFAULT_INTERVAL,
    once: bool = False,
    client_factory: Callable[[], Any] | None = None,
    report: Callable[[dict[str, Any]], None] | None = None,
    sleep: Callable[[float], None] = time.sleep,
) -> None:
    """Probe and publish, forever (or once).

    Nothing that goes wrong in a cycle ends the loop. This runs unattended on a host
    nobody is logged in to, and the alternative to carrying on is a page that silently
    stops being true; a failed cycle publishes ``unknown`` and the next one tries again
    with a fresh client.
    """
    if not once and interval < MIN_INTERVAL:
        raise LiveStatusError(f"--interval must be at least {MIN_INTERVAL:g} seconds")
    if client_factory is None:
        from appfl_bio_suite.core.endpoint import _client as client_factory

    client = None
    while True:
        started = time.monotonic()
        try:
            if client is None:
                client = client_factory()
            status = probe_targets(targets, client=client, interval=interval)
            if all(site["state"] == "unknown" for site in status["sites"].values()):
                # Every request failed: more likely the client than the federation.
                client = None
        except Exception:  # noqa: BLE001 - see docstring; the loop must survive
            logger.exception("probe cycle failed")
            client = None
            status = _unknown(targets, interval)
        try:
            write_status(status, out)
        except OSError:
            logger.exception("could not write %s", out)
            if once:
                raise
        if report is not None:
            report(status)
        if once:
            return
        sleep(max(1.0, interval - (time.monotonic() - started)))
