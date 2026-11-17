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

That identity is a *route*: an endpoint in the probe list may name the service-account
client to ask as, and the environment variable holding its secret, in the same shape as
the coordinator's ``--site-credentials`` manifest. The secret itself is never in the
probe list. Note what a route hands the prober: a single-user endpoint's owning client can
submit tasks to it, not only read its status, so give the prober its own secret on that
client -- revocable without disturbing the coordinator's.

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
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import yaml

from appfl_bio_suite.core.config import Federation

__all__ = [
    "LiveStatusError",
    "Route",
    "load_probe_list",
    "missing_secrets",
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


_ENV_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


class LiveStatusError(ValueError):
    """The probe list is missing or malformed."""


@dataclass(frozen=True)
class Route:
    """Ask about one endpoint as this service-account client instead of the default."""

    client_id: str
    client_secret_env: str

    def client(self):
        from appfl_bio_suite.core.site_credentials import SiteCredential

        return SiteCredential(self.client_id, self.client_secret_env, "status").client()


def is_placeholder(uuid: str) -> bool:
    return _PLACEHOLDER.match(uuid.lower()) is not None


def load_probe_list(path: Path | str) -> tuple[dict[str, list[str]], dict[str, Route]]:
    """Read a probe list. Returns (site id -> endpoint UUIDs, endpoint UUID -> route).

    YAML or JSON::

        {"sites":  {"<site id>": ["<endpoint uuid>", ...]},
         "routes": {"<endpoint uuid>": {"client_id": "...", "client_secret_env": "NAME"}}}

    Site ids are the ones the *published* map uses, which need not be the ones in
    ``local/federation.yaml``. ``routes`` is optional; an endpoint without one is asked
    about as the prober's default identity.
    """
    path = Path(path)
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise LiveStatusError(f"cannot read probe targets {path}: {exc}") from exc
    if not isinstance(raw, dict) or set(raw) - {"sites", "routes"}:
        raise LiveStatusError(f"{path}: expected `sites:` and optionally `routes:`, nothing else")
    sites = raw.get("sites")
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

    probed = {uuid for uuids in targets.values() for uuid in uuids}
    routes: dict[str, Route] = {}
    raw_routes = raw.get("routes") or {}
    if not isinstance(raw_routes, dict):
        raise LiveStatusError(f"{path}: `routes:` must map endpoint UUIDs to routes")
    for uuid, entry in raw_routes.items():
        where = f"{path}: route for {uuid}"
        if str(uuid).lower() not in probed:
            raise LiveStatusError(f"{where}: that endpoint is not in `sites:`")
        if not isinstance(entry, dict) or set(entry) != {"client_id", "client_secret_env"}:
            # Refused rather than ignored: the likeliest extra field is the secret itself.
            raise LiveStatusError(
                f"{where}: needs exactly client_id and client_secret_env. The secret goes "
                "in that environment variable, never in this file."
            )
        client_id, env = str(entry["client_id"]).lower(), str(entry["client_secret_env"])
        if not _UUID.match(client_id) or not _ENV_NAME.match(env):
            raise LiveStatusError(f"{where}: client_id must be a UUID, client_secret_env a name")
        routes[str(uuid).lower()] = Route(client_id, env)
    return targets, routes


def load_targets(path: Path | str) -> dict[str, list[str]]:
    """The ``sites`` half of a probe list."""
    return load_probe_list(path)[0]


def missing_secrets(routes: Mapping[str, Route], environ: Mapping[str, str] | None = None):
    """Environment variables the routes name that are not set, sorted."""
    environ = os.environ if environ is None else environ
    return sorted(
        {r.client_secret_env for r in routes.values() if not environ.get(r.client_secret_env)}
    )


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


def _check(targets: dict[str, list[str]], client_for: Callable[[str], Any]) -> dict[str, Any]:
    """Endpoint UUID -> EndpointStatus, one request each. Never raises."""
    from appfl_bio_suite.core.endpoint import EndpointStatus, get_status

    results = {}
    for site, uuids in targets.items():
        for uuid in uuids:
            if uuid in results:
                continue
            try:
                client = client_for(uuid)
            except Exception as exc:  # noqa: BLE001 - one bad route must not end the cycle
                results[uuid] = EndpointStatus(
                    uuid=uuid, label=site, status="error", error=f"{type(exc).__name__}: {exc}"
                )
                continue
            results[uuid] = get_status(uuid, label=site, client=client)
            if results[uuid].error:
                # To the log, never to the status file: the message quotes the UUID.
                logger.warning("%s: %s", site, results[uuid].error)
    return results


def _summarize(
    targets: dict[str, list[str]], results: dict[str, Any], interval: float, now=None
) -> dict[str, Any]:
    sites: dict[str, dict[str, Any]] = {}
    for site, uuids in targets.items():
        mine = [results[uuid] for uuid in uuids]
        online = sum(r.online for r in mine)
        errors = sum(r.error is not None for r in mine)
        sites[site] = {
            "state": _site_state(online, errors, len(mine)),
            "endpoints": len(mine),
            "online": online,
        }
    return {
        "schema": 1,
        "checked_at": (now or datetime.now(UTC)).isoformat(timespec="seconds"),
        "interval_sec": interval,
        "sites": sites,
    }


def probe_targets(
    targets: dict[str, list[str]],
    *,
    client: Any = None,
    interval: float = DEFAULT_INTERVAL,
    now: datetime | None = None,
) -> dict[str, Any]:
    """Ask Globus Compute about every endpoint once, as one identity. Never raises."""
    return _summarize(targets, _check(targets, lambda _uuid: client), interval, now)


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
    routes: Mapping[str, Route] | None = None,
    interval: float = DEFAULT_INTERVAL,
    once: bool = False,
    client_factory: Callable[[], Any] | None = None,
    route_factory: Callable[[Route], Any] | None = None,
    report: Callable[[dict[str, Any]], None] | None = None,
    sleep: Callable[[float], None] = time.sleep,
) -> None:
    """Probe and publish, forever (or once).

    Nothing that goes wrong in a cycle ends the loop. This runs unattended on a host
    nobody is logged in to, and the alternative to carrying on is a page that silently
    stops being true; a failed cycle publishes ``unknown`` and the next one tries again.

    One client per identity -- the default, and each route's -- reused across cycles. An
    identity whose every request failed is rebuilt next cycle: that is more likely its
    credential than every endpoint it asks about.
    """
    if not once and interval < MIN_INTERVAL:
        raise LiveStatusError(f"--interval must be at least {MIN_INTERVAL:g} seconds")
    if client_factory is None:
        from appfl_bio_suite.core.endpoint import _client as client_factory
    if route_factory is None:
        route_factory = Route.client
    routes = dict(routes or {})

    def identity(uuid: str) -> Route | None:
        return routes.get(uuid)

    clients: dict[Route | None, Any] = {}

    def client_for(uuid: str):
        key = identity(uuid)
        if key not in clients:
            clients[key] = client_factory() if key is None else route_factory(key)
        return clients[key]

    while True:
        started = time.monotonic()
        try:
            results = _check(targets, client_for)
            status = _summarize(targets, results, interval)
            for key in {identity(uuid) for uuid in results}:
                asked = [r for uuid, r in results.items() if identity(uuid) == key]
                if all(r.error is not None for r in asked):
                    clients.pop(key, None)
        except Exception:  # noqa: BLE001 - see docstring; the loop must survive
            logger.exception("probe cycle failed")
            clients.clear()
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
