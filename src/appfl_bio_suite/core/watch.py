"""Publish the federation as a map, and stream live runs onto it.

Built on `hivewatch <https://github.com/APPFL/hivewatch>`_, APPFL's monitoring toolkit.
hivewatch supplies the event schema, the map metadata artifact, the SSE server and the
Leaflet viewer; this module supplies the one thing hivewatch cannot know, which is who is
in *this* federation and what they are each doing.

TWO VIEWS, ONE MAP
------------------
**The network.** A federation exists before any run does. Partners are recruited, they
stand up endpoints, they unpack their data -- and until somebody launches something,
there is nothing to look at. So the network view is built from ``federation.yaml`` alone:
every site, where it is, which experiments it runs, what data it holds, how many samples
it declared, and whether its endpoint is answering right now. That is the page everyone
in the federation can open, and it is populated the day the config names them.

**The runs.** A run streams onto the same map through :class:`RunWatcher`: per-round,
per-site updates alongside the network view in the same runs directory, so one page shows
both the standing federation and what it is doing at this moment.

WHY THE GEOGRAPHY IS DECLARED, NOT DETECTED
-------------------------------------------
hivewatch's own APPFL example has each client call ``ipinfo.io`` and report the answer
upward. That is a reasonable default for a laptop demo and wrong for this suite in three
separate ways: a partner's compute node has no outbound web access, so the call fails; if
it succeeded it would return the institution's border router, not the cluster; and it
adds an outbound call from inside a partner's security boundary for the sole purpose of
drawing a dot -- which is a conversation with their security office that no dot is worth.

The coordinator already knows where their partners are. It is written down once, in
``federation.yaml``, next to everything else about that site. See
:class:`appfl_bio_suite.core.config.Location`.

WHAT IS DELIBERATELY NOT PUBLISHED
----------------------------------
Endpoint UUIDs are truncated to a fingerprint by default. They are not secrets -- a
partner's endpoint refuses tasks from anyone their identity mapping does not name, so
knowing one grants nothing -- but "grants nothing" is a claim about their configuration,
not ours, and the map is meant to be a page you can hand out. ``--include-endpoint-uuids``
opts in for an internal deployment.

Service accounts, ``data_dir`` and ``output_dir`` paths, and the coordinator's identity
never reach the map at all. Those describe the inside of a partner's cluster.
"""

from __future__ import annotations

import base64
import json
import logging
import queue
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from appfl_bio_suite.core.config import Federation, Site
from appfl_bio_suite.core.experiments import get_spec

__all__ = [
    "WatchError",
    "DEFAULT_RUNS_DIR",
    "DEFAULT_PORT",
    "NETWORK_RUN_ID",
    "Participation",
    "participations",
    "unplaced_sites",
    "unplaced_report",
    "build_network_events",
    "build_network_metadata",
    "metadata_from_events",
    "write_network_run",
    "export_site",
    "map_server",
    "write_sidecar",
    "watchable",
    "WATCHED_METRICS",
    "RunWatcher",
]

logger = logging.getLogger(__name__)

# Where run artifacts live. Under local/ because a run's metadata carries endpoint
# fingerprints and site names, and local/ is gitignored; `watch export` is the deliberate
# act that moves a copy somewhere publishable.
DEFAULT_RUNS_DIR = Path("local/watch/runs")

DEFAULT_PORT = 7070

# The network view is a hivewatch run like any other, which is what lets it sit in the
# same runs directory, appear in the same run list, and be served by the same server. It
# has a fixed id so that rebuilding it replaces it instead of accumulating copies.
NETWORK_RUN_ID = "network"

_PROTOCOL = "Globus Compute / APPFL"

_INSTALL_HINT = (
    "hivewatch is not installed.\n"
    "\n"
    "It is a coordinator-only extra -- the map is drawn on your machine from your\n"
    "federation config, and nothing about it reaches a partner:\n"
    "\n"
    "    pip install 'appfl-bio-suite[watch]' -c constraints.txt"
)


class WatchError(RuntimeError):
    """The map cannot be built: hivewatch is missing, or the federation has nothing to draw."""


def require_hivewatch():
    """Import hivewatch, or raise with the install line.

    Kept behind a function so that importing this module -- which ``launch`` does on every
    run -- never depends on an optional package.
    """
    try:
        import hivewatch
    except ImportError as exc:  # pragma: no cover - depends on install extras
        raise WatchError(_INSTALL_HINT) from exc
    return hivewatch


# ---------------------------------------------------------------------------
# what each site is doing
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Participation:
    """One site's involvement in one experiment, reduced to what the map shows."""

    experiment: str
    client_id: str
    endpoint_uuid: str
    data_kind: str
    samples: int | None = None

    @property
    def label(self) -> str:
        """The value of this experiment's chip on the site's card."""
        count = f"{self.samples:,} samples" if self.samples else "samples not declared"
        return f"{count} · {self.data_kind}" if self.data_kind else count


def participations(
    federation: Federation, experiment: str | None = None
) -> dict[str, list[Participation]]:
    """Site id -> what that site does, across every enabled experiment.

    Disabled experiments are skipped: ``enabled: false`` means the federation is not
    running it, and a map that shows it anyway is describing a federation that does not
    exist. Sites appear in declaration order, which is the order the config author chose.
    """
    names = [experiment] if experiment else list(federation.enabled_experiments())
    out: dict[str, list[Participation]] = {}

    for name in names:
        exp = federation.experiments.get(name)
        if exp is None or not exp.enabled:
            continue
        # An experiment declared in federation.yaml but unknown to this build has no
        # registry entry and therefore no data_kind. Draw it with what we do know rather
        # than refusing the whole map for one unrecognized name.
        try:
            data_kind = get_spec(name).data_kind
        except KeyError:
            data_kind = ""

        for entry in exp.sites:
            # FLamby names its sample count differently because it means a different
            # thing -- a training split of a public dataset, not a distributed bundle.
            samples = entry.expected_samples
            if samples is None:
                samples = entry.expected_train_samples
            out.setdefault(entry.site, []).append(
                Participation(
                    experiment=name,
                    client_id=entry.client_id,
                    endpoint_uuid=entry.endpoint_uuid,
                    data_kind=data_kind,
                    samples=samples,
                )
            )
    return out


def unplaced_sites(federation: Federation, experiment: str | None = None) -> list[Site]:
    """Participating sites with no coordinates. They cannot be drawn."""
    active = participations(federation, experiment)
    return [s for s in federation.sites if s.id in active and s.location is None]


def unplaced_report(federation: Federation, experiment: str | None = None) -> str:
    """A message naming the sites that cannot be drawn, and the YAML that fixes it."""
    missing = unplaced_sites(federation, experiment)
    if not missing:
        return ""
    lines = [
        f"{len(missing)} participating site(s) have no coordinates and are not on the map:",
        "",
    ]
    for site in missing:
        where = f" ({site.country})" if site.country else ""
        lines.append(f"  {site.id:16} {site.name}{where}")
    lines += [
        "",
        f"Add a `location` block to each under `sites:` in {federation._where()}:",
        "",
        "  - id: " + missing[0].id,
        "    name: " + missing[0].name,
        "    location:",
        "      lat: 41.8781",
        "      lng: -87.6298",
        "      city: Chicago",
        "",
        "Coordinates are declared rather than looked up -- see core/config.py::Location.",
        "Nothing else about a run depends on them; an unplaced site still trains.",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# the network view
# ---------------------------------------------------------------------------


def _fingerprint(uuid: str) -> str:
    """Enough of a UUID to match it against a config you already hold, and no more."""
    return f"{uuid[:8]}…{uuid[-4:]}"


def _site_client(
    site: Site,
    parts: list[Participation],
    *,
    status: str = "active",
    include_endpoint_uuids: bool = False,
) -> dict[str, Any]:
    """One site, as a hivewatch client update.

    ONE MARKER PER INSTITUTION, not one per experiment. A site hosting both GWAS and
    fine-mapping is one building with two endpoints in it, and hivewatch keys markers by
    ``client_id`` -- so two entries at identical coordinates would draw two dots on top of
    each other and report the federation as twice its real size. The experiments become
    fields on the one marker instead, which is also the more useful reading: the question
    the map answers is "who is in this federation and what do they hold".

    Every key that is not in hivewatch's own vocabulary is rendered as a labelled chip on
    the site's card, in insertion order. Keys prefixed with ``_`` are stored in the
    artifact but hidden from the viewer.
    """
    client: dict[str, Any] = {
        "client_id": site.id,
        "num_samples": sum(p.samples or 0 for p in parts),
        "status": status,
        "lat": site.location.lat if site.location else None,
        "lng": site.location.lng if site.location else None,
        "city": site.location.city if site.location else None,
        "country": site.country,
        # Chips, in the order they should read.
        "institution": site.name,
        "experiments": ", ".join(p.experiment for p in parts),
        "endpoints": len(parts),
        "scheduler": site.scheduler,
    }

    # One chip per experiment: "18,032 samples · Genotypes and phenotypes, PLINK bfile".
    # Keyed by the experiment name so the chip's own label says which experiment it is.
    for part in parts:
        client[part.experiment] = part.label

    for part in parts:
        if include_endpoint_uuids:
            # Opted in: a visible chip, for an internal deployment where the map is the
            # convenient place to look one up.
            client[f"endpoint · {part.experiment}"] = part.endpoint_uuid
        else:
            # Default: a fingerprint, stored under a `_`-prefixed key so hivewatch's
            # viewer keeps it out of the card. Enough to confirm which endpoint a site is
            # using against a config you already hold; not enough to be an identifier.
            client[f"_endpoint_{part.experiment}"] = _fingerprint(part.endpoint_uuid)

    return client


def _server_block(federation: Federation) -> dict[str, Any]:
    """The coordinator, as hivewatch's server metadata. Drawn as the hub."""
    coordinator = federation.coordinator
    server: dict[str, Any] = {
        "protocol": _PROTOCOL,
        "org": coordinator.organization or "Coordinator",
        "country": None,
        "city": None,
        "lat": None,
        "lng": None,
    }
    if coordinator.location is not None:
        server["lat"] = coordinator.location.lat
        server["lng"] = coordinator.location.lng
        server["city"] = coordinator.location.city
    # The coordinator has no `country` field of its own -- it is not a `Site`. If they
    # also train, their endpoint is declared and the country is unknown either way, so
    # this is left unset rather than guessed.
    return server


def build_network_events(
    federation: Federation,
    *,
    experiment: str | None = None,
    statuses: dict[str, str] | None = None,
    include_endpoint_uuids: bool = False,
    generated_at: datetime | None = None,
) -> list[dict[str, Any]]:
    """The network view as a hivewatch event stream.

    Events rather than a metadata blob because events are hivewatch's actual storage
    contract: the ``.jsonl`` is what its server lists runs from and what its viewer falls
    back to, and the ``.map.json`` is derived from it by hivewatch's own assembler. Going
    through events means this module never reimplements that schema and cannot drift from
    it.

    ``statuses`` maps site id to a hivewatch client status (``active``, ``idle``,
    ``dropped``, ``failed``). Omitted, every site is ``active`` -- the map then states
    membership, not liveness, which is the honest reading when nothing has been probed.
    """
    active = participations(federation, experiment)
    if not active:
        scope = f" for experiment '{experiment}'" if experiment else ""
        raise WatchError(
            f"no enabled experiment has any participating sites{scope} in "
            f"{federation._where()}. There is nothing to draw.\n"
            "Add sites under `experiments.<name>.sites`, or check `enabled:`."
        )

    now = (generated_at or datetime.now(UTC)).isoformat()
    statuses = statuses or {}

    clients = [
        _site_client(
            site,
            active[site.id],
            status=statuses.get(site.id, "active"),
            include_endpoint_uuids=include_endpoint_uuids,
        )
        for site in federation.sites
        if site.id in active
    ]

    scope = experiment or ", ".join(sorted({p.experiment for ps in active.values() for p in ps}))
    config = {
        "view": "federation network",
        "experiments": scope,
        "sites": len(clients),
        "generated_at": now,
        "source": str(federation.source_path) if federation.source_path else "federation config",
    }

    return [
        {
            "event_type": "init",
            "run_id": NETWORK_RUN_ID,
            "algorithm": "federation network",
            "config": config,
            "started_at": now,
        },
        {
            "event_type": "server_metadata",
            "run_id": NETWORK_RUN_ID,
            "timestamp": now,
            "server": _server_block(federation),
        },
        # A single round holding every site. The viewer's round slider is a time axis for
        # a training run; the network has no time axis, so it gets one frame.
        {
            "event_type": "round_end",
            "run_id": NETWORK_RUN_ID,
            "round": 1,
            "timestamp": now,
            "round_metrics": {
                "global_accuracy": None,
                "global_loss": None,
                "num_selected": len(clients),
                "num_completed": sum(1 for c in clients if c["status"] == "active"),
                "num_stragglers": 0,
                "round_duration_sec": None,
                "gradient_divergence": None,
            },
            "clients": clients,
        },
        {"event_type": "finished", "run_id": NETWORK_RUN_ID, "timestamp": now},
    ]


def build_network_metadata(federation: Federation, **kwargs) -> dict[str, Any]:
    """The network view as hivewatch map metadata, assembled by hivewatch itself."""
    return metadata_from_events(build_network_events(federation, **kwargs))


def metadata_from_events(events: list[dict[str, Any]]) -> dict[str, Any]:
    """Assemble map metadata from already-built events, using hivewatch's own assembler.

    Deliberately not reimplemented here. The ``.map.json`` schema is hivewatch's, and a
    second implementation of it in this repository would be one that drifts.
    """
    require_hivewatch()
    from hivewatch.map import build_map_metadata_from_events

    return build_map_metadata_from_events(events)


def write_network_run(
    federation: Federation,
    runs_dir: Path | str = DEFAULT_RUNS_DIR,
    *,
    catalog_path: Path | str | None = None,
    **kwargs,
) -> tuple[Path, Path]:
    """Write the network view into a runs directory. Returns (jsonl, map.json).

    Both artifacts, because hivewatch's server lists runs by globbing ``*.jsonl`` and
    serves detail from ``*.map.json``. A network view with only the metadata file would
    load if you asked for it by name and would not appear in the run list.
    """
    events, metadata = _network_artifacts(federation, catalog_path, **kwargs)

    runs_dir = Path(runs_dir)
    runs_dir.mkdir(parents=True, exist_ok=True)
    jsonl = runs_dir / f"{NETWORK_RUN_ID}.jsonl"
    mapjson = runs_dir / f"{NETWORK_RUN_ID}.map.json"

    # Truncating rather than appending: this file is a snapshot of the config as it is
    # now, not a history. Appending would leave the previous federation's sites in the
    # run list, which is precisely the stale-map failure the fixed run id avoids.
    jsonl.write_text("".join(json.dumps(e) + "\n" for e in events), encoding="utf-8")
    mapjson.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return jsonl, mapjson


def _network_artifacts(federation: Federation, catalog_path=None, **kwargs) -> tuple[list, dict]:
    from appfl_bio_suite.core.watch_catalog import add_partners, load_catalog

    try:
        catalog = load_catalog(catalog_path)
    except ValueError as exc:
        raise WatchError(str(exc)) from exc
    events = build_network_events(federation, **kwargs)
    add_partners(events, catalog, kwargs.get("experiment"))
    metadata = metadata_from_events(events)
    results = [
        group
        for group in catalog["results"]
        if not kwargs.get("experiment") or group["experiment"] == kwargs["experiment"]
    ]
    if results:
        metadata["results"] = results
    return events, metadata


# ---------------------------------------------------------------------------
# the publishable copy
# ---------------------------------------------------------------------------

# The extension is embedded into a fresh COPY of the pinned upstream viewer for both
# live serving and export. Upstream owns events, playback and Leaflet; our extension owns
# the federation UI and initialization. Keeping assets inline preserves the three-file
# static export, including the globe's geography and rendering dependencies.
_VIEWER_ASSETS = Path(__file__).parent / "watch_assets"
_UPSTREAM_STARTUP = "initDraggableLogPanel();\nif (METADATA_URL || EVENTS_URL)"

_INDEX = """<!doctype html>
<meta charset="utf-8">
<title>{title}</title>
<meta http-equiv="refresh" content="0; url=map.html?metadata_url=network.map.json">
<style>
  body {{ font-family: system-ui, sans-serif; margin: 4rem auto; max-width: 34rem;
         padding: 0 1.5rem; line-height: 1.6; }}
</style>
<h1>{title}</h1>
<p>Opening the map&hellip;
   <a href="map.html?metadata_url=network.map.json">continue</a> if it does not.</p>
"""


def _viewer_html() -> Path:
    """The viewer shipped inside the installed hivewatch package."""
    hivewatch = require_hivewatch()
    path = Path(hivewatch.__file__).resolve().parent / "map" / "hivewatch_map.html"
    if not path.is_file():  # pragma: no cover - would mean a broken hivewatch install
        raise WatchError(
            f"hivewatch is installed but its viewer is missing ({path}). Reinstall it:\n"
            "    pip install --force-reinstall hivewatch"
        )
    return path


def _patched_viewer() -> str:
    """Build our viewer without modifying the installed hivewatch package.

    Suppress upstream startup before installing the extension: it must own initial
    loading so network snapshots render immediately and live connection events cannot
    clear the selected view. Refuse an unexpected upstream layout rather than publishing
    a page whose event handlers or geography overrides only partly apply.
    """
    html = _viewer_html().read_text(encoding="utf-8")
    startup = html.find(_UPSTREAM_STARTUP)
    script_end = html.find("</script>", startup)
    if (
        startup == -1
        or html.count(_UPSTREAM_STARTUP) != 1
        or script_end == -1
        or "</head>" not in html
        or "</body>" not in html
    ):
        raise WatchError(
            "The installed hivewatch viewer has an unsupported layout. "
            "This UI requires hivewatch==0.2.1; reinstall the pinned [watch] extra."
        )
    html = (
        html[:startup]
        + "// appfl-bio-suite: viewer.js initializes the federation view.\n"
        + html[script_end:]
    )
    css = (_VIEWER_ASSETS / "viewer.css").read_text(encoding="utf-8")
    html = html.replace("</head>", f'<style id="bio-viewer-style">\n{css}\n</style>\n</head>', 1)
    scripts = ["<!-- appfl-bio-suite: federation viewer and globe -->"]
    logo = base64.b64encode((_VIEWER_ASSETS / "suite-logo.png").read_bytes()).decode("ascii")
    scripts.append(
        '<script type="application/json" id="bio-branding-data">'
        + json.dumps({"suite_logo": "data:image/png;base64," + logo})
        + "</script>"
    )
    notices = (_VIEWER_ASSETS / "THIRD_PARTY.md").read_text(encoding="utf-8")
    scripts.append(f'<script type="text/plain" id="bio-asset-notices">\n{notices}\n</script>')
    for name in ("vendor.js", "land.js", "globe.js", "results.js", "viewer.js"):
        source = (_VIEWER_ASSETS / name).read_text(encoding="utf-8")
        scripts.append(f'<script data-bio-asset="{name}">\n{source}\n</script>')
    return html.replace("</body>", "\n".join(scripts) + "\n</body>", 1)


@contextmanager
def map_server(
    *, host: str = "0.0.0.0", port: int = DEFAULT_PORT, runs_dir: Path | str = DEFAULT_RUNS_DIR
) -> Iterator[Any]:
    """Keep the same viewer used by export available for the server's lifetime.

    The temporary file belongs to this process; concurrent servers and exports never
    overwrite each other's viewer or the installed package. The context always stops
    the server before removing its HTML, including interruption and startup failures.
    """
    require_hivewatch()
    from hivewatch.map import MapServer

    class DirectoryMapServer(MapServer):
        def publish(self, payload: dict) -> None:
            # Hivewatch 0.2.1 sets _live_run_id for every published event, which makes
            # /runs return only that run. This viewer browses a directory: broadcasting
            # new activity must leave the network and all recorded runs discoverable.
            # Preserve upstream's bounded-queue fanout without selecting a server run.
            message = f"data: {json.dumps(payload)}\n\n"
            with self._lock:
                dead = []
                for subscriber in self._subscribers:
                    try:
                        subscriber.put_nowait(message)
                    except queue.Full:
                        dead.append(subscriber)
                for subscriber in dead:
                    self._subscribers.remove(subscriber)

    page = _patched_viewer()
    with TemporaryDirectory(prefix="appfl-bio-watch-") as directory:
        path = Path(directory) / "map.html"
        path.write_text(page, encoding="utf-8")
        server = DirectoryMapServer(
            host=host, port=port, runs_dir=str(runs_dir), map_path=str(path), watch=True
        )
        try:
            yield server
        finally:
            server.stop()


def export_site(
    federation: Federation,
    out_dir: Path | str,
    *,
    title: str = "APPFL federation network",
    catalog_path: Path | str | None = None,
    **kwargs,
) -> Path:
    """Write a self-contained static site for the network view. Returns the directory.

    THIS is the thing everyone in the federation can actually open. ``watch serve`` binds
    a port on the coordinator's login node, and a partner three time zones away cannot
    reach it -- HPC login nodes do not accept inbound connections from the internet, and
    should not. A directory of three static files can be dropped on GitHub Pages, an
    institutional web host, or any object store, and needs no server of ours at all.

    The map itself still calls out for its base tiles and for Leaflet, so a viewer needs
    ordinary web access. Nothing federation-specific leaves the page: the site data is in
    a JSON file served from alongside it.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    _, metadata = _network_artifacts(federation, catalog_path, **kwargs)

    (out_dir / "network.map.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    (out_dir / "map.html").write_text(_patched_viewer(), encoding="utf-8")
    (out_dir / "index.html").write_text(_INDEX.format(title=title), encoding="utf-8")
    return out_dir


# ---------------------------------------------------------------------------
# live runs
# ---------------------------------------------------------------------------


def write_sidecar(
    federation: Federation,
    experiment: str,
    path: Path | str,
    *,
    runs_dir: Path | str = DEFAULT_RUNS_DIR,
    port: int = DEFAULT_PORT,
    serve: bool = False,
) -> Path:
    """Write everything a driver needs to put its run on the map.

    The drivers take an APPFL server config and an APPFL client config and nothing else.
    That is deliberate -- they are handed a resolved run, not a federation -- and it is
    worth keeping: it is what makes a generated config the complete description of what
    will happen. So the map's extra knowledge travels in a sidecar written beside the
    generated configs, rather than by teaching the drivers to read federation.yaml.

    Sites are indexed by BOTH endpoint UUID and client id, because the two drivers key
    results differently: the Globus Compute driver by endpoint, the serial driver by
    client id. One lookup table serves both.
    """
    active = participations(federation, experiment)
    by_site = {s.id: s for s in federation.sites}

    clients: dict[str, dict[str, Any]] = {}
    for site_id, parts in active.items():
        site = by_site[site_id]
        for part in parts:
            entry = {
                "client_id": part.client_id,
                "institution": site.name,
                "experiment": part.experiment,
                "data": part.data_kind,
                "lat": site.location.lat if site.location else None,
                "lng": site.location.lng if site.location else None,
                "city": site.location.city if site.location else None,
                "country": site.country,
                "num_samples": part.samples,
            }
            clients[part.endpoint_uuid] = entry
            clients[part.client_id] = entry

    payload = {
        "experiment": experiment,
        "runs_dir": str(runs_dir),
        "port": port,
        "serve": serve,
        "server": _server_block(federation),
        "clients": clients,
    }
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path


# Metrics a site returns that mean something on a map. Anything else a trainer sends back
# -- arrays, per-variant tables, file paths -- is a result, not telemetry, and forwarding
# it into a page that may get published would be a leak dressed up as a feature. The
# allowlist is the enforcement; the drivers do not get to decide.
WATCHED_METRICS = (
    "local_accuracy",
    "local_loss",
    "num_samples",
    "gradient_norm",
    "train_time_sec",
    "bytes_sent",
)


def watchable(metadata: Any) -> dict[str, Any]:
    """The subset of a client's returned metadata that belongs on the map.

    Numeric scalars only. ``bool`` is excluded despite being an ``int`` -- a boolean
    rendered as a chip reading "1" is worse than no chip.
    """
    try:
        items = list(dict(metadata).items())
    except (TypeError, ValueError):
        return {}
    return {
        key: value
        for key, value in items
        if key in WATCHED_METRICS and isinstance(value, int | float) and not isinstance(value, bool)
    }


class RunWatcher:
    """Streams one run onto the map. Every method is safe to call unconditionally.

    NOTHING HERE MAY BREAK A RUN. A federated run costs a scheduler queue wait at three
    institutions; losing one because a monitoring library raised would be an absurd trade.
    So :meth:`load` returns an inert watcher when watching is off or hivewatch is absent,
    and every call is wrapped: a failure disables the watcher, logs once, and the run
    carries on without it.
    """

    def __init__(self, sidecar: dict[str, Any] | None = None):
        self._sidecar = sidecar or {}
        self._clients: dict[str, dict[str, Any]] = self._sidecar.get("clients", {})
        self._run = None
        self._live = False
        self._warned = False

    # -- construction -------------------------------------------------------

    @classmethod
    def load(cls, path: str | Path | None) -> RunWatcher:
        """Build a watcher from a sidecar path. ``None``, or an unreadable file, is inert."""
        if not path:
            return cls()
        try:
            return cls(json.loads(Path(path).read_text(encoding="utf-8")))
        except (OSError, ValueError) as exc:
            logger.warning("[watch] ignoring sidecar %s: %s", path, exc)
            return cls()

    @property
    def enabled(self) -> bool:
        return self._live

    # -- lifecycle ----------------------------------------------------------

    def start(self, algorithm: str = "APPFL", config: dict[str, Any] | None = None) -> None:
        if not self._sidecar:
            return
        try:
            import hivewatch
            from hivewatch.emitters import SSEEmitter
        except ImportError:
            logger.warning(
                "[watch] --watch was requested but hivewatch is not installed; "
                "the run will proceed without the map. %s",
                _INSTALL_HINT.splitlines()[-1].strip(),
            )
            return

        try:
            runs_dir = Path(self._sidecar.get("runs_dir", DEFAULT_RUNS_DIR))
            runs_dir.mkdir(parents=True, exist_ok=True)
            self._run = hivewatch.init(
                algorithm=algorithm,
                config=config or {},
                emitters=[
                    SSEEmitter(
                        runs_dir=str(runs_dir),
                        port=int(self._sidecar.get("port", DEFAULT_PORT)),
                        serve_map=bool(self._sidecar.get("serve", False)),
                    )
                ],
                verbose=False,
            )
            self._live = True
            self._run.set_server_metadata(**self._sidecar.get("server", {}))
        except Exception as exc:  # noqa: BLE001 - monitoring must never stop a run
            self._fail(exc)

    def round_start(self, round_no: int) -> None:
        self._guard(lambda: self._run.round_start(round_no))

    def client_update(self, key: str, round_no: int, **metrics: Any) -> None:
        """Log one site's result. ``key`` is an endpoint UUID or a client id.

        The site's identity and geography come from the sidecar, not from the worker.
        A partner returns statistics; where they are is something the coordinator already
        knew, and merging the two here is what keeps the worker payload free of anything
        that is not a result.
        """

        def _log() -> None:
            known = dict(self._clients.get(key, {}))
            client_id = known.pop("client_id", key)
            # A sample count the site actually reported beats the one the config
            # predicted; that difference is worth seeing on the map.
            known.update({k: v for k, v in metrics.items() if v is not None})
            self._run.log_client_update(client_id=client_id, round=round_no, **known)

        self._guard(_log)

    def round_end(self, round_no: int, **summary: Any) -> None:
        self._guard(lambda: self._run.log_round(round_no, **summary))

    def finish(self) -> None:
        self._guard(lambda: self._run.finish())
        self._live = False

    # -- internals ----------------------------------------------------------

    def _guard(self, call) -> None:
        if not self._live or self._run is None:
            return
        try:
            call()
        except Exception as exc:  # noqa: BLE001 - monitoring must never stop a run
            self._fail(exc)

    def _fail(self, exc: BaseException) -> None:
        self._live = False
        if not self._warned:
            self._warned = True
            logger.warning(
                "[watch] map disabled for this run: %s: %s. The run is unaffected.",
                type(exc).__name__,
                exc,
            )
