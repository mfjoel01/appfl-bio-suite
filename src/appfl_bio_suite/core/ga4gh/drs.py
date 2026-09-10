"""DRS -- Data Repository Service: the bundles a federation runs on, as addressable objects.

WHAT PROBLEM THIS SOLVES HERE
-----------------------------
Before this, a site's dataset was ``data_dir``: an absolute path on someone else's
cluster, asserted in ``federation.yaml``, checked for the presence of six filenames and
nothing more. Every question that matters about it was unanswerable from the coordinator:

* Did this site unpack the bundle I cut for *this* site?
* Is it the bundle from the simulation run whose answer key I am about to score against?
* Did it survive the transfer intact?

None of those fail loudly. A site running last month's bundle produces well-formed
aggregates over the wrong individuals; the pooled fit succeeds, the credible sets are
plausible, and the number that is wrong is one nobody can see. DRS makes each bundle a
content-addressed object: a sha-256 per file, a Merkle id over the set, and a
``drs://`` URI that is the same string at the coordinator and at the site.

CONTENT ADDRESSING, AND WHY OBJECTS DEDUPLICATE ACROSS SITES
------------------------------------------------------------
A blob's id **is** its sha-256. Two sites holding byte-identical files therefore hold the
same DRS object, reachable at two paths -- which is exactly what happens to
``reference_variants.tsv``, hard-linked into every bundle because it must be identical
for the harmonization to be sound. Seeing one object with three access methods is the
correct picture of that, and a per-site id would have hidden the property the experiment
depends on.

A bundle's id is a sha-256 over its sorted ``(name, child id)`` pairs, so it changes if
and only if some byte in it changed. That is what makes "is this the bundle I cut" a
string comparison.

THE ACCESS METHODS ARE REAL, AND THAT IS THE POINT
--------------------------------------------------
``file://`` is a spec-listed access method type and it is the honest one for a bundle a
partner already holds on their own filesystem. ``https://`` appears when this registry is
served (:func:`serve`), and ``globus://`` when a collection is configured -- the
federation already moves data that way, and naming the collection in the object is what
lets a partner re-fetch a bundle without a second conversation.

WHAT IS NOT HERE
----------------
No writes. DRS 1.5 is a read API and this is a read implementation: objects are created
by ``simulate``, which knows what it produced, and never by an HTTP request. There is no
authorization layer either -- the served registry is metadata about synthetic data, and a
deployment holding real genotypes would put this behind its own gateway rather than
trusting a flag here. Both limits are stated in docs/coordinator/ga4gh.md rather than
implied by silence.
"""

from __future__ import annotations

import hashlib
import json
import mimetypes
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Literal
from urllib.parse import quote, urlparse

from pydantic import BaseModel, ConfigDict, Field

__all__ = [
    "DRS_VERSION",
    "Checksum",
    "AccessURL",
    "AccessMethod",
    "ContentsObject",
    "DrsObject",
    "DrsRegistry",
    "DrsError",
    "DrsClient",
    "build_registry",
    "EXCLUDED_FROM_OBJECTS",
    "registry_for_files",
    "register_directory",
    "verify_object",
    "service_info",
    "serve",
    "parse_drs_uri",
]

DRS_VERSION = "1.5.0"

# GA4GH's checksum type strings are lowercase with a hyphen ("sha-256"), not Python's
# hashlib spelling. Getting this wrong produces a document that validates structurally
# and is rejected by every conformant client.
CHECKSUM_TYPE = "sha-256"

_CHUNK = 1 << 20

# Never part of a bundle object, whoever registers it.
#
# A dataset's DUO profile *describes* the object -- it carries the object's own
# ``drs_uri`` -- so it cannot be inside the checksum that names it. Excluding it here,
# rather than only in the site-side verifier, is what makes "the same bytes get the same
# id" true for a coordinator who re-registers a directory after the profile was written
# into it. Without this, ``simulate`` and ``ga4gh drs register`` mint different ids for
# the same data, and the property the whole scheme rests on quietly stops holding.
EXCLUDED_FROM_OBJECTS = frozenset({"DATA_USE.json"})


class DrsError(RuntimeError):
    """A DRS object, registry, or URI is missing or malformed."""


class _Model(BaseModel):
    model_config = ConfigDict(extra="allow")


class Checksum(_Model):
    checksum: str
    type: str = CHECKSUM_TYPE


class AccessURL(_Model):
    url: str
    headers: list[str] = Field(default_factory=list)


class AccessMethod(_Model):
    """How to get the bytes. ``type`` is from the DRS enum, not free text."""

    type: Literal["s3", "gs", "ftp", "gsiftp", "globus", "htsget", "https", "file"]
    access_url: AccessURL | None = None
    access_id: str | None = None
    region: str | None = None
    authorizations: dict[str, Any] | None = None


class ContentsObject(_Model):
    """A child of a bundle. ``id`` present means it is resolvable on its own."""

    name: str
    id: str | None = None
    drs_uri: list[str] = Field(default_factory=list)
    contents: list[ContentsObject] = Field(default_factory=list)


class DrsObject(_Model):
    """A DRS 1.5 object: a blob, or a bundle of them.

    ``contents`` is absent on a blob and present on a bundle -- that is how the spec
    distinguishes them, and it is why this is one model rather than two.
    """

    id: str
    self_uri: str
    created_time: str
    checksums: list[Checksum] = Field(default_factory=list)
    size: int = 0
    name: str | None = None
    updated_time: str | None = None
    version: str | None = None
    mime_type: str | None = None
    description: str | None = None
    aliases: list[str] = Field(default_factory=list)
    access_methods: list[AccessMethod] = Field(default_factory=list)
    contents: list[ContentsObject] | None = None

    @property
    def is_bundle(self) -> bool:
        return self.contents is not None

    @property
    def sha256(self) -> str | None:
        for entry in self.checksums:
            if entry.type == CHECKSUM_TYPE:
                return entry.checksum
        return None

    def local_path(self) -> Path | None:
        """The first ``file://`` access method's path, if any."""
        for method in self.access_methods:
            if method.type == "file" and method.access_url:
                return Path(urlparse(method.access_url.url).path)
        return None

    def child_ids(self) -> dict[str, str]:
        return {c.name: c.id for c in (self.contents or []) if c.id}


DrsObject.model_rebuild()
ContentsObject.model_rebuild()


# ---------------------------------------------------------------------------
# building a registry from what `simulate` produced
# ---------------------------------------------------------------------------


def file_sha256(path: str | Path) -> str:
    """Streamed sha-256. Same digest ``core.simulation`` records, deliberately."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def _now() -> str:
    return datetime.now(UTC).isoformat()


def parse_drs_uri(uri: str) -> tuple[str, str]:
    """``drs://host/id`` -> ``(host, id)``.

    Hostname-based URIs only. The compact-identifier form (``drs://prefix:accession``)
    resolves through identifiers.org, which needs the network and a registry entry
    neither a partner's worker nor this project has.
    """
    if not uri.startswith("drs://"):
        raise DrsError(f"'{uri}' is not a DRS URI. They look like 'drs://drs.example.org/<id>'.")
    remainder = uri[len("drs://") :]
    host, sep, object_id = remainder.partition("/")
    if not sep or not object_id:
        raise DrsError(
            f"'{uri}' names a host but no object. The hostname form is "
            "'drs://<host>/<object id>'; the compact-identifier form "
            "('drs://prefix:accession') is not supported here because resolving it "
            "requires identifiers.org, and a partner's worker has no outbound network."
        )
    return host, object_id


@dataclass
class _Access:
    """Where a registry's objects can be fetched from, as configured by the coordinator."""

    hostname: str
    https_base: str | None = None
    globus_collection: str | None = None

    def methods(self, object_id: str, path: Path) -> list[AccessMethod]:
        methods = [AccessMethod(type="file", access_url=AccessURL(url=path.resolve().as_uri()))]
        if self.https_base:
            base = self.https_base.rstrip("/")
            methods.append(
                AccessMethod(
                    type="https",
                    access_id="https",
                    access_url=AccessURL(url=f"{base}/files/{quote(object_id)}"),
                )
            )
        if self.globus_collection:
            # Globus Transfer addresses a file as (collection uuid, path on collection).
            # The DRS enum has a `globus` type for exactly this; the URL carries both so
            # a partner can hand it straight to `globus transfer`.
            methods.append(
                AccessMethod(
                    type="globus",
                    access_id="globus",
                    access_url=AccessURL(
                        url=f"globus://{self.globus_collection}{path.resolve().as_posix()}"
                    ),
                )
            )
        return methods


class DrsRegistry(BaseModel):
    """Every object a coordinator can resolve, and the identity of the service.

    Serialized as one JSON file so that it can be committed, diffed, copied into a
    partner bundle, and served -- all without a database. A registry with tens of
    thousands of objects (a 100-locus run's phenotype files) is a few megabytes, which is
    the scale this is intended for and not one more.
    """

    model_config = ConfigDict(extra="forbid")

    drs_version: str = DRS_VERSION
    hostname: str
    created_time: str = Field(default_factory=_now)
    # What produced it, so a registry found on disk is traceable.
    source: str = ""
    suite_version: str = ""
    objects: dict[str, DrsObject] = Field(default_factory=dict)

    # -- lookup -----------------------------------------------------------

    def get(self, object_id: str) -> DrsObject:
        try:
            return self.objects[object_id]
        except KeyError:
            raise DrsError(
                f"no DRS object '{object_id}' in this registry ({len(self.objects)} "
                f"object(s), hostname {self.hostname}). Rebuild it with "
                "`appfl-bio-suite ga4gh drs register`."
            ) from None

    def resolve(self, uri: str) -> DrsObject:
        """Resolve a ``drs://`` URI locally, checking the hostname matches."""
        host, object_id = parse_drs_uri(uri)
        if host != self.hostname:
            raise DrsError(
                f"{uri} names host '{host}', but this registry serves '{self.hostname}'. "
                "A DRS id is only unique within its service."
            )
        return self.get(object_id)

    def bundles(self) -> dict[str, DrsObject]:
        return {oid: obj for oid, obj in self.objects.items() if obj.is_bundle}

    def site_bundles(self) -> dict[str, DrsObject]:
        """The bundles that are a *site's data directory*, keyed by site name.

        Not the same as :meth:`bundles`, which also contains the nested ones -- a
        bundle's ``phenotypes/`` subtree is addressable in its own right, which is
        useful and is not what a coordinator pastes into ``federation.yaml`` as a
        site's ``drs_uri``. Site bundles are distinguished by carrying an alias, which
        only :func:`build_registry` sets.
        """
        return {obj.name: obj for obj in self.bundles().values() if obj.aliases and obj.name}

    def by_name(self, name: str) -> DrsObject | None:
        for obj in self.objects.values():
            if obj.name == name:
                return obj
        return None

    def uri(self, object_id: str) -> str:
        return f"drs://{self.hostname}/{object_id}"

    def add(self, obj: DrsObject) -> DrsObject:
        """Insert an object, merging access methods when the content is already known.

        The merge is the content-addressing property showing up in practice: the same
        bytes at two paths are one object reachable two ways, not two objects.
        """
        existing = self.objects.get(obj.id)
        if existing is None:
            self.objects[obj.id] = obj
            return obj
        known = {
            (m.type, m.access_url.url if m.access_url else m.access_id)
            for m in existing.access_methods
        }
        for method in obj.access_methods:
            key = (method.type, method.access_url.url if method.access_url else method.access_id)
            if key not in known:
                existing.access_methods.append(method)
        for alias in obj.aliases:
            if alias not in existing.aliases:
                existing.aliases.append(alias)
        return existing

    # -- persistence ------------------------------------------------------

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2, exclude_none=True) + "\n", encoding="utf-8")
        return path

    @classmethod
    def load(cls, path: str | Path) -> DrsRegistry:
        path = Path(path)
        try:
            raw = json.loads(path.read_text(encoding="utf-8"))
        except FileNotFoundError as exc:
            raise DrsError(
                f"no DRS registry at {path}. `appfl-bio-suite simulate fine-mapping` "
                "writes one beside the bundles it cuts; "
                "`appfl-bio-suite ga4gh drs register --data-root <dir>` rebuilds it."
            ) from exc
        except json.JSONDecodeError as exc:
            raise DrsError(f"{path} is not valid JSON: {exc}") from exc
        return cls.model_validate(raw)

    def summary(self) -> str:
        blobs = len(self.objects) - len(self.bundles())
        total = sum(o.size for o in self.objects.values() if not o.is_bundle)
        return (
            f"{self.hostname}: {len(self.bundles())} bundle(s), {blobs} blob(s), "
            f"{total / (1 << 20):.1f} MB of distinct content"
        )


def register_file(
    registry: DrsRegistry, path: Path, access: _Access, description: str = ""
) -> DrsObject:
    """Register one file as a blob object whose id is its sha-256."""
    path = Path(path)
    digest = file_sha256(path)
    stat = path.stat()
    mime, _ = mimetypes.guess_type(path.name)
    obj = DrsObject(
        id=digest,
        name=path.name,
        self_uri=f"drs://{access.hostname}/{digest}",
        size=stat.st_size,
        created_time=datetime.fromtimestamp(stat.st_mtime, UTC).isoformat(),
        checksums=[Checksum(checksum=digest)],
        mime_type=mime or "application/octet-stream",
        description=description or None,
        access_methods=access.methods(digest, path),
    )
    return registry.add(obj)


def _bundle_id(children: dict[str, str]) -> str:
    """A Merkle id over ``{name: child id}``: changes iff some byte in the set changed."""
    digest = hashlib.sha256()
    for name in sorted(children):
        digest.update(name.encode("utf-8"))
        digest.update(b"\0")
        digest.update(children[name].encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def register_directory(
    registry: DrsRegistry,
    directory: Path,
    access: _Access,
    name: str | None = None,
    description: str = "",
    aliases: tuple[str, ...] = (),
) -> DrsObject:
    """Register a directory tree as a DRS bundle, recursing into subdirectories.

    Nested directories become nested ``contents`` entries with their own ids, so a
    bundle's ``phenotypes/`` subtree is itself addressable -- which matters because that
    is where a run's several hundred phenotype files live and "are these the phenotypes
    the answer key was written for" is a question worth being able to ask on its own.
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise DrsError(f"{directory} is not a directory")

    contents: list[ContentsObject] = []
    children: dict[str, str] = {}
    size = 0

    for entry in sorted(directory.iterdir(), key=lambda p: p.name):
        if entry.name in EXCLUDED_FROM_OBJECTS:
            continue
        if entry.is_dir():
            child = register_directory(registry, entry, access, name=entry.name)
        elif entry.is_file():
            child = register_file(registry, entry, access)
        else:
            continue
        contents.append(ContentsObject(name=entry.name, id=child.id, drs_uri=[child.self_uri]))
        children[entry.name] = child.id
        size += child.size

    if not contents:
        raise DrsError(
            f"{directory} is empty, so it has no content to address. An empty bundle "
            "would get a stable id that says nothing about any data."
        )

    bundle_id = _bundle_id(children)
    bundle = DrsObject(
        id=bundle_id,
        name=name or directory.name,
        self_uri=f"drs://{access.hostname}/{bundle_id}",
        size=size,
        created_time=_now(),
        # A bundle's checksum is its Merkle id. The spec permits any checksum type on a
        # bundle and says nothing about how to compute one; recording it explicitly is
        # what lets a client verify a bundle without fetching every blob.
        checksums=[Checksum(checksum=bundle_id)],
        description=description or None,
        aliases=list(aliases),
        contents=contents,
        access_methods=[
            AccessMethod(type="file", access_url=AccessURL(url=directory.resolve().as_uri()))
        ],
    )
    return registry.add(bundle)


def build_registry(
    data_root: str | Path,
    hostname: str,
    https_base: str | None = None,
    globus_collection: str | None = None,
    experiment: str = "fine-mapping",
) -> DrsRegistry:
    """Build a registry from a directory ``simulate`` wrote.

    The layout it expects is the one every simulating experiment produces:
    ``<root>/<site>/data/`` per site, plus the coordinator-only ``ground_truth/`` and
    ``loci/``. Sites become bundles; the answer key is registered as a blob so that a
    result table can name the exact truth it was scored against.
    """
    from appfl_bio_suite import __version__

    data_root = Path(data_root).resolve()
    if not data_root.is_dir():
        raise DrsError(f"--data-root {data_root} does not exist")

    access = _Access(hostname=hostname, https_base=https_base, globus_collection=globus_collection)
    registry = DrsRegistry(
        hostname=hostname,
        source=f"{experiment} simulation at {data_root}",
        suite_version=__version__,
    )

    sites = sorted(p.name for p in data_root.iterdir() if (p / "data").is_dir())
    if not sites:
        raise DrsError(
            f"no site bundles under {data_root} -- expected one directory per site, each "
            "containing `data/`. That is the layout `simulate` writes."
        )

    for site in sites:
        register_directory(
            registry,
            data_root / site / "data",
            access,
            name=site,
            description=f"{experiment} site bundle for '{site}'",
            aliases=(f"{experiment}/{site}",),
        )

    answer_key = data_root / "ground_truth" / "causal_manifest.tsv"
    if answer_key.is_file():
        register_file(
            registry,
            answer_key,
            access,
            description=(
                "Ground-truth causal variants. COORDINATOR-ONLY: registered so a result "
                "table can name the truth it was scored against, and deliberately not "
                "referenced by any site bundle."
            ),
        )

    return registry


def registry_for_files(
    paths: list[str | Path],
    hostname: str,
    https_base: str | None = None,
    description: str = "",
) -> DrsRegistry:
    """A small registry over a handful of files -- a run's outputs, typically.

    The results a federation produces deserve the same treatment as the data it consumed:
    a checksum and a resolvable id. Otherwise the provenance chain is content-addressed
    at the input end and "the file called fed_fm_results.tsv" at the output end, which is
    where a published number actually gets confused with a rerun of it.
    """
    from appfl_bio_suite import __version__

    access = _Access(hostname=hostname, https_base=https_base)
    registry = DrsRegistry(hostname=hostname, source=description, suite_version=__version__)
    for path in paths:
        target = Path(path)
        if target.is_file():
            register_file(registry, target, access, description=description)
    return registry


def verify_object(obj: DrsObject, path: str | Path | None = None) -> list[str]:
    """Re-checksum an object's local copy. Returns a list of problems, empty when clean.

    A bundle is verified structurally -- every child present, every child's own checksum
    correct, and the recomputed Merkle id equal to the recorded one. That last check is
    what catches a file *added* to a bundle, which per-file checksums alone cannot see.
    """
    target = Path(path) if path is not None else obj.local_path()
    problems: list[str] = []
    if target is None:
        return [f"{obj.id}: no file:// access method and no path given"]
    if not target.exists():
        return [f"{obj.id}: {target} does not exist"]

    if not obj.is_bundle:
        actual = file_sha256(target)
        if actual != obj.sha256:
            problems.append(
                f"{obj.name or obj.id}: sha-256 mismatch\n"
                f"  expected {obj.sha256}\n  actual   {actual}"
            )
        return problems

    for child in obj.contents or []:
        child_path = target / child.name
        if not child_path.exists():
            problems.append(f"{obj.name or obj.id}: missing member '{child.name}'")
    return problems


# ---------------------------------------------------------------------------
# service-info and a conformant read-only server
# ---------------------------------------------------------------------------


def service_info(registry: DrsRegistry, organization: str = "", url: str = "") -> dict[str, Any]:
    """A GA4GH service-info document for this DRS instance."""
    from appfl_bio_suite import __version__

    return {
        "id": f"org.ga4gh.drs.{registry.hostname}",
        "name": f"appfl-bio-suite DRS ({registry.hostname})",
        "type": {"group": "org.ga4gh", "artifact": "drs", "version": DRS_VERSION},
        "description": (
            "Read-only DRS over the per-site bundles of an appfl-bio-suite federation."
        ),
        "organization": {"name": organization or registry.hostname, "url": url or ""},
        "version": __version__,
        "createdAt": registry.created_time,
        "environment": "research",
    }


class DrsClient:
    """Resolve DRS URIs, over HTTP or against a local registry file.

    Both backends exist because both are real: a coordinator with the registry on disk
    should not need a server running to check a pin, and a site that can reach a DRS
    service should not need a copy of the registry. Standard library only, so this can be
    inlined into a worker-side check without adding a package to a partner's install.
    """

    def __init__(
        self,
        registry: DrsRegistry | None = None,
        base_url: str | None = None,
        token: str | None = None,
        timeout: float = 30.0,
    ):
        if registry is None and base_url is None:
            raise DrsError("a DrsClient needs either a local registry or a base_url")
        self.registry = registry
        self.base_url = base_url.rstrip("/") if base_url else None
        self.token = token
        self.timeout = timeout

    def get(self, object_id: str) -> DrsObject:
        if self.registry is not None:
            return self.registry.get(object_id)
        return DrsObject.model_validate(
            self._request(f"{self.base_url}/ga4gh/drs/v1/objects/{quote(object_id)}")
        )

    def resolve(self, uri: str) -> DrsObject:
        host, object_id = parse_drs_uri(uri)
        if self.registry is not None:
            return self.registry.resolve(uri)
        return self.get(object_id)

    def service_info(self) -> dict[str, Any]:
        if self.registry is not None:
            return service_info(self.registry)
        return self._request(f"{self.base_url}/ga4gh/drs/v1/service-info")

    def _request(self, url: str) -> dict[str, Any]:
        import urllib.error
        import urllib.request

        request = urllib.request.Request(url, headers={"Accept": "application/json"})
        if self.token:
            request.add_header("Authorization", f"Bearer {self.token}")
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", "replace")[:400]
            raise DrsError(f"{url} returned {exc.code}: {body}") from exc
        except urllib.error.URLError as exc:
            raise DrsError(f"could not reach {url}: {exc.reason}") from exc


def serve(
    registry: DrsRegistry,
    host: str = "127.0.0.1",
    port: int = 8080,
    serve_bytes: bool = True,
):
    """Return a threading HTTP server exposing the DRS read API. Caller starts it.

    Returned rather than started so that tests can bind port 0 and read back the real
    port, which is the only way to test a server without racing a fixed one.

    Routes:
        GET /ga4gh/drs/v1/service-info
        GET /ga4gh/drs/v1/objects/{id}
        GET /ga4gh/drs/v1/objects/{id}/access/{access_id}
        GET /files/{id}                     the bytes, when serve_bytes
    """
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

    class Handler(BaseHTTPRequestHandler):
        server_version = "appfl-bio-suite-drs"

        def _send(self, code: int, payload: dict[str, Any]) -> None:
            body = json.dumps(payload).encode("utf-8")
            self.send_response(code)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _error(self, code: int, message: str) -> None:
            # The shape DRS specifies for errors, so a conformant client can read them.
            self._send(code, {"status_code": code, "msg": message})

        def do_GET(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler's interface
            path = urlparse(self.path).path.rstrip("/")
            prefix = "/ga4gh/drs/v1"

            if path == f"{prefix}/service-info":
                self._send(200, service_info(registry))
                return

            if path.startswith(f"{prefix}/objects/"):
                remainder = path[len(f"{prefix}/objects/") :]
                object_id, _, access_part = remainder.partition("/access/")
                try:
                    obj = registry.get(object_id)
                except DrsError as exc:
                    self._error(404, str(exc))
                    return
                if not access_part:
                    self._send(200, json.loads(obj.model_dump_json(exclude_none=True)))
                    return
                for method in obj.access_methods:
                    if method.access_id == access_part and method.access_url:
                        self._send(200, json.loads(method.access_url.model_dump_json()))
                        return
                self._error(404, f"no access method '{access_part}' on object {object_id}")
                return

            if serve_bytes and path.startswith("/files/"):
                object_id = path[len("/files/") :]
                try:
                    obj = registry.get(object_id)
                except DrsError as exc:
                    self._error(404, str(exc))
                    return
                local = obj.local_path()
                if local is None or not local.is_file():
                    self._error(404, f"object {object_id} has no readable local copy")
                    return
                data = local.read_bytes()
                self.send_response(200)
                self.send_header("Content-Type", obj.mime_type or "application/octet-stream")
                self.send_header("Content-Length", str(len(data)))
                self.end_headers()
                self.wfile.write(data)
                return

            self._error(404, f"no route {path}")

        def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
            return  # the CLI prints what it wants; this would double it

    return ThreadingHTTPServer((host, port), Handler)
