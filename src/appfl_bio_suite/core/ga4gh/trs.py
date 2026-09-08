"""TRS -- Tool Registry Service: which code, exactly, computed a site's aggregates.

THE GAP THIS FILLS
------------------
APPFL ships a trainer by reading ``trainer.py`` on the driver and sending its *source
text* to a worker. That is a good design -- it is why a partner installs two packages and
no suite -- but it means the code a site runs is whatever the coordinator's checkout
happened to contain at dispatch time. Nothing in the payload, the results table, or the
run manifest says which version that was. Two runs a month apart, one of them against an
edited trainer, are indistinguishable afterwards.

TRS closes that: the site stage is a registered tool with a version id, a descriptor, a
per-file checksum, and a container image. ``federation.yaml`` pins the tool
(:class:`ToolPin`), ``preflight`` verifies the pin against what is installed, and the
verified pin is written into every result table's provenance. "Which tool version
produced this number" becomes a lookup instead of an archaeology exercise.

WHY A CWL DESCRIPTOR EXISTS AT ALL
----------------------------------
A TRS tool that describes nothing runnable is a metadata exercise. The descriptor here
wraps a real command -- ``appfl-bio-suite site-stage`` -- which is the same site
computation the APPFL trainer performs, reachable from a shell. Writing the descriptor is
what forced that entry point to exist (``experiments/fine_mapping/site_stage.py``), and
it is what a TES task's executor invokes. So the three specifications are wired together
rather than merely co-present: TRS names the tool and pins its image, TES runs that
image, and the tool reads DRS-addressed inputs.

WHAT THIS IS NOT
----------------
Not a TRS *server* with a database. :func:`write_registry` emits a static tree of JSON
documents at the spec's paths, servable by any static file host -- which is the whole of
what a read-only TRS needs to be, and is what a coordinator can actually publish next to
their results. Dockstore registration is supported the way Dockstore expects it, through
a generated ``.dockstore.yml``; this module does not call Dockstore's API.

The container image is *named* here, not built here. A coordinator builds and pushes it
from the generated ``Containerfile``; until they do, the TES path has nothing to run and
says so. Pretending otherwise would be the one dishonest thing in this file.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any
from urllib.parse import quote

from pydantic import BaseModel, ConfigDict, Field

__all__ = [
    "TRS_VERSION",
    "Checksum",
    "ImageData",
    "FileWrapper",
    "ToolFile",
    "ToolVersion",
    "ToolClass",
    "Tool",
    "ToolPin",
    "TrsError",
    "TrsClient",
    "build_tool",
    "write_registry",
    "dockstore_yaml",
    "descriptor_checksum",
    "verify_pin",
    "service_info",
    "serve",
]

TRS_VERSION = "2.0.1"
CHECKSUM_TYPE = "sha-256"

# The TRS tool class for a runnable command-line tool. Dockstore uses these ids.
COMMAND_LINE_TOOL = {
    "id": "0",
    "name": "CommandLineTool",
    "description": "CWL described CommandLineTool",
}


class TrsError(RuntimeError):
    """A tool document, descriptor, or pin is missing or does not verify."""


class _Model(BaseModel):
    model_config = ConfigDict(extra="allow")


class Checksum(_Model):
    checksum: str
    type: str = CHECKSUM_TYPE


class ImageData(_Model):
    """A container image, ideally pinned by digest rather than by tag.

    A tag is a moving target: ``:0.1.0`` can be repushed, and a federation that agreed on
    a tag has agreed on nothing durable. ``checksum`` carries the digest when the
    coordinator has pushed the image and knows it.
    """

    image_name: str
    image_type: str = "Docker"
    registry_host: str | None = None
    size: int | None = None
    updated: str | None = None
    checksum: list[Checksum] = Field(default_factory=list)

    @property
    def digest(self) -> str | None:
        for entry in self.checksum:
            if entry.type in ("sha-256", "sha256"):
                return entry.checksum
        return None


class FileWrapper(_Model):
    """A descriptor's content, as TRS returns it from ``.../descriptor``."""

    content: str | None = None
    checksum: list[Checksum] = Field(default_factory=list)
    url: str | None = None
    image_type: str | None = None


class ToolFile(_Model):
    """One file belonging to a tool version, with the role it plays."""

    path: str
    # PRIMARY_DESCRIPTOR | SECONDARY_DESCRIPTOR | TEST_FILE | CONTAINERFILE | OTHER
    file_type: str
    checksum: list[Checksum] = Field(default_factory=list)


class ToolVersion(_Model):
    id: str
    name: str
    url: str = ""
    author: list[str] = Field(default_factory=list)
    descriptor_type: list[str] = Field(default_factory=lambda: ["CWL"])
    containerfile: bool = False
    images: list[ImageData] = Field(default_factory=list)
    included_apps: list[str] = Field(default_factory=list)
    is_production: bool = False
    meta_version: str | None = None
    signed: bool = False
    verified: bool = False
    verified_source: list[str] = Field(default_factory=list)


class ToolClass(_Model):
    id: str
    name: str
    description: str = ""


class Tool(_Model):
    id: str
    url: str = ""
    name: str = ""
    description: str = ""
    organization: str = ""
    aliases: list[str] = Field(default_factory=list)
    toolclass: ToolClass = Field(default_factory=lambda: ToolClass(**COMMAND_LINE_TOOL))
    has_checker: bool = False
    checker_url: str | None = None
    meta_version: str | None = None
    versions: list[ToolVersion] = Field(default_factory=list)

    def version(self, version_id: str) -> ToolVersion:
        for entry in self.versions:
            if entry.id == version_id:
                return entry
        known = ", ".join(v.id for v in self.versions) or "none"
        raise TrsError(f"tool {self.id} has no version '{version_id}'. Known: {known}.")


class ToolPin(BaseModel):
    """What ``federation.yaml`` records so a run can prove which tool it ran.

    Not a TRS type -- TRS has no notion of "the version I insist on". The three fields
    are the three things that can drift independently: the tool, its version, and the
    bytes of that version's files.
    """

    model_config = ConfigDict(extra="forbid")

    id: str
    version: str
    # Merkle checksum over the version's files. Optional so a coordinator can adopt the
    # pin in two steps -- name the tool first, freeze its content once they have run
    # `ga4gh trs publish` -- but preflight warns while it is empty, because a pin without
    # a checksum pins a label.
    descriptor_checksum: str | None = None
    image: str | None = None
    registry_url: str | None = None

    def render(self) -> str:
        out = f"{self.id}@{self.version}"
        if self.descriptor_checksum:
            out += f" ({CHECKSUM_TYPE}:{self.descriptor_checksum[:16]}...)"
        return out


# ---------------------------------------------------------------------------
# building the tool from what is actually installed
# ---------------------------------------------------------------------------


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def repository_slug() -> str:
    """This install's repository, from the one constant that already names it.

    Read from ``appfl_bio_suite.REPO_URL`` rather than written again here, because a TRS
    id names *whose* tool this is. A coordinator publishing from their own fork changes
    that constant once -- as they already must, for the install spec partners are given
    -- and their tool ids follow. A second literal here would silently attribute their
    runs to somebody else's repository, which is the one thing a tool registry exists to
    prevent.
    """
    from appfl_bio_suite import REPO_URL

    return REPO_URL.removeprefix("https://").removeprefix("http://").rstrip("/")


def default_tool_id(experiment: str, repo: str | None = None) -> str:
    """Dockstore's id shape for a tool registered from a git repository.

    Dockstore addresses a registered tool as ``#workflow/<repo path>/<name>`` for
    workflows and by image name for containers. A CommandLineTool registered through
    ``.dockstore.yml`` gets the former shape, so that is the default.
    """
    return f"#workflow/{repo or repository_slug()}/{experiment}-site-stage"


def build_tool(
    experiment: str = "fine-mapping",
    version: str | None = None,
    tool_id: str | None = None,
    image: str | None = None,
    image_digest: str | None = None,
    registry_url: str | None = None,
    organization: str = "appfl-bio-suite",
) -> tuple[Tool, dict[str, str]]:
    """Build the TRS ``Tool`` for an experiment's site stage, and its descriptor files.

    Returns ``(tool, files)`` where ``files`` maps a path to its content. The content is
    generated from the installed package -- the CWL from the site stage's real interface,
    the ``OTHER`` entries from the source of the modules APPFL actually ships -- so the
    checksums describe this install and not a document that once matched it.
    """
    from appfl_bio_suite import __version__
    from appfl_bio_suite.core.experiments import get_spec

    spec = get_spec(experiment)
    version = version or __version__
    tool_id = tool_id or default_tool_id(experiment)

    files: dict[str, str] = {
        f"{experiment}-site-stage.cwl": cwl_descriptor(experiment, image),
        # The Containerfile is part of the tool, not a side artifact: it is the only
        # statement of what the pinned image is supposed to contain.
        "Containerfile": containerfile(version),
    }
    # The two modules whose source APPFL ships to a worker. Their checksums are the
    # closest thing this suite has to "the code that ran", and putting them in the tool
    # version is what makes the TRS pin cover the computation rather than the wrapper.
    for module in spec.shipped_modules:
        path = spec.package_path / f"{module}.py"
        files[f"src/{spec.package}/{module}.py"] = path.read_text(encoding="utf-8")

    tool_files = [
        ToolFile(
            path=path,
            file_type=_file_type(path, experiment),
            checksum=[Checksum(checksum=_sha256_text(content))],
        )
        for path, content in sorted(files.items())
    ]

    # No image is registered until a coordinator has built and pushed one. Inventing a
    # plausible name here would produce a tool version that looks runnable and is not --
    # a TES task naming an image nobody published fails at the executor, minutes in,
    # rather than at the pin.
    images = (
        [
            ImageData(
                image_name=image,
                registry_host=image.split("/")[0] if "/" in image else None,
                checksum=[Checksum(checksum=image_digest)] if image_digest else [],
            )
        ]
        if image
        else []
    )

    tool = Tool(
        id=tool_id,
        url=f"{(registry_url or '').rstrip('/')}/tools/{quote(tool_id, safe='')}"
        if registry_url
        else "",
        name=f"{spec.title} -- site stage",
        description=(
            f"{spec.summary}\n\nThis tool is the SITE half: it reads one site's bundle "
            "and emits the aggregates that leave that site. It never sees another "
            "site's data and performs no fitting."
        ),
        organization=organization,
        aliases=[f"appfl-bio-suite:{experiment}-site-stage:{version}"],
        has_checker=False,
        meta_version=version,
        versions=[
            ToolVersion(
                id=version,
                name=version,
                url=f"{(registry_url or '').rstrip('/')}/tools/{quote(tool_id, safe='')}"
                f"/versions/{quote(version, safe='')}"
                if registry_url
                else "",
                descriptor_type=["CWL"],
                containerfile=True,
                images=images,
                # `verified` in TRS means a third party ran the tool's tests and attests
                # to it. Nobody has, so it is false. Setting it true because our own CI
                # passes would be asserting somebody else's opinion.
                verified=False,
                is_production=False,
                meta_version=version,
            )
        ],
    )
    tool.model_extra["files"] = [json.loads(f.model_dump_json()) for f in tool_files]
    return tool, files


def _file_type(path: str, experiment: str) -> str:
    if path.endswith(".cwl"):
        return "PRIMARY_DESCRIPTOR"
    if path == "Containerfile":
        return "CONTAINERFILE"
    return "OTHER"


def descriptor_checksum(files: dict[str, str]) -> str:
    """Merkle sha-256 over ``{path: content}``. The thing a :class:`ToolPin` freezes.

    Over paths as well as contents, so that adding a file changes the checksum. A digest
    over concatenated contents alone would not notice a file being added and another
    removed with the same total bytes -- unlikely, but this is the value the whole pin
    rests on.
    """
    digest = hashlib.sha256()
    for path in sorted(files):
        digest.update(path.encode("utf-8"))
        digest.update(b"\0")
        digest.update(_sha256_text(files[path]).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def verify_pin(pin: ToolPin, experiment: str = "fine-mapping") -> list[str]:
    """Check a pin against the installed package. Returns problems; empty means clean."""
    problems: list[str] = []
    try:
        tool, files = build_tool(
            experiment=experiment, version=pin.version, tool_id=pin.id, image=pin.image
        )
    except Exception as exc:  # noqa: BLE001 - report, never raise, from a check
        return [f"could not build the tool document for {pin.id}@{pin.version}: {exc}"]

    if pin.descriptor_checksum:
        actual = descriptor_checksum(files)
        if actual != pin.descriptor_checksum:
            problems.append(
                f"tool {pin.id}@{pin.version} does not match its pin.\n"
                f"  pinned  {pin.descriptor_checksum}\n"
                f"  actual  {actual}\n"
                "The installed site stage differs from the version this federation "
                "agreed on. Either you changed the code (re-publish and re-pin) or this "
                "checkout is not the one the pin was written for."
            )
    try:
        tool.version(pin.version)
    except TrsError as exc:
        problems.append(str(exc))
    return problems


# ---------------------------------------------------------------------------
# the generated descriptors
# ---------------------------------------------------------------------------


# Shown in a generated descriptor when no image has been published yet. A CWL runner
# refuses it, which is the correct outcome and a far better one than a pull that 404s
# after the workflow has already started.
_UNPUBLISHED_IMAGE = "IMAGE-NOT-PUBLISHED-see-Containerfile"


def cwl_descriptor(experiment: str = "fine-mapping", image: str = "") -> str:
    """A CWL v1.2 CommandLineTool wrapping ``appfl-bio-suite site-stage``.

    Inputs are a directory (the site's bundle, DRS-addressed when a TES engine stages it)
    and a JSON run config; the output is the aggregate payload the coordinator pools. The
    interface is deliberately two inputs and one output directory: anything richer would
    have to be kept in step with the site stage's Python signature by hand.
    """
    return f"""\
#!/usr/bin/env cwl-runner
# GENERATED by `appfl-bio-suite ga4gh trs publish`. Do not edit by hand: it is checksummed
# into the TRS tool version, and an edit here is a pin mismatch at every site.
cwlVersion: v1.2
class: CommandLineTool
id: {experiment}-site-stage
label: {experiment} site stage
doc: |
  Compute one site's aggregates for the federated {experiment} experiment.

  Reads the site's own bundle and emits only the aggregate payload the coordinator
  pools. No genotype row and no phenotype value is written to the outputs.

  The site's data use terms (DUO, in the bundle's DATA_USE.json) are enforced by the
  tool itself against the run config's data use request. A task whose request is not
  permitted exits non-zero having read no genotypes.

baseCommand: [appfl-bio-suite, site-stage]

requirements:
  DockerRequirement:
    dockerPull: {image or _UNPUBLISHED_IMAGE}
  NetworkAccess:
    networkAccess: false
  InlineJavascriptRequirement: {{}}

inputs:
  bundle:
    type: Directory
    doc: The site's own bundle -- the directory holding site_genotypes.*, the manifests,
      and phenotypes/. Staged from a DRS URI by the execution engine.
    inputBinding:
      prefix: --data-dir
  run_config:
    type: File
    doc: JSON run config -- client id, locus shard, ancestry columns, and the data use
      request this run declares.
    inputBinding:
      prefix: --config
  output_name:
    type: string
    default: aggregates
    inputBinding:
      prefix: --out-name

outputs:
  aggregates:
    type: File
    doc: The site aggregates. Second moments only.
    outputBinding:
      glob: $(inputs.output_name).npz
  provenance:
    type: File
    doc: GA4GH provenance -- the DRS object read, the TRS tool run, the DUO decision.
    outputBinding:
      glob: ga4gh_provenance.json

successCodes: [0]
"""


def containerfile(version: str) -> str:
    """A Containerfile for the site-stage image a TES executor runs.

    Installs the partner extra, not ``[all]``: the container runs the site half, and an
    image carrying the simulation and fitting toolchain would ship SuSiEx and PLINK to
    every site for no reason.
    """
    return f"""\
# GENERATED by `appfl-bio-suite ga4gh trs publish`.
#
# Build and push this before running the TES path; the TRS tool version names the image
# and a TES executor cannot invent it:
#
#   podman build -t <your registry>/appfl-bio-suite:{version} -f Containerfile .
#   podman push  <your registry>/appfl-bio-suite:{version}
#   appfl-bio-suite ga4gh trs publish --image <your registry>/appfl-bio-suite:{version} \\
#       --image-digest <the sha256 push reported>
#
# Pin the digest, not the tag. A tag can be repushed; a federation that agreed on a tag
# has agreed on nothing durable.
FROM python:3.12-slim

# The partner extra only. This image runs the SITE stage: it forms second moments and
# returns them. It has no reason to carry SuSiEx, PLINK, matplotlib, or the simulation
# dependencies, all of which are coordinator-side.
RUN pip install --no-cache-dir "appfl-bio-suite[finemapping]=={version}"

# Same thread caps the endpoint workers get. A container scheduled onto a busy node with
# an unbounded BLAS thread pool fails during numpy's import, before any of this runs.
ENV OMP_NUM_THREADS=1 \\
    OPENBLAS_NUM_THREADS=1 \\
    MKL_NUM_THREADS=1 \\
    NUMEXPR_NUM_THREADS=1 \\
    PYTHONNOUSERSITE=1

# No entrypoint. The TES executor and the CWL descriptor both name the command
# explicitly, and an entrypoint that quietly prepends arguments is the single most
# common reason a task runs something other than what its task document says.
CMD ["appfl-bio-suite", "--help"]
"""


def dockstore_yaml(
    experiment: str = "fine-mapping", descriptor_path: str | None = None
) -> str:
    """A ``.dockstore.yml`` registering the site stage with Dockstore.

    Dockstore is a TRS implementation, so registering here is how this tool becomes
    resolvable by TRS id to anyone outside this federation -- which is the only reason a
    tool registry beats a git tag.
    """
    descriptor_path = descriptor_path or f"/ga4gh/trs/{experiment}-site-stage.cwl"
    return f"""\
# Dockstore registration for the {experiment} site stage.
#
# GENERATED by `appfl-bio-suite ga4gh trs publish`; regenerate rather than edit.
#
# Dockstore reads this from the default branch of the GitHub repository once the
# Dockstore GitHub App is installed on it. The tool then resolves by TRS id at
# https://dockstore.org/api/ga4gh/trs/v2/tools/... , which is what makes the pin in
# federation.yaml meaningful to someone who is not in this federation.
version: 1.2

tools:
  - subclass: CWL
    name: {experiment}-site-stage
    primaryDescriptorPath: {descriptor_path}
    readMePath: /docs/coordinator/ga4gh.md
    topic: >-
      Site half of the federated {experiment} experiment: reads one institution's own
      cohort and returns aggregate statistics only.
    authors:
      - name: appfl-bio-suite contributors
    publish: true
"""


# ---------------------------------------------------------------------------
# a static TRS tree, and a client for a real one
# ---------------------------------------------------------------------------


def service_info(registry_url: str = "") -> dict[str, Any]:
    from appfl_bio_suite import __version__

    return {
        "id": "org.ga4gh.trs.appfl-bio-suite",
        "name": "appfl-bio-suite tool registry",
        "type": {"group": "org.ga4gh", "artifact": "trs", "version": TRS_VERSION},
        "description": "Static TRS index of the site-stage tools of an appfl-bio-suite federation.",
        "organization": {"name": "appfl-bio-suite", "url": registry_url},
        "version": __version__,
        "environment": "research",
    }


def write_registry(
    out_dir: str | Path,
    experiment: str = "fine-mapping",
    version: str | None = None,
    tool_id: str | None = None,
    image: str | None = None,
    image_digest: str | None = None,
    registry_url: str | None = None,
) -> tuple[Path, ToolPin]:
    """Write a servable TRS tree plus the descriptors. Returns ``(dir, pin)``.

    The tree mirrors the TRS 2.0.1 paths so that pointing any static file server at
    ``<out>/`` yields a conformant read-only registry:

        ga4gh/trs/v2/service-info
        ga4gh/trs/v2/toolClasses
        ga4gh/trs/v2/tools
        ga4gh/trs/v2/tools/<url-encoded id>
        ga4gh/trs/v2/tools/<id>/versions
        ga4gh/trs/v2/tools/<id>/versions/<version>
        ga4gh/trs/v2/tools/<id>/versions/<version>/CWL/descriptor
        ga4gh/trs/v2/tools/<id>/versions/<version>/CWL/files

    A collection and a member cannot both live at ``.../tools`` on a filesystem, so every
    collection document is written as ``<collection>/index.json``. That is one line of
    static-host configuration (nginx: ``index index.json;``) and it is what :func:`serve`
    resolves automatically, so the tree is usable immediately either way.
    """
    out_dir = Path(out_dir)
    tool, files = build_tool(
        experiment=experiment,
        version=version,
        tool_id=tool_id,
        image=image,
        image_digest=image_digest,
        registry_url=registry_url,
    )
    checksum = descriptor_checksum(files)
    tool_version = tool.versions[0]

    base = out_dir / "ga4gh" / "trs" / "v2"
    encoded = quote(tool.id, safe="")
    version_dir = base / "tools" / encoded / "versions" / quote(tool_version.id, safe="")

    def dump(path: Path, payload: Any) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    tool_json = json.loads(tool.model_dump_json(exclude_none=True))
    version_json = json.loads(tool_version.model_dump_json())
    dump(base / "service-info" / "index.json", service_info(registry_url or ""))
    dump(base / "toolClasses" / "index.json", [COMMAND_LINE_TOOL])
    dump(base / "tools" / "index.json", [tool_json])
    dump(base / "tools" / encoded / "index.json", tool_json)
    dump(base / "tools" / encoded / "versions" / "index.json", [version_json])
    dump(version_dir / "index.json", version_json)
    dump(version_dir / "CWL" / "files" / "index.json", tool_json.get("files", []))

    primary = f"{experiment}-site-stage.cwl"
    dump(
        version_dir / "CWL" / "descriptor" / "index.json",
        json.loads(
            FileWrapper(
                content=files[primary],
                checksum=[Checksum(checksum=_sha256_text(files[primary]))],
                url=f"{(registry_url or '').rstrip('/')}/{primary}" if registry_url else None,
            ).model_dump_json(exclude_none=True)
        ),
    )
    for path, content in files.items():
        target = out_dir / "descriptors" / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(content, encoding="utf-8")

    # The pin a coordinator copies into federation.yaml. Emitted as a file rather than
    # only printed, because a checksum retyped by hand is a checksum entered wrong.
    pin = ToolPin(
        id=tool.id,
        version=tool_version.id,
        descriptor_checksum=checksum,
        image=tool_version.images[0].image_name if tool_version.images else None,
        registry_url=registry_url,
    )
    dump(out_dir / "tool_pin.json", json.loads(pin.model_dump_json(exclude_none=True)))
    (out_dir / ".dockstore.yml").write_text(dockstore_yaml(experiment), encoding="utf-8")
    return out_dir, pin


class TrsClient:
    """Read a TRS registry -- Dockstore, a static tree served over HTTP, or a local one.

    Standard library only, for the same reason the DRS client is: whatever can run on a
    worker must not require a package a partner did not install.
    """

    def __init__(self, base_url: str, token: str | None = None, timeout: float = 30.0):
        self.base_url = base_url.rstrip("/")
        self.token = token
        self.timeout = timeout

    def _get(self, path: str) -> Any:
        import urllib.error
        import urllib.request

        url = f"{self.base_url}{path}"
        request = urllib.request.Request(url, headers={"Accept": "application/json"})
        if self.token:
            request.add_header("Authorization", f"Bearer {self.token}")
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", "replace")[:400]
            raise TrsError(f"{url} returned {exc.code}: {body}") from exc
        except urllib.error.URLError as exc:
            raise TrsError(f"could not reach {url}: {exc.reason}") from exc

    def service_info(self) -> dict[str, Any]:
        return self._get("/ga4gh/trs/v2/service-info")

    def tool(self, tool_id: str) -> Tool:
        return Tool.model_validate(self._get(f"/ga4gh/trs/v2/tools/{quote(tool_id, safe='')}"))

    def version(self, tool_id: str, version_id: str) -> ToolVersion:
        return ToolVersion.model_validate(
            self._get(
                f"/ga4gh/trs/v2/tools/{quote(tool_id, safe='')}"
                f"/versions/{quote(version_id, safe='')}"
            )
        )

    def descriptor(
        self, tool_id: str, version_id: str, descriptor_type: str = "CWL"
    ) -> FileWrapper:
        return FileWrapper.model_validate(
            self._get(
                f"/ga4gh/trs/v2/tools/{quote(tool_id, safe='')}"
                f"/versions/{quote(version_id, safe='')}/{descriptor_type}/descriptor"
            )
        )

    def files(self, tool_id: str, version_id: str, descriptor_type: str = "CWL") -> list[ToolFile]:
        payload = self._get(
            f"/ga4gh/trs/v2/tools/{quote(tool_id, safe='')}"
            f"/versions/{quote(version_id, safe='')}/{descriptor_type}/files"
        )
        return [ToolFile.model_validate(entry) for entry in payload]

    def image_for(self, tool_id: str, version_id: str) -> str:
        """The container image a TES executor should run for this tool version."""
        version = self.version(tool_id, version_id)
        if not version.images:
            raise TrsError(
                f"{tool_id}@{version_id} registers no container image, so there is "
                "nothing for a TES executor to run. Publish one: see the generated "
                "Containerfile."
            )
        image = version.images[0]
        digest = image.digest
        return f"{image.image_name}@sha256:{digest}" if digest else image.image_name


def serve(tree: str | Path, host: str = "127.0.0.1", port: int = 8080):
    """Serve a tree written by :func:`write_registry` as a read-only TRS. Caller starts it.

    Exists so the generated tree is usable without configuring a web server: it resolves
    a spec URL to the file on disk, falling back to that path's ``index.json`` for
    collection endpoints. Returned unstarted for the same reason the DRS server is --
    a test binds port 0 and reads back what it got.
    """
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
    from urllib.parse import unquote, urlparse

    root = Path(tree).resolve()

    class Handler(BaseHTTPRequestHandler):
        server_version = "appfl-bio-suite-trs"

        def do_GET(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler's interface
            path = unquote(urlparse(self.path).path).strip("/")
            # A TRS id contains '/' and is URL-encoded in the request; unquoting it back
            # into the filesystem path is exactly right here, because write_registry
            # encoded it into a single directory name.
            candidate = root / urlparse(self.path).path.strip("/")
            for target in (candidate, candidate / "index.json"):
                if target.is_file() and root in target.resolve().parents:
                    body = target.read_bytes()
                    self.send_response(200)
                    self.send_header("Content-Type", "application/json")
                    self.send_header("Content-Length", str(len(body)))
                    self.end_headers()
                    self.wfile.write(body)
                    return
            body = json.dumps({"code": 404, "message": f"no route /{path}"}).encode("utf-8")
            self.send_response(404)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
            return

    return ThreadingHTTPServer((host, port), Handler)
