"""Check that a run will work, before starting it.

The design rule here is the one thing worth reading before adding a check:

    **A check fails hard only if it would be wrong everywhere.**

Version skew across the federation is wrong everywhere -- that is a hard failure. Running
on the wrong login node is wrong *on one particular cluster*, so it is a warning, and one
the coordinator opts into via ``coordinator.host_check`` in federation.yaml.

That distinction is what keeps this suite usable by someone who is not its author. The
tooling it grew out of hard-failed on a hostname pattern belonging to one specific HPC
center; on anyone else's machine that check is not a safety net, it is a wall. Every
cluster-specific fact in this project is now either configuration or documentation.

What it preserves from the original tooling is the *spirit*: fail early, fail loudly, and
name the fix in the message rather than making someone go read a document.
"""

from __future__ import annotations

import os
import platform
import re
import socket
import sys
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path

__all__ = ["Level", "Check", "PreflightReport", "run_preflight", "CHECK_GROUPS"]

CHECK_GROUPS = ("env", "pins", "config", "configs", "data", "ga4gh", "endpoints", "all")

# Packages whose version must be identical at every site or payloads stop deserializing.
# These come from constraints.txt at runtime; the list here is what to look at.
FEDERATION_CRITICAL = (
    "globus-compute-sdk",
    "globus-compute-endpoint",
    "globus-sdk",
    "dill",
    "appfl",
    "torch",
    "numpy",
)


class Level(StrEnum):
    OK = "ok"
    WARN = "warn"
    FAIL = "fail"
    SKIP = "skip"


@dataclass
class Check:
    name: str
    level: Level
    detail: str = ""
    fix: str = ""

    def render(self) -> str:
        mark = {
            Level.OK: "  ok  ",
            Level.WARN: " warn ",
            Level.FAIL: " FAIL ",
            Level.SKIP: " skip ",
        }
        out = [f"[{mark[self.level]}] {self.name}"]
        if self.detail:
            out.extend(f"           {line}" for line in self.detail.splitlines())
        if self.fix and self.level in (Level.WARN, Level.FAIL):
            out.append("           ->")
            out.extend(f"           -> {line}" for line in self.fix.splitlines())
        return "\n".join(out)


@dataclass
class PreflightReport:
    checks: list[Check] = field(default_factory=list)

    def add(self, name: str, level: Level, detail: str = "", fix: str = "") -> None:
        self.checks.append(Check(name=name, level=level, detail=detail, fix=fix))

    @property
    def failed(self) -> list[Check]:
        return [c for c in self.checks if c.level is Level.FAIL]

    @property
    def warned(self) -> list[Check]:
        return [c for c in self.checks if c.level is Level.WARN]

    @property
    def ok(self) -> bool:
        return not self.failed

    def render(self) -> str:
        lines = [c.render() for c in self.checks]
        lines.append("")
        if self.failed:
            lines.append(f"{len(self.failed)} check(s) FAILED. Fix these before launching.")
        elif self.warned:
            lines.append(
                f"All checks passed with {len(self.warned)} warning(s). "
                "Warnings are things that are wrong on some clusters and fine on others -- "
                "read them, then proceed if they do not apply to yours."
            )
        else:
            lines.append("All checks passed.")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# environment
# ---------------------------------------------------------------------------


def _check_python(report: PreflightReport) -> None:
    major, minor, *_ = sys.version_info
    version = platform.python_version()
    if (major, minor) == (3, 12):
        report.add(
            "python minor version",
            Level.OK,
            f"{version} (federation requires 3.12.x; patch level need not match)",
        )
    else:
        report.add(
            "python minor version",
            Level.FAIL,
            f"running {version}; the federation is standardized on 3.12.x",
            "Dill ships function bytecode, and bytecode is only portable within a minor\n"
            "version. A 3.12 driver against 3.10 workers fails at decode, before any of\n"
            "your code runs (SystemError: unknown opcode).\n"
            "Create the environment with: conda create -n appfl_env python=3.12",
        )


def _check_user_site(report: PreflightReport) -> None:
    """A ~/.local install can shadow the environment's own packages, invisibly.

    This is nastier than it sounds: `pip list` reports the same version either way, so
    the only way to tell which copy is live is to look at `module.__file__`. It has
    produced a silent driver/worker version skew on this project, where the driver used
    a package from ~/.local and the workers -- whose worker_init sets PYTHONNOUSERSITE --
    used the environment's.
    """
    if os.environ.get("PYTHONNOUSERSITE"):
        report.add("user-site packages", Level.OK, "PYTHONNOUSERSITE is set")
        return

    import site

    user_site = getattr(site, "getusersitepackages", lambda: None)()
    if not user_site or not Path(str(user_site)).is_dir():
        report.add("user-site packages", Level.OK, "no user-site directory present")
        return

    shadowed = [
        name
        for name in ("torch", "numpy", "globus_compute_sdk", "appfl")
        if _resolves_under(name, str(user_site))
    ]

    if shadowed:
        report.add(
            "user-site packages",
            Level.FAIL,
            f"these are loading from {user_site}: {', '.join(shadowed)}",
            "A ~/.local copy is shadowing the environment's version. Workers run with\n"
            "PYTHONNOUSERSITE set and will use a DIFFERENT copy -- a silent skew that\n"
            "`pip list` cannot reveal, because it reports the same version either way.\n"
            "Check with: python -c 'import torch; print(torch.__file__)'\n"
            "Fix with:   PYTHONNOUSERSITE=1 python -m pip install --no-deps \\\n"
            "                --ignore-installed <package>",
        )
    else:
        report.add(
            "user-site packages",
            Level.WARN,
            f"user-site is enabled ({user_site}) but nothing critical is shadowed",
            "Set PYTHONNOUSERSITE=1 to keep it that way. A later `pip install --user`\n"
            "would otherwise shadow an environment package without changing `pip list`.",
        )


def _resolves_under(module: str, prefix: str) -> bool:
    """Would importing ``module`` load it from under ``prefix``?

    Resolved with ``find_spec`` rather than by inspecting ``sys.modules``. Consulting
    ``sys.modules`` alone only sees what something has already imported, and preflight
    runs before anything has touched torch, numpy or appfl -- so for those three the check
    reported "nothing shadowed" no matter what was installed in ~/.local. ``find_spec`` on
    a TOP-LEVEL name answers what the next import would resolve to, without executing it.
    """
    import importlib.util

    loaded = sys.modules.get(module)
    origin = getattr(loaded, "__file__", None) if loaded is not None else None
    if origin is None:
        try:
            spec = importlib.util.find_spec(module)
        except (ImportError, AttributeError, ValueError):
            return False
        origin = getattr(spec, "origin", None) if spec is not None else None
        if origin is None and spec is not None and spec.submodule_search_locations:
            origin = next(iter(spec.submodule_search_locations), None)
    return bool(origin) and str(origin).startswith(prefix)


def _check_thread_caps(report: PreflightReport) -> None:
    """OpenBLAS/OpenMP spawn one thread per core on every numpy import.

    On a busy shared login node that is enough to exhaust the per-user thread budget,
    and the failure is not obviously about threads: 'RuntimeError: can't start new
    thread', or the interchange simply failing to start.
    """
    caps = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
    missing = [c for c in caps if os.environ.get(c) != "1"]
    if not missing:
        report.add("thread caps", Level.OK, "OPENBLAS/OMP/MKL_NUM_THREADS=1")
    else:
        report.add(
            "thread caps",
            Level.WARN,
            f"not set to 1: {', '.join(missing)}",
            "On a shared login node, numpy's import-time thread fan-out can exhaust the\n"
            "per-user thread limit, surfacing as 'RuntimeError: can't start new thread'.\n"
            "Export before launching:\n"
            "    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\n"
            "The same three lines belong in every endpoint's worker_init.",
        )


def _check_tmpdir(report: PreflightReport) -> None:
    """A long $TMPDIR overflows the AF_UNIX socket path limit and workers never register."""
    tmpdir = os.environ.get("TMPDIR", "/tmp")
    # The limit is ~108 bytes for the whole socket path; parsl appends a good deal to
    # $TMPDIR, so anything long is a problem well before the raw limit.
    if len(tmpdir) > 60:
        report.add(
            "TMPDIR length",
            Level.WARN,
            f"TMPDIR is {len(tmpdir)} chars: {tmpdir}",
            "parsl's worker pool opens a SyncManager AF_UNIX socket under $TMPDIR, and\n"
            "the path limit is ~108 characters. Schedulers often set $TMPDIR to a long\n"
            "path, which surfaces as 'workers failed to register' or\n"
            "'OSError: AF_UNIX path too long'.\n"
            "Add `export TMPDIR=/tmp` to worker_init on any endpoint that hits this.",
        )
    else:
        report.add("TMPDIR length", Level.OK, f"{tmpdir} ({len(tmpdir)} chars)")


def _check_binaries(report: PreflightReport, experiment: str | None) -> None:
    """Executables an experiment shells out to, which pip cannot install.

    Only fine-mapping has any. It is worth a check rather than a stack trace because of
    *when* the failure lands: the SuSiEx binary is not touched until after every site has
    computed and transferred its second moments, which is the expensive part of the run.
    Discovering it is missing at that point wastes the whole exchange.

    Hard only when an experiment was named -- which includes every launch, since
    ``launch()`` preflights the experiment it is about to run. A bare ``preflight`` scans
    every experiment, and there a missing SuSiEx means "you could not run fine-mapping",
    not "your environment is broken": a coordinator running only the GWAS experiment is
    entitled to a clean report. Same distinction the host check draws.
    """
    from appfl_bio_suite.core.experiments import REGISTRY

    targeted = bool(experiment and experiment in REGISTRY)
    specs = [REGISTRY[experiment]] if targeted else [s for s in REGISTRY.values() if s.implemented]
    wanted = {name: spec.name for spec in specs for name in spec.required_binaries}
    if not wanted:
        report.add("external binaries", Level.SKIP, "no experiment here needs one")
        return

    import shutil as _shutil

    missing = []
    found = []
    for binary, owner in sorted(wanted.items()):
        path = _shutil.which(binary) or _vendored(binary)
        if path:
            found.append(f"{binary}: {path}")
        else:
            missing.append(f"{binary} (needed by {owner})")

    if missing:
        report.add(
            "external binaries",
            Level.FAIL if targeted else Level.WARN,
            "not found on PATH or in vendor/bin: "
            + ", ".join(missing)
            + ("" if targeted else "\n(only matters if you run that experiment)"),
            "These are C++ command-line tools the experiment shells out to; they are not\n"
            "pip-installable and no extra provides them. Vendor static builds with:\n"
            "    ./scripts/fine-mapping/install_susiex.sh\n"
            "    ./scripts/fine-mapping/install_plink.sh\n"
            "Both need outbound network access once, and write into vendor/bin/.\n"
            "They are COORDINATOR-side only -- no partner is asked to install either.",
        )
    else:
        report.add("external binaries", Level.OK, "\n".join(found))


def _vendored(binary: str) -> str | None:
    """``vendor/bin/<binary>`` relative to the working directory, if it is executable."""
    candidate = Path("vendor") / "bin" / binary
    return str(candidate.resolve()) if candidate.is_file() else None


def _check_host(report: PreflightReport, federation) -> None:
    """Optional, coordinator-declared, and a warning by default.

    Some clusters really do have login nodes that cannot run the driver -- a shared
    filesystem makes the wrong host look correct right up until an import fails. But
    that is a fact about one deployment, not about federated learning, so the pattern
    comes from federation.yaml and defaults to advisory.
    """
    pattern = getattr(federation.coordinator, "host_check", None) if federation else None
    if not pattern:
        report.add("driver host", Level.SKIP, "no coordinator.host_check configured")
        return

    host = socket.gethostname()
    if re.search(pattern, host):
        report.add("driver host", Level.OK, f"{host} matches {pattern!r}")
        return

    enforce = getattr(federation.coordinator, "host_check_enforce", False)
    report.add(
        "driver host",
        Level.FAIL if enforce else Level.WARN,
        f"hostname {host!r} does not match coordinator.host_check {pattern!r}",
        "Your federation config restricts which host the driver may run on. If this\n"
        "host is fine, update or remove `coordinator.host_check` in federation.yaml.\n"
        "If it is not, move to a matching host before launching -- on some clusters a\n"
        "shared filesystem makes the wrong node look correct until an import fails.",
    )


# ---------------------------------------------------------------------------
# pins
# ---------------------------------------------------------------------------


def _parse_constraints(path: Path) -> dict[str, str]:
    pins: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if not line or "==" not in line:
            continue
        name, _, version = line.partition("==")
        pins[name.strip().lower().replace("_", "-")] = version.strip()
    return pins


def _installed_versions() -> dict[str, str]:
    from importlib.metadata import distributions

    out = {}
    for dist in distributions():
        # `metadata` is None for a dist-info directory with no METADATA file, which is
        # what a partially-removed package leaves behind. Reading `.get` off it took down
        # the whole pin check with an AttributeError that named neither the package nor
        # the cause -- on exactly the kind of half-cleaned environment the check is for.
        metadata = dist.metadata
        if metadata is None:
            continue
        name = (metadata.get("Name") or "").lower().replace("_", "-")
        if name:
            out[name] = dist.version
    return out


def _check_pins(report: PreflightReport, constraints: Path | None) -> None:
    if constraints is None or not constraints.is_file():
        report.add(
            "version pins",
            Level.SKIP,
            "constraints.txt not found (are you running from a source checkout?)",
        )
        return

    pinned = _parse_constraints(constraints)
    installed = _installed_versions()

    mismatched, missing = [], []
    for name, want in sorted(pinned.items()):
        have = installed.get(name)
        if have is None:
            # Only complain about something absent if it is federation-critical; the
            # extras are legitimately not all installed.
            if name in FEDERATION_CRITICAL:
                missing.append(name)
        elif have != want:
            mismatched.append(f"{name}: pinned {want}, installed {have}")

    if not mismatched and not missing:
        report.add("version pins", Level.OK, f"{len(pinned)} pins match constraints.txt")
        return

    detail = []
    if mismatched:
        detail.append("version mismatch:")
        detail.extend(f"  {m}" for m in mismatched)
    if missing:
        detail.append("not installed:")
        detail.extend(f"  {m}" for m in missing)

    critical = [m for m in mismatched if m.split(":")[0] in FEDERATION_CRITICAL] or missing
    report.add(
        "version pins",
        Level.FAIL if critical else Level.WARN,
        "\n".join(detail),
        "Every site in the federation must run the same versions of these. Skew does\n"
        "not fail at install time -- it fails as a deserialization error partway into a\n"
        "run on someone else's cluster.\n"
        "Reinstall with the constraints applied:\n"
        "    pip install -e '.[all]' -c constraints.txt\n"
        "If a utility moved a pin (globus-cli is the usual culprit), install one-off\n"
        "tools in a throwaway venv instead of this environment.",
    )


def _check_appfl_importable(report: PreflightReport) -> None:
    """The version pair this suite ships needs a shim; confirm it is doing its job."""
    from appfl_bio_suite.core import compat

    try:
        import importlib

        importlib.import_module("appfl.comm.globus_compute")
    except Exception as exc:
        report.add(
            "appfl globus-compute driver imports",
            Level.FAIL,
            f"{type(exc).__name__}: {exc}",
            "This is the production dispatch path for every experiment in the suite.\n"
            "If the error names globus_compute_sdk.sdk.login_manager, the compat shim\n"
            "did not apply -- import appfl_bio_suite before importing appfl.\n"
            "See src/appfl_bio_suite/core/compat.py.",
        )
        return

    note = "compat shim active" if compat.shim_is_needed() else "no shim needed"
    report.add("appfl globus-compute driver imports", Level.OK, note)


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------


def _check_data(report: PreflightReport, federation, experiment: str | None) -> None:
    """Check per-site data, for the sites whose data this machine can actually see.

    Most of the time that is none of them, and saying so is the point. ``data_dir`` is an
    absolute path on a partner's cluster; from here it is not merely absent, it is
    unknowable. A check that treated "not found locally" as a failure would fail on every
    correctly-configured federation -- and one that quietly passed would be worse. This
    group used to be declared, run nothing, and report "All checks passed", which is the
    most misleading thing a preflight can do.

    So each site gets an honest verdict: validated where the path resolves here (a
    loopback run, a coordinator who also hosts a site, a shared filesystem), and skipped
    with a pointer to the partner-side verification where it does not.
    """
    names = [experiment] if experiment else list(federation.enabled_experiments())
    for name in names:
        try:
            exp = federation.experiment(name)
        except Exception as exc:  # noqa: BLE001 - reported, not raised
            report.add(f"data [{name}]", Level.FAIL, str(exc))
            continue

        if not exp.sites:
            report.add(f"data [{name}]", Level.SKIP, "no sites declared for this experiment")
            continue

        required = _required_input_files(name)
        if required is None:
            # FLamby partners download a public dataset into FLamby's own configured
            # location, so nothing in federation.yaml names a path to check.
            report.add(
                f"data [{name}]",
                Level.SKIP,
                "partners obtain this dataset themselves; nothing here names a path",
                "Each site confirms their own copy with the verification step in "
                f"docs/partner/experiments/{name}.md.",
            )
            continue

        for entry in exp.sites:
            _check_site_data(report, name, entry, *required)


def _required_input_files(experiment: str):
    """``(REQUIRED_FILES, PLINK_STEM)`` for an experiment, or None if it names no path.

    Read off the experiment's own shipped loader rather than listed here. That module is
    what actually fails when a file is missing, so taking the list from anywhere else
    guarantees the two drift and the preflight passes on a bundle the loader will reject.
    """
    from appfl_bio_suite.core.experiments import REGISTRY

    spec = REGISTRY.get(experiment)
    if spec is None or not spec.implemented:
        return None

    import importlib

    try:
        loader = importlib.import_module(f"appfl_bio_suite.experiments.{spec.package}.dataset")
    except ImportError:
        return None

    required = getattr(loader, "REQUIRED_FILES", None)
    if not required:
        return None
    return required, getattr(loader, "PLINK_STEM", None)


def _check_site_data(
    report: PreflightReport, experiment: str, entry, required_files, plink_stem
) -> None:
    label = f"data {experiment}/{entry.client_id}"
    if not entry.data_dir:
        report.add(label, Level.FAIL, "no data_dir declared for this site")
        return

    data_dir = Path(entry.data_dir)
    if not data_dir.is_dir():
        report.add(
            label,
            Level.SKIP,
            f"{data_dir} is not visible from here",
            "Expected -- this path is on the partner's cluster. Have them confirm the\n"
            "bundle unpacked correctly using the verification step in their setup guide.",
        )
        return

    missing = [f for f in required_files if not (data_dir / f).is_file()]
    if missing:
        report.add(
            label,
            Level.FAIL,
            f"{data_dir} is missing: {', '.join(missing)}",
            "The required files must sit DIRECTLY in data_dir, not in a subdirectory.\n"
            "A bundle unpacked one level too deep is the usual cause.",
        )
        return

    if plink_stem is None:
        report.add(label, Level.OK, f"{data_dir}\nall required files present")
        return

    fam = data_dir / f"{plink_stem}.fam"
    try:
        n_samples = sum(1 for line in fam.read_text(encoding="utf-8").splitlines() if line.strip())
    except OSError as exc:
        report.add(label, Level.WARN, f"{data_dir}\ncould not read {fam.name}: {exc}")
        return

    detail = f"{data_dir}\n{n_samples:,} samples"
    if entry.expected_samples is not None and n_samples != entry.expected_samples:
        report.add(
            label,
            Level.FAIL,
            f"{detail}, but federation.yaml declares {entry.expected_samples:,}",
            "This site holds a different dataset from the one you recorded. Usually a\n"
            "bundle from an earlier simulation run, or two sites' bundles swapped.\n"
            "Compare their run_manifest.json checksums before launching.",
        )
        return

    report.add(label, Level.OK, detail)


# ---------------------------------------------------------------------------
# endpoints
# ---------------------------------------------------------------------------


def _check_endpoints(report: PreflightReport, federation, experiment: str | None) -> None:
    from appfl_bio_suite.core.endpoint import get_status

    names = [experiment] if experiment else list(federation.enabled_experiments())
    for name in names:
        try:
            exp = federation.experiment(name)
        except Exception as exc:
            report.add(f"endpoints [{name}]", Level.FAIL, str(exc))
            continue

        if not exp.sites:
            report.add(
                f"endpoints [{name}]",
                Level.SKIP,
                "no sites declared for this experiment yet",
            )
            continue

        for entry in exp.sites:
            status = get_status(entry.endpoint_uuid, label=entry.client_id)
            if status.online:
                report.add(f"endpoint {name}/{entry.client_id}", Level.OK, status.uuid)
            else:
                report.add(
                    f"endpoint {name}/{entry.client_id}",
                    Level.FAIL,
                    f"{status.uuid} reports '{status.status}'"
                    + (f" ({status.error})" if status.error else ""),
                    "The cloud status is authoritative -- local `globus-compute-endpoint\n"
                    "list` reads a pidfile that may record a PID from another login node.\n"
                    "If the partner believes it is running, have them confirm which\n"
                    "account started it and re-read the UUID: endpoint IDs are\n"
                    "per-account, so an endpoint restarted under a different user is a\n"
                    "different endpoint.",
                )


# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# GA4GH
#
# Offline, all of it. The point of a preflight check is to be cheaper than the failure it
# prevents, and every one of these prevents a failure that costs a scheduler queue wait --
# a site refusing on consent grounds, a bundle that is not the one this run is about, a
# tool version that is not the one the results will claim.
#
# TES reachability is deliberately NOT here: it needs the network, so it belongs with the
# endpoint checks, where the contract is already "this group talks to things".
# ---------------------------------------------------------------------------


def _check_ga4gh(report: PreflightReport, federation, experiment: str | None) -> None:
    names = [experiment] if experiment else list(federation.enabled_experiments())
    configured = [
        name
        for name in names
        if getattr(federation.experiments.get(name), "ga4gh", None) is not None
    ]
    if not configured:
        report.add(
            "ga4gh",
            Level.SKIP,
            "no experiment declares a `ga4gh` block",
            "The four standards are opt-in. See docs/coordinator/ga4gh.md for what each\n"
            "one buys and what it costs to adopt.",
        )
        return

    _check_duo_snapshot(report)
    for name in configured:
        _check_data_use(report, federation, name)
        _check_drs(report, federation, name)
        _check_trs(report, federation, name)


def _check_duo_snapshot(report: PreflightReport) -> None:
    """The vendored ontology must load. Everything else here depends on it."""
    from appfl_bio_suite.core.ga4gh.duo import DuoError, ontology

    try:
        snapshot = ontology()
    except DuoError as exc:
        report.add("duo ontology", Level.FAIL, str(exc))
        return
    report.add(
        "duo ontology",
        Level.OK,
        f"{len(snapshot['terms'])} terms from {snapshot['version_iri']}",
    )


def _check_data_use(report: PreflightReport, federation, experiment: str) -> None:
    from appfl_bio_suite.core.ga4gh.resolve import resolve_data_use

    block = federation.experiment(experiment).ga4gh
    request = federation.data_use_request(experiment)
    name = f"data use [{experiment}]"

    if request is None:
        report.add(
            name,
            Level.WARN,
            "no `ga4gh.data_use_request` declared",
            "Any site whose bundle carries a DATA_USE.json will refuse the task: a\n"
            "dataset with declared terms cannot be used by a study that declares\n"
            "nothing about itself. Declare the study's purposes under\n"
            f"experiments.{experiment}.ga4gh.data_use_request.",
        )
        return

    decisions = resolve_data_use(federation, experiment)
    blocking = [d for d in decisions if d.blocking]
    unchecked = [d for d in decisions if d.status == "no-profile"]
    detail = "\n".join(d.render() for d in decisions)

    if blocking:
        report.add(
            name,
            Level.FAIL if block.enforce_data_use else Level.WARN,
            detail,
            "Those sites' terms do not permit this study. Either the declared purposes\n"
            "are wrong for what you are actually doing, or those sites should not be in\n"
            "this run. An 'undetermined' is a term a program cannot evaluate -- read it,\n"
            "then record that you did by adding the term id to `acknowledged`.",
        )
        return
    if unchecked:
        report.add(
            name,
            Level.WARN,
            detail,
            f"{len(unchecked)} site(s) have no `data_use_profile` recorded here, so their\n"
            "terms cannot be checked before dispatch. Ask them for their DATA_USE.json;\n"
            "their own worker enforces it either way.",
        )
        return
    report.add(name, Level.OK, detail)


def _check_drs(report: PreflightReport, federation, experiment: str) -> None:
    from appfl_bio_suite.core.ga4gh.drs import DrsError
    from appfl_bio_suite.core.ga4gh.resolve import load_registry

    name = f"drs [{experiment}]"
    sites = federation.experiment(experiment).sites
    with_uri = [s for s in sites if s.drs_uri]

    if not with_uri:
        report.add(
            name,
            Level.SKIP,
            "no site names a `drs_uri`",
            "Without one, a site's data is identified by a filesystem path and nothing\n"
            "checks that it holds the bundle you cut for it. `simulate` writes a\n"
            "registry; `ga4gh drs register` builds one from an existing directory.",
        )
        return

    try:
        registry = load_registry(federation)
    except (DrsError, FileNotFoundError) as exc:
        report.add(name, Level.FAIL, str(exc))
        return
    if registry is None:
        report.add(
            name,
            Level.FAIL,
            f"{len(with_uri)} site(s) name a drs_uri but no `ga4gh.drs.registry` is set",
        )
        return

    problems = []
    for entry in with_uri:
        try:
            obj = registry.resolve(entry.drs_uri)
        except DrsError as exc:
            problems.append(f"{entry.client_id}: {exc}")
            continue
        if not obj.is_bundle:
            problems.append(
                f"{entry.client_id}: {entry.drs_uri} is a single file, not a bundle. A "
                "site's data_dir is a directory of several files."
            )

    if problems:
        report.add(
            name,
            Level.FAIL,
            "\n".join(problems),
            "Rebuild the registry against the data you actually distributed:\n"
            "    appfl-bio-suite ga4gh drs register --data-root <simulation output>",
        )
        return

    missing = [s.client_id for s in sites if not s.drs_uri]
    level = Level.WARN if missing else Level.OK
    detail = f"{len(with_uri)} bundle(s) resolve in {registry.hostname}"
    if missing:
        detail += f"; no drs_uri for {', '.join(missing)}"
    report.add(name, level, detail)


def _check_trs(report: PreflightReport, federation, experiment: str) -> None:
    from appfl_bio_suite.core.ga4gh.resolve import resolve_tool

    name = f"trs pin [{experiment}]"
    pin = federation.tool_pin(experiment)
    if pin is None:
        report.add(
            name,
            Level.SKIP,
            "no `ga4gh.tool` pin declared",
            "Without a pin, the code a site runs is whatever this checkout happens to\n"
            "contain, and the results record no version. `ga4gh trs publish` emits one.",
        )
        return

    _, problems = resolve_tool(federation, experiment)
    if problems:
        report.add(
            name,
            Level.FAIL,
            "\n".join(problems),
            "The installed site stage is not the version this federation pinned. If you\n"
            "changed it deliberately, re-publish and re-pin:\n"
            "    appfl-bio-suite ga4gh trs publish --out local/trs",
        )
        return
    if not pin.descriptor_checksum:
        report.add(
            name,
            Level.WARN,
            f"{pin.render()} -- pinned by name and version only",
            "A pin without a descriptor checksum pins a label. Freeze it:\n"
            "    appfl-bio-suite ga4gh trs publish --out local/trs\n"
            "then copy `descriptor_checksum` from local/trs/tool_pin.json.",
        )
        return
    report.add(name, Level.OK, f"{pin.render()} matches this install")


def _check_tes(report: PreflightReport, federation, experiment: str | None) -> None:
    """Reachability only. Runs with the endpoint checks, because it uses the network."""
    from appfl_bio_suite.core.ga4gh.tes import TesClient, TesError

    service = federation.ga4gh.tes if federation.ga4gh else None
    if service is None or not service.url:
        return
    try:
        info = TesClient(service.url).service_info()
    except TesError as exc:
        report.add(
            "tes service",
            Level.WARN,
            str(exc).splitlines()[0],
            "Only the TES execution path needs this; a Globus Compute run is unaffected.",
        )
        return
    kind = info.get("type", {})
    report.add(
        "tes service",
        Level.OK,
        f"{service.url}: {info.get('name', '?')} "
        f"({kind.get('artifact', '?')} {kind.get('version', '?')})",
    )


def run_preflight(
    federation=None,
    experiment: str | None = None,
    check: str = "all",
    constraints: Path | None = None,
    repo_root: Path | None = None,
) -> PreflightReport:
    """Run the requested checks and return a report.

    ``check`` selects a group so that fast, offline checks can run in CI without
    contacting any endpoint.
    """
    report = PreflightReport()
    root = repo_root or Path.cwd()
    if constraints is None:
        candidate = root / "constraints.txt"
        constraints = candidate if candidate.is_file() else None

    want = (
        {check}
        if check != "all"
        else {"env", "pins", "config", "data", "ga4gh", "endpoints"}
    )

    if "env" in want:
        _check_python(report)
        _check_user_site(report)
        _check_thread_caps(report)
        _check_tmpdir(report)
        _check_binaries(report, experiment)
        if federation is not None:
            _check_host(report, federation)

    if "pins" in want:
        _check_pins(report, constraints)
        _check_appfl_importable(report)

    if {"config", "configs"} & want:
        from appfl_bio_suite.core.experiments import check_shipped_configs

        check_shipped_configs(report, experiment)

    if "data" in want:
        if federation is not None:
            _check_data(report, federation, experiment)
        else:
            report.add("data", Level.SKIP, "no federation config loaded")

    if "ga4gh" in want:
        if federation is not None:
            _check_ga4gh(report, federation, experiment)
        else:
            report.add("ga4gh", Level.SKIP, "no federation config loaded")

    if "endpoints" in want and federation is not None:
        _check_endpoints(report, federation, experiment)
        _check_tes(report, federation, experiment)
    elif "endpoints" in want:
        report.add("endpoints", Level.SKIP, "no federation config loaded")

    return report
