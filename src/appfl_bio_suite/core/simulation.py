"""The shared simulation contract: orchestration, provenance, and determinism.

WHAT BELONGS HERE, AND WHAT DOES NOT
------------------------------------
This module holds only what generalizes across experiments: the pipeline shape
(``simulate -> partition -> materialize``), the site-manifest schema, the run manifest and
its provenance recorder, seed management, and the determinism guarantees.

The science stays in the experiment. Phenotype simulation, polygenic scoring and
PLINK-level bundling live in ``experiments/gwas/simulation/``; credible-set simulation
will live in the fine-mapping package. Those two share a *shape*, not an implementation,
and inventing a common abstraction over domain code that exists in one instance would be
guessing at the second one's requirements.

SIMULATION IS COORDINATOR-SIDE ONLY
-----------------------------------
Everything else in this suite eventually reaches a partner cluster. Simulation never
does. It runs on the coordinator's hardware and produces per-site bundles that are then
distributed. Two consequences:

* The shipped-module import rules do not apply here. Simulation code may freely import
  from ``appfl_bio_suite``, use heavy dependencies, and assume a full environment.
* Its dependencies must never land on a partner. That is why they live in a separate
  ``gwas-sim`` extra rather than in the ``gwas`` extra a partner installs.

NOT EVERY EXPERIMENT SIMULATES
------------------------------
FLamby partners download a public dataset; there is no simulation stage at all. The
contract therefore has to tolerate an experiment with no simulation, rather than assuming
every experiment has one and giving FLamby an empty stub to satisfy an interface.

WHY THE RUN MANIFEST IS NOT OPTIONAL
------------------------------------
The pipeline this replaces emitted a summary CSV of diagnostics but no record of *how*
the outputs were produced: no seeds, no input checksums, no package versions, no commit.
For a repository that supports a published result, that means the numbers cannot be
regenerated with confidence -- not because the code is wrong, but because nobody can
prove which code and which inputs produced them. Recording it costs a few hundred
milliseconds at the end of a run that takes minutes.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Protocol

__all__ = [
    "SiteAllocation",
    "SimulationScenario",
    "RunManifest",
    "SimulationPipeline",
    "file_checksum",
    "checksum_tree",
    "record_provenance",
    "verify_against_manifest",
    "MANIFEST_FILENAME",
]

MANIFEST_FILENAME = "run_manifest.json"

# Read in chunks so a 6 GB genotype file does not have to be resident to be checksummed.
_CHUNK = 1 << 20


@dataclass
class SiteAllocation:
    """One site's share of a simulated cohort.

    Sizes live in a versioned scenario config rather than as literals in the splitting
    code, so that "five sites with this skew" and "two sites, balanced" are configuration
    changes rather than code edits. The original bundler had the five sizes inline, which
    made every alternative split a source modification.
    """

    site_id: str
    n_samples: int
    # Optional: pin a site to specific sample indices instead of taking the next slice.
    # Used for reproducing an existing split exactly.
    sample_indices: list[int] | None = None


@dataclass
class SimulationScenario:
    """A named, versioned description of one simulated dataset.

    Everything that affects the output belongs here, so that the manifest can record a
    single object and the run is reproducible from it.
    """

    name: str
    description: str = ""
    seed: int = 42
    sites: list[SiteAllocation] = field(default_factory=list)
    # Experiment-specific parameters, opaque to this layer.
    parameters: dict[str, Any] = field(default_factory=dict)

    @property
    def total_samples(self) -> int:
        return sum(site.n_samples for site in self.sites)

    @property
    def site_count(self) -> int:
        return len(self.sites)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> SimulationScenario:
        sites = [
            SiteAllocation(
                site_id=str(entry["site_id"]),
                n_samples=int(entry["n_samples"]),
                sample_indices=entry.get("sample_indices"),
            )
            for entry in data.get("sites", [])
        ]
        return cls(
            name=data["name"],
            description=data.get("description", ""),
            seed=int(data.get("seed", 42)),
            sites=sites,
            parameters=data.get("parameters", {}) or {},
        )

    @classmethod
    def load(cls, path: str | Path) -> SimulationScenario:
        import yaml

        path = Path(path)
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise ValueError(f"{path} must contain a YAML mapping")
        return cls.from_dict(data)

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description,
            "seed": self.seed,
            "sites": [asdict(s) for s in self.sites],
            "parameters": self.parameters,
        }


@dataclass
class RunManifest:
    """Everything needed to reproduce, or to prove you reproduced, one simulation run."""

    experiment: str
    scenario: dict[str, Any]
    seeds: dict[str, int]
    suite_version: str
    suite_commit: str | None
    python_version: str
    platform: str
    packages: dict[str, str]
    inputs: dict[str, str]
    outputs: dict[str, str]
    started_at: str
    finished_at: str
    notes: str = ""

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, sort_keys=True) + "\n"

    def write(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.to_json(), encoding="utf-8")
        return path

    @classmethod
    def load(cls, path: str | Path) -> RunManifest:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls(**data)

    def compare_outputs(self, other: RunManifest) -> dict[str, tuple[str, str]]:
        """Return {path: (expected, actual)} for every output that differs."""
        differences: dict[str, tuple[str, str]] = {}
        for name, digest in self.outputs.items():
            actual = other.outputs.get(name, "<missing>")
            if actual != digest:
                differences[name] = (digest, actual)
        for name in other.outputs:
            if name not in self.outputs:
                differences[name] = ("<not in reference>", other.outputs[name])
        return differences


def file_checksum(path: str | Path) -> str:
    """SHA-256 of a file, streamed."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def checksum_tree(root: str | Path, patterns: tuple[str, ...] = ("**/*",)) -> dict[str, str]:
    """Checksum every matching file under ``root``, keyed by relative POSIX path.

    Sorted, and relative, so the result is comparable across machines and directory
    layouts -- an absolute path would make two correct runs look different.

    Patterns are anchored at ``root``, which is the whole point of letting a caller pass
    them: they name the subset of the tree the run is responsible for. Matching with
    ``rglob`` instead silently prepends ``**/``, so an anchored pattern like ``anl/**/*``
    also matches ``processed/anl/**`` and ``ground_truth/phenotypes/anl/**``. That is not
    a cosmetic difference -- it pulled PLINK's ``.log`` files into the manifest, and those
    record the absolute output path and the hostname, so they cannot match on a rerun
    anywhere else. Use ``**/*`` to mean "the whole tree"; do not reach for ``rglob``.

    The run manifest itself is excluded. It is written into the output directory after
    the outputs are checksummed, so including it would make every run report its own
    manifest as an unexpected extra file on re-verification -- and no manifest can
    contain its own checksum anyway.
    """
    root = Path(root)
    seen: dict[str, str] = {}
    for pattern in patterns:
        for path in sorted(root.glob(pattern)):
            if path.is_file() and path.name != MANIFEST_FILENAME:
                seen[path.relative_to(root).as_posix()] = file_checksum(path)
    return seen


def _git_commit() -> str | None:
    """The suite's commit, if we are in a checkout and git is available.

    Returned with a '-dirty' suffix when the tree has uncommitted changes, because a
    manifest claiming a commit that does not describe the code that ran is worse than no
    commit at all.
    """
    try:
        root = Path(__file__).resolve().parent.parent.parent.parent
        commit = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if commit.returncode != 0:
            return None
        sha = commit.stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if dirty.returncode == 0 and dirty.stdout.strip():
            return f"{sha}-dirty"
        return sha
    except (OSError, subprocess.SubprocessError):
        return None


def _package_versions(names: tuple[str, ...]) -> dict[str, str]:
    from importlib.metadata import PackageNotFoundError, version

    out: dict[str, str] = {}
    for name in names:
        try:
            out[name] = version(name)
        except PackageNotFoundError:
            continue
    return out


# Packages whose version can change numerical output. Recorded on every run so that a
# result that stops reproducing can be traced to the thing that moved.
PROVENANCE_PACKAGES = (
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "pandas-plink",
    "xarray",
    "appfl",
    "appfl-bio-suite",
)


def record_provenance(
    experiment: str,
    scenario: SimulationScenario,
    seeds: dict[str, int],
    inputs: dict[str, str],
    output_root: str | Path,
    started_at: datetime,
    output_patterns: tuple[str, ...] = ("*",),
    notes: str = "",
) -> RunManifest:
    """Build the manifest for a completed simulation run."""
    from appfl_bio_suite import __version__

    return RunManifest(
        experiment=experiment,
        scenario=scenario.to_dict(),
        seeds=seeds,
        suite_version=__version__,
        suite_commit=_git_commit(),
        python_version=platform.python_version(),
        platform=f"{platform.system()} {platform.machine()}",
        packages=_package_versions(PROVENANCE_PACKAGES),
        inputs=inputs,
        outputs=checksum_tree(output_root, output_patterns),
        started_at=started_at.astimezone(UTC).isoformat(),
        finished_at=datetime.now(UTC).isoformat(),
        notes=notes,
    )


def verify_against_manifest(
    manifest_path: str | Path, output_root: str | Path, output_patterns: tuple[str, ...] = ("*",)
) -> tuple[bool, str]:
    """Re-checksum outputs and compare against a recorded manifest.

    Backs ``simulate <experiment> --verify <manifest>``. The point is to be able to
    answer "did this rerun produce the same data?" mechanically rather than by
    inspecting summary statistics and hoping.
    """
    reference = RunManifest.load(manifest_path)
    actual = checksum_tree(output_root, output_patterns)

    differences: list[str] = []
    for name, digest in sorted(reference.outputs.items()):
        got = actual.get(name)
        if got is None:
            differences.append(f"  MISSING  {name}")
        elif got != digest:
            differences.append(f"  DIFFERS  {name}")
            differences.append(f"           expected {digest[:16]}...")
            differences.append(f"           actual   {got[:16]}...")
    for name in sorted(actual):
        if name not in reference.outputs:
            differences.append(f"  EXTRA    {name}")

    if not differences:
        return True, (
            f"All {len(reference.outputs)} output(s) match the manifest.\n"
            f"Reference run: {reference.started_at} "
            f"(suite {reference.suite_version}, commit {reference.suite_commit})"
        )

    report = [
        f"{len(differences)} difference(s) against {manifest_path}:",
        *differences,
        "",
        "Reproducibility depends on the seeds, the input data, and the versions of the",
        "numerical packages. Compare the manifest's `packages` and `inputs` blocks",
        "against this environment -- a moved numpy or pandas-plink is the usual cause.",
    ]
    return False, "\n".join(report)


class SimulationPipeline(Protocol):
    """The shape every experiment's simulation implements.

    A Protocol rather than a base class: the stages have nothing in common
    implementation-wise, only in sequence, and inheritance would imply shared behaviour
    that does not exist.
    """

    experiment: str

    def simulate(self, scenario: SimulationScenario, out_dir: Path) -> dict[str, str]:
        """Generate the pooled dataset. Returns {label: checksum} for the inputs used."""
        ...

    def partition(self, scenario: SimulationScenario, out_dir: Path) -> None:
        """Split the pooled dataset into per-site datasets."""
        ...

    def materialize(self, scenario: SimulationScenario, out_dir: Path) -> None:
        """Write each site's dataset in the exact layout its loader expects."""
        ...


def run_pipeline(
    pipeline: SimulationPipeline,
    scenario: SimulationScenario,
    out_dir: str | Path,
    notes: str = "",
) -> RunManifest:
    """Run the three stages in order and write the manifest.

    The manifest is written last and always, so its presence means the run completed.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    started = datetime.now(UTC)

    inputs = pipeline.simulate(scenario, out_dir)
    pipeline.partition(scenario, out_dir)
    pipeline.materialize(scenario, out_dir)

    manifest = record_provenance(
        experiment=pipeline.experiment,
        scenario=scenario,
        seeds={"scenario": scenario.seed},
        inputs=inputs,
        output_root=out_dir,
        started_at=started,
        notes=notes,
    )
    manifest.write(out_dir / MANIFEST_FILENAME)
    return manifest
