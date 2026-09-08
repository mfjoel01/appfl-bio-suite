"""Unit tests for core.simulation's manifest: what it checksums, and what --verify says.

These exist because the pattern anchoring below was wrong for the whole life of the
feature and nothing noticed. ``checksum_tree`` matched with ``rglob``, which silently
prepends ``**/``, so the fine-mapping run manifest's ``<site>/**/*`` patterns also
swept up ``processed/<site>/**`` and ``ground_truth/phenotypes/<site>/**`` -- 30 extra
files in a 62-file manifest, including PLINK ``.log`` files that record the absolute
output path and the hostname and therefore cannot match on any rerun elsewhere.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from appfl_bio_suite.core.experiments import repo_root
from appfl_bio_suite.core.simulation import checksum_tree, verify_against_manifest


def _tree(root: Path) -> None:
    """A layout with the same shape as a fine-mapping run directory."""
    for rel in (
        "anl/data/site_genotypes.bed",           # the bundle -- in scope
        "anl/data/site_manifest.tsv",
        "processed/anl/anl_chr1.bed",            # simulation scratch -- out of scope
        "processed/anl/anl_chr1.log",            # ...and not reproducible anywhere else
        "ground_truth/causal_manifest.tsv",      # the answer key -- in scope
        "ground_truth/phenotypes/anl/L0.pheno",  # out of scope
        "loci/selected_loci.tsv",
    ):
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(rel)


_PATTERNS = ("anl/**/*", "ground_truth/causal_manifest.tsv", "loci/selected_loci.tsv")


def test_patterns_are_anchored_at_root(tmp_path: Path) -> None:
    """`anl/**/*` means the anl/ directory, not every directory named anl anywhere."""
    _tree(tmp_path)
    got = set(checksum_tree(tmp_path, _PATTERNS))
    assert got == {
        "anl/data/site_genotypes.bed",
        "anl/data/site_manifest.tsv",
        "ground_truth/causal_manifest.tsv",
        "loci/selected_loci.tsv",
    }
    assert not [name for name in got if name.startswith("processed/")]
    assert not [name for name in got if name.endswith(".log")]


def test_default_pattern_still_walks_the_whole_tree(tmp_path: Path) -> None:
    """The default has to stay recursive -- the GWAS manifest relies on it."""
    _tree(tmp_path)
    assert len(checksum_tree(tmp_path)) == 7


def test_verify_is_clean_when_nothing_in_scope_changed(tmp_path: Path) -> None:
    _tree(tmp_path)
    reference = checksum_tree(tmp_path, _PATTERNS)
    ok, report = verify_against_manifest(_manifest(tmp_path, reference), tmp_path, _PATTERNS)
    assert ok, report


def test_a_later_stage_writing_beside_the_scratch_is_not_drift(tmp_path: Path) -> None:
    """The regression: qc writes processed/<site>/<site>_kinship.* after the manifest.

    Those are a downstream stage's output, not the simulation's, and reporting them as
    unexpected extras made `--verify` cry wolf on a correct run.
    """
    _tree(tmp_path)
    reference = checksum_tree(tmp_path, _PATTERNS)
    manifest = _manifest(tmp_path, reference)
    (tmp_path / "processed/anl/anl_kinship.kin0").write_text("kinship")
    (tmp_path / "processed/anl/anl_kinship.log").write_text("kinship log")

    ok, report = verify_against_manifest(manifest, tmp_path, _PATTERNS)
    assert ok, report


def test_verify_still_detects_real_drift(tmp_path: Path) -> None:
    """Narrowing the scan must not make the check vacuous."""
    _tree(tmp_path)
    reference = checksum_tree(tmp_path, _PATTERNS)
    manifest = _manifest(tmp_path, reference)

    (tmp_path / "anl/data/site_manifest.tsv").write_text("tampered")
    (tmp_path / "anl/data/site_genotypes.bed").unlink()
    (tmp_path / "anl/data/UNEXPECTED.txt").write_text("new")

    ok, report = verify_against_manifest(manifest, tmp_path, _PATTERNS)
    assert not ok
    assert "DIFFERS  anl/data/site_manifest.tsv" in report
    assert "MISSING  anl/data/site_genotypes.bed" in report
    assert "EXTRA    anl/data/UNEXPECTED.txt" in report


def _manifest(root: Path, outputs: dict[str, str]) -> Path:
    """A minimal RunManifest on disk carrying `outputs`."""
    import json

    path = root / "reference_manifest.json"
    path.write_text(
        json.dumps(
            {
                "experiment": "fine-mapping",
                "scenario": {},
                "seeds": {},
                "suite_version": "0.0.0",
                "suite_commit": None,
                "python_version": "3.12.0",
                "platform": "test",
                "packages": {},
                "inputs": {},
                "outputs": outputs,
                "started_at": "2026-01-01T00:00:00+00:00",
                "finished_at": "2026-01-01T00:00:01+00:00",
                "notes": "",
            }
        )
    )
    return path



# ---------------------------------------------------------------------------
# the two paths that build a data package must build the same one
# ---------------------------------------------------------------------------
#
# `simulate` runs steps 0-5 in one process. The published scenario is far too large for
# that -- 79 loci x 15 architectures x 10 replicates over 150,000 individuals -- so it is
# built stage by stage, and `run_stage.py bundle` is the last of those stages.
#
# `bundle` used to stop after assembling the bundles: DRS registration and the DUO
# profiles lived only inside `simulate`. So the one scenario that ships produced bundles
# with no DRS registry and no consent terms, while the tiny scenario used for smoke tests
# produced both -- and the gap was invisible precisely because every test used the tiny
# one. These pin the two together.

needs_plink = pytest.mark.skipif(
    not (repo_root() / "vendor" / "bin" / "plink").exists(),
    reason="PLINK is not vendored; run scripts/fine-mapping/install_plink.sh",
)


@pytest.fixture(scope="module")
def two_packages(tmp_path_factory):
    """Build ci-tiny with `simulate`, then rebuild the last stage over the same inputs."""
    from appfl_bio_suite.experiments.fine_mapping.simulation import (
        load_scenario,
        run_simulation,
    )
    from appfl_bio_suite.experiments.fine_mapping.simulation.stages import bundle_main

    root = tmp_path_factory.mktemp("convergence")
    whole, staged = root / "whole", root / "staged"
    scenario = load_scenario("ci-tiny")
    run_simulation(scenario, whole)

    # The stage path is handed a run directory that already holds the earlier stages'
    # output, which is exactly what the sharded PBS job hands it.
    staged.mkdir()
    for name in ("raw", "processed", "loci", "ground_truth"):
        (staged / name).symlink_to(whole / name)
    bundle_main(["--scenario", "ci-tiny", "--out", str(staged)])
    return whole, staged


@needs_plink
def test_the_stage_path_registers_drs_and_duo_like_simulate(two_packages) -> None:
    """The regression: bundles that no coordinator could verify and no site could consent to."""
    import json

    whole, staged = two_packages
    for package in (whole, staged):
        assert (package / "drs_registry.json").is_file(), f"no DRS registry in {package.name}"
        for site in ("anl", "covenant", "mbzuai"):
            assert (package / site / "data" / "DATA_USE.json").is_file(), (
                f"no consent terms in {package.name}/{site}"
            )

    # Ids are content-addressed, so identical bundles must mint identical ids no matter
    # which path built them. This is the strongest statement of "the same package".
    def ids(package):
        registry = json.loads((package / "drs_registry.json").read_text())
        return {o["name"]: o["id"] for o in registry["objects"].values() if o.get("name")}

    assert ids(whole) == ids(staged)


@needs_plink
def test_both_paths_checksum_the_same_files(two_packages) -> None:
    """Including DATA_USE.json, which the stage path used to write after the manifest.

    Order is the whole point: the profile carries its bundle's drs_uri, so the DRS object
    must exist first -- and if the manifest is taken before the profiles are written,
    every bundle ends up holding a file the manifest does not list, which `--verify`
    reports as three unexpected extras on a package that is perfectly correct.
    """
    from appfl_bio_suite.core.simulation import RunManifest

    whole, staged = two_packages
    a = RunManifest.load(whole / "run_manifest.json").outputs
    b = RunManifest.load(staged / "run_manifest.json").outputs

    assert set(a) == set(b)
    assert [name for name in b if name.endswith("DATA_USE.json")], (
        "the manifest does not cover the consent profiles, so it was written before them"
    )
    assert {k: v for k, v in a.items()} == {k: v for k, v in b.items()}


@needs_plink
def test_validation_writes_beside_the_run_not_into_the_caller(two_packages, tmp_path) -> None:
    """`--out` used to default to the literal "reports/validation", a cwd-relative path.

    The Polaris job cds to the suite checkout before launching stages, so every
    production validation report was written into the source tree -- while the job's own
    completion banner pointed at ${FM_RUN_DIR}/reports/validation/SUMMARY.txt, which was
    never created. Running from anywhere must put the report beside the run.
    """
    import os
    import subprocess
    import sys

    whole, _ = two_packages
    elsewhere = tmp_path / "some-other-cwd"
    elsewhere.mkdir()

    env = {**os.environ, "PYTHONPATH": str(repo_root() / "src")}
    result = subprocess.run(
        [
            sys.executable,
            str(repo_root() / "scripts" / "fine-mapping" / "run_stage.py"),
            "validation",
            "--config",
            str(whole / "pipeline_config.yaml"),
        ],
        cwd=elsewhere,
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr[-2000:]
    assert (whole / "reports" / "validation" / "SUMMARY.txt").is_file(), (
        "validation did not write beside the run it validated"
    )
    assert not (elsewhere / "reports").exists(), (
        f"validation leaked output into the working directory: "
        f"{sorted(p.name for p in elsewhere.iterdir())}"
    )
