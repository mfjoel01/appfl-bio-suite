"""The GWAS simulation must be output-identical to the pipeline it replaces.

This is the acceptance criterion for the whole simulation migration. The code is
scientifically load-bearing -- it reproduces a published simulation design -- so a
refactor is only correct if it changes nothing observable.

Two levels of check, deliberately:

``test_port_matches_committed_golden``
    Compares against checksums captured from the pre-migration scripts and committed to
    the repository. Runs anywhere, forever, including after the original source trees are
    retired. This is the one CI runs.

``test_port_matches_freshly_run_legacy``
    Re-runs the actual legacy scripts and compares. Stronger, because it cannot drift
    from reality, but only possible while those trees still exist. Skipped otherwise.

If the golden file ever needs regenerating, that is a decision to change published
numerical behaviour, not a routine fix. See tests/legacy_harness.py.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from appfl_bio_suite.core.simulation import file_checksum

GOLDEN_PATH = Path(__file__).parent / "fixtures" / "gwas_legacy_golden.json"


@pytest.fixture(scope="module")
def golden() -> dict:
    assert GOLDEN_PATH.is_file(), (
        f"{GOLDEN_PATH} is missing. It is the reference the simulation port is held "
        "against and must be committed."
    )
    return json.loads(GOLDEN_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def ported_run(tmp_path_factory, golden) -> Path:
    """Run the ported pipeline against the same fixture the golden file used."""
    import sys

    sys.path.insert(0, str(Path(__file__).parent / "fixtures"))
    from make_fixture import build

    from appfl_bio_suite.experiments.gwas.simulation.bundler import bundle_sites
    from appfl_bio_suite.experiments.gwas.simulation.phenotypes import simulate_phenotypes

    root = tmp_path_factory.mktemp("ported")
    fixture_cfg = golden["fixture"]
    site_sizes = golden["site_sizes"]

    info = build(
        root,
        n_samples=fixture_cfg["n_samples"],
        n_variants=fixture_cfg["n_variants"],
        seed=fixture_cfg["seed"],
    )
    simulate_phenotypes(
        plink_prefix=info["plink_prefix"],
        pgs_t2d_path=info["pgs_t2d"],
        pgs_bmi_path=info["pgs_bmi"],
        out_dir=root / "out",
    )
    bundle_sites(
        plink_prefix=info["plink_prefix"],
        pooled_dir=root / "out",
        sites_dir=root / "sites",
        site_sizes=site_sizes,
        split_seed=42,
    )
    return root


def _resolve(root: Path, key: str) -> Path:
    if key.startswith("pooled/"):
        return root / "out" / key.split("/", 1)[1]
    site, name = key.split("/", 1)
    return root / "sites" / site / "data" / name


def test_golden_file_is_well_formed(golden):
    assert golden["outputs"], "golden file records no outputs"
    assert golden["site_sizes"], "golden file records no site split"
    # Full variant set, not the smoke fraction -- a golden captured at 0.02 scaling would
    # leave most of the pipeline uncharacterized.
    assert golden["data_sim_scaling"] == 1.0


@pytest.mark.parametrize(
    "key", sorted(json.loads(GOLDEN_PATH.read_text())["outputs"]) if GOLDEN_PATH.is_file() else []
)
def test_port_matches_committed_golden(ported_run: Path, golden: dict, key: str):
    """Every output must be byte-identical to the pre-migration pipeline."""
    path = _resolve(ported_run, key)
    assert path.is_file(), f"the port did not produce {key}"
    assert file_checksum(path) == golden["outputs"][key], (
        f"{key} differs from the pre-migration pipeline.\n"
        "\n"
        "This code backs published figures, so any difference is a defect unless you "
        "are deliberately changing the science.\n"
        "\n"
        "Most likely causes, in order:\n"
        "  1. An RNG call was reordered, or numpy's legacy global RNG "
        "(np.random.seed + np.random.shuffle in bundler.py) was replaced with "
        "default_rng -- a different algorithm producing a different split.\n"
        "  2. ndarray.std() (ddof=0) replaced with pandas .std() (ddof=1).\n"
        "  3. The PGS chunk size changed, altering float accumulation order.\n"
        "  4. The groupby/join/filter order in _prep_pgs changed."
    )


def test_no_unexpected_outputs(ported_run: Path, golden: dict):
    """The port must not quietly add or drop files a downstream consumer might rely on."""
    produced = {f"pooled/{p.name}" for p in (ported_run / "out").glob("*") if p.is_file()}
    for site_dir in sorted((ported_run / "sites").glob("Site*")):
        produced |= {
            f"{site_dir.name}/{p.name}" for p in (site_dir / "data").glob("*") if p.is_file()
        }
    expected = set(golden["outputs"])
    assert produced == expected, (
        f"output set changed.\n  only in port:   {sorted(produced - expected)}\n"
        f"  only in golden: {sorted(expected - produced)}"
    )


@pytest.mark.slow
@pytest.mark.needs_legacy_tree
def test_port_matches_freshly_run_legacy(tmp_path, golden):
    """Re-run the real legacy scripts and compare. Stronger, but needs the old trees."""
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    sys.path.insert(0, str(Path(__file__).parent / "fixtures"))
    from legacy_harness import collect_outputs, legacy_tree_available, run_legacy_pipeline

    if not legacy_tree_available():
        pytest.skip(
            "legacy GWAS tree not configured: set APPFL_BIO_SUITE_LEGACY_GWAS to an "
            "APPFL_GWAS/APPFL checkout (skips whether the trees exist or not)"
        )

    site_sizes = golden["site_sizes"]
    local = run_legacy_pipeline(tmp_path, site_sizes)
    fresh = collect_outputs(local, n_sites=len(site_sizes))

    assert fresh == golden["outputs"], (
        "the committed golden file no longer matches a fresh run of the legacy scripts. "
        "Either the legacy tree changed, or the environment's numerical packages moved."
    )
