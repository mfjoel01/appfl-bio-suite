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

import numpy as np
import pandas as pd
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
        "  4. The BLAS thread pin around the accumulation in _compute_pgs was removed. "
        "Only the PGS-derived files fail in that case -- pgs_scores, phenotypes_* and "
        "the summaries -- while covariates and genotypes pass, and it reproduces only "
        "on machines with enough cores. test_pgs_accumulation_does_not_follow_the_"
        "machine is the direct check.\n"
        "  5. The groupby/join/filter order in _prep_pgs changed.\n"
        "\n"
        "If none of those changed, compare this environment against constraints.txt "
        "before touching the golden file -- see docs/coordinator/releasing.md."
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


def test_pgs_accumulation_does_not_follow_the_machine(golden):
    """The polygenic score must not depend on how many cores the machine has.

    ``_compute_pgs`` sums each chunk with a BLAS matrix-vector product. Given eight or
    more threads OpenBLAS splits that reduction and adds the partial sums in a different
    order, which moves the low-order bits of every score and, through them, every
    phenotype and summary file downstream. Nothing raises; the outputs are simply not
    the published ones.

    Not hypothetical: this is what made the committed golden file stop matching on a
    larger runner while every line of the pipeline was unchanged. phenotypes.py pins the
    accumulation to one BLAS thread, and this holds it to that.

    Eight is the threshold measured for a chunk this wide, so the comparison is made
    there rather than at the machine's core count -- that way the test still exercises
    the split on a small runner, where asking for the core count would ask for two.
    """
    import xarray as xr
    from threadpoolctl import threadpool_limits

    from appfl_bio_suite.experiments.gwas.simulation.phenotypes import _compute_pgs
    from appfl_bio_suite.experiments.gwas.simulation.schema import PhenotypeParams

    fixture = golden["fixture"]
    n_samples, n_variants = fixture["n_samples"], fixture["n_variants"]

    # Dosages and weights shaped like the real pipeline's, including the missingness
    # that sends _compute_pgs through its mean-imputation path. The chunk size is the
    # production default, because chunk width is what decides whether BLAS splits at all.
    rng = np.random.default_rng(fixture["seed"])
    dosage = rng.integers(0, 3, size=(n_samples, n_variants)).astype(np.float64)
    dosage[rng.random(dosage.shape) < 0.01] = np.nan
    G = xr.DataArray(dosage, dims=("sample", "variant"))

    snps = [f"rs{i}" for i in range(n_variants)]
    pgs_df = pd.DataFrame(
        {
            "SNP": snps,
            "effect_allele": np.where(rng.random(n_variants) < 0.5, "A", "G"),
            "A1": "A",
            "effect_weight": rng.normal(size=n_variants),
        }
    )
    snp_to_pos = {snp: i for i, snp in enumerate(snps)}
    chunk_size = PhenotypeParams().pgs_chunk_size

    with threadpool_limits(limits=1, user_api="blas"):
        serial = _compute_pgs(pgs_df, G, snp_to_pos, chunk_size)
    with threadpool_limits(limits=8, user_api="blas"):
        parallel = _compute_pgs(pgs_df, G, snp_to_pos, chunk_size)

    assert np.array_equal(serial, parallel), (
        "the polygenic score changed when BLAS was allowed 8 threads instead of 1. "
        "phenotypes.py is supposed to pin that accumulation to a single thread; if that "
        "pin was removed, these numbers now depend on the machine and the committed "
        "golden file will fail wherever the core count differs."
    )
