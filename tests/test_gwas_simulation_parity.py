"""The GWAS simulation must reproduce its committed checksums, on any machine.

This is the acceptance criterion for the whole simulation migration. The code is
scientifically load-bearing -- it reproduces a published simulation design -- so a
refactor is only correct if it changes nothing observable.

Three levels of check, deliberately:

``test_port_matches_committed_golden``
    Compares every output against checksums committed to the repository. Runs anywhere,
    forever, including after the original source trees are retired. This is the one CI
    runs.

``test_golden_holds_on_an_older_cpu``
    Runs the same pipeline in a fresh interpreter made to compute like an older x86-64
    CPU, and holds it to the same checksums. A result that follows the hardware passes
    the first check on some machines and fails it on others. This check fails on every
    machine instead.

``test_port_matches_freshly_run_legacy``
    Re-runs the actual legacy scripts and compares. Stronger, because it cannot drift
    from reality, but only possible while those trees still exist. Skipped otherwise.

If the golden file ever needs regenerating, that is a decision to change published
numerical behaviour, not a routine fix. See tests/legacy_harness.py.
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from appfl_bio_suite.core.simulation import file_checksum

FIXTURES = Path(__file__).parent / "fixtures"
GOLDEN_PATH = FIXTURES / "gwas_legacy_golden.json"
PIPELINE_SCRIPT = FIXTURES / "golden_pipeline.py"


def _pipeline():
    """tests/fixtures/golden_pipeline.py: the run the golden file describes, defined once."""
    if str(FIXTURES) not in sys.path:
        sys.path.insert(0, str(FIXTURES))
    import golden_pipeline

    return golden_pipeline


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
    return _pipeline().run(tmp_path_factory.mktemp("ported"), golden)


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
    """Every output must be byte-identical to the committed checksum."""
    path = _pipeline().resolve(ported_run, key)
    assert path.is_file(), f"the port did not produce {key}"
    assert file_checksum(path) == golden["outputs"][key], (
        f"{key} differs from the committed golden checksum.\n"
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
        "  4. The polygenic score's summation in _compute_pgs was reordered, or handed "
        "to BLAS (@, np.dot, np.einsum). Only the PGS-derived files fail in that case -- "
        "pgs_scores, phenotypes_* and the summary -- while covariates and genotypes "
        "pass. Through BLAS, the failure comes and goes with the CI host, and "
        "test_golden_holds_on_an_older_cpu and test_pgs_is_the_reference_blas_product "
        "fail on every machine.\n"
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


@pytest.mark.skipif(
    platform.machine().lower() not in {"x86_64", "amd64"},
    reason="the emulated CPU is an x86-64 one; other architectures have other knobs",
)
def test_golden_holds_on_an_older_cpu(golden, tmp_path):
    """Every output must come out the same whatever CPU computes it.

    The committed checksums are only worth something if they hold on whatever machine
    runs the tests, and the CPU can leak into results without any code changing. OpenBLAS
    selects its kernels by CPU model, each rounding differently, and numpy dispatches its
    loops to the widest vector instructions the CPU has. A result that depends on either
    passes the golden check on some hosts and fails it on others.

    Not hypothetical. The polygenic score used to be a BLAS matrix-vector product, and
    the first golden checksums matched only what OpenBLAS's Haswell-family kernel
    produces -- the kernel on this project's AMD nodes. CI passed on hosts that select
    that kernel and failed the 13 PGS-derived files on the others, so re-running the same
    commit could go either way. Ten of those 13 failing files reproduce byte for byte
    with OPENBLAS_CORETYPE=Nehalem on an AMD node.

    So the pipeline runs again in a fresh interpreter made to compute like an older CPU:
    OpenBLAS's Nehalem kernel, and numpy's vectorized loops switched off down to its
    baseline. Both are read at load time, hence the separate process. Nehalem rather
    than the oldest kernel, Prescott, because Prescott happens to add in the reference
    order the checksums are defined by, so a score routed back through BLAS would still
    pass under it. Under Nehalem it fails on every machine, and so does anything that
    follows the vector width numpy dispatches to on this one.
    """
    pipeline = _pipeline()
    here = pipeline.cpu_state()
    persona = {"OPENBLAS_CORETYPE": "Nehalem"}
    if here["numpy_dispatch"]:
        persona["NPY_DISABLE_CPU_FEATURES"] = " ".join(here["numpy_dispatch"])

    result = subprocess.run(
        [sys.executable, str(PIPELINE_SCRIPT), str(tmp_path)],
        env={**os.environ, **persona},
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert result.returncode == 0, f"the pipeline failed under {persona}:\n{result.stderr}"
    report = json.loads(result.stdout.strip().splitlines()[-1])

    # A comparison between two identical configurations would prove nothing.
    if report["cpu"] == here:
        pytest.skip(
            f"{persona} changed nothing here (numpy and BLAS report {here} either way), "
            "so there is no second configuration to compare against"
        )

    differ = sorted(k for k, v in golden["outputs"].items() if report["outputs"].get(k) != v)
    assert not differ, (
        f"{len(differ)} output(s) changed when the same pipeline ran on a different CPU:\n"
        + "".join(f"  {key}\n" for key in differ)
        + f"\nthis machine: {here}\nemulated:     {report['cpu']}\n\n"
        "Something on that path now depends on the hardware. If only the PGS-derived "
        "files are listed, the score is going through BLAS again: _compute_pgs must sum "
        "its terms itself, in a fixed order. Do not regenerate the golden file to make "
        "this pass. A checksum captured on one machine would just fail on another."
    )


def _score_inputs(golden: dict):
    """Dosages and weights shaped like the real pipeline's, as _compute_pgs takes them.

    Includes the missingness that sends _compute_pgs through its mean-imputation path,
    and effect alleles on both strands so that the dosage flip runs too.
    """
    import xarray as xr

    fixture = golden["fixture"]
    n_samples, n_variants = fixture["n_samples"], fixture["n_variants"]

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
    return pgs_df, G, snp_to_pos, dosage


def _reference_pgs(pgs_df: pd.DataFrame, dosage: np.ndarray, chunk_size: int) -> np.ndarray:
    """What _compute_pgs must compute, in plain Python floats, one operation per term.

    Deliberately literal and slow, so the order of operations is not open to
    interpretation. Python floats are IEEE-754 doubles, rounded like numpy's. Every
    variant in ``pgs_df`` is genotyped here, in order, as _score_inputs builds them.
    """
    n_samples, n_variants = dosage.shape
    flip = list(pgs_df["effect_allele"].str.upper() == pgs_df["A1"].str.upper())
    weights = [float(w) for w in pgs_df["effect_weight"]]

    columns = []
    for j in range(n_variants):
        column = [float(d) for d in dosage[:, j]]
        if flip[j]:
            column = [2.0 - d for d in column]
        observed = [d for d in column if d == d]  # NaN is the one value unequal to itself
        mean = sum(observed) / len(observed)
        columns.append([d if d == d else mean for d in column])

    scores = [0.0] * n_samples
    for start in range(0, n_variants, chunk_size):
        # One chunk's product, as the reference BLAS defines DGEMV: from zero, a variant
        # at a time, one multiplication and one addition per term.
        product = [0.0] * n_samples
        for j in range(start, min(start + chunk_size, n_variants)):
            for i in range(n_samples):
                product[i] = product[i] + columns[j][i] * weights[j]
        for i in range(n_samples):
            scores[i] = scores[i] + product[i]
    return np.array(scores)


def test_pgs_is_the_reference_blas_product(golden):
    """_compute_pgs is chunk @ w the way the reference BLAS defines it, chunk by chunk.

    This is the definition the golden checksums rest on, checked against the same
    arithmetic written out term by term. It also covers what the golden run cannot: the
    fixture's scores fit in one chunk, while a chunk size of 7 makes each partial product
    restart from zero 43 times.
    """
    from appfl_bio_suite.experiments.gwas.simulation.phenotypes import _compute_pgs

    pgs_df, G, snp_to_pos, dosage = _score_inputs(golden)
    for chunk_size in (7, 2000):
        got = _compute_pgs(pgs_df, G, snp_to_pos, chunk_size)
        want = _reference_pgs(pgs_df, dosage, chunk_size)
        assert got.tobytes() == want.tobytes(), (
            f"at chunk size {chunk_size}, _compute_pgs no longer adds the terms in the "
            "reference order. If it now hands the product to BLAS, the scores follow the "
            "CPU again; if the order changed on purpose, that changes published numbers."
        )


def test_pgs_is_the_dot_product_to_rounding(golden):
    """The score is the dot product it claims to be, to within rounding.

    Checked against a plain matrix-vector product of the same flipped, imputed dosages,
    so a dropped, doubled or misaligned variant fails here, whatever the golden file and
    the reference-order check say.

    The bound is per sample, relative to the sum of that sample's terms' magnitudes. That
    is the scale the rounding error of any summation order is measured against, and
    near-zero scores from cancellation would make a relative tolerance meaningless.
    """
    from appfl_bio_suite.experiments.gwas.simulation.phenotypes import _compute_pgs
    from appfl_bio_suite.experiments.gwas.simulation.schema import PhenotypeParams

    pgs_df, G, snp_to_pos, dosage = _score_inputs(golden)
    score = _compute_pgs(pgs_df, G, snp_to_pos, PhenotypeParams().pgs_chunk_size)

    d = dosage.copy()
    flip = (pgs_df["effect_allele"] == pgs_df["A1"]).to_numpy()
    d[:, flip] = 2.0 - d[:, flip]
    d = np.where(np.isnan(d), np.nanmean(d, axis=0), d)
    w = pgs_df["effect_weight"].to_numpy()

    magnitude = np.abs(d * w).sum(axis=1)
    error = np.abs(score - d @ w)
    assert np.all(error <= 1e-12 * magnitude), (
        f"the score is not the dot product it should be: worst error "
        f"{np.max(error / magnitude):.3g} of the terms' magnitude, where rounding "
        "accounts for about 1e-14"
    )


@pytest.mark.slow
@pytest.mark.needs_legacy_tree
@pytest.mark.skipif(
    platform.machine().lower() not in {"x86_64", "amd64"},
    reason="the legacy scripts' BLAS kernel can only be pinned on x86-64",
)
def test_port_matches_freshly_run_legacy(tmp_path, golden):
    """Re-run the real legacy scripts and compare. Stronger, but needs the old trees.

    The harness pins OpenBLAS to its reference-order kernel, without which the legacy
    polygenic scores, and every file computed from them, would follow this machine's
    CPU. See "WHY THE BLAS KERNEL IS PINNED" in tests/legacy_harness.py.
    """
    sys.path.insert(0, str(Path(__file__).parent))
    sys.path.insert(0, str(FIXTURES))
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
