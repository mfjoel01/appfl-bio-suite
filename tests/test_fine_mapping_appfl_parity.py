"""The APPFL path must reproduce the standalone driver's results, column for column.

WHAT THIS PROTECTS
------------------
``tests/test_fine_mapping_federated.py`` already proves the hard thing: that the federated
statistics equal the centralized ones (algo.md Corollary 1), running the vendored code
directly. That test knows nothing about APPFL.

This one closes the remaining gap. Between the vendored code and a real federated run sit
two pieces of this repository's own machinery -- the shipped trainer that computes a
site's aggregates, and the aggregator that decodes a wire payload back into them. Either
could lose or corrupt a number without anything else noticing: a mis-keyed tensor, a
silent dtype narrowing, a variant list that arrives in a different order. The credible
sets would still look like credible sets.

So this runs both paths over one simulated package and asserts the result tables match on
every statistic column, at zero tolerance. ``runtime_s`` is wall clock and is excluded;
everything else -- credible-set sizes, capture, PIPs, purities, per-ancestry minimum
p-values -- must be identical.

If this fails while ``test_fine_mapping_shipped_parity.py`` passes, the aggregate
computation is fine and the fault is in transport or reassembly: look at the payload keys
and at ``FineMappingAggregator._decode``.

COST
----
Marked slow, and skipped without PLINK and SuSiEx on PATH or in ``vendor/bin``. It
simulates a real package and runs SuSiEx twice over every instance, which is tens of
seconds -- worth it, because this is the only test that exercises the actual production
path end to end.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import pandas as pd
import pytest

from appfl_bio_suite.experiments.fine_mapping.dataset import get_dataset
from appfl_bio_suite.experiments.fine_mapping.simulation import (
    load_scenario,
    run_simulation,
)
from appfl_bio_suite.experiments.fine_mapping.trainer import SiteFineMappingTrainer

# Statistics only. runtime_s is wall clock; `error` is empty string on one path and NaN on
# the other, and carries no result.
_EXCLUDED_COLUMNS = {"runtime_s", "error"}

# Must match what the aggregator is constructed with below, or the two paths would be
# fine-mapping with different settings and the comparison would mean nothing.
LEVEL, PVAL_THRESH, MAF = 0.95, 1e-5, 0.001

# What the shipped ci-tiny scenario's three profiles require between them: biomedical
# research (so `covenant`'s HMB permission is satisfied), by a non-commercial
# not-for-profit that agrees to publish. Same request `run --driver serial` synthesizes.
DATA_USE_REQUEST = {
    "requester": "parity-test@localhost",
    "purposes": ["DUO:0000038"],
    "non_commercial": True,
    "not_for_profit_organisation": True,
    "publication_agreed": True,
}

pytestmark = pytest.mark.slow


def _binary(name: str) -> str | None:
    found = shutil.which(name)
    if found:
        return found
    vendored = Path("vendor") / "bin" / name
    return str(vendored.resolve()) if vendored.is_file() else None


@pytest.fixture(scope="module")
def binaries():
    resolved = {name: _binary(name) for name in ("SuSiEx", "plink")}
    missing = [name for name, path in resolved.items() if path is None]
    if missing:
        pytest.skip(
            f"fine-mapping binaries unavailable: {missing}. "
            "Run scripts/fine-mapping/install_susiex.sh and install_plink.sh."
        )
    return resolved


@pytest.fixture(scope="module")
def package(tmp_path_factory: pytest.TempPathFactory, binaries):
    """Simulate the ci-tiny scenario once, and share it between both paths."""
    out = tmp_path_factory.mktemp("fm-appfl-parity") / "data"
    scenario = load_scenario("ci-tiny")
    run_simulation(scenario, out)
    return out, scenario


@pytest.fixture(scope="module")
def standalone(package, binaries):
    """Run the vendored driver, exactly as the standalone repository does."""
    from appfl_bio_suite.experiments.fine_mapping.fedfm.fed_fine_mapping import (
        run_fed_fine_mapping,
    )

    out, scenario = package
    cfg = scenario.to_pipeline_config(out)
    run_fed_fine_mapping(cfg, n_workers=1, level=LEVEL, pval_thresh=PVAL_THRESH, maf=MAF)
    results = cfg.resolved_path("reports_dir") / "fed_fine_mapping" / "fed_fm_results.tsv"
    assert results.is_file(), "the standalone driver wrote no results"
    return pd.read_csv(results, sep="\t")


@pytest.fixture(scope="module")
def through_appfl(package, binaries, tmp_path_factory):
    """Run the shipped trainer at each site, then the aggregator, as APPFL would.

    Driven directly rather than through ``run --driver serial``, so a failure points at
    this repository's trainer or aggregator rather than at APPFL's scheduler. The serial
    driver does no more than this loop.
    """
    from appfl_bio_suite.experiments.fine_mapping.aggregator import FineMappingAggregator

    out, scenario = package
    output_dir = tmp_path_factory.mktemp("fm-appfl-out")
    logger = logging.getLogger("appfl-parity")

    # The simulated bundles now carry DUO terms, so the loader needs the study that is
    # asking -- exactly as a real run supplies it. Passing one here is not test
    # scaffolding: a bundle with declared terms REFUSES a run that declares nothing, and
    # a parity test that bypassed that would be testing a path no run takes.
    local_models = {}
    for site in scenario.site_ids:
        dataset, _ = get_dataset(out / site / "data", site, data_use_request=DATA_USE_REQUEST)
        trainer = SiteFineMappingTrainer(
            train_dataset=dataset,
            train_configs={"trainer_output_dirname": str(output_dir / site)},
            logger=logger,
            client_id=site,
        )
        trainer.train()
        local_models[site] = trainer.get_parameters()

    aggregator = FineMappingAggregator(
        aggregator_configs={
            "level": LEVEL,
            "pval_thresh": PVAL_THRESH,
            "maf": MAF,
            "n_workers": 1,
            "make_figures": False,
            "output_dir": str(output_dir),
            "causal_manifest": str(out / "ground_truth" / "causal_manifest.tsv"),
            "susiex_binary": binaries["SuSiEx"],
            "plink_binary": binaries["plink"],
        },
        logger=logger,
    )
    aggregator.aggregate(local_models)
    return pd.read_csv(output_dir / "data" / "fed_fm_results.tsv", sep="\t")


def _aligned(a: pd.DataFrame, b: pd.DataFrame):
    key = ["locus_id", "architecture_id", "replicate"]
    return (
        a.sort_values(key).reset_index(drop=True),
        b.sort_values(key).reset_index(drop=True),
    )


def test_both_paths_fine_map_the_same_instances(standalone, through_appfl):
    appfl, stand = _aligned(through_appfl, standalone)
    assert len(appfl) == len(stand) > 0
    assert list(appfl["locus_id"]) == list(stand["locus_id"])
    assert list(appfl["architecture_id"]) == list(stand["architecture_id"])


def test_both_paths_report_the_same_columns(standalone, through_appfl):
    """A column present on one side only is a silently narrower result table."""
    only_appfl = set(through_appfl.columns) - set(standalone.columns) - _EXCLUDED_COLUMNS
    only_stand = set(standalone.columns) - set(through_appfl.columns) - _EXCLUDED_COLUMNS
    assert not only_appfl, f"columns only in the APPFL results: {sorted(only_appfl)}"
    assert not only_stand, f"columns only in the standalone results: {sorted(only_stand)}"


def test_every_statistic_is_identical(standalone, through_appfl):
    """Zero tolerance. The two paths run the same coordinator code on the same numbers,
    so anything other than equality means a number was lost in transport."""
    import numpy as np

    appfl, stand = _aligned(through_appfl, standalone)
    compared = [c for c in appfl.columns if c in stand.columns and c not in _EXCLUDED_COLUMNS]
    assert len(compared) > 15, "suspiciously few columns compared; did the schema change?"

    mismatched = []
    for column in compared:
        left, right = appfl[column], stand[column]
        if pd.api.types.is_numeric_dtype(left) and pd.api.types.is_numeric_dtype(right):
            same = np.array_equal(
                left.to_numpy(dtype=float), right.to_numpy(dtype=float), equal_nan=True
            )
        else:
            same = left.astype(str).equals(right.astype(str))
        if not same:
            mismatched.append(f"{column}: appfl={left.tolist()} standalone={right.tolist()}")

    assert not mismatched, "the APPFL path and the standalone driver disagree:\n" + "\n".join(
        mismatched
    )


def test_the_run_actually_produced_credible_sets(through_appfl):
    """Guard the fixture. Two paths that both find nothing agree trivially.

    ci-tiny uses a deliberately enormous effect size precisely so that credible sets come
    back; if that stops being true the parity assertions above go vacuous while still
    passing.
    """
    assert through_appfl["n_credible_sets"].sum() > 0, (
        "no credible sets were produced, so the parity comparison proves nothing. "
        "Check the ci-tiny scenario's h2 and the aggregator's pval_thresh."
    )
    assert through_appfl["n_causal"].max() > 0, (
        "no ground truth was loaded, so the scoring columns are all empty and excluded "
        "from the comparison by being identical NaN. Check causal_manifest."
    )


def test_the_aggregator_writes_its_site_summary(through_appfl, package, tmp_path_factory):
    """The per-site record of what each partner sent. It is the only place the uplink
    size and the harmonization count are retained after a run."""
    out, scenario = package
    # Written beside the results by the same aggregate() call the fixture made.
    summaries = list(Path(tmp_path_factory.getbasetemp()).rglob("fed_fm_site_summary.csv"))
    assert summaries, "no fed_fm_site_summary.csv was written"
    summary = pd.read_csv(summaries[0])
    assert set(summary["CLIENT_ID"]) == set(scenario.site_ids)
    assert (summary["N_GENO_BLOCKS"] > 0).all()
