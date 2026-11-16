"""The coordinator-side statistics: pooling, and the QQ plot's quantiles.

Both bugs these cover produced a plausible-looking wrong answer rather than an error,
which is the category that most needs a test. Neither is caught by the loopback run: the
ci-tiny cohort happens to have no monomorphic variants, and it is far below the
downsampling threshold.
"""

from __future__ import annotations

import numpy as np
import pytest

from appfl_bio_suite.experiments.gwas.aggregator import MetaAnalysisAggregator
from appfl_bio_suite.experiments.gwas.plotting import plot_qq


@pytest.fixture
def aggregator(tmp_path):
    """A bare aggregator. Only `_meta_analyze` is exercised, so no APPFL run is needed."""
    return MetaAnalysisAggregator(aggregator_configs={"output_dir": str(tmp_path)})


def _variant_frame(n: int):
    import pandas as pd

    return pd.DataFrame(
        {
            "SNP": [f"rs{i}" for i in range(n)],
            "CHR": [1] * n,
            "BP": list(range(1, n + 1)),
        }
    )


def test_one_site_missing_a_variant_does_not_wipe_it_out_for_everyone(aggregator):
    """A variant monomorphic at ONE site must still be pooled from the others.

    That site's trainer returns beta=NaN, se=NaN. The aggregator gives it weight 0 --
    but `0.0 * NaN` is NaN, so summing the weighted betas used to propagate the NaN into
    the pooled estimate, and the variant silently vanished from the hits table and the
    plots. It is the exact failure the zero-weighting exists to prevent.
    """
    beta = np.array([[0.5, 0.2], [np.nan, 0.3]])
    se = np.array([[0.1, 0.1], [np.nan, 0.1]])

    out = aggregator._meta_analyze(beta, se, np.array([0.2, 0.3]), "BMI", 1000, _variant_frame(2))

    assert np.isfinite(out["BETA"]).all(), "a site with no estimate must not poison the pool"
    # Variant 0: only site 1 contributes, so the pooled estimate is that site's.
    assert out["BETA"].iloc[0] == pytest.approx(0.5)
    # Variant 1: both contribute with equal precision.
    assert out["BETA"].iloc[1] == pytest.approx(0.25)
    assert np.isfinite(out["P"]).all()


def test_a_variant_missing_everywhere_is_nan_not_a_fabricated_zero(aggregator):
    """No site has an estimate: the honest answer is NaN, not a pooled 0."""
    beta = np.array([[np.nan], [np.nan]])
    se = np.array([[np.nan], [np.nan]])

    out = aggregator._meta_analyze(beta, se, np.array([0.0]), "BMI", 1000, _variant_frame(1))

    assert np.isnan(out["BETA"].iloc[0])
    assert np.isnan(out["SE"].iloc[0])


def test_qq_expected_quantiles_use_the_full_test_count_when_downsampling(tmp_path):
    """Downsampling must not move the null line.

    Recomputing expected quantiles over the truncated array ranks the kept p-values as if
    they were the whole experiment, lifting the entire curve off the diagonal -- about
    0.6 log units on perfectly null data, which reads as genome-wide inflation that is
    not there.
    """
    rng = np.random.default_rng(0)
    n_tests, max_points = 200_000, 20_000
    p_values = rng.uniform(size=n_tests)

    captured = {}

    def capture(expected, observed, *args, **kwargs):
        captured["expected"] = np.asarray(expected)
        captured["observed"] = np.asarray(observed)

    import matplotlib.axes

    original = matplotlib.axes.Axes.scatter
    matplotlib.axes.Axes.scatter = lambda self, x, y, **kw: capture(x, y)
    try:
        plot_qq(p_values, "BMI", max_points, tmp_path / "qq.png")
    finally:
        matplotlib.axes.Axes.scatter = original

    expected, observed = captured["expected"], captured["observed"]
    assert expected.size == max_points, "the plot should carry only the retained points"

    # Uniform p-values are null by construction, so the cloud must sit on the diagonal.
    assert abs(float(np.mean(observed - expected))) < 0.05

    # The most significant p-value's expected quantile is set by the FULL test count.
    assert float(expected.max()) == pytest.approx(-np.log10(0.5 / n_tests))


# ---------------------------------------------------------------------------
# the variant axis is one axis
# ---------------------------------------------------------------------------


def _meta_tensor(frame):
    import json

    import torch

    blob = json.dumps({column: frame[column].tolist() for column in frame}).encode()
    return torch.frombuffer(bytearray(blob), dtype=torch.uint8)


@pytest.fixture
def pooling(aggregator):
    """The aggregator with a logger, as APPFL constructs it, for calls through `aggregate`."""
    import logging

    aggregator.logger = logging.getLogger("test-gwas-aggregator")
    return aggregator


def _axis(n=4):
    import pandas as pd

    return pd.DataFrame(
        {
            "CHR": [1] * n,
            "SNP": [f"rs{i}" for i in range(n)],
            "BP": [100 * (i + 1) for i in range(n)],
            "EA": ["A"] * n,
            "NEA": ["G"] * n,
        }
    )


def test_sites_coding_the_effect_allele_differently_are_not_pooled(pooling):
    """Equal lengths stack without complaint; the effects then cancel instead of adding."""
    swapped = _axis()
    swapped.loc[2, ["EA", "NEA"]] = ["G", "A"]
    models = {
        "Site1": {"variant_meta": _meta_tensor(_axis())},
        "Site2": {"variant_meta": _meta_tensor(swapped)},
    }
    with pytest.raises(ValueError) as raised:
        pooling.aggregate(models)
    message = str(raised.value)
    assert "Site2" in message and "row 2" in message
    assert "effect allele swapped" in message and "opposite signs" in message


def test_sites_with_the_same_variants_in_another_order_are_not_pooled(pooling):
    reordered = _axis().iloc[[0, 2, 1, 3]].reset_index(drop=True)
    models = {
        "Site1": {"variant_meta": _meta_tensor(_axis())},
        "Site2": {"variant_meta": _meta_tensor(reordered)},
    }
    with pytest.raises(ValueError, match="disagree at 2 of 4 variants"):
        pooling.aggregate(models)


def test_a_site_with_a_different_number_of_variants_is_named(pooling):
    models = {
        "Site1": {"variant_meta": _meta_tensor(_axis(4))},
        "Site2": {"variant_meta": _meta_tensor(_axis(3))},
    }
    with pytest.raises(ValueError, match="Site2 analyzed 3 variants but Site1 analyzed 4"):
        pooling.aggregate(models)
