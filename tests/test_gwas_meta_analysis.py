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
