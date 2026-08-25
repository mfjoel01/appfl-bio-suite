"""Inverse-variance fixed-effect meta-analysis of per-site GWAS summary statistics.

COORDINATOR-SIDE. APPFL's ServerAgent loads ``aggregator_path`` locally and never ships
it, so unlike ``trainer.py`` and ``dataset.py`` this module may import freely from the
suite and depend on matplotlib.

WHAT IT DOES
------------
Each site returns per-variant effect sizes and standard errors from its own cohort. This
combines them into pooled estimates by weighting each site's estimate by its precision
(the inverse of its squared standard error), which is the standard fixed-effect
meta-analysis:

    w_i    = 1 / SE_i^2
    beta   = sum(w_i * beta_i) / sum(w_i)
    SE     = sqrt(1 / sum(w_i))

Larger and less noisy sites therefore count for more, automatically and without anyone
configuring a weight. Under the assumption that all sites estimate the same underlying
effect, this recovers what a pooled analysis of the combined cohort would have found --
which is the claim the experiment exists to test, and why the central baseline is worth
computing alongside it.

Single round. Sites do not train iteratively; they each compute once and report.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
from appfl.algorithm.aggregator import BaseAggregator  # noqa: E402
from scipy.stats import norm  # noqa: E402

from appfl_bio_suite.experiments.gwas.plotting import (  # noqa: E402
    plot_manhattan,
    plot_qq,
    write_hits_table,
)

DEFAULT_HIT_P_THRESHOLD = 5e-8


def _to_numpy(value):
    if torch.is_tensor(value):
        return value.detach().cpu().numpy()
    return np.asarray(value)


class MetaAnalysisAggregator(BaseAggregator):
    """Pools per-site summary statistics; writes tables and plots."""

    def __init__(self, model=None, aggregator_configs=None, logger=None):
        self.model = model
        self.logger = logger
        self.aggregator_configs = aggregator_configs or {}
        self.hit_threshold = float(
            self.aggregator_configs.get("hit_p_threshold", DEFAULT_HIT_P_THRESHOLD)
        )
        self.qq_max_points = int(self.aggregator_configs.get("qq_max_points", 250000))
        self.plot_per_site = bool(self.aggregator_configs.get("plot_per_site", True))

        self.output_dir = Path(
            self.aggregator_configs.get("output_dir", "local/output/gwas")
        ).resolve()
        self.data_dir = self.output_dir / "data"
        self.graphs_dir = self.output_dir / "graphs"
        self.logs_dir = self.output_dir / "logs"
        for directory in (self.data_dir, self.graphs_dir, self.logs_dir):
            directory.mkdir(parents=True, exist_ok=True)

        self.global_state = {"meta_ready": torch.tensor([0], dtype=torch.int64)}

    def get_parameters(self, **kwargs):
        return self.global_state

    def aggregate(self, local_models, **kwargs):
        client_ids = list(local_models.keys())
        self.logger.info(f"meta-analysis across {len(client_ids)} site(s)")

        # The variant axis comes from the first site's metadata. Every site must be on
        # the same variant set for meta-analysis to mean anything; the shape check below
        # is what enforces it.
        meta_bytes = _to_numpy(local_models[client_ids[0]]["variant_meta"]).tobytes()
        variant_df = pd.DataFrame(json.loads(meta_bytes.decode("utf-8")))

        def stack(key):
            return np.vstack([_to_numpy(local_models[cid][key]) for cid in client_ids])

        def scalars(key, dtype=np.float64):
            return np.array(
                [float(_to_numpy(local_models[cid][key])[0]) for cid in client_ids], dtype=dtype
            )

        bmi_beta_stack = stack("bmi_beta")
        bmi_se_stack = stack("bmi_se")
        t2d_beta_stack = stack("t2d_beta")
        t2d_se_stack = stack("t2d_se")
        maf_stack = stack("maf")
        gwas_n = scalars("gwas_n")
        eval_n = scalars("eval_n")
        bmi_r2 = scalars("local_bmi_r2")
        t2d_auc = scalars("local_t2d_auc")

        n_variants = len(variant_df)
        if bmi_beta_stack.shape[1] != n_variants:
            raise ValueError(
                f"site payloads carry {bmi_beta_stack.shape[1]} variants but the variant "
                f"metadata from {client_ids[0]} describes {n_variants}.\n"
                "Every site must analyze the SAME variant set. The usual cause is sites "
                "running with different variant_scaling values, or one site having been "
                "sent a bundle built from a different simulation run."
            )

        total_n = int(gwas_n.sum())
        # Weighted by cohort size: a site's allele frequency should count in proportion
        # to how many people it observed.
        meta_maf = np.average(maf_stack, axis=0, weights=gwas_n)

        bmi_df = self._meta_analyze(
            bmi_beta_stack, bmi_se_stack, meta_maf, "BMI", total_n, variant_df
        )
        t2d_df = self._meta_analyze(
            t2d_beta_stack, t2d_se_stack, meta_maf, "T2D", total_n, variant_df
        )

        bmi_df.to_csv(self.data_dir / "appfl_meta_gwas_bmi.csv.gz", index=False)
        t2d_df.to_csv(self.data_dir / "appfl_meta_gwas_t2d.csv.gz", index=False)
        write_hits_table(
            bmi_df, t2d_df, self.hit_threshold, self.data_dir / "appfl_meta_gwas_hits.csv"
        )
        self._write_site_metrics(client_ids, gwas_n, eval_n, bmi_r2, t2d_auc)

        for trait, frame in (("bmi", bmi_df), ("t2d", t2d_df)):
            plot_manhattan(
                frame, trait.upper(), self.hit_threshold,
                self.graphs_dir / f"appfl_meta_gwas_{trait}_manhattan.png", label="meta",
            )
            plot_qq(
                frame["P"].to_numpy(dtype=np.float64), trait.upper(), self.qq_max_points,
                self.graphs_dir / f"appfl_meta_gwas_{trait}_qq.png",
            )

        if self.plot_per_site:
            self._plot_per_site(
                client_ids, variant_df, bmi_beta_stack, bmi_se_stack,
                t2d_beta_stack, t2d_se_stack, maf_stack, gwas_n,
            )

        self.global_state = {
            "meta_ready": torch.tensor([1], dtype=torch.int64),
            "num_clients": torch.tensor([len(client_ids)], dtype=torch.int64),
            "num_variants": torch.tensor([n_variants], dtype=torch.int64),
            "bmi_hits": torch.tensor(
                [int((bmi_df["P"] < self.hit_threshold).sum())], dtype=torch.int64
            ),
            "t2d_hits": torch.tensor(
                [int((t2d_df["P"] < self.hit_threshold).sum())], dtype=torch.int64
            ),
            "weighted_local_bmi_r2": torch.tensor(
                [float(np.average(bmi_r2, weights=eval_n))], dtype=torch.float64
            ),
            "weighted_local_t2d_auc": torch.tensor(
                [float(np.average(t2d_auc, weights=eval_n))], dtype=torch.float64
            ),
        }
        self.logger.info(f"meta-analysis written to {self.output_dir}/{{data,graphs}}")
        return self.global_state

    def _meta_analyze(self, beta_stack, se_stack, maf, trait, total_n, variant_df):
        """Inverse-variance fixed-effect pooling, per variant."""
        with np.errstate(divide="ignore", invalid="ignore"):
            # A site contributes zero weight where its estimate is missing or its
            # standard error is non-positive -- which happens for monomorphic variants
            # at that site. Zero weight rather than NaN so one bad site does not wipe
            # out a variant for everyone.
            usable = np.isfinite(se_stack) & (se_stack > 0) & np.isfinite(beta_stack)
            weights = np.where(usable, 1.0 / np.square(se_stack), 0.0)
        weight_sum = weights.sum(axis=0)

        # Zeroing the weight is not enough on its own. A monomorphic variant comes back
        # as beta=NaN, and `0.0 * NaN` is NaN rather than 0, so summing the products
        # would propagate that site's NaN into the pooled estimate for every OTHER site
        # -- the variant then vanishes from the hits table and the plots, which is
        # precisely what the zero weight above exists to prevent. Select the discarded
        # terms out instead of relying on multiplication by zero to erase them.
        contributions = np.where(usable, weights * beta_stack, 0.0)

        beta = np.divide(
            contributions.sum(axis=0), weight_sum,
            out=np.full(weight_sum.shape, np.nan, dtype=np.float64), where=weight_sum > 0,
        )
        se = np.sqrt(
            np.divide(
                1.0, weight_sum,
                out=np.full(weight_sum.shape, np.nan, dtype=np.float64), where=weight_sum > 0,
            )
        )
        stat = np.divide(
            beta, se,
            out=np.zeros(weight_sum.shape, dtype=np.float64),
            where=np.isfinite(se) & (se > 0),
        )
        p_value = np.clip(2.0 * norm.sf(np.abs(stat)), np.finfo(np.float64).tiny, 1.0)

        out = variant_df.copy()
        out["TRAIT"] = trait
        out["BETA"] = beta
        out["SE"] = se
        out["STAT"] = stat
        if trait == "T2D":
            out["OR"] = np.exp(np.clip(beta, -50, 50))
        out["P"] = p_value
        out["MAF"] = maf
        out["N_META"] = total_n
        return out

    def _write_site_metrics(self, client_ids, gwas_n, eval_n, bmi_r2, t2d_auc):
        pd.DataFrame(
            {
                "CLIENT_ID": client_ids,
                "GWAS_N": gwas_n.astype(int),
                "EVAL_N": eval_n.astype(int),
                "LOCAL_BMI_R2": bmi_r2,
                "LOCAL_T2D_AUROC": t2d_auc,
            }
        ).to_csv(self.data_dir / "appfl_site_pgs_metrics.csv", index=False)

        pd.DataFrame(
            [
                {
                    "NUM_CLIENTS": len(client_ids),
                    "TOTAL_GWAS_N": int(gwas_n.sum()),
                    "TOTAL_EVAL_N": int(eval_n.sum()),
                    "WEIGHTED_LOCAL_BMI_R2": float(np.average(bmi_r2, weights=eval_n)),
                    "WEIGHTED_LOCAL_T2D_AUROC": float(np.average(t2d_auc, weights=eval_n)),
                }
            ]
        ).to_csv(self.data_dir / "appfl_meta_summary.csv", index=False)

    def _plot_per_site(
        self, client_ids, variant_df, bmi_beta, bmi_se, t2d_beta, t2d_se, maf, gwas_n
    ):
        """Render each site's own Manhattan and QQ plots, here rather than at the site.

        Everything needed is in the payload. Producing them here means matplotlib is not
        a partner dependency, and it puts every site's diagnostics somewhere the
        coordinator can actually look at them -- previously they were written to the
        partner's filesystem and never seen.
        """
        per_site_dir = self.graphs_dir / "per_site"
        per_site_dir.mkdir(parents=True, exist_ok=True)

        for i, client_id in enumerate(client_ids):
            for trait, beta_row, se_row in (
                ("BMI", bmi_beta[i], bmi_se[i]),
                ("T2D", t2d_beta[i], t2d_se[i]),
            ):
                with np.errstate(divide="ignore", invalid="ignore"):
                    stat = np.divide(
                        beta_row, se_row,
                        out=np.zeros_like(beta_row),
                        where=np.isfinite(se_row) & (se_row > 0),
                    )
                p_value = np.clip(2.0 * norm.sf(np.abs(stat)), np.finfo(np.float64).tiny, 1.0)

                frame = variant_df.copy()
                frame["BETA"] = beta_row
                frame["SE"] = se_row
                frame["P"] = p_value
                frame["MAF"] = maf[i]
                frame["N"] = int(gwas_n[i])

                lower = trait.lower()
                plot_manhattan(
                    frame, trait, self.hit_threshold,
                    per_site_dir / f"{client_id}_{lower}_manhattan.png", label=str(client_id),
                )
                plot_qq(
                    p_value, trait, self.qq_max_points,
                    per_site_dir / f"{client_id}_{lower}_qq.png",
                )
