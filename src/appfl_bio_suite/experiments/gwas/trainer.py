"""Per-site GWAS and polygenic-score evaluation. SHIPPED TO WORKERS.

=============================================================================
THIS FILE'S SOURCE IS SENT OVER THE WIRE AND EXECUTED ON A PARTNER'S CLUSTER.
=============================================================================

Self-contained by design: standard library, numpy/pandas/scipy/sklearn, pandas-plink,
torch, and ``appfl`` only. No imports from ``appfl_bio_suite``, no sibling imports.
Enforced by tests/test_shipped_modules.py.

The original version imported two helper modules by bare name, which is why every partner
had to add a ``PYTHONPATH`` line to their endpoint's ``worker_init`` -- and a missing one
was the largest single source of partner-side breakage on this project. Making this file
stand alone deleted that line from the setup guide.

WHAT EACH SITE COMPUTES, AND WHAT LEAVES
----------------------------------------
Runs a complete local GWAS over its own cohort, then returns **only** per-variant summary
statistics: effect size, standard error, minor-allele frequency, sample counts, and two
scalar polygenic-score metrics. Genotypes and phenotypes never leave the site, and the
payload is a few megabytes regardless of cohort size.

Single round. ``_has_run`` guards it: this is summary-statistic federated learning, not
iterative training, so there is exactly one exchange.

  BMI (continuous)  -- residualize phenotype and genotypes on covariates, then per
                       variant beta = G'y / ||G||^2 with a two-sided t test.
  T2D (binary)      -- fit a covariate-only null, then a logistic score test per variant.
                       A score test rather than per-variant model fitting because fitting
                       240,000 logistic regressions is not tractable and the score test
                       is the standard, and near-equivalent, alternative.

PLOTTING IS DELIBERATELY NOT DONE HERE
--------------------------------------
The original wrote Manhattan, QQ, ROC and scatter plots at each site. That has been moved
to the coordinator, which can produce identical plots from the returned summary
statistics -- it receives per-variant beta, se and MAF for every site.

Three reasons, all of which favour the partner:

* It removes matplotlib from the partner's dependency surface entirely.
* Per-site plots written to a partner's filesystem are diagnostics for the coordinator,
  and they were landing where the coordinator could not see them.
* It removes ~100 lines that would otherwise have to be duplicated between this shipped
  module and the coordinator-side aggregator, since neither may import the other.

Numerical outputs are unchanged: the per-site CSVs are still written locally, and the
returned payload is identical.
"""

import json
import os

# Set before numpy is imported, because BLAS reads these at load time. Without them a
# busy shared node can refuse to spawn threads during import ("can't start new thread").
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
from appfl.algorithm.trainer.base_trainer import BaseTrainer  # noqa: E402
from pandas_plink import read_plink1_bin  # noqa: E402
from scipy.stats import chi2, t  # noqa: E402
from sklearn.metrics import r2_score, roc_auc_score  # noqa: E402

DEFAULT_HIT_P_THRESHOLD = 5e-8

# PLINK codes the non-autosomes numerically; some sources use letters.
CHROM_MAP = {"X": 23, "Y": 24, "XY": 25, "MT": 26, "M": 26}


def _normalize_chr(chrom):
    return chrom.astype(str).str.strip().str.upper().replace(CHROM_MAP).astype(int)


def _linear_regression(use_cuml=False, **kwargs):
    """Least squares. cuML is opt-in and lazily imported so it is never required."""
    if use_cuml:
        from cuml.linear_model import LinearRegression
    else:
        from sklearn.linear_model import LinearRegression
    return LinearRegression(**kwargs)


def _fit_binary_model(X, y, use_cuml=False):
    """Unpenalized logistic regression.

    scikit-learn has spelled "no penalty" two different ways across versions -- None in
    current releases, the string "none" in older ones -- and passing the wrong one raises
    rather than falling back. Trying both keeps this working across the range of versions
    a partner might have, which matters because we do not control their environment as
    tightly as our own.
    """
    if use_cuml:
        from cuml.linear_model import LogisticRegression

        # cuML rejects penalty=None; a very large C is the equivalent.
        model = LogisticRegression(fit_intercept=False, C=1e10, max_iter=1000, tol=1e-8)
        model.fit(X, y)
        return model

    from sklearn.linear_model import LogisticRegression

    for penalty in (None, "none"):
        try:
            model = LogisticRegression(
                fit_intercept=False, penalty=penalty, solver="lbfgs", max_iter=1000, tol=1e-8
            )
            model.fit(X, y)
            return model
        except (TypeError, ValueError):
            continue
    raise RuntimeError(
        "could not fit an unpenalized LogisticRegression. Check the scikit-learn "
        "version in the environment the endpoint's worker_init activates."
    )


def _apply_variant_scaling(variant_df, G=None, fraction=1.0, seed=42):
    """Deterministically subsample variants at analysis time.

    Distinct seed from the simulation-time scaling so the two compound independently.
    """
    if fraction >= 1.0:
        return variant_df, G
    n_total = len(variant_df)
    n_keep = max(1, round(n_total * fraction))
    rng = np.random.default_rng(seed)
    chosen = np.sort(rng.choice(n_total, size=n_keep, replace=False))
    filtered = variant_df.iloc[chosen].reset_index(drop=True)
    return filtered, (G.isel(variant=chosen) if G is not None else None)


def _load_genotype_chunk(G, start, end):
    """Load a variant block, mean-impute missing calls, and return dosages plus MAF.

    Chunked because the full matrix does not fit in a worker's memory at real cohort
    sizes, and `single-threaded` because the thread caps above exist for a reason -- dask
    spawning its own pool would defeat them.
    """
    chunk = (
        G.isel(variant=slice(start, end))
        .compute(scheduler="single-threaded")
        .values.astype(np.float64, copy=False)
    )
    observed = np.isfinite(chunk)
    counts = observed.sum(axis=0)
    means = np.divide(
        np.nansum(chunk, axis=0),
        counts,
        out=np.zeros(chunk.shape[1], dtype=np.float64),
        where=counts > 0,
    )
    if not observed.all():
        row_idx, col_idx = np.where(~observed)
        chunk[row_idx, col_idx] = means[col_idx]
    maf = np.minimum(means / 2.0, 1.0 - means / 2.0)
    return chunk, maf


def _write_hits_table(bmi_df, t2d_df, hit_threshold, out_path):
    """Write genome-wide significant variants, or the top 100 if there are none.

    Falling back to the top 100 rather than writing an empty file is deliberate: an empty
    hits table is indistinguishable from a failed run, and the top variants are what
    someone looks at first when a run produces no hits.
    """
    tables = []
    for trait_df in (bmi_df, t2d_df):
        hits = trait_df.loc[trait_df["P"] < hit_threshold].copy()
        if hits.empty:
            hits = trait_df.nsmallest(100, "P").copy()
            hits["HIT_SET"] = "top100"
        else:
            hits = hits.sort_values("P").copy()
            hits["HIT_SET"] = f"p<{hit_threshold:g}"
        tables.append(hits)
    pd.concat(tables, ignore_index=True).to_csv(out_path, index=False)


class SiteGWASTrainer(BaseTrainer):
    """Runs one site's local GWAS and returns summary statistics."""

    def __init__(
        self,
        model=None,
        loss_fn=None,
        metric=None,
        train_dataset=None,
        val_dataset=None,
        train_configs=None,
        logger=None,
        client_id=None,
        **kwargs,
    ):
        super().__init__(
            model=model,
            loss_fn=loss_fn,
            metric=metric,
            train_dataset=train_dataset,
            val_dataset=val_dataset,
            train_configs=train_configs,
            logger=logger,
            client_id=client_id,
            **kwargs,
        )
        self.client_id = str(client_id) if client_id is not None else str(train_dataset.site_id)
        self.chunk_size = int(self.train_configs.get("gwas_chunk_size", 256))
        self.hit_threshold = float(
            self.train_configs.get("hit_p_threshold", DEFAULT_HIT_P_THRESHOLD)
        )
        self.variant_scaling = float(self.train_configs.get("variant_scaling", 1.0))
        self.variant_scaling_seed = int(self.train_configs.get("variant_scaling_seed", 42))
        self.use_cuml = bool(self.train_configs.get("use_cuml", False))

        default_out = str(self.train_dataset.data_dir.parent / "output")
        self.output_dir = Path(
            self.train_configs.get(
                "trainer_output_dirname",
                self.train_configs.get("logging_output_dirname", default_out),
            )
        ).resolve()
        self.data_dir = self.output_dir / "data"
        self.logs_dir = self.output_dir / "logs"
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)

        self.model_state = {}
        self.global_state = {}
        self._has_run = False

    def load_parameters(self, params):
        self.global_state = params if isinstance(params, dict) else {}

    def get_parameters(self):
        return self.model_state

    def train(self, **kwargs):
        # Single round. This is summary-statistic FL: there is one exchange, and a second
        # invocation would recompute identical results at full cost.
        if self._has_run:
            return

        self.logger.info(f"{self.client_id}: loading site data from {self.train_dataset.data_dir}")
        G = read_plink1_bin(str(self.train_dataset.plink_bed), verbose=False, ref="a1")
        (
            sample_df,
            variant_df,
            X_cov,
            y_bmi_gwas,
            y_t2d_gwas,
            y_bmi_eval,
            y_t2d_eval,
        ) = self._load_site_tables(G)

        n_total = len(variant_df)
        variant_df, G = _apply_variant_scaling(
            variant_df, G, self.variant_scaling, self.variant_scaling_seed
        )
        if len(variant_df) < n_total:
            self.logger.info(
                f"{self.client_id}: variant_scaling -> {len(variant_df)} / {n_total} variants"
            )

        bmi_df = self._run_bmi_gwas(G, variant_df, X_cov, y_bmi_gwas)
        t2d_df = self._run_t2d_gwas(G, variant_df, X_cov, y_t2d_gwas)
        pgs_metrics = self._run_local_pgs(
            G, sample_df[["FID", "IID"]], X_cov, y_bmi_eval, y_t2d_eval, bmi_df, t2d_df
        )

        bmi_df.to_csv(self.data_dir / f"{self.client_id}_local_gwas_bmi.csv.gz", index=False)
        t2d_df.to_csv(self.data_dir / f"{self.client_id}_local_gwas_t2d.csv.gz", index=False)
        _write_hits_table(
            bmi_df,
            t2d_df,
            self.hit_threshold,
            self.data_dir / f"{self.client_id}_local_gwas_hits.csv",
        )
        pgs_metrics.to_csv(self.data_dir / f"{self.client_id}_local_pgs_metrics.csv", index=False)

        # How the server learns the variant axis. Sent as JSON bytes because the payload
        # is a tensor dict, and this is the only non-numeric field.
        variant_meta = json.dumps(
            {
                "CHR": variant_df["CHR"].tolist(),
                "SNP": variant_df["SNP"].tolist(),
                "BP": variant_df["BP"].tolist(),
                "EA": variant_df["EA"].tolist(),
                "NEA": variant_df["NEA"].tolist(),
            }
        ).encode("utf-8")

        # THE ONLY THING THAT LEAVES THIS SITE.
        self.model_state = {
            "bmi_beta": torch.from_numpy(bmi_df["BETA"].to_numpy(dtype=np.float64)),
            "bmi_se": torch.from_numpy(bmi_df["SE"].to_numpy(dtype=np.float64)),
            "t2d_beta": torch.from_numpy(t2d_df["BETA"].to_numpy(dtype=np.float64)),
            "t2d_se": torch.from_numpy(t2d_df["SE"].to_numpy(dtype=np.float64)),
            "maf": torch.from_numpy(bmi_df["MAF"].to_numpy(dtype=np.float64)),
            "gwas_n": torch.tensor([len(y_bmi_gwas)], dtype=torch.int64),
            "eval_n": torch.tensor([len(y_bmi_eval)], dtype=torch.int64),
            "local_bmi_r2": torch.tensor([float(pgs_metrics.loc[0, "VALUE"])], dtype=torch.float64),
            "local_t2d_auc": torch.tensor(
                [float(pgs_metrics.loc[1, "VALUE"])], dtype=torch.float64
            ),
            "variant_meta": torch.frombuffer(bytearray(variant_meta), dtype=torch.uint8),
        }

        self._has_run = True
        self.round += 1
        self.logger.info(f"{self.client_id}: local GWAS and PGS complete.")

    def _load_site_tables(self, G):
        """Merge phenotypes and covariates onto the genotype sample order."""
        fam = pd.DataFrame({"FID": G.fid.values.astype(str), "IID": G.iid.values.astype(str)})
        pheno_gwas = pd.read_csv(self.train_dataset.pheno_gwas)
        pheno_eval = pd.read_csv(self.train_dataset.pheno_eval)
        cov = pd.read_csv(self.train_dataset.covariates)

        # validate="one_to_one" is load-bearing: a duplicated FID/IID would silently
        # expand the merge and pair genotypes with the wrong phenotypes.
        sample_df = (
            fam.merge(pheno_gwas, on=["FID", "IID"], how="left", validate="one_to_one")
            .rename(columns={"T2D": "T2D_gwas", "BMI": "BMI_gwas"})
            .merge(pheno_eval, on=["FID", "IID"], how="left", validate="one_to_one")
            .rename(columns={"T2D": "T2D_eval", "BMI": "BMI_eval"})
            .merge(cov, on=["FID", "IID"], how="left", validate="one_to_one")
        )

        required = ["BMI_gwas", "T2D_gwas", "BMI_eval", "T2D_eval", "age", "sex"]
        missing = sample_df[required].isna().sum()
        if missing.any():
            raise ValueError(
                f"{self.client_id}: some samples have no phenotype or covariate values: "
                f"{missing.to_dict()}\n"
                "Every individual in the .fam must appear in all three CSVs. This "
                "usually means the files came from different bundles."
            )

        age = sample_df["age"].to_numpy(dtype=np.float64)
        age_sd = age.std(ddof=0)
        if age_sd == 0:
            raise ValueError(f"{self.client_id}: age has zero variance; cannot standardize.")
        age_z = (age - age.mean()) / age_sd

        X_cov = np.column_stack(
            [
                np.ones(len(sample_df), dtype=np.float64),
                age_z,
                sample_df["sex"].to_numpy(dtype=np.float64),
            ]
        )
        variant_df = pd.DataFrame(
            {
                "CHR": _normalize_chr(pd.Series(G.chrom.values)),
                "SNP": pd.Series(G.snp.values).astype(str),
                "BP": pd.Series(G.pos.values).astype(np.int64),
                "EA": pd.Series(G.a1.values).astype(str),
                "NEA": pd.Series(G.a0.values).astype(str),
            }
        )
        return (
            sample_df,
            variant_df,
            X_cov,
            sample_df["BMI_gwas"].to_numpy(dtype=np.float64),
            sample_df["T2D_gwas"].to_numpy(dtype=np.float64),
            sample_df["BMI_eval"].to_numpy(dtype=np.float64),
            sample_df["T2D_eval"].to_numpy(dtype=np.float64),
        )

    def _run_bmi_gwas(self, G, variant_df, X_cov, y_bmi):
        """Covariate-adjusted linear association, by residualization.

        Residualizing both phenotype and genotypes on the covariates gives the same
        per-variant coefficient as fitting the full model, at a fraction of the cost --
        the covariate fit happens once rather than once per variant.
        """
        self.logger.info(f"{self.client_id}: BMI GWAS")
        n_samples, n_cov = X_cov.shape
        df = n_samples - n_cov - 1

        y_model = _linear_regression(self.use_cuml, fit_intercept=False)
        y_model.fit(X_cov, y_bmi)
        y_res = y_bmi - y_model.predict(X_cov)

        out_chunks = []
        n_variants = len(variant_df)
        for start in range(0, n_variants, self.chunk_size):
            end = min(start + self.chunk_size, n_variants)
            if start == 0 or start % (20 * self.chunk_size) == 0:
                self.logger.info(f"{self.client_id}: BMI chunk {start + 1}-{end}/{n_variants}")

            G_chunk, maf = _load_genotype_chunk(G, start, end)
            g_model = _linear_regression(self.use_cuml, fit_intercept=False)
            g_model.fit(X_cov, G_chunk)
            G_res = G_chunk - g_model.predict(X_cov)

            ss_g = np.einsum("ij,ij->j", G_res, G_res)
            beta = np.divide(
                G_res.T @ y_res,
                ss_g,
                out=np.full(end - start, np.nan, dtype=np.float64),
                where=ss_g > 0,
            )
            resid = y_res[:, None] - G_res * beta[None, :]
            sigma2 = np.einsum("ij,ij->j", resid, resid) / df
            se = np.sqrt(
                np.divide(
                    sigma2,
                    ss_g,
                    out=np.full(end - start, np.nan, dtype=np.float64),
                    where=ss_g > 0,
                )
            )
            stat = np.divide(
                beta,
                se,
                out=np.zeros(end - start, dtype=np.float64),
                where=np.isfinite(se) & (se > 0),
            )
            # Clipped away from zero: a p-value of exactly 0 becomes infinity under
            # -log10 and breaks every downstream plot and threshold comparison.
            p_value = np.clip(2.0 * t.sf(np.abs(stat), df=df), np.finfo(np.float64).tiny, 1.0)

            chunk_df = variant_df.iloc[start:end].copy()
            chunk_df["TRAIT"] = "BMI"
            chunk_df["BETA"] = beta
            chunk_df["SE"] = se
            chunk_df["STAT"] = stat
            chunk_df["P"] = p_value
            chunk_df["MAF"] = maf
            chunk_df["N"] = n_samples
            out_chunks.append(chunk_df)

        return pd.concat(out_chunks, ignore_index=True)

    def _run_t2d_gwas(self, G, variant_df, X_cov, y_t2d):
        """Logistic score test against a covariate-only null.

        Fitting a separate logistic regression per variant is not tractable at genome
        scale. The score test evaluates each variant against one null fit, which is the
        standard approach and near-equivalent at the effect sizes that matter.
        """
        self.logger.info(f"{self.client_id}: T2D GWAS")
        null_model = _fit_binary_model(X_cov, y_t2d, self.use_cuml)
        # Clipped away from 0 and 1: the working weights mu*(1-mu) go to zero otherwise
        # and the information matrix becomes singular.
        mu = np.clip(null_model.predict_proba(X_cov)[:, 1], 1e-8, 1.0 - 1e-8)
        w = mu * (1.0 - mu)
        sqrt_w = np.sqrt(w)
        X_cov_w = X_cov * sqrt_w[:, None]
        score_resid = y_t2d - mu
        n_samples = X_cov.shape[0]

        out_chunks = []
        n_variants = len(variant_df)
        for start in range(0, n_variants, self.chunk_size):
            end = min(start + self.chunk_size, n_variants)
            if start == 0 or start % (20 * self.chunk_size) == 0:
                self.logger.info(f"{self.client_id}: T2D chunk {start + 1}-{end}/{n_variants}")

            G_chunk, maf = _load_genotype_chunk(G, start, end)
            # Weighted residualization: the projection is fit in the weighted space,
            # then applied to the unweighted genotypes.
            g_model = _linear_regression(self.use_cuml, fit_intercept=False)
            g_model.fit(X_cov_w, G_chunk * sqrt_w[:, None])
            G_res = G_chunk - g_model.predict(X_cov)

            info = np.einsum("ij,i,ij->j", G_res, w, G_res)
            beta = np.divide(
                G_res.T @ score_resid,
                info,
                out=np.full(end - start, np.nan, dtype=np.float64),
                where=info > 0,
            )
            se = np.sqrt(
                np.divide(
                    1.0,
                    info,
                    out=np.full(end - start, np.nan, dtype=np.float64),
                    where=info > 0,
                )
            )
            stat = np.divide(
                beta,
                se,
                out=np.zeros(end - start, dtype=np.float64),
                where=np.isfinite(se) & (se > 0),
            )
            p_value = np.clip(chi2.sf(stat * stat, df=1), np.finfo(np.float64).tiny, 1.0)

            chunk_df = variant_df.iloc[start:end].copy()
            chunk_df["TRAIT"] = "T2D"
            chunk_df["BETA"] = beta
            chunk_df["SE"] = se
            chunk_df["STAT"] = stat
            # Clipped before exp() so an extreme beta cannot overflow to inf.
            chunk_df["OR"] = np.exp(np.clip(beta, -50, 50))
            chunk_df["P"] = p_value
            chunk_df["MAF"] = maf
            chunk_df["N"] = n_samples
            chunk_df["TEST"] = "score"
            out_chunks.append(chunk_df)

        return pd.concat(out_chunks, ignore_index=True)

    def _run_local_pgs(self, G, sample_df, X_cov, y_bmi_eval, y_t2d_eval, bmi_df, t2d_df):
        """Score this site's cohort with its own GWAS betas and evaluate.

        Evaluated against an INDEPENDENT phenotype replicate, not the one the betas were
        estimated from. Scoring against the training phenotypes would report a badly
        inflated R2 and AUROC.
        """
        self.logger.info(f"{self.client_id}: local PGS")
        beta_matrix = np.column_stack(
            [
                bmi_df["BETA"].fillna(0.0).to_numpy(dtype=np.float64),
                t2d_df["BETA"].fillna(0.0).to_numpy(dtype=np.float64),
            ]
        )
        scores = np.zeros((G.sizes["sample"], 2), dtype=np.float64)
        n_variants = G.sizes["variant"]

        for start in range(0, n_variants, self.chunk_size):
            end = min(start + self.chunk_size, n_variants)
            weights = beta_matrix[start:end, :]
            if not np.any(weights):
                continue
            G_chunk, _ = _load_genotype_chunk(G, start, end)
            scores += G_chunk @ weights

        score_sd = scores.std(axis=0, ddof=0)
        if np.any(score_sd == 0):
            raise ValueError(
                f"{self.client_id}: a polygenic score has zero variance, so it carries "
                "no information. Usually means every effect estimate was NaN -- check "
                "that the genotypes and phenotypes actually correspond."
            )
        scores = (scores - scores.mean(axis=0)[None, :]) / score_sd[None, :]
        pgs_bmi, pgs_t2d = scores[:, 0], scores[:, 1]

        X_bmi = np.column_stack([X_cov, pgs_bmi])
        bmi_model = _linear_regression(self.use_cuml, fit_intercept=False)
        bmi_model.fit(X_bmi, y_bmi_eval)
        y_bmi_pred = bmi_model.predict(X_bmi)
        bmi_r2 = float(r2_score(y_bmi_eval, y_bmi_pred))

        X_t2d = np.column_stack([X_cov, pgs_t2d])
        t2d_model = _fit_binary_model(X_t2d, y_t2d_eval, self.use_cuml)
        y_t2d_prob = t2d_model.predict_proba(X_t2d)[:, 1]
        t2d_auc = float(roc_auc_score(y_t2d_eval, y_t2d_prob))

        pgs_df = sample_df.copy()
        pgs_df["local_pgs_BMI"] = pgs_bmi
        pgs_df["local_pgs_T2D"] = pgs_t2d
        pgs_df["BMI_observed"] = y_bmi_eval
        pgs_df["BMI_predicted"] = y_bmi_pred
        pgs_df["T2D_observed"] = y_t2d_eval
        pgs_df["T2D_probability"] = y_t2d_prob
        pgs_df.to_csv(self.data_dir / f"{self.client_id}_local_pgs_scores.csv", index=False)

        # Row order is the payload contract: the server reads VALUE at row 0 for BMI R2
        # and row 1 for T2D AUROC.
        return pd.DataFrame(
            [
                {"CLIENT_ID": self.client_id, "TRAIT": "BMI", "METRIC": "R2", "VALUE": bmi_r2},
                {"CLIENT_ID": self.client_id, "TRAIT": "T2D", "METRIC": "AUROC", "VALUE": t2d_auc},
            ]
        )
