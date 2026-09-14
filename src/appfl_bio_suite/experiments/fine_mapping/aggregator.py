"""Pool site aggregates, standardize once, and fine-map with SuSiEx.

COORDINATOR-SIDE. APPFL's ServerAgent loads ``aggregator_path`` locally and never ships
it, so unlike ``trainer.py`` and ``dataset.py`` this module may import freely from the
suite -- and it does, deliberately.

IT DELEGATES THE STATISTICS RATHER THAN REIMPLEMENTING THEM
------------------------------------------------------------
Every number this produces comes out of ``fedfm/fed_fine_mapping.py``, the module
vendored verbatim from the standalone repository: :func:`pool_geno` sums the sites' raw
moments, :func:`build_columns` standardizes them once and writes each ancestry's LD panel,
:func:`pool_pheno` sums the phenotype moments, and :func:`fed_finemap_instance` runs
SuSiEx and scores the credible sets.

This module's whole job is the part APPFL adds: decode the wire payload back into the
aggregate objects that code already accepts, group them, and write the results out. That
split is what makes the two paths comparable -- ``appfl-bio-suite run fine-mapping`` and
``python -m ...fedfm.fed_fine_mapping`` execute the same coordinator statistics over the
same aggregates, so a difference between their outputs is a transport bug and can be
nothing else. ``tests/test_fine_mapping_appfl_parity.py`` asserts they agree.

THE ORDER OF OPERATIONS IS THE WHOLE DESIGN
-------------------------------------------
Sum first, standardize second. Sites ship raw ``X'X`` and ``1'X`` on the 0/1/2 dosage
scale precisely so that the correlation structure is formed against the *pooled* moments.
Standardizing at the sites and averaging the results looks equivalent and is not: each
site would centre against its own column means, deleting the between-site allele-frequency
variation before the coordinator could see it. algo.md A.10 works the example where that
mistake yields 0.866 against a truth of 0.816.

Single round. Sites compute once and report; there is no iterative exchange to converge.

GA4GH: THE TWO HALVES OF THE PROVENANCE MEET HERE
-------------------------------------------------
The coordinator knows what it dispatched -- the data use request, the DRS object it
expected each site to hold, the TRS tool pin it sent. Each site knows what it actually
did, and says so in its manifest. This is the only place both are in one process, so
this is where they are compared and where the merged record is written
(``ga4gh_provenance.json``, beside the results).

A DRS mismatch is fatal here rather than a warning. If a site computed over an object
other than the one this run's provenance claims, every number in the results table is
attributed to the wrong data, and no downstream reader could tell.
"""

from __future__ import annotations

import json
import os
import shutil
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
from appfl.algorithm.aggregator import BaseAggregator  # noqa: E402

from appfl_bio_suite.experiments.fine_mapping.fedfm.fed_fine_mapping import (  # noqa: E402
    GenoAggregate,
    PhenoAggregate,
    build_columns,
    fed_finemap_instance,
    pool_pheno,
)
from appfl_bio_suite.experiments.fine_mapping.fedfm.fine_mapping import (  # noqa: E402
    _resolve_binary,
)
from appfl_bio_suite.experiments.fine_mapping.trainer import (  # noqa: E402
    KEY_GRAM,
    KEY_MANIFEST,
    KEY_SCALARS,
    KEY_SEP,
    KEY_USUM,
    KEY_XTY,
)

__all__ = ["FineMappingAggregator", "split_instance_id"]

# Defaults mirror the standalone driver's argparse, so a run launched through APPFL and a
# run launched from the command line fine-map with the same settings unless told
# otherwise.
DEFAULT_LEVEL = 0.95
DEFAULT_PVAL_THRESH = 1e-5
DEFAULT_MAF = 0.005


def _plain(value):
    """OmegaConf node -> plain dict/list. A no-op for anything already plain."""
    try:
        from omegaconf import OmegaConf

        if OmegaConf.is_config(value):
            return OmegaConf.to_container(value, resolve=True)
    except ImportError:  # pragma: no cover - omegaconf ships with appfl
        pass
    return value


def _to_numpy(value):
    if torch.is_tensor(value):
        return value.detach().cpu().numpy()
    return np.asarray(value)


def split_instance_id(instance: str, locus_id: str) -> tuple[str, int]:
    """``L0000_ncsl1_h2-0.0005_rg1_rep3`` -> ``("ncsl1_h2-0.0005_rg1", 3)``.

    The instance id is the join of three things the simulation knows and the wire format
    does not carry separately. Splitting on the last ``_rep`` rather than on every
    underscore is what keeps architecture ids free to contain them, which they do.
    """
    if not instance.startswith(f"{locus_id}_"):
        raise ValueError(f"instance {instance!r} does not belong to locus {locus_id!r}")
    remainder = instance[len(locus_id) + 1 :]
    arch_id, sep, rep = remainder.rpartition("_rep")
    if not sep or not rep.isdigit():
        raise ValueError(
            f"instance {instance!r} does not end in '_rep<N>'. Instance ids are "
            "'<locus_id>_<architecture_id>_rep<replicate>'."
        )
    return arch_id, int(rep)


class FineMappingAggregator(BaseAggregator):
    """Reassembles per-site aggregates into pooled ancestry columns and fine-maps them."""

    def __init__(self, model=None, aggregator_configs=None, logger=None):
        self.model = model
        self.logger = logger
        self.aggregator_configs = aggregator_configs or {}
        cfg = self.aggregator_configs

        from appfl_bio_suite.experiments.fine_mapping.inference import DEFAULT_OPTIONS

        self.inference_options = {k: cfg.get(k, v) for k, v in DEFAULT_OPTIONS.items()}
        self.population_order = list(
            cfg.get("population_order", ["EUR", "AFR", "AMR", "EAS", "CSA", "MID"])
        )
        self.level = float(cfg.get("level", DEFAULT_LEVEL))
        self.pval_thresh = float(cfg.get("pval_thresh", DEFAULT_PVAL_THRESH))
        self.maf = float(cfg.get("maf", DEFAULT_MAF))
        self.n_workers = int(cfg.get("n_workers", 1))
        self.keep_work = bool(cfg.get("keep_work", False))
        self.make_figures = bool(cfg.get("make_figures", True))

        # What this coordinator dispatched, from federation.yaml via launch.py. Absent
        # for a federation not using GA4GH, in which case the checks below no-op and no
        # provenance file is written.
        #
        # Forced to plain containers: APPFL hands the aggregator an OmegaConf node, and a
        # DictConfig that reads exactly like a dict is not JSON-serializable -- which
        # surfaces at the very end of a run, while writing the provenance, after the
        # expensive part is done.
        self.ga4gh = _plain(cfg.get("ga4gh", {}) or {})

        self.output_dir = Path(cfg.get("output_dir", "local/output/fine-mapping")).resolve()
        self.data_dir = self.output_dir / "data"
        self.graphs_dir = self.output_dir / "graphs"
        self.logs_dir = self.output_dir / "logs"
        self.work_root = Path(cfg.get("work_dir", self.output_dir / "work")).resolve()
        for directory in (self.data_dir, self.graphs_dir, self.logs_dir, self.work_root):
            directory.mkdir(parents=True, exist_ok=True)

        # Ground truth. Coordinator-only by design -- the causal manifest is the answer
        # key and is deliberately not distributed with the site bundles. Without it the
        # run still produces credible sets; it just cannot score them, which is the right
        # behaviour for a real federation where nobody knows the truth.
        truth = cfg.get("causal_manifest", None)
        self.causal_manifest = Path(truth).resolve() if truth else None

        # SuSiEx is a C++ CLI we shell out to. PLINK is demanded by its argument
        # validator (main.cpp:274) but never invoked once precomputed LD files are
        # present, which on this path they always are -- a coordinator holds no genotypes
        # and so has no reason to have PLINK at all.
        self.vendor_root = Path(cfg.get("vendor_root", ".")).resolve()
        self.susiex = str(cfg.get("susiex_binary") or _resolve_binary("SuSiEx", self.vendor_root))
        try:
            self.plink = str(cfg.get("plink_binary") or _resolve_binary("plink", self.vendor_root))
        except FileNotFoundError:
            self.plink = "plink"

        self.global_state = {"fine_mapping_ready": torch.tensor([0], dtype=torch.int64)}

    def get_parameters(self, **kwargs):
        return self.global_state

    # -- decoding ----------------------------------------------------------

    def _decode(self, client_ids, local_models):
        """Turn the flat tensor payloads back into the aggregate objects.

        Returns ``(geno_by_locus, pheno_by_instance, loci, site_manifests)``, keyed so
        that the pooling below never has to search.
        """
        geno_by_locus: dict[str, list[GenoAggregate]] = {}
        pheno_by_instance: dict[str, dict[str, list[PhenoAggregate]]] = {}
        loci: dict[str, pd.Series] = {}
        manifests: dict[str, dict] = {}

        for cid in client_ids:
            payload = local_models[cid]
            if KEY_MANIFEST not in payload:
                raise ValueError(
                    f"site {cid} returned no '{KEY_MANIFEST}' block. Its payload cannot "
                    "be interpreted -- the site is running a different trainer version."
                )
            manifest = json.loads(_to_numpy(payload[KEY_MANIFEST]).tobytes().decode("utf-8"))
            manifests[cid] = manifest

            for entry in manifest["loci"]:
                # Every site must agree on a locus's window, or they are not fine-mapping
                # the same region. Checked rather than assumed because the windows come
                # from each site's own bundle.
                as_series = pd.Series(
                    {
                        "locus_id": entry["locus_id"],
                        "chrom": int(entry["chrom"]),
                        "start_bp": int(entry["start_bp"]),
                        "end_bp": int(entry["end_bp"]),
                        "stratum": entry.get("stratum", ""),
                    }
                )
                known = loci.get(entry["locus_id"])
                if known is not None and not known.drop("stratum").equals(
                    as_series.drop("stratum")
                ):
                    raise ValueError(
                        f"site {cid} defines locus {entry['locus_id']} as "
                        f"chr{as_series['chrom']}:{as_series['start_bp']}-"
                        f"{as_series['end_bp']}, which disagrees with another site.\n"
                        "Every site's selected_loci.tsv must come from the same "
                        "simulation run."
                    )
                loci[entry["locus_id"]] = as_series

            for block in manifest["geno_blocks"]:
                locus_id, pop = block["locus_id"], block["pop"]
                variants = pd.DataFrame(block["variants"])
                variants["bp"] = variants["bp"].astype(np.int64)
                key = f"{KEY_GRAM}{KEY_SEP}{locus_id}{KEY_SEP}{pop}"
                # Upcast unconditionally: a site may have shipped the Gram as float32
                # where that was provably lossless, and the pooling sums in float64.
                G = _to_numpy(payload[key]).astype(np.float64, copy=False)
                u = _to_numpy(payload[f"{KEY_USUM}{KEY_SEP}{locus_id}{KEY_SEP}{pop}"])
                geno_by_locus.setdefault(locus_id, []).append(
                    GenoAggregate(
                        site=str(cid),
                        pop=pop,
                        locus_id=locus_id,
                        variants=variants,
                        G=G,
                        u=np.asarray(u, dtype=np.float64),
                        n=int(block["n"]),
                        n_incomplete_variants=int(block["n_incomplete_variants"]),
                    )
                )

            # snp_ids for a phenotype block are the geno block's, by construction at the
            # site: the same rows and the same surviving columns produced both.
            snp_ids_of = {
                (b["locus_id"], b["pop"]): np.asarray(b["variants"]["snp_id"], dtype=object)
                for b in manifest["geno_blocks"]
            }
            for block in manifest["pheno_blocks"]:
                inst, pop = block["instance"], block["pop"]
                c = _to_numpy(payload[f"{KEY_XTY}{KEY_SEP}{inst}{KEY_SEP}{pop}"])
                q, w, n = _to_numpy(payload[f"{KEY_SCALARS}{KEY_SEP}{inst}{KEY_SEP}{pop}"])
                pheno_by_instance.setdefault(inst, {}).setdefault(pop, []).append(
                    PhenoAggregate(
                        site=str(cid),
                        pop=pop,
                        instance=inst,
                        snp_ids=snp_ids_of[(block["locus_id"], pop)],
                        c=np.asarray(c, dtype=np.float64),
                        q=float(q),
                        w=float(w),
                        n=int(n),
                    )
                )

        return geno_by_locus, pheno_by_instance, loci, manifests

    def _truth_map(self) -> dict[tuple[str, str, int], list[str]]:
        """Ground-truth causal variants, keyed by (locus, architecture, replicate)."""
        if self.causal_manifest is None:
            self.logger.warning(
                "no `causal_manifest` configured: credible sets will be produced but not "
                "scored against the truth, so coverage and PIP columns will be empty. "
                "That is correct for a real federation and wrong for a benchmark -- point "
                "aggregator_kwargs.causal_manifest at the simulation's "
                "causal_manifest.tsv if this is one."
            )
            return {}
        if not self.causal_manifest.is_file():
            raise FileNotFoundError(
                f"causal_manifest {self.causal_manifest} does not exist. It is written by "
                "`appfl-bio-suite simulate fine-mapping` into the run's ground_truth/ "
                "directory."
            )
        manifest = pd.read_csv(self.causal_manifest, sep="\t")
        return {
            (r.locus_id, r.architecture_id, int(r.replicate)): str(r.causal_snp_ids).split(",")
            for r in manifest.itertuples(index=False)
        }

    # -- the exchange ------------------------------------------------------

    def aggregate(self, local_models, **kwargs):
        started = time.time()
        client_ids = list(local_models.keys())
        self.logger.info(f"fine-mapping across {len(client_ids)} site(s)")

        geno_by_locus, pheno_by_instance, loci, manifests = self._decode(client_ids, local_models)
        self._check_ga4gh(client_ids, manifests)
        truth_map = self._truth_map()

        # The ancestry columns to build, in a stable order. Union across sites, ordered by
        # the first site's declared order so a rerun writes its columns the same way --
        # SuSiEx's output is per-column and an unstable order makes two runs' files
        # gratuitously different.
        pops: list[str] = []
        for cid in client_ids:
            for pop in manifests[cid]["pops"]:
                if pop not in pops:
                    pops.append(pop)

        pops = [p for p in self.population_order if p in pops] + sorted(
            set(pops) - set(self.population_order)
        )
        self._log_uplink(client_ids, manifests)

        rows: list[dict] = []
        for locus_id in sorted(geno_by_locus):
            locus = loci[locus_id]
            blocks = geno_by_locus[locus_id]
            locus_dir = self.work_root / f"{locus_id}_ld"
            locus_dir.mkdir(parents=True, exist_ok=True)

            # Sum, then standardize once, then write each ancestry's LD panel. Genotype
            # only, so it is built once per locus and shared by every instance there.
            columns = build_columns(blocks, pops, locus_dir, self.maf)
            for col in columns:
                if col.n_incomplete_variants:
                    self.logger.warning(
                        f"  {locus_id}/{col.pop}: {col.n_incomplete_variants} window "
                        "variant(s) dropped as not fully observed at some site; the "
                        "centralized comparator applies the same complete-variant policy"
                    )
            self.logger.info(
                f"  {locus_id}: {len(blocks)} block(s) from "
                f"{len({b.site for b in blocks})} site(s) -> "
                + ", ".join(f"{c.pop}(n={c.n}, M={len(c.variants)})" for c in columns)
            )

            instances = sorted(
                inst for inst in pheno_by_instance if inst.startswith(f"{locus_id}_")
            )
            jobs = []
            for inst in instances:
                arch_id, rep = split_instance_id(inst, locus_id)
                pooled = {
                    col.pop: pool_pheno(
                        pheno_by_instance[inst][col.pop],
                        col.variants["snp_id"].to_numpy(),
                    )
                    for col in columns
                    if col.pop in pheno_by_instance[inst]
                }
                missing = [c.pop for c in columns if c.pop not in pooled]
                if missing:
                    raise ValueError(
                        f"instance {inst} has no phenotype aggregates for ancestry "
                        f"column(s) {missing}, but genotype aggregates exist for them. "
                        "A site returned one half of a block and not the other."
                    )
                jobs.append((arch_id, rep, truth_map.get((locus_id, arch_id, rep), []), pooled))

            rows.extend(self._run_instances(locus, jobs, columns))
            if not self.keep_work:
                shutil.rmtree(locus_dir, ignore_errors=True)

        results = pd.DataFrame(rows)
        self._write_outputs(results, client_ids, manifests)

        elapsed = time.time() - started
        self.global_state = {
            "fine_mapping_ready": torch.tensor([1], dtype=torch.int64),
            "num_clients": torch.tensor([len(client_ids)], dtype=torch.int64),
            "num_loci": torch.tensor([len(geno_by_locus)], dtype=torch.int64),
            "num_instances": torch.tensor([len(results)], dtype=torch.int64),
            "power_any_causal": torch.tensor(
                [float(results["any_causal_captured"].mean()) if len(results) else np.nan],
                dtype=torch.float64,
            ),
            "mean_credible_sets": torch.tensor(
                [float(results["n_credible_sets"].mean()) if len(results) else np.nan],
                dtype=torch.float64,
            ),
            "elapsed_s": torch.tensor([elapsed], dtype=torch.float64),
        }
        self.logger.info(
            f"fine-mapping complete in {elapsed:.1f}s -> {self.output_dir}/{{data,graphs}}"
        )
        return self.global_state

    def _run_instances(self, locus, jobs, columns) -> list[dict]:
        """Fine-map every instance at one locus, in parallel where asked.

        ``fed_finemap_instance`` catches its own failures and records them in the result
        row, so one instance that SuSiEx chokes on does not lose the other 149.
        """
        if not jobs:
            return []

        def call(arch_id, rep, truth, pooled):
            return fed_finemap_instance(
                locus,
                arch_id,
                rep,
                truth,
                columns,
                pooled,
                self.work_root,
                self.susiex,
                self.plink,
                self.level,
                self.pval_thresh,
                self.keep_work,
                self.inference_options,
            )

        if self.n_workers <= 1:
            return [call(*job) for job in jobs]

        from joblib import Parallel, delayed

        return list(
            Parallel(n_jobs=self.n_workers, backend="loky")(
                delayed(fed_finemap_instance)(
                    locus,
                    arch_id,
                    rep,
                    truth,
                    columns,
                    pooled,
                    self.work_root,
                    self.susiex,
                    self.plink,
                    self.level,
                    self.pval_thresh,
                    self.keep_work,
                    self.inference_options,
                )
                for arch_id, rep, truth, pooled in jobs
            )
        )

    def _log_uplink(self, client_ids, manifests) -> None:
        """Report what each site actually sent. The payload size is this experiment's
        headline cost, so it belongs in the log rather than in a reader's head."""
        for cid in client_ids:
            manifest = manifests[cid]
            self.logger.info(
                f"  {cid}: n={manifest['sample_size']:,}, "
                f"ancestries {sorted(manifest['composition'])}, "
                f"{len(manifest['geno_blocks'])} genotype block(s), "
                f"{len(manifest['pheno_blocks'])} phenotype block(s), "
                f"{manifest['n_flipped']} variant(s) recoded to reference allele order"
            )

    # -- GA4GH -------------------------------------------------------------

    def _check_ga4gh(self, client_ids, manifests) -> None:
        """Compare what each site reports doing against what this run dispatched.

        Three comparisons, in decreasing severity:

        * **DRS object identity.** A site that verified its bundle against a different
          object than this run's provenance names has computed over different data than
          the results will claim. Fatal.
        * **Data use.** A site whose decision was not ``permitted`` should never have
          reached this point -- enforcement is on by default -- so arriving here means it
          was explicitly disabled. The results are still produced, and the fact is
          recorded and logged, because a result computed under a waived consent check
          must not look identical to one computed under a satisfied check.
        * **Tool pin.** A site that echoes a different pin than the coordinator holds was
          dispatched by something else, which usually means two drivers ran against one
          federation.
        """
        expected_objects = (self.ga4gh.get("drs") or {}).get("objects", {}) or {}
        expected_tool = self.ga4gh.get("tool") or {}

        for cid in client_ids:
            reported = manifests[cid].get("ga4gh") or {}
            if not reported:
                if self.ga4gh:
                    self.logger.warning(
                        f"  {cid}: returned no GA4GH provenance. It is running a trainer "
                        "from before this was recorded; its results cannot be attributed "
                        "to a tool version or a data object."
                    )
                continue

            drs = reported.get("drs") or {}
            expected_uri = expected_objects.get(str(cid))
            actual_uri = drs.get("self_uri") or ""
            if expected_uri and actual_uri and expected_uri != actual_uri:
                raise ValueError(
                    f"site {cid} computed over DRS object {actual_uri}, but this run "
                    f"expects {expected_uri}.\n"
                    "Every number this site contributed would be attributed to data it "
                    "did not read. Fix federation.yaml's `drs_uri` for this site, or "
                    "re-send the bundle this run is about."
                )

            data_use = reported.get("data_use") or {}
            status = data_use.get("status") or data_use.get("outcome")
            if status in ("denied", "undetermined"):
                self.logger.warning(
                    f"  {cid}: DATA USE {status.upper()} -- this site's terms do not "
                    "permit this study, and enforcement was disabled. Its contribution "
                    "is in these results and the decision is recorded in "
                    "ga4gh_provenance.json."
                )
            elif status == "no-profile":
                self.logger.info(f"  {cid}: no data use profile in its bundle; nothing checked")
            elif status == "permitted":
                self.logger.info(
                    f"  {cid}: data use permitted "
                    f"({len(data_use.get('reasons', []))} DUO term(s) satisfied)"
                )

            tool = reported.get("tool") or {}
            if expected_tool and tool and tool.get("id") != expected_tool.get("id"):
                self.logger.warning(
                    f"  {cid}: ran tool {tool.get('id')}@{tool.get('version')}, but this "
                    f"run pins {expected_tool.get('id')}@{expected_tool.get('version')}. "
                    "Two drivers against one federation is the usual cause."
                )

    def _write_ga4gh_provenance(self, results_path: Path, client_ids, manifests) -> Path | None:
        """Merge dispatched intent with site attestation, and address the outputs.

        Written even when a site reported nothing, because "this site attested to
        nothing" is itself the fact a reader needs. Skipped entirely when the run carries
        no GA4GH configuration at all.
        """
        if not self.ga4gh:
            return None

        from appfl_bio_suite.core.ga4gh.drs import registry_for_files

        record = {
            "run": {
                "experiment": "fine-mapping",
                "finished_at": pd.Timestamp.utcnow().isoformat(),
                "sites": len(client_ids),
            },
            "dispatched": self.ga4gh,
            "attested": {
                str(cid): (manifests[cid].get("ga4gh") or {"status": "not-reported"})
                for cid in client_ids
            },
        }

        hostname = (self.ga4gh.get("drs") or {}).get("hostname")
        if hostname:
            # The two provenance files are about the outputs rather than among them, and
            # neither can contain its own checksum. Including them would also make a
            # rerun's registry list the previous run's leftovers -- the same trap
            # core/simulation.py's checksum_tree avoids by excluding the run manifest.
            outputs = sorted(
                path
                for path in self.data_dir.glob("*")
                if path.is_file() and path.name not in ("drs_outputs.json", "ga4gh_provenance.json")
            )
            registry = registry_for_files(outputs, hostname, description="fine-mapping results")
            registry.save(self.data_dir / "drs_outputs.json")
            record["outputs"] = {
                obj.name: obj.self_uri for obj in registry.objects.values() if not obj.is_bundle
            }

        path = self.data_dir / "ga4gh_provenance.json"
        path.write_text(json.dumps(record, indent=2, sort_keys=False) + "\n", encoding="utf-8")
        self.logger.info(f"wrote GA4GH provenance -> {path}")
        return path

    def _write_outputs(self, results: pd.DataFrame, client_ids, manifests) -> None:
        results_path = self.data_dir / "fed_fm_results.tsv"
        results.to_csv(results_path, sep="\t", index=False)
        self.logger.info(f"wrote {len(results)} result row(s) -> {results_path}")

        pd.DataFrame(
            [
                {
                    "CLIENT_ID": cid,
                    "N_SAMPLES": manifests[cid]["sample_size"],
                    "ANCESTRIES": ",".join(sorted(manifests[cid]["composition"])),
                    "N_GENO_BLOCKS": len(manifests[cid]["geno_blocks"]),
                    "N_PHENO_BLOCKS": len(manifests[cid]["pheno_blocks"]),
                    "N_VARIANTS_FLIPPED": manifests[cid]["n_flipped"],
                }
                for cid in client_ids
            ]
        ).to_csv(self.data_dir / "fed_fm_site_summary.csv", index=False)

        if results.empty:
            self.logger.warning("no instances were fine-mapped; skipping rollups")
            self._write_ga4gh_provenance(results_path, client_ids, manifests)
            return

        self._write_rollups(results)
        # After the rollups, so the output objects it addresses include them. A
        # provenance record that names half the outputs is worse than none: it reads as
        # complete.
        self._write_ga4gh_provenance(results_path, client_ids, manifests)
        if self.make_figures:
            from appfl_bio_suite.experiments.fine_mapping.plotting import write_figures

            write_figures(results, self.graphs_dir, self.logger)

    def _write_rollups(self, results: pd.DataFrame) -> None:
        """Aggregate power and credible-set metrics, the same way the standalone path
        does, so the two tables line up column for column."""
        meta = results["architecture_id"].map(_parse_architecture).apply(pd.Series)
        frame = pd.concat([results, meta], axis=1)

        def rollup(group: pd.DataFrame) -> pd.Series:
            return pd.Series(
                {
                    "n_instances": len(group),
                    "power_any_causal": group["any_causal_captured"].mean(),
                    "mean_n_cs": group["n_credible_sets"].mean(),
                    "median_best_cs_size": group["best_cs_size"].median(),
                    "mean_causal_pip": group["causal_pip_mean"].mean(),
                }
            )

        for by, name in (
            (["ncsl", "h2_target", "rg"], "by_architecture"),
            (["stratum", "rg"], "by_stratum_rg"),
        ):
            cols = [c for c in by if c in frame.columns and frame[c].notna().any()]
            if not cols:
                continue
            path = self.data_dir / f"fed_fm_rollup_{name}.tsv"
            (
                frame.groupby(cols)
                .apply(rollup, include_groups=False)
                .reset_index()
                .to_csv(path, sep="\t", index=False)
            )
            self.logger.info(f"wrote rollup -> {path}")


def _parse_architecture(architecture_id: str) -> dict:
    """``ncsl2_h2-0.001_rg0.7`` -> its three parameters.

    Parsed from the id rather than looked up, because the id is the only thing the site
    payload carries and the coordinator should not need the simulation's architecture
    grid in hand to roll up its own results. Unparseable ids yield NaNs and drop out of
    the grouped rollups rather than failing the run.
    """
    out = {"ncsl": np.nan, "h2_target": np.nan, "rg": np.nan}
    for part in str(architecture_id).split("_"):
        try:
            if part.startswith("ncsl"):
                out["ncsl"] = int(part[4:])
            elif part.startswith("h2-"):
                out["h2_target"] = float(part[3:])
            elif part.startswith("rg"):
                out["rg"] = float(part[2:])
        except ValueError:
            continue
    return out
