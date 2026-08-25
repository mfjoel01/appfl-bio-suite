"""The fine-mapping SITE stage. SHIPPED TO WORKERS.

=============================================================================
THIS FILE'S SOURCE IS SENT OVER THE WIRE AND EXECUTED ON A PARTNER'S CLUSTER.
=============================================================================

Self-contained by design: standard library, numpy/pandas/scipy, torch, and ``appfl``
only. No imports from ``appfl_bio_suite``, no sibling imports, no relative imports.
Enforced by tests/test_shipped_modules.py.

WHAT THIS COMPUTES, AND WHAT LEAVES
-----------------------------------
For each locus in this run's shard, and for each superpopulation this site holds, it
reads its own genotype window once and emits the raw second-moment aggregates of
``algo.md`` §7, on the **unstandardized** 0/1/2 dosage scale::

    G = X'X   (M x M)      u = 1'X   (M)      n           <- once per (locus, ancestry)
    c = X'y   (M)          q = 1'y            w = y'y     <- once per instance

Nothing else crosses the site boundary. No genotype row, no phenotype value, and in
particular **no locally standardized LD matrix**: standardizing here and combining the
results at the coordinator is biased, because each site would centre against its own
column means and delete the between-site allele-frequency variation before the
coordinator could see it (algo.md §9). Standardization is the coordinator's job and doing
it here is the one mistake this whole protocol is shaped to prevent.

The split between the two groups above is an economy of transmission, not of information.
Only ``c``, ``q`` and ``w`` depend on the phenotype, so the O(M^2) half is emitted once
per locus and the O(M) half once per instance; together they are exactly the six objects
of §7.

WHY THIS CODE EXISTS TWICE
--------------------------
``fedfm/fed_fine_mapping.py::site_aggregates`` is the same computation, vendored verbatim
from the standalone repository, and the coordinator-side pooling calls into it. This
module cannot: shipped source may not import from the suite, so the site stage has to
stand alone here.

The duplication is deliberate and it is tested rather than trusted --
``tests/test_fine_mapping_shipped_parity.py`` runs both implementations over the same
fixture and asserts ``G``, ``u``, ``c``, ``q``, ``w`` and ``n`` agree exactly, not
approximately. If you change the statistics in one, that test fails until you change the
other.

ALLELE HARMONIZATION IS NOT OPTIONAL HERE
-----------------------------------------
Theorem 1 condition (i) wants one ordered, allele-harmonized variant list across sites,
and it does not come free. Per-site filesets cut without ``--keep-allele-order`` have A1
set to each site's own minor allele, and on real chr1 data about 6% of variants end up
coded oppositely at two sites. A flipped site contributes ``2 - x`` where the others
contribute ``x``, which sums in silently and negates every off-diagonal that variant
touches. So this recodes against the bundle's ``reference_variants.tsv`` before forming a
single moment, and reports how many it had to flip; the coordinator re-checks rather than
trusting the report.

MISSING GENOTYPES
-----------------
Theorem 1 condition (ii) needs one complete-case cohort per ancestry. A window variant not
fully observed within a block is dropped from that block and counted, because no single
``(G, u, n)`` triple can express per-variant complete-case handling. On HAPNEST chr1 this
costs about four variants in two thousand. Where the count is nonzero the federated and
centralized paths fine-map slightly different variant sets and the exactness guarantee no
longer bites -- which is why the number is reported rather than silently absorbed.
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

# ---------------------------------------------------------------------------
# The uplink format. The aggregator imports these names -- it is coordinator-side and may
# import freely -- so the wire contract is single-sourced here rather than described in
# two places that can drift.
# ---------------------------------------------------------------------------

# Payload keys are "<kind>::<scope>::<pop>". A flat tensor dict is what APPFL's
# communicators serialize, so structure has to live in the key and in the manifest.
KEY_SEP = "::"
KEY_GRAM = "G"        # G::<locus_id>::<pop>    -- (M, M) genotype Gram
KEY_USUM = "u"        # u::<locus_id>::<pop>    -- (M,)   column sums
KEY_XTY = "c"         # c::<instance>::<pop>    -- (M,)   X'y
KEY_SCALARS = "s"     # s::<instance>::<pop>    -- (3,)   (q, w, n)
KEY_MANIFEST = "manifest"

# ``X'X`` entries are sums of products of dosages in {0,1,2}, so they are integers
# bounded by 4n. float32 holds every integer below 2^24 exactly, so at any realistic
# cohort size the Gram survives a float32 round trip bit-for-bit -- and it is the
# dominant term in the uplink, so halving it halves the run's network cost. Guarded
# rather than assumed: the trainer verifies the bound before downcasting and falls back
# to float64 otherwise. See `uplink_gram_dtype` in the client config.
_EXACT_FLOAT32_INT = 1 << 24

# Rows read per pass when accumulating X'y in float64 out of a float32 dosage block: big
# enough to keep BLAS busy, small enough that the upcast temporary stays ~100 MB at
# locus-scale M. Matches fedfm/fed_fine_mapping.py, and must, for bit-identical sums:
# floating-point addition is not associative, so the chunk size is part of the answer.
_ROW_CHUNK = 8192

# Same bound, for choosing whether the Gram can be formed in single precision. Above it,
# partial sums could exceed float32's exact-integer range and we pay for float64.
_EXACT_GRAM_MAX_N = (1 << 24) // 4

# PLINK 1 .bed two-bit code -> dosage of the A1 allele.
# 00 = hom A1 (2), 01 = missing, 10 = het (1), 11 = hom A2 (0).
_BED_CODE_TO_DOSAGE = np.array([2, np.nan, 1, 0], dtype=np.float32)

_BIM_COLUMNS = ["chrom", "snp_id", "cm", "bp", "a1", "a2"]
_FAM_COLUMNS = ["FID", "IID", "father", "mother", "sex", "pheno"]


def block_key(kind, scope, pop):
    """Build a payload key. The aggregator splits on :data:`KEY_SEP` to invert it."""
    return f"{kind}{KEY_SEP}{scope}{KEY_SEP}{pop}"


# ---------------------------------------------------------------------------
# PLINK readers -- inlined because this module may not import its siblings
# ---------------------------------------------------------------------------


def read_bim(path):
    """Read a PLINK .bim. Row index is the variant's position inside the .bed."""
    return pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=_BIM_COLUMNS,
        dtype={"chrom": str, "snp_id": str, "cm": float, "bp": np.int64,
               "a1": str, "a2": str},
        engine="python",
    )


def read_fam(path):
    return pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=_FAM_COLUMNS,
        dtype={"FID": str, "IID": str, "father": str, "mother": str,
               "sex": str, "pheno": str},
        engine="python",
    )


def read_bed_variants(bed_prefix, variant_indices, n_samples):
    """Read a subset of variants from a PLINK 1 .bed as an (n_samples, M) dosage matrix.

    Seeks straight to each requested variant rather than streaming the file, which is
    what makes a locus window cheap to read out of a chromosome-scale fileset. Missing
    calls come back as NaN. float32, because a dosage is one of four values and the
    matrix is the largest thing on the worker.

    Only variant-major (mode byte 0x01) files are supported. PLINK 1.9+ always writes
    variant-major; this is the universal format.
    """
    bed_path = Path(str(bed_prefix) + ".bed")
    if not bed_path.exists():
        raise FileNotFoundError(f"No .bed file at {bed_path}")

    bytes_per_variant = (n_samples + 3) // 4
    variant_indices = np.asarray(variant_indices, dtype=np.int64)
    out = np.empty((n_samples, len(variant_indices)), dtype=np.float32)

    with bed_path.open("rb") as handle:
        magic = handle.read(3)
        if len(magic) != 3:
            raise ValueError(f"Truncated .bed header in {bed_path}")
        if magic[0] != 0x6C or magic[1] != 0x1B:
            raise ValueError(f"Bad .bed magic in {bed_path}: {magic!r}")
        if magic[2] != 0x01:
            raise ValueError(
                f"Sample-major .bed not supported ({bed_path}); "
                "re-encode with `plink --bfile <x> --make-bed --out <x>`"
            )

        for out_col, vidx in enumerate(variant_indices):
            handle.seek(3 + int(vidx) * bytes_per_variant)
            block = handle.read(bytes_per_variant)
            if len(block) != bytes_per_variant:
                raise ValueError(
                    f"Truncated variant {vidx} in {bed_path}: expected "
                    f"{bytes_per_variant} bytes, got {len(block)}"
                )
            arr = np.frombuffer(block, dtype=np.uint8)
            codes = np.empty(bytes_per_variant * 4, dtype=np.uint8)
            codes[0::4] = arr & 0x03
            codes[1::4] = (arr >> 2) & 0x03
            codes[2::4] = (arr >> 4) & 0x03
            codes[3::4] = (arr >> 6) & 0x03
            out[:, out_col] = _BED_CODE_TO_DOSAGE[codes[:n_samples]]

    return out


# ---------------------------------------------------------------------------
# The moments themselves
# ---------------------------------------------------------------------------


def columns_with_missing(X, chunk=_ROW_CHUNK):
    """Boolean mask of columns carrying at least one missing dosage.

    Chunked so the intermediate boolean array stays small beside a locus-scale block.
    """
    bad = np.zeros(X.shape[1], dtype=bool)
    for start in range(0, X.shape[0], chunk):
        bad |= np.isnan(X[start:start + chunk]).any(axis=0)
    return bad


def gram(X):
    """``X'X`` as exact integers in float64.

    Formed in single precision while ``4n`` stays under float32's exact-integer ceiling:
    every partial sum a BLAS sgemm produces is itself an integer below that bound, so the
    result is exact and costs half the memory. The integrality assertion is not a
    formality -- it is what catches a block that still carries missing values or is not
    on the raw 0/1/2 scale, either of which would corrupt every downstream moment
    silently.
    """
    n = X.shape[0]
    work = X if n <= _EXACT_GRAM_MAX_N else X.astype(np.float64)
    G = np.asarray(work.T @ work, dtype=np.float64)
    if not np.array_equal(G, np.rint(G)):
        raise RuntimeError(
            "X'X is not integral -- the dosage block is not raw 0/1/2 counts, "
            "or it still carries missing values"
        )
    return G


def xty(X, Y, chunk=_ROW_CHUNK):
    """``X'Y`` accumulated in float64 out of a float32 dosage block.

    In row chunks rather than by upcasting the whole block: single precision would cost
    ~1e-5 relative on a 50,000-row sum, which is coarser than the six significant digits
    the summary statistics are written at.
    """
    out = np.zeros((X.shape[1], Y.shape[1]), dtype=np.float64)
    for start in range(0, X.shape[0], chunk):
        stop = start + chunk
        out += X[start:stop].astype(np.float64).T @ Y[start:stop]
    return out


def locus_window_variants(bim, chrom, start_bp, end_bp):
    """In-window variants of a ``.bim``, in ``.bed`` order, index preserved.

    The preserved index is the variant's row inside the matching ``.bed``, which is what
    :func:`read_bed_variants` wants.
    """
    mask = (bim["chrom"].astype(str) == str(chrom)) & bim["bp"].between(
        int(start_bp), int(end_bp)
    )
    return bim.loc[mask]


def harmonize_bim(bim, reference, site_id):
    """Recode this site's ``.bim`` to the reference allele order; return (bim, flip).

    ``flip`` marks the variants this site has A1/A2 the other way round from the agreed
    reference list. Their dosages become ``2 - x`` before any moment is formed, and the
    returned bim carries the canonical alleles so that what the site emits is already
    harmonized. A variant whose alleles the reference does not list at all is a hard
    error rather than a flip: that is not a coding difference, it is a different variant.
    """
    ref = reference.reindex(bim["snp_id"].to_numpy())
    if ref["a1"].isna().any():
        n_absent = int(ref["a1"].isna().sum())
        example = bim.loc[ref["a1"].isna().to_numpy(), "snp_id"].to_numpy()[:3]
        raise RuntimeError(
            f"{site_id} holds {n_absent} variant(s) absent from reference_variants.tsv, "
            f"e.g. {list(example)}.\n"
            "The reference variant list is agreed before any data is touched and must "
            "cover every variant a site holds; a site with extra variants was cut from a "
            "different fileset than the one the list describes."
        )
    a1, a2 = bim["a1"].to_numpy(), bim["a2"].to_numpy()
    r1, r2 = ref["a1"].to_numpy(), ref["a2"].to_numpy()
    same = (a1 == r1) & (a2 == r2)
    flip = (a1 == r2) & (a2 == r1)
    if not (same | flip).all():
        bad = bim.loc[~(same | flip), "snp_id"].to_numpy()
        raise RuntimeError(
            f"{site_id} codes {len(bad)} variant(s) with alleles the reference does not "
            f"list at all, e.g. {bad[0]}: not a strand or order swap."
        )
    return bim.assign(a1=r1, a2=r2), flip


class SiteFineMappingTrainer(BaseTrainer):
    """Computes this site's aggregates for the run's loci and returns only those."""

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

        # Which loci this exchange covers. Sharding is by locus because the uplink is
        # O(M^2) per (locus, ancestry) and a full 100-locus sweep is hundreds of
        # megabytes per site. Every site must shard identically or the coordinator would
        # be pooling different loci -- which is why these come from the server config
        # rather than from anything site-local.
        self.locus_n_shards = int(self.train_configs.get("locus_n_shards", 1))
        self.locus_shard_index = int(self.train_configs.get("locus_shard_index", 0))
        self.locus_limit = self.train_configs.get("locus_limit", None)
        self.instance_limit = self.train_configs.get("instance_limit", None)

        # Ancestry columns the coordinator asked for. A site simply skips any it does not
        # hold; the coordinator pools whoever contributed.
        self.pops = list(self.train_configs.get("pops", []) or [])

        self.gram_dtype = str(self.train_configs.get("uplink_gram_dtype", "float64"))
        if self.gram_dtype not in ("float64", "float32"):
            raise ValueError(
                f"uplink_gram_dtype must be 'float64' or 'float32', got "
                f"{self.gram_dtype!r}"
            )

        default_out = str(self.train_dataset.data_dir.parent / "output")
        self.output_dir = Path(
            self.train_configs.get(
                "trainer_output_dirname",
                self.train_configs.get("logging_output_dirname", default_out),
            )
        ).resolve()
        self.logs_dir = self.output_dir / "logs"
        self.logs_dir.mkdir(parents=True, exist_ok=True)

        self.model_state = {}
        self.global_state = {}
        self._has_run = False

    def load_parameters(self, params):
        self.global_state = params if isinstance(params, dict) else {}

    def get_parameters(self):
        return self.model_state

    # -- helpers ----------------------------------------------------------

    def _selected_loci(self):
        """The loci this run covers, after sharding and any limit.

        Sharded by ``index % n_shards`` over the bundle's own locus order, which every
        site shares because they were all cut from one ``selected_loci.tsv``.
        """
        loci = list(self.train_dataset.loci)
        if self.locus_n_shards > 1:
            if not 0 <= self.locus_shard_index < self.locus_n_shards:
                raise ValueError(
                    f"locus_shard_index must satisfy 0 <= index < locus_n_shards, got "
                    f"{self.locus_shard_index} of {self.locus_n_shards}"
                )
            loci = [
                locus for i, locus in enumerate(loci)
                if i % self.locus_n_shards == self.locus_shard_index
            ]
        if self.locus_limit is not None:
            loci = loci[: int(self.locus_limit)]
        if not loci:
            raise RuntimeError(
                f"{self.client_id}: no loci selected. Shard "
                f"{self.locus_shard_index}/{self.locus_n_shards} of "
                f"{len(self.train_dataset.loci)} locus/loci is empty."
            )
        return loci

    def _instances_for(self, locus_id):
        """Instances belonging to one locus.

        Instance ids are ``<locus_id>_<architecture_id>_rep<k>``, which is the naming the
        simulation writes and the causal manifest keys on. Matching by prefix rather than
        by parsing keeps this module ignorant of how an architecture id is spelled.
        """
        prefix = f"{locus_id}_"
        found = [i for i in self.train_dataset.instances if i.startswith(prefix)]
        if self.instance_limit is not None:
            found = found[: int(self.instance_limit)]
        return found

    def _phenotypes(self, instances, iids):
        """Phenotypes for ``instances``, aligned to this site's fam order.

        A partially phenotyped block is refused rather than dropped: it would make ``n``
        differ between ``G`` and ``c``, and no single ``(G, c, u, n, q, w)`` tuple can
        express that (Theorem 1 condition (ii)).
        """
        Y = np.empty((len(iids), len(instances)), dtype=np.float64)
        for col, inst in enumerate(instances):
            path = self.train_dataset.phenotype_dir / f"{inst}.pheno"
            pheno = pd.read_csv(path, sep="\t", dtype={"FID": str, "IID": str})
            y = pheno.set_index("IID")["y"].reindex(iids.values).to_numpy(dtype=np.float64)
            if np.isnan(y).any():
                raise RuntimeError(
                    f"{self.client_id}: {path.name} does not cover "
                    f"{int(np.isnan(y).sum())} of its {len(iids)} individuals"
                )
            Y[:, col] = y
        return Y

    def _pack_gram(self, G):
        """Downcast the Gram for transmission when that is provably lossless."""
        if self.gram_dtype == "float32" and np.abs(G).max(initial=0.0) < _EXACT_FLOAT32_INT:
            return torch.from_numpy(G.astype(np.float32))
        return torch.from_numpy(G)

    # -- the exchange ------------------------------------------------------

    def train(self, **kwargs):
        # Single round. This is aggregate federated learning: there is one exchange, and
        # a second invocation would recompute identical results at full cost.
        if self._has_run:
            return

        dataset = self.train_dataset
        self.logger.info(f"{self.client_id}: reading site index from {dataset.data_dir}")

        bim = read_bim(dataset.plink_bim)
        fam = read_fam(dataset.plink_fam)
        manifest = pd.read_csv(
            dataset.manifest_path, sep="\t",
            dtype={"FID": str, "IID": str, "superpopulation": str},
        )
        reference = pd.read_csv(
            dataset.reference_path, sep="\t",
            dtype={"snp_id": str, "chrom": str, "bp": np.int64, "a1": str, "a2": str},
        )
        if reference["snp_id"].duplicated().any():
            raise RuntimeError(
                f"{self.client_id}: reference_variants.tsv has duplicate variant ids"
            )
        reference = reference.set_index("snp_id")[["chrom", "bp", "a1", "a2"]]

        # Harmonize once, chromosome-wide, not once per locus: the .bim is half a million
        # rows and re-checking it per locus would cost more than every Gram in the run.
        bim, flip_all = harmonize_bim(bim, reference, self.client_id)
        n_flipped = int(flip_all.sum())
        if n_flipped:
            self.logger.info(
                f"{self.client_id}: {n_flipped} of {len(bim)} variants recoded to the "
                "reference allele order"
            )

        pops = self.pops or list(dataset.composition)
        row_of = pd.Series(np.arange(len(fam)), index=fam["IID"].values)
        loci = self._selected_loci()

        payload = {}
        geno_meta = []
        pheno_meta = []
        all_instances = []

        for locus in loci:
            locus_id = str(locus["locus_id"])
            instances = self._instances_for(locus_id)
            if not instances:
                self.logger.warning(
                    f"{self.client_id}: locus {locus_id} has no phenotype instances; "
                    "emitting its genotype moments only"
                )
            all_instances.extend(instances)

            variants = locus_window_variants(
                bim, locus["chrom"], locus["start_bp"], locus["end_bp"]
            )
            if variants.empty:
                raise RuntimeError(
                    f"{self.client_id}: no variants in the {locus_id} window "
                    f"(chr{locus['chrom']}:{locus['start_bp']}-{locus['end_bp']})"
                )

            Y = (
                self._phenotypes(instances, fam["IID"])
                if instances
                else np.empty((len(fam), 0), dtype=np.float64)
            )

            X = read_bed_variants(dataset.plink_prefix, variants.index.to_numpy(), len(fam))
            # Recode the variants this site has the other way round, so that every site's
            # column j counts the same allele and the Grams are addable.
            flip = flip_all[variants.index.to_numpy()]
            if flip.any():
                X[:, flip] = 2.0 - X[:, flip]
            variants = variants[["chrom", "snp_id", "bp", "a1", "a2"]].reset_index(drop=True)

            try:
                for pop in pops:
                    members = manifest.loc[
                        manifest["superpopulation"] == pop, "IID"
                    ]
                    if members.empty:
                        continue  # this site holds none of this ancestry
                    rows = np.sort(
                        row_of.reindex(members.values).to_numpy(dtype=np.int64)
                    )
                    Xp = X[rows]
                    keep = ~columns_with_missing(Xp)
                    if not keep.all():
                        Xp = Xp[:, keep]
                    n = int(Xp.shape[0])
                    block_vars = (
                        variants if keep.all()
                        else variants.loc[keep].reset_index(drop=True)
                    )

                    payload[block_key(KEY_GRAM, locus_id, pop)] = self._pack_gram(gram(Xp))
                    payload[block_key(KEY_USUM, locus_id, pop)] = torch.from_numpy(
                        Xp.sum(axis=0, dtype=np.float64)
                    )
                    geno_meta.append({
                        "locus_id": locus_id,
                        "pop": pop,
                        "n": n,
                        "n_incomplete_variants": int((~keep).sum()),
                        "variants": {
                            "chrom": block_vars["chrom"].astype(str).tolist(),
                            "snp_id": block_vars["snp_id"].astype(str).tolist(),
                            "bp": block_vars["bp"].astype(np.int64).tolist(),
                            "a1": block_vars["a1"].astype(str).tolist(),
                            "a2": block_vars["a2"].astype(str).tolist(),
                        },
                    })

                    if instances:
                        Yp = Y[rows]
                        C = xty(Xp, Yp)
                        q = Yp.sum(axis=0)
                        w = np.einsum("ij,ij->j", Yp, Yp)
                        for col, inst in enumerate(instances):
                            payload[block_key(KEY_XTY, inst, pop)] = torch.from_numpy(
                                np.ascontiguousarray(C[:, col])
                            )
                            payload[block_key(KEY_SCALARS, inst, pop)] = torch.tensor(
                                [float(q[col]), float(w[col]), float(n)],
                                dtype=torch.float64,
                            )
                            pheno_meta.append({
                                "instance": inst,
                                "locus_id": locus_id,
                                "pop": pop,
                                "n": n,
                            })
                        del Yp
                    del Xp
            finally:
                del X

            self.logger.info(
                f"{self.client_id}: locus {locus_id} -- "
                f"{sum(1 for m in geno_meta if m['locus_id'] == locus_id)} ancestry "
                f"block(s), {len(instances)} instance(s)"
            )

        manifest_blob = json.dumps({
            "client_id": self.client_id,
            "sample_size": int(len(fam)),
            "composition": dataset.composition,
            "pops": pops,
            "loci": [
                {
                    "locus_id": str(locus["locus_id"]),
                    "chrom": str(locus["chrom"]),
                    "start_bp": int(locus["start_bp"]),
                    "end_bp": int(locus["end_bp"]),
                    "stratum": str(locus.get("stratum", "")),
                }
                for locus in loci
            ],
            "instances": all_instances,
            "geno_blocks": geno_meta,
            "pheno_blocks": pheno_meta,
            "n_flipped": n_flipped,
            "gram_dtype": self.gram_dtype,
            "locus_shard": [self.locus_shard_index, self.locus_n_shards],
        }).encode("utf-8")

        # THE ONLY THING THAT LEAVES THIS SITE.
        payload[KEY_MANIFEST] = torch.frombuffer(bytearray(manifest_blob), dtype=torch.uint8)
        self.model_state = payload

        uplink_mb = sum(
            t.numel() * t.element_size() for t in payload.values()
        ) / (1024 * 1024)
        self.logger.info(
            f"{self.client_id}: {len(geno_meta)} genotype block(s), "
            f"{len(pheno_meta)} phenotype block(s), {uplink_mb:.1f} MB uplink"
        )

        self._has_run = True
        self.round += 1
