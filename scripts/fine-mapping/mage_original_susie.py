"""Re-fit the original MAGE panel with released FastQTL slopes and pinned susieR."""

import argparse
import json
import struct
import subprocess
import tempfile
import urllib.request
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import read_bed_variants
from appfl_bio_suite.experiments.fine_mapping.real_data.mage import (
    MAGE_DATA_ZIP_URL,
    _parse_zip_central_directory,
)
from appfl_bio_suite.experiments.fine_mapping.real_data.vcf_slice import (
    http_content_length,
    http_get_bytes,
)


def fetch_nominal(path):
    """Stream the two requested archive members while bounding memory usage."""
    url = MAGE_DATA_ZIP_URL
    members = _parse_zip_central_directory(url, http_content_length(url))
    path.parent.mkdir(parents=True, exist_ok=True)
    for dest in (path, Path(str(path) + ".tbi")):
        name = "eQTL_FastQTL_results.nominal_pass.allAssociations.MAGE.v1.0.txt.gz"
        if dest != path:
            name += ".tbi"
        (entry,) = [v for k, v in members.items() if Path(k).name == name]
        method, offset, size, expected = entry
        if dest.exists() and dest.stat().st_size == expected:
            continue
        header = http_get_bytes(url, offset, offset + 29)
        if header[:4] != b"PK\x03\x04" or method != 8:
            raise ValueError(f"Unexpected ZIP header for {name}")
        start = offset + 30 + sum(struct.unpack("<HH", header[26:30]))
        request = urllib.request.Request(
            url, headers={"Range": f"bytes={start}-{start + size - 1}"}
        )
        decoder = zlib.decompressobj(-15)
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as tmp:
            temp_path = Path(tmp.name)
            try:
                with urllib.request.urlopen(request, timeout=120) as response:
                    if response.status != 206:
                        raise ValueError(f"Expected partial response for {name}")
                    while chunk := response.read(4 << 20):
                        tmp.write(decoder.decompress(chunk))
                if not decoder.eof or tmp.tell() != expected:
                    raise ValueError(f"Incomplete ZIP member {name}")
                tmp.flush()
                temp_path.replace(dest)
            finally:
                temp_path.unlink(missing_ok=True)


def prepare(locus, base, meta, tabix, inputs):
    """Write aligned ALT dosages and the authors' z scores in VCF order."""
    prefix = base / "cache/plink/with_swap" / f"{locus.locus_id}_1000000"
    bim = pd.read_csv(f"{prefix}.bim", sep="\t", header=None)
    samples = pd.read_csv(f"{prefix}.fam", sep=r"\s+", header=None)[1].tolist()
    associations = {}
    for line in tabix.fetch(f"chr{locus.chrom}", locus.start_bp - 1, locus.end_bp + 1):
        fields = line.split("\t")
        if fields[6] == locus.ensembl_id:
            snp = fields[4]
            if snp in associations:
                raise ValueError(f"Duplicate nominal association for {snp}")
            associations[snp] = float(fields[13]) / float(fields[14])
    if not associations:
        raise ValueError(f"No nominal results for {locus.ensembl_id}")
    selected = bim[bim[1].isin(associations)]
    missing = associations.keys() - set(selected[1])
    if missing or len(selected) != len(associations) or len(selected) != locus.variants_tested:
        raise ValueError(
            f"{locus.locus_id}: nominal={len(associations)}, aligned={len(selected)}, "
            f"published={locus.variants_tested}, missing examples={sorted(missing)[:5]}"
        )
    dosage = read_bed_variants(prefix, selected.index.to_numpy(), len(samples))
    alt = selected[1].str.split(":").str[-1]
    flip = selected[4] != alt
    if not (selected.loc[flip, 5] == alt[flip]).all() or not np.isfinite(dosage).all():
        raise ValueError("Invalid dosage or allele orientation")
    dosage[:, flip] = 2 - dosage[:, flip]
    inputs.mkdir(parents=True, exist_ok=True)
    dosage.T.astype("uint8").tofile(inputs / "dosage.bin")
    pd.DataFrame({"snp_id": selected[1], "z": selected[1].map(associations)}).to_csv(
        inputs / "z.tsv", sep="\t", index=False
    )
    with np.load(meta / "expression_panel.npz", allow_pickle=True) as expr:
        row = expr["locus_ids"].tolist().index(locus.locus_id)
        y = pd.Series(expr["expr_raw"][row], index=expr["samples"])
    cov = pd.read_csv(meta / "covariates.tsv", sep="\t").set_index("sample_id").loc[samples]
    cov.insert(0, "expression", y.loc[samples].to_numpy())
    cov.to_csv(inputs / "samples.tsv", sep="\t")


def set_diagnostics(panel, meta, out):
    """Find best matches in both directions, retaining set labels, sizes and overlap."""
    published = pd.read_csv(meta / "published_cs_variants.tsv", sep="\t")
    approx = Path("local/output/mage-replication/arms/cen_susie_fwl_731_1mb")
    rows = []
    for locus in panel.itertuples():
        gene = published[published.ensemblID == locus.ensembl_id]
        reference = {k: set(v.variant_kgpID) for k, v in gene.groupby("variantCredibleSet")}
        for method, root in (("susieR_RSS", out), ("pooled_sufficient_statistics", approx)):
            (path,) = (root / "data/instances").glob(f"{locus.locus_id}_*/credible_sets.tsv")
            table = pd.read_csv(path, sep="\t")
            observed = {str(k): set(v.snp_id) for k, v in table.groupby("component")}
            for source, targets, candidates in (
                ("published", reference, observed),
                ("reported", observed, reference),
            ):
                for label, members in targets.items():
                    scores = {k: len(members & v) / len(members | v) for k, v in candidates.items()}
                    best = max(scores, key=scores.get) if scores else None
                    rows.append(
                        dict(
                            locus=locus.locus_id,
                            method=method,
                            source=source,
                            component=label,
                            size=len(members),
                            best_component=best,
                            best_jaccard=scores.get(best, 0.0),
                        )
                    )
    pd.DataFrame(rows).to_csv(out / "set_matches.tsv", sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--nominal", type=Path, required=True, help="Indexed MAGE all-associations BGZF"
    )
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--locus", help="Fit one manifest locus")
    parser.add_argument("--report", action="store_true")
    parser.add_argument("--download", action="store_true", help="Stage the nominal file and index")
    args = parser.parse_args()
    if args.download:
        fetch_nominal(args.nominal)
        return
    base = Path("local/data/mage")
    root = base / "federation_with_swap_fwl_1mb"
    meta = root / "meta"
    out = Path("local/output/mage-original-susie")
    panel = pd.read_csv(meta / "panel_manifest.tsv", sep="\t")
    if args.report:
        from appfl_bio_suite.experiments.fine_mapping.real_data.benchmark import (
            score_against_published,
        )

        instances = out / "data/instances"
        pub = pd.read_csv(meta / "published_cs_variants.tsv", sep="\t")
        checks = []
        for locus in panel.itertuples():
            inst = instances / f"{locus.locus_id}_mage"
            fit = pd.read_csv(inst / "fit.tsv", sep="\t").iloc[0]
            pips = pd.read_csv(inst / "variants.tsv", sep="\t").set_index("snp_id").pip
            ref = pub[pub.ensemblID == locus.ensembl_id].set_index("variant_kgpID")
            diff = pips.loc[ref.index] - ref.variantPIP
            checks.append(
                dict(
                    locus_id=locus.locus_id,
                    **fit.to_dict(),
                    max_abs_diff_published_pip=diff.abs().max(),
                )
            )
        set_diagnostics(panel, meta, out)
        summary, per_locus = score_against_published(out, root)
        per_locus = per_locus.merge(pd.DataFrame(checks), on="locus_id", validate="one_to_one")
        per_locus.to_csv(out / "replication_loci.tsv", sep="\t", index=False)
        (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        return
    if args.locus:
        if args.locus not in set(panel.locus_id):
            raise ValueError(f"Unknown locus {args.locus!r}")
        panel = panel[panel.locus_id == args.locus]
    with pysam.TabixFile(str(args.nominal)) as tabix:
        for locus in panel.itertuples():
            inputs = out / "inputs" / locus.locus_id
            dest = out / "data/instances" / f"{locus.locus_id}_mage"
            prepare(locus, base, meta, tabix, inputs)
            subprocess.run(
                [args.rscript, str(Path(__file__).with_suffix(".R")), str(inputs), str(dest)],
                check=True,
            )


if __name__ == "__main__":
    main()
