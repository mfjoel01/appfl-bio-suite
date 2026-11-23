"""Run the prespecified MAGE parity comparison from the repository root."""

import argparse
from pathlib import Path

import pandas as pd

from appfl_bio_suite.experiments.fine_mapping.real_data.benchmark import (
    compare_arm_pairs,
    run_centralized_arm,
    run_federated_arm,
)
from appfl_bio_suite.experiments.fine_mapping.real_data.federation_builder import (
    build_mage_federation,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("build", "fit", "all"), default="all")
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    panel = pd.read_csv("docs/experiments/fine-mapping/mage/heldout_panel.tsv", sep="\t")
    root = Path("local/data/mage/federation_heldout")
    out = Path("local/output/mage-heldout")
    if args.stage in ("build", "all"):
        build_mage_federation(root, phenotype="raw", panel_gene_ids=panel.ensembl_id.tolist())
    if args.stage == "build":
        return
    built = pd.read_csv(root / "meta/panel_manifest.tsv", sep="\t")
    if built.ensembl_id.tolist() != panel.ensembl_id.tolist():
        raise ValueError("Federation differs from the prespecified panel")
    loci = built[["locus_id", "gene_symbol"]]
    tables = []
    for method in ("susie", "multisusie"):
        settings = dict(
            backend="multisusie",
            single_population=method == "susie",
            tol=1e-5,
            n_workers=args.workers,
        )
        pooled, fed = out / f"pooled_{method}", out / f"federated_{method}"
        run_centralized_arm(root, pooled, pooled.name, fwl=True, **settings)
        run_federated_arm(root, fed, fed.name, covariate_adjustment="fwl", **settings)
        parity, _ = compare_arm_pairs(pooled, fed, method)
        table = loci.merge(parity, on="locus_id", how="left", validate="one_to_one")
        table["comparison"] = method
        for label, directory in (("pooled", pooled), ("federated", fed)):
            fits = pd.read_csv(directory / "data/fed_fm_results.tsv", sep="\t")
            status = fits.set_index("locus_id")["fit_status"]
            table[f"{label}_status"] = table.locus_id.map(status)
        tables.append(table)
        pd.concat(tables).to_csv(out / "heldout_parity.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
