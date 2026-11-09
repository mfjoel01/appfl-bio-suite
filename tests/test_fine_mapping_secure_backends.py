"""Secure and plaintext site payloads must fit identically with either backend."""

import json
import logging
import shutil

import numpy as np
import pandas as pd
import pytest


@pytest.mark.parametrize("method", ["susiex", "multisusie"])
def test_secure_backend_matches_plaintext(tmp_path, tmp_path_factory, monkeypatch, method):
    # This fixture includes allele flips and two ancestries with unequal site sizes.
    from test_fine_mapping_shipped_parity import POPS
    from test_fine_mapping_shipped_parity import package as make_package

    from appfl_bio_suite.experiments.fine_mapping.aggregator import FineMappingAggregator
    from appfl_bio_suite.experiments.fine_mapping.dataset import get_dataset
    from appfl_bio_suite.experiments.fine_mapping.trainer import SiteFineMappingTrainer

    susiex = shutil.which("SuSiEx")
    if method == "susiex" and not susiex:
        pytest.skip("SuSiEx binary required")
    if method == "multisusie":
        pytest.importorskip("MultiSuSiE")
    root, _, _ = make_package.__wrapped__(tmp_path_factory)
    sites = ["covenant", "mbzuai"]
    config = {"run_id": "test-" + method, "participants": sites}
    payloads = {False: {}, True: {}}
    for site in sites:
        dataset, _ = get_dataset(root / site / "data", site)
        for secure in (False, True):
            settings = {"trainer_output_dirname": str(tmp_path / site), "pops": POPS}
            if secure:
                path = tmp_path / f"{site}.json"
                path.write_text(
                    json.dumps(
                        {
                            **config,
                            "client_id": site,
                            "keys": {peer: "a1" * 32 for peer in sites if peer != site},
                        }
                    )
                )
                monkeypatch.setenv("APPFL_FINEMAP_SECAGG_KEYS", str(path))
                settings["secure_aggregation"] = config
            trainer = SiteFineMappingTrainer(
                train_dataset=dataset,
                train_configs=settings,
                client_id=site,
                logger=logging.getLogger("backend-parity"),
            )
            trainer.train()
            payloads[secure][site] = trainer.get_parameters()
    tables = []
    for secure in (False, True):
        output = tmp_path / str(secure)
        cfg = {
            "backend": method,
            "output_dir": str(output),
            "make_figures": False,
            "n_workers": 1,
            "pval_thresh": 1.0,
            "maf": 0.001,
            "n_signals": 2,
            "susiex_binary": susiex,
            "plink_binary": shutil.which("plink") or "plink",
        }
        if secure:
            cfg["secure_aggregation"] = config
        agg = FineMappingAggregator(
            aggregator_configs=cfg, logger=logging.getLogger("backend-parity")
        )
        agg.aggregate(payloads[secure])
        table = pd.read_csv(output / "data/fed_fm_results.tsv", sep="\t")
        assert len(table) == 2 and table.error.fillna("").eq("").all()
        assert table.converged.all()
        assert len(list((output / "data/instances").glob("*/variants.tsv"))) == 2
        if secure:
            assert (output / "data/secure_aggregation.json").is_file()
        tables.append(table)
    for plain_path in (tmp_path / "False/data/instances").glob("*/variants.tsv"):
        secure_path = tmp_path / "True/data/instances" / plain_path.parent.name / plain_path.name
        a = pd.read_csv(plain_path, sep="\t").set_index("snp_id").sort_index()
        b = pd.read_csv(secure_path, sep="\t").set_index("snp_id").sort_index()
        assert a.index.equals(b.index)
        np.testing.assert_allclose(a.pip, b.pip, atol=1e-5, rtol=1e-4)
    np.testing.assert_array_equal(tables[0].n_credible_sets, tables[1].n_credible_sets)
    np.testing.assert_allclose(tables[0].top_pip, tables[1].top_pip, atol=1e-5, rtol=1e-4)
