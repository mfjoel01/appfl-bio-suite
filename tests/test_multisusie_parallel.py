"""Parallel checkpoint preparation must preserve the serial numerical result."""

from pathlib import Path

import numpy as np
import pytest
from test_fine_mapping_federated import _build_package


def test_parallel_locus_matches_serial_moments_and_posterior(tmp_path, monkeypatch):
    pytest.importorskip("MultiSuSiE")
    monkeypatch.syspath_prepend(str(Path(__file__).parents[1] / "scripts/fine-mapping"))
    import parallel_multisusie_locus as m

    config = _build_package(tmp_path / "package")
    serial = tmp_path / "serial"
    serial.mkdir()
    config = m.p.freeze_config(serial, "benchmark", "federation", config)
    task = m.p.make_tasks(
        serial, "benchmark", [dict(name="federation", config=str(config), central_check=False)]
    )[0]
    m.p.work_item(serial, task)
    parallel = tmp_path / "parallel"
    parallel.mkdir()
    m.complete_locus(parallel, task, 2)
    files = list((serial / "benchmark/federation/L0000").glob("*/multisusie/fit.json"))
    assert len(files) == 2
    for original in files:
        result = parallel / original.relative_to(serial)
        assert m.p.read(original)["moments_sha256"] == m.p.read(result)["moments_sha256"]
        with (
            np.load(original.parent / "posterior.npz") as a,
            np.load(result.parent / "posterior.npz") as b,
        ):
            assert set(a.files) == set(b.files)
            for key in a.files:
                np.testing.assert_array_equal(a[key], b[key])
    assert (parallel / "benchmark/federation/L0000/accepted.json").exists()
