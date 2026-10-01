"""Offline checks for source identity and beta parsing at the GDC boundary."""

import copy

import numpy as np
import pytest

from appfl_bio_suite.experiments.methylation.gdc import PRIMARY, read_beta, select_cohort


def cohort():
    files, cases = [], []
    for label in ("M1", "M2", "M4", "M5"):
        for i in range(20):
            key = f"{label}-{i}"
            case = {
                "case_id": key,
                "submitter_id": key,
                "diagnoses": [{"fab_morphology_code": label}],
            }
            cases.append(case)
            files.append(
                {
                    "access": "open",
                    "platform": "Illumina Human Methylation 450",
                    "data_type": "Methylation Beta Value",
                    "cases": [
                        {
                            "case_id": key,
                            "submitter_id": key,
                            "samples": [
                                {
                                    "submitter_id": key + "-sample",
                                    "sample_type": PRIMARY,
                                }
                            ],
                        }
                    ],
                }
            )
    return files, cases


def test_exact_join_and_gate():
    files, cases = cohort()
    rows, counts = select_cohort(files, list(reversed(cases)))
    assert len(rows) == 80
    assert counts == dict.fromkeys(("M1", "M2", "M4", "M5"), 20)
    with pytest.raises(ValueError, match="20-patient"):
        select_cohort(files[1:], cases)
    with pytest.raises(ValueError, match="Multiple beta"):
        select_cohort(files + [files[0]], cases)


@pytest.mark.parametrize("change", ["label", "identity", "normal", "controlled"])
def test_reject_ambiguous_or_ineligible_data(change):
    files, cases = copy.deepcopy(cohort())
    if change == "label":
        cases[0]["diagnoses"].append({"fab_morphology_code": "M0"})
    elif change == "identity":
        cases[0]["submitter_id"] = "wrong"
    elif change == "normal":
        files[0]["cases"][0]["samples"][0]["sample_type"] = "Blood Derived Normal"
    else:
        files[0]["access"] = "controlled"
    with pytest.raises(ValueError):
        select_cohort(files, cases)


def test_beta_preserves_missingness_and_order(tmp_path):
    path = tmp_path / "beta.txt"
    path.write_text("cg002\t0.25\ncg001\tNA\nrs123\t0.4\n")
    probes, beta = read_beta(path)
    assert probes.tolist() == ["cg002", "cg001"]
    assert beta[0] == 0.25 and np.isnan(beta[1])
    for contents in ("cg001\t1.1\n", "cg001\tinf\n", "cg001\t0\ncg001\t1\n"):
        path.write_text(contents)
        with pytest.raises(ValueError):
            read_beta(path)


def test_download_checks_cache_and_rejects_corrupt_response(tmp_path, monkeypatch):
    import hashlib
    import io

    from appfl_bio_suite.experiments.methylation import gdc

    payload = b"cg001\t0.5\n"
    file = {
        "file_id": "12345678-1234-1234-1234-123456789012",
        "file_size": len(payload),
        "md5sum": hashlib.md5(payload).hexdigest(),
    }
    calls = []

    def response(*args, **kwargs):
        calls.append(args)
        return io.BytesIO(payload)

    monkeypatch.setattr(gdc, "urlopen", response)
    path = gdc.download(file, tmp_path)
    assert path.read_bytes() == payload
    assert gdc.download(file, tmp_path) == path
    assert len(calls) == 1
    path.write_bytes(b"corrupt")
    assert gdc.download(file, tmp_path).read_bytes() == payload
    path.unlink()
    monkeypatch.setattr(gdc, "urlopen", lambda *a, **k: io.BytesIO(b"bad"))
    monkeypatch.setattr(gdc.time, "sleep", lambda _: None)
    with pytest.raises(ValueError, match="checksum"):
        gdc.download(file, tmp_path)
    assert not path.exists()
    assert not path.with_suffix(".part").exists()


def test_build_writes_pipeline_contract_and_manifest(tmp_path, monkeypatch):
    import json

    from appfl_bio_suite.experiments.methylation import gdc
    from appfl_bio_suite.experiments.methylation.dataset import MethylationData

    files, cases = cohort()
    (tmp_path / "metadata.json").write_text(
        json.dumps(
            {
                "retrieved_utc": "fixture",
                "files": {"data": {"hits": files}},
                "cases": {"data": {"hits": cases}},
            }
        )
    )
    beta = tmp_path / "fixture.txt"
    beta.write_text("cg001\t0.2\ncg002\tNA\n")
    monkeypatch.setattr(gdc, "download", lambda file, cache: beta)
    manifest = gdc.build(tmp_path)
    data = MethylationData.load(tmp_path / "cohort.npz")
    assert data.X.shape == (80, 2)
    assert len(set(data.patient_ids)) == 80
    assert manifest["missing_fraction"] == 0.5
    assert manifest["dataset_sha256"] == gdc.digest(tmp_path / "cohort.npz")
    assert json.loads((tmp_path / "manifest.json").read_text())["samples"] == 80
    with pytest.raises(FileExistsError):
        gdc.build(tmp_path)
