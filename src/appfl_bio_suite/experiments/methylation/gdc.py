"""Download the prespecified TCGA-LAML FAB cohort from open GDC data.

Run with ``python -m appfl_bio_suite.experiments.methylation.gdc --out DIR``.
Metadata and verified raw downloads are retained for reproducibility and restart.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np

from .dataset import MethylationData

API = "https://api.gdc.cancer.gov"
FAB_CLASSES = ("M1", "M2", "M4", "M5")
PRIMARY = "Primary Blood Derived Cancer - Peripheral Blood"


def query(endpoint: str, filters: dict, fields: str) -> dict:
    url = f"{API}/{endpoint}?" + urlencode(
        {"filters": json.dumps(filters), "fields": fields, "size": 10000}
    )
    with urlopen(url, timeout=120) as response:
        result = json.load(response)
    if len(result["data"]["hits"]) != result["data"]["pagination"]["total"]:
        raise ValueError("GDC query was truncated")
    return result


def select_cohort(files: list[dict], cases: list[dict]) -> tuple[list[dict], dict]:
    """Join exact case UUIDs; reject ambiguous patients, samples or labels."""
    by_id = {c["case_id"]: c for c in cases}
    if len(by_id) != len(cases):
        raise ValueError("Duplicate clinical cases")
    rows, seen, counts = [], set(), Counter()
    for file in files:
        if (file["access"], file["platform"], file["data_type"]) != (
            "open",
            "Illumina Human Methylation 450",
            "Methylation Beta Value",
        ):
            raise ValueError("Unexpected file access, platform or type")
        if len(file["cases"]) != 1:
            raise ValueError("Ambiguous file case")
        case = file["cases"][0]
        samples = case["samples"]
        if len(samples) != 1 or samples[0]["sample_type"] != PRIMARY:
            raise ValueError("Expected one primary blood cancer sample")
        clinical = by_id[case["case_id"]]
        if clinical["submitter_id"] != case["submitter_id"]:
            raise ValueError("Patient identifier mismatch")
        if case["case_id"] in seen:
            raise ValueError("Multiple beta files for a patient")
        seen.add(case["case_id"])
        labels = {d.get("fab_morphology_code") for d in clinical.get("diagnoses", [])}
        if len(labels) != 1 or None in labels:
            raise ValueError("Missing or conflicting FAB label")
        label = labels.pop()
        counts[label] += 1
        if label in FAB_CLASSES:
            rows.append(
                {
                    "file": file,
                    "patient": case["submitter_id"],
                    "sample": samples[0]["submitter_id"],
                    "label": label,
                }
            )
    if any(counts[label] < 20 for label in FAB_CLASSES):
        raise ValueError("Prespecified four-class cohort no longer meets 20-patient gate")
    return sorted(rows, key=lambda row: row["patient"]), dict(sorted(counts.items()))


def digest(path: Path, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download(file: dict, cache: Path) -> Path:
    file_id = file["file_id"]
    if not re.fullmatch(r"[0-9a-f-]{36}", file_id):
        raise ValueError("Invalid GDC file UUID")
    path = cache / f"{file_id}.txt"

    def verified(candidate: Path) -> bool:
        return (
            candidate.is_file()
            and candidate.stat().st_size == file["file_size"]
            and digest(candidate, "md5") == file["md5sum"]
        )

    if verified(path):
        return path
    temporary = path.with_suffix(".part")
    for attempt in range(3):
        try:
            with urlopen(f"{API}/data/{file_id}", timeout=120) as response:
                with temporary.open("wb") as stream:
                    shutil.copyfileobj(response, stream)
            if not verified(temporary):
                raise ValueError("GDC download checksum/size mismatch")
            temporary.replace(path)
            return path
        except (OSError, ValueError):
            temporary.unlink(missing_ok=True)
            if attempt == 2:
                raise
            time.sleep(attempt + 1)
    raise RuntimeError("Unreachable")


def read_beta(path: Path) -> tuple[np.ndarray, np.ndarray]:
    probes, values = [], []
    with path.open() as stream:
        for line in stream:
            probe, value = line.rstrip("\n").split("\t")
            if not re.fullmatch(r"cg\d+", probe):
                continue
            beta = float("nan") if value in ("NA", "NaN", "nan") else float(value)
            probes.append(probe)
            values.append(beta)
    if not probes or len(set(probes)) != len(probes):
        raise ValueError("Empty or duplicate CpG probes")
    beta = np.asarray(values, dtype=np.float64)
    if np.isinf(beta).any() or ((beta < 0) | (beta > 1)).any():
        raise ValueError("Invalid beta value")
    return np.asarray(probes), beta.astype(np.float32)


def build(out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    target = out / "cohort.npz"
    if target.exists():
        raise FileExistsError(target)
    cache = out / "raw"
    cache.mkdir(exist_ok=True)
    metadata = out / "metadata.json"
    if metadata.exists():
        record = json.loads(metadata.read_text())
    else:

        def filters(pairs):
            return {
                "op": "and",
                "content": [{"op": "in", "content": {"field": k, "value": [v]}} for k, v in pairs],
            }

        files = query(
            "files",
            filters(
                [
                    ("cases.project.project_id", "TCGA-LAML"),
                    ("access", "open"),
                    ("platform", "Illumina Human Methylation 450"),
                    ("data_type", "Methylation Beta Value"),
                ]
            ),
            "file_id,file_name,file_size,md5sum,access,platform,data_type,"
            "cases.case_id,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type",
        )
        cases = query(
            "cases",
            filters([("project.project_id", "TCGA-LAML")]),
            "case_id,submitter_id,diagnoses.fab_morphology_code",
        )
        record = {"retrieved_utc": datetime.now(UTC).isoformat(), "files": files, "cases": cases}
        temporary = metadata.with_suffix(".tmp")
        temporary.write_text(json.dumps(record, indent=2) + "\n")
        temporary.replace(metadata)
    rows, counts = select_cohort(record["files"]["data"]["hits"], record["cases"]["data"]["hits"])
    print(f"Selected {len(rows)} primary patients; FAB counts: {counts}", flush=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        paths = list(pool.map(lambda row: download(row["file"], cache), rows))
    first_probes, first_values = read_beta(paths[0])
    matrix = np.empty((len(rows), len(first_probes)), dtype=np.float32)
    matrix[0] = first_values
    for index, path in enumerate(paths[1:], 1):
        probes, values = read_beta(path)
        if not np.array_equal(probes, first_probes):
            raise ValueError("CpG schema/order differs across files")
        matrix[index] = values
    data = MethylationData(
        matrix,
        np.array([r["label"] for r in rows]),
        first_probes,
        patient_ids=np.array([r["patient"] for r in rows]),
        sample_ids=np.array([r["sample"] for r in rows]),
    )
    temporary = out / "cohort.partial.npz"
    temporary.unlink(missing_ok=True)
    data.save(temporary)
    manifest = {
        "project": "TCGA-LAML",
        "label": "FAB morphology",
        "selected_classes": FAB_CLASSES,
        "matched_class_counts": counts,
        "samples": len(rows),
        "cpgs": len(first_probes),
        "missing_fraction": float(np.isnan(matrix).mean()),
        "dataset_sha256": digest(temporary),
        "metadata_sha256": digest(metadata),
        "loader_sha256": digest(Path(__file__)),
        "retrieved_utc": record["retrieved_utc"],
        "rows": rows,
        "probe_policy": "All cg probes in source order; no variance filtering or imputation",
        "terms": "https://gdc.cancer.gov/analyze-data/data-analysis-policies",
    }
    manifest_tmp = out / "manifest.tmp"
    manifest_tmp.write_text(json.dumps(manifest, indent=2) + "\n")
    manifest_tmp.replace(out / "manifest.json")
    temporary.replace(target)
    print(json.dumps({k: v for k, v in manifest.items() if k != "rows"}), flush=True)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    build(parser.parse_args().out)
