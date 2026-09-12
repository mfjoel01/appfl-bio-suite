"""Explicit, coordinator-owned partner and result catalogue for the network viewer.

The catalogue is presentation data. Adding an institution here never creates a compute
endpoint or enables a runnable experiment. Result files are selected explicitly; no
output directories are crawled, and their filesystem paths are never published.
"""

from __future__ import annotations

import base64
import csv
import io
import json
import math
from pathlib import Path
from urllib.parse import urlsplit

MAX_ASSET_BYTES = 16 * 1024 * 1024
PREVIEW_ROWS = 200
_MIMES = {
    ".png": "image/png",
    ".jpg": "image/jpeg",
    ".jpeg": "image/jpeg",
    ".webp": "image/webp",
    ".html": "text/html",
    ".csv": "text/csv",
    ".tsv": "text/tab-separated-values",
}


def _text(value: object, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"catalogue {field} must be a nonempty string")
    return value.strip()


def _url(value: object) -> str:
    url = _text(value, "source_url")
    parsed = urlsplit(url)
    if parsed.scheme not in {"http", "https"} or not parsed.netloc:
        raise ValueError("catalogue source_url must be an http(s) URL")
    return url


def _partner(raw: dict) -> dict:
    partner = {
        "client_id": _text(raw.get("id"), "partner id"),
        "institution": _text(raw.get("name"), "partner name"),
        "country": _text(raw.get("country"), "partner country"),
        "status": "idle",
        "num_samples": None,
        "lat": None,
        "lng": None,
        "partnership_stage": _text(raw.get("stage"), "partner stage"),
    }
    projects = raw.get("projects", [])
    if not isinstance(projects, list):
        raise ValueError("catalogue partner projects must be a list")
    partner["experiments"] = ", ".join(dict.fromkeys(_text(p, "project") for p in projects))
    for key in ("notes", "location_basis"):
        if raw.get(key):
            partner[key] = _text(raw[key], key)
    if raw.get("source_url"):
        partner["location_source"] = _url(raw["source_url"])
    location = raw.get("location")
    if location is not None:
        if not isinstance(location, dict):
            raise ValueError("catalogue location must be an object or null")
        for key, bound in (("lat", 90), ("lng", 180)):
            value = location.get(key)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or abs(value) > bound
            ):
                raise ValueError(f"catalogue partner {partner['client_id']} has invalid {key}")
            partner[key] = value
        partner["city"] = _text(location.get("city"), "city")
    contacts = raw.get("contacts", [])
    if not isinstance(contacts, list):
        raise ValueError("catalogue contacts must be a list")
    partner["contacts"] = []
    for contact in contacts:
        email = _text(contact.get("email"), "contact email")
        if "@" not in email or any(c in email for c in '\r\n<>" '):
            raise ValueError("catalogue contact email is invalid")
        partner["contacts"].append(
            {"name": _text(contact.get("name"), "contact name"), "email": email}
        )
    return partner


def _artifact(raw: dict, directory: Path) -> dict:
    path = Path(_text(raw.get("path"), "artifact path"))
    if not path.is_absolute():
        path = directory / path
    suffix = path.suffix.lower()
    if suffix not in _MIMES:
        raise ValueError(f"unsupported result file type: {suffix}")
    if path.stat().st_size > MAX_ASSET_BYTES:
        raise ValueError(f"result file exceeds the 16 MiB limit: {path.name}")
    content = path.read_bytes()
    result = {
        "title": _text(raw.get("title"), "artifact title"),
        "description": str(raw.get("description", "")),
        "filename": path.name,
        "download": f"data:{_MIMES[suffix]};base64,{base64.b64encode(content).decode('ascii')}",
    }
    if suffix in {".csv", ".tsv"}:
        reader = csv.reader(
            io.StringIO(content.decode("utf-8-sig")), delimiter="\t" if suffix == ".tsv" else ","
        )
        columns = next(reader, [])
        if not columns:
            raise ValueError(f"result table has no column header: {path.name}")
        rows = [row for row in reader if row]
        if any(len(row) != len(columns) for row in rows):
            raise ValueError(f"result table has inconsistent columns: {path.name}")
        result.update(kind="table", columns=columns, rows=rows[:PREVIEW_ROWS], total_rows=len(rows))
    else:
        result["kind"] = "report" if suffix == ".html" else "image"
    return result


def load_catalog(path: Path | str | None = None) -> dict:
    """Compile an optional JSON catalogue, embedding selected result artifacts."""
    if path is None:
        return {"partners": [], "results": []}
    path = Path(path)
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(raw, dict):
            raise ValueError("catalogue root must be an object")
        partners = [_partner(p) for p in raw.get("partners", [])]
        ids = [p["client_id"] for p in partners]
        if len(set(ids)) != len(ids):
            raise ValueError("catalogue partner ids must be unique")
        results = []
        for group in raw.get("results", []):
            artifacts = [_artifact(a, path.parent) for a in group.get("artifacts", [])]
            results.append(
                {
                    "experiment": _text(group.get("experiment"), "result experiment"),
                    "title": _text(group.get("title"), "result title"),
                    "description": str(group.get("description", "")),
                    "artifacts": artifacts,
                }
            )
        return {"partners": partners, "results": results}
    except (OSError, ValueError, TypeError, AttributeError) as exc:
        raise ValueError(f"Cannot load watch catalogue {path}: {exc}") from exc


def add_partners(events: list[dict], catalog: dict, experiment: str | None = None) -> None:
    """Merge by institution id, retaining the federation's samples and endpoint state."""
    snapshot = next(event for event in events if event["event_type"] == "round_end")
    clients = {c["client_id"]: c for c in snapshot["clients"]}
    for partner in catalog["partners"]:
        if experiment and experiment not in partner["experiments"].split(", "):
            continue
        identity = partner["client_id"]
        if identity in clients:
            current = clients[identity]
            projects = list(
                dict.fromkeys((current["experiments"] + ", " + partner["experiments"]).split(", "))
            )
            # Coordinates and samples declared by a participating site are authoritative.
            current.update(
                {
                    k: v
                    for k, v in partner.items()
                    if k
                    not in {
                        "status",
                        "num_samples",
                        "lat",
                        "lng",
                        "city",
                        "experiments",
                        "institution",
                        "country",
                        "location_basis",
                        "location_source",
                    }
                }
            )
            current["experiments"] = ", ".join(p for p in projects if p)
        else:
            clients[identity] = dict(partner)
    snapshot["clients"] = list(clients.values())
    snapshot["round_metrics"]["num_selected"] = len(clients)
    init = next(event for event in events if event["event_type"] == "init")
    init["config"]["sites"] = len(clients)
    init["config"]["experiments"] = ", ".join(
        sorted(
            {
                project
                for client in clients.values()
                for project in client.get("experiments", "").split(", ")
                if project
            }
        )
    )
