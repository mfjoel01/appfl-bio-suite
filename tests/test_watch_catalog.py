"""Publishing a catalogue merges partners without creating runnable endpoints."""

import base64
import json

import pytest
from click.testing import CliRunner

from appfl_bio_suite.core.config import load_federation
from appfl_bio_suite.core.experiments import repo_root
from appfl_bio_suite.core.watch import WatchError, export_site, write_network_run
from appfl_bio_suite.core.watch_catalog import PREVIEW_ROWS, load_catalog

pytest.importorskip("hivewatch")


@pytest.fixture
def catalog(tmp_path):
    table = tmp_path / "measurements.tsv"
    table.write_text("variant\tscore\n" + "".join(f"v{i}\t{i}\n" for i in range(205)))
    report = tmp_path / "report.html"
    report.write_text("<h1>Study report</h1>")
    raw = {
        "partners": [
            {
                "id": "partner-new",
                "name": "New Institute",
                "country": "Example country",
                "projects": ["fine-mapping"],
                "stage": "1 - Initial Planning",
                "contacts": [{"name": "Research contact", "email": "research@example.org"}],
                "location": {"lat": 0, "lng": 0, "city": "Example city"},
                "location_basis": "Representative headquarters.",
                "source_url": "https://example.org/location",
            }
        ],
        "results": [
            {
                "experiment": "fine-mapping",
                "title": "Reference study",
                "artifacts": [
                    {"title": "Scores", "path": "measurements.tsv"},
                    {"title": "Report", "path": "report.html"},
                ],
            }
        ],
    }
    path = tmp_path / "catalog.json"
    path.write_text(json.dumps(raw))
    return path, raw


def test_export_embeds_selected_results_and_preserves_full_download(catalog, tmp_path):
    path, _ = catalog
    federation = load_federation(repo_root() / "federation.yaml.example")
    original = federation.model_dump()
    out = tmp_path / "public"
    export_site(federation, out, catalog_path=path)
    assert federation.model_dump() == original
    assert {p.name for p in out.iterdir()} == {"index.html", "map.html", "network.map.json"}
    text = (out / "network.map.json").read_text()
    assert str(tmp_path) not in text
    metadata = json.loads(text)
    clients = metadata["rounds"][0]["clients"]
    new = next(c for c in clients if c["client_id"] == "partner-new")
    assert new["lat"] == new["lng"] == 0
    assert new.get("num_samples") is None
    assert new["status"] == "idle"
    assert new["contacts"][0]["email"] == "research@example.org"
    assert "endpoints" not in new
    table, report = metadata["results"][0]["artifacts"]
    assert table["total_rows"] == 205
    assert len(table["rows"]) == PREVIEW_ROWS
    assert (
        base64.b64decode(table["download"].split(",", 1)[1])
        == (path.parent / "measurements.tsv").read_bytes()
    )
    assert report["kind"] == "report"
    _, live = write_network_run(federation, tmp_path / "runs", catalog_path=path)
    assert json.loads(live.read_text())["results"] == metadata["results"]


def test_matching_partner_retains_federation_identity_and_samples(catalog, tmp_path):
    path, raw = catalog
    raw["partners"][0]["id"] = "site-north"
    path.write_text(json.dumps(raw))
    federation = load_federation(repo_root() / "federation.yaml.example")
    _, baseline = write_network_run(federation, tmp_path / "baseline")
    _, enriched = write_network_run(federation, tmp_path / "enriched", catalog_path=path)
    before = json.loads(baseline.read_text())["rounds"][0]["clients"]
    after = json.loads(enriched.read_text())["rounds"][0]["clients"]
    assert len(before) == len(after)
    old = next(c for c in before if c["client_id"] == "site-north")
    new = next(c for c in after if c["client_id"] == "site-north")
    for field in ("institution", "num_samples", "lat", "lng", "status", "endpoints"):
        assert old[field] == new[field]
    assert new["partnership_stage"] == "1 - Initial Planning"


@pytest.mark.parametrize(
    "change",
    [
        lambda raw: raw["partners"][0]["location"].update(lat=91),
        lambda raw: raw["partners"].append(raw["partners"][0]),
        lambda raw: raw["partners"][0].update(source_url="javascript:alert(1)"),
        lambda raw: raw["results"][0]["artifacts"][0].update(path="missing.tsv"),
    ],
)
def test_invalid_catalogue_fails_before_writing_export(catalog, tmp_path, change):
    path, raw = catalog
    change(raw)
    path.write_text(json.dumps(raw))
    federation = load_federation(repo_root() / "federation.yaml.example")
    with pytest.raises(WatchError, match="Cannot load watch catalogue"):
        export_site(federation, tmp_path / "public", catalog_path=path)
    assert not (tmp_path / "public" / "network.map.json").exists()


def test_catalogue_is_explicit_and_cli_can_publish_it(catalog, tmp_path):
    from appfl_bio_suite.cli import main

    assert load_catalog() == {"partners": [], "results": []}
    path, _ = catalog
    result = CliRunner().invoke(
        main,
        [
            "watch",
            "export",
            "--federation",
            str(repo_root() / "federation.yaml.example"),
            "--catalog",
            str(path),
            "--out",
            str(tmp_path / "cli-export"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "Reference study" in (tmp_path / "cli-export" / "network.map.json").read_text()
