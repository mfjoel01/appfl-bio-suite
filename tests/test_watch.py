"""The federation network map: what it draws, and what it refuses to publish.

Three things are checked, and the middle one is the reason this file is longer than the
feature warrants:

1. **The map describes the federation the config describes.** One marker per institution,
   its experiments, its data, its declared samples -- and disabled experiments absent.

2. **Nothing that describes the inside of a partner's cluster reaches the page.** The
   export is a file a coordinator hands out. Service accounts, worker-side paths, the
   coordinator's own identity and (by default) full endpoint UUIDs must not be in it.
   That is asserted against the serialized bytes, not against the code's intent, because
   the failure mode is a field somebody adds later without thinking about where it goes.

3. **Watching cannot break a run.** A federated run costs a scheduler queue wait at every
   participating institution. Every path where monitoring could raise is exercised with
   something that raises.
"""

from __future__ import annotations

import json

import pytest

from appfl_bio_suite.core.config import FederationError, load_federation
from appfl_bio_suite.core.experiments import REGISTRY, repo_root
from appfl_bio_suite.core.watch import (
    NETWORK_RUN_ID,
    RunWatcher,
    WatchError,
    build_network_events,
    export_site,
    participations,
    unplaced_report,
    unplaced_sites,
    watchable,
    write_network_run,
    write_sidecar,
)

EXAMPLE = repo_root() / "federation.yaml.example"

# hivewatch assembles the .map.json and ships the viewer. It is in [all], so CI has it;
# a partial install skips the tests that need it rather than failing the suite.
hivewatch = pytest.importorskip("hivewatch", reason="the [watch] extra is not installed")


@pytest.fixture
def federation():
    return load_federation(EXAMPLE)


def _clients(federation, **kwargs) -> dict[str, dict]:
    events = build_network_events(federation, **kwargs)
    (round_end,) = [e for e in events if e["event_type"] == "round_end"]
    return {c["client_id"]: c for c in round_end["clients"]}


# ---------------------------------------------------------------------------
# what it draws
# ---------------------------------------------------------------------------


def test_one_marker_per_institution_not_per_participation(federation):
    """A site in three experiments is one building, and gets one dot.

    hivewatch keys markers by client_id, so an entry per participation would stack three
    markers at identical coordinates and report the federation as three times its size.
    """
    clients = _clients(federation)
    assert set(clients) == {s.id for s in federation.sites}
    assert clients["site-north"]["endpoints"] == 3


def test_each_site_carries_its_experiments_data_and_samples(federation):
    north = _clients(federation)["site-north"]
    assert north["institution"] == "Northern Institute of Technology"
    assert north["country"] == "Canada"
    assert (north["lat"], north["lng"]) == (43.6532, -79.3832)

    # One chip per experiment, labelled by the experiment, reading "<n> samples · <kind>".
    assert "18,032 samples" in north["gwas"]
    assert REGISTRY["gwas"].data_kind in north["gwas"]
    assert "199 samples" in north["flamby-heart-disease"]


def test_declared_sample_counts_are_summed_into_the_popup_field(federation):
    """`num_samples` is hivewatch's own field; the per-experiment chips break it down."""
    north = _clients(federation)["site-north"]
    assert north["num_samples"] == 199 + 18032 + 50000


def test_disabled_experiments_are_not_drawn(federation):
    """`enabled: false` means the federation is not running it. Neither is the map."""
    federation.experiments["gwas"].enabled = False
    north = _clients(federation)["site-north"]
    assert "gwas" not in north
    assert "gwas" not in north["experiments"]
    assert north["endpoints"] == 2


def test_experiment_filter_narrows_to_one(federation):
    north = _clients(federation, experiment="gwas")["site-north"]
    assert north["experiments"] == "gwas"
    assert north["num_samples"] == 18032


def test_the_coordinator_is_the_hub(federation):
    events = build_network_events(federation)
    (server,) = [e for e in events if e["event_type"] == "server_metadata"]
    assert server["server"]["lat"] == 41.8781
    assert server["server"]["org"] == federation.coordinator.organization


def test_probed_status_reaches_the_marker(federation):
    clients = _clients(federation, statuses={"site-north": "failed", "site-east": "idle"})
    assert clients["site-north"]["status"] == "failed"
    assert clients["site-east"]["status"] == "idle"
    # Unprobed sites are `active`: membership, which is what an unprobed map states.
    assert clients["site-south"]["status"] == "active"


def test_a_federation_with_no_participating_sites_is_refused(federation):
    for exp in federation.experiments.values():
        exp.sites = []
    with pytest.raises(WatchError, match="nothing to draw"):
        build_network_events(federation)


# ---------------------------------------------------------------------------
# what it refuses to publish
# ---------------------------------------------------------------------------


def test_export_publishes_no_partner_internals(federation, tmp_path):
    """Asserted against the bytes, not the intent.

    The way this breaks is somebody adding a useful-looking field to the client dict
    without asking where that dict ends up. Reading the file catches that; reading the
    code does not.
    """
    export_site(federation, tmp_path)
    published = (tmp_path / "network.map.json").read_text(encoding="utf-8")

    forbidden = {
        "service account": "gwas_svc",
        "another service account": "flamby_svc",
        "worker-side data path": "/scratch/gwas_svc/data/site2",
        "worker-side output path": "/work/finemap_svc/outputs",
        "coordinator identity": federation.coordinator.identity,
        "coordinator identity UUID": federation.coordinator.identity_id,
    }
    for what, value in forbidden.items():
        assert value not in published, f"{what} ({value!r}) reached the published map"


def test_endpoint_uuids_are_fingerprinted_by_default(federation, tmp_path):
    full = federation.experiments["gwas"].sites[0].endpoint_uuid
    export_site(federation, tmp_path)
    published = (tmp_path / "network.map.json").read_text(encoding="utf-8")

    assert full not in published
    assert full[:8] in published  # the fingerprint, which is the point of having one


def test_fingerprints_are_hidden_from_the_viewer(federation):
    """A leading underscore is hivewatch's contract for 'store it, do not render it'.

    Keyed off the value rather than the name: `endpoints: 3` is a legitimate visible
    chip, and a test that matched on "endpoint" in the key would call it a leak.
    """
    north = _clients(federation)["site-north"]
    fingerprints = [k for k, v in north.items() if isinstance(v, str) and "\u2026" in v]
    assert fingerprints, "no endpoint fingerprints were emitted at all"
    assert all(k.startswith("_") for k in fingerprints)


def test_full_uuids_only_when_explicitly_opted_into(federation, tmp_path):
    full = federation.experiments["gwas"].sites[0].endpoint_uuid
    export_site(federation, tmp_path, include_endpoint_uuids=True)
    assert full in (tmp_path / "network.map.json").read_text(encoding="utf-8")


def test_only_allowlisted_numeric_metrics_leave_a_run():
    """Everything else a trainer returns is a result, and results are not telemetry."""
    forwarded = watchable(
        {
            "local_accuracy": 0.91,
            "num_samples": 400,
            "beta": [0.1, 0.2, 0.3],
            "output_path": "/scratch/gwas_svc/outputs/site2",
            "converged": True,
            "site_id": "Site2",
        }
    )
    assert forwarded == {"local_accuracy": 0.91, "num_samples": 400}


def test_watchable_survives_something_that_is_not_a_mapping():
    assert watchable(None) == {}
    assert watchable(["not", "a", "mapping"]) == {}


# ---------------------------------------------------------------------------
# coordinates
# ---------------------------------------------------------------------------


def test_a_site_without_coordinates_is_reported_not_guessed(federation):
    federation.site("site-east").location = None

    assert [s.id for s in unplaced_sites(federation)] == ["site-east"]
    report = unplaced_report(federation)
    assert "site-east" in report
    assert "location:" in report  # the YAML to paste

    # It is still listed as a participant -- it just has nowhere to be drawn.
    east = _clients(federation)["site-east"]
    assert east["lat"] is None and east["lng"] is None
    assert east["institution"] == "Eastern Technical University"


def test_a_fully_placed_federation_reports_nothing(federation):
    assert unplaced_sites(federation) == []
    assert unplaced_report(federation) == ""


@pytest.mark.parametrize(
    "lat,lng",
    [(91.0, 0.0), (-91.0, 0.0), (0.0, 181.0), (0.0, -181.0)],
)
def test_out_of_range_coordinates_are_refused(tmp_path, lat, lng):
    """Caught here, in a second, rather than as a marker in the wrong hemisphere."""
    text = EXAMPLE.read_text(encoding="utf-8").replace(
        "      lat: 43.6532\n      lng: -79.3832\n", f"      lat: {lat}\n      lng: {lng}\n"
    )
    path = tmp_path / "federation.yaml"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(FederationError, match="out of range"):
        load_federation(path)


# ---------------------------------------------------------------------------
# artifacts
# ---------------------------------------------------------------------------


def test_the_network_view_is_a_hivewatch_run_in_both_its_formats(federation, tmp_path):
    """Both files, because hivewatch lists runs by globbing .jsonl and serves .map.json.

    A network view with only the metadata would load if you asked for it by name and
    would never appear in the run list.
    """
    jsonl, mapjson = write_network_run(federation, tmp_path)
    assert jsonl.name == f"{NETWORK_RUN_ID}.jsonl"

    header = json.loads(jsonl.read_text(encoding="utf-8").splitlines()[0])
    assert header["event_type"] == "init"
    assert header["run_id"] == NETWORK_RUN_ID

    metadata = json.loads(mapjson.read_text(encoding="utf-8"))
    assert metadata["schema_version"] == 1
    assert len(metadata["rounds"]) == 1
    assert len(metadata["rounds"][0]["clients"]) == 3


def test_rebuilding_replaces_rather_than_accumulates(federation, tmp_path):
    """Otherwise a renamed site haunts the run list forever."""
    write_network_run(federation, tmp_path)
    federation.sites = federation.sites[:1]
    for exp in federation.experiments.values():
        exp.sites = [e for e in exp.sites if e.site == "site-north"]

    jsonl, mapjson = write_network_run(federation, tmp_path)
    assert len(jsonl.read_text(encoding="utf-8").splitlines()) == 4
    metadata = json.loads(mapjson.read_text(encoding="utf-8"))
    assert len(metadata["rounds"][0]["clients"]) == 1


def test_export_is_three_self_contained_files(federation, tmp_path):
    export_site(federation, tmp_path)
    assert {p.name for p in tmp_path.iterdir()} == {"index.html", "map.html", "network.map.json"}

    index = (tmp_path / "index.html").read_text(encoding="utf-8")
    assert "map.html?metadata_url=network.map.json" in index

    # The viewer is hivewatch's, plus the one documented splice, and it is still HTML.
    page = (tmp_path / "map.html").read_text(encoding="utf-8")
    assert page.count("</body>") == 1
    assert page.index("appfl-bio-suite: render") < page.index("</body>")


# ---------------------------------------------------------------------------
# the run bridge
# ---------------------------------------------------------------------------


def test_the_sidecar_is_keyed_by_both_names_a_driver_might_use(federation, tmp_path):
    """The Globus Compute driver keys results by endpoint; the serial driver by client id."""
    path = write_sidecar(federation, "gwas", tmp_path / "gwas_watch.json")
    payload = json.loads(path.read_text(encoding="utf-8"))

    entry = federation.experiments["gwas"].sites[0]
    assert payload["clients"][entry.endpoint_uuid] == payload["clients"][entry.client_id]
    assert payload["clients"][entry.client_id]["institution"] == "Northern Institute of Technology"


def test_the_sidecar_carries_geography_the_worker_never_sends(federation, tmp_path):
    """A partner returns statistics. Where they are is joined in on the coordinator."""
    path = write_sidecar(federation, "gwas", tmp_path / "gwas_watch.json")
    site = json.loads(path.read_text(encoding="utf-8"))["clients"]["Site1"]
    assert (site["lat"], site["lng"]) == (43.6532, -79.3832)


def test_no_sidecar_means_an_inert_watcher():
    watcher = RunWatcher.load(None)
    watcher.start()
    watcher.round_start(1)
    watcher.client_update("Site1", 1, local_loss=0.5)
    watcher.round_end(1)
    watcher.finish()
    assert not watcher.enabled


def test_an_unreadable_sidecar_means_an_inert_watcher(tmp_path):
    broken = tmp_path / "watch.json"
    broken.write_text("{not json", encoding="utf-8")
    assert not RunWatcher.load(broken).enabled
    assert not RunWatcher.load(tmp_path / "absent.json").enabled


def test_an_emitter_that_raises_disables_the_map_and_not_the_run(tmp_path, caplog):
    """The whole contract of this module, in one test."""

    class Exploding:
        def on_init(self, *_):
            pass

        def on_client_update(self, *_):
            raise RuntimeError("emitter fell over")

    watcher = RunWatcher({"runs_dir": str(tmp_path), "server": {}})
    watcher._run = hivewatch.init(emitters=[Exploding()], verbose=False)
    watcher._live = True

    watcher.client_update("Site1", 1, local_loss=0.5)  # must not raise
    watcher.round_end(1)
    watcher.finish()

    assert not watcher.enabled


def test_a_run_update_merges_the_sites_identity_with_its_metrics(tmp_path):
    """The map needs both halves, and only the coordinator has both."""
    captured = []

    class Recording:
        def on_init(self, *_):
            pass

        def on_client_update(self, client):
            captured.append(client)

    sidecar = {
        "runs_dir": str(tmp_path),
        "server": {},
        "clients": {
            "endpoint-uuid": {
                "client_id": "Site1",
                "institution": "Northern Institute of Technology",
                "lat": 43.6532,
                "lng": -79.3832,
                "num_samples": 18032,
            }
        },
    }
    watcher = RunWatcher(sidecar)
    watcher._run = hivewatch.init(emitters=[Recording()], verbose=False)
    watcher._live = True
    watcher.client_update("endpoint-uuid", 1, local_loss=0.42, num_samples=17998)

    (client,) = captured
    assert client.client_id == "Site1"
    assert client.lat == 43.6532
    assert client.local_loss == 0.42
    # What the site actually reported wins over what the config predicted.
    assert client.num_samples == 17998
    assert client.extra["institution"] == "Northern Institute of Technology"


# ---------------------------------------------------------------------------
# registry uniformity
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", [n for n, s in REGISTRY.items() if s.implemented])
def test_every_implemented_experiment_says_what_data_it_holds(name):
    """The map has one line per site to answer 'what is behind that dot'. Fill it in."""
    assert REGISTRY[name].data_kind, (
        f"experiment '{name}' has no data_kind. It is the phrase shown beside every site "
        "running it on the network map -- a few words naming what a participating site "
        "actually holds on disk."
    )


def test_participations_are_keyed_by_site_not_by_client_id(federation):
    """Client ids are per-experiment; a site's identity across the federation is its id."""
    active = participations(federation)
    assert set(active) == {"site-north", "site-south", "site-east"}
    assert {p.client_id for p in active["site-north"]} == {"Site1", "anl"}
