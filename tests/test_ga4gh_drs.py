"""DRS: content addressing, the registry, the served API, and what verification catches.

The tests worth having here are the ones about *identity*: that the same bytes get the
same id at two sites, that a changed byte changes a bundle's id, and that the site-side
verifier notices a bundle that is not the one it was supposed to hold. Those are the
failures DRS was added to make visible, and each of them is silent without it.
"""

from __future__ import annotations

import json
import threading
import urllib.error
import urllib.request

import pytest

from appfl_bio_suite.core.ga4gh.drs import (
    DrsClient,
    DrsError,
    DrsRegistry,
    build_registry,
    parse_drs_uri,
    registry_for_files,
    serve,
    verify_object,
)


def make_bundles(root, sites=("anl", "covenant"), payload=b"genotypes"):
    """A directory shaped like `simulate` output: <root>/<site>/data/."""
    for site in sites:
        data = root / site / "data"
        (data / "phenotypes").mkdir(parents=True)
        (data / "site_genotypes.bed").write_bytes(payload + site.encode())
        (data / "site_manifest.tsv").write_text(f"FID\tIID\tsuperpopulation\n1\t{site}1\tEUR\n")
        # Identical at every site by construction -- the harmonization depends on it.
        (data / "reference_variants.tsv").write_text(
            "snp_id\tchrom\tbp\ta1\ta2\nrs1\t1\t100\tA\tG\n"
        )
        (data / "phenotypes" / "L0_a_rep0.pheno").write_text("FID\tIID\ty\n1\t1\t0.5\n")
    truth = root / "ground_truth"
    truth.mkdir()
    (truth / "causal_manifest.tsv").write_text("locus_id\tcausal_snp_ids\nL0\trs1\n")
    return root


@pytest.fixture
def registry(tmp_path):
    make_bundles(tmp_path)
    return build_registry(tmp_path, hostname="drs.test.org")


# ---------------------------------------------------------------------------
# content addressing
# ---------------------------------------------------------------------------


def test_a_blobs_id_is_its_sha256(registry, tmp_path):
    from appfl_bio_suite.core.ga4gh.drs import file_sha256

    path = tmp_path / "anl" / "data" / "site_genotypes.bed"
    obj = registry.get(file_sha256(path))
    assert obj.name == "site_genotypes.bed"
    assert obj.sha256 == file_sha256(path)


def test_identical_files_at_two_sites_are_one_object_with_two_paths(registry):
    """reference_variants.tsv is hard-linked into every bundle and must dedupe."""
    matching = [o for o in registry.objects.values() if o.name == "reference_variants.tsv"]
    assert len(matching) == 1
    file_methods = [m for m in matching[0].access_methods if m.type == "file"]
    assert len(file_methods) == 2


def test_a_bundle_id_is_stable_across_rebuilds(tmp_path):
    make_bundles(tmp_path)
    first = build_registry(tmp_path, hostname="drs.test.org").by_name("anl")
    second = build_registry(tmp_path, hostname="drs.test.org").by_name("anl")
    assert first.id == second.id


def test_changing_one_byte_changes_the_bundle_id(tmp_path):
    """The whole point: a stale bundle is a different object, not a similar one."""
    make_bundles(tmp_path)
    before = build_registry(tmp_path, hostname="drs.test.org").by_name("anl").id
    (tmp_path / "anl" / "data" / "site_manifest.tsv").write_text(
        "FID\tIID\tsuperpopulation\n1\tanl1\tAFR\n"
    )
    after = build_registry(tmp_path, hostname="drs.test.org").by_name("anl").id
    assert before != after


def test_adding_a_file_changes_the_bundle_id(tmp_path):
    """A Merkle id over names as well as contents, so an addition is not invisible."""
    make_bundles(tmp_path)
    before = build_registry(tmp_path, hostname="drs.test.org").by_name("anl").id
    (tmp_path / "anl" / "data" / "extra.tsv").write_text("surprise\n")
    after = build_registry(tmp_path, hostname="drs.test.org").by_name("anl").id
    assert before != after


def test_the_consent_file_is_never_part_of_the_object(tmp_path):
    """A profile carries its object's drs_uri, so it cannot be inside that object's id.

    The property this protects is user-visible: `simulate` registers the bundles and
    *then* writes each profile into one, so re-registering the same directory later --
    which `ga4gh drs register` does -- must mint the same ids. Without the exclusion it
    mints different ones, and "the same bytes get the same id" quietly stops holding for
    exactly the directories this project produces.
    """
    from appfl_bio_suite.core.ga4gh.duo import DataUseProfile, write_profile

    make_bundles(tmp_path)
    before = build_registry(tmp_path, hostname="drs.test.org").by_name("anl")
    write_profile(
        DataUseProfile(dataset_id="x", permission="DUO:0000042", drs_uri=before.self_uri),
        tmp_path / "anl" / "data" / "DATA_USE.json",
    )
    after = build_registry(tmp_path, hostname="drs.test.org").by_name("anl")

    assert after.id == before.id
    assert "DATA_USE.json" not in after.child_ids()


def test_two_sites_get_different_bundle_ids(registry):
    assert registry.by_name("anl").id != registry.by_name("covenant").id


def test_the_answer_key_is_registered_but_belongs_to_no_bundle(registry):
    """Coordinator-only: addressable so results can cite it, in nobody's bundle."""
    key = next(o for o in registry.objects.values() if o.name == "causal_manifest.tsv")
    for bundle in registry.bundles().values():
        assert key.id not in bundle.child_ids().values()


# ---------------------------------------------------------------------------
# the registry
# ---------------------------------------------------------------------------


def test_registry_round_trips(registry, tmp_path):
    path = registry.save(tmp_path / "drs_registry.json")
    reloaded = DrsRegistry.load(path)
    assert set(reloaded.objects) == set(registry.objects)
    assert reloaded.hostname == registry.hostname


def test_resolving_a_uri_from_another_service_is_refused(registry):
    """A DRS id is only unique within its service; crossing them silently is worse."""
    uri = registry.uri(registry.by_name("anl").id).replace("drs.test.org", "drs.other.org")
    with pytest.raises(DrsError, match="only unique within its service"):
        registry.resolve(uri)


def test_a_missing_object_names_how_to_rebuild(registry):
    with pytest.raises(DrsError, match="ga4gh drs register"):
        registry.get("0" * 64)


def test_compact_identifiers_are_refused_with_the_reason():
    with pytest.raises(DrsError, match="identifiers.org"):
        parse_drs_uri("drs://prefix:accession")


def test_a_url_is_not_a_drs_uri():
    with pytest.raises(DrsError, match="not a DRS URI"):
        parse_drs_uri("https://drs.test.org/abc")


def test_an_empty_directory_cannot_be_addressed(tmp_path):
    (tmp_path / "empty" / "data").mkdir(parents=True)
    with pytest.raises(DrsError, match="no site bundles|empty"):
        build_registry(tmp_path, hostname="drs.test.org")


# ---------------------------------------------------------------------------
# verification
# ---------------------------------------------------------------------------


def test_verify_is_clean_on_untouched_data(registry):
    for bundle in registry.bundles().values():
        assert verify_object(bundle) == []


def test_verify_notices_a_missing_member(registry, tmp_path):
    bundle = registry.by_name("anl")
    (tmp_path / "anl" / "data" / "site_manifest.tsv").unlink()
    assert any("missing member" in p for p in verify_object(bundle))


def test_verify_notices_a_corrupted_blob(registry, tmp_path):
    path = tmp_path / "anl" / "data" / "site_genotypes.bed"
    blob = registry.get(next(iter(
        [o.id for o in registry.objects.values() if o.local_path() == path.resolve()]
    )))
    path.write_bytes(b"corrupted")
    assert any("sha-256 mismatch" in p for p in verify_object(blob))


# ---------------------------------------------------------------------------
# the served API
# ---------------------------------------------------------------------------


@pytest.fixture
def server(registry):
    # Port 0: the OS picks a free one and we read it back. A fixed port races anything
    # else on the machine, including another copy of this test suite.
    httpd = serve(registry, "127.0.0.1", 0)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_address[1]}"
    httpd.shutdown()


def _get(url):
    with urllib.request.urlopen(url, timeout=10) as response:
        return json.loads(response.read().decode("utf-8"))


def test_service_info_declares_drs(server):
    info = _get(f"{server}/ga4gh/drs/v1/service-info")
    assert info["type"] == {"group": "org.ga4gh", "artifact": "drs", "version": "1.5.0"}


def test_objects_endpoint_returns_the_bundle(server, registry):
    bundle = registry.by_name("anl")
    payload = _get(f"{server}/ga4gh/drs/v1/objects/{bundle.id}")
    assert payload["name"] == "anl"
    assert {c["name"] for c in payload["contents"]} == set(bundle.child_ids())


def test_access_endpoint_returns_the_url_for_a_served_registry(tmp_path):
    """With an https base configured, /access/{id} hands back a fetchable URL."""
    make_bundles(tmp_path)
    registry = build_registry(
        tmp_path, hostname="drs.test.org", https_base="https://drs.test.org"
    )
    httpd = serve(registry, "127.0.0.1", 0)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}"
    try:
        blob = next(o for o in registry.objects.values() if o.name == "site_manifest.tsv")
        payload = _get(f"{base}/ga4gh/drs/v1/objects/{blob.id}/access/https")
        assert payload["url"].endswith(f"/files/{blob.id}")
    finally:
        httpd.shutdown()


def test_an_access_id_that_does_not_exist_is_a_404(server, registry):
    """No https base was configured here, so the server must refuse rather than invent."""
    blob = next(o for o in registry.objects.values() if o.name == "site_manifest.tsv")
    with pytest.raises(urllib.error.HTTPError) as exc:
        _get(f"{server}/ga4gh/drs/v1/objects/{blob.id}/access/https")
    assert exc.value.code == 404


def test_unknown_object_is_a_404_in_the_shape_drs_specifies(server):
    with pytest.raises(urllib.error.HTTPError) as exc:
        _get(f"{server}/ga4gh/drs/v1/objects/{'0' * 64}")
    body = json.loads(exc.value.read().decode("utf-8"))
    assert body["status_code"] == 404 and "msg" in body


def test_bytes_are_served_and_match_the_checksum(server, registry, tmp_path):
    from appfl_bio_suite.core.ga4gh.drs import file_sha256

    blob = next(o for o in registry.objects.values() if o.name == "site_manifest.tsv")
    with urllib.request.urlopen(f"{server}/files/{blob.id}", timeout=10) as response:
        data = response.read()
    assert file_sha256(tmp_path / "anl" / "data" / "site_manifest.tsv") == blob.sha256
    assert len(data) == blob.size


def test_client_over_http_matches_the_local_registry(server, registry):
    bundle = registry.by_name("anl")
    remote = DrsClient(base_url=server).get(bundle.id)
    local = DrsClient(registry=registry).get(bundle.id)
    assert remote.id == local.id
    assert remote.child_ids() == local.child_ids()


def test_a_client_needs_a_backend():
    with pytest.raises(DrsError, match="local registry or a base_url"):
        DrsClient()


# ---------------------------------------------------------------------------
# outputs
# ---------------------------------------------------------------------------


def test_results_can_be_addressed_too(tmp_path):
    results = tmp_path / "fed_fm_results.tsv"
    results.write_text("locus_id\tn_credible_sets\nL0\t1\n")
    registry = registry_for_files([results], "drs.test.org")
    obj = registry.by_name("fed_fm_results.tsv")
    assert obj is not None and obj.self_uri.startswith("drs://drs.test.org/")
