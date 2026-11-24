"""What the launcher hands to APPFL, and what a loopback run resolves to.

Every case here is a silent failure rather than a loud one: a weighting mode that is
configured and then ignored, a federation setting that is declared and never applied, a
data path that is off by one directory. None of them raise; they just produce a different
run from the one the config describes.
"""

from __future__ import annotations

import pytest

from appfl_bio_suite.core.config import FederationError, load_federation
from appfl_bio_suite.core.experiments import repo_root
from appfl_bio_suite.core.launch import _driver_command, build_client_configs, build_server_config
from appfl_bio_suite.core.loopback import localize_for_loopback, loopback_federation

EXAMPLE = repo_root() / "federation.yaml.example"


@pytest.fixture(scope="module")
def federation():
    return load_federation(EXAMPLE)


# --- driver invocation ----------------------------------------------------


def test_sample_size_weighting_asks_the_driver_for_sample_sizes(tmp_path):
    """Without the flag APPFL falls back to 1/num_clients, silently.

    `client_weights_mode: sample_size` only takes effect if the server was told each
    client's sample size, and it is told by an extra round trip the Globus Compute driver
    makes only on request.
    """
    config = {"server_configs": {"aggregator_kwargs": {"client_weights_mode": "sample_size"}}}
    cmd = _driver_command("globus_compute", tmp_path / "s.yaml", tmp_path / "c.yaml", config)
    assert "--get-sample-size" in cmd


def test_equal_weighting_does_not_pay_for_a_round_trip_it_does_not_need(tmp_path):
    config = {"server_configs": {"aggregator_kwargs": {"client_weights_mode": "equal"}}}
    cmd = _driver_command("globus_compute", tmp_path / "s.yaml", tmp_path / "c.yaml", config)
    assert "--get-sample-size" not in cmd


# --- federation settings actually reaching the configs --------------------


def test_gwas_analysis_settings_come_from_the_federation_config(federation):
    """`variant_scaling` and `hit_p_threshold` were declared, documented, and ignored."""
    exp = federation.experiment("gwas")
    exp.variant_scaling = 0.25
    exp.hit_p_threshold = 1e-6
    try:
        clients = build_client_configs(federation, "gwas")
        server = build_server_config(federation, "gwas")
    finally:
        exp.variant_scaling = 1.0
        exp.hit_p_threshold = 5.0e-8

    for client in clients:
        assert client["train_configs"]["variant_scaling"] == 0.25
        assert client["train_configs"]["hit_p_threshold"] == 1e-6
    assert server["server_configs"]["aggregator_kwargs"]["hit_p_threshold"] == 1e-6


def test_every_site_gets_the_same_variant_scaling(federation):
    """The aggregator refuses payloads computed over different variant sets."""
    clients = build_client_configs(federation, "gwas")
    scalings = {c["train_configs"]["variant_scaling"] for c in clients}
    assert len(scalings) == 1


# --- loopback -------------------------------------------------------------


def _make_site_tree(root, names):
    for name in names:
        (root / name / "data").mkdir(parents=True)
    return root


def test_loopback_data_dir_points_at_where_the_simulator_writes(tmp_path):
    """`simulate --out X` writes X/<site>/data; the loader needs that directory itself.

    Resolving to X/<site> instead put every loopback run one level above its own data,
    and the loader failed on every site.
    """
    root = _make_site_tree(tmp_path / "sim", ["Site1", "Site2"])
    fed = loopback_federation("gwas", data_root=root)

    entries = fed.experiment("gwas").sites
    assert [e.client_id for e in entries] == ["Site1", "Site2"]
    for entry in entries:
        assert entry.data_dir == str(root / entry.client_id / "data")


def test_loopback_needs_no_federation_config(tmp_path):
    """The documented flow runs the loopback BEFORE local/federation.yaml exists."""
    root = _make_site_tree(tmp_path / "sim", ["Site1", "Site2"])
    fed = loopback_federation("gwas", data_root=root)

    clients = build_client_configs(fed, "gwas")
    assert len(clients) == 2
    assert build_server_config(fed, "gwas", "loopback")["server_configs"]["num_clients"] == 2


def test_loopback_says_what_to_do_when_the_data_is_not_there(tmp_path):
    with pytest.raises(FederationError, match="simulate gwas"):
        loopback_federation("gwas", data_root=tmp_path / "never-simulated")


def test_localizing_moves_output_off_the_partners_cluster(federation, tmp_path):
    """A real config's output_dir is a path on someone else's machine.

    Running it serially here without rewriting it means trying to create a directory
    under a partner's service-account home, which fails with a PermissionError that has
    nothing to do with what was being tested.
    """
    root = _make_site_tree(tmp_path / "sim", ["Site1", "Site2", "Site3"])
    fed = load_federation(EXAMPLE)

    before = [e.output_dir for e in fed.experiment("gwas").sites]
    assert any(p.startswith("/home/gwas_svc") for p in before)

    localize_for_loopback(fed, "gwas", data_root=root, out_root=tmp_path / "out")

    for entry in fed.experiment("gwas").sites:
        assert entry.output_dir.startswith(str(tmp_path / "out"))
        assert entry.data_dir == str(root / entry.client_id / "data")


def test_localizing_moves_the_drs_registry_with_the_objects(tmp_path):
    """A config started from the shipped example must survive the documented loopback run.

    The example names a DRS registry under local/data/ and leaves each site's `drs_uri`
    empty. Localizing pointed every site at an object in the SIMULATED data's registry
    and kept resolving them in the federation's own -- a file that, for anyone who had
    just copied the example, does not exist. `run fine-mapping --driver serial` then
    failed before starting, on the first thing the README tells a new coordinator to do.
    """
    from appfl_bio_suite.core.ga4gh.drs import build_registry
    from appfl_bio_suite.core.ga4gh.resolve import resolve_drs_object

    fed = load_federation(EXAMPLE)
    exp = fed.experiment("fine-mapping")
    root = tmp_path / "sim"
    for client_id in [entry.client_id for entry in exp.sites]:
        data = root / client_id / "data"
        data.mkdir(parents=True)
        (data / "site_manifest.tsv").write_text(f"IID\n{client_id}\n")
        (data / "DATA_USE.json").write_text("{}")
    build_registry(root, hostname="drs.loopback.test").save(root / "drs_registry.json")

    declared_study = exp.ga4gh.data_use_request
    # A value a real coordinator would have recorded, which describes some other bundle.
    exp.sites[0].drs_uri = "drs://drs.example-university.edu/" + "a" * 64
    exp.sites[0].data_use_profile = "/coordinator/copies/anl/DATA_USE.json"

    localize_for_loopback(fed, "fine-mapping", data_root=root, out_root=tmp_path / "out")

    assert fed.ga4gh.drs.hostname == "drs.loopback.test"
    assert fed.ga4gh.drs.registry == str(root / "drs_registry.json")
    # Everything about the federation that is not the simulated data is left alone.
    assert fed.ga4gh.trs.publish_dir == "local/trs"
    assert fed.ga4gh.tes.timeout_seconds == 7200
    assert exp.ga4gh.data_use_request == declared_study
    for entry in exp.sites:
        assert entry.data_use_profile == str(root / entry.client_id / "data" / "DATA_USE.json")
        resolved = resolve_drs_object(fed, entry)
        assert resolved["name"] == entry.client_id
        assert resolved["self_uri"] == entry.drs_uri


def test_localizing_without_a_simulated_registry_clears_what_it_cannot_repoint(tmp_path):
    """A real `drs_uri` left in place would be checked against a bundle it never described."""
    fed = load_federation(EXAMPLE)
    exp = fed.experiment("fine-mapping")
    root = _make_site_tree(tmp_path / "sim", [entry.client_id for entry in exp.sites])
    exp.sites[0].drs_uri = "drs://drs.example-university.edu/" + "a" * 64
    exp.sites[0].data_use_profile = "/coordinator/copies/anl/DATA_USE.json"

    localize_for_loopback(fed, "fine-mapping", data_root=root, out_root=tmp_path / "out")

    assert [entry.drs_uri for entry in exp.sites] == [None, None, None]
    assert [entry.data_use_profile for entry in exp.sites] == [None, None, None]
    # No registry travelled with the data, so the federation's own setting is untouched.
    assert fed.ga4gh.drs.hostname == "drs.example-university.edu"


def test_secure_fine_mapping_configuration_reaches_both_roles(federation):
    fed = federation.model_copy(deep=True)
    exp = fed.experiment("fine-mapping")
    config = {"run_id": "launch-test", "participants": [s.client_id for s in exp.sites]}
    exp.secure_aggregation = config
    server = build_server_config(fed, "fine-mapping")
    public = server["server_configs"]["aggregator_kwargs"]["secure_aggregation"]
    assert public["run_id"] == config["run_id"]
    assert server["client_configs"]["train_configs"]["secure_aggregation"] == public
    for client in build_client_configs(fed, "fine-mapping"):
        assert client["train_configs"]["secure_aggregation"] == config
    exp.secure_aggregation = {**config, "participants": ["wrong", "roster"]}
    with pytest.raises(ValueError, match="participants"):
        build_server_config(fed, "fine-mapping")


def test_multisusie_backend_is_coordinator_only(federation):
    fed = federation.model_copy(deep=True)
    exp = fed.experiment("fine-mapping")
    exp.fine_mapping_backend = "multisusie"
    exp.multisusie_rho = 0.6
    exp.multisusie_min_abs_corr = 0.7
    kwargs = build_server_config(fed, "fine-mapping")["server_configs"]["aggregator_kwargs"]
    assert kwargs["backend"] == "multisusie"
    assert kwargs["multisusie_rho"] == 0.6
    assert kwargs["multisusie_min_abs_corr"] == 0.7
    for client in build_client_configs(fed, "fine-mapping"):
        assert "backend" not in client["train_configs"]
