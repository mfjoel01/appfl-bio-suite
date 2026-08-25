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
