"""Every shipped config parses, validates, and resolves.

A config that does not parse fails at launch time -- after a scheduler queue wait, on
someone else's cluster. Catching it here costs milliseconds.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from appfl_bio_suite.core.config import FederationError, load_federation
from appfl_bio_suite.core.experiments import REGISTRY, repo_root

EXAMPLE = repo_root() / "federation.yaml.example"


def _shipped_configs() -> list[Path]:
    out = []
    for spec in REGISTRY.values():
        if spec.configs_path.is_dir():
            out.extend(sorted(spec.configs_path.rglob("*.yaml")))
    return out


@pytest.mark.parametrize("path", _shipped_configs(), ids=[p.name for p in _shipped_configs()])
def test_shipped_config_parses(path: Path):
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert data is not None, f"{path} is empty"
    assert isinstance(data, dict), f"{path} must be a mapping at the top level"


def test_every_implemented_experiment_has_its_required_configs():
    for name, spec in REGISTRY.items():
        if not spec.implemented:
            continue
        for required in ("server.yaml", "clients.template.yaml"):
            assert (spec.configs_path / required).is_file(), (
                f"{name} is implemented but has no {required}"
            )


def test_every_implemented_experiment_has_a_loopback_config():
    """A stranger must be able to validate the install before recruiting anyone."""
    for name, spec in REGISTRY.items():
        if not spec.implemented:
            continue
        assert (spec.configs_path / "server.loopback.yaml").is_file(), (
            f"{name} has no loopback config. Without one there is no way to check an "
            "install end to end without an external partner."
        )


def test_loopback_configs_are_synchronous():
    """A serial run executes clients one after another; async has nothing to overlap."""
    for spec in REGISTRY.values():
        if not spec.implemented:
            continue
        path = spec.configs_path / "server.loopback.yaml"
        config = yaml.safe_load(path.read_text(encoding="utf-8"))
        scheduler = config["server_configs"].get("scheduler", "SyncScheduler")
        assert "Async" not in scheduler, f"{path} uses {scheduler}"


def test_gwas_server_config_is_single_round():
    """Summary-statistic FL has exactly one exchange. More rounds recompute the same
    result at full cost, which looks like a hang rather than an error."""
    spec = REGISTRY["gwas"]
    config = yaml.safe_load((spec.configs_path / "server.yaml").read_text(encoding="utf-8"))
    assert config["server_configs"]["num_global_epochs"] == 1


# ---------------------------------------------------------------------------
# federation.yaml
# ---------------------------------------------------------------------------


def test_example_federation_loads():
    federation = load_federation(EXAMPLE)
    assert federation.coordinator.identity
    assert federation.sites
    assert federation.experiments


def test_example_federation_is_obviously_fictional():
    """It is committed and public. It must not look like anyone's real deployment."""
    text = EXAMPLE.read_text(encoding="utf-8").lower()
    assert "example" in text
    federation = load_federation(EXAMPLE)
    assert "example" in federation.coordinator.identity.lower(), (
        "the example coordinator identity should be obviously fictional"
    )


def test_example_covers_both_scheduler_kinds():
    """Both provider shapes must be exercised, or one silently rots."""
    federation = load_federation(EXAMPLE)
    schedulers = {s.scheduler for s in federation.sites}
    assert {"slurm", "pbspro"} <= schedulers


def test_every_registered_experiment_appears_in_the_example():
    """Including the unimplemented one -- that is what proves the slot works."""
    federation = load_federation(EXAMPLE)
    for name in REGISTRY:
        assert name in federation.experiments, (
            f"{name} is registered but absent from federation.yaml.example, so nobody "
            "copying the example would know it exists."
        )


def test_configs_resolve_for_every_enabled_experiment():
    """End to end: federation.yaml -> the exact configs APPFL is handed."""
    from appfl_bio_suite.core.launch import resolve_run

    federation = load_federation(EXAMPLE)
    for name, experiment in federation.enabled_experiments().items():
        run = resolve_run(federation, name)
        assert run.num_clients == len(experiment.sites)
        assert run.server_config["server_configs"]["num_clients"] == run.num_clients

        for client in run.client_configs:
            # dataset_path is resolved on the DRIVER and must be absolute -- APPFL
            # resolves it against the driver's working directory.
            assert Path(client["data_configs"]["dataset_path"]).is_absolute()
            assert Path(client["data_configs"]["dataset_path"]).is_file()
            # Worker-side paths must be absolute on the partner's own cluster.
            assert client["train_configs"]["logging_output_dirname"].startswith("/")


def test_no_config_path_is_left_unresolved():
    """A `*_path: null` that survives resolution fails on the worker, not here."""
    from appfl_bio_suite.core.launch import resolve_run

    federation = load_federation(EXAMPLE)
    for name in federation.enabled_experiments():
        run = resolve_run(federation, name)

        # `experiment` is bound as a default argument rather than closed over: a closure
        # would capture the loop variable by reference and every message would name the
        # last experiment, which is exactly the kind of misleading failure output that
        # wastes time when a test finally does fail.
        def walk(node, trail="", experiment=name):
            if isinstance(node, dict):
                for key, value in node.items():
                    if key.endswith("_path"):
                        assert value is not None, f"{experiment}: {trail}.{key} is still null"
                    walk(value, f"{trail}.{key}", experiment)
            elif isinstance(node, list):
                for item in node:
                    walk(item, trail, experiment)

        walk(run.server_config)


# ---------------------------------------------------------------------------
# validation catches real mistakes
# ---------------------------------------------------------------------------


def _example_dict() -> dict:
    return yaml.safe_load(EXAMPLE.read_text(encoding="utf-8"))


def test_duplicate_client_ids_are_rejected(tmp_path):
    from appfl_bio_suite.core.config import Federation

    data = _example_dict()
    sites = data["experiments"]["gwas"]["sites"]
    sites[1]["client_id"] = sites[0]["client_id"]
    with pytest.raises(Exception, match="duplicate client_id"):
        Federation.model_validate(data)


def test_two_sites_on_the_same_center_are_rejected():
    """Two sites on one shard means training twice on the same data."""
    from appfl_bio_suite.core.config import Federation

    data = _example_dict()
    sites = data["experiments"]["flamby-heart-disease"]["sites"]
    sites[1]["center"] = sites[0]["center"]
    with pytest.raises(Exception, match="same dataset center"):
        Federation.model_validate(data)


def test_a_center_beyond_num_clients_is_rejected():
    from appfl_bio_suite.core.config import Federation

    data = _example_dict()
    data["experiments"]["flamby-heart-disease"]["sites"][0]["center"] = 99
    with pytest.raises(Exception, match="num_clients"):
        Federation.model_validate(data)


def test_relative_worker_paths_are_rejected():
    """A relative path resolves inside the endpoint's task working directory."""
    from appfl_bio_suite.core.config import Federation

    data = _example_dict()
    data["experiments"]["gwas"]["sites"][0]["data_dir"] = "relative/path"
    with pytest.raises(Exception, match="absolute path"):
        Federation.model_validate(data)


def test_a_placeholder_uuid_is_rejected_with_a_useful_message():
    from appfl_bio_suite.core.config import Federation

    data = _example_dict()
    data["experiments"]["gwas"]["sites"][0]["endpoint_uuid"] = "REPLACE_WITH_UUID"
    with pytest.raises(Exception, match="placeholder"):
        Federation.model_validate(data)


def test_an_unknown_site_reference_is_rejected():
    from appfl_bio_suite.core.config import Federation

    data = _example_dict()
    data["experiments"]["gwas"]["sites"][0]["site"] = "site-that-does-not-exist"
    with pytest.raises(Exception, match="not in the top-level"):
        Federation.model_validate(data)


def test_a_typo_in_a_key_is_rejected():
    """Unknown keys are refused, so a typo is not silently ignored.

    Without this, a misspelled key means a default is quietly applied where the user
    thought they had set something -- a failure that only shows up on a partner's cluster.
    """
    from pydantic import ValidationError

    from appfl_bio_suite.core.config import Federation

    data = _example_dict()
    data["coordinator"]["identtiy"] = "typo@example.org"
    with pytest.raises(ValidationError, match="identtiy"):
        Federation.model_validate(data)


def test_missing_federation_file_explains_how_to_create_one():
    with pytest.raises(FederationError) as excinfo:
        load_federation("/nonexistent/federation.yaml")
    assert "no federation config" in str(excinfo.value).lower()
