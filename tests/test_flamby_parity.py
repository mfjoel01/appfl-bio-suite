"""FLamby migration parity against the pre-migration configuration.

WHAT IS AND IS NOT CLAIMED HERE
-------------------------------
**Structural parity only.** These tests assert that this repository produces
configurations semantically equivalent to the ones the pre-migration tree carried, for the
same federation: same trainer, same optimizer and learning rate, same local step count,
same batch sizes, same aggregator, same weighting mode, same per-site centre assignment.

**Numerical parity is NOT established, and cannot be.** The pre-migration tree contains no
completed run of this experiment -- its only output artifact holds a single logging header
line and nothing else. There is no baseline result to compare against.

That is recorded plainly rather than glossed, because the tempting mistake is to treat a
first green run from this repository as evidence of parity. It would be evidence that this
repository works, which is a weaker and different claim. See
docs/experiments/flamby-heart-disease/ABOUT.md#record-of-runs.

Tests needing the pre-migration tree skip unless APPFL_BIO_SUITE_LEGACY_FLAMBY points at
it. The tests that do not need it -- the invariants below -- always run, and they are the
ones that keep protecting the migration after the old trees are gone.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
import yaml

from appfl_bio_suite.core.config import load_federation
from appfl_bio_suite.core.experiments import repo_root
from appfl_bio_suite.core.launch import build_client_configs, build_server_config

EXAMPLE = repo_root() / "federation.yaml.example"
EXPERIMENT = "flamby-heart-disease"

_LEGACY_ENV = "APPFL_BIO_SUITE_LEGACY_FLAMBY"
_LEGACY = Path(os.environ.get(_LEGACY_ENV, "/nonexistent"))
_LEGACY_CONFIGS = _LEGACY / "examples" / "resources" / "config_gc" / "flamby" / "heart_disease"

needs_legacy = pytest.mark.skipif(
    not _LEGACY_CONFIGS.is_dir(),
    reason=f"set {_LEGACY_ENV} to the pre-migration APPFL_FLamby checkout",
)


@pytest.fixture(scope="module")
def federation():
    return load_federation(EXAMPLE)


@pytest.fixture(scope="module")
def server_config(federation):
    return build_server_config(federation, EXPERIMENT)


@pytest.fixture(scope="module")
def client_configs(federation):
    return build_client_configs(federation, EXPERIMENT)


# ---------------------------------------------------------------------------
# invariants -- always run, and outlive the pre-migration trees
# ---------------------------------------------------------------------------


def test_training_hyperparameters_are_unchanged(server_config):
    """The values the pre-migration configs used. Changing one changes the experiment."""
    train = server_config["client_configs"]["train_configs"]
    assert train["trainer"] == "VanillaTrainer"
    assert train["mode"] == "step"
    assert train["num_local_steps"] == 100
    assert train["optim"] == "Adam"
    assert train["optim_args"]["lr"] == 0.001
    assert train["train_batch_size"] == 4
    assert train["val_batch_size"] == 4


def test_aggregation_is_unchanged(server_config):
    server = server_config["server_configs"]
    assert server["aggregator"] == "FedAvgAggregator"
    assert server["scheduler"] == "SyncScheduler"
    assert server["scheduler_kwargs"]["same_init_model"] is True


def test_model_loss_and_metric_are_unchanged(server_config):
    assert server_config["client_configs"]["model_configs"]["model_name"] == "Baseline"
    assert server_config["client_configs"]["train_configs"]["loss_fn_name"] == "BaselineLoss"
    assert server_config["client_configs"]["train_configs"]["metric_name"] == "metric"


def test_dataset_kwargs_match_the_flamby_loader_contract(client_configs):
    """`num_clients` is the dataset's split size, not the number of sites."""
    for client in client_configs:
        kwargs = client["data_configs"]["dataset_kwargs"]
        assert kwargs["dataset"] == "HeartDisease"
        assert kwargs["num_clients"] == 4
        assert 0 <= kwargs["client_id"] < 4


def test_each_site_gets_a_distinct_centre(client_configs):
    centres = [c["data_configs"]["dataset_kwargs"]["client_id"] for c in client_configs]
    assert len(centres) == len(set(centres)), (
        f"two sites share a centre: {centres}. The federation would train twice on the "
        "same data while reporting otherwise."
    )


def test_worker_side_paths_are_absolute(client_configs):
    """Resolved on the partner's worker, inside the endpoint's task working directory."""
    for client in client_configs:
        assert client["train_configs"]["logging_output_dirname"].startswith("/")


def test_driver_side_paths_are_absolute_and_exist(client_configs, server_config):
    """Read on the driver, so they must resolve here regardless of working directory."""
    for client in client_configs:
        path = Path(client["data_configs"]["dataset_path"])
        assert path.is_absolute() and path.is_file()

    for key, node in (
        ("model_path", server_config["client_configs"]["model_configs"]),
        ("loss_fn_path", server_config["client_configs"]["train_configs"]),
        ("metric_path", server_config["client_configs"]["train_configs"]),
    ):
        assert Path(node[key]).is_absolute() and Path(node[key]).is_file()


def test_no_client_config_asks_the_partner_for_a_dataset_path(federation):
    """The pre-migration setup guide asked partners for `dataset_path`. It was wrong.

    That path is resolved on the coordinator's driver, which reads the file and ships its
    source. A partner supplying one supplies a value that is never used.
    """
    from appfl_bio_suite.core.partner import _client_config_block

    for entry in federation.experiment(EXPERIMENT).sites:
        block = _client_config_block(federation, EXPERIMENT, entry.client_id)
        assert "dataset_path" not in block["data_configs"]


# ---------------------------------------------------------------------------
# direct comparison against the pre-migration configs
# ---------------------------------------------------------------------------


def _legacy(name: str) -> dict:
    return yaml.safe_load((_LEGACY_CONFIGS / name).read_text(encoding="utf-8"))


@needs_legacy
def test_server_config_matches_the_pre_migration_one(server_config):
    legacy = _legacy("server_fedavg_2site.yaml")
    legacy_train = legacy["client_configs"]["train_configs"]
    train = server_config["client_configs"]["train_configs"]

    for key in (
        "trainer",
        "mode",
        "num_local_steps",
        "optim",
        "train_batch_size",
        "val_batch_size",
        "do_validation",
    ):
        assert train[key] == legacy_train[key], f"{key} differs from the pre-migration config"

    assert train["optim_args"]["lr"] == legacy_train["optim_args"]["lr"]
    assert server_config["server_configs"]["aggregator"] == legacy["server_configs"]["aggregator"]
    assert (
        server_config["client_configs"]["model_configs"]["model_name"]
        == legacy["client_configs"]["model_configs"]["model_name"]
    )


@needs_legacy
def test_client_config_shape_matches_the_pre_migration_one(client_configs):
    """Same keys, same dataset kwargs shape -- differing only in per-site values."""
    legacy_clients = _legacy("clients_2site.yaml")["clients"]
    legacy_sample = legacy_clients[0]
    ours = client_configs[0]

    assert set(ours) >= {"endpoint_id", "client_id", "train_configs", "data_configs"}
    assert set(ours["data_configs"]) >= set(legacy_sample["data_configs"])
    assert set(ours["data_configs"]["dataset_kwargs"]) == set(
        legacy_sample["data_configs"]["dataset_kwargs"]
    )
    assert (
        ours["data_configs"]["dataset_kwargs"]["dataset"]
        == (legacy_sample["data_configs"]["dataset_kwargs"]["dataset"])
    )


def test_the_shipped_loader_is_behaviourally_equivalent():
    """Same dispatch and the same client cap as the pre-migration loader.

    The names differ -- `get_flamby` became `get_dataset`, which is uniform across
    experiments -- and the new one raises instead of asserting, with a message that names
    the fix. Behaviour for valid input is identical.

    Not gated on the legacy tree: it asserts the NEW loader's behaviour, described here
    by reference to the old one. Gating it meant the only test of that behaviour never
    ran anywhere -- not in CI, and not for anyone without the pre-migration checkout.
    """
    from appfl_bio_suite.experiments.flamby_heart_disease import dataset

    with pytest.raises(ValueError, match="4 centers"):
        dataset.get_dataset(dataset="HeartDisease", num_clients=5, client_id=0)

    with pytest.raises(NotImplementedError, match="HeartDisease"):
        dataset.get_dataset(dataset="TcgaBrca", num_clients=6, client_id=0)


def test_numerical_parity_is_explicitly_unverified():
    """A standing reminder, in the test suite where it cannot be overlooked.

    If a baseline run of the pre-migration tree ever exists, replace this with a real
    comparison and update ABOUT.md's record of runs.
    """
    about = (repo_root() / "docs" / "experiments" / EXPERIMENT / "ABOUT.md").read_text(
        encoding="utf-8"
    )
    assert "no completed baseline" in about.lower(), (
        "ABOUT.md must keep stating that no baseline run exists for this experiment, so "
        "that a first green run from this repository is not mistaken for parity."
    )
