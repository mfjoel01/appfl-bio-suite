"""Credential isolation and the real APPFL round loop over an offline transport."""

import json
from concurrent.futures import Future
from types import SimpleNamespace
from uuid import uuid4

import pytest

from appfl_bio_suite.core.site_credentials import (
    SiteCredential,
    SiteExecutorRouter,
    load_site_credentials,
)

ENDPOINTS = ["10000000-0000-0000-0000-000000000001", "10000000-0000-0000-0000-000000000002"]
CLIENTS = ["20000000-0000-0000-0000-000000000001", "20000000-0000-0000-0000-000000000002"]


@pytest.fixture
def manifest(tmp_path, monkeypatch):
    raw = {"version": 1, "endpoints": {}}
    for i, endpoint in enumerate(ENDPOINTS):
        name = f"TEST_SITE_{i}_SECRET"
        monkeypatch.setenv(name, f"fake-secret-{i}")
        raw["endpoints"][endpoint] = {
            "client_id": CLIENTS[i],
            "client_secret_env": name,
            "expected_user": "testuser",
        }
    path = tmp_path / "routes.json"
    path.write_text(json.dumps(raw))
    return path


def test_routes_never_contain_secret_values(manifest):
    routes = load_site_credentials(manifest, ENDPOINTS)
    assert [routes[e].client_id for e in ENDPOINTS] == CLIENTS
    assert "fake-secret" not in repr(routes)
    assert "fake-secret" not in manifest.read_text()


def test_serial_simulation_rejects_site_credentials(manifest, monkeypatch):
    from click.testing import CliRunner

    from appfl_bio_suite import cli

    monkeypatch.setattr(
        cli, "get_spec", lambda _: SimpleNamespace(implemented=True, serial_runner="unused")
    )
    result = CliRunner().invoke(
        cli.main, ["run", "fine-mapping", "--site-credentials", str(manifest)]
    )
    assert result.exit_code == 1
    assert "not a serial simulation" in result.output


@pytest.mark.parametrize("mutation", ["missing", "shared", "inline", "bad_env", "unknown"])
def test_invalid_routes_fail_before_authentication(manifest, mutation):
    raw = json.loads(manifest.read_text())
    entry = raw["endpoints"][ENDPOINTS[1]]
    if mutation == "missing":
        del raw["endpoints"][ENDPOINTS[1]]
    elif mutation == "shared":
        entry["client_id"] = CLIENTS[0]
    elif mutation == "inline":
        entry["client_secret"] = "must-not-appear-in-errors"
    elif mutation == "bad_env":
        entry["client_secret_env"] = "$INVALID"
    else:
        entry["unexpected"] = True
    manifest.write_text(json.dumps(raw))
    with pytest.raises(ValueError) as error:
        load_site_credentials(manifest, ENDPOINTS)
    assert "must-not-appear" not in str(error.value)


def test_missing_environment_fails_closed(manifest, monkeypatch):
    monkeypatch.delenv("TEST_SITE_1_SECRET")
    with pytest.raises(ValueError, match="TEST_SITE_1_SECRET"):
        load_site_credentials(manifest, ENDPOINTS)


@pytest.fixture
def transport(monkeypatch):
    import globus_compute_sdk

    instances = []

    class OfflineExecutor:
        def __init__(self, endpoint_id, client):
            self.endpoint_id = endpoint_id
            self.client = client
            self.task_group_id = "inherited-other-identity"
            self.calls = []
            self.closed = False
            instances.append(self)

        def submit(self, fn, *args, **kwargs):
            assert not self.closed
            if not self.calls:
                assert self.task_group_id is None
                self.task_group_id = str(uuid4())
            self.calls.append(kwargs)
            future = Future()
            try:
                future.set_result(fn(*args, **kwargs))
            except Exception as exc:
                future.set_exception(exc)
            return future

        def shutdown(self, **kwargs):
            self.closed = True

    monkeypatch.setattr(globus_compute_sdk, "Executor", OfflineExecutor)
    monkeypatch.setattr(SiteCredential, "client", lambda self: SimpleNamespace(id=self.client_id))
    return instances


def test_router_preserves_accounts_groups_and_rejects_unknown_targets(manifest, transport):
    router = SiteExecutorRouter(load_site_credentials(manifest, ENDPOINTS))
    for endpoint in [*ENDPOINTS, ENDPOINTS[0]]:
        router.endpoint_id = endpoint
        assert router.submit(lambda: 42).result() == 42
    assert {ex.endpoint_id: ex.client.id for ex in transport} == dict(
        zip(ENDPOINTS, CLIENTS, strict=True)
    )
    assert len({ex.task_group_id for ex in transport}) == 2
    with pytest.raises(ValueError):
        router.endpoint_id = str(uuid4())
    router.shutdown(wait=False, cancel_futures=True)
    assert all(ex.closed for ex in transport)
    with pytest.raises(RuntimeError):
        router.submit(lambda: 42)


def test_shutdown_attempts_every_executor(manifest, transport):
    router = SiteExecutorRouter(load_site_credentials(manifest, ENDPOINTS))

    def fail(**kwargs):
        raise RuntimeError("transport failed")

    transport[0].shutdown = fail
    with pytest.raises(RuntimeError, match="1 site"):
        router.shutdown(wait=False)
    assert transport[1].closed


def test_real_appfl_two_site_two_round_training(manifest, transport, tmp_path, monkeypatch):
    """Real server, communicator, worker training, and aggregation; no Globus network."""
    from appfl_bio_suite.core.drivers.globus_compute_driver import main
    from appfl_bio_suite.core.site_trial import prepare_trial

    monkeypatch.setenv("APPFL_BIO_SUITE_GLOBUS_CREDENTIALS", str(manifest))
    server, clients = prepare_trial(manifest, tmp_path / "trial")
    assert (
        main(["--server-config", str(server), "--client-config", str(clients), "--get-sample-size"])
        == 0
    )
    assert len(transport) == 2
    for executor in transport:
        assert [call["task_name"] for call in executor.calls] == [
            "get_sample_size",
            "train",
            "train",
        ]
        assert "fake-secret" not in repr(executor.calls)
        assert "client_secret_env" not in repr(executor.calls)
        assert executor.closed
    assert len({ex.task_group_id for ex in transport}) == 2


def test_trial_does_not_require_secrets_or_overwrite(manifest, tmp_path, monkeypatch):
    from appfl_bio_suite.core.site_trial import prepare_trial

    monkeypatch.delenv("TEST_SITE_0_SECRET")
    monkeypatch.delenv("TEST_SITE_1_SECRET")
    output = tmp_path / "trial"
    server, clients = prepare_trial(manifest, output)
    assert "fake-secret" not in server.read_text() + clients.read_text()
    with pytest.raises(FileExistsError):
        prepare_trial(manifest, output)


@pytest.mark.parametrize(
    "failure,expected",
    [
        ("ENDPOINT_ACCESS_FORBIDDEN", 0),
        ("FUNCTION_NOT_PERMITTED", 1),
        ("timeout", 1),
        ("executed", 1),
    ],
)
def test_isolation_check_requires_specific_denial(manifest, monkeypatch, failure, expected):
    import globus_compute_sdk

    from appfl_bio_suite.core import endpoint, site_trial

    monkeypatch.setattr(SiteCredential, "client", lambda self: object())
    monkeypatch.setattr(
        endpoint,
        "smoke_test",
        lambda *a, **k: endpoint.SmokeResult(
            "endpoint", "site", True, 0, {"user": "testuser", "python": "3.12.13"}
        ),
    )

    class Denied(Exception):
        http_status = 403
        code = failure

    class Executor:
        def __init__(self, **kwargs):
            self.task_group_id = "wrong-owner"

        def submit(self, function):
            assert self.task_group_id is None
            return self

        def result(self, **kwargs):
            if failure == "executed":
                return {}
            if failure == "timeout":
                raise TimeoutError()
            raise Denied()

        def shutdown(self, **kwargs):
            pass

    monkeypatch.setattr(globus_compute_sdk, "Executor", Executor)
    assert site_trial.check_sites(manifest, timeout=1) == expected


def test_smoke_cli_routes_credentials_and_fails_wrong_local_user(manifest, monkeypatch):
    from click.testing import CliRunner

    from appfl_bio_suite import cli
    from appfl_bio_suite.core.endpoint import SmokeResult

    monkeypatch.setattr(cli, "_load", lambda _: SimpleNamespace(experiments={}))
    monkeypatch.setattr(cli, "_resolve_targets", lambda *a: [("a", ENDPOINTS[0])])
    monkeypatch.setattr(SiteCredential, "client", lambda self: self.client_id)

    def smoke(uuid, *, client, **kwargs):
        assert client == CLIENTS[0]
        return SmokeResult(uuid, "a", True, 0, {"user": "wrong-account"})

    monkeypatch.setattr("appfl_bio_suite.core.endpoint.smoke_test", smoke)
    result = CliRunner().invoke(
        cli.main, ["endpoint", "smoke", "--site-credentials", str(manifest)]
    )
    assert result.exit_code == 1
    assert "expected local user" in result.output


def test_launch_passes_manifest_only_in_coordinator_environment(manifest, monkeypatch, tmp_path):
    from appfl_bio_suite.core import launch, preflight

    monkeypatch.setattr(
        preflight, "run_preflight", lambda **k: SimpleNamespace(ok=True, render=lambda: "ok")
    )
    monkeypatch.setattr(launch, "enforce_data_use", lambda *a: (True, ""))
    run = SimpleNamespace(
        num_clients=2,
        server_config={},
        client_configs=[
            {"client_id": f"site-{i}", "endpoint_id": endpoint}
            for i, endpoint in enumerate(ENDPOINTS)
        ],
    )
    monkeypatch.setattr(launch, "resolve_run", lambda *a: run)
    monkeypatch.setattr(
        launch,
        "write_run_configs",
        lambda *a: (tmp_path / "server.yaml", tmp_path / "clients.yaml"),
    )
    calls = []
    monkeypatch.setattr(launch.subprocess, "call", lambda cmd, env: calls.append((cmd, env)) or 0)
    fed = SimpleNamespace(
        experiment=lambda _: SimpleNamespace(
            sites=[SimpleNamespace(endpoint_uuid=e) for e in ENDPOINTS]
        )
    )
    assert launch.launch(fed, "gwas", site_credentials=manifest) == 0
    assert calls[0][1]["APPFL_BIO_SUITE_GLOBUS_CREDENTIALS"] == str(manifest.resolve())
    assert "fake-secret" not in repr(calls[0][0])
    assert "credentials" not in repr(run.client_configs)
