"""Offline regression checks for the live authorization probe's controls."""

import argparse
import importlib.util
import sqlite3
from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.fixture
def probe(tmp_path, monkeypatch):
    path = Path(__file__).resolve().parents[1] / "scripts/globus_service_account_probe.py"
    spec = importlib.util.spec_from_file_location("service_account_probe", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "STATE_DIR", tmp_path)
    return module


def test_shared_secrets_have_separate_caches_and_fresh_preserves_them(probe, monkeypatch):
    import globus_sdk
    from globus_sdk.token_storage import MemoryTokenStorage

    monkeypatch.setattr(globus_sdk, "ClientApp", lambda **kwargs: SimpleNamespace(**kwargs))
    creds = {"client_id": "shared-client", "credential_id": "alpha", "secret": "test-only"}
    alpha = probe._probe_app(creds).config.token_storage
    beta = probe._probe_app({**creds, "credential_id": "beta"}).config.token_storage
    assert alpha.namespace != beta.namespace
    assert alpha.get_token_data_by_resource_server() == {}
    db = probe.STATE_DIR / "probe_tokens.db"
    with sqlite3.connect(db) as connection:
        connection.execute(
            "INSERT INTO token_storage VALUES (?, ?, ?)",
            (alpha.namespace, "fake-resource", "{}"),
        )
    assert beta.get_token_data_by_resource_server() == {}
    assert isinstance(probe._probe_app(creds, fresh=True).config.token_storage, MemoryTokenStorage)
    with sqlite3.connect(db) as connection:
        assert connection.execute("SELECT count(*) FROM token_storage").fetchone()[0] == 1
    assert db.stat().st_mode & 0o777 == 0o600


def test_client_identity_is_explicit_without_environment_switching(probe, monkeypatch):
    import globus_compute_sdk

    monkeypatch.setenv("GLOBUS_COMPUTE_CLIENT_ID", "operator-client")
    monkeypatch.setenv("GLOBUS_COMPUTE_CLIENT_SECRET", "operator-secret")
    creds = {"client_id": "site-client", "credential_id": "site-secret", "secret": "test-only"}
    monkeypatch.setattr(probe, "_load_state", lambda _: {"sites": {"alpha": creds}})
    sentinel = object()
    monkeypatch.setattr(probe, "_probe_app", lambda creds, fresh: sentinel)
    monkeypatch.setattr(globus_compute_sdk, "Client", lambda **kwargs: SimpleNamespace(**kwargs))
    client, _ = probe._client_as("probe", "alpha", fresh=True)
    assert client.app is sentinel
    assert probe.os.environ["GLOBUS_COMPUTE_CLIENT_ID"] == "operator-client"
    assert probe.os.environ["GLOBUS_COMPUTE_CLIENT_SECRET"] == "operator-secret"


def test_whoami_uses_same_app_as_compute(probe, monkeypatch, capsys):
    from globus_compute_sdk.sdk.auth import auth_client

    app = object()
    creds = {"secret": "test-only", "credential_id": "alpha"}
    monkeypatch.setattr(probe, "_client_as", lambda *a, **k: (SimpleNamespace(app=app), creds))

    def auth(**kwargs):
        assert kwargs["app"] is app
        return SimpleNamespace(userinfo=lambda: {"sub": "site-identity"})

    monkeypatch.setattr(auth_client, "ComputeAuthClient", auth)
    probe.cmd_whoami(SimpleNamespace(as_site="alpha", federation="probe", fresh=True))
    assert "site-identity" in capsys.readouterr().out


def test_timeout_is_inconclusive_and_task_group_is_reset(probe, monkeypatch, capsys):
    import globus_compute_sdk

    monkeypatch.setattr(probe, "_client_as", lambda *a, **k: (object(), {"secret": "test-only"}))

    class FakeExecutor:
        task_group_id = "another-identity-group"

        def __init__(self, **kwargs):
            pass

        def submit(self, fn):
            assert self.task_group_id is None
            return self

        def result(self, timeout):
            raise TimeoutError()

        def shutdown(self, **kwargs):
            pass

    monkeypatch.setattr(globus_compute_sdk, "Executor", FakeExecutor)
    args = SimpleNamespace(as_site="alpha", federation="probe", fresh=True, to="uuid", timeout=1)
    assert probe.cmd_submit(args) == 2
    assert "INCONCLUSIVE" in capsys.readouterr().out


def test_provision_cannot_overwrite_existing_credentials(probe, monkeypatch):
    path = probe._state_file("existing")
    path.write_text("existing state")
    monkeypatch.setattr(probe, "_auth_client", lambda: pytest.fail("must not contact Auth"))
    with pytest.raises(SystemExit, match="already exists"):
        probe.cmd_provision(argparse.Namespace(federation="existing"))
    assert path.read_text() == "existing state"
