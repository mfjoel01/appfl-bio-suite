"""Single-user endpoints, host-local workers, mixed credential routes, and large specs.

The bundle path that a cloud VM takes: no batch scheduler, no root, an endpoint that runs
as the service account under a Globus service-account client, and a coordinator that
reaches it through a --site-credentials route. Endpoint processes and Globus calls are
replaced at their boundary; everything up to that boundary is the real code.
"""

import json
import os
import stat
import subprocess
import zipfile
from pathlib import Path
from types import SimpleNamespace
from uuid import uuid4

import pytest
import yaml
from click.testing import CliRunner

from appfl_bio_suite.cli import main
from appfl_bio_suite.core import onboarding as api
from appfl_bio_suite.core.config import FederationError, load_federation
from appfl_bio_suite.core.experiments import repo_root
from appfl_bio_suite.core.partner import generate_bundle, unresolved_fields
from appfl_bio_suite.core.site_credentials import SiteCredential, load_site_credentials

SERVICE = "fedlearn_svc"
CLIENT = "30000000-0000-0000-0000-000000000003"
ENDPOINT = "40000000-0000-0000-0000-000000000004"


def _federation(tmp_path, *, scheduler="local", mode="single-user"):
    raw = yaml.safe_load((repo_root() / "federation.yaml.example").read_text())
    for site in raw["sites"]:
        if site["id"] == "site-south":
            site["scheduler"] = scheduler
            if scheduler == "local":
                for field in ("partition", "queue", "account"):
                    site.pop(field, None)
    for experiment in raw["experiments"].values():
        for entry in experiment.get("sites", []):
            if entry["site"] == "site-south":
                entry["endpoint_mode"] = mode
    path = tmp_path / "federation.yaml"
    path.write_text(yaml.safe_dump(raw, sort_keys=False))
    return load_federation(path)


@pytest.fixture
def fed(tmp_path):
    return _federation(tmp_path)


@pytest.fixture
def su_dir(fed, tmp_path):
    return generate_bundle(fed, "gwas", "site-south", tmp_path / "bundles")


@pytest.fixture
def su_bundle(su_dir):
    return api.read_bundle(su_dir)


@pytest.fixture
def su_site(su_bundle, tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    return api.SiteConfig(data_dir=data, work_dir=tmp_path / "work")


def _report(bundle, **changes):
    return {
        "bundle_digest": api.digest(bundle),
        "python": "3.12.11",
        "versions": {**bundle.versions, "appfl-bio-suite": bundle.suite_version},
        "user": bundle.service_account,
        "data_readable": True,
        "work_writable": True,
        **changes,
    }


def _fake_install(monkeypatch, bundle, report=None):
    calls = []

    def fake_run(args, **kwargs):
        calls.append(args)
        if "venv" in args:
            target = Path(args[-1]) / "bin/python"
            target.parent.mkdir(parents=True)
            target.touch()
        return SimpleNamespace(stdout="3.12\n")

    monkeypatch.setattr(api, "run", fake_run)
    monkeypatch.setattr(api, "inspect_runtime", lambda *a: report or _report(bundle))
    return calls


# -- federation and bundle ---------------------------------------------------


def test_local_site_rejects_scheduler_settings(tmp_path):
    raw = yaml.safe_load((repo_root() / "federation.yaml.example").read_text())
    raw["sites"][1]["scheduler"] = "local"  # keeps partition and account
    path = tmp_path / "federation.yaml"
    path.write_text(yaml.safe_dump(raw))
    with pytest.raises(FederationError, match="scheduler: local but sets"):
        load_federation(path)


def test_single_user_local_bundle_declares_both(fed, su_dir, su_bundle):
    assert su_bundle.scheduler == "local"
    assert su_bundle.endpoint_mode == "single-user"
    assert {"data_dir", "work_dir"} <= set(su_bundle.defaults)
    assert not {"account", "partition", "queue"} & set(su_bundle.defaults)
    assert unresolved_fields(fed, "gwas", "site-south") == []
    template = (su_dir / "user_config_template.yaml.j2").read_text()
    assert "LocalProvider" in template and "REPLACE_WITH" not in template
    guide = (su_dir / "0-onboarding.md").read_text()
    assert "single-user" in guide and "--credentials-file" in guide
    assert "As root, using the prepared runtime" not in guide
    assert "does not apply to you" in (su_dir / "1-endpoint-setup.md").read_text()


def test_default_bundle_is_still_multi_user(tmp_path):
    bundle = api.read_bundle(
        generate_bundle(
            load_federation(repo_root() / "federation.yaml.example"),
            "gwas",
            "site-north",
            tmp_path,
        )
    )
    assert bundle.endpoint_mode == "multi-user" and bundle.scheduler == "slurm"
    assert {"account", "partition"} <= set(bundle.defaults)


def test_large_specification_reads_from_its_zip(su_dir, tmp_path, monkeypatch):
    # A DRS-verified fine-mapping bundle carries a checksum per data file: megabytes.
    spec = json.loads((su_dir / "onboarding.json").read_text())
    spec["client_config"]["padding"] = "x" * 2_500_000
    (su_dir / "onboarding.json").write_text(json.dumps(spec))
    archive = tmp_path / "site.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as zipped:
        for path in su_dir.iterdir():
            zipped.write(path, path.name)
    assert api.digest(api.read_bundle(archive)) == api.digest(api.read_bundle(su_dir))
    monkeypatch.setattr(api, "MAX_SPEC_BYTES", 1_000_000)
    with pytest.raises(ValueError, match="the limit is 1,000,000"):
        api.read_bundle(archive)


# -- local settings ----------------------------------------------------------


def test_local_settings_refuse_scheduler_fields(su_bundle, tmp_path):
    local = tmp_path / "site.yaml"
    local.write_text(yaml.safe_dump({"account": "alloc", "data_dir": "/d", "work_dir": "/w"}))
    with pytest.raises(ValueError, match="no scheduler; remove account"):
        api.resolve_site(su_bundle, local)


def test_slurm_names_what_is_missing(tmp_path):
    bundle = api.read_bundle(
        generate_bundle(_federation(tmp_path, scheduler="slurm"), "gwas", "site-south", tmp_path)
    )
    local = tmp_path / "site.yaml"
    local.write_text(yaml.safe_dump({"account": "", "data_dir": "/d", "work_dir": "/w"}))
    # A blank bundle default reads as missing and is named, not a pydantic regex error.
    bundle = bundle.model_copy(update={"defaults": {"account": "", "partition": ""}})
    with pytest.raises(ValueError, match=r"missing: account, partition"):
        api.resolve_site(bundle, local)


def test_join_local_prompts_only_for_paths(su_dir, tmp_path, monkeypatch):
    calls = []

    def setup(spec, config, state, **kwargs):
        calls.append(config)
        state.mkdir(parents=True, exist_ok=True)
        return state

    monkeypatch.setattr(api, "setup_site", setup)
    data, work = tmp_path / "d", tmp_path / "w"
    result = CliRunner().invoke(
        main, ["join", str(su_dir), "--no-profile"], input=f"{data}\n{work}\n"
    )
    assert result.exit_code == 0, result.output
    assert "Scheduler account" not in result.output
    assert calls[0].account is None and calls[0].data_dir == data
    assert "Endpoint mode: single-user; scheduler: local" in result.output
    assert "--credentials-file" in result.output and "Administrator" not in result.output


# -- setup and activation ----------------------------------------------------


def test_setup_writes_endpoint_config_without_staged_data(
    su_bundle, su_site, tmp_path, monkeypatch
):
    _fake_install(monkeypatch, su_bundle, _report(su_bundle, data_readable=False))
    state = api.setup_site(su_bundle, su_site, tmp_path / "state")
    assert not (state / "identity_mapping.json").exists()
    config = yaml.safe_load((state / "endpoint_config.yaml").read_text())
    provider = config["engine"]["provider"]
    assert provider["type"] == "LocalProvider" and "account" not in provider
    assert str(state / "runtime/bin/activate") in provider["worker_init"]
    assert su_bundle.endpoint_name in config["display_name"]
    assert config["idle_heartbeats_hard"] == 5760
    assert set(api.prepared_files(su_bundle)) == {"runtime/bin/python", "endpoint_config.yaml"}
    # Data is required where a worker first reads it, not before the software exists.
    monkeypatch.setattr(api, "inspect_runtime", lambda *a: _report(su_bundle, data_readable=False))
    with pytest.raises(ValueError, match="data_dir must exist"):
        api.verify_site(state)


@pytest.fixture
def prepared(su_bundle, su_site, tmp_path, monkeypatch):
    _fake_install(monkeypatch, su_bundle)
    state = api.setup_site(su_bundle, su_site, tmp_path / "state")
    home = tmp_path / "home"
    monkeypatch.setattr(Path, "home", lambda: home)
    monkeypatch.delenv("GLOBUS_COMPUTE_USER_DIR", raising=False)
    monkeypatch.delenv("GLOBUS_COMPUTE_CLIENT_ID", raising=False)
    monkeypatch.delenv("GLOBUS_COMPUTE_CLIENT_SECRET", raising=False)
    # The real uid, so the credentials file's ownership check sees its true owner.
    uid = os.geteuid()
    if uid == 0:
        pytest.skip("single-user activation is exercised as an unprivileged account")
    monkeypatch.setattr(api.os, "geteuid", lambda: uid)
    import pwd

    monkeypatch.setattr(
        pwd, "getpwuid", lambda uid: SimpleNamespace(pw_name=su_bundle.service_account)
    )
    monkeypatch.setattr(api, "verify_site", lambda state, endpoint_id=None: {})
    return state, home


def _credentials(tmp_path, mode=0o600, text=None):
    path = tmp_path / "creds.env"
    path.write_text(
        text or f"GLOBUS_COMPUTE_CLIENT_ID={CLIENT}\nGLOBUS_COMPUTE_CLIENT_SECRET='s3cr3t'\n"
    )
    path.chmod(mode)
    return path


def _fake_endpoint(monkeypatch, home, endpoint_name, *, write_uuid=True):
    calls = []

    def run(args, **kwargs):
        calls.append((args, kwargs.get("env", {})))
        if args[1] == "configure":
            (home / ".globus_compute" / endpoint_name).mkdir(parents=True)
        if args[1] == "start" and write_uuid:
            info = home / ".globus_compute" / endpoint_name / "endpoint.json"
            info.write_text(json.dumps({"endpoint_id": ENDPOINT}))
        return SimpleNamespace(stdout="")

    monkeypatch.setattr(api, "run", run)
    monkeypatch.setattr(api.subprocess, "run", lambda args, **kwargs: calls.append((args, {})))
    return calls


def test_single_user_activation_starts_endpoint_as_service_account(
    su_bundle, prepared, tmp_path, monkeypatch
):
    state, home = prepared
    calls = _fake_endpoint(monkeypatch, home, su_bundle.endpoint_name)
    activated = api.activate_single_user(state, _credentials(tmp_path), wait=1)
    assert activated["endpoint_id"] == ENDPOINT
    (configure, env), (start, start_env) = calls
    assert configure[1:] == ["configure", "--multi-user", "false", su_bundle.endpoint_name]
    assert start[1:] == ["start", "--detach", su_bundle.endpoint_name]
    assert start_env["GLOBUS_COMPUTE_CLIENT_ID"] == CLIENT
    assert start_env["GLOBUS_COMPUTE_CLIENT_SECRET"] == "s3cr3t"  # quotes stripped
    endpoint = home / ".globus_compute" / su_bundle.endpoint_name
    assert (endpoint / "config.yaml").read_text() == (state / "endpoint_config.yaml").read_text()
    unit = activated["unit"].read_text()
    assert f"User={su_bundle.service_account}" in unit
    assert f"EnvironmentFile={tmp_path / 'creds.env'}" in unit
    assert "s3cr3t" not in unit and "--detach" not in unit
    # A second activation restarts the endpoint it manages instead of refusing it.
    calls.clear()
    api.activate_single_user(state, _credentials(tmp_path), wait=1)
    assert [c[0][1] for c in calls] == ["stop", "start"]


@pytest.mark.parametrize(
    "change, message",
    [
        ("root", "not root"),
        ("other_user", "this shell is someone"),
        ("no_credentials", "personal Globus login"),
        ("group_readable", "readable only by its owner"),
        ("bad_key", "only GLOBUS_COMPUTE_CLIENT_ID"),
        ("unmanaged", "unmanaged endpoint"),
    ],
)
def test_single_user_activation_refuses(
    su_bundle, prepared, tmp_path, monkeypatch, change, message
):
    state, home = prepared
    calls = _fake_endpoint(monkeypatch, home, su_bundle.endpoint_name)
    credentials = _credentials(tmp_path)
    if change == "root":
        monkeypatch.setattr(api.os, "geteuid", lambda: 0)
    elif change == "other_user":
        import pwd

        monkeypatch.setattr(pwd, "getpwuid", lambda uid: SimpleNamespace(pw_name="someone"))
    elif change == "no_credentials":
        credentials = None
    elif change == "group_readable":
        credentials = _credentials(tmp_path, mode=0o640)
    elif change == "bad_key":
        credentials = _credentials(tmp_path, text="AWS_SECRET=hunter2\n")
    elif change == "unmanaged":
        (home / ".globus_compute" / su_bundle.endpoint_name).mkdir(parents=True)
    with pytest.raises(ValueError, match=message) as caught:
        api.activate_single_user(state, credentials, wait=1)
    assert "hunter2" not in str(caught.value) and "s3cr3t" not in str(caught.value)
    assert calls == []


def test_multi_user_activation_refuses_single_user_bundle(prepared):
    state, _home = prepared
    with pytest.raises(ValueError, match="single-user endpoint"):
        api.activate_site(state)


def test_activate_cli_prints_uuid_and_systemd_steps(su_bundle, prepared, tmp_path, monkeypatch):
    state, home = prepared
    _fake_endpoint(monkeypatch, home, su_bundle.endpoint_name)
    monkeypatch.setattr(api, "load_state", lambda s: (su_bundle, None))
    result = CliRunner().invoke(
        main,
        [
            "site",
            "activate",
            "--state",
            str(state),
            "--credentials-file",
            str(_credentials(tmp_path)),
        ],
    )
    assert result.exit_code == 0, result.output
    assert f"Endpoint UUID: {ENDPOINT}" in result.output
    assert "systemctl daemon-reload" in result.output and "s3cr3t" not in result.output


# -- registration ------------------------------------------------------------


def _receipt(bundle, site, tmp_path):
    receipt = api.Receipt(
        bundle_digest=api.digest(bundle),
        experiment=bundle.experiment,
        client_id=bundle.client_id,
        endpoint_id=ENDPOINT,
        site=site,
        runtime=_report(bundle),
        state_dir=tmp_path / "state",
    )
    path = tmp_path / "registration.json"
    api.write_json(path, receipt.model_dump(mode="json"))
    return path


def _executor(monkeypatch, report):
    import globus_compute_sdk

    calls = []

    class Executor:
        def __init__(self, **kwargs):
            calls.append(kwargs)

        def submit(self, function, *args):
            return SimpleNamespace(result=lambda timeout: report)

        def shutdown(self, **kwargs):
            pass

    monkeypatch.setattr(globus_compute_sdk, "Executor", Executor)
    return calls


def test_single_user_registration_needs_a_route(su_bundle, su_site, tmp_path, monkeypatch):
    calls = _executor(monkeypatch, _report(su_bundle))
    with pytest.raises(ValueError, match="Pass --site-credentials"):
        api.register_site(su_bundle, _receipt(su_bundle, su_site, tmp_path), tmp_path / "o", 5)
    assert calls == []


def test_single_user_registration_probes_as_the_site_client(
    su_bundle, su_site, tmp_path, monkeypatch
):
    calls = _executor(monkeypatch, _report(su_bundle))
    client = object()
    credential = SimpleNamespace(expected_user=su_bundle.service_account, client=lambda: client)
    out = tmp_path / "ready.json"
    result = api.register_site(
        su_bundle, _receipt(su_bundle, su_site, tmp_path), out, 5, credential
    )
    assert result["status"] == "ready" and out.exists()
    assert calls == [{"endpoint_id": ENDPOINT, "client": client}]

    wrong = SimpleNamespace(expected_user="ubuntu", client=lambda: client)
    with pytest.raises(ValueError, match="expects user ubuntu"):
        api.register_site(
            su_bundle, _receipt(su_bundle, su_site, tmp_path), tmp_path / "x", 5, wrong
        )


# -- mixed credential routes -------------------------------------------------


def _routes(tmp_path, entries):
    path = tmp_path / "routes.json"
    path.write_text(json.dumps({"version": 1, "endpoints": entries}))
    return path


def test_coordinator_login_route_is_explicit_and_needs_no_secret(tmp_path, monkeypatch):
    mep, single = str(uuid4()), str(uuid4())
    monkeypatch.setenv("SITE_SECRET", "x")
    path = _routes(
        tmp_path,
        {
            mep: {"coordinator_login": True, "expected_user": "mfjoel01"},
            single: {
                "client_id": CLIENT,
                "client_secret_env": "SITE_SECRET",
                "expected_user": SERVICE,
            },
        },
    )
    routes = load_site_credentials(path, [mep, single])
    assert routes[mep].coordinator_login and routes[mep].client_id is None
    assert not routes[single].coordinator_login
    # Still no fallback: an endpoint missing from the manifest is refused.
    with pytest.raises(ValueError, match="Missing site credentials"):
        load_site_credentials(path, [mep, str(uuid4())])


@pytest.mark.parametrize(
    "entry",
    [
        {"coordinator_login": False, "expected_user": "u"},
        {"coordinator_login": "yes", "expected_user": "u"},
        {"coordinator_login": True},
        {"coordinator_login": True, "expected_user": "u", "client_id": CLIENT},
    ],
)
def test_malformed_coordinator_login_routes_are_refused(tmp_path, entry):
    endpoint = str(uuid4())
    with pytest.raises(ValueError, match="coordinator_login: true"):
        load_site_credentials(_routes(tmp_path, {endpoint: entry}), [endpoint])


def test_coordinator_login_uses_the_default_client(monkeypatch):
    import globus_compute_sdk

    made = []
    monkeypatch.setattr(globus_compute_sdk, "Client", lambda **kwargs: made.append(kwargs) or "c")
    assert SiteCredential(None, None, "u", coordinator_login=True).client() == "c"
    assert set(made[0]) == {"code_serialization_strategy"}


def test_isolation_trial_refuses_coordinator_login(tmp_path, monkeypatch):
    from appfl_bio_suite.core.site_trial import check_sites

    monkeypatch.setenv("SITE_SECRET", "x")
    path = _routes(
        tmp_path,
        {
            str(uuid4()): {"coordinator_login": True, "expected_user": "u"},
            str(uuid4()): {
                "client_id": CLIENT,
                "client_secret_env": "SITE_SECRET",
                "expected_user": "u",
            },
        },
    )
    with pytest.raises(ValueError, match="no site identity to isolate"):
        check_sites(path)


def test_register_cli_routes_through_the_manifest(
    su_dir, su_bundle, su_site, tmp_path, monkeypatch
):
    seen = []
    monkeypatch.setattr(
        api,
        "register_site",
        lambda spec, receipt, out, timeout, credential=None: (
            seen.append(credential) or {"client_id": spec.client_id, "endpoint_id": ENDPOINT}
        ),
    )
    monkeypatch.setattr(api, "updated_federation", lambda *a: {"schema_version": 1})
    monkeypatch.setenv("SITE_SECRET", "x")
    routes = _routes(
        tmp_path,
        {
            ENDPOINT: {
                "client_id": CLIENT,
                "client_secret_env": "SITE_SECRET",
                "expected_user": SERVICE,
            }
        },
    )
    federation = tmp_path / "federation.yaml"
    result = CliRunner().invoke(
        main,
        [
            "experiment", "register", str(_receipt(su_bundle, su_site, tmp_path)),
            "--bundle", str(su_dir), "--federation", str(federation),
            "--updated-federation", str(tmp_path / "next.yaml"),
            "--out", str(tmp_path / "ready.json"), "--site-credentials", str(routes),
        ],
    )  # fmt: skip
    assert result.exit_code == 0, result.output
    assert seen[0].client_id == CLIENT


# -- runs ---------------------------------------------------------------------


def test_spooling_is_enabled_only_for_globus_compute_fine_mapping():
    from appfl_bio_suite.core.launch import enable_result_spool
    from appfl_bio_suite.core.result_spool import INLINE_LIMIT_BYTES

    clients = [{"train_configs": {}}, {"train_configs": {"result_spool_bytes": 5}}]
    enable_result_spool("fine-mapping", clients)
    assert [c["train_configs"]["result_spool_bytes"] for c in clients] == [INLINE_LIMIT_BYTES, 5]
    gwas = [{"train_configs": {}}]
    enable_result_spool("gwas", gwas)
    assert gwas == [{"train_configs": {}}]


def test_driver_fetches_spooled_results_from_the_right_endpoint(tmp_path, monkeypatch):
    import torch

    from appfl_bio_suite.core import result_spool
    from appfl_bio_suite.core.drivers.globus_compute_driver import _fetch_if_spooled

    inline = {"manifest": torch.zeros(3, dtype=torch.uint8)}
    log = SimpleNamespace(info=lambda m: None, warning=lambda m: None)
    assert _fetch_if_spooled(inline, ENDPOINT, "site", None, {}, tmp_path, log) is inline

    info = {
        "format": result_spool.SPOOL_FORMAT,
        "path": "/work/spool/site-1.pt",
        "size": 10,
        "sha256": "0" * 64,
        "client_id": "site",
    }
    blob = torch.frombuffer(bytearray(json.dumps(info).encode()), dtype=torch.uint8)
    made, removed = [], []
    executor = SimpleNamespace(
        submit=lambda fn, *a: removed.append(a) or SimpleNamespace(result=lambda timeout: True)
    )
    monkeypatch.setattr(
        result_spool, "spool_executor", lambda e, c: made.append((e, c)) or executor
    )
    monkeypatch.setattr(result_spool, "fetch_spooled", lambda i, submit, d, log: {"ok": i["path"]})
    route = SimpleNamespace(client_id=CLIENT)
    executors = {}
    out = _fetch_if_spooled(
        {"result_spool": blob}, ENDPOINT, "site", {ENDPOINT: route}, executors, tmp_path, log
    )
    assert out == {"ok": "/work/spool/site-1.pt"}
    assert made == [(ENDPOINT, route)] and removed == [("/work/spool/site-1.pt",)]


# -- provenance ---------------------------------------------------------------


def test_installed_wheel_never_reports_an_enclosing_checkout(tmp_path):
    from appfl_bio_suite.core.simulation import _git_commit

    repo = tmp_path / "someone-elses-repo"
    site_packages = repo / ".venv/lib/python3.12"
    site_packages.mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    assert _git_commit(site_packages) is None
    (site_packages / "pyproject.toml").write_text("")  # still not the repository root
    assert _git_commit(site_packages) is None
    assert _git_commit(repo_root()) not in (None, "")


def test_credentials_file_mode_check_uses_owner_bits_only(tmp_path):
    path = _credentials(tmp_path, mode=0o600)
    assert api.read_credentials_file(path)["GLOBUS_COMPUTE_CLIENT_ID"] == CLIENT
    path.chmod(0o604)
    assert stat.S_IMODE(os.stat(path).st_mode) & 0o077
    with pytest.raises(ValueError, match="readable only by its owner"):
        api.read_credentials_file(path)
