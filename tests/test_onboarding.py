"""Onboarding contracts: reject drift before mutation, and prove remote readiness."""

import json
import zipfile
from pathlib import Path
from types import SimpleNamespace
from uuid import uuid4

import pytest
import yaml
from click.testing import CliRunner

from appfl_bio_suite.cli import main
from appfl_bio_suite.core import onboarding as api
from appfl_bio_suite.core.config import load_federation
from appfl_bio_suite.core.experiments import repo_root
from appfl_bio_suite.core.partner import generate_bundle


@pytest.fixture
def bundle_dir(tmp_path):
    return generate_bundle(
        load_federation(repo_root() / "federation.yaml.example"),
        "gwas",
        "site-north",
        tmp_path / "bundles",
    )


@pytest.fixture
def bundle(bundle_dir):
    return api.read_bundle(bundle_dir)


@pytest.fixture
def local(tmp_path, bundle):
    data = tmp_path / "data"
    data.mkdir()
    values = {**bundle.defaults, "data_dir": str(data), "work_dir": str(tmp_path / "work")}
    path = tmp_path / "site.yaml"
    path.write_text(yaml.safe_dump(values))
    return path


def good_report(bundle):
    return {
        "bundle_digest": api.digest(bundle),
        "python": "3.12.11",
        "versions": {**bundle.versions, "appfl-bio-suite": bundle.suite_version},
        "user": bundle.service_account,
        "data_readable": True,
        "work_writable": True,
    }


def make_receipt(bundle, site):
    return api.Receipt(
        bundle_digest=api.digest(bundle),
        experiment=bundle.experiment,
        client_id=bundle.client_id,
        endpoint_id=uuid4(),
        site=site,
        runtime=good_report(bundle),
        state_dir=site.work_dir.parent / "state",
    )


def test_zip_and_directory_have_same_revision(bundle_dir, bundle, tmp_path):
    archive = tmp_path / "site.zip"
    with zipfile.ZipFile(archive, "w") as stream:
        stream.write(bundle_dir / "onboarding.json", "onboarding.json")
        stream.writestr("../../never-extract", "bad")
    assert api.digest(api.read_bundle(archive)) == api.digest(bundle)
    assert not (tmp_path / "never-extract").exists()


def test_duplicate_specs_refused(bundle_dir, tmp_path):
    archive = tmp_path / "site.zip"
    with zipfile.ZipFile(archive, "w") as stream:
        for name in ["one/onboarding.json", "two/onboarding.json"]:
            stream.write(bundle_dir / "onboarding.json", name)
    with pytest.raises(ValueError, match="exactly one"):
        api.read_bundle(archive)


@pytest.mark.parametrize(
    "change", [{"schema_version": 2}, {"suite_version": "99.0"}, {"unexpected": True}]
)
def test_invalid_bundle_refused(bundle_dir, change):
    path = bundle_dir / "onboarding.json"
    raw = json.loads(path.read_text())
    raw.update(change)
    path.write_text(json.dumps(raw))
    with pytest.raises(ValueError):
        api.read_bundle(path)


@pytest.mark.parametrize(
    "change",
    [
        {"account": ""},
        {"account": "a\nmalicious"},
        {"data_dir": "relative"},
        {"service_account": "root"},
        {"typo": "value"},
    ],
)
def test_local_allowlist(bundle, local, change):
    raw = yaml.safe_load(local.read_text())
    raw.update(change)
    local.write_text(yaml.safe_dump(raw))
    with pytest.raises(ValueError):
        api.resolve_site(bundle, local)


def test_scheduler_mismatch(bundle, local):
    raw = yaml.safe_load(local.read_text())
    other = "queue" if bundle.scheduler == "slurm" else "partition"
    raw[other] = "wrong"
    local.write_text(yaml.safe_dump(raw))
    with pytest.raises(ValueError, match="do not set"):
        api.resolve_site(bundle, local)


def test_setup_resumes_and_rejects_changed_inputs(bundle, local, tmp_path, monkeypatch):
    site = api.resolve_site(bundle, local)
    calls = []

    def fake_run(args, **kwargs):
        calls.append(args)
        if "venv" in args:
            target = Path(args[-1]) / "bin/python"
            target.parent.mkdir(parents=True)
            target.touch()
        return SimpleNamespace(stdout="3.12\n")

    monkeypatch.setattr(api, "run", fake_run)
    monkeypatch.setattr(api, "inspect_runtime", lambda *a: good_report(bundle))
    state = tmp_path / "state"
    api.setup_site(bundle, site, state)
    api.setup_site(bundle, site, state)
    assert sum("venv" in call for call in calls) == 1
    assert (state / "identity_mapping.json").exists()
    config = yaml.safe_load((state / "user_config_template.yaml.j2").read_text())
    assert config["engine"]["provider"]["account"] == site.account
    assert "REPLACE_WITH" not in (state / "user_config_template.yaml.j2").read_text()
    calls.clear()
    with pytest.raises(ValueError, match="different inputs"):
        api.setup_site(bundle, site.model_copy(update={"account": "new"}), state)
    assert calls == []


def test_setup_does_not_claim_unrelated_directory(bundle, local, tmp_path):
    state = tmp_path / "other"
    state.mkdir()
    (state / "important").write_text("keep")
    with pytest.raises(ValueError, match="not empty"):
        api.setup_site(bundle, api.resolve_site(bundle, local), state)
    assert (state / "important").read_text() == "keep"


@pytest.mark.parametrize(
    "change",
    [
        {"user": "root"},
        {"python": "3.11.9"},
        {"data_readable": False},
        {"work_writable": False},
        {"versions": {}},
    ],
)
def test_runtime_failures_block_readiness(bundle, change):
    with pytest.raises(ValueError):
        api.check_runtime(bundle, {**good_report(bundle), **change})


def test_cpu_torch_local_version_accepted(bundle):
    report = good_report(bundle)
    report["versions"]["torch"] += "+cpu"
    api.check_runtime(bundle, report)


def fake_executor(monkeypatch, report=None, error=None):
    import globus_compute_sdk

    calls = []

    class Executor:
        def __init__(self, **kwargs):
            calls.append(kwargs)

        def submit(self, function, *args):
            calls.append((function, args))
            return self

        def result(self, timeout):
            if error:
                raise error
            return report

        def shutdown(self, **kwargs):
            calls.append(kwargs)

    monkeypatch.setattr(globus_compute_sdk, "Executor", Executor)
    return calls


def test_registration_runs_real_probe_contract(bundle, local, tmp_path, monkeypatch):
    site = api.resolve_site(bundle, local)
    receipt = make_receipt(bundle, site)
    receipt_path = tmp_path / "receipt.json"
    api.write_json(receipt_path, receipt.model_dump(mode="json"))
    calls = fake_executor(monkeypatch, good_report(bundle))
    out = tmp_path / "ready.json"
    assert api.register_site(bundle, receipt_path, out, 5)["status"] == "ready"
    assert calls[1][0] is api.runtime_report
    assert calls[-1] == {"wait": False, "cancel_futures": True}
    assert out.exists()


def test_stale_receipt_never_submits(bundle, local, tmp_path, monkeypatch):
    receipt = make_receipt(bundle, api.resolve_site(bundle, local))
    receipt.bundle_digest = "stale"
    path = tmp_path / "receipt.json"
    api.write_json(path, receipt.model_dump(mode="json"))
    calls = fake_executor(monkeypatch, good_report(bundle))
    with pytest.raises(ValueError, match="revision"):
        api.register_site(bundle, path, tmp_path / "ready.json", 5)
    assert not calls


@pytest.mark.parametrize("failure", ["timeout", "wrong_user"])
def test_remote_failure_never_writes_ready(bundle, local, tmp_path, monkeypatch, failure):
    receipt = make_receipt(bundle, api.resolve_site(bundle, local))
    path = tmp_path / "receipt.json"
    api.write_json(path, receipt.model_dump(mode="json"))
    report = {**good_report(bundle), "user": "wrong"}
    calls = fake_executor(
        monkeypatch, report, TimeoutError("queued") if failure == "timeout" else None
    )
    out = tmp_path / "ready.json"
    with pytest.raises((ValueError, TimeoutError)):
        api.register_site(bundle, path, out, 1)
    assert not out.exists()
    assert calls[-1]["wait"] is False


def test_registration_syncs_local_paths(bundle, local):
    receipt = make_receipt(bundle, api.resolve_site(bundle, local))
    raw = api.updated_federation(bundle, receipt, repo_root() / "federation.yaml.example")
    updated = next(
        x for x in raw["experiments"]["gwas"]["sites"] if x["client_id"] == bundle.client_id
    )
    assert updated["endpoint_uuid"] == str(receipt.endpoint_id)
    assert updated["data_dir"] == str(receipt.site.data_dir)
    assert updated["output_dir"] == str(receipt.site.work_dir)


def test_package_without_endpoint_ids(tmp_path):
    from appfl_bio_suite.core.launch import build_client_configs

    source = repo_root() / "federation.yaml.example"
    raw = yaml.safe_load(source.read_text())
    for entry in raw["experiments"]["gwas"]["sites"]:
        entry.pop("endpoint_uuid")
    path = tmp_path / "federation.yaml"
    path.write_text(yaml.safe_dump(raw))
    fed = load_federation(path)
    with pytest.raises(ValueError, match="register endpoints"):
        build_client_configs(fed, "gwas")
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "experiment",
            "package",
            "gwas",
            "--federation",
            str(path),
            "--out",
            str(tmp_path / "bundles"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert len(list((tmp_path / "bundles").glob("*.zip"))) == 3


def test_cli_init_preserves_local_file(bundle_dir, tmp_path):
    local = tmp_path / "site.yaml"
    local.write_text("keep")
    result = CliRunner().invoke(main, ["site", "init", str(bundle_dir), "--out", str(local)])
    assert result.exit_code != 0
    assert local.read_text() == "keep"


def test_wheel_transport_and_corruption(bundle_dir, tmp_path):
    wheel = tmp_path / "appfl_bio_suite-0.2.1-py3-none-any.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr(
            "appfl_bio_suite-0.2.1.dist-info/METADATA", "Name: appfl-bio-suite\nVersion: 0.2.1\n"
        )
    api.attach_wheel(bundle_dir, wheel)
    bundle = api.read_bundle(bundle_dir)
    assert api.read_wheel(bundle_dir, bundle) == wheel.read_bytes()
    (bundle_dir / wheel.name).write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="checksum"):
        api.read_wheel(bundle_dir, bundle)


def test_wrong_wheel_rejected_before_copy(bundle_dir, tmp_path):
    wheel = tmp_path / "appfl_bio_suite-99-py3-none-any.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr(
            "appfl_bio_suite-99.dist-info/METADATA", "Name: appfl-bio-suite\nVersion: 99\n"
        )
    with pytest.raises(ValueError, match="does not match"):
        api.attach_wheel(bundle_dir, wheel)
    assert not (bundle_dir / wheel.name).exists()


def test_changed_experiment_refuses_registration(bundle, local, tmp_path):
    receipt = make_receipt(bundle, api.resolve_site(bundle, local))
    raw = yaml.safe_load((repo_root() / "federation.yaml.example").read_text())
    raw["experiments"]["gwas"]["variant_scaling"] = 0.5
    source = tmp_path / "changed.yaml"
    source.write_text(yaml.safe_dump(raw))
    with pytest.raises(ValueError, match="settings changed"):
        api.updated_federation(bundle, receipt, source)


def make_state(tmp_path, bundle, site):
    state = tmp_path / "state"
    api.write_json(
        state / "state.json",
        {
            "bundle_digest": api.digest(bundle),
            "site": site.model_dump(mode="json"),
            "install_source": bundle.install_spec,
        },
    )
    api.write_json(state / "onboarding.json", bundle.model_dump())
    api.write_json(state / "identity_mapping.json", [{"test": "mapping"}])
    (state / "user_config_template.yaml.j2").write_text(
        api.worker_config(bundle, site, state / "runtime")
    )
    return state


def test_activation_refuses_unprivileged_user(bundle, local, tmp_path, monkeypatch):
    state = make_state(tmp_path, bundle, api.resolve_site(bundle, local))
    monkeypatch.setattr(api.os, "geteuid", lambda: 1000)
    with pytest.raises(ValueError, match="requires root"):
        api.activate_site(state)


def test_activation_preserves_generated_mapping_path_and_endpoint(
    bundle, local, tmp_path, monkeypatch
):
    import pwd

    state = make_state(tmp_path, bundle, api.resolve_site(bundle, local))
    monkeypatch.setattr(api.os, "geteuid", lambda: 0)
    monkeypatch.setattr(pwd, "getpwnam", lambda name: object())
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    endpoint = tmp_path / ".globus_compute" / bundle.endpoint_name
    calls = []

    def run(args, **kwargs):
        calls.append(args)
        if "configure" in args:
            endpoint.mkdir(parents=True)
            (endpoint / "config.yaml").write_text(
                "identity_mapping_config_path: generated-mapping.json\n"
            )
            (endpoint / "endpoint.json").write_text('{"endpoint_id": "unchanged"}')
        return SimpleNamespace(stdout="")

    monkeypatch.setattr(api, "run", run)
    api.activate_site(state)
    api.activate_site(state)
    assert sum("configure" in call for call in calls) == 1
    assert (endpoint / "generated-mapping.json").read_text() == (
        state / "identity_mapping.json"
    ).read_text()
    assert json.loads((endpoint / "endpoint.json").read_text())["endpoint_id"] == "unchanged"
    assert calls[0][0] == "runuser"


def test_activation_refuses_existing_endpoint(bundle, local, tmp_path, monkeypatch):
    import pwd

    state = make_state(tmp_path, bundle, api.resolve_site(bundle, local))
    monkeypatch.setattr(api.os, "geteuid", lambda: 0)
    monkeypatch.setattr(pwd, "getpwnam", lambda name: object())
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    calls = []
    monkeypatch.setattr(api, "run", lambda args, **kwargs: calls.append(args))
    endpoint = tmp_path / ".globus_compute" / bundle.endpoint_name
    endpoint.mkdir(parents=True)
    (endpoint / "config.yaml").write_text("keep")
    with pytest.raises(ValueError, match="unmanaged endpoint"):
        api.activate_site(state)
    assert (endpoint / "config.yaml").read_text() == "keep"
    assert not any("configure" in call for call in calls)


def test_no_receipt_on_failed_local_verification(bundle, local, tmp_path, monkeypatch):
    state = make_state(tmp_path, bundle, api.resolve_site(bundle, local))
    monkeypatch.setattr(
        api, "inspect_runtime", lambda *args: {**good_report(bundle), "user": "bad"}
    )
    with pytest.raises(ValueError, match="worker must run"):
        api.verify_site(state, str(uuid4()))
    assert not (state / "registration.json").exists()


def test_register_cli_writes_ready_and_updated_config(
    bundle_dir, bundle, local, tmp_path, monkeypatch
):
    receipt = make_receipt(bundle, api.resolve_site(bundle, local))
    receipt_path = tmp_path / "registration.json"
    api.write_json(receipt_path, receipt.model_dump(mode="json"))
    fake_executor(monkeypatch, good_report(bundle))
    out = tmp_path / "ready.json"
    updated = tmp_path / "registered.yaml"
    result = CliRunner().invoke(
        main,
        [
            "experiment",
            "register",
            str(receipt_path),
            "--bundle",
            str(bundle_dir),
            "--federation",
            str(repo_root() / "federation.yaml.example"),
            "--out",
            str(out),
            "--updated-federation",
            str(updated),
        ],
    )
    assert result.exit_code == 0, result.output
    entry = load_federation(updated).experiment_site("gwas", "site-north")
    assert entry.endpoint_uuid == str(receipt.endpoint_id)
    assert entry.output_dir == str(receipt.site.work_dir)
    # Registering another partner must not invalidate already-issued bundles.
    assert api.experiment_digest(load_federation(updated), "gwas") == bundle.experiment_digest


def test_membership_change_invalidates_bundle(bundle, local, tmp_path):
    receipt = make_receipt(bundle, api.resolve_site(bundle, local))
    raw = yaml.safe_load((repo_root() / "federation.yaml.example").read_text())
    raw["experiments"]["gwas"]["sites"].pop()
    source = tmp_path / "changed.yaml"
    source.write_text(yaml.safe_dump(raw))
    with pytest.raises(ValueError, match="settings changed"):
        api.updated_federation(bundle, receipt, source)


def test_probe_executes_in_subprocess_and_reads_installed_revision(bundle, local, tmp_path):
    import sys

    site = api.resolve_site(bundle, local)
    site.work_dir.mkdir()
    state = make_state(tmp_path, bundle, site)
    python = state / "runtime/bin/python"
    python.parent.mkdir(parents=True)
    python.symlink_to(sys.executable)
    report = api.inspect_runtime(python, bundle, site)
    assert report["bundle_digest"] == api.digest(bundle)
    assert report["data_readable"] is True
    assert report["work_writable"] is True
    (state / "onboarding.json").unlink()
    report = api.inspect_runtime(python, bundle, site)
    assert report["bundle_digest"] == "missing"


def test_remote_bundle_drift_blocks_readiness(bundle, local, tmp_path, monkeypatch):
    receipt = make_receipt(bundle, api.resolve_site(bundle, local))
    path = tmp_path / "receipt.json"
    api.write_json(path, receipt.model_dump(mode="json"))
    fake_executor(monkeypatch, {**good_report(bundle), "bundle_digest": "changed"})
    out = tmp_path / "ready.json"
    with pytest.raises(ValueError, match="bundle revision"):
        api.register_site(bundle, path, out, 5)
    assert not out.exists()


@pytest.fixture
def join_setup(monkeypatch):
    calls = []

    def setup(bundle, site, state, **kwargs):
        calls.append((bundle, site, state, kwargs))
        state.mkdir(parents=True, exist_ok=True)
        api.write_json(state / "onboarding.json", bundle.model_dump())
        api.write_json(
            state / "state.json",
            {
                "bundle_digest": api.digest(bundle),
                "site": site.model_dump(mode="json"),
                "install_source": kwargs.get("install_source")
                or bundle.wheel_sha256
                or bundle.install_spec,
            },
        )
        for name in ("runtime/bin/python", "identity_mapping.json", "user_config_template.yaml.j2"):
            target = state / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.touch()
        return state

    monkeypatch.setattr(api, "setup_site", setup)
    return calls


def test_join_uses_prepared_local_file_without_prompts(bundle_dir, local, tmp_path, join_setup):
    profile = tmp_path / "cluster.yaml"
    result = CliRunner().invoke(
        main, ["join", str(bundle_dir), "--local", str(local), "--profile", str(profile)]
    )
    assert result.exit_code == 0, result.output
    assert "Scheduler account:" not in result.output
    assert "site activate" in result.output
    assert "--endpoint-id UUID" in result.output
    assert "registration.json" not in result.output
    values = yaml.safe_load(profile.read_text())
    assert "data_dir" not in values and "work_dir" not in values
    assert len(join_setup) == 1
    assert (join_setup[0][2] / "site.yaml").exists()


def test_join_guided_prompts_and_reprompts_bad_path(bundle_dir, tmp_path, join_setup):
    answers = f"my_account\nmy_queue\nrelative\n{tmp_path / 'data'}\n{tmp_path / 'work'}\n"
    result = CliRunner().invoke(main, ["join", str(bundle_dir), "--no-profile"], input=answers)
    assert result.exit_code == 0, result.output
    assert "Invalid data_dir" in result.output
    assert join_setup[0][1].account == "my_account"
    assert join_setup[0][1].data_dir == tmp_path / "data"


def test_join_reuses_cluster_but_prompts_for_experiment_paths(
    bundle_dir, bundle, tmp_path, join_setup
):
    profile = tmp_path / "cluster.yaml"
    slot = "partition" if bundle.scheduler == "slurm" else "queue"
    profile.write_text(
        yaml.safe_dump({"scheduler": bundle.scheduler, "account": "remembered", slot: "cpu"})
    )
    result = CliRunner().invoke(
        main,
        ["join", str(bundle_dir), "--profile", str(profile)],
        input=f"{tmp_path / 'data'}\n{tmp_path / 'work'}\n",
    )
    assert result.exit_code == 0, result.output
    assert "Scheduler account:" not in result.output
    assert "Data directory (absolute path)" in result.output
    assert join_setup[0][1].account == "remembered"


def test_join_rejects_wrong_cluster_profile(bundle_dir, bundle, tmp_path, join_setup):
    profile = tmp_path / "cluster.yaml"
    profile.write_text(
        yaml.safe_dump(
            {"scheduler": "pbspro" if bundle.scheduler == "slurm" else "slurm", "account": "other"}
        )
    )
    result = CliRunner().invoke(main, ["join", str(bundle_dir), "--profile", str(profile)])
    assert result.exit_code != 0
    assert "different scheduler" in result.output
    assert not join_setup


def test_join_noninteractive_missing_fields_fails_before_install(
    bundle_dir, local, tmp_path, join_setup
):
    values = yaml.safe_load(local.read_text())
    values["account"] = ""
    local.write_text(yaml.safe_dump(values))
    result = CliRunner().invoke(
        main, ["join", str(bundle_dir), "--local", str(local), "--non-interactive", "--no-profile"]
    )
    assert result.exit_code != 0
    assert not join_setup


def test_join_resume_registration_skips_install_and_prompts(
    bundle_dir, bundle, local, tmp_path, monkeypatch, join_setup
):
    state = tmp_path / "state with spaces"
    runner = CliRunner()
    first = runner.invoke(
        main,
        ["join", str(bundle_dir), "--local", str(local), "--state", str(state), "--no-profile"],
    )
    assert first.exit_code == 0, first.output
    verified = []
    monkeypatch.setattr(api, "verify_site", lambda *args: verified.append(args))
    endpoint_id = str(uuid4())
    result = runner.invoke(
        main,
        [
            "join",
            str(bundle_dir),
            "--state",
            str(state),
            "--endpoint-id",
            endpoint_id,
            "--no-profile",
        ],
    )
    assert result.exit_code == 0, result.output
    assert len(join_setup) == 1
    assert verified == [(state, endpoint_id)]
    assert "registration.json" in result.output
    assert "remote probe is still required" in result.output
    assert f"'{state}'" in first.output


def test_join_rejects_state_from_different_bundle(bundle_dir, bundle, local, tmp_path, join_setup):
    state = make_state(tmp_path, bundle, api.resolve_site(bundle, local))
    raw = json.loads((bundle_dir / "onboarding.json").read_text())
    raw["client_id"] = "AnotherSite"
    (bundle_dir / "onboarding.json").write_text(json.dumps(raw))
    result = CliRunner().invoke(
        main, ["join", str(bundle_dir), "--state", str(state), "--no-profile"]
    )
    assert result.exit_code != 0
    assert "another bundle revision" in result.output
    assert not join_setup


def test_join_does_not_enable_unfinished_single_user_mode(bundle_dir, join_setup):
    result = CliRunner().invoke(main, ["join", str(bundle_dir), "--endpoint-mode", "single-user"])
    assert result.exit_code != 0
    assert "Invalid value" in result.output
    assert not join_setup


def test_join_cancel_does_not_install_or_save_profile(
    bundle_dir, tmp_path, join_setup, monkeypatch
):
    import click

    def abort(*args, **kwargs):
        raise click.Abort()

    monkeypatch.setattr(click, "prompt", abort)
    profile = tmp_path / "cluster.yaml"
    result = CliRunner().invoke(
        main, ["join", str(bundle_dir), "--profile", str(profile)], input=""
    )
    assert result.exit_code != 0
    assert not join_setup and not profile.exists()
