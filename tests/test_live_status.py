"""Live endpoint status: what the prober publishes, and what it must never say.

Two properties carry the weight. The status file is served from a public host, so no
endpoint UUID may reach it -- asserted against the serialized bytes. And a probe that
could not be made is ``unknown``, never ``offline``: one expired credential must not
become a statement that every partner is down.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from appfl_bio_suite.cli import main
from appfl_bio_suite.core.config import load_federation
from appfl_bio_suite.core.experiments import repo_root
from appfl_bio_suite.core.live_status import (
    STATUS_FILENAME,
    LiveStatusError,
    Route,
    load_probe_list,
    load_targets,
    missing_secrets,
    probe_targets,
    run,
    targets_from_federation,
)

UP = "11111111-aaaa-4bbb-8ccc-000000000001"
DOWN = "11111111-aaaa-4bbb-8ccc-000000000002"
DENIED = "11111111-aaaa-4bbb-8ccc-000000000003"


class FakeClient:
    """Globus Compute, reduced to the one call the prober makes."""

    def __init__(self, broken: bool = False):
        self.broken = broken
        self.calls: list[str] = []

    def get_endpoint_status(self, uuid: str) -> dict:
        self.calls.append(uuid)
        if self.broken or uuid == DENIED:
            raise RuntimeError(f"403 ENDPOINT_ACCESS_FORBIDDEN {uuid}")
        return {"status": "online" if uuid == UP else "offline"}


def test_site_states():
    targets = {
        "up": [UP],
        "down": [DOWN],
        "half": [UP, DOWN],
        "denied": [DENIED],
        "down-and-denied": [DOWN, DENIED],
    }
    sites = probe_targets(targets, client=FakeClient())["sites"]
    assert {site: entry["state"] for site, entry in sites.items()} == {
        "up": "online",
        "down": "offline",
        "half": "degraded",
        "denied": "unknown",
        "down-and-denied": "offline",
    }
    assert sites["half"] == {"state": "degraded", "endpoints": 2, "online": 1}


def test_status_file_never_carries_a_uuid(tmp_path):
    out = tmp_path / STATUS_FILENAME
    run({"a": [UP], "b": [DENIED]}, out, once=True, client_factory=FakeClient)
    text = out.read_text(encoding="utf-8")
    for uuid in (UP, DENIED):
        assert uuid not in text and uuid[:8] not in text
    assert json.loads(text)["sites"]["a"]["state"] == "online"
    assert not list(tmp_path.glob(".*tmp")), "the scratch file was renamed, not left behind"


def test_a_broken_prober_reports_unknown_not_offline(tmp_path):
    out = tmp_path / STATUS_FILENAME

    def no_client():
        raise RuntimeError("invalid_client")

    run({"a": [UP], "b": [DOWN]}, out, once=True, client_factory=no_client)
    states = {entry["state"] for entry in json.loads(out.read_text())["sites"].values()}
    assert states == {"unknown"}


def test_loop_survives_failures_and_replaces_a_dead_client(tmp_path):
    out = tmp_path / STATUS_FILENAME
    made: list[FakeClient] = []
    seen: list[str] = []

    def factory():
        made.append(FakeClient(broken=len(made) == 0))
        return made[-1]

    class Stop(Exception):
        pass

    def sleep(_seconds: float) -> None:
        if len(seen) == 2:
            raise Stop

    with pytest.raises(Stop):
        run(
            {"a": [UP]},
            out,
            interval=60,
            client_factory=factory,
            report=lambda status: seen.append(status["sites"]["a"]["state"]),
            sleep=sleep,
        )
    assert seen == ["unknown", "online"]
    assert len(made) == 2, "a client whose every request failed is not reused"


OWNED = "11111111-aaaa-4bbb-8ccc-000000000004"
SITE_CLIENT = "22222222-aaaa-4bbb-8ccc-000000000001"


class OwnerOnly(FakeClient):
    """A single-user endpoint: its status answers only the service account that owns it."""

    def __init__(self, identity: str | None):
        super().__init__()
        self.identity = identity

    def get_endpoint_status(self, uuid: str) -> dict:
        if uuid != OWNED:
            return super().get_endpoint_status(uuid)
        if self.identity != SITE_CLIENT:
            raise RuntimeError("403 ENDPOINT_ACCESS_FORBIDDEN")
        return {"status": "online"}


def test_a_route_asks_as_the_owning_account(tmp_path):
    out = tmp_path / STATUS_FILENAME
    made: list[str | None] = []

    def default():
        made.append(None)
        return OwnerOnly(None)

    def routed(route: Route):
        made.append(route.client_id)
        return OwnerOnly(route.client_id)

    targets = {"anl": [UP], "mbzuai": [OWNED]}
    run(targets, out, once=True, client_factory=default)
    assert json.loads(out.read_text())["sites"]["mbzuai"]["state"] == "unknown"

    made.clear()
    routes = {OWNED: Route(SITE_CLIENT, "SITE_SECRET")}
    run(targets, out, routes=routes, once=True, client_factory=default, route_factory=routed)
    sites = json.loads(out.read_text())["sites"]
    assert sites["mbzuai"]["state"] == "online" and sites["anl"]["state"] == "online"
    assert len(made) == 2 and set(made) == {None, SITE_CLIENT}, "one client per identity"


def test_a_route_that_cannot_authenticate_is_unknown_for_its_site_only(tmp_path):
    out = tmp_path / STATUS_FILENAME

    def broken(route: Route):
        raise ValueError(f"Set coordinator environment variable {route.client_secret_env}.")

    run(
        {"anl": [UP], "mbzuai": [OWNED]},
        out,
        routes={OWNED: Route(SITE_CLIENT, "SITE_SECRET")},
        once=True,
        client_factory=FakeClient,
        route_factory=broken,
    )
    sites = json.loads(out.read_text())["sites"]
    assert sites == {
        "anl": {"state": "online", "endpoints": 1, "online": 1},
        "mbzuai": {"state": "unknown", "endpoints": 1, "online": 0},
    }


def test_probe_list_routes(tmp_path):
    path = tmp_path / "targets.json"
    route = {"client_id": SITE_CLIENT, "client_secret_env": "SITE_SECRET"}
    path.write_text(json.dumps({"sites": {"mbzuai": [OWNED]}, "routes": {OWNED: route}}))
    targets, routes = load_probe_list(path)
    assert targets == {"mbzuai": [OWNED]}
    assert routes == {OWNED: Route(SITE_CLIENT, "SITE_SECRET")}
    assert missing_secrets(routes, {}) == ["SITE_SECRET"]
    assert missing_secrets(routes, {"SITE_SECRET": "x"}) == []

    for bad in (
        {**route, "client_secret": "s3"},  # the secret itself, inline
        {"client_id": "nope", "client_secret_env": "SITE_SECRET"},
        {"client_id": SITE_CLIENT, "client_secret_env": "not a name"},
    ):
        path.write_text(json.dumps({"sites": {"mbzuai": [OWNED]}, "routes": {OWNED: bad}}))
        with pytest.raises(LiveStatusError):
            load_probe_list(path)
    path.write_text(json.dumps({"sites": {"anl": [UP]}, "routes": {OWNED: route}}))
    with pytest.raises(LiveStatusError, match="not in `sites:`"):
        load_probe_list(path)


def test_cli_targets_adds_endpoints_and_copies_routes_not_secrets(tmp_path, monkeypatch):
    monkeypatch.setenv("SITE_SECRET", "do-not-copy-me")
    manifest = tmp_path / "site-credentials.json"
    manifest.write_text(
        json.dumps(
            {
                "version": 1,
                "endpoints": {
                    OWNED: {
                        "client_id": SITE_CLIENT,
                        "client_secret_env": "SITE_SECRET",
                        "expected_user": "finemap_svc",
                    }
                },
            }
        )
    )
    out = tmp_path / "targets.json"
    result = CliRunner().invoke(
        main,
        [
            "watch", "targets",
            "--federation", str(repo_root() / "website" / "federation.yaml"),
            "--endpoint", f"mbzuai={OWNED}",
            "--endpoint", f"covenant={UP}",
            "--endpoint", f"anl={DOWN}",
            "--skip", "anl",
            "--site-credentials", str(manifest),
            "--out", str(out),
        ],
    )  # fmt: skip
    assert result.exit_code == 0, result.output
    text = out.read_text()
    assert "do-not-copy-me" not in text and "expected_user" not in text
    document = json.loads(text)
    assert document["sites"] == {"mbzuai": [OWNED], "covenant": [UP]}
    assert document["routes"] == {
        OWNED: {"client_id": SITE_CLIENT, "client_secret_env": "SITE_SECRET"}
    }
    assert "as the account in $SITE_SECRET" in result.output
    assert "anl: left out (--skip)" in result.output
    assert (out.stat().st_mode & 0o777) == 0o600


def test_cli_live_refuses_a_route_whose_secret_is_unset(tmp_path, monkeypatch):
    monkeypatch.delenv("SITE_SECRET", raising=False)
    path = tmp_path / "targets.json"
    route = {"client_id": SITE_CLIENT, "client_secret_env": "SITE_SECRET"}
    path.write_text(json.dumps({"sites": {"mbzuai": [OWNED]}, "routes": {OWNED: route}}))
    result = CliRunner().invoke(
        main, ["watch", "live", "--targets", str(path), "--out", str(tmp_path), "--once"]
    )
    assert result.exit_code != 0 and "SITE_SECRET" in result.output


def test_interval_has_a_floor(tmp_path):
    with pytest.raises(LiveStatusError, match="at least"):
        run({"a": [UP]}, tmp_path / "s.json", interval=1, client_factory=FakeClient)


@pytest.mark.parametrize(
    "body",
    ["{}", '{"sites": {}}', '{"sites": {"a": []}}', '{"sites": {"a": ["not-a-uuid"]}}', "["],
)
def test_malformed_targets_are_refused(tmp_path, body):
    path = tmp_path / "targets.json"
    path.write_text(body)
    with pytest.raises(LiveStatusError):
        load_targets(path)


def test_targets_from_federation_skips_placeholders():
    """The example federation is fictional throughout: nothing in it is probed."""
    example = load_federation(repo_root() / "federation.yaml.example")
    targets, notes = targets_from_federation(example)
    assert targets == {}
    assert notes and all("placeholder" in note or "no endpoint" in note for note in notes)


def test_targets_from_federation_renames_and_dedupes(tmp_path):
    text = (repo_root() / "website" / "federation.yaml").read_text(encoding="utf-8")
    path = tmp_path / "federation.yaml"
    path.write_text(text.replace("bbbbbbbb-0000-0000-0000-000000000004", UP))
    targets, _ = targets_from_federation(load_federation(path), {"cdac": "cdac-pune"})
    assert targets == {"cdac-pune": [UP]}


def test_cli_probes_once_into_an_export_directory(tmp_path, monkeypatch):
    monkeypatch.setattr("appfl_bio_suite.core.endpoint._client", FakeClient)
    targets = tmp_path / "targets.json"
    targets.write_text(json.dumps({"sites": {"anl": [UP], "cdac": [DOWN]}}))
    result = CliRunner().invoke(
        main, ["watch", "live", "--targets", str(targets), "--out", str(tmp_path), "--once"]
    )
    assert result.exit_code == 0, result.output
    assert "anl=online  cdac=offline" in result.output
    assert json.loads((tmp_path / STATUS_FILENAME).read_text())["sites"]["cdac"]["online"] == 0


def test_export_names_the_status_file_only_when_asked(tmp_path):
    pytest.importorskip("hivewatch", reason="the [watch] extra is not installed")
    from appfl_bio_suite.core.watch import export_site

    federation = load_federation(repo_root() / "website" / "federation.yaml")
    plain = json.loads((export_site(federation, tmp_path / "a") / "network.map.json").read_text())
    live = json.loads(
        (
            export_site(federation, tmp_path / "b", live_status=STATUS_FILENAME)
            / "network.map.json"
        ).read_text()
    )
    assert "live_status" not in plain["config"]
    assert live["config"]["live_status"] == STATUS_FILENAME


# ---------------------------------------------------------------------------
# scripts/aws/live_map.py -- the parts that run before any AWS call
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def live_map():
    path = repo_root() / "scripts" / "aws" / "live_map.py"
    spec = importlib.util.spec_from_file_location("live_map", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SHA = "0123456789abcdef0123456789abcdef01234567"


def test_user_data_carries_no_secret_and_only_validated_values(live_map):
    script = live_map.render_user_data(
        sha=SHA, host="203-0-113-10.sslip.io", region="us-east-1", name="live-map", interval=60
    )
    assert f"git fetch -q --depth 1 origin {SHA}" in script
    assert "LIVE_MAP_HOST=203-0-113-10.sslip.io" in script
    assert "secret" not in script.lower()
    for bad in (
        dict(sha="main; rm -rf /"),
        dict(host="x.io; curl evil"),
        dict(name="Bad Name"),
        dict(region="us-east-1 && true"),
        dict(interval=1),
    ):
        args = dict(sha=SHA, host="a.example.org", region="us-east-1", name="live-map", interval=60)
        with pytest.raises(live_map.DeployError):
            live_map.render_user_data(**{**args, **bad})


def test_cloudfront_fronts_the_instance_over_https_without_caching(live_map):
    origin = live_map.ec2_public_dns("52.22.193.198", "us-east-1")
    assert origin == "ec2-52-22-193-198.compute-1.amazonaws.com"
    assert live_map.ec2_public_dns("3.4.5.6", "us-west-2") == (
        "ec2-3-4-5-6.us-west-2.compute.amazonaws.com"
    )
    config = live_map.distribution_config(
        name="live-map", origin=origin, cache_policy="C", request_policy="R", reference="x"
    )
    behaviour = config["DefaultCacheBehavior"]
    assert config["Comment"] == "live-map", "destroy and status find it by this"
    assert behaviour["ViewerProtocolPolicy"] == "redirect-to-https"
    assert (behaviour["CachePolicyId"], behaviour["OriginRequestPolicyId"]) == ("C", "R")
    assert config["Origins"]["Items"][0]["DomainName"] == origin
    assert config["Origins"]["Items"][0]["CustomOriginConfig"]["OriginProtocolPolicy"] == (
        "http-only"
    )


def test_update_resets_the_served_host_name(live_map):
    commands = live_map.update_commands(SHA)
    assert any("LIVE_MAP_HOST=/' /etc/live-map.env" in c for c in commands)
    assert any(
        "LIVE_MAP_HOST=map.example.org/'" in c
        for c in live_map.update_commands(SHA, "map.example.org")
    )
    with pytest.raises(live_map.DeployError):
        live_map.update_commands(SHA, "x/; rm -rf /")


def test_targets_are_checked_against_the_published_map(live_map, tmp_path):
    path = tmp_path / "targets.json"
    path.write_text(json.dumps({"sites": {"cdac": [UP], "nz": [DOWN]}}))
    body, warnings, secrets = live_map.read_targets(path, live_map.public_site_ids())
    assert json.loads(body) == {"sites": {"cdac": [UP], "nz": [DOWN]}}
    assert len(warnings) == 1 and "'nz'" in warnings[0]
    assert secrets == []
    with pytest.raises(live_map.DeployError, match="watch targets"):
        live_map.read_targets(tmp_path / "missing.json")


def test_routed_targets_carry_secret_names_and_never_secrets(live_map, tmp_path):
    path = tmp_path / "targets.json"
    route = {"client_id": SITE_CLIENT, "client_secret_env": "SITE_SECRET"}
    path.write_text(json.dumps({"sites": {"mbzuai": [OWNED]}, "routes": {OWNED: route}}))
    body, _, secrets = live_map.read_targets(path)
    assert json.loads(body)["routes"] == {OWNED: route}
    assert secrets == ["SITE_SECRET"]
    with pytest.raises(live_map.DeployError, match="SITE_SECRET"):
        live_map.read_site_secrets(secrets, {})
    assert live_map.read_site_secrets(secrets, {"SITE_SECRET": "x"}) == {"SITE_SECRET": "x"}

    path.write_text(
        json.dumps({"sites": {"mbzuai": [OWNED]}, "routes": {OWNED: {**route, "secret": "x"}}})
    )
    with pytest.raises(live_map.DeployError, match="never the secret"):
        live_map.read_targets(path)


def test_update_keeps_the_stored_probe_list_without_its_secrets(live_map, tmp_path, monkeypatch):
    from types import SimpleNamespace

    monkeypatch.delenv("SITE_SECRET", raising=False)
    path = tmp_path / "targets.json"
    route = {"client_id": SITE_CLIENT, "client_secret_env": "SITE_SECRET"}
    path.write_text(json.dumps({"sites": {"mbzuai": [OWNED]}, "routes": {OWNED: route}}))
    args = SimpleNamespace(targets=path, globus_client=tmp_path / "none.json")
    stored = SimpleNamespace(has_parameters=lambda *keys: True)
    monkeypatch.setattr(live_map, "public_site_ids", lambda: {"mbzuai"})

    targets, globus, secrets = live_map.gather_inputs(args, stored, required=False)
    assert (targets, globus, secrets) == (None, None, None)
    with pytest.raises(live_map.DeployError, match="SITE_SECRET"):
        live_map.gather_inputs(args, stored, required=True)

    monkeypatch.setenv("SITE_SECRET", "x")
    targets, _, secrets = live_map.gather_inputs(args, stored, required=False)
    assert targets is not None and secrets == {"SITE_SECRET": "x"}


def test_globus_client_from_file_or_environment(live_map, tmp_path):
    path = tmp_path / "client.json"
    assert live_map.read_globus_client(path, {}) is None
    path.write_text(json.dumps({"client_id": UP, "client_secret": "s3"}))
    assert live_map.read_globus_client(path, {}) == (UP, "s3")
    environ = {"LIVE_MAP_GLOBUS_CLIENT_ID": DOWN, "LIVE_MAP_GLOBUS_CLIENT_SECRET": "env"}
    assert live_map.read_globus_client(path, environ) == (DOWN, "env")
    path.write_text(json.dumps({"client_id": "nope", "client_secret": "s3"}))
    with pytest.raises(live_map.DeployError):
        live_map.read_globus_client(path, {})


def test_instance_scripts_are_present_and_executable():
    setup = Path(repo_root()) / "scripts" / "aws" / "instance_setup.sh"
    assert setup.stat().st_mode & 0o111
    text = setup.read_text(encoding="utf-8")
    assert "SITE_LIVE_STATUS=1" in text and "live_map_prober.py" in text
