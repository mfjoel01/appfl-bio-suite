"""Two things the command line owes whoever runs it mid-onboarding.

1. **A site with no endpoint yet is a state, not a crash.** `endpoint_uuid` is optional
   so that a site can be added to an experiment, packaged, and sent its bundle before it
   has stood anything up. Every site is in that state between `experiment package` and
   `experiment register` -- which is exactly when a coordinator opens the map, checks
   endpoint status and runs preflight. Each of those used to fail on the `None`.

2. **A problem with the federation config is reported as its message.** Those messages
   are written to name the fix. Raised from a lookup made after the file had loaded, they
   arrived as a traceback instead.
"""

from __future__ import annotations

import pytest
from click.testing import CliRunner

from appfl_bio_suite.cli import main
from appfl_bio_suite.core.config import load_federation
from appfl_bio_suite.core.launch import resolve_run
from appfl_bio_suite.core.preflight import Level, PreflightReport, _check_endpoints

REGISTERED = "cccccccc-0000-0000-0000-000000000001"

FEDERATION = f"""\
schema_version: 1
coordinator:
  identity: a.researcher@example-university.edu
  site: example-university
sites:
  - id: example-university
    name: Example University
    country: United States
    location: {{lat: 41.8781, lng: -87.6298, city: Chicago}}
    scheduler: pbspro
    queue: debug
    account: EXAMPLE-ALLOC
  - id: site-north
    name: Northern Institute of Technology
    country: Canada
    location: {{lat: 43.6532, lng: -79.3832, city: Toronto}}
    scheduler: slurm
    partition: compute
    account: fedlearn
experiments:
  gwas:
    service_account: gwas_svc
    endpoint_name: gwas-ep
    sites:
      - site: example-university
        client_id: Site1
        endpoint_uuid: {REGISTERED}
        data_dir: /home/gwas_svc/data/site1
        output_dir: /home/gwas_svc/outputs
        expected_samples: 18032
      # Packaged and sent its bundle; has not registered an endpoint yet.
      - site: site-north
        client_id: Site2
        data_dir: /scratch/gwas_svc/data/site2
        output_dir: /scratch/gwas_svc/outputs
        expected_samples: 9237
"""


@pytest.fixture
def pending(tmp_path):
    path = tmp_path / "federation.yaml"
    path.write_text(FEDERATION, encoding="utf-8")
    return path


def invoke(*args: str):
    result = CliRunner().invoke(main, [str(a) for a in args])
    # A ClickException exits through SystemExit; anything else is the traceback a user
    # would have seen.
    assert result.exception is None or isinstance(result.exception, SystemExit), (
        f"uncaught {type(result.exception).__name__}: {result.exception}"
    )
    return result


# ---------------------------------------------------------------------------
# a site that has not registered an endpoint yet
# ---------------------------------------------------------------------------


def test_endpoint_status_reports_a_pending_site_without_probing_it(pending, monkeypatch):
    probed = []

    def fake_status(uuid, label="", **_):
        from appfl_bio_suite.core.endpoint import EndpointStatus

        probed.append(uuid)
        return EndpointStatus(uuid=uuid, label=label, status="online")

    monkeypatch.setattr("appfl_bio_suite.core.endpoint.get_status", fake_status)
    result = invoke("endpoint", "status", "--federation", pending)
    assert probed == [REGISTERED]
    assert "gwas/Site2" in result.output and "no endpoint registered yet" in result.output
    # One site is not up, so the command is not a success -- but it is an answer.
    assert result.exit_code == 1


def test_endpoint_smoke_does_not_submit_to_a_pending_site(pending, monkeypatch):
    from appfl_bio_suite.core.endpoint import SmokeResult

    submitted = []

    def fake_smoke(uuid, label="", timeout=0, **_):
        submitted.append(uuid)
        return SmokeResult(uuid=uuid, label=label, ok=True, elapsed=0.1, payload={"user": "x"})

    monkeypatch.setattr("appfl_bio_suite.core.endpoint.smoke_test", fake_smoke)
    result = invoke("endpoint", "smoke", "--federation", pending)
    assert submitted == [REGISTERED]
    assert "no endpoint registered yet" in result.output
    assert result.exit_code == 1


def test_preflight_names_the_pending_site_instead_of_asking_about_none(pending, monkeypatch):
    from appfl_bio_suite.core.endpoint import EndpointStatus

    asked = []

    def fake_status(uuid, label="", **_):
        asked.append(uuid)
        return EndpointStatus(uuid=uuid, label=label, status="online")

    monkeypatch.setattr("appfl_bio_suite.core.endpoint.get_status", fake_status)
    report = PreflightReport()
    _check_endpoints(report, load_federation(pending), "gwas")
    assert asked == [REGISTERED]
    by_name = {check.name: check for check in report.checks}
    assert by_name["endpoint gwas/Site1"].level is Level.OK
    assert by_name["endpoint gwas/Site2"].level is Level.FAIL
    assert "no endpoint registered yet" in by_name["endpoint gwas/Site2"].detail
    assert "experiment\nregister" in by_name["endpoint gwas/Site2"].fix


def test_a_real_run_still_refuses_a_pending_site_by_name(pending, tmp_path):
    result = invoke("run", "gwas", "--dry-run", "--federation", pending, "--out-dir", tmp_path)
    assert result.exit_code == 1
    assert "Sites must register endpoints before launch: Site2" in result.output


def test_a_serial_run_does_not_need_endpoints(pending):
    """It contacts none: the loopback driver drops `endpoint_id` before building a client."""
    fed = load_federation(pending)
    with pytest.raises(ValueError, match="must register endpoints"):
        resolve_run(fed, "gwas")
    run = resolve_run(fed, "gwas", variant="loopback", allow_unregistered=True)
    assert [c["endpoint_id"] for c in run.client_configs] == [REGISTERED, None]


def test_site_credentials_are_only_required_for_registered_endpoints(pending, tmp_path):
    """A pending site has nothing to route to; it used to fail the manifest check on None."""
    import json

    from appfl_bio_suite.core.site_credentials import load_site_credentials

    manifest = tmp_path / "routes.json"
    manifest.write_text(
        json.dumps(
            {
                "version": 1,
                "endpoints": {REGISTERED: {"coordinator_login": True, "expected_user": "svc"}},
            }
        )
    )
    assert list(load_site_credentials(manifest, [REGISTERED, None])) == [REGISTERED]


# ---------------------------------------------------------------------------
# config problems are messages
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("args", "expected"),
    [
        (
            ["identity", "validate", "mapping.json", "--experiment", "fine-mapping"],
            "experiment 'fine-mapping' is not declared",
        ),
        (["partner-bundle", "gwas", "--site", "nowhere"], "does not participate"),
        (["ga4gh", "duo", "check", "--experiment", "fine-mapping"], "is not declared"),
        (["ga4gh", "trs", "verify", "--experiment", "fine-mapping"], "is not declared"),
        (["ga4gh", "tes", "task", "--experiment", "gwas"], "no `ga4gh.tes` block"),
        (["ga4gh", "drs", "resolve", "drs://host/abc"], "no `ga4gh.drs` block"),
        (["ga4gh", "drs", "serve"], "no `ga4gh.drs` block"),
        (["run", "fine-mapping", "--backend", "multisusie", "--dry-run"], "is not declared"),
    ],
)
def test_a_lookup_the_config_cannot_answer_is_reported_as_its_message(pending, args, expected):
    result = invoke(*args, "--federation", pending)
    assert result.exit_code == 1
    assert result.output.startswith("Error: ") or "\nError: " in result.output
    assert expected in result.output


def test_packaging_an_undeclared_experiment_is_a_message(pending, tmp_path):
    result = invoke(
        "experiment", "package", "fine-mapping", "--federation", pending, "--out", tmp_path
    )
    assert result.exit_code == 1
    assert "experiment 'fine-mapping' is not declared" in result.output


def test_an_unknown_config_variant_lists_the_ones_that_exist(pending, tmp_path):
    registered = "        client_id: Site2\n        endpoint_uuid: " + REGISTERED[:-1] + "2\n"
    pending.write_text(pending.read_text().replace("        client_id: Site2\n", registered))
    result = invoke(
        "run", "gwas", "--config", "nosuch", "--dry-run", "--federation", pending,
        "--out-dir", tmp_path,
    )  # fmt: skip
    assert result.exit_code == 1
    assert "no config template 'server.nosuch.yaml'" in result.output
    assert "server.loopback.yaml" in result.output


def test_preflight_env_for_an_experiment_the_config_does_not_declare(pending):
    """A question about the environment. It used to die looking the experiment up."""
    result = invoke(
        "preflight", "--experiment", "fine-mapping", "--check", "env", "--federation", pending
    )
    assert "python minor version" in result.output
    assert "external binaries" in result.output


def test_a_drs_block_with_no_registry_says_so(pending):
    pending.write_text(pending.read_text() + "ga4gh:\n  drs:\n    hostname: drs.example.org\n")
    result = invoke("ga4gh", "drs", "resolve", "drs://drs.example.org/abc", "--federation", pending)
    assert result.exit_code == 1
    assert "names no `registry`" in result.output


def test_an_experiment_with_no_sites_names_itself_in_the_fix(pending):
    """The second line of this message was never interpolated: it said `{experiment}`."""
    text = pending.read_text()
    pending.write_text(text[: text.index("    sites:\n      - site: example-university")])
    with pytest.raises(RuntimeError) as raised:
        resolve_run(load_federation(pending), "gwas")
    assert "experiments.gwas.sites" in str(raised.value)
    assert "{experiment}" not in str(raised.value)
