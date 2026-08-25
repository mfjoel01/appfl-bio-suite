"""``appfl-bio-suite`` -- one command for every operational task.

Replaces a scattering of shell scripts and one-off Python files that had been
reimplemented independently in two projects: three endpoint probes, two identity-mapping
validators, two driver launchers, and a status script. They did the same things with
different hardcoded UUIDs, and a fix to one never reached the other.

Everything here reads from ``federation.yaml``. No command contains a coordinator
identity, an endpoint UUID, or a site name, which is what lets a different coordinator
run the same commands against their own federation with no code change.
"""

from __future__ import annotations

import sys
from pathlib import Path

import click

from appfl_bio_suite import __version__
from appfl_bio_suite.core.config import FederationError, load_federation
from appfl_bio_suite.core.experiments import REGISTRY, experiment_names, get_spec

_FEDERATION_HELP = (
    "Path to federation.yaml. Defaults to local/federation.yaml, then ./federation.yaml, "
    "or $APPFL_BIO_SUITE_FEDERATION."
)


def _load(path: str | None, required: bool = True):
    try:
        return load_federation(path)
    except FederationError as exc:
        if not required:
            return None
        raise click.ClickException(str(exc)) from exc


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="appfl-bio-suite")
def main() -> None:
    """Federated learning experiments in computational biology.

    New here? Read docs/coordinator/new-federation.md, then:

    \b
      cp federation.yaml.example local/federation.yaml
      $EDITOR local/federation.yaml
      appfl-bio-suite preflight
    """


# ---------------------------------------------------------------------------
# preflight
# ---------------------------------------------------------------------------


@main.command()
@click.option("--experiment", type=click.Choice(experiment_names()), default=None)
@click.option(
    "--check",
    type=click.Choice(["env", "pins", "configs", "data", "endpoints", "all"]),
    default="all",
    help="Run one group of checks. Everything except 'endpoints' and 'all' is offline.",
)
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def preflight(experiment: str | None, check: str, federation_path: str | None) -> None:
    """Check the environment, pins, configs, and endpoints before a run.

    Fails hard only on things that are wrong everywhere -- version skew, a broken config.
    Anything that depends on which cluster you are on is a warning, because a coordinator
    elsewhere has different constraints.
    """
    from appfl_bio_suite.core.preflight import run_preflight

    # Endpoint checks need a federation; the rest are useful without one, which matters
    # when someone is validating a fresh install before they have written their config.
    needs_federation = check in ("endpoints", "all")
    fed = _load(federation_path, required=needs_federation)

    report = run_preflight(federation=fed, experiment=experiment, check=check)
    click.echo(report.render())
    sys.exit(0 if report.ok else 1)


# ---------------------------------------------------------------------------
# endpoint
# ---------------------------------------------------------------------------


@main.group()
def endpoint() -> None:
    """Inspect and test Globus Compute endpoints."""


@endpoint.command("status")
@click.argument("site", required=False)
@click.option("--experiment", type=click.Choice(experiment_names()), default=None)
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def endpoint_status(site: str | None, experiment: str | None, federation_path: str | None) -> None:
    """Report whether endpoints are online. SITE may be a site id or a client id.

    The cloud status is authoritative. Local `globus-compute-endpoint list` reads a
    pidfile that, on a shared filesystem, can record a PID from a different login node.
    """
    from appfl_bio_suite.core.endpoint import get_status

    fed = _load(federation_path)
    targets = _resolve_targets(fed, experiment, site)
    if not targets:
        raise click.ClickException("no endpoints matched. Check --experiment and SITE.")

    any_down = False
    for label, uuid in targets:
        status = get_status(uuid, label=label)
        click.echo(status.render())
        any_down |= not status.online
    sys.exit(1 if any_down else 0)


@endpoint.command("smoke")
@click.argument("site", required=False)
@click.option("--experiment", type=click.Choice(experiment_names()), default=None)
@click.option("--timeout", default=300.0, show_default=True, help="Seconds to wait.")
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def endpoint_smoke(
    site: str | None, experiment: str | None, timeout: float, federation_path: str | None
) -> None:
    """Round-trip a trivial task through an endpoint and report where it ran.

    This is the real test, not `status`. It exercises authentication, identity mapping,
    the scheduler, the environment, and the return path -- every one of which has failed
    independently while the endpoint reported 'online'.

    The 'user' in the result is the load-bearing field: on a partner's multi-user
    endpoint it must be the experiment's service account. If it is the account that
    started the endpoint, identity mapping was skipped.

    Runs no experiment code and touches no data, so it is safe at any time.
    """
    from appfl_bio_suite.core.endpoint import smoke_test

    fed = _load(federation_path)
    targets = _resolve_targets(fed, experiment, site)
    if not targets:
        raise click.ClickException("no endpoints matched. Check --experiment and SITE.")

    expected = {}
    for exp in fed.experiments.values():
        for entry in exp.sites:
            expected[entry.endpoint_uuid] = exp.service_account

    failed = False
    for label, uuid in targets:
        click.echo(f"submitting to {label} ({uuid}) ...")
        result = smoke_test(uuid, label=label, timeout=timeout)
        click.echo(result.render())

        if result.ok:
            want = expected.get(uuid)
            got = (result.payload or {}).get("user")
            if want and got and got != want:
                click.echo(
                    f"      WARNING: ran as '{got}', but this experiment's service\n"
                    f"      account is '{want}'. If this is a multi-user endpoint, the\n"
                    f"      identity mapping is not taking effect -- most likely the\n"
                    f"      endpoint was started by an unprivileged user."
                )
        failed |= not result.ok
        click.echo("")
    sys.exit(1 if failed else 0)


def _resolve_targets(fed, experiment: str | None, site: str | None) -> list[tuple[str, str]]:
    """Turn --experiment/SITE into (label, uuid) pairs."""
    names = [experiment] if experiment else list(fed.enabled_experiments())
    targets: list[tuple[str, str]] = []
    for name in names:
        try:
            exp = fed.experiment(name)
        except FederationError:
            continue
        for entry in exp.sites:
            if site and site not in (entry.site, entry.client_id):
                continue
            targets.append((f"{name}/{entry.client_id}", entry.endpoint_uuid))

    # A coordinator that also trains has its own endpoint; include it when unfiltered.
    if not site and fed.coordinator.endpoint:
        targets.append(
            (
                f"coordinator/{fed.coordinator.endpoint.name}",
                fed.coordinator.endpoint.uuid,
            )
        )
    return targets


# ---------------------------------------------------------------------------
# identity
# ---------------------------------------------------------------------------


@main.group()
def identity() -> None:
    """Work with Globus identity mappings."""


@identity.command("validate")
@click.argument("mapping_file", type=click.Path(path_type=Path))
@click.option("--experiment", type=click.Choice(experiment_names()), default=None)
@click.option("--identity", "identity_str", default=None, help="Override the coordinator identity.")
@click.option("--expect", default=None, help="Override the expected local account.")
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def identity_validate(
    mapping_file: Path,
    experiment: str | None,
    identity_str: str | None,
    expect: str | None,
    federation_path: str | None,
) -> None:
    """Check that a MEP's identity-mapping file will accept us. Offline, ~1 second.

    Pass the exact path `identity_mapping_config_path` names in the endpoint's
    config.yaml. Editing a different copy of the file changes nothing, and after
    switching to a privileged user that path may still point into a personal home
    directory.

    Catches the failure that has cost the most time on this project: `^` and `$` anchors
    in `match`. That field is not a full regex -- the mapper escapes anchors into literal
    characters and then anchors the pattern itself, so an anchored pattern can never
    match, and every submission fails with a 422.

    \b
    A partner can run this themselves, in the endpoint's environment:
        appfl-bio-suite identity validate <the path config.yaml names> \\
            --identity <coordinator identity> --expect <service account>
    """
    from appfl_bio_suite.core.identity import validate_mapping_file

    identity_id = None
    if identity_str is None or expect is None:
        fed = _load(federation_path)
        identity_str = identity_str or fed.coordinator.identity
        identity_id = fed.coordinator.identity_id
        if expect is None:
            if experiment is None:
                raise click.ClickException(
                    "specify --experiment (to take the service account from "
                    "federation.yaml) or pass --expect explicitly."
                )
            expect = fed.experiment(experiment).service_account

    report = validate_mapping_file(
        mapping_file, identity=identity_str, expected_output=expect, identity_id=identity_id
    )
    click.echo(report.render())
    sys.exit(0 if report.passed else 1)


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------


@main.command("run")
@click.argument("experiment", type=click.Choice(experiment_names()))
@click.option("--config", "variant", default="default", show_default=True,
              help="Server config variant, e.g. 'fedcompass' or 'loopback'.")
@click.option(
    "--driver",
    type=click.Choice(["globus_compute", "serial"]),
    default="globus_compute",
    show_default=True,
    help="'serial' is the loopback path: all sites in one process, no Globus.",
)
@click.option("--out-dir", type=click.Path(path_type=Path), default=None,
              help="Where to write the resolved configs.")
@click.option("--dry-run", is_flag=True, help="Resolve and write configs, but do not launch.")
@click.option("--data-root", type=click.Path(path_type=Path), default=None,
              help="Loopback only: the --out directory `simulate` wrote per-site data to.")
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def run_cmd(
    experiment: str,
    variant: str,
    driver: str,
    out_dir: Path | None,
    dry_run: bool,
    data_root: Path | None,
    federation_path: str | None,
) -> None:
    """Run an experiment.

    Resolves federation.yaml into APPFL server and client configs, runs a preflight, and
    launches. --dry-run stops after writing the configs, which is the fastest way to see
    exactly what a partner's endpoint will be sent.

    \b
    --driver serial is the loopback path and needs no federation config at all: it
    synthesizes one describing simulated sites on this machine, so an install can be
    proved before anything has been filled in. If you do have one it is used, with its
    worker-side paths repointed here -- they name directories on partner clusters.
    """
    # Planned-but-unimplemented experiments are accepted and then refused by name, the
    # same way `simulate` refuses them. Leaving them out of the choice list instead gave
    # a generic "is not one of" usage error that named neither the reason nor the doc.
    if not get_spec(experiment).implemented:
        raise click.ClickException(
            f"'{experiment}' is planned but not implemented. "
            f"See docs/experiments/{experiment}/ABOUT.md."
        )

    from appfl_bio_suite.core.launch import launch
    from appfl_bio_suite.core.loopback import localize_for_loopback, loopback_federation

    loopback = driver == "serial"
    fed = _load(federation_path, required=not loopback)

    if loopback:
        try:
            if fed is None:
                fed = loopback_federation(experiment, data_root)
            else:
                localize_for_loopback(fed, experiment, data_root)
        except FederationError as exc:
            raise click.ClickException(str(exc)) from exc
    elif data_root is not None:
        raise click.ClickException(
            "--data-root only applies to a loopback run. It repoints every site at "
            "simulated data on THIS machine, which is not what a real federation wants: "
            "each site's data_dir is a path on its own cluster. Add --driver serial, or "
            "drop --data-root."
        )

    code = launch(
        federation=fed,
        experiment=experiment,
        variant=variant,
        out_dir=out_dir,
        dry_run=dry_run,
        driver=driver,
    )
    sys.exit(code)


# ---------------------------------------------------------------------------
# simulate
# ---------------------------------------------------------------------------


@main.command("simulate")
@click.argument("experiment", type=click.Choice(experiment_names()))
@click.option("--scenario", default=None, help="Scenario name or path to a scenario YAML.")
@click.option("--out", "out_dir", type=click.Path(path_type=Path), default=None,
              help="Output directory for the simulated per-site data.")
@click.option("--verify", "verify_manifest", type=click.Path(path_type=Path), default=None,
              help="Re-checksum outputs against a recorded run manifest and report drift.")
@click.option("--list-scenarios", is_flag=True, help="List available scenarios and exit.")
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def simulate_cmd(
    experiment: str,
    scenario: str | None,
    out_dir: Path | None,
    verify_manifest: Path | None,
    list_scenarios: bool,
    federation_path: str | None,
) -> None:
    """Generate synthetic data and split it into per-site datasets.

    Coordinator-side only. This never runs on a partner cluster, and its dependencies
    live in a separate extra so they never land on one.

    Every run writes a manifest recording the scenario, every seed, the suite commit,
    input and output checksums, and package versions. --verify re-checksums an output
    directory against one, which answers "did this rerun produce the same data?"
    mechanically rather than by eyeballing summary statistics.
    """
    spec = get_spec(experiment)
    if not spec.has_simulation:
        raise click.ClickException(
            f"'{experiment}' has no simulation stage -- partners obtain the data "
            f"themselves. See docs/experiments/{experiment}/DATA.md."
        )
    if not spec.implemented:
        raise click.ClickException(
            f"'{experiment}' is planned but not implemented. "
            f"See docs/experiments/{experiment}/ABOUT.md."
        )

    # `simulation_scenario` in federation.yaml is the scenario that produced the data
    # this federation runs on. Defaulting to it here is what makes it a fact about the
    # federation rather than a note: rerunning `simulate` without arguments regenerates
    # the data the config describes, instead of asking which scenario that was.
    if scenario is None and not list_scenarios and verify_manifest is None:
        fed = _load(federation_path, required=False)
        declared = None
        if fed is not None:
            declared = getattr(fed.experiments.get(experiment), "simulation_scenario", None)
        if declared:
            scenario = declared
            click.echo(f"scenario '{scenario}', from {fed.source_path}")

    # Dispatched from the registry rather than imported by name. This line used to read
    # `from appfl_bio_suite.experiments.gwas.simulation import cli_entry`, which meant
    # `simulate fine-mapping` would have run the GWAS pipeline -- the exact class of bug
    # the registry exists to prevent, sitting inside the command that consults it.
    import importlib

    module = f"appfl_bio_suite.experiments.{spec.package}.simulation"
    try:
        cli_entry = importlib.import_module(module).cli_entry
    except (ImportError, AttributeError) as exc:
        raise click.ClickException(
            f"experiment '{experiment}' declares has_simulation=True but {module} does "
            f"not provide a cli_entry(): {exc}"
        ) from exc

    sys.exit(
        cli_entry(
            scenario=scenario,
            out_dir=out_dir,
            verify_manifest=verify_manifest,
            list_scenarios=list_scenarios,
        )
    )


# ---------------------------------------------------------------------------
# partner-bundle
# ---------------------------------------------------------------------------


@main.command("partner-bundle")
@click.argument("experiment", type=click.Choice(experiment_names(implemented_only=True)))
@click.option("--site", required=True, help="Site id or client id from federation.yaml.")
@click.option("--out", "out_dir", type=click.Path(path_type=Path), default=None,
              help="Where to write the bundle. Defaults to local/partner_bundles/.")
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def partner_bundle_cmd(
    experiment: str, site: str, out_dir: Path | None, federation_path: str | None
) -> None:
    """Generate one partner's complete, fully-filled-in setup bundle.

    Emits their config with every value already resolved -- their client id, their data
    assignment, YOUR identity, their expected sample count -- plus their two setup
    documents, with no placeholders left for them to interpret.

    Every question partners have asked has been about which placeholder applied to them.
    Each partner runs one experiment at one data center, so there is no reason to make
    them choose.

    Documents are copied only from docs/partner/. Nothing outside that tree can end up in
    a bundle, which makes "did I accidentally send them my internal notes" a structural
    impossibility rather than a judgment call.
    """
    from appfl_bio_suite.core.partner import generate_bundle

    fed = _load(federation_path)
    destination = generate_bundle(fed, experiment, site, out_dir)
    click.echo(f"\nBundle written to: {destination}")
    click.echo("\nContents:")
    for path in sorted(destination.rglob("*")):
        if path.is_file():
            click.echo(f"  {path.relative_to(destination)}")


# ---------------------------------------------------------------------------
# experiments
# ---------------------------------------------------------------------------


@main.command("experiments")
def experiments_cmd() -> None:
    """List the experiments this build knows about."""
    for name, spec in REGISTRY.items():
        state = "implemented" if spec.implemented else "PLANNED, not implemented"
        sim = "generates its own data" if spec.has_simulation else "partners obtain the data"
        click.echo(f"\n{name}")
        click.echo(f"  {spec.title}")
        click.echo(f"  status: {state}   |   data: {sim}")
        click.echo(f"  install: pip install 'appfl-bio-suite[{','.join(spec.partner_extras)}]'")
        for line in _wrap(spec.summary, 74):
            click.echo(f"  {line}")


def _wrap(text: str, width: int) -> list[str]:
    import textwrap

    return textwrap.wrap(text, width=width)


if __name__ == "__main__":
    main()
