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

from appfl_bio_suite import __version__, install_spec
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
    type=click.Choice(["env", "pins", "configs", "data", "ga4gh", "endpoints", "all"]),
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
    needs_federation = check in ("endpoints", "ga4gh", "all")
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
@click.option(
    "--config",
    "variant",
    default="default",
    show_default=True,
    help="Server config variant, e.g. 'fedcompass' or 'loopback'.",
)
@click.option(
    "--driver",
    type=click.Choice(["globus_compute", "serial", "tes"]),
    default="globus_compute",
    show_default=True,
    help="'serial' is the loopback path: all sites in one process, no Globus. 'tes' "
    "dispatches the site stage as GA4GH Task Execution Service tasks; it needs a "
    "published container image and a TES service per site.",
)
@click.option(
    "--out-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Where to write the resolved configs.",
)
@click.option("--dry-run", is_flag=True, help="Resolve and write configs, but do not launch.")
@click.option(
    "--watch",
    is_flag=True,
    help="Stream this run onto the federation map (`watch serve`). Never fatal: if "
    "hivewatch is missing the run proceeds without it.",
)
@click.option(
    "--data-root",
    type=click.Path(path_type=Path),
    default=None,
    help="Loopback only: the --out directory `simulate` wrote per-site data to.",
)
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def run_cmd(
    experiment: str,
    variant: str,
    driver: str,
    out_dir: Path | None,
    dry_run: bool,
    watch: bool,
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
        watch=watch,
    )
    sys.exit(code)


# ---------------------------------------------------------------------------
# watch
# ---------------------------------------------------------------------------

_WATCH_OPTIONS = [
    click.option(
        "--catalog",
        "catalog_path",
        type=click.Path(path_type=Path, exists=True),
        default=None,
        help="Optional JSON partner roster and selected result artifacts.",
    ),
    click.option(
        "--experiment",
        type=click.Choice(experiment_names()),
        default=None,
        help="Draw only this experiment's sites. Default: every enabled experiment.",
    ),
    click.option(
        "--probe",
        is_flag=True,
        help="Ask Globus Compute whether each endpoint is up, and colour the map by the "
        "answer. Costs one API call per endpoint; without it the map states membership, "
        "not liveness.",
    ),
    click.option(
        "--include-endpoint-uuids",
        is_flag=True,
        help="Publish full endpoint UUIDs instead of fingerprints. For an internal "
        "deployment; see core/watch.py on why the default is the other way.",
    ),
    click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP),
]


def _watch_options(command):
    for option in reversed(_WATCH_OPTIONS):
        command = option(command)
    return command


@main.group()
def watch() -> None:
    """The federation network map: who is in it, where, and what they hold.

    Built on hivewatch (https://github.com/APPFL/hivewatch), which supplies the event
    schema, the map metadata format, and the viewer. What this adds is the federation:
    hivewatch knows how to draw a client, and federation.yaml knows who the clients are.

    \b
    A federation exists before any run does, so the map does too. It is drawn from
    federation.yaml alone -- no run required, nothing to launch first.

    \b
      appfl-bio-suite watch build            refresh it from federation.yaml
      appfl-bio-suite watch serve            open it locally, live
      appfl-bio-suite watch export --out DIR a static copy anyone can host

    Needs the coordinator-only extra:  pip install 'appfl-bio-suite[watch]'
    """


def _statuses(fed, experiment: str | None, probe: bool) -> dict[str, str] | None:
    """Site id -> hivewatch status, by asking Globus Compute. None when not probing."""
    if not probe:
        return None

    from appfl_bio_suite.core.endpoint import get_status
    from appfl_bio_suite.core.watch import participations

    out: dict[str, str] = {}
    for site_id, parts in participations(fed, experiment).items():
        online = 0
        for part in parts:
            label = f"{part.experiment}/{part.client_id}"
            status = get_status(part.endpoint_uuid, label=label)
            click.echo(status.render())
            online += int(status.online)
        # A site with two endpoints and one of them down is neither up nor down. It gets
        # `idle`, which the viewer draws differently from both -- a site half-online is
        # exactly the state worth being able to see at a glance.
        out[site_id] = "active" if online == len(parts) else "idle" if online else "failed"
    click.echo("")
    return out


@watch.command("build")
@_watch_options
@click.option(
    "--runs-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Where run artifacts live. Defaults to local/watch/runs.",
)
def watch_build(
    catalog_path: Path | None,
    experiment: str | None,
    probe: bool,
    include_endpoint_uuids: bool,
    federation_path: str | None,
    runs_dir: Path | None,
) -> None:
    """Refresh the network view from federation.yaml.

    Writes it as an ordinary hivewatch run under a fixed id, so it sits in the same runs
    directory as your real runs and appears alongside them in the viewer's run list.
    Rerun it whenever federation.yaml changes; it replaces rather than accumulates.
    """
    from appfl_bio_suite.core.watch import DEFAULT_RUNS_DIR, WatchError, unplaced_report
    from appfl_bio_suite.core.watch import write_network_run as write_run

    fed = _load(federation_path)
    try:
        _, mapjson = write_run(
            fed,
            runs_dir or DEFAULT_RUNS_DIR,
            catalog_path=catalog_path,
            experiment=experiment,
            statuses=_statuses(fed, experiment, probe),
            include_endpoint_uuids=include_endpoint_uuids,
        )
    except WatchError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"network view: {mapjson}")
    report = unplaced_report(fed, experiment)
    if report:
        click.echo("")
        click.echo(report)


@watch.command("serve")
@click.option("--host", default="0.0.0.0", show_default=True)
@click.option("--port", default=None, type=int, help="Defaults to 7070.")
@click.option(
    "--runs-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Where run artifacts live. Defaults to local/watch/runs.",
)
def watch_serve(host: str, port: int | None, runs_dir: Path | None) -> None:
    """Serve the map and every run in the runs directory, live.

    Watches the directory, so a run launched with `run --watch` appears here as it
    happens without restarting anything.

    \b
    THIS IS THE LOCAL VIEW. It binds a port on the machine you run it on, which on an
    HPC login node is reachable by you and by nobody else -- partners cannot open it, and
    the cluster is right not to let them. `watch export` is what you hand out.
    """
    from appfl_bio_suite.core.watch import DEFAULT_PORT, DEFAULT_RUNS_DIR, WatchError, map_server

    runs_dir = runs_dir or DEFAULT_RUNS_DIR
    port = port or DEFAULT_PORT
    runs_dir.mkdir(parents=True, exist_ok=True)

    if not any(runs_dir.glob("*.jsonl")):
        click.echo(f"{runs_dir} is empty -- run `appfl-bio-suite watch build` first.\n")

    try:
        with map_server(host=host, port=port, runs_dir=runs_dir) as server:
            server.start()
            click.echo(f"map:  http://localhost:{port}")
            click.echo(f"runs: {runs_dir}")
            click.echo("Ctrl-C to stop.")
            server.serve_forever()
    except KeyboardInterrupt:
        pass
    except WatchError as exc:
        raise click.ClickException(str(exc)) from exc


@watch.command("export")
@_watch_options
@click.option(
    "--out",
    "out_dir",
    type=click.Path(path_type=Path),
    required=True,
    help="Directory to write the static site into.",
)
@click.option("--title", default="APPFL federation network", show_default=True)
def watch_export(
    catalog_path: Path | None,
    experiment: str | None,
    probe: bool,
    include_endpoint_uuids: bool,
    federation_path: str | None,
    out_dir: Path,
    title: str,
) -> None:
    """Write a static copy of the network map that anyone can host.

    Three files: the viewer, the data, and an index. No server of ours, no Python at the
    other end. Drop it on GitHub Pages, an institutional web host, or any object store,
    and every coordinator and partner in the federation has the same page.

    \b
    CHECK WHAT YOU ARE PUBLISHING. Site names, countries, coordinates, sample counts and
    experiment names are all in network.map.json by design -- that is the map. Service
    accounts, cluster paths, your Globus identity and (unless you pass
    --include-endpoint-uuids) full endpoint UUIDs are not.
    """
    from appfl_bio_suite.core.watch import WatchError, export_site, unplaced_report

    fed = _load(federation_path)
    try:
        destination = export_site(
            fed,
            out_dir,
            catalog_path=catalog_path,
            title=title,
            experiment=experiment,
            statuses=_statuses(fed, experiment, probe),
            include_endpoint_uuids=include_endpoint_uuids,
        )
    except WatchError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Static site written to: {destination}")
    for path in sorted(destination.iterdir()):
        click.echo(f"  {path.name}")
    click.echo("")
    click.echo("Preview it locally (a plain file:// open will not fetch the data):")
    click.echo(f"  python -m http.server -d {destination} 8000")

    report = unplaced_report(fed, experiment)
    if report:
        click.echo("")
        click.echo(report)


# ---------------------------------------------------------------------------
# simulate
# ---------------------------------------------------------------------------


@main.command("simulate")
@click.argument("experiment", type=click.Choice(experiment_names()))
@click.option("--scenario", default=None, help="Scenario name or path to a scenario YAML.")
@click.option(
    "--out",
    "out_dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Output directory for the simulated per-site data.",
)
@click.option(
    "--verify",
    "verify_manifest",
    type=click.Path(path_type=Path),
    default=None,
    help="Re-checksum outputs against a recorded run manifest and report drift.",
)
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
@click.option(
    "--out",
    "out_dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Where to write the bundle. Defaults to local/partner_bundles/.",
)
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
# ga4gh
# ---------------------------------------------------------------------------


@main.group()
def ga4gh() -> None:
    """GA4GH standards: data use terms, data objects, tools, and task execution.

    \b
      duo   Data Use Ontology       may this site's data be used for this study?
      drs   Data Repository Service exactly which bytes did a site compute over?
      trs   Tool Registry Service   exactly which tool version computed them?
      tes   Task Execution Service  run that tool, there, on those bytes.

    Every one of these is opt-in per experiment, and a federation using none of them runs
    exactly as it did before. See docs/coordinator/ga4gh.md.
    """


# -- duo --------------------------------------------------------------------


@ga4gh.group("duo")
def ga4gh_duo() -> None:
    """Data use terms: what a dataset permits, and what this study intends."""


@ga4gh_duo.command("check")
@click.option("--experiment", type=click.Choice(experiment_names()), default="fine-mapping")
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def ga4gh_duo_check(experiment: str, federation_path: str | None) -> None:
    """Match every site's DUO terms against this experiment's data use request.

    The same evaluation each site's worker performs on its own copy, run offline against
    the copies you hold. A site that refuses here will refuse there -- so this is how you
    find out before spending a scheduler allocation to be told.

    Exits non-zero if any site's terms do not permit the study.
    """
    from appfl_bio_suite.core.ga4gh.resolve import resolve_data_use

    fed = _load(federation_path)
    request = fed.data_use_request(experiment)
    if request is None:
        raise click.ClickException(
            f"experiment '{experiment}' declares no `ga4gh.data_use_request` in "
            f"{fed.source_path}. A study has to say what it is before its access to "
            "anything can be evaluated."
        )
    click.echo(request.render())
    click.echo("")

    decisions = resolve_data_use(fed, experiment)
    if not decisions:
        raise click.ClickException(f"experiment '{experiment}' has no participating sites.")
    for decision in decisions:
        click.echo(decision.render())

    blocking = [d for d in decisions if d.blocking]
    click.echo("")
    if blocking:
        click.echo(f"{len(blocking)} site(s) do not permit this study.")
    else:
        click.echo("Every site with a recorded profile permits this study.")
    sys.exit(1 if blocking else 0)


@ga4gh_duo.command("show")
@click.argument("profile", type=click.Path(path_type=Path))
def ga4gh_duo_show(profile: Path) -> None:
    """Read a DATA_USE.json and print its terms in full."""
    from appfl_bio_suite.core.ga4gh.duo import DuoError, describe_term, load_profile

    try:
        loaded = load_profile(profile)
    except DuoError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(loaded.render())
    click.echo("")
    click.echo(describe_term(loaded.permission))
    for entry in loaded.modifiers:
        click.echo(describe_term(entry.id))
    if loaded.description:
        click.echo("")
        for line in _wrap(loaded.description, 78):
            click.echo(line)


@ga4gh_duo.command("terms")
@click.option(
    "--kind",
    type=click.Choice(["permission", "modifier", "purpose", "all"]),
    default="all",
    show_default=True,
)
def ga4gh_duo_terms(kind: str) -> None:
    """List the DUO terms this build understands, from the vendored release."""
    from appfl_bio_suite.core.ga4gh.duo import MODIFIERS, PERMISSIONS, PURPOSES, ontology

    snapshot = ontology()
    click.echo(f"{snapshot['title']}  {snapshot['version_iri']}\n")
    groups = {
        "permission": ("data use permissions (a dataset declares exactly one)", PERMISSIONS),
        "modifier": ("data use modifiers (extra conditions)", MODIFIERS),
        "purpose": ("research purposes (a STUDY declares these)", PURPOSES),
    }
    for key, (title, table) in groups.items():
        if kind not in ("all", key):
            continue
        click.echo(title)
        for term_id, label in table.items():
            click.echo(f"  {term_id}  {label}")
        click.echo("")


# -- drs --------------------------------------------------------------------


@ga4gh.group("drs")
def ga4gh_drs() -> None:
    """Data objects: content-addressed bundles, and where to get them."""


@ga4gh_drs.command("register")
@click.option(
    "--data-root",
    type=click.Path(path_type=Path),
    required=True,
    help="A directory `simulate` wrote: one <site>/data/ per site.",
)
@click.option("--hostname", default=None, help="DRS hostname. Defaults to ga4gh.drs.hostname.")
@click.option("--https-base", default=None, help="Base URL where this registry is served.")
@click.option("--globus-collection", default=None, help="Collection UUID holding the bundles.")
@click.option("--out", "out_path", type=click.Path(path_type=Path), default=None)
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def ga4gh_drs_register(
    data_root: Path,
    hostname: str | None,
    https_base: str | None,
    globus_collection: str | None,
    out_path: Path | None,
    federation_path: str | None,
) -> None:
    """Build a DRS registry over a directory of per-site bundles.

    Ids are content-addressed: a blob's id is its sha-256 and a bundle's is a Merkle hash
    over its members. So two coordinators registering byte-identical data mint identical
    ids, and a bundle whose id changed is a bundle whose content changed.
    """
    from appfl_bio_suite.core.ga4gh.drs import DrsError, build_registry

    fed = _load(federation_path, required=False)
    service = fed.ga4gh.drs if fed and fed.ga4gh else None
    hostname = hostname or (service.hostname if service else None)
    if not hostname:
        raise click.ClickException(
            "no --hostname, and no `ga4gh.drs.hostname` in a federation config. A DRS id "
            "is only unique within its service, so the hostname is part of the identity "
            "-- 'drs.your-org.example', not a URL."
        )

    try:
        registry = build_registry(
            data_root,
            hostname=hostname,
            https_base=https_base or (service.https_base if service else None),
            globus_collection=globus_collection or (service.globus_collection if service else None),
        )
    except DrsError as exc:
        raise click.ClickException(str(exc)) from exc

    destination = registry.save(out_path or (Path(data_root) / "drs_registry.json"))
    click.echo(registry.summary())
    click.echo("")
    # Site bundles only. The nested ones -- each bundle's phenotypes/ subtree -- are
    # addressable too, and listing them here would bury the three lines a coordinator
    # actually has to copy.
    for name, obj in sorted(registry.site_bundles().items()):
        click.echo(f"  {name:12} {obj.self_uri}")
    click.echo(f"\nregistry: {destination}")
    click.echo("\nPut these in federation.yaml as each site's `drs_uri`, and the registry")
    click.echo("path as `ga4gh.drs.registry`.")


@ga4gh_drs.command("resolve")
@click.argument("uri")
@click.option("--registry", "registry_path", type=click.Path(path_type=Path), default=None)
@click.option("--url", "base_url", default=None, help="Resolve over HTTP against a DRS service.")
@click.option("--verify", is_flag=True, help="Re-checksum the local copy, if there is one.")
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def ga4gh_drs_resolve(
    uri: str,
    registry_path: Path | None,
    base_url: str | None,
    verify: bool,
    federation_path: str | None,
) -> None:
    """Resolve a drs:// URI (or a bare id) and print the object."""
    from appfl_bio_suite.core.ga4gh.drs import DrsClient, DrsError, DrsRegistry, verify_object

    registry = None
    if base_url is None:
        if registry_path is None:
            fed = _load(federation_path)
            registry_path = Path(fed.drs_service().registry or "")
        try:
            registry = DrsRegistry.load(registry_path)
        except DrsError as exc:
            raise click.ClickException(str(exc)) from exc

    client = DrsClient(registry=registry, base_url=base_url)
    try:
        obj = client.resolve(uri) if uri.startswith("drs://") else client.get(uri)
    except DrsError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"{obj.name or '(unnamed)'}  {obj.self_uri}")
    click.echo(f"  {'bundle' if obj.is_bundle else 'blob'}, {obj.size:,} bytes")
    for entry in obj.checksums:
        click.echo(f"  {entry.type}: {entry.checksum}")
    for method in obj.access_methods:
        url = method.access_url.url if method.access_url else "(via access_id)"
        click.echo(f"  {method.type}: {url}")
    for child in obj.contents or []:
        click.echo(f"    {child.name}  {child.id[:16]}...")

    if verify:
        problems = verify_object(obj)
        click.echo("")
        click.echo("\n".join(problems) if problems else "verified: content matches")
        sys.exit(1 if problems else 0)


@ga4gh_drs.command("serve")
@click.option("--registry", "registry_path", type=click.Path(path_type=Path), default=None)
@click.option("--host", default="127.0.0.1", show_default=True)
@click.option("--port", default=8080, show_default=True)
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def ga4gh_drs_serve(
    registry_path: Path | None, host: str, port: int, federation_path: str | None
) -> None:
    """Serve the registry as a read-only DRS 1.5 service.

    \b
    NO AUTHORIZATION. Anything this serves is readable by anyone who can reach the port,
    so bind it to localhost or an internal interface. It exists to make a registry
    resolvable during development and for a deployment that puts its own gateway in
    front; it is not a data access control.
    """
    from appfl_bio_suite.core.ga4gh.drs import DrsError, DrsRegistry, serve

    if registry_path is None:
        fed = _load(federation_path)
        registry_path = Path(fed.drs_service().registry or "")
    try:
        registry = DrsRegistry.load(registry_path)
    except DrsError as exc:
        raise click.ClickException(str(exc)) from exc

    server = serve(registry, host=host, port=port)
    click.echo(f"DRS {registry.hostname}: http://{host}:{port}/ga4gh/drs/v1/")
    click.echo(f"  {registry.summary()}")
    click.echo("  no authorization -- do not expose this to a network you do not control")
    click.echo("Ctrl-C to stop.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        server.shutdown()


# -- trs --------------------------------------------------------------------


@ga4gh.group("trs")
def ga4gh_trs() -> None:
    """The site stage as a registered, versioned, checksummed tool."""


@ga4gh_trs.command("publish")
@click.option("--experiment", type=click.Choice(experiment_names()), default="fine-mapping")
@click.option("--out", "out_dir", type=click.Path(path_type=Path), default=Path("local/trs"))
@click.option("--image", default=None, help="Container image the executor runs.")
@click.option("--image-digest", default=None, help="Its sha256 digest. Pin this, not the tag.")
@click.option("--tool-id", default=None, help="TRS id. Defaults to a Dockstore-shaped id.")
@click.option("--registry-url", default=None, help="Where this tree will be served.")
@click.option("--version", "version", default=None, help="Defaults to the suite version.")
def ga4gh_trs_publish(
    experiment: str,
    out_dir: Path,
    image: str | None,
    image_digest: str | None,
    tool_id: str | None,
    registry_url: str | None,
    version: str | None,
) -> None:
    """Write a servable TRS tree, the CWL descriptor, a Containerfile, and the pin.

    The pin is the output that matters: copy it into federation.yaml under
    `experiments.<name>.ga4gh.tool` and every run afterwards records -- and verifies --
    which tool version produced its numbers.
    """
    from appfl_bio_suite.core.ga4gh.trs import write_registry

    destination, pin = write_registry(
        out_dir,
        experiment=experiment,
        version=version,
        tool_id=tool_id,
        image=image,
        image_digest=image_digest,
        registry_url=registry_url,
    )
    click.echo(f"TRS tree: {destination}")
    for path in sorted(destination.rglob("*")):
        if path.is_file():
            click.echo(f"  {path.relative_to(destination)}")
    click.echo("")
    click.echo("Pin for federation.yaml:")
    click.echo("      tool:")
    click.echo(f'        id: "{pin.id}"')
    click.echo(f'        version: "{pin.version}"')
    click.echo(f"        descriptor_checksum: {pin.descriptor_checksum}")
    if pin.image:
        click.echo(f"        image: {pin.image}")

    # Three states, and they call for different things. No image at all is the state a
    # first publish is in and is not a mistake; a tag is a moving target; a digest is
    # what a federation should agree on.
    if pin.image is None:
        click.echo("")
        click.echo("NO CONTAINER IMAGE is registered, so the TES path has nothing to run.")
        click.echo("(The Globus Compute path is unaffected -- it ships source, not images.)")
        click.echo("Build and push one from the generated Containerfile, then:")
        click.echo(
            f"  appfl-bio-suite ga4gh trs publish --out {out_dir} \\\n"
            "      --image <registry>/appfl-bio-suite:<version> --image-digest <sha256>"
        )
    elif not image_digest:
        click.echo("")
        click.echo("The image is pinned by TAG. A tag can be repushed, so a federation that")
        click.echo("agreed on one has agreed on nothing durable. Re-run with --image-digest")
        click.echo("(the sha256 the push reported) to pin the content instead.")


@ga4gh_trs.command("verify")
@click.option("--experiment", type=click.Choice(experiment_names()), default="fine-mapping")
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def ga4gh_trs_verify(experiment: str, federation_path: str | None) -> None:
    """Check this federation's tool pin against the installed package."""
    from appfl_bio_suite.core.ga4gh.resolve import resolve_tool

    fed = _load(federation_path)
    pin = fed.tool_pin(experiment)
    if pin is None:
        raise click.ClickException(
            f"experiment '{experiment}' declares no `ga4gh.tool` pin. Create one with "
            "`appfl-bio-suite ga4gh trs publish`."
        )
    _, problems = resolve_tool(fed, experiment)
    if problems:
        for problem in problems:
            click.echo(problem)
        sys.exit(1)
    click.echo(f"{pin.render()} matches this install.")


@ga4gh_trs.command("show")
@click.argument("tool_id", required=False)
@click.option("--url", "base_url", default=None, help="A TRS registry, e.g. Dockstore's API.")
@click.option("--version", "version_id", default=None)
@click.option("--experiment", type=click.Choice(experiment_names()), default="fine-mapping")
def ga4gh_trs_show(
    tool_id: str | None, base_url: str | None, version_id: str | None, experiment: str
) -> None:
    """Print a tool document -- from a registry with --url, or from this install."""
    from appfl_bio_suite.core.ga4gh.trs import TrsClient, TrsError, build_tool

    if base_url:
        if not tool_id:
            raise click.ClickException("pass the TRS id of the tool to fetch.")
        try:
            client = TrsClient(base_url)
            tool = client.tool(tool_id)
            version = client.version(tool_id, version_id or tool.versions[0].id)
        except TrsError as exc:
            raise click.ClickException(str(exc)) from exc
    else:
        tool, _ = build_tool(experiment=experiment, version=version_id, tool_id=tool_id)
        version = tool.versions[0]

    click.echo(f"{tool.id}\n  {tool.name}")
    click.echo(f"  class: {tool.toolclass.name}   organization: {tool.organization}")
    click.echo(f"  version: {version.id}  descriptors: {', '.join(version.descriptor_type)}")
    for image in version.images:
        digest = image.digest
        suffix = f"  sha256:{digest[:16]}..." if digest else ""
        click.echo(f"  image: {image.image_name}{suffix}")
    for entry in (tool.model_extra or {}).get("files", []):
        checksum = (entry.get("checksum") or [{}])[0].get("checksum", "")
        click.echo(f"    {entry['file_type']:20} {entry['path']}  {checksum[:12]}")


# -- tes --------------------------------------------------------------------


@ga4gh.group("tes")
def ga4gh_tes() -> None:
    """Run the site stage as GA4GH tasks, at sites that speak TES."""


@ga4gh_tes.command("task")
@click.option("--experiment", type=click.Choice(experiment_names()), default="fine-mapping")
@click.option("--site", default=None, help="One site. Default: every participating site.")
@click.option("--out", "out_dir", type=click.Path(path_type=Path), default=None)
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def ga4gh_tes_task(
    experiment: str, site: str | None, out_dir: Path | None, federation_path: str | None
) -> None:
    """Write the TES task document for each site, without submitting anything.

    What to send a partner who asks what will actually run on their cluster: the image,
    the command, the inputs by DRS URI, the resource request, and the data use decision
    that authorized it. Nothing is hidden in a framework.
    """
    from appfl_bio_suite.core.ga4gh.tes import build_site_task, write_task

    fed = _load(federation_path)
    exp = fed.experiment(experiment)
    tes = fed.tes_service()
    pin = fed.tool_pin(experiment)
    image = (pin.image if pin else None) or "<no image pinned>"

    entries = [e for e in exp.sites if site is None or site in (e.site, e.client_id)]
    if not entries:
        raise click.ClickException(f"no site matched '{site}' in experiment '{experiment}'.")

    out_dir = out_dir or Path("local/tes")
    from appfl_bio_suite.core.launch import build_client_configs

    clients = {c["client_id"]: c for c in build_client_configs(fed, experiment)}

    for entry in entries:
        client = clients[entry.client_id]
        task = build_site_task(
            client_id=entry.client_id,
            experiment=experiment,
            image=image,
            bundle_url=entry.drs_uri or "<no drs_uri for this site>",
            outputs_url=entry.tes_outputs_url or "<no tes_outputs_url for this site>",
            run_config={
                "client_id": entry.client_id,
                "train_configs": client.get("train_configs", {}),
                "dataset_kwargs": {
                    k: v
                    for k, v in (client.get("data_configs", {}).get("dataset_kwargs", {})).items()
                    if k != "data_dir"
                },
            },
            cpu_cores=tes.cpu_cores,
            ram_gb=tes.ram_gb,
            disk_gb=tes.disk_gb,
            preemptible=tes.preemptible,
            tags={"trs_id": pin.id if pin else "", "trs_version": pin.version if pin else ""},
        )
        path = write_task(task, out_dir / f"{entry.client_id}.task.json")
        click.echo(task.render())
        click.echo(f"  -> {path}\n")


@ga4gh_tes.command("service-info")
@click.option("--url", "base_url", default=None, help="Defaults to ga4gh.tes.url.")
@click.option("--federation", "federation_path", default=None, help=_FEDERATION_HELP)
def ga4gh_tes_service_info(base_url: str | None, federation_path: str | None) -> None:
    """Ask a TES service what it is. The cheapest proof that the URL and token work."""
    from appfl_bio_suite.core.ga4gh.tes import TesClient, TesError

    if base_url is None:
        fed = _load(federation_path)
        base_url = fed.tes_service().url
        if not base_url:
            raise click.ClickException("no `ga4gh.tes.url` in the federation config.")
    try:
        info = TesClient(base_url).service_info()
    except TesError as exc:
        raise click.ClickException(str(exc)) from exc
    for key in ("id", "name", "type", "version", "environment"):
        if key in info:
            click.echo(f"  {key}: {info[key]}")


# ---------------------------------------------------------------------------
# site-stage
# ---------------------------------------------------------------------------


@main.command(
    "site-stage",
    context_settings={"ignore_unknown_options": True},
)
@click.argument("experiment", type=click.Choice(experiment_names()), default="fine-mapping")
@click.option("--config", "config_path", required=True, type=click.Path(path_type=Path))
@click.option("--data-dir", required=True, type=click.Path(path_type=Path))
@click.option("--out", "out_dir", type=click.Path(path_type=Path), default=Path("."))
@click.option("--out-name", default="aggregates", show_default=True)
def site_stage_cmd(
    experiment: str, config_path: Path, data_dir: Path, out_dir: Path, out_name: str
) -> None:
    """Compute ONE site's aggregates from its own bundle, and write them to a file.

    The site half of the experiment, as a command. This is what a TES executor runs and
    what the CWL descriptor wraps, and it is the same dataset and trainer classes APPFL
    ships to a Globus Compute worker -- not a reimplementation.

    \b
    It enforces the bundle's own data use terms before reading anything. A study the
    terms do not permit exits 77 with nothing read.
    """
    import importlib

    spec = get_spec(experiment)
    module_name = f"appfl_bio_suite.experiments.{spec.package}.site_stage"
    try:
        module = importlib.import_module(module_name)
    except ImportError as exc:
        raise click.ClickException(
            f"experiment '{experiment}' has no command-line site stage ({module_name}).\n"
            "Only fine-mapping has one: FLamby's site stage is a multi-round training "
            "loop with no single-shot form, and the GWAS one has not been given a CLI."
        ) from exc

    sys.exit(
        module.main(
            [
                experiment,
                "--config",
                str(config_path),
                "--data-dir",
                str(data_dir),
                "--out",
                str(out_dir),
                "--out-name",
                out_name,
            ]
        )
    )


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
        click.echo(f"  install: pip install '{install_spec(','.join(spec.partner_extras))}'")
        for line in _wrap(spec.summary, 74):
            click.echo(f"  {line}")


def _wrap(text: str, width: int) -> list[str]:
    import textwrap

    return textwrap.wrap(text, width=width)


if __name__ == "__main__":
    main()
