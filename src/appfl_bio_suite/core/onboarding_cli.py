"""CLI for coordinator packaging and partner setup."""

from __future__ import annotations

import subprocess
import sys
import zipfile
from functools import wraps
from pathlib import Path

import click

from appfl_bio_suite.core import onboarding as api


def errors(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        try:
            return function(*args, **kwargs)
        except (
            ValueError,
            OSError,
            KeyError,
            RuntimeError,
            subprocess.CalledProcessError,
            zipfile.BadZipFile,
            TimeoutError,
        ) as exc:
            raise click.ClickException(str(exc)) from exc

    return wrapped


@click.group()
def site():
    """Set up a partner from a coordinator bundle and a small local YAML file."""


@site.command("init")
@click.argument("bundle", type=click.Path(exists=True, path_type=Path))
@click.option("--out", type=click.Path(path_type=Path), default="site.yaml")
@errors
def init(bundle, out):
    """Write the local settings to edit; never overwrite an existing file."""
    import yaml

    spec = api.read_bundle(bundle)
    with out.open("x") as handle:
        handle.write(yaml.safe_dump(spec.defaults, sort_keys=False))
    click.echo(f"Edit {out}: scheduler account/queue and local data/work directories.")


@site.command("setup")
@click.argument("bundle", type=click.Path(exists=True, path_type=Path))
@click.option("--local", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--state", required=True, type=click.Path(path_type=Path))
@click.option("--python", default=sys.executable, help="Python 3.12 executable for the runtime.")
@click.option(
    "--install-source", help="Override release requirement with a wheel/path plus extras."
)
@errors
def setup(bundle, local, state, python, install_source):
    """Install the pinned runtime and generate endpoint files; safe to rerun."""
    spec = api.read_bundle(bundle)
    config = api.resolve_site(spec, local)
    destination = api.setup_site(
        spec,
        config,
        state,
        python=python,
        install_source=install_source,
        wheel=api.read_wheel(bundle, spec) if not install_source else None,
    )
    click.echo(f"Runtime and configuration prepared in {destination}")
    click.echo(f"Administrator: create/authorize POSIX account {spec.service_account}, then run:")
    click.echo(f"  appfl-bio-suite site activate --state {destination}")


@site.command("activate")
@click.option("--state", required=True, type=click.Path(exists=True, path_type=Path))
@errors
def activate(state):
    """As root, configure and start the multi-user endpoint using the prepared runtime."""
    endpoint = api.activate_site(state)
    click.echo(f"Endpoint configured at {endpoint}")
    click.echo("Read its UUID with globus-compute-endpoint list as this same administrator.")
    click.echo("As the service account, run site verify --state PATH --endpoint-id UUID.")


@site.command("verify")
@click.option("--state", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--endpoint-id", help="UUID reported by the administrator's endpoint.")
@click.option(
    "--runtime-only", is_flag=True, help="Check runtime/access without producing a receipt."
)
@errors
def verify(state, endpoint_id, runtime_only):
    """As the service account, verify local runtime/access and produce registration.json."""
    if not runtime_only and not endpoint_id:
        raise ValueError("Supply --endpoint-id UUID, or use --runtime-only")
    if runtime_only and endpoint_id:
        raise ValueError("Use either --runtime-only or --endpoint-id")
    api.verify_site(state, endpoint_id)
    click.echo("Local checks passed.")
    if endpoint_id:
        click.echo(f"Send {state.resolve() / 'registration.json'} to the coordinator.")
        click.echo("Remote worker verification is still required before this site is ready.")


@click.group()
def experiment():
    """Package onboarding bundles and verify partner registrations."""


@experiment.command("package")
@click.argument("experiment_name")
@click.option("--federation", required=True, type=click.Path(exists=True, path_type=Path))
@click.option(
    "--wheel",
    type=click.Path(exists=True, path_type=Path),
    help="Include a matching suite wheel so partners need no repository checkout.",
)
@click.option("--site", "sites", multiple=True, help="Site/client ID; omit to package every site.")
@click.option("--out", required=True, type=click.Path(path_type=Path))
@errors
def package(experiment_name, federation, sites, out, wheel):
    """Generate one distributable ZIP per site using an existing federation.yaml."""
    from appfl_bio_suite.core.config import load_federation
    from appfl_bio_suite.core.partner import generate_bundle

    fed = load_federation(federation)
    entries = fed.experiment(experiment_name).sites
    selected = sites or tuple(entry.client_id for entry in entries)
    for key in selected:
        destination = generate_bundle(fed, experiment_name, key, out)
        if wheel:
            api.attach_wheel(destination, wheel)
        # Archive only the generated manifest's allowlist, never stale/extra local files.
        import json

        manifest = json.loads((destination / "bundle_manifest.json").read_text())
        archive = Path(str(destination) + ".zip")
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as zipped:
            for name in [*manifest["files"], "bundle_manifest.json"]:
                zipped.write(destination / name, name)
        click.echo(archive)


@experiment.command("register")
@click.argument("receipt", type=click.Path(exists=True, path_type=Path))
@click.option("--bundle", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--out", required=True, type=click.Path(path_type=Path))
@click.option("--federation", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--updated-federation", required=True, type=click.Path(path_type=Path))
@click.option("--timeout", default=300.0, type=click.FloatRange(min=1))
@errors
def register(receipt, bundle, out, timeout, federation, updated_federation):
    """Check a receipt against its bundle and verify a real remote worker before readiness."""
    import yaml

    spec = api.read_bundle(bundle)
    if updated_federation.exists() or out.exists():
        raise ValueError(
            "Choose new output files; existing registration/config files are preserved"
        )
    if updated_federation.resolve() == out.resolve():
        raise ValueError("Registration and federation outputs must be different files")
    registration = api.Receipt.model_validate_json(receipt.read_text())
    document = api.updated_federation(spec, registration, federation)
    result = api.register_site(spec, receipt, out, timeout)
    updated_federation.parent.mkdir(parents=True, exist_ok=True)
    updated_federation.write_text(yaml.safe_dump(document, sort_keys=False))
    click.echo(f"Ready: {result['client_id']} ({result['endpoint_id']}); saved {out}")
