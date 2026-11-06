"""Guided partner onboarding; endpoint activation remains an explicit backend step.

The current backend is MEP. The single-user Globus pattern is being developed separately;
joining must not silently select it or infer it from the invoking user's privileges.
"""

from __future__ import annotations

import json
import shlex
import sys
from pathlib import Path
from uuid import UUID

import click
import yaml
from pydantic import ValidationError

from appfl_bio_suite.core import onboarding as api
from appfl_bio_suite.core.onboarding_cli import errors


def read_mapping(path: Path) -> dict:
    values = yaml.safe_load(path.read_text())
    if not isinstance(values, dict):
        raise ValueError(f"{path} must contain a YAML mapping")
    return values


def cluster_settings(path: Path, scheduler: str) -> dict:
    if not path.exists():
        return {}
    values = read_mapping(path)
    allowed = {"scheduler", "account", "partition", "queue", "interface"}
    if set(values) - allowed:
        raise ValueError(f"Unknown cluster profile fields: {sorted(set(values) - allowed)}")
    if values.pop("scheduler", None) != scheduler:
        raise ValueError(f"Profile {path} uses a different scheduler; select another --profile")
    return values


def prompt_value(field: str, default):
    labels = {
        "account": "Scheduler account",
        "partition": "Slurm partition",
        "queue": "PBS queue",
        "data_dir": "Data directory (absolute path)",
        "work_dir": "Work directory (absolute path)",
    }
    while True:
        value = click.prompt(labels[field], default=str(default) if default else None)
        # Validate this individual answer immediately, so only this question repeats.
        try:
            candidate = {
                "account": "unused",
                "data_dir": "/unused",
                "work_dir": "/unused",
                field: value,
            }
            api.SiteConfig.model_validate(candidate)
            return value
        except ValidationError as exc:
            click.echo(f"Invalid {field}: {exc.errors()[0]['msg']}", err=True)


def collect_settings(bundle, profile_values, local_values, *, interactive):
    supplied = {**profile_values, **local_values}
    unknown = set(supplied) - set(api.SiteConfig.model_fields)
    if unknown:
        raise ValueError(f"Unknown local fields: {sorted(unknown)}")
    values = {**bundle.defaults, **supplied}
    slot = "partition" if bundle.scheduler == "slurm" else "queue"
    for field in ("account", slot, "data_dir", "work_dir"):
        # An explicit file/profile is authoritative. Coordinator defaults are suggestions
        # during the guided flow, rather than assumptions about a partner's filesystem.
        if interactive and not supplied.get(field):
            values[field] = prompt_value(field, values.get(field))
    site = api.SiteConfig.model_validate(values)
    other = "queue" if slot == "partition" else "partition"
    if not getattr(site, slot) or getattr(site, other):
        raise ValueError(f"{bundle.scheduler} needs {slot}; do not set {other}")
    return site


def command_line(*args) -> str:
    return shlex.join(str(arg) for arg in args)


def mep_handoff(state, bundle_path, spec):
    executable = state / "runtime/bin/appfl-bio-suite"
    click.echo(
        f"Prepared. Endpoint mode: multi-user (MEP); worker account: {spec.service_account}."
    )
    click.echo("Administrator: ensure that account has scheduler and directory access, then run:")
    click.echo(command_line(executable, "site", "activate", "--state", state))
    click.echo("As that same administrator, obtain the endpoint UUID:")
    click.echo(command_line(state / "runtime/bin/globus-compute-endpoint", "list"))
    click.echo(f"Then, as {spec.service_account}, finish registration:")
    click.echo(
        command_line(
            executable, "join", bundle_path.resolve(), "--state", state, "--endpoint-id", "UUID"
        )
    )


@click.command("join")
@click.argument("bundle", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--local",
    type=click.Path(exists=True, path_type=Path),
    help="Local YAML prepared before joining; supplied values are not prompted again.",
)
@click.option(
    "--profile",
    type=click.Path(path_type=Path),
    help="Reusable cluster YAML; defaults to the app config directory's cluster.yaml.",
)
@click.option("--no-profile", is_flag=True, help="Do not read or save a reusable cluster profile.")
@click.option(
    "--state",
    type=click.Path(path_type=Path),
    help="Shared installation directory; defaults to WORK_DIR/.biosim/BUNDLE_REVISION.",
)
@click.option("--python", default=sys.executable, help="Python 3.12 for the isolated runtime.")
@click.option("--install-source", help="Optional wheel/path override, with experiment extras.")
@click.option("--non-interactive", is_flag=True, help="Validate supplied values without prompting.")
@click.option(
    "--endpoint-id", type=click.UUID, help="Finish verification after endpoint activation."
)
@click.option(
    "--endpoint-mode",
    type=click.Choice(["multi-user"]),
    default="multi-user",
    help="Current supported mode. Single-user support will follow separate validation.",
)
@errors
def join(
    bundle,
    local,
    profile,
    no_profile,
    state,
    python,
    install_source,
    non_interactive,
    endpoint_id: UUID | None,
    endpoint_mode,
):
    """Join an experiment from its ZIP, optionally using settings prepared beforehand."""
    if profile and no_profile:
        raise ValueError("Use either --profile or --no-profile")
    spec = api.read_bundle(bundle)
    click.echo(
        f"Joining {spec.experiment} / {spec.client_id}; coordinator: {spec.coordinator_identity}"
    )
    local_values = read_mapping(local) if local else {}
    saved = None
    if state and (state / "state.json").exists():
        previous, previous_site = api.load_state(state)
        if api.digest(previous) != api.digest(spec):
            raise ValueError("State belongs to another bundle revision; choose a new --state")
        local_values = {**previous_site.model_dump(mode="json"), **local_values}
        saved = json.loads((state / "state.json").read_text())
        if install_source is None:
            previous_source = saved["install_source"]
            if previous_source != (spec.wheel_sha256 or spec.install_spec):
                install_source = previous_source
    profile_path = profile or Path(click.get_app_dir("appfl-bio-suite")) / "cluster.yaml"
    # A resumed join uses its recorded settings, independent of other cluster profiles.
    profile_values = {} if no_profile or saved else cluster_settings(profile_path, spec.scheduler)
    config = collect_settings(spec, profile_values, local_values, interactive=not non_interactive)
    state = (state or config.work_dir / ".biosim" / api.digest(spec)[:16]).resolve()
    if saved and config.model_dump(mode="json") != saved["site"]:
        raise ValueError("Local settings changed; choose a new --state directory")
    if saved and install_source and install_source != saved["install_source"]:
        raise ValueError("Install source changed; choose a new --state directory")
    # A prepared installation can finish registration without another network install.
    prepared = saved and all(
        (state / name).is_file()
        for name in ("runtime/bin/python", "identity_mapping.json", "user_config_template.yaml.j2")
    )
    if not (prepared and endpoint_id):
        api.setup_site(
            spec,
            config,
            state,
            python=python,
            install_source=install_source,
            wheel=api.read_wheel(bundle, spec) if not install_source else None,
        )
    (state / "site.yaml").write_text(
        yaml.safe_dump(config.model_dump(mode="json"), sort_keys=False)
    )
    if not no_profile and not saved:
        slot = "partition" if spec.scheduler == "slurm" else "queue"
        values = {
            "scheduler": spec.scheduler,
            "account": config.account,
            slot: getattr(config, slot),
        }
        if config.interface:
            values["interface"] = config.interface
        profile_path.parent.mkdir(parents=True, exist_ok=True)
        profile_path.write_text(yaml.safe_dump(values, sort_keys=False))
        click.echo(
            f"Cluster settings saved to {profile_path}. Data/work paths stay with this experiment."
        )
    click.echo(f"Local settings: {state / 'site.yaml'}")
    if endpoint_id:
        api.verify_site(state, str(endpoint_id))
        click.echo(f"Locally verified. Send {state / 'registration.json'} to the coordinator.")
        click.echo("The coordinator's remote probe is still required before readiness.")
    else:
        mep_handoff(state, bundle, spec)
