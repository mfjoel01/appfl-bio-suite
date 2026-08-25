"""Build APPFL configs from a federation config, and launch runs.

Two responsibilities, deliberately separable:

1. **Resolve** a federation config plus an experiment's shipped templates into the exact
   server and client config APPFL expects. This is pure and testable, which matters --
   it is what the structural-parity test compares against the pre-migration tree.
2. **Launch** a run, after a preflight, with the failure modes named up front.

WHERE PATHS RESOLVE -- the distinction that has caused the most partner confusion
---------------------------------------------------------------------------------
APPFL's Globus Compute communicator reads dataset and trainer files **on the coordinating
driver**, inlines their source, and ships the text to the worker::

    with open(client_config.data_configs.dataset_path) as file:
        client_config.data_configs.dataset_source = file.read()
    del client_config.data_configs.dataset_path

So ``dataset_path`` and ``trainer_path`` are paths on *your* machine, for every client,
including partners on other continents. A partner is never asked for one, and a setup
guide that asks for one is wrong.

The opposite holds for ``data_dir``, ``output_dir``, and every logging path: those are
used on the worker, so they are absolute paths on the partner's own cluster. Relative
paths there resolve inside the endpoint's task working directory, which is not where
anyone expects.

This module keeps the two straight so that neither the config author nor the partner has
to.
"""

from __future__ import annotations

import os
import shlex
import subprocess
import sys
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from appfl_bio_suite.core.config import Federation
from appfl_bio_suite.core.experiments import ExperimentSpec, get_spec

__all__ = [
    "ResolvedRun",
    "build_client_configs",
    "build_server_config",
    "resolve_run",
    "write_run_configs",
    "launch",
    "PlaceholderError",
]

# Tokens that mean a config was templated and never filled in. Launching with one of
# these present wastes a scheduler queue wait to discover a typo.
_PLACEHOLDER_TOKENS = ("REPLACE_WITH", "<your", "<YOUR", "TODO_FILL", "CHANGEME")


class PlaceholderError(RuntimeError):
    """A config still contains an unfilled placeholder."""


@dataclass
class ResolvedRun:
    """A run, fully resolved, ready to hand to APPFL."""

    experiment: str
    server_config: dict[str, Any]
    client_configs: list[dict[str, Any]]
    spec: ExperimentSpec

    @property
    def num_clients(self) -> int:
        return len(self.client_configs)


def _load_template(spec: ExperimentSpec, name: str) -> dict[str, Any]:
    path = spec.configs_path / name
    if not path.is_file():
        available = sorted(p.name for p in spec.configs_path.glob("*.yaml"))
        raise FileNotFoundError(
            f"no config template '{name}' for experiment '{spec.name}'.\n"
            f"Looked in: {spec.configs_path}\n"
            f"Available: {', '.join(available) or 'none'}"
        )
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _shipped_module_path(spec: ExperimentSpec, module: str) -> str:
    """Absolute path, on the DRIVER, to a module whose source gets shipped.

    Absolute because APPFL resolves it relative to the driver's working directory, and
    depending on where a launcher happened to `cd` is exactly the kind of fragility that
    turns into "works for me".
    """
    path = spec.package_path / f"{module}.py"
    if not path.is_file():
        raise FileNotFoundError(
            f"experiment '{spec.name}' declares shipped module '{module}', but "
            f"{path} does not exist. Either create it or remove it from the registry's "
            "shipped_modules."
        )
    return str(path)


def build_client_configs(
    federation: Federation, experiment: str, spec: ExperimentSpec | None = None
) -> list[dict[str, Any]]:
    """One APPFL client config per participating site, fully resolved.

    Every value comes from federation.yaml. Nothing here is a literal, which is what lets
    a different coordinator run the same command and get correct configs for their own
    federation.
    """
    spec = spec or get_spec(experiment)
    exp = federation.experiment(experiment)
    template = _load_template(spec, "clients.template.yaml")
    defaults = template.get("client_defaults", {}) or {}

    configs: list[dict[str, Any]] = []
    for entry in exp.sites:
        # Deep copy, not dict(): a shallow copy leaves nested values (optim_args, and
        # anything else a template nests) shared between clients. PyYAML then emits them
        # as anchors and aliases -- `&id001` / `*id001` -- which is valid YAML but reads
        # as noise in a config a partner is asked to check, and makes a hand edit to one
        # client silently affect another.
        train = deepcopy(defaults.get("train_configs", {}))
        # Worker-side: absolute paths on the partner's own cluster.
        train["logging_output_dirname"] = entry.output_dir
        train["logging_output_filename"] = f"appfl_{entry.client_id}_{experiment}"

        data_kwargs = deepcopy(defaults.get("dataset_kwargs", {}))
        if experiment == "flamby-heart-disease":
            data_kwargs["dataset"] = exp.dataset
            data_kwargs["num_clients"] = exp.num_clients
            data_kwargs["client_id"] = entry.center
        elif experiment == "gwas":
            data_kwargs["data_dir"] = entry.data_dir
            data_kwargs["site_id"] = entry.client_id
            train["trainer_output_dirname"] = entry.output_dir
            # Federation-wide analysis settings, applied here rather than left to the
            # template. Both MUST be identical at every site: the aggregator refuses to
            # combine payloads computed over different variant sets, and two sites
            # disagreeing about significance would produce two different hits tables for
            # one federation. Setting them from the one config that describes the whole
            # federation is what makes "identical at every site" structural instead of a
            # thing the coordinator has to remember.
            if exp.variant_scaling is not None:
                train["variant_scaling"] = exp.variant_scaling
            if exp.hit_p_threshold is not None:
                train["hit_p_threshold"] = exp.hit_p_threshold
        elif experiment == "fine-mapping":
            data_kwargs["data_dir"] = entry.data_dir
            data_kwargs["site_id"] = entry.client_id
            train["trainer_output_dirname"] = entry.output_dir
            # Same argument as the GWAS settings above, with sharper teeth. The
            # coordinator pools one locus's second moments across sites; two sites on
            # different locus shards would each contribute to loci the other skipped,
            # and every affected locus would be fine-mapped on a smaller cohort than its
            # reported `n`. Nothing downstream can detect that, so the only defence is
            # for the values to come from the one config describing the whole federation.
            for field, value in (
                ("locus_n_shards", exp.locus_n_shards),
                ("locus_shard_index", exp.locus_shard_index),
                ("locus_limit", exp.locus_limit),
                ("instance_limit", exp.instance_limit),
                ("uplink_gram_dtype", exp.uplink_gram_dtype),
            ):
                if value is not None:
                    train[field] = value
            if exp.pops:
                train["pops"] = list(exp.pops)
        else:
            # A registered experiment with no branch here would get a client config
            # carrying no data assignment at all, and would fail on the worker with a
            # loader error that named nothing useful. Failing here names the fix.
            raise NotImplementedError(
                f"build_client_configs has no branch for experiment '{experiment}'.\n"
                "Add one that maps its federation.yaml fields onto train_configs and "
                "dataset_kwargs. Registering an experiment is not enough: the registry "
                "says what exists, this says how its per-site settings reach APPFL."
            )

        client: dict[str, Any] = {
            "endpoint_id": entry.endpoint_uuid,
            "client_id": entry.client_id,
            "train_configs": train,
            "data_configs": {
                # DRIVER-side path. APPFL reads this file here and ships its source.
                "dataset_path": _shipped_module_path(spec, "dataset"),
                "dataset_name": defaults.get("dataset_name", "get_dataset"),
                "dataset_kwargs": data_kwargs,
            },
        }
        configs.append(client)
    return configs


def build_server_config(
    federation: Federation,
    experiment: str,
    variant: str = "default",
    spec: ExperimentSpec | None = None,
) -> dict[str, Any]:
    """The APPFL server config, with federation-derived values applied over a template."""
    spec = spec or get_spec(experiment)
    exp = federation.experiment(experiment)
    filename = "server.yaml" if variant == "default" else f"server.{variant}.yaml"
    config = _load_template(spec, filename)

    server = config.setdefault("server_configs", {})
    server["num_clients"] = len(exp.sites)

    if exp.rounds is not None:
        server["num_global_epochs"] = exp.rounds
    if exp.client_weights_mode is not None:
        server.setdefault("aggregator_kwargs", {})["client_weights_mode"] = (
            exp.client_weights_mode
        )
    # The pooled hits table must use the same threshold the sites did, so it comes from
    # the same place their client configs get it from.
    if exp.hit_p_threshold is not None:
        server.setdefault("aggregator_kwargs", {})["hit_p_threshold"] = exp.hit_p_threshold

    # Fine-mapping's aggregator settings. The answer key in particular belongs in
    # federation.yaml rather than in the shipped server config: whether a run can be
    # scored against ground truth is a fact about that federation's data, not about the
    # experiment. A real federation leaves it unset and gets credible sets it cannot
    # score, which is correct.
    for field, key in (
        ("causal_manifest", "causal_manifest"),
        ("susiex_binary", "susiex_binary"),
        ("credible_set_level", "level"),
        ("pval_thresh", "pval_thresh"),
        ("maf", "maf"),
    ):
        value = getattr(exp, field, None)
        if value is not None:
            server.setdefault("aggregator_kwargs", {})[key] = value

    # Fill in every `*_path` the experiment declares. These are DRIVER-side absolute
    # paths into the installed package: APPFL reads each file here and ships its source.
    # Absolute rather than relative because APPFL resolves them against the driver's
    # working directory, so a relative path silently depends on where the process was
    # started from.
    _resolve_config_paths(config, spec)
    return config


# Where each `*_path` key lives in the config tree. Server-level keys sit under
# `server_configs`; the rest are client-side.
_PATH_KEY_LOCATION = {
    "model_path": ("client_configs", "model_configs"),
    "loss_fn_path": ("client_configs", "train_configs"),
    "metric_path": ("client_configs", "train_configs"),
    "trainer_path": ("client_configs", "train_configs"),
    "aggregator_path": ("server_configs",),
}


def _resolve_config_paths(config: dict[str, Any], spec: ExperimentSpec) -> None:
    """Replace declared `*_path: null` entries with absolute package paths."""
    for key, module in spec.config_path_modules.items():
        location = _PATH_KEY_LOCATION.get(key)
        if location is None:
            raise KeyError(
                f"experiment '{spec.name}' declares config_path_modules['{key}'], but "
                f"launch.py does not know where '{key}' belongs in the config tree. "
                "Add it to _PATH_KEY_LOCATION."
            )
        node = config
        for part in location:
            node = node.setdefault(part, {})
        node[key] = _shipped_module_path(spec, module)


def resolve_run(
    federation: Federation, experiment: str, variant: str = "default"
) -> ResolvedRun:
    spec = get_spec(experiment)
    exp = federation.experiment(experiment)
    if not exp.sites:
        raise RuntimeError(
            f"experiment '{experiment}' has no sites in {federation._where()}.\n"
            "Add at least one under experiments.{experiment}.sites before running."
        )
    return ResolvedRun(
        experiment=experiment,
        server_config=build_server_config(federation, experiment, variant, spec),
        client_configs=build_client_configs(federation, experiment, spec),
        spec=spec,
    )


def check_no_placeholders(config: Any, where: str) -> None:
    """Refuse to launch with an unfilled placeholder.

    Inherited from a launcher script that did this with grep, and worth keeping: the
    cost of catching it here is nothing, and the cost of catching it after a scheduler
    queue wait is a wasted allocation and a confusing failure.

    Comments are not consulted, which the original grep-based version got wrong -- a
    comment mentioning the placeholder token tripped the guard.
    """
    text = yaml.safe_dump(config)
    for token in _PLACEHOLDER_TOKENS:
        if token in text:
            raise PlaceholderError(
                f"{where} still contains the placeholder '{token}'.\n"
                "Fill in the real value in your federation.yaml before launching."
            )


def write_run_configs(run: ResolvedRun, out_dir: Path) -> tuple[Path, Path]:
    """Write resolved configs to disk and return (server_path, client_path)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    server_path = out_dir / f"{run.experiment}_server.yaml"
    client_path = out_dir / f"{run.experiment}_clients.yaml"

    check_no_placeholders(run.server_config, str(server_path))
    check_no_placeholders(run.client_configs, str(client_path))

    server_path.write_text(yaml.safe_dump(run.server_config, sort_keys=False), encoding="utf-8")
    client_path.write_text(
        yaml.safe_dump({"clients": run.client_configs}, sort_keys=False), encoding="utf-8"
    )
    return server_path, client_path


def launch(
    federation: Federation,
    experiment: str,
    variant: str = "default",
    out_dir: Path | None = None,
    dry_run: bool = False,
    driver: str = "globus_compute",
) -> int:
    """Resolve, write, and run. Returns the driver's exit code.

    Preflight runs first and a hard failure stops the launch. That ordering is the
    valuable part of the shell launchers this replaces: it is much cheaper to be told
    the environment is wrong now than to find out from a worker three time zones away.
    """
    from appfl_bio_suite.core.preflight import run_preflight

    report = run_preflight(federation=federation, experiment=experiment, check="env")
    print(report.render(), file=sys.stderr)
    if not report.ok:
        print(
            "\nPreflight failed; not launching. Re-run with `preflight` for detail.",
            file=sys.stderr,
        )
        return 1

    run = resolve_run(federation, experiment, variant)
    out_dir = out_dir or Path("local/configs/generated")
    server_path, client_path = write_run_configs(run, out_dir)

    print(f"\nResolved {run.num_clients} client(s) for '{experiment}':", file=sys.stderr)
    for client in run.client_configs:
        print(f"  {client['client_id']:10} -> {client['endpoint_id']}", file=sys.stderr)
    print(f"\nserver config: {server_path}\nclient config: {client_path}", file=sys.stderr)

    if dry_run:
        print("\n--dry-run: stopping before launch.", file=sys.stderr)
        return 0

    cmd = _driver_command(driver, server_path, client_path, run.server_config)
    print(f"\nlaunching: {shlex.join(cmd)}\n", file=sys.stderr)

    env = dict(os.environ)
    # The same caps the workers get. Without them a busy shared login node can refuse to
    # spawn threads during numpy import, and the driver dies before it dispatches.
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")

    return subprocess.call(cmd, env=env)


def _driver_command(
    driver: str,
    server_path: Path,
    client_path: Path,
    server_config: dict[str, Any] | None = None,
) -> list[str]:
    """Build the APPFL driver invocation.

    ``serial`` is the loopback path: coordinator and simulated sites in one process, no
    Globus Compute and no scheduler. It is what lets someone validate an install before
    recruiting a single partner, and it is the only federated path CI can run.
    """
    if driver == "globus_compute":
        cmd = [
            sys.executable,
            "-m",
            "appfl_bio_suite.core.drivers.globus_compute_driver",
            "--server-config",
            str(server_path),
            "--client-config",
            str(client_path),
        ]
        # APPFL's aggregators weight by sample size only if the server was TOLD each
        # client's sample size, and it is told by an extra round trip the driver makes
        # only on request. Without the flag the aggregator falls back to 1/num_clients
        # silently -- no warning, no error, just a differently-weighted result than the
        # one federation.yaml asked for. The serial driver collects sizes unconditionally
        # (they cost nothing in-process), so this only applies here.
        if _wants_sample_size_weighting(server_config):
            cmd.append("--get-sample-size")
        return cmd
    if driver == "serial":
        return [
            sys.executable,
            "-m",
            "appfl_bio_suite.core.drivers.serial_driver",
            "--server-config",
            str(server_path),
            "--client-config",
            str(client_path),
        ]
    raise ValueError(f"unknown driver '{driver}'. Known: globus_compute, serial.")


def _wants_sample_size_weighting(server_config: dict[str, Any] | None) -> bool:
    """True when the resolved server config asks for sample-size client weighting."""
    if not server_config:
        return False
    aggregator_kwargs = (server_config.get("server_configs") or {}).get("aggregator_kwargs") or {}
    return aggregator_kwargs.get("client_weights_mode") == "sample_size"
