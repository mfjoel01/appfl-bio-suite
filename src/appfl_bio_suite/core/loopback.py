"""Make a serial (loopback) run work on one machine, with or without a federation.yaml.

WHY THIS EXISTS
---------------
The loopback run is the first thing the README, the new-federation guide and CI all ask
someone to do: *prove your install works, before involving anyone*. That instruction comes
BEFORE the step that creates ``local/federation.yaml``, and deliberately so -- describing a
federation is the thing you do once the install is known good.

So a loopback run cannot require a federation config. It also cannot reuse one unmodified
when there is one: every path in a federation config is a path on somebody else's cluster.
``data_dir`` is where a partner unpacked their bundle and ``output_dir`` is where their
service account can write. Running that config serially on the coordinator's laptop means
trying to ``mkdir /home/gwas_svc/outputs`` locally, which fails with a ``PermissionError``
that has nothing to do with what was being tested.

This module answers both halves: synthesize a federation when none exists, and rewrite the
worker-side paths of one that does so they point at this machine.

Nothing here ever touches a real run. It is reached only from ``run --driver serial``.

THE LOOPBACK RUN EXERCISES THE GA4GH LAYER TOO
----------------------------------------------
``simulate`` writes a DUO profile into every bundle and registers each one with DRS, so a
loopback run that ignored both would be proving less than the real path does -- and would
fail, because a bundle that declares terms is refused by a run that declares no study.

So a synthesized federation carries a data use request, points each site at its own
profile, and names each site's DRS object. The request is truthful about what it is: a
self-run over synthetic data by the person who generated it.

It also means the shipped scenario's *differing* per-site terms are matched for real on
every CI run. Adding population-origins research to the purposes below makes the
``covenant`` site refuse, which is exactly the behaviour that is worth having proven
before a real partner's terms are in play.
"""

from __future__ import annotations

from pathlib import Path

from appfl_bio_suite.core.config import Federation, FederationError

__all__ = ["loopback_federation", "localize_for_loopback", "default_output_root"]

# Structurally required by the schema and meaningless here: the serial driver strips
# `endpoint_id` before constructing a ClientAgent, and nothing in a loopback run contacts
# Globus at all.
_PLACEHOLDER_UUID = "00000000-0000-0000-0000-{:012d}"

# Not a real identity, and it does not need to be -- no partner authorizes anything for a
# run that never leaves this process. Kept obviously fake so it cannot be mistaken for a
# coordinator's own identity if it turns up in a generated config.
_LOOPBACK_IDENTITY = "loopback@localhost"

_DEFAULT_SITE_COUNT = 2


def default_output_root(experiment: str) -> Path:
    """Where a loopback run writes per-site output. Under ``local/``, which is gitignored."""
    return (Path("local/output") / f"{experiment}-loopback" / "sites").resolve()


def _site_output_dir(out_root: Path, client_id: str) -> str:
    return str((Path(out_root) / client_id).resolve())


def _discover_simulated_sites(experiment: str, data_root: Path) -> list[str]:
    """Site directories the simulator wrote under ``--out``.

    Every simulating experiment writes ``<out>/<site>/data/`` and hands the loader the
    ``data`` directory itself, because each one requires its input files to sit directly
    in ``data_dir``. That shared layout is why this does not need to know which
    experiment it is discovering for -- only which one to name in the error.
    """
    if not data_root.is_dir():
        raise FederationError(
            f"--data-root {data_root} does not exist.\n"
            "Generate the per-site data first:\n"
            f"    appfl-bio-suite simulate {experiment} --scenario ci-tiny "
            f"--out {data_root}"
        )
    sites = sorted(child.name for child in data_root.iterdir() if (child / "data").is_dir())
    if not sites:
        raise FederationError(
            f"no simulated sites under {data_root}.\n"
            "Expected one directory per site, each containing a `data/` subdirectory --\n"
            "the layout `simulate` writes. Generate it with:\n"
            f"    appfl-bio-suite simulate {experiment} --scenario ci-tiny "
            f"--out {data_root}"
        )
    return sites


def loopback_federation(
    experiment: str, data_root: Path | None = None, out_root: Path | None = None
) -> Federation:
    """Build a complete, valid federation describing simulated sites on this machine.

    Used when a serial run finds no federation config, so that the documented
    "prove your install works" flow needs nothing to have been filled in first.
    """
    out_root = Path(out_root) if out_root else default_output_root(experiment)

    if experiment in ("gwas", "fine-mapping"):
        if data_root is None:
            raise FederationError(
                f"a loopback {experiment} run needs --data-root: it is where `simulate` "
                "wrote the\nper-site data.\n"
                f"    appfl-bio-suite simulate {experiment} --scenario ci-tiny "
                f"--out /tmp/{experiment}-ci\n"
                f"    appfl-bio-suite run {experiment} --config loopback --driver serial "
                f"\\\n        --data-root /tmp/{experiment}-ci"
            )
        data_root = Path(data_root).resolve()
        client_ids = _discover_simulated_sites(experiment, data_root)
        extra: dict = {}
        if experiment == "fine-mapping":
            extra["ga4gh"] = _loopback_ga4gh(data_root)
            # A loopback run is the one case where the coordinator legitimately holds the
            # answer key: it simulated the data seconds ago. Wiring it up automatically is
            # what makes the loopback run prove something -- without it the run completes
            # and reports credible sets it cannot score, which looks identical to success
            # whether or not the statistics are right.
            answer_key = data_root / "ground_truth" / "causal_manifest.tsv"
            if answer_key.is_file():
                extra["causal_manifest"] = str(answer_key)
            # Same argument as the answer key: the coordinator simulated this data, so it
            # knows which MAF filter the scenario declared. Taking it from there rather
            # than from the shipped loopback template is what keeps the loopback result
            # comparable to `run_stage.py centralized`, which reads the same field.
            maf = _scenario_maf(data_root)
            if maf is not None:
                extra["maf"] = maf
        per_site = [
            {"data_dir": str(data_root / cid / "data"), **_loopback_site_ga4gh(data_root, cid)}
            for cid in client_ids
        ]
    elif experiment == "flamby-heart-disease":
        client_ids = [f"Site{i}" for i in range(1, _DEFAULT_SITE_COUNT + 1)]
        # Distinct centers: two simulated sites training on the same shard would be
        # training twice on the same data, which is not the thing being proved.
        extra = {"dataset": "HeartDisease", "num_clients": 4, "client_weights_mode": "equal"}
        per_site = [{"center": i} for i in range(len(client_ids))]
    else:
        raise FederationError(
            f"no loopback federation is defined for experiment '{experiment}'. "
            "Pass --federation with a config that declares its sites."
        )

    sites = [{"id": cid.lower(), "name": f"Loopback {cid}"} for cid in client_ids]
    experiment_sites = [
        {
            "site": site["id"],
            "client_id": cid,
            "endpoint_uuid": _PLACEHOLDER_UUID.format(index + 1),
            "output_dir": _site_output_dir(out_root, cid),
            **fields,
        }
        for index, (cid, site, fields) in enumerate(zip(client_ids, sites, per_site, strict=True))
    ]

    services = _loopback_drs_service(data_root) if experiment == "fine-mapping" else None

    return Federation.model_validate(
        {
            "schema_version": 1,
            "coordinator": {"identity": _LOOPBACK_IDENTITY},
            **({"ga4gh": services} if services else {}),
            "sites": sites,
            "experiments": {
                experiment: {
                    "enabled": True,
                    "service_account": "loopback",
                    "endpoint_name": f"{experiment}-loopback",
                    "sites": experiment_sites,
                    **extra,
                }
            },
        }
    )


# The study a loopback run declares. Truthful, which matters more here than it looks:
# these attestations are matched against the shipped scenario's per-site terms on every
# CI run, and a request that claimed something false to get past a check would make the
# check meaningless in exactly the place it is cheapest to keep honest.
#
# `DUO:0000038 genetic research` is a subclass of biomedical research, so it satisfies the
# `covenant` site's HMB permission. Adding a non-biomedical purpose here makes that site
# refuse -- which is the point of it being declared rather than assumed.
_LOOPBACK_DATA_USE = {
    "requester": _LOOPBACK_IDENTITY,
    "project": "loopback install validation",
    "institution": "loopback",
    "purposes": ["DUO:0000038"],
    "non_commercial": True,
    "not_for_profit_organisation": True,
    "publication_agreed": True,
    "ethics_approval": "not applicable: synthetic data, no human subjects",
}


def _scenario_maf(data_root: Path) -> float | None:
    """The pooled-MAF filter the simulated scenario declared, if it declared one.

    Read straight out of ``pipeline_config.yaml`` rather than re-validating the whole
    scenario: this runs before APPFL starts, a malformed field should not take the run
    down, and the aggregator's own default is the right fallback.
    """
    config_path = Path(data_root) / "pipeline_config.yaml"
    if not config_path.is_file():
        return None
    import yaml

    try:
        block = (yaml.safe_load(config_path.read_text()) or {}).get("fine_mapping") or {}
        maf = block.get("maf")
        return float(maf) if maf is not None else None
    except (yaml.YAMLError, OSError, TypeError, ValueError):
        return None


def _loopback_drs_service(data_root: Path) -> dict | None:
    """The `ga4gh.drs` block for a simulated data directory, if it has a registry."""
    registry = Path(data_root) / "drs_registry.json"
    if not registry.is_file():
        return None
    import json

    try:
        hostname = json.loads(registry.read_text(encoding="utf-8"))["hostname"]
    except (OSError, ValueError, KeyError):
        return None
    return {"drs": {"hostname": hostname, "registry": str(registry)}}


def _loopback_ga4gh(data_root: Path) -> dict:
    """The experiment-level GA4GH block for a loopback run."""
    return {
        "data_use_request": dict(_LOOPBACK_DATA_USE),
        "enforce_data_use": True,
        # `metadata` rather than `full`: the .bed was written seconds ago by the process
        # that is about to read it, so hashing it proves nothing and costs the whole
        # point of a fast smoke test.
        "verify_bundles": "metadata",
    }


def _loopback_site_ga4gh(data_root: Path, client_id: str) -> dict:
    """A simulated site's own GA4GH fields: its profile, and its DRS object."""
    out: dict = {}
    profile = Path(data_root) / client_id / "data" / "DATA_USE.json"
    if profile.is_file():
        out["data_use_profile"] = str(profile)

    registry_path = Path(data_root) / "drs_registry.json"
    if registry_path.is_file():
        from appfl_bio_suite.core.ga4gh.drs import DrsError, DrsRegistry

        try:
            registry = DrsRegistry.load(registry_path)
            obj = registry.by_name(client_id)
            if obj is not None:
                out["drs_uri"] = obj.self_uri
        except DrsError:
            # A registry that will not load is a reason to run without DRS verification,
            # never a reason not to run: the loopback run's job is to prove the install.
            pass
    return out


def localize_for_loopback(
    federation: Federation,
    experiment: str,
    data_root: Path | None = None,
    out_root: Path | None = None,
) -> None:
    """Repoint a real federation's worker-side paths at this machine, in place.

    Both kinds of path have to move, not just the data:

    * ``data_dir`` -- only when ``--data-root`` is given, and it resolves to
      ``<root>/<client_id>/data`` because that is the layout ``simulate`` writes and the
      loader requires the six input files to sit directly in ``data_dir``.
    * ``output_dir`` -- always. It is an absolute path on a partner's cluster, and a
      serial run executes here, so leaving it alone means the run dies trying to create
      a directory under someone else's home.
    """
    out_root = Path(out_root) if out_root else default_output_root(experiment)
    exp = federation.experiment(experiment)
    for entry in exp.sites:
        if data_root is not None:
            local_data = (Path(data_root) / entry.client_id / "data").resolve()
            object.__setattr__(entry, "data_dir", str(local_data))
        object.__setattr__(entry, "output_dir", _site_output_dir(out_root, entry.client_id))

    # The GA4GH fields point at the simulated data too, for the same reason the paths do:
    # a real federation's `drs_uri` names an object in the coordinator's registry and its
    # `data_use_profile` a copy of a partner's terms, neither of which describes the
    # bundle sitting in --data-root. Only filled in where the coordinator left them unset.
    if experiment == "fine-mapping" and data_root is not None:
        _localize_ga4gh(federation, exp, Path(data_root))

    # Same reasoning as in loopback_federation: locally-simulated data comes with its
    # answer key, and a loopback run that cannot score itself proves only that nothing
    # crashed. Only filled in when the file is actually there, and never overriding a
    # path the coordinator set deliberately.
    if experiment == "fine-mapping" and data_root is not None and exp.causal_manifest is None:
        answer_key = (Path(data_root) / "ground_truth" / "causal_manifest.tsv").resolve()
        if answer_key.is_file():
            object.__setattr__(exp, "causal_manifest", str(answer_key))

    # And the filter the simulated scenario declared, so a loopback run driven from a real
    # federation.yaml still matches `run_stage.py centralized` on the same data. Never
    # overrides a value the coordinator set deliberately.
    if experiment == "fine-mapping" and data_root is not None and exp.maf is None:
        maf = _scenario_maf(Path(data_root))
        if maf is not None:
            object.__setattr__(exp, "maf", maf)


def _localize_ga4gh(federation: Federation, exp, data_root: Path) -> None:
    """Repoint an existing federation's GA4GH fields at locally simulated bundles."""
    from appfl_bio_suite.core.config import ExperimentGA4GH, GA4GHServices

    if exp.ga4gh is None:
        object.__setattr__(exp, "ga4gh", ExperimentGA4GH.model_validate(_loopback_ga4gh(data_root)))

    services = _loopback_drs_service(data_root)
    if services is not None and federation.ga4gh is None:
        object.__setattr__(federation, "ga4gh", GA4GHServices.model_validate(services))

    for entry in exp.sites:
        for key, value in _loopback_site_ga4gh(data_root, entry.client_id).items():
            object.__setattr__(entry, key, value)
