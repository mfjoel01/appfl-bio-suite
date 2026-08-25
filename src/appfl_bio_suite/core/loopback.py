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
            # A loopback run is the one case where the coordinator legitimately holds the
            # answer key: it simulated the data seconds ago. Wiring it up automatically is
            # what makes the loopback run prove something -- without it the run completes
            # and reports credible sets it cannot score, which looks identical to success
            # whether or not the statistics are right.
            answer_key = data_root / "ground_truth" / "causal_manifest.tsv"
            if answer_key.is_file():
                extra["causal_manifest"] = str(answer_key)
        per_site = [{"data_dir": str(data_root / cid / "data")} for cid in client_ids]
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

    return Federation.model_validate(
        {
            "schema_version": 1,
            "coordinator": {"identity": _LOOPBACK_IDENTITY},
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

    # Same reasoning as in loopback_federation: locally-simulated data comes with its
    # answer key, and a loopback run that cannot score itself proves only that nothing
    # crashed. Only filled in when the file is actually there, and never overriding a
    # path the coordinator set deliberately.
    if experiment == "fine-mapping" and data_root is not None and exp.causal_manifest is None:
        answer_key = (Path(data_root) / "ground_truth" / "causal_manifest.tsv").resolve()
        if answer_key.is_file():
            object.__setattr__(exp, "causal_manifest", str(answer_key))
