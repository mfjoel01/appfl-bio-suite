"""The simulation split into separately-runnable stages, for sharded HPC runs.

WHY THIS EXISTS ALONGSIDE ``run_simulation``
---------------------------------------------
``appfl-bio-suite simulate fine-mapping`` runs the whole pipeline in one process. That is
the right shape for ``ci-tiny``, for a moderate scenario, and for anyone finding out
whether the thing works -- one command, one directory, one manifest.

It is the wrong shape for the published scenario. That run simulates 100 loci x 15
architectures x 10 replicates over 150,000 individuals; two of its stages are
embarrassingly parallel over loci and the upstream pipeline shards them across compute
nodes for that reason. A single process would take days.

So the same pipeline is also reachable one stage at a time. The vendored stages already
support ``--n-shards``/``--shard-index``/``--merge`` -- that is their own design, not
something added here -- and what was missing was the first and last steps, which belong
to the suite rather than to upstream:

    prepare   materialize the genotype pool and write pipeline_config.yaml, so the
              sharded stages have a config to point at
    bundle    rearrange the finished outputs into per-site bundles and check disjointness

Between them, the sharded run is::

    run_stage.py prepare        --scenario three-site-hapnest --out $RUN
    run_stage.py sampling       --config $RUN/pipeline_config.yaml
    run_stage.py locus-selection --config $RUN/pipeline_config.yaml   # shardable
    run_stage.py phenotypes     --config $RUN/pipeline_config.yaml \\
        --shard-index $i --n-shards $N                                # shardable
    run_stage.py phenotypes     --config $RUN/pipeline_config.yaml --merge --n-shards $N
    run_stage.py bundle         --scenario three-site-hapnest --out $RUN

WHAT YOU GIVE UP BY SHARDING
-----------------------------
The run manifest. ``run_simulation`` writes it because it is the one process that saw the
whole run and can checksum its outputs against the inputs and seeds that produced them.
A stage-by-stage run has no such vantage point, and a manifest assembled after the fact
from whatever happens to be on disk would assert a provenance nobody verified.

``bundle`` therefore writes the manifest only when it can honestly do so -- it re-derives
the scenario, re-checksums the inputs, and records that the run was sharded, with the
shard count in the notes. What it cannot record is that every shard ran the code this
manifest names, so a sharded run's manifest is weaker evidence than a single-process
one's, and says so.
"""

from __future__ import annotations

import argparse
import logging
import sys
from datetime import UTC, datetime
from pathlib import Path

__all__ = ["prepare_main", "bundle_main"]

log = logging.getLogger(__name__)


def _scenario_and_out(args) -> tuple:
    from appfl_bio_suite.experiments.fine_mapping.simulation import load_scenario

    return load_scenario(args.scenario), Path(args.out).resolve()


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--scenario",
        required=True,
        help="scenario name or path, as `simulate --list-scenarios` reports",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="the run directory; the same value every stage of the run is given",
    )


def prepare_main(argv=None) -> int:
    """Materialize the genotype pool and write the resolved pipeline config."""
    from appfl_bio_suite.experiments.fine_mapping.simulation import _write_pipeline_config
    from appfl_bio_suite.experiments.fine_mapping.simulation.cohort import materialize_pool

    parser = argparse.ArgumentParser(
        prog="run_stage.py prepare",
        description=prepare_main.__doc__,
    )
    _add_common(parser)
    args = parser.parse_args(list(argv) if argv is not None else None)

    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s %(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    scenario, out_dir = _scenario_and_out(args)
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = scenario.to_pipeline_config(out_dir)
    prefix = materialize_pool(scenario.cohort, scenario.chromosome, scenario.pool_dir(out_dir))
    destination = _write_pipeline_config(scenario, cfg, out_dir)

    print(f"pool:   {prefix}")
    print(f"config: {destination}")
    print("\nNext, with that config:")
    print(f"  run_stage.py sampling --config {destination}")
    return 0


def bundle_main(argv=None) -> int:
    """Rearrange the finished pipeline outputs into per-site bundles."""
    from appfl_bio_suite.core.simulation import (
        MANIFEST_FILENAME,
        SimulationScenario,
        SiteAllocation,
        file_checksum,
        record_provenance,
    )
    from appfl_bio_suite.experiments.fine_mapping.simulation.bundler import bundle_sites
    from appfl_bio_suite.experiments.fine_mapping.simulation.cohort import POOL_STEM

    parser = argparse.ArgumentParser(
        prog="run_stage.py bundle",
        description=bundle_main.__doc__,
    )
    _add_common(parser)
    parser.add_argument(
        "--n-shards",
        type=int,
        default=1,
        help="how many shards produced this run; recorded in the manifest",
    )
    parser.add_argument(
        "--no-manifest",
        action="store_true",
        help="skip the run manifest (checksumming a large package is not free)",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s %(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    started = datetime.now(UTC)
    scenario, out_dir = _scenario_and_out(args)
    cfg = scenario.to_pipeline_config(out_dir)

    site_ids = list(cfg.sites)
    pool_prefix = scenario.pool_dir(out_dir) / POOL_STEM(scenario.chromosome)
    bundles = bundle_sites(cfg, site_ids, pool_prefix, out_dir)

    for site, path in bundles.items():
        print(f"{site}: {path}")

    if args.no_manifest:
        print("\nno run manifest written (--no-manifest)")
        return 0

    inputs = {
        f"pool{ext}": file_checksum(Path(str(pool_prefix) + ext))
        for ext in (".bed", ".bim", ".fam")
    }
    inputs["population_manifest"] = file_checksum(cfg.resolved_path("population_manifest"))

    manifest = record_provenance(
        experiment="fine-mapping",
        scenario=SimulationScenario(
            name=scenario.name,
            description=scenario.description,
            seed=int(scenario.pipeline.get("master_seed", 0)),
            sites=[
                SiteAllocation(site_id=site, n_samples=int(cfg.sites[site].n)) for site in site_ids
            ],
            parameters={"cohort": scenario.cohort.to_dict(), "pipeline": scenario.pipeline},
        ),
        seeds=scenario.seeds(),
        inputs=inputs,
        output_root=out_dir,
        started_at=started,
        output_patterns=tuple(f"{site}/**/*" for site in site_ids)
        + ("ground_truth/causal_manifest.tsv", "loci/selected_loci.tsv"),
        notes=(
            f"SHARDED RUN across {args.n_shards} shard(s), assembled by `run_stage.py "
            "bundle`. Weaker evidence than a single-process `simulate` manifest: it "
            "records the scenario, the seeds and the output checksums, but nothing here "
            "witnessed each shard, so it cannot attest that every shard ran this code. "
            f"Cohort provider: {scenario.cohort.provider}."
        ),
    )
    manifest.write(out_dir / MANIFEST_FILENAME)
    print(f"\nmanifest: {out_dir / MANIFEST_FILENAME}")
    print(f"answer key: {out_dir / 'ground_truth' / 'causal_manifest.tsv'} (NOT bundled)")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(prepare_main())
