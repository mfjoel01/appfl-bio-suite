"""Fine-mapping data simulation: pool -> site cohorts -> loci -> phenotypes -> bundles.

COORDINATOR-SIDE ONLY. None of this ever runs on a partner cluster, which is why its
dependencies live in the separate ``finemapping-sim`` extra and why the shipped-module
import rules do not apply here. This code may freely import from ``appfl_bio_suite``.

Four stages, run in order. The middle three are the vendored upstream pipeline, called
unmodified; the first and last are what the suite adds so that the chain runs from
nothing to distributable bundles::

    cohort.py                 step 0 -- obtain or generate the genotype pool
    fedfm/sampling.py         step 1 -- disjoint per-site cohorts, no individual reused
    fedfm/locus_selection.py  step 2 -- LD-divergence scoring, stratified locus choice
    fedfm/phenotype_sim.py    step 3 -- causal variants, effect sizes, phenotypes
    bundler.py                step 4 -- rearrange into one self-contained bundle per site

Step 0 exists for the same reason the GWAS experiment's does: without it, running this
experiment requires a 135 GB download, and a repository that cannot be exercised without
one is a repository nobody checks. See cohort.py for what the substitute is and, more
importantly, what it is not.

Every run writes a manifest recording the scenario, every seed, the suite commit, input
and output checksums, and package versions -- so a result can be traced back to the data
and the code that produced it.

THIS STAGE NEEDS PLINK
----------------------
Step 1 shells out to ``plink --keep --make-bed`` to cut each site's fileset, exactly as
upstream does. Unlike the GWAS simulation, which is pure numpy, this one has a real
external toolchain dependency -- and so does the experiment itself, which shells out to
SuSiEx. ``scripts/fine-mapping/install_plink.sh`` and ``install_susiex.sh`` vendor static
builds of both. Preflight checks for them before a run rather than after a queue wait.
"""

from __future__ import annotations

import logging
import sys
from datetime import UTC, datetime
from pathlib import Path

__all__ = [
    "run_simulation",
    "list_scenarios",
    "load_scenario",
    "cli_entry",
    "SCENARIO_DIR",
]

log = logging.getLogger(__name__)

SCENARIO_DIR = Path(__file__).resolve().parent.parent / "configs" / "simulation"


def list_scenarios() -> dict[str, Path]:
    """Every shipped scenario, by name."""
    if not SCENARIO_DIR.is_dir():
        return {}
    return {p.stem: p for p in sorted(SCENARIO_DIR.glob("*.yaml"))}


def load_scenario(name_or_path: str):
    """Load a scenario by shipped name or by path."""
    from appfl_bio_suite.experiments.fine_mapping.simulation.schema import (
        FineMappingScenario,
    )

    available = list_scenarios()
    if name_or_path in available:
        scenario = FineMappingScenario.load(available[name_or_path])
    else:
        path = Path(name_or_path)
        if not path.is_file():
            raise FileNotFoundError(
                f"no scenario '{name_or_path}'.\n"
                f"Shipped scenarios: {', '.join(available) or 'none'}\n"
                f"(from {SCENARIO_DIR})\n"
                "You can also pass a path to your own scenario YAML."
            )
        scenario = FineMappingScenario.load(path)
    scenario.validate()
    return scenario


def run_simulation(scenario, out_dir: Path, notes: str = ""):
    """Run steps 0 through 4 and write the run manifest. Returns the manifest."""
    from appfl_bio_suite.core.simulation import (
        MANIFEST_FILENAME,
        SimulationScenario,
        SiteAllocation,
        file_checksum,
        record_provenance,
    )
    from appfl_bio_suite.experiments.fine_mapping.fedfm.locus_selection import (
        run_locus_selection,
    )
    from appfl_bio_suite.experiments.fine_mapping.fedfm.phenotype_sim import (
        run_phenotype_sim,
    )
    from appfl_bio_suite.experiments.fine_mapping.fedfm.sampling import run_sampling
    from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import setup_logging
    from appfl_bio_suite.experiments.fine_mapping.simulation.bundler import bundle_sites
    from appfl_bio_suite.experiments.fine_mapping.simulation.cohort import materialize_pool

    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    started = datetime.now(UTC)

    cfg = scenario.to_pipeline_config(out_dir)
    # The vendored modules log through their own logger. Wiring it up here means their
    # progress lands in the run's own directory alongside everything else, rather than in
    # whatever `logs/` happened to be relative to the working directory.
    setup_logging(log_dir=cfg.resolved_path("logs_dir"))

    # Write the resolved pipeline config beside the data, before anything runs.
    #
    # Two audiences. The vendored stages that this suite does NOT wrap -- the centralized
    # SuSiEx baseline, the validation harness, QC -- take a config path on the command
    # line and have no idea what a scenario is; without this file, running them against a
    # simulated package would mean hand-writing one and getting a path wrong. And a
    # reader wanting to know what parameters actually produced a directory can read one
    # file instead of reconstructing scenario plus injected paths.
    #
    # Written first rather than last so it is present even if a later stage fails, which
    # is exactly when someone needs to re-run one stage by hand.
    _write_pipeline_config(scenario, cfg, out_dir)

    # -- step 0 --------------------------------------------------------
    log.info("step 0: genotype pool (provider=%s)", scenario.cohort.provider)
    pool_prefix = materialize_pool(
        scenario.cohort, scenario.chromosome, scenario.pool_dir(out_dir)
    )

    # -- step 1 --------------------------------------------------------
    log.info("step 1: per-site cohorts")
    run_sampling(cfg, extract_plink=True)

    # -- step 2 --------------------------------------------------------
    log.info("step 2: locus selection")
    run_locus_selection(cfg)

    # -- step 3 --------------------------------------------------------
    log.info("step 3: causal architectures and phenotypes")
    run_phenotype_sim(cfg)

    # -- step 4 --------------------------------------------------------
    log.info("step 4: per-site bundles")
    site_ids = list(cfg.sites)
    bundles = bundle_sites(cfg, site_ids, pool_prefix, out_dir)

    # -- manifest ------------------------------------------------------
    core_scenario = SimulationScenario(
        name=scenario.name,
        description=scenario.description,
        seed=int(scenario.pipeline.get("master_seed", 0)),
        sites=[
            SiteAllocation(site_id=site, n_samples=int(cfg.sites[site].n))
            for site in site_ids
        ],
        parameters={"cohort": scenario.cohort.to_dict(), "pipeline": scenario.pipeline},
    )

    inputs = {
        f"pool{ext}": file_checksum(Path(str(pool_prefix) + ext))
        for ext in (".bed", ".bim", ".fam")
    }
    inputs["population_manifest"] = file_checksum(cfg.resolved_path("population_manifest"))

    # Only the bundles and the answer key are checksummed, not every intermediate. The
    # pipeline's scratch (per-site ID files, candidate windows, LD work) is large,
    # uninteresting, and in places not bit-reproducible across PLINK builds -- including
    # it would make `--verify` report drift that says nothing about the data.
    manifest = record_provenance(
        experiment="fine-mapping",
        scenario=core_scenario,
        seeds=scenario.seeds(),
        inputs=inputs,
        output_root=out_dir,
        started_at=started,
        output_patterns=tuple(f"{site}/**/*" for site in site_ids)
        + ("ground_truth/causal_manifest.tsv", "loci/selected_loci.tsv"),
        notes=notes
        or (
            "Cohort provider: "
            + scenario.cohort.provider
            + (
                ". If 'synthetic_hapnest', this pool is a documented substitute for "
                "HAPNEST and does NOT reproduce the published figures -- see "
                "docs/experiments/fine-mapping/DATA.md."
            )
        ),
    )
    manifest.write(out_dir / MANIFEST_FILENAME)
    log.info("manifest: %s", out_dir / MANIFEST_FILENAME)

    log.info(
        "bundles ready: %s",
        ", ".join(f"{site} -> {path}" for site, path in bundles.items()),
    )
    return manifest


PIPELINE_CONFIG_FILENAME = "pipeline_config.yaml"


def causal_manifest_path(out_dir: Path) -> Path:
    """Where the answer key lands. Coordinator-side; never in a bundle."""
    return Path(out_dir).resolve() / "ground_truth" / "causal_manifest.tsv"


def pipeline_config_path(out_dir: Path) -> Path:
    """The resolved upstream config for a run, readable by the vendored stages."""
    return Path(out_dir).resolve() / PIPELINE_CONFIG_FILENAME


def _write_pipeline_config(scenario, cfg, out_dir: Path) -> Path:
    """Serialize the scenario's `pipeline:` block with this run's paths filled in."""
    import yaml

    destination = pipeline_config_path(out_dir)
    body = {"paths": scenario.paths_for(out_dir, scenario.pool_dir(out_dir))}
    body.update(scenario.pipeline)

    header = (
        f"# Resolved pipeline config for the '{scenario.name}' run in this directory.\n"
        "#\n"
        "# GENERATED -- edit the scenario, not this file. It is the scenario's `pipeline:`\n"
        "# block with `paths:` filled in from --out, which is exactly what the vendored\n"
        "# stages expect on their --config argument:\n"
        "#\n"
        "#   python scripts/fine-mapping/run_stage.py centralized \\\n"
        f"#       --config {destination}\n"
        "#\n"
        "# Those stages -- the centralized SuSiEx baseline, validation, QC -- are the\n"
        "# comparators the federated run is judged against. They run on the coordinator\n"
        "# and have no APPFL involvement.\n"
    )
    destination.write_text(
        header + yaml.safe_dump(body, sort_keys=False), encoding="utf-8"
    )
    log.info("pipeline config: %s", destination)
    return destination


def cli_entry(
    scenario: str | None = None,
    out_dir: Path | None = None,
    verify_manifest: Path | None = None,
    list_scenarios_flag: bool = False,
    **kwargs,
) -> int:
    """Back ``appfl-bio-suite simulate fine-mapping``. Returns a process exit code."""
    from appfl_bio_suite.core.simulation import verify_against_manifest

    logging.basicConfig(
        level=logging.INFO, format="[%(asctime)s %(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    if list_scenarios_flag or kwargs.get("list_scenarios"):
        available = list_scenarios()
        if not available:
            print(f"no scenarios found in {SCENARIO_DIR}")
            return 1
        print("Available fine-mapping simulation scenarios:\n")
        for name in available:
            try:
                loaded = load_scenario(name)
                print(f"  {name}")
                print(f"    {loaded.description.strip()}")
                print(
                    f"    {loaded.site_count} site(s), {loaded.total_samples:,} "
                    f"individuals, pool provider '{loaded.cohort.provider}'"
                )
            except Exception as exc:  # noqa: BLE001 - listing must not die on one bad file
                print(f"  {name}  (could not load: {exc})")
            print()
        return 0

    if verify_manifest is not None:
        target = Path(out_dir) if out_dir else Path(verify_manifest).parent
        reference = Path(verify_manifest)
        from appfl_bio_suite.core.simulation import RunManifest

        patterns = _output_patterns(RunManifest.load(reference))
        ok, report = verify_against_manifest(reference, target, patterns)
        print(report)
        return 0 if ok else 1

    if scenario is None:
        print(
            "specify --scenario. Available:\n  "
            + "\n  ".join(list_scenarios())
            + "\n\nOr --list-scenarios for detail.",
            file=sys.stderr,
        )
        return 2

    loaded = load_scenario(scenario)
    destination = (
        Path(out_dir) if out_dir else Path("local/data") / f"fine-mapping-{loaded.name}"
    )
    manifest = run_simulation(loaded, destination)

    print(f"\nSimulation complete: {destination}")
    print(f"  scenario   {loaded.name}")
    print(f"  sites      {loaded.site_count} ({', '.join(loaded.site_ids)})")
    print(f"  enrolled   {loaded.total_samples:,} individuals")
    print(f"  answer key {causal_manifest_path(destination)}")
    print(f"  manifest   {destination / 'run_manifest.json'}")
    print(f"  commit     {manifest.suite_commit}")
    print("\nThe answer key stays here. It is NOT in any site bundle -- point the")
    print("aggregator at it with aggregator_kwargs.causal_manifest to score a run.")
    print("\nProve the whole path on this machine, without partners:")
    print(
        f"  appfl-bio-suite run fine-mapping --config loopback --driver serial "
        f"--data-root {destination}"
    )
    print("\nVerify a later rerun against this run with:")
    print(
        f"  appfl-bio-suite simulate fine-mapping --verify "
        f"{destination / 'run_manifest.json'}"
    )
    return 0


def _output_patterns(manifest) -> tuple[str, ...]:
    """Reconstruct the patterns a manifest's outputs were collected under.

    ``--verify`` has to re-checksum the same subset the run recorded, or every
    intermediate the run deliberately skipped is reported as an unexpected extra file.
    The site list comes from the manifest itself, so verifying a three-site run and a
    two-site run each look at the right directories.
    """
    sites = [entry["site_id"] for entry in manifest.scenario.get("sites", [])]
    return tuple(f"{site}/**/*" for site in sites) + (
        "ground_truth/causal_manifest.tsv",
        "loci/selected_loci.tsv",
    )
