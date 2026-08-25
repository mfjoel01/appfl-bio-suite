"""GWAS data simulation: cohort -> phenotypes -> per-site bundles.

COORDINATOR-SIDE ONLY. None of this ever runs on a partner cluster, which is why its
dependencies live in the separate ``gwas-sim`` extra and why the shipped-module import
rules do not apply here. This code may freely import from ``appfl_bio_suite``.

Three stages, run in order:

    cohort.py       step 0 -- obtain or generate the input genotypes
    phenotypes.py   step 1 -- simulate covariates, phenotypes, polygenic scores
    bundler.py      step 2 -- split into disjoint per-site datasets

Steps 1 and 2 are numerically frozen and verified byte-identical against the
pre-migration pipeline. Step 0 is new: it exists because the original cohort is not
regenerable, and without it a stranger could not run this experiment at all. See
cohort.py and docs/experiments/gwas/DATA.md.

Every run writes a manifest recording the scenario, every seed, the suite commit, input
and output checksums, and package versions -- so a result can be traced back to the data
and the code that produced it.
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
    from appfl_bio_suite.experiments.gwas.simulation.schema import GwasScenario

    available = list_scenarios()
    if name_or_path in available:
        scenario = GwasScenario.load(available[name_or_path])
    else:
        path = Path(name_or_path)
        if not path.is_file():
            raise FileNotFoundError(
                f"no scenario '{name_or_path}'.\n"
                f"Shipped scenarios: {', '.join(available) or 'none'}\n"
                f"(from {SCENARIO_DIR})\n"
                "You can also pass a path to your own scenario YAML."
            )
        scenario = GwasScenario.load(path)
    scenario.validate()
    return scenario


def run_simulation(scenario, out_dir: Path, notes: str = ""):
    """Run steps 0 through 2 and write the run manifest. Returns the manifest."""
    from appfl_bio_suite.core.simulation import (
        MANIFEST_FILENAME,
        SimulationScenario,
        SiteAllocation,
        file_checksum,
        record_provenance,
    )
    from appfl_bio_suite.experiments.gwas.simulation.bundler import bundle_sites, verify_disjoint
    from appfl_bio_suite.experiments.gwas.simulation.cohort import materialize_cohort

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    started = datetime.now(UTC)

    # -- step 0 --------------------------------------------------------
    log.info("step 0: cohort (provider=%s)", scenario.cohort.provider)
    plink_prefix = materialize_cohort(scenario.cohort, out_dir)

    pgs_t2d, pgs_bmi = _resolve_pgs_files(scenario, plink_prefix, out_dir)

    # -- step 1 --------------------------------------------------------
    log.info("step 1: phenotypes")
    from appfl_bio_suite.experiments.gwas.simulation.phenotypes import simulate_phenotypes

    pooled_dir = out_dir / "pooled"
    simulate_phenotypes(
        plink_prefix=plink_prefix,
        pgs_t2d_path=pgs_t2d,
        pgs_bmi_path=pgs_bmi,
        out_dir=pooled_dir,
        params=scenario.phenotypes,
    )

    # -- step 2 --------------------------------------------------------
    log.info("step 2: per-site bundles")
    site_ids = [f"Site{i}" for i in range(1, scenario.site_count + 1)]
    site_dirs = bundle_sites(
        plink_prefix=plink_prefix,
        pooled_dir=pooled_dir,
        sites_dir=out_dir,
        site_sizes=list(scenario.sites),
        split_seed=scenario.split_seed,
        site_ids=site_ids,
        data_sim_scaling=scenario.phenotypes.data_sim_scaling,
        data_sim_scaling_seed=scenario.phenotypes.data_sim_scaling_seed,
    )

    overlaps = verify_disjoint(site_dirs)
    if overlaps:
        raise RuntimeError(
            f"sites share samples, which breaks the premise of the federation: {overlaps}. "
            "Every individual must appear at exactly one site, or the meta-analysis "
            "double-counts them while still reporting the full sample size."
        )

    # -- manifest ------------------------------------------------------
    core_scenario = SimulationScenario(
        name=scenario.name,
        description=scenario.description,
        seed=scenario.split_seed,
        sites=[
            SiteAllocation(site_id=sid, n_samples=size)
            for sid, size in zip(site_ids, scenario.sites, strict=True)
        ],
        parameters={
            "cohort": scenario.cohort.to_dict(),
            "phenotypes": scenario.phenotypes.to_dict(),
        },
    )

    inputs = {
        f"cohort{ext}": file_checksum(Path(str(plink_prefix) + ext))
        for ext in (".bed", ".bim", ".fam")
    }
    inputs["pgs_t2d"] = file_checksum(pgs_t2d)
    inputs["pgs_bmi"] = file_checksum(pgs_bmi)

    manifest = record_provenance(
        experiment="gwas",
        scenario=core_scenario,
        seeds={
            "cohort": scenario.cohort.seed,
            "covariates": scenario.phenotypes.seed_covariates,
            "gwas_replicate": scenario.phenotypes.seed_gwas_replicate,
            "pgs_replicate": scenario.phenotypes.seed_pgs_replicate,
            "data_sim_scaling": scenario.phenotypes.data_sim_scaling_seed,
            "site_split": scenario.split_seed,
        },
        inputs=inputs,
        output_root=out_dir,
        started_at=started,
        notes=notes
        or (
            "Cohort provider: "
            + scenario.cohort.provider
            + ". If 'synthetic_ld', this cohort is a documented substitute and does NOT "
            "reproduce the originally published figures -- see "
            "docs/experiments/gwas/DATA.md."
        ),
    )
    manifest.write(out_dir / MANIFEST_FILENAME)
    log.info("manifest: %s", out_dir / MANIFEST_FILENAME)
    return manifest


def _resolve_pgs_files(scenario, plink_prefix: Path, out_dir: Path) -> tuple[Path, Path]:
    """Locate the two PGS Catalog scoring files, or synthesize matched substitutes.

    Real runs use the published PGS Catalog files. When they are absent -- CI, or someone
    exercising the pipeline before obtaining them -- matched substitutes are generated
    against this cohort's own variants, so the pipeline is runnable end to end. The
    manifest records their checksums either way, so a run using substitutes is never
    mistaken for one using the real weights.
    """
    params = scenario.cohort.to_dict()
    t2d = params.get("pgs_t2d_path")
    bmi = params.get("pgs_bmi_path")
    if t2d and bmi:
        return Path(t2d), Path(bmi)

    beside = Path(str(plink_prefix)).parent
    candidates = {
        "t2d": beside / "PGS003443_hmPOS_GRCh37.txt.gz",
        "bmi": beside / "PGS004994_hmPOS_GRCh37.txt.gz",
    }
    if all(p.is_file() for p in candidates.values()):
        log.info("using PGS Catalog files found beside the cohort")
        return candidates["t2d"], candidates["bmi"]

    log.warning(
        "PGS Catalog scoring files not found beside the cohort; generating matched "
        "substitutes. Results will not correspond to any published score. "
        "See docs/experiments/gwas/DATA.md for how to obtain the real files."
    )
    return _synthesize_pgs(plink_prefix, out_dir / "pgs", scenario.cohort.seed)


def _synthesize_pgs(plink_prefix: Path, out_dir: Path, seed: int) -> tuple[Path, Path]:
    """Generate PGS-Catalog-shaped weight files matched to a cohort's variants."""
    import gzip

    import numpy as np
    import pandas as pd

    out_dir.mkdir(parents=True, exist_ok=True)
    bim = pd.read_csv(
        Path(str(plink_prefix) + ".bim"),
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
        dtype=str,
    )
    rng = np.random.default_rng(seed + 991)
    # Partial overlap, as real catalogs have -- the pipeline's overlap step is
    # load-bearing and a fully-overlapping set would not exercise it.
    keep = np.sort(rng.choice(len(bim), max(1, int(len(bim) * 0.6)), replace=False))
    subset = bim.iloc[keep]

    paths = []
    for pgs_id, scale in (("PGS003443", 0.35), ("PGS004994", 0.30)):
        path = out_dir / f"{pgs_id}_hmPOS_GRCh37.txt.gz"
        weights = rng.normal(0, scale, len(subset))
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write("###PGS CATALOG SCORING FILE - SUBSTITUTE, NOT A PUBLISHED SCORE\n")
            handle.write(f"#pgs_id={pgs_id}_substitute\n")
            handle.write("#genome_build=GRCh37\n")
            handle.write("rsID\thm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n")
            for i, (_, row) in enumerate(subset.iterrows()):
                effect, other = ("G", "A") if i % 2 == 0 else ("A", "G")
                handle.write(
                    f"{row.SNP}\t{row.CHR}\t{row.BP}\t{effect}\t{other}\t{weights[i]:.6f}\n"
                )
        paths.append(path)
    return paths[0], paths[1]


def cli_entry(
    scenario: str | None = None,
    out_dir: Path | None = None,
    verify_manifest: Path | None = None,
    list_scenarios_flag: bool = False,
    **kwargs,
) -> int:
    """Back ``appfl-bio-suite simulate gwas``. Returns a process exit code."""
    from appfl_bio_suite.core.simulation import verify_against_manifest

    logging.basicConfig(
        level=logging.INFO, format="[%(asctime)s %(levelname)s] %(message)s", datefmt="%H:%M:%S"
    )

    if list_scenarios_flag or kwargs.get("list_scenarios"):
        available = list_scenarios()
        if not available:
            print(f"no scenarios found in {SCENARIO_DIR}")
            return 1
        print("Available GWAS simulation scenarios:\n")
        for name in available:
            try:
                loaded = load_scenario(name)
                print(f"  {name}")
                print(f"    {loaded.description}")
                print(
                    f"    {loaded.site_count} site(s), {loaded.total_samples:,} samples, "
                    f"cohort provider '{loaded.cohort.provider}'"
                )
            except Exception as exc:  # noqa: BLE001 - listing must not die on one bad file
                print(f"  {name}  (could not load: {exc})")
            print()
        return 0

    if verify_manifest is not None:
        target = Path(out_dir) if out_dir else Path(verify_manifest).parent
        ok, report = verify_against_manifest(verify_manifest, target)
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
    destination = Path(out_dir) if out_dir else Path("local/data") / f"gwas-{loaded.name}"
    manifest = run_simulation(loaded, destination)

    print(f"\nSimulation complete: {destination}")
    print(f"  scenario   {loaded.name}")
    print(f"  sites      {loaded.site_count} ({', '.join(str(s) for s in loaded.sites)})")
    print(f"  samples    {loaded.total_samples:,}")
    print(f"  manifest   {destination / 'run_manifest.json'}")
    print(f"  commit     {manifest.suite_commit}")
    print("\nVerify a later rerun against this run with:")
    print(f"  appfl-bio-suite simulate gwas --verify {destination / 'run_manifest.json'}")
    return 0
