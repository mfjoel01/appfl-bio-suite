"""The experiment registry -- one place that knows what experiments exist.

Every experiment-shaped question in the suite resolves here: where its configs live,
which of its modules get shipped to workers, which docs it must have, whether it has a
simulation stage. Nothing else hardcodes an experiment name.

That is what makes adding the third experiment a fill-in rather than a restructure. To
add one you write a :class:`ExperimentSpec` entry, create the directories it names, and
write the four documents; the CLI, preflight, bundle generator, and test suite all pick
it up without modification. If adding an experiment ever requires editing anything other
than this registry and the experiment's own package, that is a design bug.

WHAT THE THIRD EXPERIMENT ACTUALLY COST
---------------------------------------
``fine-mapping`` was declared here with ``implemented=False`` while the suite shipped two
experiments, precisely so that the claim above could be tested rather than asserted.
Implementing it flipped one flag and filled in one package -- and turned up four places
that had quietly hardcoded ``gwas`` where they meant "an experiment with this property":

* ``cli.py``'s ``simulate`` imported the GWAS simulation package by name.
* ``loopback.py`` knew how to synthesize a federation for two named experiments.
* ``preflight.py``'s data check tested ``if name != "gwas"``.
* ``launch.py``'s client builder had a per-experiment branch with no fallback.

The first three are now dispatched from the registry or from the experiment's own
package. The fourth stays an explicit per-experiment branch, because mapping
``federation.yaml`` onto APPFL's config tree really is per-experiment -- but it now
raises on an experiment it does not know instead of silently emitting a client config
with no data assignment, which is how that branch would otherwise fail: on a worker,
minutes in.

None of the four required a restructure, which is the outcome the empty slot existed to
predict. But "no special cases were needed" would have been the wrong summary, and the
seams are worth naming so the fourth experiment does not have to rediscover them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "ExperimentSpec",
    "REGISTRY",
    "get_spec",
    "implemented_experiments",
    "experiment_names",
    "package_root",
    "repo_root",
    "check_shipped_configs",
]


@dataclass(frozen=True)
class ExperimentSpec:
    """Everything the shared machinery needs to know about one experiment."""

    # Canonical hyphenated name. Used on the CLI, in federation.yaml, and as the
    # directory name under docs/experiments/ and docs/partner/experiments/.
    name: str

    # Python package name under appfl_bio_suite.experiments. Underscored, because it is
    # an identifier; the two differ and conflating them is a common slip.
    package: str

    title: str
    summary: str

    # Modules whose SOURCE is read on the coordinator and shipped to workers as text.
    # These are subject to the self-containment rule: they may not import anything the
    # partner has not installed, and in particular may not import from appfl_bio_suite.
    # tests/test_shipped_modules.py enforces this against exactly this list.
    shipped_modules: tuple[str, ...] = ()

    # Does this experiment generate its own data? FLamby does not -- partners download a
    # public dataset -- so the pipeline must tolerate an experiment with no simulation
    # stage at all rather than assuming every experiment has one.
    has_simulation: bool = False

    # False means the slots exist and the docs explain the plan, but there is no code.
    implemented: bool = True

    # Extras a coordinator needs, and extras a partner needs. They differ: simulation
    # dependencies are coordinator-only and must never land on a partner.
    partner_extras: tuple[str, ...] = ()
    coordinator_extras: tuple[str, ...] = ()

    # External executables the COORDINATOR must be able to find, beyond anything pip can
    # install. Only fine-mapping has any: it shells out to the SuSiEx C++ CLI, and its
    # simulation shells out to PLINK. A missing binary is a hard failure discovered
    # partway through a run otherwise -- after the site aggregates have been computed and
    # transferred, which is the expensive part.
    #
    # Deliberately not a partner concern. Partners compute second moments in numpy; the
    # toolchain lives entirely on the coordinator.
    required_binaries: tuple[str, ...] = ()

    # Which shipped module backs each `*_path` key in the configs.
    #
    # Shipped configs carry `model_path: null` and so on, and these are filled in with
    # absolute paths at run time. Absolute because APPFL resolves them relative to the
    # driver's working directory, and depending on where a launcher happened to `cd` is
    # how a config becomes "works on my machine".
    #
    # Declared rather than inferred from the key name: loss and metric live alongside the
    # model rather than in files of their own, and a naming convention that has to encode
    # that is a convention waiting to be broken.
    config_path_modules: dict[str, str] = field(default_factory=dict)

    @property
    def package_path(self) -> Path:
        return package_root() / "experiments" / self.package

    @property
    def configs_path(self) -> Path:
        return self.package_path / "configs"

    @property
    def docs_path(self) -> Path:
        return repo_root() / "docs" / "experiments" / self.name

    @property
    def partner_doc_path(self) -> Path:
        return repo_root() / "docs" / "partner" / "experiments" / f"{self.name}.md"

    def shipped_module_paths(self) -> list[Path]:
        return [self.package_path / f"{m}.py" for m in self.shipped_modules]


def package_root() -> Path:
    """Directory of the installed ``appfl_bio_suite`` package."""
    return Path(__file__).resolve().parent.parent


def repo_root() -> Path:
    """Repository root when running from a source checkout.

    Docs are not packaged, so anything that reads them only works from a checkout. That
    is fine -- the audiences for docs (coordinators, reviewers, partners) all have one --
    but callers must handle the directory being absent in an installed-only environment.
    """
    return package_root().parent.parent


REGISTRY: dict[str, ExperimentSpec] = {
    "flamby-heart-disease": ExperimentSpec(
        name="flamby-heart-disease",
        package="flamby_heart_disease",
        title="FLamby Fed-Heart-Disease",
        summary=(
            "Multi-round FedAvg on a small public tabular benchmark with natural "
            "hospital-level client splits. The scientific claim is about the "
            "infrastructure, not the model: every failure is an infrastructure failure, "
            "never a compute failure."
        ),
        shipped_modules=("dataset", "model"),
        has_simulation=False,
        implemented=True,
        partner_extras=("flamby",),
        coordinator_extras=("flamby",),
        config_path_modules={
            "model_path": "model",
            # Loss and metric live in model.py alongside the network they belong to.
            "loss_fn_path": "model",
            "metric_path": "model",
        },
    ),
    "gwas": ExperimentSpec(
        name="gwas",
        package="gwas",
        title="Federated GWAS with polygenic-score evaluation",
        summary=(
            "Single-round summary-statistic federated learning. Each site runs a "
            "complete local GWAS and returns per-variant effect sizes, standard errors "
            "and allele frequencies; the server combines them by inverse-variance "
            "fixed-effect meta-analysis. Genotypes never leave the site."
        ),
        # NOT "aggregator": APPFL's ServerAgent loads aggregator_path locally and never
        # ships it. It runs on the coordinator, so it may import freely from this package.
        # Verified against appfl/agent/server.py -- model_path, loss_fn_path, metric_path
        # and trainer_path are all read server-side and shipped as *_source, while
        # aggregator_path is not.
        shipped_modules=("dataset", "trainer"),
        has_simulation=True,
        implemented=True,
        partner_extras=("gwas",),
        coordinator_extras=("gwas", "gwas-sim"),
        config_path_modules={
            "trainer_path": "trainer",
            "aggregator_path": "aggregator",
        },
    ),
    "fine-mapping": ExperimentSpec(
        name="fine-mapping",
        package="fine_mapping",
        title="Federated cross-ancestry statistical fine-mapping",
        summary=(
            "Single-round aggregate federated learning. Each site returns raw "
            "second-moment aggregates over its own individuals at a locus; the server "
            "sums them, standardizes once against the pooled moments, and runs SuSiEx. "
            "Unlike meta-analysis this does not approximate a pooled fit -- it equals "
            "it, and the test suite holds it to that."
        ),
        # Same two as GWAS, and for the same reason: APPFL's ServerAgent reads
        # aggregator_path locally and never ships it, while dataset_path and
        # trainer_path are read on the driver and shipped as source.
        shipped_modules=("dataset", "trainer"),
        has_simulation=True,
        implemented=True,
        partner_extras=("finemapping",),
        coordinator_extras=("finemapping", "finemapping-sim"),
        # SuSiEx does the fitting; PLINK cuts the per-site filesets during simulation.
        # scripts/fine-mapping/install_susiex.sh and install_plink.sh vendor both.
        required_binaries=("SuSiEx", "plink"),
        config_path_modules={
            "trainer_path": "trainer",
            "aggregator_path": "aggregator",
        },
    ),
}


def experiment_names(implemented_only: bool = False) -> list[str]:
    return [n for n, s in REGISTRY.items() if s.implemented or not implemented_only]


def implemented_experiments() -> dict[str, ExperimentSpec]:
    return {n: s for n, s in REGISTRY.items() if s.implemented}


def get_spec(name: str) -> ExperimentSpec:
    """Look up an experiment, with a message that lists the alternatives."""
    try:
        return REGISTRY[name]
    except KeyError:
        known = ", ".join(REGISTRY)
        # Underscores for hyphens is the single most common mistake here, because the
        # package name really does use them.
        suggestion = name.replace("_", "-")
        hint = f" Did you mean '{suggestion}'?" if suggestion in REGISTRY else ""
        raise KeyError(f"unknown experiment '{name}'. Known: {known}.{hint}") from None


def check_shipped_configs(report, experiment: str | None = None) -> None:
    """Assert every shipped config file parses. Used by preflight and by CI.

    A config that does not parse is a failure that costs a scheduler queue wait to
    discover if it is found at launch time instead of here.
    """
    import yaml

    from appfl_bio_suite.core.preflight import Level

    specs = (
        [get_spec(experiment)]
        if experiment
        else [s for s in REGISTRY.values() if s.implemented]
    )

    for spec in specs:
        configs = spec.configs_path
        if not configs.is_dir():
            level = Level.SKIP if not spec.implemented else Level.FAIL
            report.add(
                f"configs [{spec.name}]",
                level,
                f"no configs directory at {configs}",
                "Every implemented experiment ships its configs inside its package.",
            )
            continue

        files = sorted(configs.rglob("*.yaml"))
        bad = []
        for path in files:
            try:
                yaml.safe_load(path.read_text(encoding="utf-8"))
            except yaml.YAMLError as exc:
                bad.append(f"{path.relative_to(configs)}: {exc}")

        if bad:
            report.add(
                f"configs [{spec.name}]",
                Level.FAIL,
                "\n".join(bad),
                "Fix the YAML. A config that does not parse fails at launch, after a\n"
                "scheduler queue wait, rather than here.",
            )
        elif not files:
            report.add(
                f"configs [{spec.name}]", Level.WARN, f"no .yaml files under {configs}"
            )
        else:
            report.add(
                f"configs [{spec.name}]", Level.OK, f"{len(files)} config file(s) parse"
            )
