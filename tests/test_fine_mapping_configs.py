"""Fine-mapping's configuration, scenarios, and the vendoring contract.

Three groups, each guarding something the port could quietly break.

1. **The vendoring contract.** ``fedfm/`` claims to be byte-identical to the standalone
   repository. That claim is only useful if it is checkable, so it is checked.

2. **The duplicated site list.** ``schema.py`` repeats ``sampling.py``'s ``SITE_ORDER``
   rather than importing it. Two copies of a list is exactly the arrangement that drifts,
   and drift here is silent: a scenario naming a site the sampler skips produces a
   smaller cohort than asked for, with no error.

3. **Scenario validation.** Every refusal in ``FineMappingScenario.validate`` exists
   because the alternative is a failure much later and much less legible. Each one is
   tested, because a validator nobody exercises is a validator that stops working.
"""

from __future__ import annotations

import filecmp
import os
from pathlib import Path

import pytest
import yaml

from appfl_bio_suite.core.experiments import REGISTRY, get_spec
from appfl_bio_suite.experiments.fine_mapping.simulation import (
    SCENARIO_DIR,
    list_scenarios,
    load_scenario,
)
from appfl_bio_suite.experiments.fine_mapping.simulation.schema import (
    KNOWN_SITE_IDS,
    FineMappingScenario,
)

SPEC = get_spec("fine-mapping")

# The eight modules copied verbatim. Listed here rather than globbed so that a new file
# appearing in fedfm/ -- which would not be vendored, and would break the contract --
# fails this test rather than being silently accepted.
VENDORED_MODULES = (
    "utils.py",
    "sampling.py",
    "locus_selection.py",
    "phenotype_sim.py",
    "qc.py",
    "validation.py",
    "fine_mapping.py",
    "fed_fine_mapping.py",
)

# Where the standalone repository lives, if it is available. Overridable, because it is
# outside this repo and its location is a fact about one machine.
_UPSTREAM = Path(
    os.environ.get(
        "FEDFM_UPSTREAM",
        Path(__file__).resolve().parents[3] / "fedfm-simulation",
    )
)


# ---------------------------------------------------------------------------
# the vendoring contract
# ---------------------------------------------------------------------------


def test_the_vendored_package_contains_exactly_the_declared_modules():
    """A ninth file in fedfm/ is either an unvendored edit or an undeclared addition."""
    package = SPEC.package_path / "fedfm"
    present = {p.name for p in package.glob("*.py")} - {"__init__.py"}
    assert present == set(VENDORED_MODULES), (
        f"fedfm/ holds {sorted(present)}, expected {sorted(VENDORED_MODULES)}.\n"
        "Everything in that package is a verbatim copy from the standalone repository. "
        "New code belongs in the experiment package above it, not here -- see "
        "fedfm/__init__.py."
    )


@pytest.mark.skipif(
    not (_UPSTREAM / "src").is_dir(),
    reason=f"the standalone repository is not at {_UPSTREAM}; set FEDFM_UPSTREAM",
)
@pytest.mark.parametrize("module", VENDORED_MODULES)
def test_vendored_module_is_byte_identical_to_upstream(module: str):
    """The whole value of vendoring verbatim is that this can be checked with cmp.

    A failure here does not necessarily mean something is wrong -- upstream may simply
    have moved on. It means the copy is stale and the claim in fedfm/__init__.py no
    longer holds, so either re-copy or amend the claim.
    """
    ours = SPEC.package_path / "fedfm" / module
    theirs = _UPSTREAM / "src" / module
    assert filecmp.cmp(ours, theirs, shallow=False), (
        f"{module} differs from {theirs}.\n"
        "fedfm/ is a verbatim vendoring; a local edit there breaks the guarantee that "
        "the code producing the numbers is the code the proofs are about. Make the "
        "change upstream and re-copy, or move it into the experiment package."
    )


def test_the_shipped_trainer_does_not_import_the_vendored_package():
    """The site stage exists twice on purpose; this is what keeps the second copy honest.

    tests/test_shipped_modules.py already forbids importing ``appfl_bio_suite`` at all.
    This states the specific consequence for this experiment, so that a failure explains
    why the duplication is there instead of looking like an arbitrary rule.
    """
    import ast

    tree = ast.parse((SPEC.package_path / "trainer.py").read_text(encoding="utf-8"))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
    # Substring matching would trip on the module docstring, which explains at length
    # why this import is forbidden -- prose about a rule is not a violation of it.
    assert not any("fedfm" in name for name in imported), (
        "trainer.py references the vendored package. Its source is shipped to partner "
        "workers, which do not install this suite, so it must stand alone -- that is why "
        "the site stage is implemented twice and why "
        "tests/test_fine_mapping_shipped_parity.py compares the two at zero tolerance."
    )


# ---------------------------------------------------------------------------
# the duplicated site list
# ---------------------------------------------------------------------------


def test_known_site_ids_matches_the_samplers_site_order():
    from appfl_bio_suite.experiments.fine_mapping.fedfm.sampling import SITE_ORDER

    assert KNOWN_SITE_IDS == SITE_ORDER, (
        f"schema.KNOWN_SITE_IDS is {KNOWN_SITE_IDS} but sampling.SITE_ORDER is "
        f"{SITE_ORDER}.\n"
        "They are duplicated so that schema.py stays importable without the scientific "
        "stack. Drift between them is silent in the worst way: the sampler skips any "
        "site it does not enumerate, so a scenario would produce a cohort smaller than "
        "it asked for with nothing raised."
    )


# ---------------------------------------------------------------------------
# the shipped scenarios
# ---------------------------------------------------------------------------


def test_scenarios_are_shipped():
    available = list_scenarios()
    assert available, f"no scenarios under {SCENARIO_DIR}"
    assert "ci-tiny" in available, "the smoke scenario is what CI and a new user run first"


@pytest.mark.parametrize("name", sorted(list_scenarios()))
def test_every_scenario_loads_and_validates(name: str):
    scenario = load_scenario(name)  # validate() runs inside
    assert scenario.site_ids
    assert scenario.total_samples > 0


@pytest.mark.parametrize("name", sorted(list_scenarios()))
def test_every_scenario_builds_a_valid_pipeline_config(name: str, tmp_path, monkeypatch):
    """The `pipeline:` block is handed to the vendored pydantic model unmodified.

    A scenario that parses as YAML but not as a SimulationConfig fails at the very start
    of a run that someone may have queued for hours.
    """
    scenario = load_scenario(name)
    if scenario.cohort.provider == "hapnest":
        # A staged pool is not present in CI; point at a directory so the config can be
        # built. Whether the pool is really there is materialize_pool's business.
        monkeypatch.setenv("FEDFM_HAPNEST_DIR", str(tmp_path / "staged"))
    cfg = scenario.to_pipeline_config(tmp_path)
    assert cfg.chromosome == scenario.chromosome
    assert set(cfg.sites) == set(scenario.site_ids)
    # Paths must be absolute and under the run directory, or two runs collide.
    for key in ("processed_dir", "loci_dir", "ground_truth_dir"):
        assert Path(cfg.resolved_path(key)).is_absolute()
        assert str(cfg.resolved_path(key)).startswith(str(tmp_path.resolve()))


def test_ci_tiny_needs_no_download():
    """The smoke scenario has to run on a fresh checkout, or nobody runs it."""
    scenario = load_scenario("ci-tiny")
    assert scenario.cohort.provider == "synthetic_hapnest"
    assert scenario.cohort.hapnest_dir is None


def test_the_published_scenario_uses_real_data_and_says_so():
    scenario = load_scenario("three-site-hapnest")
    assert scenario.cohort.provider == "hapnest"
    assert scenario.total_samples == 150_000
    assert set(scenario.site_ids) == set(KNOWN_SITE_IDS)


def test_the_published_scenario_refuses_to_run_without_a_staged_pool(tmp_path, monkeypatch):
    """`hapnest_dir: null` is the shipped state. It must fail loudly rather than
    silently substituting the synthetic pool -- a published figure produced from a
    substitute is the exact confusion the two providers exist to prevent.

    It must ALSO still load and list, which is why the check lives here and not in
    validate(): a coordinator has to be able to read the scenario before deciding to
    spend 135 GB on it.
    """
    from appfl_bio_suite.experiments.fine_mapping.simulation.cohort import materialize_pool

    monkeypatch.delenv("FEDFM_HAPNEST_DIR", raising=False)
    scenario = load_scenario("three-site-hapnest")   # loads fine

    with pytest.raises(ValueError, match="FEDFM_HAPNEST_DIR"):
        scenario.pool_dir(tmp_path)
    with pytest.raises(FileNotFoundError, match="FEDFM_HAPNEST_DIR"):
        materialize_pool(scenario.cohort, scenario.chromosome, tmp_path)


def test_the_staged_pool_can_come_from_the_environment(tmp_path, monkeypatch):
    """So that a committed scenario never carries one machine's path."""
    monkeypatch.setenv("FEDFM_HAPNEST_DIR", str(tmp_path))
    scenario = load_scenario("three-site-hapnest")
    assert scenario.pool_dir(tmp_path / "unused") == tmp_path.resolve()


# ---------------------------------------------------------------------------
# validation refuses what it claims to refuse
# ---------------------------------------------------------------------------


def _ci_tiny_dict() -> dict:
    return yaml.safe_load((SCENARIO_DIR / "ci-tiny.yaml").read_text(encoding="utf-8"))


def test_an_unknown_site_is_refused_with_the_reason():
    data = _ci_tiny_dict()
    data["pipeline"]["sites"]["harvard"] = data["pipeline"]["sites"].pop("anl")
    with pytest.raises(ValueError, match="silently drop"):
        FineMappingScenario.from_dict(data).validate()


def test_a_pool_too_small_for_the_sites_is_refused():
    data = _ci_tiny_dict()
    data["cohort"]["n_individuals"] = 100
    with pytest.raises(ValueError, match="disjoint"):
        FineMappingScenario.from_dict(data).validate()


def test_a_superpopulation_the_pool_cannot_supply_is_refused():
    """Caught here rather than inside the sampler, which raises without naming the fix."""
    data = _ci_tiny_dict()
    data["cohort"]["superpopulation_weights"]["MID"] = 0.001
    with pytest.raises(ValueError, match="MID"):
        FineMappingScenario.from_dict(data).validate()


def test_an_unknown_provider_is_refused():
    data = _ci_tiny_dict()
    data["cohort"]["provider"] = "1000genomes"
    with pytest.raises(ValueError, match="unknown cohort provider"):
        FineMappingScenario.from_dict(data).validate()


def test_a_scenario_missing_a_pipeline_block_is_refused():
    data = _ci_tiny_dict()
    del data["pipeline"]["architecture"]
    with pytest.raises(ValueError, match="pipeline.architecture"):
        FineMappingScenario.from_dict(data).validate()


# ---------------------------------------------------------------------------
# the APPFL configs
# ---------------------------------------------------------------------------


def test_server_config_is_single_round():
    """Aggregate FL has exactly one exchange. More rounds recompute identical second
    moments at full cost, which looks like a hang rather than an error."""
    config = yaml.safe_load(
        (SPEC.configs_path / "server.yaml").read_text(encoding="utf-8")
    )
    assert config["server_configs"]["num_global_epochs"] == 1


def test_the_answer_key_is_not_shipped_in_the_server_config():
    """`causal_manifest` must default to null. A committed path to someone's ground truth
    would make a fresh federation either crash or, worse, score against the wrong file."""
    for name in ("server.yaml", "server.loopback.yaml"):
        config = yaml.safe_load((SPEC.configs_path / name).read_text(encoding="utf-8"))
        kwargs = config["server_configs"]["aggregator_kwargs"]
        assert kwargs["causal_manifest"] is None, f"{name} hardcodes a causal_manifest"


def test_locus_sharding_defaults_are_coherent():
    """An out-of-range shard index selects no loci and reports success on an empty table."""
    for name in ("server.yaml", "server.loopback.yaml", "clients.template.yaml"):
        config = yaml.safe_load((SPEC.configs_path / name).read_text(encoding="utf-8"))
        block = (
            config.get("client_configs", {}).get("train_configs")
            or config["client_defaults"]["train_configs"]
        )
        assert block["locus_n_shards"] >= 1
        assert 0 <= block["locus_shard_index"] < block["locus_n_shards"]


def test_the_registry_declares_the_binaries_the_experiment_shells_out_to():
    """pip cannot install SuSiEx. Preflight only knows to look because the spec says so."""
    assert set(SPEC.required_binaries) == {"SuSiEx", "plink"}
    # And no other experiment claims to need one, which is what keeps the check quiet for
    # a coordinator running only GWAS or FLamby.
    for name, spec in REGISTRY.items():
        if name != "fine-mapping":
            assert spec.required_binaries == ()
