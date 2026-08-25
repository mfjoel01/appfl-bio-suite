"""Enforce the import rule for modules whose source is shipped to partner workers.

THE RULE
--------
APPFL reads a shipped module on the coordinating driver, inlines its source text, and
sends it to a partner's worker for execution. That worker has whatever the partner
installed and nothing else.

So a shipped module may import only:

* the standard library
* the scientific packages a partner installs with the experiment's extra
* ``appfl`` itself, which a partner needs anyway

and specifically may NOT import from ``appfl_bio_suite``, because a partner does not
install the suite. Nor may it import a sibling module by bare name -- the shipped source
runs from a temporary working directory where its siblings do not exist.

WHY THIS IS A TEST AND NOT A CONVENTION
---------------------------------------
The failure mode is remote, delayed, and misleading. A bad import surfaces as a
``ModuleNotFoundError`` on someone else's cluster, minutes into a run, after a scheduler
queue wait, and it reads like an infrastructure problem rather than an import problem.

This is not hypothetical. The GWAS trainer originally did::

    from gwas_config import ...
    from gwas_plot_utils import ...

Those two lines are the entire reason every GWAS partner had to add a ``PYTHONPATH``
entry to their endpoint's ``worker_init``, and a missing or wrong one was recorded as the
single largest source of partner-side breakage on that project. Removing the sibling
imports removed that line from the partner setup guide.

Catching it here costs milliseconds.
"""

from __future__ import annotations

import ast
import sys

import pytest

from appfl_bio_suite.core.experiments import REGISTRY, ExperimentSpec

# Third-party packages a partner installs as part of an experiment's extra. Keep this
# tight: every addition is something a partner must have for the experiment to run, and
# widening it to fix a failing test is almost always the wrong fix.
ALLOWED_THIRD_PARTY = {
    "numpy",
    "pandas",
    "scipy",
    "sklearn",
    "torch",
    "torchvision",
    "matplotlib",
    "pandas_plink",
    "flamby",
    "appfl",
}

# Packages a shipped module may reference ONLY from inside a function body, never at
# module scope. A module-level import of one of these would fail on any worker that does
# not have it; a function-level import guarded by configuration cannot, because the
# function is not called unless the feature is switched on.
#
# cuml is the case this exists for: GPU regression is opt-in, defaults to off, and is
# imported inside `if use_cuml:`. A partner who never enables it never needs the package.
# Keeping the capability costs nothing as long as the import stays deferred -- which is
# exactly what test_optional_imports_are_deferred checks.
OPTIONAL_THIRD_PARTY = {"cuml"}

FORBIDDEN_ROOTS = {"appfl_bio_suite"}


def _shipped_modules() -> list[tuple[str, str, object]]:
    out = []
    for spec in REGISTRY.values():
        if not spec.implemented:
            continue
        for module in spec.shipped_modules:
            out.append((spec.name, module, spec))
    return out


def _module_scope_imports(tree: ast.AST) -> set[str]:
    """Root package names imported at MODULE scope only.

    These run at import time on the worker, so every one of them must be present or the
    module fails to load before any code executes.
    """
    roots: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.Import):
            for alias in node.names:
                roots.add(alias.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom) and node.module:
            roots.add(node.module.split(".")[0])
    return roots


def _top_level_imports(tree: ast.AST) -> set[str]:
    """Every root package name imported anywhere in the module, at any nesting depth.

    Function-level imports count. Deferring a heavy import inside a function is good
    practice and this suite does it deliberately, but it does not change whether the
    package must be installed on the worker -- it only changes when the failure happens.
    """
    roots: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                roots.add(alias.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            if node.level and node.level > 0:
                # A relative import cannot resolve in shipped source, which executes
                # standalone with no package context.
                roots.add(f"<relative import, level {node.level}>")
            elif node.module:
                roots.add(node.module.split(".")[0])
    return roots


@pytest.mark.parametrize(
    ("experiment", "module", "spec"),
    _shipped_modules(),
    ids=[f"{e}:{m}" for e, m, _ in _shipped_modules()],
)
def test_shipped_module_imports_are_worker_safe(
    experiment: str, module: str, spec: ExperimentSpec
):
    path = spec.package_path / f"{module}.py"
    assert path.is_file(), (
        f"{spec.name} declares shipped module '{module}' but {path} does not exist"
    )

    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imports = _top_level_imports(tree)

    stdlib = sys.stdlib_module_names
    offenders = sorted(
        name
        for name in imports
        if name not in stdlib
        and name not in ALLOWED_THIRD_PARTY
        and name not in OPTIONAL_THIRD_PARTY
    )

    assert not offenders, (
        f"{path.name} imports {offenders}, which a partner's worker may not have.\n"
        "\n"
        "This module's SOURCE is shipped over the wire and executed on a partner's\n"
        "cluster. It may import the standard library, appfl, and the scientific\n"
        "packages the experiment's extra installs -- nothing else.\n"
        "\n"
        "If the import is genuinely needed, move the code that needs it into a\n"
        "coordinator-side module, or inline it. Widening ALLOWED_THIRD_PARTY means\n"
        "committing every partner to installing that package."
    )


@pytest.mark.parametrize(
    ("experiment", "module", "spec"),
    _shipped_modules(),
    ids=[f"{e}:{m}" for e, m, _ in _shipped_modules()],
)
def test_shipped_module_does_not_import_the_suite(
    experiment: str, module: str, spec: ExperimentSpec
):
    """The specific rule, stated separately so a failure names the actual problem."""
    path = spec.package_path / f"{module}.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imports = _top_level_imports(tree)

    bad = sorted(imports & FORBIDDEN_ROOTS)
    assert not bad, (
        f"{path.name} imports {bad}. A partner installs the experiment's extra, not\n"
        "appfl_bio_suite, so this import fails on their worker.\n"
        "\n"
        "Shipped modules are self-contained by design. That is what lets a partner run\n"
        "the experiment without installing this package, and it is what removed the\n"
        "PYTHONPATH line from the partner setup guide."
    )


@pytest.mark.parametrize(
    ("experiment", "module", "spec"),
    _shipped_modules(),
    ids=[f"{e}:{m}" for e, m, _ in _shipped_modules()],
)
def test_shipped_module_has_no_relative_imports(
    experiment: str, module: str, spec: ExperimentSpec
):
    """Relative imports cannot resolve in shipped source; it runs with no package."""
    path = spec.package_path / f"{module}.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    relative = [n for n in _top_level_imports(tree) if n.startswith("<relative")]
    assert not relative, (
        f"{path.name} uses a relative import. Shipped source is executed standalone in a\n"
        "temporary working directory with no package context, so `from . import x` has\n"
        "nothing to resolve against."
    )


@pytest.mark.parametrize(
    ("experiment", "module", "spec"),
    _shipped_modules(),
    ids=[f"{e}:{m}" for e, m, _ in _shipped_modules()],
)
def test_shipped_module_parses_standalone(experiment: str, module: str, spec: ExperimentSpec):
    """The source must compile on its own -- that is exactly how the worker gets it."""
    path = spec.package_path / f"{module}.py"
    source = path.read_text(encoding="utf-8")
    compile(source, str(path), "exec")


def test_every_declared_shipped_module_exists():
    """A registry entry naming a file that does not exist fails at launch, not here."""
    missing = []
    for spec in REGISTRY.values():
        if not spec.implemented:
            continue
        for path in spec.shipped_module_paths():
            if not path.is_file():
                missing.append(str(path))
    assert not missing, f"declared shipped modules do not exist: {missing}"


@pytest.mark.parametrize(
    ("experiment", "module", "spec"),
    _shipped_modules(),
    ids=[f"{e}:{m}" for e, m, _ in _shipped_modules()],
)
def test_optional_imports_are_deferred(experiment: str, module: str, spec: ExperimentSpec):
    """An optional package must never be imported at module scope.

    Optional means "a partner may not have this". Importing it at module scope makes it
    mandatory in practice, because the module fails to load on any worker without it --
    before a single line of the experiment runs, and regardless of whether the feature it
    supports was ever switched on.

    Deferring it into the function that needs it keeps the capability available to
    whoever has the package, and free for everyone else.
    """
    path = spec.package_path / f"{module}.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    eager = _module_scope_imports(tree) & OPTIONAL_THIRD_PARTY
    assert not eager, (
        f"{path.name} imports {sorted(eager)} at module scope. These are optional "
        "packages -- move the import inside the function that uses it, so a worker "
        "without the package can still load this module."
    )


def test_partner_extras_do_not_include_coordinator_only_packages():
    """A partner's install must not pull in what only the coordinator needs.

    This is a packaging fact rather than a documentation promise: keeping simulation and
    plotting dependencies out of the extra a partner installs is what makes
    "simulation is coordinator-side only" true rather than merely stated.
    """
    import re

    from appfl_bio_suite.core.experiments import repo_root

    pyproject = (repo_root() / "pyproject.toml").read_text(encoding="utf-8")

    coordinator_only = {"matplotlib", "pandas-plink", "xarray", "dask", "joblib"}
    # pandas-plink IS needed by the GWAS trainer to read genotypes, so it is legitimately
    # in the partner extra. The rest are not.
    #
    # joblib is the fine-mapping case: it parallelizes SuSiEx across instances, which
    # happens entirely at the coordinator. The fine-mapping partner extra is deliberately
    # two packages, and this is what keeps it that way -- note that the FINE-MAPPING
    # trainer does not import pandas-plink either, because it reads the .bed itself.
    coordinator_only.discard("pandas-plink")

    for extra in ("flamby", "gwas", "finemapping"):
        match = re.search(rf"^{extra} = \[(.*?)\]", pyproject, re.S | re.M)
        assert match, f"could not find the [{extra}] extra in pyproject.toml"
        body = match.group(1)
        for package in sorted(coordinator_only):
            assert package not in body, (
                f"[{extra}] is installed by partners and must not include {package!r}.\n"
                "Coordinator-only dependencies belong in the gwas-sim extra."
            )
