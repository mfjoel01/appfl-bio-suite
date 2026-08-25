"""Preflight must not report success for work it did not do.

A false green here is worse than no check at all: the whole point of the command is that
someone trusts it instead of going and looking.
"""

from __future__ import annotations

import pytest

from appfl_bio_suite.core.config import Federation, load_federation
from appfl_bio_suite.core.experiments import repo_root
from appfl_bio_suite.core.preflight import Level, run_preflight

EXAMPLE = repo_root() / "federation.yaml.example"


@pytest.fixture(scope="module")
def federation():
    return load_federation(EXAMPLE)


def test_every_declared_check_group_actually_runs_something(federation):
    """`--check data` was a declared group that ran nothing and printed 'All checks passed'."""
    for group in ("env", "pins", "configs", "data"):
        report = run_preflight(federation=federation, check=group)
        assert report.checks, f"--check {group} ran no checks but still reported a verdict"


def test_data_checks_do_not_claim_to_have_verified_a_partners_cluster(federation):
    """A path on someone else's machine is unknowable from here -- and must say so."""
    report = run_preflight(federation=federation, check="data")
    gwas = [c for c in report.checks if c.name.startswith("data gwas/")]
    assert len(gwas) == len(federation.experiment("gwas").sites)
    assert all(c.level is Level.SKIP for c in gwas)
    assert all("not visible from here" in c.detail for c in gwas)


def test_data_checks_validate_a_site_whose_data_is_visible(tmp_path, federation):
    from appfl_bio_suite.experiments.gwas.dataset import PLINK_STEM, REQUIRED_FILES

    site = tmp_path / "Site1" / "data"
    site.mkdir(parents=True)
    for name in REQUIRED_FILES:
        (site / name).write_text("")
    (site / f"{PLINK_STEM}.fam").write_text("F1 I1 0 0 1 -9\nF2 I2 0 0 2 -9\n")

    entry = federation.experiment("gwas").sites[0]
    original_dir, original_n = entry.data_dir, entry.expected_samples
    object.__setattr__(entry, "data_dir", str(site))
    object.__setattr__(entry, "expected_samples", 2)
    try:
        ok = run_preflight(federation=federation, check="data")
        object.__setattr__(entry, "expected_samples", 99)
        wrong = run_preflight(federation=federation, check="data")
    finally:
        object.__setattr__(entry, "data_dir", original_dir)
        object.__setattr__(entry, "expected_samples", original_n)

    assert [c for c in ok.checks if c.name == "data gwas/Site1"][0].level is Level.OK
    mismatch = [c for c in wrong.checks if c.name == "data gwas/Site1"][0]
    assert mismatch.level is Level.FAIL, "a site holding the wrong dataset must not pass"


def test_duplicate_site_ids_are_named(tmp_path):
    """Detection worked; the message reported an empty list and named nobody."""
    with pytest.raises(Exception, match=r"duplicate site id\(s\): \['dup'\]"):
        Federation.model_validate(
            {
                "schema_version": 1,
                "coordinator": {"identity": "a.researcher@example-university.edu"},
                "sites": [
                    {"id": "dup", "name": "One"},
                    {"id": "dup", "name": "Two"},
                    {"id": "other", "name": "Three"},
                ],
            }
        )
