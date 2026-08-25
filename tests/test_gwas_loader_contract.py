"""The simulated output layout and the client loader must agree.

WHY THIS IS A TEST
------------------
The bundler writes per-site data; the loader reads it. Those two lived in different
directories, maintained by hand, and nothing checked that they agreed. A change to either
side surfaced as a FileNotFoundError on a partner's cluster, minutes into a run, after a
scheduler queue wait -- which is both the slowest and the most embarrassing place to
discover a filename typo.

Now a mismatch fails in CI instead.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from appfl_bio_suite.experiments.gwas import dataset as loader
from appfl_bio_suite.experiments.gwas.simulation import bundler

sys.path.insert(0, str(Path(__file__).parent / "fixtures"))


def test_loader_requirements_are_a_subset_of_what_the_bundler_writes():
    """Every file the loader demands must actually be produced."""
    written = set(bundler.SITE_DATA_FILES)
    required = set(loader.REQUIRED_FILES)
    missing = required - written
    assert not missing, (
        f"the loader requires {sorted(missing)}, which the bundler does not write.\n"
        "Either the bundler stopped producing them or the loader started asking for "
        "files that were never part of a site bundle. Sites already given data would "
        "fail on their next run."
    )


def test_bundler_and_loader_agree_on_the_required_set():
    """The two modules' own declarations of 'required' must not drift apart."""
    assert set(bundler.REQUIRED_BY_LOADER) == set(loader.REQUIRED_FILES), (
        "bundler.REQUIRED_BY_LOADER and dataset.REQUIRED_FILES disagree.\n"
        f"  bundler: {sorted(bundler.REQUIRED_BY_LOADER)}\n"
        f"  loader:  {sorted(loader.REQUIRED_FILES)}\n"
        "These describe the same contract from opposite sides and must match."
    )


def test_plink_stem_matches():
    """A stem mismatch renames every genotype file and breaks every existing bundle."""
    assert bundler.PLINK_STEM == loader.PLINK_STEM


@pytest.fixture(scope="module")
def simulated(tmp_path_factory) -> Path:
    """Actually run the pipeline, so the contract is checked against real output."""
    from appfl_bio_suite.experiments.gwas.simulation import load_scenario, run_simulation

    out = tmp_path_factory.mktemp("gwas-contract")
    run_simulation(load_scenario("ci-tiny"), out)
    return out


def test_simulated_output_loads(simulated: Path):
    """The end-to-end check: what simulation writes, the loader must open."""
    for site_id in ("Site1", "Site2"):
        data_dir = simulated / site_id / "data"
        assert data_dir.is_dir(), f"simulation produced no data directory for {site_id}"

        dataset, val = loader.get_dataset(data_dir=str(data_dir), site_id=site_id)
        assert val is None, "GWAS is single-round; there is no validation split"
        assert len(dataset) > 0, f"{site_id} loaded zero samples"
        assert dataset.site_id == site_id


def test_simulated_sites_have_the_expected_sample_counts(simulated: Path):
    """Site sizes must match the scenario, or the split silently lost people."""
    from appfl_bio_suite.experiments.gwas.simulation import load_scenario

    scenario = load_scenario("ci-tiny")
    for index, expected in enumerate(scenario.sites, start=1):
        site_id = f"Site{index}"
        dataset, _ = loader.get_dataset(data_dir=str(simulated / site_id / "data"), site_id=site_id)
        assert len(dataset) == expected, (
            f"{site_id} has {len(dataset)} samples, scenario declares {expected}"
        )


def test_loader_error_names_the_missing_files(tmp_path):
    """A partner reading this message must be able to act on it without asking."""
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError) as excinfo:
        loader.get_dataset(data_dir=str(empty), site_id="Site1")

    message = str(excinfo.value)
    assert "missing" in message
    for name in loader.REQUIRED_FILES:
        assert name in message, f"the error does not name the missing file {name}"


def test_loader_rejects_a_nonexistent_directory(tmp_path):
    """The most common partner mistake: a path that is right on the wrong machine."""
    with pytest.raises(FileNotFoundError, match="does not exist"):
        loader.get_dataset(data_dir=str(tmp_path / "nope"), site_id="Site1")
