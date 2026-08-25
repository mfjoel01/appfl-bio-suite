"""Prove the repository works for a coordinator who is not its author.

This is the primary acceptance criterion, and the one most likely to be quietly missed --
because everything works on the original author's machine whether or not it holds.

The check: load a second, entirely fictional federation with a different coordinator,
different institutions, a different site count, a different scheduler mix and different
service accounts, then generate configs and a full partner bundle from it. No code change,
no edit to any shipped file.

If any of these fail, something in the suite is carrying an assumption about one
particular federation.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from appfl_bio_suite.core.config import load_federation
from appfl_bio_suite.core.experiments import repo_root

SECOND = Path(__file__).parent / "fixtures" / "second_federation.yaml"
EXAMPLE = repo_root() / "federation.yaml.example"


@pytest.fixture(scope="module")
def federation():
    return load_federation(SECOND)


def test_it_loads(federation):
    assert federation.coordinator.identity == "dr.k.tanaka@marine-genomics.example.org"
    assert federation.coordinator.organization == "Marine Genomics Consortium"
    assert len(federation.sites) == 2


def test_it_shares_nothing_with_the_shipped_example(federation):
    """A test that accidentally matched the example would prove nothing."""
    example = load_federation(EXAMPLE)
    assert federation.coordinator.identity != example.coordinator.identity
    assert {s.id for s in federation.sites}.isdisjoint({s.id for s in example.sites})
    assert (
        federation.experiment("gwas").service_account
        != example.experiment("gwas").service_account
    )


def test_a_coordinator_who_only_dispatches_is_supported(federation):
    """Not every coordinator is also a training site."""
    assert federation.coordinator.endpoint is None


def test_configs_resolve(federation):
    from appfl_bio_suite.core.launch import resolve_run

    run = resolve_run(federation, "gwas")
    assert run.num_clients == 2
    assert run.server_config["server_configs"]["num_clients"] == 2

    client_ids = {c["client_id"] for c in run.client_configs}
    assert client_ids == {"Harbour", "Reef"}

    for client in run.client_configs:
        kwargs = client["data_configs"]["dataset_kwargs"]
        assert kwargs["data_dir"].startswith("/")
        assert kwargs["site_id"] in client_ids


@pytest.mark.parametrize("site", ["harbour-lab", "reef-station"])
def test_partner_bundle_generates(federation, tmp_path, site):
    from appfl_bio_suite.core.partner import generate_bundle

    destination = generate_bundle(federation, "gwas", site, tmp_path)

    for expected in (
        "1-endpoint-setup.md",
        "2-gwas.md",
        "your-config-block.yaml",
        "example_identity_mapping_config.json",
        "user_config_template.yaml.j2",
        "bundle_manifest.json",
    ):
        assert (destination / expected).is_file(), f"bundle is missing {expected}"


def test_the_bundle_carries_the_second_coordinators_identity(federation, tmp_path):
    """The whole point: the identity comes from the config, never from a literal."""
    from appfl_bio_suite.core.partner import generate_bundle

    destination = generate_bundle(federation, "gwas", "harbour-lab", tmp_path)

    mapping = json.loads(
        (destination / "example_identity_mapping_config.json").read_text(encoding="utf-8")
    )
    rule = mapping[0]["mappings"][0]

    assert rule["match"] == r"dr\.k\.tanaka@marine-genomics\.example\.org"
    assert rule["output"] == "marine_gwas"

    # The anchors that break everything must not be there.
    assert not rule["match"].startswith("^")
    assert not rule["match"].endswith("$")

    # And nothing from the shipped example may have leaked in.
    example = load_federation(EXAMPLE)
    assert example.coordinator.identity not in json.dumps(mapping)


def test_the_generated_mapping_actually_resolves(federation, tmp_path):
    """Validate the generated file with the real Globus library, end to end."""
    pytest.importorskip("globus_identity_mapping")
    from appfl_bio_suite.core.identity import validate_mapping_file
    from appfl_bio_suite.core.partner import generate_bundle

    destination = generate_bundle(federation, "gwas", "harbour-lab", tmp_path)
    report = validate_mapping_file(
        destination / "example_identity_mapping_config.json",
        identity=federation.coordinator.identity,
        expected_output="marine_gwas",
        identity_id=federation.coordinator.identity_id,
    )
    assert report.passed, report.render()


def test_the_endpoint_template_matches_each_sites_scheduler(federation, tmp_path):
    """A PBS site and a SLURM site must get different provider blocks."""
    from appfl_bio_suite.core.partner import generate_bundle

    pbs = (generate_bundle(federation, "gwas", "harbour-lab", tmp_path / "a")
           / "user_config_template.yaml.j2").read_text(encoding="utf-8")
    assert "PBSProProvider" in pbs
    assert "queue: research" in pbs
    assert "MpiExecLauncher" in pbs
    assert "SlurmProvider" not in pbs

    slurm = (generate_bundle(federation, "gwas", "reef-station", tmp_path / "b")
             / "user_config_template.yaml.j2").read_text(encoding="utf-8")
    assert "SlurmProvider" in slurm
    assert "partition: batch" in slurm
    assert "SrunLauncher" in slurm
    assert "PBSProProvider" not in slurm
    # This site declares an interface; the template must include it.
    assert "ifname: ib0" in slurm


def test_the_config_block_carries_this_sites_own_paths(federation, tmp_path):
    from appfl_bio_suite.core.partner import generate_bundle

    destination = generate_bundle(federation, "gwas", "reef-station", tmp_path)
    block = yaml.safe_load(
        (destination / "your-config-block.yaml").read_text(encoding="utf-8")
    )[0]

    assert block["client_id"] == "Reef"
    assert block["data_configs"]["dataset_kwargs"]["data_dir"] == "/scratch/marine_gwas/reef"
    # dataset_path is resolved on the coordinator's driver; a partner is never given one.
    assert "dataset_path" not in block["data_configs"]


def test_the_partner_docs_render_with_no_placeholders(federation, tmp_path):
    import re

    from appfl_bio_suite.core.partner import generate_bundle

    destination = generate_bundle(federation, "gwas", "harbour-lab", tmp_path)
    for doc in destination.glob("*.md"):
        text = doc.read_text(encoding="utf-8")
        assert not re.findall(r"\{\{[^}]+\}\}", text), f"{doc.name} has unrendered variables"
        assert "dr.k.tanaka@marine-genomics.example.org" in text or doc.name.startswith("2-")


def test_a_disabled_experiment_is_not_offered(federation):
    enabled = federation.enabled_experiments()
    assert "gwas" in enabled
    assert "flamby-heart-disease" not in enabled
    assert "fine-mapping" not in enabled
