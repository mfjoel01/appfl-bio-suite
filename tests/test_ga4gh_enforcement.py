"""What actually happens to a run when the GA4GH layer says no.

Three places can refuse, and they are not redundant:

* the **site's loader** (``dataset.py``), in the site's own process, before any genotype
  is opened. The only one that is a control rather than a courtesy.
* the **coordinator's launch gate**, which stops a run that every site would refuse
  anyway, saving three scheduler queue waits to be told so.
* **preflight**, which says the same thing without launching anything at all.

These tests pin the behaviour of all three, including the cases where refusing would be
wrong -- a site with no declared terms, and a run with enforcement deliberately waived.
"""

from __future__ import annotations

import json

import pytest

from appfl_bio_suite.core.config import Federation
from appfl_bio_suite.core.ga4gh.drs import build_registry
from appfl_bio_suite.core.ga4gh.duo import DataUseProfile, write_profile
from appfl_bio_suite.core.launch import build_client_configs, enforce_data_use
from appfl_bio_suite.experiments.fine_mapping.dataset import SiteFineMappingDataset

# A study that is genetic research, non-commercial, by a not-for-profit, publishing.
STUDY = {
    "requester": "me@example.org",
    "project": "fedfm",
    "institution": "Example University",
    "purposes": ["DUO:0000038"],
    "non_commercial": True,
    "not_for_profit_organisation": True,
    "publication_agreed": True,
}


def make_bundle(root, terms: dict | None = None):
    """A minimal but complete fine-mapping bundle, optionally with DUO terms."""
    data = root / "data"
    (data / "phenotypes").mkdir(parents=True)
    (data / "site_genotypes.bed").write_bytes(bytes([0x6C, 0x1B, 0x01]))
    (data / "site_genotypes.bim").write_text("1\trs1\t0\t100\tA\tG\n")
    (data / "site_genotypes.fam").write_text("1 IND1 0 0 0 -9\n")
    (data / "site_manifest.tsv").write_text("FID\tIID\tsuperpopulation\n1\tIND1\tEUR\n")
    (data / "reference_variants.tsv").write_text("snp_id\tchrom\tbp\ta1\ta2\nrs1\t1\t100\tA\tG\n")
    (data / "selected_loci.tsv").write_text("locus_id\tchrom\tstart_bp\tend_bp\nL0\t1\t1\t1000\n")
    (data / "phenotypes" / "L0_arch_rep0.pheno").write_text("FID\tIID\ty\n1\tIND1\t0.5\n")
    if terms is not None:
        write_profile(
            DataUseProfile(dataset_id="test/site", **terms), data / "DATA_USE.json"
        )
    return data


# ---------------------------------------------------------------------------
# the site
# ---------------------------------------------------------------------------


def test_a_permitted_study_loads(tmp_path):
    data = make_bundle(tmp_path, {"permission": "DUO:0000006"})
    dataset = SiteFineMappingDataset(data, "site", data_use_request=STUDY)
    assert dataset.data_use["outcome"] == "permitted"
    assert dataset.sample_size == 1


def test_a_refused_study_never_reads_the_data(tmp_path):
    """The refusal has to happen before the loader touches the cohort, not after."""
    data = make_bundle(tmp_path, {"permission": "DUO:0000011"})  # ancestry research only
    with pytest.raises(PermissionError) as exc:
        SiteFineMappingDataset(data, "site", data_use_request=STUDY)
    assert "do not permit this study" in str(exc.value)
    assert "No data was read" in str(exc.value)
    assert "DUO:0000011" in str(exc.value)


def test_a_bundle_with_terms_refuses_a_run_that_declares_nothing(tmp_path):
    """A consent code cannot be satisfied by a study that says nothing about itself."""
    data = make_bundle(tmp_path, {"permission": "DUO:0000042"})
    with pytest.raises(PermissionError, match="declared no\ndata use request|no data use request"):
        SiteFineMappingDataset(data, "site")


def test_a_bundle_without_terms_is_not_refused(tmp_path):
    """A site that declared nothing has refused nothing. Backwards compatible on purpose."""
    data = make_bundle(tmp_path)
    dataset = SiteFineMappingDataset(data, "site", data_use_request=STUDY)
    assert dataset.data_use["status"] == "no-profile"


def test_waiving_enforcement_still_records_the_refusal(tmp_path):
    """A result computed under a waived check must not look like a permitted one."""
    data = make_bundle(tmp_path, {"permission": "DUO:0000011"})
    dataset = SiteFineMappingDataset(
        data, "site", data_use_request=STUDY, enforce_data_use=False
    )
    assert dataset.data_use["outcome"] == "denied"
    assert dataset.data_use["enforced"] is False


def test_an_unreadable_profile_is_fatal_rather_than_ignored(tmp_path):
    data = make_bundle(tmp_path, {"permission": "DUO:0000042"})
    (data / "DATA_USE.json").write_text("{ not json", encoding="utf-8")
    with pytest.raises(ValueError, match="must not be computed over"):
        SiteFineMappingDataset(data, "site", data_use_request=STUDY)


# ---------------------------------------------------------------------------
# DRS verification at the site
# ---------------------------------------------------------------------------


def drs_object_for(root, site="s1"):
    registry = build_registry(root, hostname="drs.test.org")
    obj = registry.by_name(site)
    return {
        "id": obj.id,
        "name": obj.name,
        "self_uri": obj.self_uri,
        "size": obj.size,
        "checksums": [c.model_dump() for c in obj.checksums],
        "contents": [{"name": c.name, "id": c.id} for c in obj.contents or []],
    }


def test_a_matching_bundle_verifies(tmp_path):
    make_bundle(tmp_path / "s1")
    obj = drs_object_for(tmp_path)
    dataset = SiteFineMappingDataset(tmp_path / "s1" / "data", "s1", drs_object=obj)
    assert dataset.drs["problems"] == []
    assert dataset.drs["checked"] >= 3


def test_a_modified_bundle_is_refused(tmp_path):
    """The failure DRS was added for: aggregates over data the results do not describe."""
    make_bundle(tmp_path / "s1")
    obj = drs_object_for(tmp_path)
    (tmp_path / "s1" / "data" / "site_manifest.tsv").write_text(
        "FID\tIID\tsuperpopulation\n1\tIND1\tAFR\n"
    )
    with pytest.raises(ValueError, match="sha-256 mismatch"):
        SiteFineMappingDataset(tmp_path / "s1" / "data", "s1", drs_object=obj)


def test_an_extra_file_is_refused(tmp_path):
    make_bundle(tmp_path / "s1")
    obj = drs_object_for(tmp_path)
    (tmp_path / "s1" / "data" / "leftover.tsv").write_text("from another run\n")
    with pytest.raises(ValueError, match="not in the DRS object"):
        SiteFineMappingDataset(tmp_path / "s1" / "data", "s1", drs_object=obj)


def test_the_consent_file_is_exempt_from_the_object(tmp_path):
    """DATA_USE.json describes the object, so it cannot be inside its own checksum."""
    make_bundle(tmp_path / "s1")
    obj = drs_object_for(tmp_path)
    write_profile(
        DataUseProfile(dataset_id="test/site", permission="DUO:0000042"),
        tmp_path / "s1" / "data" / "DATA_USE.json",
    )
    dataset = SiteFineMappingDataset(
        tmp_path / "s1" / "data", "s1", drs_object=obj, data_use_request=STUDY
    )
    assert dataset.drs["problems"] == []


def test_metadata_mode_skips_the_genotype_matrix(tmp_path):
    """The only file large enough for hashing to cost anything, and the default skips it."""
    make_bundle(tmp_path / "s1")
    obj = drs_object_for(tmp_path)
    dataset = SiteFineMappingDataset(
        tmp_path / "s1" / "data", "s1", drs_object=obj, verify_bundles="metadata"
    )
    assert any("site_genotypes.bed" in s for s in dataset.drs["skipped"])


def test_full_mode_checks_it(tmp_path):
    make_bundle(tmp_path / "s1")
    obj = drs_object_for(tmp_path)
    dataset = SiteFineMappingDataset(
        tmp_path / "s1" / "data", "s1", drs_object=obj, verify_bundles="full"
    )
    assert not any("site_genotypes.bed" in s for s in dataset.drs["skipped"])


def test_off_checks_nothing(tmp_path):
    """Tampering with the genotype matrix -- which the loader never opens -- goes
    unnoticed with verification off, and is caught with it on `full`. That contrast is
    the whole argument for the setting having three values rather than a boolean."""
    make_bundle(tmp_path / "s1")
    obj = drs_object_for(tmp_path)
    (tmp_path / "s1" / "data" / "site_genotypes.bed").write_bytes(b"tampered")

    dataset = SiteFineMappingDataset(
        tmp_path / "s1" / "data", "s1", drs_object=obj, verify_bundles="off"
    )
    assert dataset.drs["checked"] == 0

    with pytest.raises(ValueError, match="sha-256 mismatch"):
        SiteFineMappingDataset(
            tmp_path / "s1" / "data", "s1", drs_object=obj, verify_bundles="full"
        )


# ---------------------------------------------------------------------------
# the coordinator
# ---------------------------------------------------------------------------


def federation_for(tmp_path, terms: dict, request: dict | None = None) -> Federation:
    data = make_bundle(tmp_path / "s1", terms)
    ga4gh: dict = {"enforce_data_use": True, "verify_bundles": "metadata"}
    if request is not None:
        ga4gh["data_use_request"] = request
    return Federation.model_validate(
        {
            "schema_version": 1,
            "coordinator": {"identity": "me@example.org", "organization": "Example University"},
            "sites": [{"id": "s1", "name": "Site One"}],
            "experiments": {
                "fine-mapping": {
                    "service_account": "svc",
                    "endpoint_name": "ep",
                    "ga4gh": ga4gh,
                    "sites": [
                        {
                            "site": "s1",
                            "client_id": "S1",
                            "endpoint_uuid": "dddddddd-0000-0000-0000-000000000001",
                            "output_dir": "/out",
                            "data_dir": str(data),
                            "data_use_profile": str(data / "DATA_USE.json"),
                        }
                    ],
                }
            },
        }
    )


def test_the_launch_gate_passes_a_permitted_study(tmp_path):
    fed = federation_for(tmp_path, {"permission": "DUO:0000006"}, STUDY)
    ok, report = enforce_data_use(fed, "fine-mapping")
    assert ok and "permitted" in report


def test_the_launch_gate_stops_a_refused_study(tmp_path):
    fed = federation_for(tmp_path, {"permission": "DUO:0000011"}, STUDY)
    ok, report = enforce_data_use(fed, "fine-mapping")
    assert not ok
    assert "Not launching" in report


def test_the_launch_gate_warns_but_proceeds_when_enforcement_is_off(tmp_path):
    fed = federation_for(tmp_path, {"permission": "DUO:0000011"}, STUDY)
    fed.experiment("fine-mapping").ga4gh.enforce_data_use = False
    ok, report = enforce_data_use(fed, "fine-mapping")
    assert ok
    assert "WARNING" in report and "own workers will still refuse" in report


def test_the_request_reaches_the_site_as_a_request_not_a_verdict(tmp_path):
    """The site decides. Sending a coordinator-computed verdict would make it a formality."""
    fed = federation_for(tmp_path, {"permission": "DUO:0000042"}, STUDY)
    kwargs = build_client_configs(fed, "fine-mapping")[0]["data_configs"]["dataset_kwargs"]
    assert kwargs["data_use_request"]["purposes"] == ["DUO:0000038"]
    assert "outcome" not in kwargs["data_use_request"]
    assert kwargs["enforce_data_use"] is True
    assert kwargs["verify_bundles"] == "metadata"


def test_the_requester_defaults_to_the_coordinator_identity(tmp_path):
    """The same identity partners already authorize, so DUO:0000026 checks something real."""
    request = dict(STUDY)
    request.pop("requester")
    fed = federation_for(tmp_path, {"permission": "DUO:0000042"}, request)
    assert fed.data_use_request("fine-mapping").requester == "me@example.org"


def test_a_federation_without_ga4gh_gets_untouched_client_configs(tmp_path):
    """Adding four specifications to a working experiment had to change nothing by default."""
    fed = federation_for(tmp_path, {"permission": "DUO:0000042"}, STUDY)
    fed.experiment("fine-mapping").ga4gh = None
    kwargs = build_client_configs(fed, "fine-mapping")[0]["data_configs"]["dataset_kwargs"]
    assert set(kwargs) == {"data_dir", "site_id"}


# ---------------------------------------------------------------------------
# preflight
# ---------------------------------------------------------------------------


def test_preflight_fails_on_a_refused_study(tmp_path):
    from appfl_bio_suite.core.preflight import Level, run_preflight

    fed = federation_for(tmp_path, {"permission": "DUO:0000011"}, STUDY)
    report = run_preflight(federation=fed, experiment="fine-mapping", check="ga4gh")
    failures = [c for c in report.failed if c.name.startswith("data use")]
    assert failures and failures[0].level is Level.FAIL
    assert "acknowledged" in failures[0].fix


def test_preflight_warns_when_the_study_is_undeclared(tmp_path):
    from appfl_bio_suite.core.preflight import run_preflight

    fed = federation_for(tmp_path, {"permission": "DUO:0000042"}, None)
    report = run_preflight(federation=fed, experiment="fine-mapping", check="ga4gh")
    warned = [c for c in report.warned if c.name.startswith("data use")]
    assert warned, "a study that declares nothing must not pass silently"
    assert "no `ga4gh.data_use_request` declared" in warned[0].detail
    assert "will refuse the task" in warned[0].fix


def test_preflight_skips_cleanly_when_ga4gh_is_not_configured(tmp_path):
    from appfl_bio_suite.core.preflight import Level, run_preflight

    fed = federation_for(tmp_path, {"permission": "DUO:0000042"}, STUDY)
    fed.experiment("fine-mapping").ga4gh = None
    report = run_preflight(federation=fed, experiment="fine-mapping", check="ga4gh")
    assert [c for c in report.checks if c.name == "ga4gh" and c.level is Level.SKIP]


def test_preflight_reports_the_ontology_release(tmp_path):
    from appfl_bio_suite.core.preflight import run_preflight

    fed = federation_for(tmp_path, {"permission": "DUO:0000042"}, STUDY)
    report = run_preflight(federation=fed, experiment="fine-mapping", check="ga4gh")
    ontology = [c for c in report.checks if c.name == "duo ontology"]
    assert ontology and "duo/releases/" in ontology[0].detail


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------


def test_the_trainer_carries_the_decision_into_the_payload(tmp_path):
    """The site attests to what it did; the coordinator can only record what it sent."""
    from appfl_bio_suite.experiments.fine_mapping.trainer import SiteFineMappingTrainer

    data = make_bundle(tmp_path / "s1", {"permission": "DUO:0000042"})
    dataset = SiteFineMappingDataset(data, "S1", data_use_request=STUDY)
    trainer = SiteFineMappingTrainer(
        train_dataset=dataset,
        train_configs={
            "trainer_output_dirname": str(tmp_path / "out"),
            "ga4gh_tool": {"id": "#workflow/x", "version": "0.1.0"},
        },
        logger=_QuietLogger(),
        client_id="S1",
    )
    provenance = trainer._provenance()
    assert provenance["data_use"]["outcome"] == "permitted"
    assert provenance["tool"]["id"] == "#workflow/x"
    assert provenance["drs"]["mode"] in ("off", "metadata")


class _QuietLogger:
    def info(self, *args, **kwargs):
        pass

    def warning(self, *args, **kwargs):
        pass


def test_run_provenance_names_every_specification_version(tmp_path):
    from appfl_bio_suite.core.ga4gh.resolve import run_provenance

    fed = federation_for(tmp_path, {"permission": "DUO:0000042"}, STUDY)
    record = run_provenance(fed, "fine-mapping")
    assert set(record["specifications"]) == {"drs", "tes", "trs", "duo"}
    assert record["data_use_request"]["purposes"] == ["DUO:0000038"]
    assert json.dumps(record)  # must be serializable; it is written to disk verbatim


def test_the_launch_gate_stops_a_run_that_declares_no_study(tmp_path):
    """Not the coordinator second-guessing the site: this site's worker refuses a
    study-less run unconditionally, so dispatching spends an allocation to learn nothing."""
    fed = federation_for(tmp_path, {"permission": "DUO:0000042"}, None)
    ok, report = enforce_data_use(fed, "fine-mapping")
    assert not ok
    assert "no-request" in report
