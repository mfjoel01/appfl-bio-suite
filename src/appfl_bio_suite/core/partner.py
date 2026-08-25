"""Generate a partner's complete setup bundle: config, docs, and mapping, all filled in.

WHY GENERATED RATHER THAN TEMPLATED
-----------------------------------
The guide this replaces handed partners a template full of ``<angle brackets>`` and
``# leave as-is`` annotations. Essentially every question a partner asked was about which
placeholder applied to them -- and one partner spent a week on a config whose placeholder
they had filled in correctly, because the real error was elsewhere and the placeholders
made it impossible to tell at a glance what was intentional.

Each partner runs exactly one experiment at one data center. There is no reason for them
to choose anything. The bundle contains their config with every value already resolved:
their client id, their data assignment, their expected sample count, their scheduler's
provider block, and the *coordinator's* identity read from federation.yaml.

THE CLOSED-SET GUARANTEE
------------------------
Documents are copied from ``docs/partner/`` and nowhere else. That directory is a closed
set: everything in it is safe to send, and nothing outside it can reach a partner.

This is a structural property, not a convention. :func:`generate_bundle` resolves every
source path and refuses to read anything outside that tree, and
``tests/test_docs.py::test_bundle_reads_only_partner_docs`` asserts it. "Did I
accidentally include my internal notes" stops being a judgment call made under time
pressure and becomes impossible.

It matters because the coordinator's own notes legitimately contain things a partner must
not see -- other partners' names, the fact that a dataset is simulated, debugging
narrative about a third site. Keeping those out by remembering to is not a control.
"""

from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from datetime import UTC
from pathlib import Path
from typing import Any

import yaml

from appfl_bio_suite.core.config import Federation
from appfl_bio_suite.core.experiments import get_spec, repo_root
from appfl_bio_suite.core.identity import build_mapping_document

__all__ = ["generate_bundle", "BundleContents", "partner_docs_root", "render_context"]


class BundleSecurityError(RuntimeError):
    """An attempt to read bundle content from outside docs/partner/."""


@dataclass
class BundleContents:
    """What a partner receives."""

    directory: Path
    files: list[Path]


def partner_docs_root() -> Path:
    """The only directory bundle documents may come from."""
    return repo_root() / "docs" / "partner"


def _assert_inside_partner_docs(path: Path) -> Path:
    """Refuse any path outside docs/partner/, symlinks included."""
    root = partner_docs_root().resolve()
    resolved = path.resolve()
    if not resolved.is_relative_to(root):
        raise BundleSecurityError(
            f"refusing to include {resolved}: bundle content may only come from "
            f"{root}.\n"
            "This is deliberate. Everything under docs/partner/ is safe to send; "
            "nothing outside it is. If a partner needs this file, move it into "
            "docs/partner/ and make sure it names no other site and no internal detail."
        )
    return resolved


def render_context(federation: Federation, experiment: str, site_key: str) -> dict[str, Any]:
    """Every value a partner-facing template may reference.

    Deliberately flat and fully resolved: a template that has to compute anything is a
    template that can compute it wrong.
    """
    spec = get_spec(experiment)
    exp = federation.experiment(experiment)
    entry = federation.experiment_site(experiment, site_key)
    site = entry.resolved_site
    coordinator = federation.coordinator

    return {
        # -- who is asking ------------------------------------------------
        # From federation.yaml, never a literal. A different coordinator running this
        # command gets correct bundles for their own federation with no code change.
        "coordinator_identity": coordinator.identity,
        "coordinator_identity_id": coordinator.identity_id or "",
        "coordinator_identity_match": coordinator.identity_regex,
        "coordinator_organization": coordinator.organization or "the coordinating site",
        "coordinator_contact": coordinator.contact or "",
        # -- the experiment -----------------------------------------------
        "experiment": experiment,
        "experiment_title": spec.title,
        "experiment_package": spec.package,
        "partner_extras": ",".join(spec.partner_extras),
        "service_account": exp.service_account,
        "endpoint_name": exp.endpoint_name,
        # -- this partner -------------------------------------------------
        "site_id": site.id,
        "site_name": site.name,
        "client_id": entry.client_id,
        "output_dir": entry.output_dir,
        # Where the partner puts their data. GWAS sites are told exactly where, because
        # the coordinator must know the path to configure the run. FLamby sites choose
        # for themselves -- FLamby resolves the location internally and never tells us --
        # so a concrete default is derived instead. It is a suggestion they can override,
        # but they are not asked to invent one: a placeholder here is a decision handed
        # to someone with less context than we have.
        "data_dir": entry.data_dir or f"{entry.output_dir.rstrip('/')}/data",
        "center": entry.center if entry.center is not None else "",
        "expected_train_samples": entry.expected_train_samples or "",
        "expected_samples": entry.expected_samples or "",
        # -- their cluster ------------------------------------------------
        "scheduler": site.scheduler,
        "provider_type": site.provider_type,
        "launcher_type": site.launcher_type,
        "scheduler_slot": site.scheduler_slot,
        "scheduler_slot_value": site.scheduler_slot_value or "",
        "account": site.account or "",
        "conda_module": site.conda_module or "",
        "interface": site.interface or "",
        # -- versions -----------------------------------------------------
        "globus_compute_version": "4.9.0",
        "python_version": "3.12",
    }


def _render(template_text: str, context: dict[str, Any]) -> str:
    from jinja2 import StrictUndefined, Template

    # StrictUndefined so a template referencing a variable that does not exist fails
    # here, loudly, rather than silently rendering an empty string into a partner's
    # config -- which is exactly the class of bug this generator exists to prevent.
    return Template(template_text, undefined=StrictUndefined, keep_trailing_newline=True).render(
        **context
    )


def _client_config_block(federation: Federation, experiment: str, site_key: str) -> dict:
    """The partner's own client entry, exactly as it will appear in the run config.

    Note what is NOT here: `dataset_path`. That path is resolved on the coordinator's
    driver, which reads the file and ships its source to the worker. Asking a partner for
    it -- as an earlier version of the setup guide did -- is asking for a value that will
    never be used, and inviting them to create a file that will never be read.
    """
    from appfl_bio_suite.core.launch import build_client_configs

    entry = federation.experiment_site(experiment, site_key)
    for client in build_client_configs(federation, experiment):
        if client["client_id"] == entry.client_id:
            block = json.loads(json.dumps(client))  # deep copy
            block["data_configs"].pop("dataset_path", None)
            return block
    raise RuntimeError(f"could not resolve a client config for site '{site_key}'")


def generate_bundle(
    federation: Federation,
    experiment: str,
    site_key: str,
    out_dir: Path | None = None,
) -> Path:
    """Write one partner's bundle and return its directory."""
    entry = federation.experiment_site(experiment, site_key)
    context = render_context(federation, experiment, site_key)

    out_dir = out_dir or (Path("local") / "partner_bundles")
    destination = out_dir / f"{experiment}_{entry.client_id}"
    destination.mkdir(parents=True, exist_ok=True)

    docs_root = partner_docs_root()
    if not docs_root.is_dir():
        raise RuntimeError(
            f"{docs_root} does not exist. Partner documents are not packaged, so bundle "
            "generation only works from a source checkout of the repository."
        )

    written: list[Path] = []

    # -- 1. the two setup documents, rendered ------------------------------
    # Part 1 is shared and identical for every experiment; Part 2 is this experiment's.
    doc_sources = [
        (docs_root / "endpoint-setup.md", "1-endpoint-setup.md"),
        (docs_root / "experiments" / f"{experiment}.md", f"2-{experiment}.md"),
    ]
    for source, name in doc_sources:
        _assert_inside_partner_docs(source)
        if not source.is_file():
            raise RuntimeError(
                f"missing partner document {source}. Every experiment must have one; "
                "see docs/README.md for the four document roles."
            )
        rendered = _render(source.read_text(encoding="utf-8"), context)
        _check_no_placeholders_remain(rendered, source)
        target = destination / name
        target.write_text(rendered, encoding="utf-8")
        written.append(target)

    # -- 2. their client config block --------------------------------------
    config_path = destination / "your-config-block.yaml"
    block = _client_config_block(federation, experiment, site_key)
    config_path.write_text(
        "# Your site's entry in the coordinator's client config.\n"
        "# Every value is already filled in -- send it back as-is unless a path is wrong\n"
        "# for your cluster. Paths here are on YOUR cluster and must be readable and\n"
        f"# writable by the service account '{context['service_account']}'.\n"
        "#\n"
        "# There is deliberately no `dataset_path`: that file is read on the\n"
        "# coordinator's machine and its contents are shipped to your worker, so it is\n"
        "# never a path on your cluster.\n\n" + yaml.safe_dump([block], sort_keys=False),
        encoding="utf-8",
    )
    written.append(config_path)

    # -- 3. their identity mapping, ready to install -----------------------
    mapping_path = destination / "example_identity_mapping_config.json"
    mapping = build_mapping_document(
        identity_match=context["coordinator_identity_match"],
        output_account=context["service_account"],
    )
    mapping_path.write_text(json.dumps(mapping, indent=2) + "\n", encoding="utf-8")
    written.append(mapping_path)

    # -- 4. their endpoint template, with their scheduler filled in --------
    template_source = repo_root() / "partners" / "templates" / "user_config_template.yaml.j2"
    if template_source.is_file():
        target = destination / "user_config_template.yaml.j2"
        target.write_text(
            _render(template_source.read_text(encoding="utf-8"), context), encoding="utf-8"
        )
        written.append(target)

    # -- 5. the validator, so they can check their own mapping -------------
    validator = repo_root() / "partners" / "templates" / "validate_mapping.py"
    if validator.is_file():
        target = destination / "validate_mapping.py"
        shutil.copy2(validator, target)
        written.append(target)

    # -- 6. a manifest, so a stale bundle is detectable --------------------
    # A partner working from an outdated copy has already cost this project real time.
    # Recording what generated this one means the question is answerable.
    _write_bundle_manifest(destination, federation, experiment, entry, written)

    return destination


def _check_no_placeholders_remain(rendered: str, source: Path) -> None:
    """A rendered partner document must contain no placeholders. That is the point."""
    import re

    leftovers = set(re.findall(r"\{\{[^}]+\}\}", rendered))
    leftovers |= set(re.findall(r"<(?:your|YOUR|a )[^>\n]{2,40}>", rendered))
    if leftovers:
        raise RuntimeError(
            f"{source} still contains placeholders after rendering: "
            f"{sorted(leftovers)}\n"
            "A generated bundle must be complete -- every placeholder left in is a "
            "question the partner has to ask. Add the missing value to federation.yaml "
            "or to render_context()."
        )


def _write_bundle_manifest(
    destination: Path, federation: Federation, experiment: str, entry, written: list[Path]
) -> None:
    from datetime import datetime

    from appfl_bio_suite import __version__
    from appfl_bio_suite.core.simulation import _git_commit

    manifest = {
        "experiment": experiment,
        "client_id": entry.client_id,
        "site": entry.site,
        "generated_at": datetime.now(UTC).isoformat(),
        "suite_version": __version__,
        "suite_commit": _git_commit(),
        "coordinator_identity": federation.coordinator.identity,
        "files": sorted(p.name for p in written),
    }
    (destination / "bundle_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
