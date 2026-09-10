"""Turn a federation config into the GA4GH facts a run needs, in one place.

Every consumer of these standards asks the same three questions -- may this site
participate, which object is it supposed to hold, and which tool version are we sending
-- and each of them is a small piece of assembly over ``federation.yaml``, a registry
file, and the installed package. Doing that assembly in ``launch.py``, again in
``preflight.py``, and again in the CLI is how the three come to disagree.

So it lives here, and the three answers are:

    resolve_data_use(fed, experiment)   -> a decision per site, and why
    resolve_drs_object(fed, entry)      -> the object that site must hold, trimmed
    resolve_tool(fed, experiment)       -> the pin, verified against this install

WHAT "TRIMMED" MEANS, AND WHY IT MATTERS
----------------------------------------
The DRS object a coordinator holds carries ``file://`` access methods pointing at the
coordinator's own filesystem. Shipping those to a partner's worker would be shipping
paths that do not exist there, and inviting a confused bug report about a directory on
somebody else's machine. What the worker needs is the identity -- the id, the self URI,
and the per-member checksums -- so that is all it gets.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

__all__ = [
    "SiteDataUse",
    "resolve_data_use",
    "resolve_drs_object",
    "resolve_tool",
    "load_registry",
    "run_provenance",
]


@dataclass
class SiteDataUse:
    """One site's data use decision, as the coordinator can compute it offline."""

    site: str
    client_id: str
    # "permitted" | "denied" | "undetermined" | "no-profile" | "no-request"
    status: str
    detail: str = ""
    decision: dict[str, Any] | None = None

    @property
    def blocking(self) -> bool:
        """True when this site must not be dispatched to.

        ``no-profile`` is not blocking: a site that has declared no terms has not refused
        anything, and its own worker will reach the same conclusion.

        ``no-request`` IS blocking, and that is not the coordinator second-guessing the
        site. This site has declared terms and the run has declared no study, which its
        worker refuses unconditionally -- so dispatching would spend a scheduler
        allocation to be told something already known here.
        """
        return self.status in ("denied", "undetermined", "no-request")

    def render(self) -> str:
        head = f"{self.client_id:12} {self.status}"
        if self.decision:
            lines = [head]
            for reason in self.decision.get("reasons", []):
                if reason["outcome"] != "permitted":
                    lines.append(
                        f"               {reason['term_id']} "
                        f"{reason['term_label']}: {reason['detail']}"
                    )
            return "\n".join(lines)
        if not self.detail:
            return head
        # Kept whole rather than truncated to its first line. These details are validation
        # errors naming a file and a field, and the first line of one of those is
        # reliably the half that says nothing.
        first, *rest = [line for line in self.detail.splitlines() if line.strip()]
        return "\n".join([f"{head}  {first}", *(f"               {line}" for line in rest)])


def resolve_data_use(federation, experiment: str) -> list[SiteDataUse]:
    """Evaluate every participating site's DUO profile against the run's request.

    Uses the coordinator's *copy* of each profile (``data_use_profile`` in
    federation.yaml). The authoritative copy is the one in the site's own bundle, which
    its worker reads -- this is the offline preview, and the two can legitimately differ
    if a site has updated its terms since the coordinator was sent a copy. That is worth
    knowing about, which is why the decision travels back in the payload and is compared
    on arrival rather than assumed to match.
    """
    from appfl_bio_suite.core.ga4gh.duo import DuoError, evaluate, load_profile, term_label

    exp = federation.experiment(experiment)
    request = federation.data_use_request(experiment)
    out: list[SiteDataUse] = []

    for entry in exp.sites:
        if entry.data_use_profile is None:
            out.append(
                SiteDataUse(
                    site=entry.site,
                    client_id=entry.client_id,
                    status="no-profile",
                    detail=(
                        "no `data_use_profile` recorded for this site, so its terms "
                        "cannot be checked here. Its own worker still enforces whatever "
                        "is in its bundle."
                    ),
                )
            )
            continue
        try:
            profile = load_profile(entry.data_use_profile)
        except DuoError as exc:
            out.append(
                SiteDataUse(
                    site=entry.site,
                    client_id=entry.client_id,
                    status="undetermined",
                    detail=str(exc),
                )
            )
            continue

        if request is None:
            out.append(
                SiteDataUse(
                    site=entry.site,
                    client_id=entry.client_id,
                    status="no-request",
                    detail=(
                        f"this site declares terms ({profile.permission}, "
                        f"{term_label(profile.permission)}) but the experiment declares no "
                        "`ga4gh.data_use_request`. Its worker will refuse the task."
                    ),
                )
            )
            continue

        decision = evaluate(profile.to_dict(), request.to_dict())
        out.append(
            SiteDataUse(
                site=entry.site,
                client_id=entry.client_id,
                status=decision["outcome"],
                decision=decision,
            )
        )
    return out


def load_registry(federation, required: bool = True):
    """Load the DRS registry this federation names, or None when it names none."""
    from appfl_bio_suite.core.ga4gh.drs import DrsRegistry

    service = federation.ga4gh.drs if federation.ga4gh else None
    if service is None or not service.registry:
        if required:
            raise FileNotFoundError(
                "no `ga4gh.drs.registry` in the federation config. Build one with "
                "`appfl-bio-suite ga4gh drs register --data-root <simulation output>`."
            )
        return None
    return DrsRegistry.load(Path(service.registry))


def resolve_drs_object(federation, entry, registry=None) -> dict[str, Any] | None:
    """The DRS object a site must hold, trimmed to what its worker can use.

    ``None`` when the site names no ``drs_uri``, which is the state a federation not
    using DRS stays in permanently and is not an error.
    """
    if not entry.drs_uri:
        return None
    if registry is None:
        try:
            registry = load_registry(federation)
        except FileNotFoundError as exc:
            # Refused rather than skipped. A site naming an object nobody can resolve is
            # a half-configured federation, and quietly dropping the verification would
            # produce a run that claims DRS provenance while checking nothing.
            raise FileNotFoundError(
                f"site '{entry.site}' declares drs_uri {entry.drs_uri}, but {exc}"
            ) from exc
    if registry is None:
        return None
    obj = registry.resolve(entry.drs_uri)
    return {
        "id": obj.id,
        "name": obj.name,
        "self_uri": obj.self_uri,
        "size": obj.size,
        "checksums": [c.model_dump() for c in obj.checksums],
        # Members, by name and content id. This is the whole of what a worker verifies
        # against; access methods are deliberately dropped (see the module docstring).
        "contents": [{"name": child.name, "id": child.id} for child in (obj.contents or [])],
    }


def resolve_tool(federation, experiment: str, verify: bool = True) -> tuple[dict | None, list[str]]:
    """The TRS pin for this run, and any problems verifying it against this install.

    Returns ``(pin dict or None, problems)``. Problems are returned rather than raised so
    that preflight can report all of them and a launch can decide -- a pin mismatch is
    fatal to a run that claims reproducibility and merely interesting to a developer
    iterating on the trainer.
    """
    from appfl_bio_suite.core.ga4gh.trs import verify_pin

    pin = federation.tool_pin(experiment)
    if pin is None:
        return None, []
    problems = verify_pin(pin, experiment) if verify else []
    return pin.model_dump(exclude_none=True), problems


def run_provenance(federation, experiment: str) -> dict[str, Any]:
    """The GA4GH facts a coordinator can state about a run before it starts.

    Written beside the results by the aggregator, and merged there with what each site
    reports having actually done.
    """
    from appfl_bio_suite import __version__
    from appfl_bio_suite.core.ga4gh import drs as drs_module
    from appfl_bio_suite.core.ga4gh import duo as duo_module
    from appfl_bio_suite.core.ga4gh import tes as tes_module
    from appfl_bio_suite.core.ga4gh import trs as trs_module

    request = federation.data_use_request(experiment)
    pin, problems = resolve_tool(federation, experiment)
    service = federation.ga4gh.drs if federation.ga4gh else None

    return {
        "suite_version": __version__,
        "specifications": {
            "drs": drs_module.DRS_VERSION,
            "tes": tes_module.TES_VERSION,
            "trs": trs_module.TRS_VERSION,
            "duo": duo_module.ontology().get("version_iri", ""),
        },
        "data_use_request": request.to_dict() if request else None,
        "tool": pin,
        "tool_pin_problems": problems,
        "drs": (
            {
                "hostname": service.hostname,
                "registry": service.registry,
                "objects": {
                    entry.client_id: entry.drs_uri
                    for entry in federation.experiment(experiment).sites
                    if entry.drs_uri
                },
            }
            if service
            else None
        ),
    }
