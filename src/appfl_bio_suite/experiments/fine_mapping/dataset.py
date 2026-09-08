"""Fine-mapping site data loader. SHIPPED TO WORKERS -- self-contained by design.

=============================================================================
THIS FILE'S SOURCE IS SENT OVER THE WIRE AND EXECUTED ON A PARTNER'S CLUSTER.
=============================================================================

It may import only the standard library and packages the partner installed with
``appfl-bio-suite[finemapping]``. It must NOT import from ``appfl_bio_suite``, and it must
not import a sibling module by bare name -- the shipped source runs from a temporary
working directory where its siblings do not exist.

Enforced by tests/test_shipped_modules.py. The rule exists because the alternative was
tried on the GWAS experiment: sibling imports meant every partner had to add a
``PYTHONPATH`` entry to their endpoint's ``worker_init``, and a missing or wrong one was
the single largest source of partner-side breakage on that project.

THE GA4GH CHECKS RUN HERE, AND THAT IS THE WHOLE POINT
------------------------------------------------------
Two of them, in this order, before a single genotype is read:

1. **DUO.** If the bundle carries a ``DATA_USE.json``, this site's own consent code is
   matched against the data use request the coordinator dispatched. A study the terms do
   not permit is refused *here*, in the site's process, on the site's hardware -- which is
   the only version of consent enforcement that means anything. The coordinator runs the
   same check offline before dispatch (``preflight --check ga4gh``), but that one is a
   courtesy that saves a queue wait, not a control.

2. **DRS.** The bundle is re-checksummed against the DRS object the coordinator says this
   site should hold. It catches the failure nothing downstream can: a site running a
   *different* bundle -- last month's, or another site's -- produces well-formed
   aggregates over the wrong individuals, pools without complaint, and is wrong in no
   visible way.

The DUO matcher below is a copy of ``core/ga4gh/duo.py::evaluate``. It has to be: this
file's source is shipped to a worker that does not have the suite installed. The copy is
held to the original by ``tests/test_ga4gh_shipped_parity.py``, which runs both over the
same case table and asserts identical decisions -- the same arrangement, for the same
reason, as the site-aggregate computation in ``trainer.py``.

WHAT THIS LOADER DOES NOT DO
----------------------------
It does not read the genotypes. A site's chromosome fileset is tens of gigabytes and the
trainer reads only the few thousand variants inside each locus window, seeking straight
to them in the ``.bed``. Loading the matrix here would exhaust a worker before any
statistics were computed.

What it does instead is validate -- eagerly, at construction, before a scheduler
allocation has been spent -- that every file the trainer will need is present and mutually
consistent: the manifest covers the fam, the phenotypes cover the manifest, and the loci
are inside the chromosome. Each of those has failed on a real bundle, and each is far
cheaper to discover here than three minutes into a run on someone else's cluster.
"""

import csv
import hashlib
import json
from datetime import date
from pathlib import Path

# The stem of the PLINK1 triple inside a site's bundle. Fixed rather than derived from
# the client id: a bundle should be inspectable without knowing which site it was cut
# for, and a stem that encodes the site name is a stem that goes stale when a site is
# renamed in federation.yaml.
PLINK_STEM = "site_genotypes"

# Files that must sit DIRECTLY in data_dir. A bundle unpacked one level too deep is the
# most common partner-side mistake, and naming the files in the error is what makes that
# diagnosable without a round trip.
REQUIRED_FILES = (
    f"{PLINK_STEM}.bed",
    f"{PLINK_STEM}.bim",
    f"{PLINK_STEM}.fam",
    "site_manifest.tsv",
    "reference_variants.tsv",
    "selected_loci.tsv",
)

# Phenotypes live one level down because there is one file per instance and a real run
# has hundreds -- 100 loci x 15 architectures x 10 replicates. Flattening them into
# data_dir would make the directory unlistable.
PHENOTYPE_DIRNAME = "phenotypes"


def _read_tsv(path, required_columns):
    """Read a small TSV into a list of dicts, checking its header names.

    Deliberately ``csv`` rather than pandas. These files are metadata -- a few thousand
    rows at most -- and the loader's job is to fail with a clear message before the
    trainer imports anything heavy. A header typo caught here names the column; the same
    typo caught inside pandas surfaces as a KeyError three functions deeper.
    """
    path = Path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [c for c in required_columns if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(
                f"{path} is missing column(s) {missing}.\n"
                f"  found: {reader.fieldnames}\n"
                f"  expected at least: {list(required_columns)}"
            )
        return list(reader)


# ---------------------------------------------------------------------------
# GA4GH DUO -- this site's consent code, enforced in this site's process
#
# A copy of core/ga4gh/duo.py::evaluate and its rule table, inlined because shipped
# source may not import the suite. tests/test_ga4gh_shipped_parity.py runs both over the
# same case table and asserts they agree; if you change a rule in one, that test fails
# until you change it in the other.
# ---------------------------------------------------------------------------

DATA_USE_FILENAME = "DATA_USE.json"

PERMITTED = "permitted"
DENIED = "denied"
UNDETERMINED = "undetermined"

# Purposes that are biomedical research (DUO:0000037 and its subclasses). What satisfies
# HMB, and nothing else does.
BIOMEDICAL_PURPOSES = frozenset({"DUO:0000037", "DUO:0000038", "DUO:0000039", "DUO:0000040"})
# Population and ancestry research. What POA permits, and what DUO:0000044 prohibits.
ANCESTRY_PURPOSES = frozenset({"DUO:0000032", "DUO:0000033"})

DUO_LABELS = {
    "DUO:0000004": "no restriction",
    "DUO:0000006": "health or medical or biomedical research",
    "DUO:0000007": "disease specific research",
    "DUO:0000011": "population origins or ancestry research only",
    "DUO:0000012": "research specific restrictions",
    "DUO:0000015": "no general methods research",
    "DUO:0000016": "genetic studies only",
    "DUO:0000018": "not for profit, non commercial use only",
    "DUO:0000019": "publication required",
    "DUO:0000020": "collaboration required",
    "DUO:0000021": "ethics approval required",
    "DUO:0000022": "geographical restriction",
    "DUO:0000024": "publication moratorium",
    "DUO:0000025": "time limit on use",
    "DUO:0000026": "user specific restriction",
    "DUO:0000027": "project specific restriction",
    "DUO:0000028": "institution specific restriction",
    "DUO:0000029": "return to database or resource",
    "DUO:0000031": "method development",
    "DUO:0000032": "population research",
    "DUO:0000033": "ancestry research",
    "DUO:0000034": "age category research",
    "DUO:0000035": "gender category research",
    "DUO:0000036": "research control",
    "DUO:0000037": "biomedical research",
    "DUO:0000038": "genetic research",
    "DUO:0000039": "drug development research",
    "DUO:0000040": "disease category research",
    "DUO:0000042": "general research use",
    "DUO:0000043": "clinical care use",
    "DUO:0000044": "population origins or ancestry research prohibited",
    "DUO:0000045": "not for profit organisation use only",
    "DUO:0000046": "non-commercial use only",
}

PERMISSION_IDS = frozenset(
    {"DUO:0000004", "DUO:0000042", "DUO:0000006", "DUO:0000007", "DUO:0000011"}
)


def evaluate_data_use(profile, request, as_of=None):
    """Match a DUO profile against a data use request. Returns a decision dict."""
    today = as_of or date.today().isoformat()
    reasons = []

    permission = str(profile.get("permission", "")).strip()
    purposes = [str(p) for p in request.get("purposes", [])]
    acknowledged = {str(a) for a in request.get("acknowledged", [])}

    def add(term_id, outcome, detail):
        reasons.append(
            {
                "term_id": term_id,
                "term_label": DUO_LABELS.get(term_id, ""),
                "outcome": outcome,
                "detail": detail,
            }
        )

    if permission not in PERMISSION_IDS:
        add(permission or "(none)", DENIED, "not a DUO data use permission")
    elif permission == "DUO:0000004":
        add(permission, PERMITTED, "no restriction on use")
    elif permission == "DUO:0000042":
        add(permission, PERMITTED, "general research use permits any research purpose")
    elif permission == "DUO:0000006":
        outside = [p for p in purposes if p not in BIOMEDICAL_PURPOSES]
        if outside:
            add(
                permission,
                DENIED,
                "permits health/medical/biomedical research only, but the study declares "
                + ", ".join(f"{p} ({DUO_LABELS.get(p, '?')})" for p in outside),
            )
        else:
            add(permission, PERMITTED, "every declared purpose is biomedical research")
    elif permission == "DUO:0000007":
        allowed = [str(v).strip().lower() for v in profile.get("permission_value", [])]
        disease = str(request.get("disease") or "").strip().lower()
        if not disease:
            add(
                permission,
                DENIED,
                "disease specific research: the study declares no disease. Permitted "
                "for: " + (", ".join(allowed) or "(unspecified)"),
            )
        elif disease not in allowed:
            add(
                permission,
                DENIED,
                f"disease specific research for {', '.join(allowed) or '(unspecified)'}; "
                f"the study declares '{disease}'",
            )
        else:
            add(permission, PERMITTED, f"disease specific research, and the study is on {disease}")
    elif permission == "DUO:0000011":
        outside = [p for p in purposes if p not in ANCESTRY_PURPOSES]
        if outside:
            add(
                permission,
                DENIED,
                "permits population origins or ancestry research only, but the study "
                "declares " + ", ".join(f"{p} ({DUO_LABELS.get(p, '?')})" for p in outside),
            )
        else:
            add(permission, PERMITTED, "the study is population/ancestry research")

    for entry in profile.get("modifiers", []):
        term_id = str(entry.get("id", "")).strip()
        values = [str(v) for v in entry.get("value", [])]
        lowered = [v.strip().lower() for v in values]
        outcome, detail = _evaluate_modifier(term_id, values, lowered, request, purposes, today)
        if outcome == UNDETERMINED and term_id in acknowledged:
            outcome = PERMITTED
            detail = f"acknowledged by the requester: {detail}"
        add(term_id, outcome, detail)

    outcomes = {r["outcome"] for r in reasons}
    if DENIED in outcomes:
        overall = DENIED
    elif UNDETERMINED in outcomes:
        overall = UNDETERMINED
    else:
        overall = PERMITTED

    return {
        "outcome": overall,
        "dataset_id": str(profile.get("dataset_id", "")),
        "requester": str(request.get("requester", "")),
        "reasons": reasons,
        "evaluated_at": today,
        "duo_version": str(profile.get("duo_version", "")),
    }


def _evaluate_modifier(term_id, values, lowered, request, purposes, today):
    """One modifier against the request. Returns ``(outcome, detail)``."""
    ask = request.get

    if term_id == "DUO:0000015":  # no general methods research
        if "DUO:0000031" in purposes:
            return DENIED, "the study declares method development (DUO:0000031)"
        return PERMITTED, "the study declares no methods development"

    if term_id == "DUO:0000016":  # genetic studies only
        outside = [p for p in purposes if p != "DUO:0000038"]
        if outside:
            return DENIED, "genetic studies only; the study also declares " + ", ".join(outside)
        return PERMITTED, "the study is genetic research"

    if term_id == "DUO:0000044":  # population/ancestry research prohibited
        offending = [p for p in purposes if p in ANCESTRY_PURPOSES]
        if offending:
            return DENIED, "population/ancestry research is prohibited: " + ", ".join(offending)
        return PERMITTED, "the study declares no population or ancestry research"

    if term_id in ("DUO:0000018", "DUO:0000046"):  # non-commercial
        if not ask("non_commercial", False):
            return DENIED, "the requester does not attest to non-commercial use"
        if term_id == "DUO:0000018" and not ask("not_for_profit_organisation", False):
            return DENIED, "the requester does not attest to being a not-for-profit organisation"
        return PERMITTED, "attested"

    if term_id == "DUO:0000045":  # not for profit organisation use only
        if not ask("not_for_profit_organisation", False):
            return DENIED, "the requester does not attest to being a not-for-profit organisation"
        return PERMITTED, "attested"

    if term_id == "DUO:0000019":  # publication required
        if not ask("publication_agreed", False):
            return DENIED, "the requester has not agreed to publish results"
        return PERMITTED, "attested"

    if term_id == "DUO:0000020":  # collaboration required
        if not ask("collaboration_agreed", False):
            return DENIED, "the requester has not agreed to collaborate" + (
                f" with {', '.join(values)}" if values else ""
            )
        return PERMITTED, "attested" + (f" ({', '.join(values)})" if values else "")

    if term_id == "DUO:0000021":  # ethics approval required
        approval = str(ask("ethics_approval") or "").strip()
        if not approval:
            return DENIED, "no ethics/IRB approval reference was supplied"
        return PERMITTED, f"approval {approval}"

    if term_id == "DUO:0000022":  # geographical restriction
        where = str(ask("geography") or "").strip().lower()
        if not values:
            return UNDETERMINED, "restricted geographically, but no region was specified"
        if not where:
            return DENIED, "restricted to " + ", ".join(values) + "; the request names no region"
        if where not in lowered:
            return DENIED, f"restricted to {', '.join(values)}; the request names '{where}'"
        return PERMITTED, f"within {where}"

    if term_id == "DUO:0000024":  # publication moratorium
        if not ask("moratorium_accepted", False):
            return DENIED, "the requester has not accepted the publication moratorium" + (
                f" (until {values[0]})" if values else ""
            )
        return PERMITTED, "accepted" + (f" until {values[0]}" if values else "")

    if term_id == "DUO:0000025":  # time limit on use
        if not values:
            return UNDETERMINED, "a time limit applies, but no date was specified"
        limit = values[0].strip()
        if today > limit:
            return DENIED, f"use expired on {limit} (evaluating as of {today})"
        return PERMITTED, f"within the limit, which ends {limit}"

    if term_id == "DUO:0000026":  # user specific restriction
        who = str(ask("requester") or "").strip().lower()
        if not values:
            return UNDETERMINED, "restricted to approved users, but none were listed"
        if who not in lowered:
            return DENIED, f"'{who}' is not among the approved users"
        return PERMITTED, f"'{who}' is an approved user"

    if term_id == "DUO:0000027":  # project specific restriction
        project = str(ask("project") or "").strip().lower()
        if not values:
            return UNDETERMINED, "restricted to approved projects, but none were listed"
        if project not in lowered:
            return DENIED, (
                f"restricted to {', '.join(values)}; the request names "
                f"'{project or '(no project)'}'"
            )
        return PERMITTED, f"project '{project}' is approved"

    if term_id == "DUO:0000028":  # institution specific restriction
        institution = str(ask("institution") or "").strip().lower()
        if not values:
            return UNDETERMINED, "restricted to approved institutions, but none were listed"
        if institution not in lowered:
            return DENIED, (
                f"restricted to {', '.join(values)}; the request names "
                f"'{institution or '(no institution)'}'"
            )
        return PERMITTED, f"institution '{institution}' is approved"

    if term_id == "DUO:0000029":  # return to database or resource
        if not ask("return_to_resource_agreed", False):
            return DENIED, "the requester has not agreed to return derived data to the resource"
        return PERMITTED, "attested"

    if term_id == "DUO:0000043":  # clinical care use
        return PERMITTED, "clinical care use is additionally permitted; no effect on research use"

    if term_id == "DUO:0000012":  # research specific restrictions
        return UNDETERMINED, (
            "free-text research restriction"
            + (f": {'; '.join(values)}" if values else "")
            + ". A person must read this and record the decision in `acknowledged`."
        )

    return UNDETERMINED, (
        "this matcher has no rule for that term, so it cannot be treated as satisfied"
    )


def render_decision(decision):
    """The decision as a block of text, for a log or an exception message."""
    mark = {PERMITTED: "ok", DENIED: "DENY", UNDETERMINED: "????"}
    lines = [
        f"{decision['outcome'].upper()}: {decision['dataset_id'] or 'this dataset'} "
        f"for {decision['requester'] or '(no requester declared)'}"
    ]
    for reason in decision["reasons"]:
        lines.append(
            f"  [{mark.get(reason['outcome'], '?'):>4}] {reason['term_id']} "
            f"{reason['term_label']}: {reason['detail']}"
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# GA4GH DRS -- is this the bundle the coordinator cut for this site?
# ---------------------------------------------------------------------------

_HASH_CHUNK = 1 << 20

# Excluded from `metadata` verification: the genotype matrix, which is the only file
# large enough for hashing to cost real time. Everything else in a bundle is metadata,
# phenotypes, or the variant list, and hashing all of them is seconds.
_LARGE_FILES = (f"{PLINK_STEM}.bed",)


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(_HASH_CHUNK)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def verify_against_drs(data_dir, drs_object, mode="metadata"):
    """Re-checksum a bundle against its DRS object. Returns a report dict.

    ``drs_object`` is the coordinator's DRS record for this site's bundle, inlined into
    the client config -- the worker resolves nothing over the network, because a partner's
    compute node has no outbound access and a check that needs one is a check that gets
    switched off.

    Both directions are checked. A missing member is a truncated or wrong bundle; an
    *extra* file means this directory is not the object it claims to be, which per-file
    checksums alone would never notice.
    """
    report = {
        "mode": mode,
        "object_id": (drs_object or {}).get("id", ""),
        "self_uri": (drs_object or {}).get("self_uri", ""),
        "checked": 0,
        "skipped": [],
        "problems": [],
    }
    if not drs_object or mode == "off":
        report["mode"] = "off"
        return report

    contents = {entry["name"]: entry for entry in drs_object.get("contents", []) or []}
    if not contents:
        report["problems"].append(
            "the DRS object for this bundle lists no contents, so there is nothing to "
            "verify against. Rebuild the registry with `ga4gh drs register`."
        )
        return report

    present = {p.name for p in Path(data_dir).iterdir()}
    # DATA_USE.json is never part of the object: the profile describes the bundle and
    # carries its drs_uri, so it cannot be inside the checksum that names it. The
    # registrar excludes it for the same reason (core/ga4gh/drs.py), which is what keeps
    # a re-registered directory minting the same id.
    extra = sorted(present - set(contents) - {DATA_USE_FILENAME})
    missing = sorted(set(contents) - present)
    if missing:
        report["problems"].append(f"missing from this bundle: {', '.join(missing)}")
    if extra:
        report["problems"].append(
            f"present in this bundle but not in the DRS object: {', '.join(extra)}"
        )

    for name, entry in sorted(contents.items()):
        path = Path(data_dir) / name
        if not path.exists():
            continue
        if path.is_dir():
            # Nested bundles (phenotypes/) carry their own Merkle id; verifying them
            # would need the child object, which is not inlined. Named as skipped rather
            # than silently passed.
            report["skipped"].append(f"{name}/ (nested DRS bundle)")
            continue
        if mode == "metadata" and name in _LARGE_FILES:
            report["skipped"].append(f"{name} (mode=metadata)")
            continue
        expected = entry.get("id") or ""
        actual = file_sha256(path)
        report["checked"] += 1
        if expected and actual != expected:
            report["problems"].append(
                f"{name}: sha-256 mismatch\n"
                f"    expected {expected}\n"
                f"    actual   {actual}"
            )
    return report


class SiteFineMappingDataset:
    """One site's cohort for fine-mapping: paths, ancestry composition, and loci.

    Validation is strict and happens at construction, for the reason above. Everything it
    reads is small; the genotype matrix is never touched.
    """

    def __init__(
        self,
        data_dir,
        site_id,
        instances=None,
        data_use_request=None,
        enforce_data_use=True,
        drs_object=None,
        verify_bundles="metadata",
    ):
        self.data_dir = Path(data_dir).resolve()
        self.site_id = str(site_id)

        if not self.data_dir.is_dir():
            raise FileNotFoundError(
                f"{self.site_id}: data_dir does not exist: {self.data_dir}\n"
                "\n"
                "This path is resolved on YOUR cluster, by the worker. It must be an "
                "absolute path readable by the service account the endpoint maps tasks "
                "to -- not a path on the coordinator's machine."
            )

        missing = [name for name in REQUIRED_FILES if not (self.data_dir / name).is_file()]
        if missing:
            present = sorted(p.name for p in self.data_dir.iterdir())
            raise FileNotFoundError(
                f"{self.site_id}: missing required input files in {self.data_dir}\n"
                "  missing: " + ", ".join(missing) + "\n"
                "  present: " + (", ".join(present) if present else "(directory is empty)") + "\n\n"
                "Unpack the bundle you were sent so that these files sit directly in "
                "data_dir, not in a nested subdirectory."
            )

        # -- GA4GH, before any data is read -----------------------------------
        #
        # Order matters and is not cosmetic. Consent first: a study this site's terms do
        # not permit must not get as far as opening a phenotype file. Identity second:
        # there is no point checksumming a bundle nobody was allowed to compute over.
        self.data_use = self._check_data_use(data_use_request, enforce_data_use)
        self.drs = self._verify_bundle(drs_object, verify_bundles)

        self.plink_prefix = self.data_dir / PLINK_STEM
        self.plink_bed = self.data_dir / f"{PLINK_STEM}.bed"
        self.plink_bim = self.data_dir / f"{PLINK_STEM}.bim"
        self.plink_fam = self.data_dir / f"{PLINK_STEM}.fam"
        self.manifest_path = self.data_dir / "site_manifest.tsv"
        self.reference_path = self.data_dir / "reference_variants.tsv"
        self.loci_path = self.data_dir / "selected_loci.tsv"
        self.phenotype_dir = self.data_dir / PHENOTYPE_DIRNAME

        # -- the fam: one line per individual, in .bed column order ------------
        with self.plink_fam.open("r", encoding="utf-8") as handle:
            fam_ids = [line.split()[1] for line in handle if line.strip()]
        self.sample_size = len(fam_ids)
        if self.sample_size == 0:
            raise ValueError(
                f"{self.site_id}: {self.plink_fam} lists no samples. The bundle is "
                "empty or truncated -- check the transfer."
            )

        # -- the manifest: which of those individuals are which ancestry ------
        manifest = _read_tsv(self.manifest_path, ("FID", "IID", "superpopulation"))
        fam_set = set(fam_ids)
        manifest_ids = {row["IID"] for row in manifest}
        unknown = manifest_ids - fam_set
        if unknown:
            raise ValueError(
                f"{self.site_id}: site_manifest.tsv names {len(unknown)} individual(s) "
                f"that are not in {self.plink_fam.name}, e.g. {sorted(unknown)[:3]}.\n"
                "The manifest and the fileset must describe the same cohort. This "
                "usually means two bundles' files were mixed."
            )
        unlabelled = fam_set - manifest_ids
        if unlabelled:
            raise ValueError(
                f"{self.site_id}: {len(unlabelled)} individual(s) in "
                f"{self.plink_fam.name} have no ancestry label in site_manifest.tsv, "
                f"e.g. {sorted(unlabelled)[:3]}.\n"
                "Every individual must be assigned to a superpopulation: the site stage "
                "partitions its cohort by ancestry and an unlabelled individual would be "
                "silently dropped from every aggregate."
            )

        # Ancestries this site actually holds, and how many of each. Reported to the
        # coordinator so it can check its column plan against reality rather than
        # assuming the composition it recorded when the bundle was cut.
        composition = {}
        for row in manifest:
            pop = row["superpopulation"]
            composition[pop] = composition.get(pop, 0) + 1
        self.composition = dict(sorted(composition.items()))

        # -- the loci to fine-map ---------------------------------------------
        self.loci = _read_tsv(self.loci_path, ("locus_id", "chrom", "start_bp", "end_bp"))
        if not self.loci:
            raise ValueError(f"{self.site_id}: selected_loci.tsv lists no loci")

        # -- the instances (locus x architecture x replicate) ------------------
        # A site holds one .pheno per instance. The coordinator decides which subset a
        # run covers; absent an explicit list, everything present is offered.
        if not self.phenotype_dir.is_dir():
            raise FileNotFoundError(
                f"{self.site_id}: no {PHENOTYPE_DIRNAME}/ directory in {self.data_dir}.\n"
                "It holds one <instance>.pheno per (locus, architecture, replicate)."
            )
        available = sorted(p.stem for p in self.phenotype_dir.glob("*.pheno"))
        if not available:
            raise FileNotFoundError(
                f"{self.site_id}: {self.phenotype_dir} contains no .pheno files."
            )
        self.available_instances = available

        if instances is None:
            self.instances = available
        else:
            requested = [str(i) for i in instances]
            absent = [i for i in requested if i not in set(available)]
            if absent:
                raise FileNotFoundError(
                    f"{self.site_id}: the coordinator asked for {len(absent)} instance(s) "
                    f"this site has no phenotype for, e.g. {absent[:3]}.\n"
                    "Every site must hold the same instance set -- the coordinator pools "
                    "one instance's aggregates across all of them, and a site missing an "
                    "instance would silently shrink that instance's cohort."
                )
            self.instances = requested

    # -- GA4GH ------------------------------------------------------------

    def _check_data_use(self, request, enforce):
        """Match this bundle's DUO terms against the run's request, or explain why not.

        Four cases, and each behaves differently on purpose:

        * terms present, request present -- evaluate. A refusal raises PermissionError.
        * terms present, NO request      -- refuse. A dataset that declares its terms
          must not be computed over by a study that declared nothing about itself; that
          is the whole content of a consent code.
        * NO terms, request present      -- proceed, and say so. Nothing to check
          against; the run records that this site declared none.
        * neither                        -- proceed silently. A federation not using DUO
          is unaffected by any of this.
        """
        profile_path = self.data_dir / DATA_USE_FILENAME
        has_profile = profile_path.is_file()

        if not has_profile:
            return {
                "status": "no-profile",
                "detail": (
                    f"this bundle carries no {DATA_USE_FILENAME}, so no data use terms "
                    "were checked"
                ),
            }

        try:
            profile = json.loads(profile_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise ValueError(
                f"{self.site_id}: {profile_path} is present but unreadable ({exc}).\n"
                "It declares this dataset's use terms, and a bundle whose terms cannot "
                "be read must not be computed over."
            ) from exc

        if not request:
            raise PermissionError(
                f"{self.site_id}: this bundle declares data use terms in "
                f"{DATA_USE_FILENAME}, but the run that dispatched this task declared no "
                "data use request.\n\n"
                + _profile_summary(profile)
                + "\n\nThe coordinator must declare the study's purpose under "
                "`experiments.<name>.ga4gh.data_use_request` in federation.yaml. A "
                "dataset with terms cannot be used by a study that states nothing about "
                "itself."
            )

        decision = evaluate_data_use(profile, request)
        decision["status"] = decision["outcome"]
        decision["enforced"] = bool(enforce)
        decision["profile_path"] = str(profile_path)

        if decision["outcome"] != PERMITTED:
            if enforce:
                raise PermissionError(
                    f"{self.site_id}: this site's data use terms do not permit this "
                    f"study.\n\n{render_decision(decision)}\n\n"
                    "No data was read. Resolve this with the data steward for this "
                    "dataset -- either the study's declared purpose is wrong, or this "
                    "site should not be in this run."
                )
            # Not enforcing is a coordinator's choice about their own dry run; the
            # decision still travels in the payload, so a result produced this way is
            # not mistakable for a permitted one.
            decision["detail"] = "not permitted, and enforcement was disabled for this run"
        return decision

    def _verify_bundle(self, drs_object, mode):
        """Re-checksum this bundle against its DRS record. Mismatches are fatal."""
        report = verify_against_drs(self.data_dir, drs_object, mode or "off")
        if report["problems"]:
            raise ValueError(
                f"{self.site_id}: this directory is not the DRS object the coordinator "
                f"expects ({report['self_uri'] or report['object_id'] or 'unknown'}).\n\n"
                + "\n".join(f"  {p}" for p in report["problems"])
                + "\n\nThe usual cause is a bundle from a different simulation run, or "
                "a transfer that did not finish. Re-fetch the bundle you were sent for "
                "this run; do not repair it file by file."
            )
        return report

    def __len__(self):
        return self.sample_size

    def __repr__(self):
        return (
            f"SiteFineMappingDataset({self.site_id!r}, n={self.sample_size}, "
            f"pops={list(self.composition)}, loci={len(self.loci)}, "
            f"instances={len(self.instances)}, dir={self.data_dir})"
        )


def _profile_summary(profile):
    """A one-block rendering of a profile, for an error a partner will read."""
    lines = [
        f"  dataset: {profile.get('dataset_id', '(unnamed)')}",
        f"  permission: {profile.get('permission', '?')} "
        f"{DUO_LABELS.get(profile.get('permission', ''), '')}",
    ]
    for entry in profile.get("modifiers", []) or []:
        value = ", ".join(str(v) for v in entry.get("value", []) or [])
        lines.append(
            f"  modifier: {entry.get('id')} {DUO_LABELS.get(entry.get('id'), '')}"
            + (f" ({value})" if value else "")
        )
    return "\n".join(lines)


def get_dataset(
    data_dir,
    site_id,
    instances=None,
    data_use_request=None,
    enforce_data_use=True,
    drs_object=None,
    verify_bundles="metadata",
    **kwargs,
):
    """APPFL dataset entry point.

    Returns ``(train_dataset, val_dataset)``. There is no validation split, for the same
    reason the GWAS experiment has none: this is single-round aggregate federated
    learning, not iterative training, so there is no global model to score against a
    held-out set each round. Fine-mapping accuracy is measured against the simulated
    ground truth at the coordinator, which is the only place that knows it.
    """
    return (
        SiteFineMappingDataset(
            data_dir=data_dir,
            site_id=site_id,
            instances=instances,
            data_use_request=data_use_request,
            enforce_data_use=enforce_data_use,
            drs_object=drs_object,
            verify_bundles=verify_bundles,
        ),
        None,
    )
