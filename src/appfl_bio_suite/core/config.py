"""Load, validate, and resolve ``federation.yaml``.

This is the single source of truth for who the coordinator is, who the partners are, and
what each of them is doing. Everything else in the suite -- config generation, partner
bundles, preflight, the CLI -- reads from here and never from a literal.

That is what makes the suite reusable by someone other than its author: standing up a new
federation is filling in one file, not editing code. If any component in this package
needs a coordinator identity, an endpoint UUID, or a site name, it takes it from a
:class:`Federation` instance.

Validation is deliberately strict and the messages name the fix. A malformed federation
config surfaces here, in a second, rather than as a 422 from a partner's endpoint three
time zones away.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from appfl_bio_suite.core.ga4gh.duo import DataUseRequest
from appfl_bio_suite.core.ga4gh.trs import ToolPin

__all__ = [
    "Federation",
    "Coordinator",
    "GA4GHServices",
    "DrsService",
    "TesService",
    "TrsService",
    "ExperimentGA4GH",
    "CoordinatorEndpoint",
    "Location",
    "Site",
    "Experiment",
    "ExperimentSite",
    "FederationError",
    "load_federation",
    "find_federation_file",
    "EXPERIMENT_NAMES",
]

SCHEMA_VERSION = 1

# Canonical experiment identifiers. Directory names under docs/experiments/ and
# docs/partner/experiments/ must match these, which tests assert.
EXPERIMENT_NAMES = ("flamby-heart-disease", "gwas", "fine-mapping")

# Where load_federation() looks, in order. The real config lives in local/ because it
# carries live endpoint UUIDs; see local/README.md.
_SEARCH_PATH = (
    Path("local/federation.yaml"),
    Path("federation.yaml"),
)

_ENV_VAR = "APPFL_BIO_SUITE_FEDERATION"

_UUID_RE = re.compile(
    r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-"
    r"[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"
)


class FederationError(Exception):
    """A federation config is missing, malformed, or internally inconsistent."""


# The ONLY backslash escapes the Globus expression mapper accepts, read straight out of
# globus_identity_mapping/mappers/expression.py::_compile_match:
#
#     bad_escape = re.search(r"(\\[^.?*|()\\0-9])", match)
#     if bad_escape: raise InvalidMappingError(...)
#
# So `\.` is fine and `\-` is a hard error -- the endpoint rejects the whole mapping
# document. This is why the identity cannot simply be passed through `re.escape`, which
# escapes `-` and would ship every partner a mapping file their endpoint refuses to load.
_VALID_ESCAPES = ".?*|()\\"

# Characters the mapper escapes into literals *for* you, if you leave them bare:
#
#     match = re.sub(r"(?<!\\)([\^\$\+\{\}\[\]])", r"\\\1", match)
#
# These must be left ALONE. Escaping them yourself trips the bad-escape check above,
# while leaving them bare gets the literal match you wanted. An identity containing `+`
# (`user+tag@example.org`) therefore works correctly with no special handling.
_MAPPER_ESCAPES_FOR_YOU = "^$+{}[]"


def escape_identity_for_mapping(identity: str) -> str:
    """Render an identity for a partner's identity-mapping ``match`` field.

    Escapes exactly the characters the mapper permits escaping, leaves alone the ones it
    escapes for you, and adds no anchors. See :attr:`Coordinator.identity_regex` for why
    the anchors matter and :func:`compile_match_expression` to check the result against
    the real library.
    """
    identity = identity.strip()

    # An identity that starts with '^' or ends with '$' is almost always a regex someone
    # pasted instead of the identity itself. The mapper would literalize the anchors and
    # the mapping would silently never match -- which is precisely the failure this whole
    # code path exists to prevent, so it is worth refusing loudly.
    if identity.startswith("^") or identity.endswith("$"):
        raise FederationError(
            f"coordinator.identity is {identity!r}, which looks like a regex rather than "
            "an identity. Use the bare identity (e.g. 'you@example.org'). The mapper "
            "anchors the pattern itself, and a '^' or '$' you write becomes a literal "
            "character that can never match."
        )

    out = []
    for char in identity:
        if char in _VALID_ESCAPES:
            out.append("\\" + char)
        else:
            # Includes _MAPPER_ESCAPES_FOR_YOU, which must stay bare, and ordinary
            # characters like '-', '@', '_' and '/' which need no escaping at all.
            out.append(char)
    return "".join(out)


def compile_match_expression(match: str) -> str:
    """Compile a ``match`` value through the real Globus mapper and return the pattern.

    Uses the same library the endpoint uses, so this is authoritative rather than a
    reimplementation of its rules. Raises :class:`FederationError` if the mapper rejects
    the expression.
    """
    try:
        from globus_identity_mapping.mappers.expression import ExpressionIdentityMapping
    except ImportError as exc:  # pragma: no cover - depends on install extras
        raise FederationError(
            "globus-identity-mapping is not installed. It ships with "
            "globus-compute-endpoint; run this in the same environment as the endpoint."
        ) from exc

    try:
        return ExpressionIdentityMapping._compile_match(match).pattern
    except Exception as exc:
        raise FederationError(f"the Globus mapper rejects match={match!r}: {exc}") from exc


class _Strict(BaseModel):
    """Reject unknown keys.

    A typo in a YAML key is otherwise silently ignored, and the resulting behaviour --
    a default quietly applied where the user thought they had set something -- is
    exactly the kind of failure that only shows up on a partner's cluster.
    """

    model_config = ConfigDict(extra="forbid")


class Location(_Strict):
    """Where a participant physically is, for the network map.

    DECLARED, NOT DETECTED. The obvious alternative -- resolve each participant's IP
    through a geolocation service at run time -- does not work here and should not be
    made to. A partner's compute node has no outbound web access, its public IP belongs
    to the institution's border router rather than to the machine, and asking a site to
    call an external service so that a coordinator can draw a dot is a data-governance
    conversation nobody wants to have for a dot. The coordinator already knows where
    their partners are; writing it down once is both cheaper and more accurate.

    There is no country-centroid fallback either. A site with no coordinates is reported
    as unplaced rather than drawn somewhere plausible -- see
    :func:`appfl_bio_suite.core.watch.unplaced_sites`. Inventing a position would put a
    partner's name on a map at a location they never gave, which is worse than a gap.
    """

    lat: float
    lng: float

    # City is display-only. The map shows "city, country" under each marker; country
    # comes from the site itself, which already declares one.
    city: str | None = None

    @field_validator("lat")
    @classmethod
    def _lat_range(cls, value: float) -> float:
        if not -90.0 <= value <= 90.0:
            raise ValueError(f"latitude {value} is out of range (-90 to 90).")
        return value

    @field_validator("lng")
    @classmethod
    def _lng_range(cls, value: float) -> float:
        if not -180.0 <= value <= 180.0:
            raise ValueError(
                f"longitude {value} is out of range (-180 to 180). Note the order: "
                "`lat` first, then `lng`. Swapping them is the usual cause -- a "
                "longitude in the latitude slot is often still in range and produces a "
                "marker in the wrong hemisphere with no error at all."
            )
        return value


class CoordinatorEndpoint(_Strict):
    """The coordinator's own compute endpoint, when they also train."""

    name: str
    uuid: str
    scheduler: Literal["pbspro", "slurm", "none"] = "none"
    queue: str | None = None
    partition: str | None = None
    account: str | None = None

    @field_validator("uuid")
    @classmethod
    def _uuid_shape(cls, value: str) -> str:
        return _require_uuid(value, "coordinator.endpoint.uuid")


class Coordinator(_Strict):
    """The identity that dispatches, and the machine it dispatches from."""

    identity: str
    identity_id: str | None = None
    organization: str | None = None
    contact: str | None = None
    # Where the driver runs. Drawn as the hub of the network map; every partner marker
    # is joined back to it. Optional -- omitting it costs the map its centre, nothing else.
    location: Location | None = None
    endpoint: CoordinatorEndpoint | None = None
    host_check: str | None = None
    host_check_enforce: bool = False

    @field_validator("identity")
    @classmethod
    def _identity_present(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError(
                "coordinator.identity is required -- it is the string every partner "
                "authorizes in their endpoint's identity mapping. Find it with "
                "`globus whoami`."
            )
        if value.startswith("<") or "example" in value and "example-university" not in value:
            # Only a nudge; the shipped example deliberately uses example-university.edu.
            pass
        return value

    @field_validator("identity_id")
    @classmethod
    def _identity_id_shape(cls, value: str | None) -> str | None:
        if value is None:
            return None
        return _require_uuid(value, "coordinator.identity_id")

    @field_validator("host_check")
    @classmethod
    def _host_check_compiles(cls, value: str | None) -> str | None:
        if value is None:
            return None
        try:
            re.compile(value)
        except re.error as exc:
            raise ValueError(f"coordinator.host_check is not a valid regex: {exc}") from exc
        return value

    @property
    def identity_regex(self) -> str:
        """The identity as it must appear in a partner's mapping ``match`` field.

        Two things make this trickier than ``re.escape``, and both have bitten this
        project:

        1. **No anchors.** The Globus expression mapper escapes ``^`` and ``$`` into
           *literal* characters and then wraps the pattern in its own ``^...$``. Writing
           ``^me@example\\.org$`` therefore compiles to a pattern matching only the
           literal 17-character string, and every task submission fails with
           ``422 ... Identity failed to map to a local user name``. Omitting the anchors
           does not loosen the match; the mapper anchors it for you.

        2. **Minimal escaping.** ``re.escape`` also escapes ``-``, producing
           ``a\\-b@example\\.org``. That is harmless to the regex engine but it looks
           wrong to a partner reading their own config, and a partner who "fixes" it by
           hand is a partner debugging a mapping failure. Only characters that actually
           need escaping are escaped.

        See core/identity.py, which validates a real mapping file against this.
        """
        return escape_identity_for_mapping(self.identity)


class Site(_Strict):
    """An institution in the federation. Metadata only; participation is per-experiment."""

    id: str
    name: str
    country: str | None = None
    # Where this institution is. Only the network map reads it; nothing about a run
    # depends on it. See :class:`Location` for why it is declared rather than detected.
    location: Location | None = None
    scheduler: Literal["pbspro", "slurm"] = "slurm"
    queue: str | None = None
    partition: str | None = None
    account: str | None = None
    conda_module: str | None = None
    interface: str | None = None

    @model_validator(mode="after")
    def _scheduler_slot(self) -> Site:
        """SLURM sites use `partition`, PBS sites use `queue`. Catch the mix-up early."""
        if self.scheduler == "slurm" and self.queue and not self.partition:
            raise ValueError(
                f"site '{self.id}' is scheduler: slurm but sets `queue`. "
                "SLURM sites use `partition`. (PBS sites use `queue`.)"
            )
        if self.scheduler == "pbspro" and self.partition and not self.queue:
            raise ValueError(
                f"site '{self.id}' is scheduler: pbspro but sets `partition`. "
                "PBS sites use `queue`. (SLURM sites use `partition`.)"
            )
        return self

    @property
    def scheduler_slot(self) -> str:
        """The provider key this site's scheduler expects: 'queue' or 'partition'."""
        return "queue" if self.scheduler == "pbspro" else "partition"

    @property
    def scheduler_slot_value(self) -> str | None:
        return self.queue if self.scheduler == "pbspro" else self.partition

    @property
    def provider_type(self) -> str:
        """The parsl provider class name for this scheduler."""
        return "PBSProProvider" if self.scheduler == "pbspro" else "SlurmProvider"

    @property
    def launcher_type(self) -> str:
        return "MpiExecLauncher" if self.scheduler == "pbspro" else "SrunLauncher"


class ExperimentSite(_Strict):
    """One site's participation in one experiment.

    Everything the bundle generator needs to emit that partner's fully-resolved config
    and setup docs, with no placeholder left for them to interpret.
    """

    site: str
    client_id: str
    endpoint_uuid: str
    output_dir: str

    # FLamby
    center: int | None = None
    expected_train_samples: int | None = None

    # GWAS and fine-mapping. Both distribute per-site bundles the coordinator generated,
    # so both name a directory on the partner's cluster and a sample count to check it
    # against. Shared rather than duplicated because they mean the same thing.
    data_dir: str | None = None
    expected_samples: int | None = None

    # -- GA4GH ------------------------------------------------------------
    #
    # Which DRS object this site's `data_dir` is supposed to contain. The coordinator
    # resolves it against the registry and ships the object's checksums to the worker,
    # which re-checks them before computing. That is what makes "did this site unpack the
    # bundle I cut for THIS site" answerable -- the failure it catches is a site running
    # last month's bundle, which produces well-formed aggregates over the wrong people
    # and is invisible in every downstream number.
    drs_uri: str | None = None

    # The coordinator's copy of this site's DUO profile, for the offline check. The
    # authoritative copy is the one inside the site's own bundle, which its worker reads;
    # this is what lets preflight give the same answer before a queue wait.
    data_use_profile: str | None = None

    # TES path only. A per-site service URL (falls back to the federation-wide one) and
    # the URL its outputs are staged to -- which the driver must be able to read, since
    # that is where the aggregates come back from.
    tes_url: str | None = None
    tes_outputs_url: str | None = None

    @field_validator("endpoint_uuid")
    @classmethod
    def _uuid_shape(cls, value: str) -> str:
        return _require_uuid(value, "endpoint_uuid")

    @field_validator("output_dir", "data_dir")
    @classmethod
    def _absolute(cls, value: str | None) -> str | None:
        """Worker-side paths must be absolute.

        A Globus Compute worker runs in the endpoint's task working directory, not in a
        repo root, so a relative path resolves somewhere nobody intended. This is a
        classic mid-run failure: the task dispatches fine and fails minutes later.
        """
        if value is None:
            return None
        if not value.startswith("/"):
            raise ValueError(
                f"'{value}' must be an absolute path. It is resolved on the partner's "
                "worker, which runs in the endpoint's task working directory -- a "
                "relative path will not resolve where you expect."
            )
        return value

    # `resolved_site` is attached by Federation after cross-referencing `sites`.
    @property
    def resolved_site(self) -> Site:
        site = getattr(self, "_resolved_site", None)
        if site is None:
            raise FederationError(
                f"site '{self.site}' was not resolved -- load the config through "
                "load_federation() rather than constructing models directly."
            )
        return site


class DrsService(_Strict):
    """Where this federation's data objects are named and, optionally, served.

    ``hostname`` is the authority half of every ``drs://`` URI the federation uses. It is
    required as soon as DRS is used at all, because a DRS id is only unique within its
    service -- two coordinators who both wrote ``drs://localhost/<id>`` would have
    produced URIs that collide in provenance and resolve to different bytes.
    """

    hostname: str
    # The registry `simulate` wrote, or `ga4gh drs register` built. Coordinator-side path.
    registry: str | None = None
    # Base URL where this registry is served, when it is. Adds an `https` access method
    # to every object.
    https_base: str | None = None
    # A Globus collection UUID holding the bundles, adding a `globus` access method --
    # which is how a partner actually re-fetches one in this federation.
    globus_collection: str | None = None

    @field_validator("hostname")
    @classmethod
    def _hostname_shape(cls, value: str) -> str:
        value = value.strip().rstrip("/")
        if "://" in value:
            raise ValueError(
                f"drs.hostname is '{value}', which is a URL. It is the authority half of "
                "a DRS URI -- 'drs.example.org', not 'https://drs.example.org'. Put the "
                "URL in `https_base` if you serve the registry."
            )
        if not value:
            raise ValueError("drs.hostname must not be empty")
        return value


class TrsService(_Strict):
    """The tool registry this federation's pins resolve against."""

    # Dockstore's is https://dockstore.org/api . Optional: a pin verifies against the
    # installed package with no registry at all, which is what preflight does offline.
    registry_url: str | None = None
    # Where `ga4gh trs publish` writes the static tree.
    publish_dir: str | None = None


class TesService(_Strict):
    """A TES service, for the TES execution path."""

    url: str | None = None
    poll_seconds: float = 15.0
    # A site stage at production scale is minutes to tens of minutes. Two hours is a
    # ceiling that catches a wedged task without cancelling a slow one.
    timeout_seconds: float = 7200.0
    cpu_cores: int = 4
    ram_gb: float = 16.0
    disk_gb: float = 64.0
    preemptible: bool = False


class GA4GHServices(_Strict):
    """Federation-wide GA4GH service configuration.

    Every field is optional. A federation that uses none of this runs exactly as it did
    before, which is the property that let the four standards be added to a working
    system rather than replacing it.
    """

    drs: DrsService | None = None
    trs: TrsService | None = None
    tes: TesService | None = None


class ExperimentGA4GH(_Strict):
    """One experiment's GA4GH settings: what it may do, and with which tool.

    WHICH EXPERIMENTS THIS ACTUALLY DOES SOMETHING FOR
    ---------------------------------------------------
    The coordinator-side half works for any experiment: the data use gate, the preflight
    checks, and the tool pin all read this block and nothing experiment-specific.

    The site-side half -- a worker refusing a study its terms do not permit, and
    re-checksumming its bundle before computing -- is implemented in the fine-mapping
    loader only. Declaring this block on another experiment therefore gets you the checks
    a coordinator can make and none of the ones that matter, which is worth knowing
    before relying on it. Adding them elsewhere is per-experiment work in that
    experiment's shipped loader; see ``experiments/fine_mapping/dataset.py``.
    """

    # The study, in DUO terms. Evaluated against every participating site's profile.
    # One request for the whole experiment, deliberately: a coordinator who could vary
    # the declared purpose per site to get past a refusal would have a gate that gates
    # nothing.
    data_use_request: DataUseRequest | None = None

    # Refuse to launch when a site's terms do not permit the request. Default on: a
    # federation that declares data use terms and then dispatches anyway has written
    # documentation, not a control. Turning it off leaves the SITE-side check in place --
    # that one is not the coordinator's to disable.
    enforce_data_use: bool = True

    # The TRS pin. Verified against the installed package by preflight, recorded in the
    # results, and used to choose the container image on the TES path.
    tool: ToolPin | None = None

    # How much of a site's bundle its worker re-checksums against the DRS record before
    # computing.
    #
    #   off       trust the filesystem
    #   metadata  every file except the genotype .bed  (default; seconds)
    #   full      everything, .bed included            (minutes, at chromosome scale)
    #
    # `metadata` is the useful default because it catches the failures that actually
    # happen -- a stale bundle, a half-finished transfer, two runs' files mixed -- for a
    # cost nobody notices. `full` is for the run whose result gets published.
    verify_bundles: Literal["off", "metadata", "full"] = "metadata"


class Experiment(_Strict):
    """One experiment's federation-wide settings and its participating sites."""

    enabled: bool = True
    service_account: str
    endpoint_name: str
    sites: list[ExperimentSite] = Field(default_factory=list)

    # FLamby
    dataset: str | None = None
    num_clients: int | None = None
    rounds: int | None = None
    client_weights_mode: Literal["equal", "sample_size"] | None = None

    # GWAS and fine-mapping: which simulation scenario produced the distributed data.
    simulation_scenario: str | None = None

    # GA4GH: the study's data use request, the tool pin, and how much of each bundle a
    # site re-checksums. Absent means this experiment uses none of it and behaves exactly
    # as it did before -- see :class:`ExperimentGA4GH`.
    ga4gh: ExperimentGA4GH | None = None

    # GWAS
    variant_scaling: float | None = None
    hit_p_threshold: float | None = None

    # Fine-mapping.
    #
    # The sharding fields are federation-wide for a load-bearing reason: a site's uplink
    # is O(M^2) per locus, so a production run is split into several launches, and two
    # sites disagreeing about which split they are in would silently fine-map several
    # loci on a smaller cohort than the results claim. Declaring them here is what makes
    # "the same at every site" a property of the config rather than of the coordinator's
    # memory.
    locus_n_shards: int | None = None
    locus_shard_index: int | None = None
    locus_limit: int | None = None
    instance_limit: int | None = None
    pops: list[str] | None = None
    uplink_gram_dtype: Literal["float64", "float32"] | None = None

    # Coordinator-side. The ground-truth causal variants, for scoring a simulated run;
    # a real federation has no such file and leaves this unset.
    causal_manifest: str | None = None
    susiex_binary: str | None = None
    credible_set_level: float | None = None
    pval_thresh: float | None = None
    maf: float | None = None

    @model_validator(mode="after")
    def _unique_client_ids(self) -> Experiment:
        seen: dict[str, int] = {}
        for site in self.sites:
            if site.client_id in seen:
                raise ValueError(
                    f"duplicate client_id '{site.client_id}'. Each participating site "
                    "needs a distinct one -- APPFL uses it to key results, and two sites "
                    "sharing it means one silently overwrites the other."
                )
            seen[site.client_id] = 1

        centers = [s.center for s in self.sites if s.center is not None]
        if len(centers) != len(set(centers)):
            raise ValueError(
                "two sites are assigned the same dataset center. Each site must train "
                "on a distinct shard, or the federation is training twice on the same "
                "data and the result is not what it claims to be."
            )
        if self.num_clients is not None:
            over = [c for c in centers if c >= self.num_clients]
            if over:
                raise ValueError(
                    f"center(s) {over} are >= num_clients ({self.num_clients}). "
                    "The dataset loader asserts on this and will fail on the worker."
                )
        return self

    @model_validator(mode="after")
    def _locus_shard_is_coherent(self) -> Experiment:
        """An out-of-range shard index silently produces an empty run.

        The trainer would select no loci and return a payload with no genotype blocks;
        the aggregator would write an empty results table and report success. Refusing it
        here is the only place the mistake is visible.
        """
        shards = self.locus_n_shards
        index = self.locus_shard_index
        if shards is not None and shards < 1:
            raise ValueError(f"locus_n_shards must be at least 1, got {shards}")
        if index is not None:
            limit = shards if shards is not None else 1
            if not 0 <= index < limit:
                raise ValueError(
                    f"locus_shard_index {index} is out of range for locus_n_shards "
                    f"{limit}. Valid indices are 0 to {limit - 1}. A run with an "
                    "out-of-range index selects no loci and reports success on an empty "
                    "result table."
                )
        return self

    @property
    def site_count(self) -> int:
        return len(self.sites)


class Federation(_Strict):
    """A complete, validated, cross-referenced federation configuration."""

    schema_version: int = SCHEMA_VERSION
    coordinator: Coordinator
    sites: list[Site] = Field(default_factory=list)
    experiments: dict[str, Experiment] = Field(default_factory=dict)

    # Federation-wide GA4GH services. Which experiments use them is per-experiment; where
    # they live is not, because a DRS hostname or a TES URL that differed per experiment
    # would be two services described as one.
    ga4gh: GA4GHServices | None = None

    # Set by load_federation() so error messages can name the file.
    source_path: Path | None = None

    @field_validator("schema_version")
    @classmethod
    def _known_schema(cls, value: int) -> int:
        if value != SCHEMA_VERSION:
            raise ValueError(
                f"schema_version {value} is not supported (this build understands "
                f"{SCHEMA_VERSION}). Compare your file against federation.yaml.example."
            )
        return value

    @model_validator(mode="after")
    def _cross_reference(self) -> Federation:
        by_id = {s.id: s for s in self.sites}
        if len(by_id) != len(self.sites):
            # Counted over the declared ids, not over `by_id` -- that dict has already
            # collapsed the duplicates, so counting its keys always yields 1 and the
            # message named an empty list while correctly refusing the config.
            declared = [s.id for s in self.sites]
            dupes = sorted({sid for sid in declared if declared.count(sid) > 1})
            raise ValueError(f"duplicate site id(s): {dupes}")

        for exp_name, experiment in self.experiments.items():
            if exp_name not in EXPERIMENT_NAMES:
                raise ValueError(
                    f"unknown experiment '{exp_name}'. Known: {', '.join(EXPERIMENT_NAMES)}. "
                    "Adding a new one means adding it to the package, not just this file."
                )
            for entry in experiment.sites:
                site = by_id.get(entry.site)
                if site is None:
                    raise ValueError(
                        f"experiment '{exp_name}' references site '{entry.site}', which "
                        f"is not in the top-level `sites` list "
                        f"(have: {', '.join(sorted(by_id)) or 'none'})."
                    )
                object.__setattr__(entry, "_resolved_site", site)
        return self

    # -- lookups ------------------------------------------------------------

    def experiment(self, name: str) -> Experiment:
        """Get one experiment, with a message that lists the alternatives."""
        try:
            return self.experiments[name]
        except KeyError:
            known = ", ".join(sorted(self.experiments)) or "none"
            raise FederationError(
                f"experiment '{name}' is not declared in {self._where()}. Declared: {known}."
            ) from None

    def site(self, site_id: str) -> Site:
        for site in self.sites:
            if site.id == site_id:
                return site
        known = ", ".join(s.id for s in self.sites) or "none"
        raise FederationError(
            f"site '{site_id}' is not declared in {self._where()}. Declared: {known}."
        )

    def experiment_site(self, experiment: str, site_id: str) -> ExperimentSite:
        """Find a site's participation in an experiment.

        Accepts either the site id (`site-north`) or the client id (`Site1`), because
        both are natural things to type and getting told 'no such site' for the one you
        happened to pick is pure friction.
        """
        exp = self.experiment(experiment)
        for entry in exp.sites:
            if site_id in (entry.site, entry.client_id):
                return entry
        known = ", ".join(f"{e.site} ({e.client_id})" for e in exp.sites) or "none"
        raise FederationError(
            f"site '{site_id}' does not participate in experiment '{experiment}'. "
            f"Participating: {known}."
        )

    def enabled_experiments(self) -> dict[str, Experiment]:
        return {n: e for n, e in self.experiments.items() if e.enabled}

    # -- GA4GH ------------------------------------------------------------

    def drs_service(self) -> DrsService:
        """The DRS service, or raise naming what to add.

        Raises rather than returning None because every caller is already inside a code
        path that needs DRS -- resolving a site's ``drs_uri``, building a registry -- and
        an Optional here would only move the same message into five call sites.
        """
        service = (self.ga4gh.drs if self.ga4gh else None)
        if service is None:
            raise FederationError(
                f"no `ga4gh.drs` block in {self._where()}, but something asked for a DRS "
                "object. Add:\n"
                "    ga4gh:\n"
                "      drs:\n"
                "        hostname: drs.your-org.example\n"
                "        registry: local/data/<run>/drs_registry.json\n"
                "See docs/coordinator/ga4gh.md."
            )
        return service

    def tes_service(self) -> TesService:
        service = (self.ga4gh.tes if self.ga4gh else None)
        if service is None:
            raise FederationError(
                f"no `ga4gh.tes` block in {self._where()}, but a TES run was requested. "
                "Add `ga4gh.tes.url`, or run with --driver globus_compute."
            )
        return service

    def data_use_request(self, experiment: str) -> DataUseRequest | None:
        """This experiment's data use request, with the requester defaulted.

        The requester defaults to ``coordinator.identity`` -- the same string partners
        already authorize in their identity mapping. That is what makes DUO:0000026
        (user specific restriction) a check against a real, verified identity rather than
        a name somebody typed into a request form.
        """
        block = self.experiment(experiment).ga4gh
        if block is None or block.data_use_request is None:
            return None
        request = block.data_use_request
        if not request.requester:
            request = request.model_copy(update={"requester": self.coordinator.identity})
        if request.institution is None and self.coordinator.organization:
            request = request.model_copy(update={"institution": self.coordinator.organization})
        return request

    def tool_pin(self, experiment: str) -> ToolPin | None:
        block = self.experiment(experiment).ga4gh
        return block.tool if block else None

    def _where(self) -> str:
        return str(self.source_path) if self.source_path else "the federation config"


def _require_uuid(value: str, field: str) -> str:
    value = value.strip()
    if not _UUID_RE.match(value):
        hint = ""
        if value.startswith("<") or "REPLACE" in value.upper():
            hint = (
                " This still looks like a placeholder -- fill in the real value the "
                "partner sent you."
            )
        raise ValueError(
            f"{field}: '{value}' is not a UUID.{hint} A partner reads theirs from "
            "`globus-compute-endpoint list`, run as the same privileged user that "
            "started the endpoint (endpoint IDs are per-account)."
        )
    return value.lower()


def find_federation_file(explicit: str | Path | None = None) -> Path:
    """Locate the federation config, or raise with instructions for creating one."""
    if explicit is not None:
        path = Path(explicit).expanduser()
        if not path.is_file():
            raise FederationError(f"no federation config at {path}")
        return path

    env = os.environ.get(_ENV_VAR)
    if env:
        path = Path(env).expanduser()
        if not path.is_file():
            raise FederationError(f"{_ENV_VAR} points at {path}, which does not exist.")
        return path

    for candidate in _SEARCH_PATH:
        if candidate.is_file():
            return candidate

    searched = ", ".join(str(p) for p in _SEARCH_PATH)
    raise FederationError(
        "No federation config found.\n"
        "\n"
        "Create one from the worked example:\n"
        "    cp federation.yaml.example local/federation.yaml\n"
        "    $EDITOR local/federation.yaml\n"
        "\n"
        f"Searched: {searched}\n"
        f"You can also set {_ENV_VAR} or pass --federation.\n"
        "\n"
        "New here? Start with docs/coordinator/new-federation.md."
    )


def load_federation(path: str | Path | None = None) -> Federation:
    """Load and validate a federation config.

    Raises :class:`FederationError` with a message that names the fix -- this runs before
    anything touches the network, so it is the cheapest place to catch a mistake.
    """
    resolved = find_federation_file(path)

    try:
        raw: Any = yaml.safe_load(resolved.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise FederationError(f"{resolved} is not valid YAML:\n{exc}") from exc

    if raw is None:
        raise FederationError(f"{resolved} is empty. Start from federation.yaml.example.")
    if not isinstance(raw, dict):
        raise FederationError(
            f"{resolved} must be a YAML mapping at the top level, got {type(raw).__name__}."
        )

    try:
        federation = Federation.model_validate(raw)
    except Exception as exc:
        raise FederationError(f"{resolved} is not a valid federation config:\n\n{exc}") from exc

    object.__setattr__(federation, "source_path", resolved)
    return federation
