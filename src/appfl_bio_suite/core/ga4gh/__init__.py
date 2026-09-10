"""GA4GH standards, as the connective tissue of a federated fine-mapping run.

Four GA4GH specifications are implemented here, and each one closes a hole this suite
had. They are not decoration on a working system; the run is assembled out of them.

    DUO   Data Use Ontology            May this site's data be used for this study?
    DRS   Data Repository Service      Exactly which bytes did this site compute over?
    TRS   Tool Registry Service        Exactly which tool version computed them?
    TES   Task Execution Service       Run that tool, there, on those bytes.

Read in that order they compose into one sentence: *a named tool version, pinned by
checksum, runs over content-addressed data whose use terms permit the study that
dispatched it.* Every clause of that sentence used to be an assumption.

WHAT EACH ONE REPLACED
----------------------
**DUO** replaced nothing, which is the point. A site's willingness to participate lived
in an email thread. Now it lives in ``DATA_USE.json`` inside the site's own bundle, the
site's worker refuses a study its terms do not permit (``dataset.py``), and the
coordinator sees the same refusal offline, before dispatch (``preflight --check ga4gh``).
A site enforcing its own consent code, in its own process, is the only version of this
that means anything -- the coordinator-side check is a courtesy that saves a queue wait.

**DRS** replaced "the path I told them to unpack it into". A bundle is now a content-
addressed object with a sha-256 per file and a Merkle id over the set, so "did site B run
on the bundle I cut for site B" is answerable rather than assumed. The failure it catches
is real and silent: two simulation runs differing in seed produce bundles that validate,
load, and pool into results that are wrong in no visible way.

**TRS** replaced "whatever ``pip install`` resolved that morning". The site stage is a
registered tool with a version and a descriptor checksum; the pin is verified at launch
and recorded in the results. This is the standard the other three lean on -- a TES task
without a TRS pin runs *something*, and a DRS object without a tool that produced it is
provenance with a hole in the middle.

**TES** replaced nothing yet, and is the one addition that is genuinely optional. Globus
Compute remains the transport this federation runs on. TES is here because a partner who
already operates a TES service (Funnel, TESK, or a cloud vendor's) can join without
standing up an endpoint at all, and because expressing the site stage as a TES task is
what forced it to have a real command-line form (``site_stage.py``) instead of existing
only as an APPFL trainer object.

WHERE THE SEAMS ARE
-------------------
    duo.py         terms, profiles, requests, and the matcher. Pure functions.
    drs.py         DRS 1.5 objects, a registry built from a simulation run, a client,
                   and a conformant read-only server.
    trs.py         TRS 2.0.1 tool + version documents, descriptor generation
                   (CWL and Dockstore), a client, and pin verification.
    tes.py         TES 1.1 task documents and a client.
    data/duo.json  the vendored DUO release. Not fetched at run time: a partner's
                   worker has no outbound network, and an ontology that changes under a
                   running federation is a consent decision changing silently.

None of it adds a dependency. The models are pydantic (already required) and the clients
are ``urllib`` (standard library), because every one of these has to be usable from a
worker that has the experiment's two-package extra and nothing else.

THE ONE RULE THAT CONSTRAINS EVERYTHING HERE
--------------------------------------------
``experiments/fine_mapping/dataset.py`` is shipped to partner workers and may not import
this package -- so the DUO matcher exists twice, here and inlined there, exactly as the
site-aggregate computation does. ``tests/test_ga4gh_shipped_parity.py`` runs both over
the same case table and asserts identical decisions. Duplication that is tested is a
cost; duplication that is trusted is a bug waiting for a partner to find.
"""

from __future__ import annotations

__all__ = ["drs", "duo", "tes", "trs"]
