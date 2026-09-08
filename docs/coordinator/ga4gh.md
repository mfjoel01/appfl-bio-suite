# GA4GH standards in the fine-mapping experiment

Four GA4GH specifications are implemented, and each closes a hole this federation had.
They are not decoration on a working system; the fine-mapping run is assembled out of
them.

| | | answers |
| --- | --- | --- |
| **DUO** | Data Use Ontology | May this site's data be used for this study? |
| **DRS** | Data Repository Service | Exactly which bytes did a site compute over? |
| **TRS** | Tool Registry Service | Exactly which tool version computed them? |
| **TES** | Task Execution Service | Run that tool, there, on those bytes. |

Read in that order they compose into one sentence: *a named tool version, pinned by
checksum, runs over content-addressed data whose use terms permit the study that
dispatched it.* Every clause of that sentence used to be an assumption.

**All of it is opt-in per experiment.** A `federation.yaml` with no `ga4gh` blocks
produces byte-identical client configs to the ones it produced before any of this
existed, and `tests/test_ga4gh_enforcement.py` asserts that.

---

## What each one replaced

**DUO replaced nothing, which is the point.** A site's willingness to participate lived in
an email thread. Now it lives in `DATA_USE.json` *inside the site's own bundle*, the
site's worker refuses a study its terms do not permit, and the coordinator sees the same
refusal offline before dispatch. A site enforcing its own consent code, in its own
process, is the only version of this that means anything — the coordinator-side check is
a courtesy that saves a queue wait.

**DRS replaced "the path I told them to unpack it into".** A bundle is now a
content-addressed object with a sha-256 per file and a Merkle id over the set, so *did
site B run on the bundle I cut for site B* is answerable rather than assumed. The failure
it catches is real and silent: a site running last month's bundle produces well-formed
aggregates over the wrong individuals, pools without complaint, and is wrong in no visible
way.

**TRS replaced "whatever `pip install` resolved that morning".** APPFL ships a trainer by
reading `trainer.py` on the driver and sending its *source text* to a worker — a good
design, and the reason a partner installs two packages and no suite. But it means the code
a site runs is whatever the coordinator's checkout contained at dispatch time, and nothing
in the payload or the results said which version that was. Now the site stage is a
registered tool with a version, a descriptor checksum covering the shipped modules, and a
container image; the pin is verified at launch and recorded in the results.

**TES replaced nothing yet, and is the one genuinely optional addition.** Globus Compute
remains this federation's transport. TES is here because a partner who already operates a
TES service (Funnel, TESK, a cloud vendor's) can join without standing up an endpoint —
and because expressing the site stage as a TES task forced it to have a real command-line
form (`appfl-bio-suite site-stage`) instead of existing only as an APPFL trainer object.

---

## DUO — data use

### What a site declares

Every bundle carries a `DATA_USE.json`:

```json
{
  "dataset_id": "fine-mapping/covenant",
  "drs_uri": "drs://drs.example.org/620f6100…",
  "permission": "DUO:0000006",
  "modifiers": [{"id": "DUO:0000018"}],
  "steward": "Data Access Committee, Example Institute",
  "duo_version": "http://purl.obolibrary.org/obo/duo/releases/2021-02-23/duo.owl"
}
```

Exactly one **permission** (`DUO:0000001`'s five subclasses) and any number of
**modifiers** (`DUO:0000017`'s). `appfl-bio-suite ga4gh duo terms` lists what this build
understands.

The simulation writes one into every generated bundle, with terms from the scenario's
`data_use:` block. The three shipped sites declare *different* terms on purpose — a
scenario where every site says "general research use" exercises the consent machinery
about as well as a test suite with one passing assertion.

### What the study declares

In `federation.yaml`, once, for the whole experiment:

```yaml
experiments:
  fine-mapping:
    ga4gh:
      data_use_request:
        # requester defaults to coordinator.identity -- the same string partners
        # already authorize in their identity mapping, which is what makes
        # DUO:0000026 (user specific restriction) a check against a verified identity.
        project: cross-ancestry fine-mapping
        purposes: [DUO:0000038]        # genetic research
        non_commercial: true
        not_for_profit_organisation: true
        publication_agreed: true
        ethics_approval: IRB-2026-0114
```

Purposes are **research purpose** terms (`DUO:0000031`–`DUO:0000040`), not permissions.
Confusing the two is the mistake this vocabulary invites, so a permission in the
`purposes` list is refused by name.

One request for the whole experiment, deliberately: a coordinator who could vary the
declared purpose per site to get past a refusal would have a gate that gates nothing.

### The three outcomes

`permitted`, `denied`, and `undetermined`. The third is not a hedge — it is the honest
result for terms carrying free text a program cannot evaluate (`DUO:0000012`, *research
specific restrictions*, whose value might be "no use in studies of X"). It blocks dispatch
exactly as `denied` does, and is cleared only by a human recording that they read it:

```yaml
        acknowledged: [DUO:0000012]
```

Silently treating an unevaluable restriction as satisfied is the failure this design
exists to prevent. Acknowledging cannot clear a `denied` — only an `undetermined`.

### Where it is enforced

| | |
| --- | --- |
| `appfl-bio-suite ga4gh duo check` | offline, against your copies of the profiles |
| `appfl-bio-suite preflight --check ga4gh` | the same, with everything else |
| `appfl-bio-suite run …` | the launch gate; a refusal stops the run before dispatch |
| the site's worker | **the control.** Refuses before opening a genotype file |

Four graduated behaviours at the site, each deliberate:

| bundle has terms | run declares a study | what happens |
| --- | --- | --- |
| yes | yes | evaluated; a refusal raises `PermissionError` and reads nothing |
| yes | no | **refused.** A dataset with terms cannot be used by a study that states nothing about itself |
| no | yes | proceeds, recording that this site declared none |
| no | no | proceeds silently — a federation not using DUO is unaffected |

`enforce_data_use: false` waives the *coordinator's* gate for a dry run. It does not
disable the site's check — that one is not the coordinator's to disable — and a decision
that was overridden is recorded in the provenance, so a result computed under a waived
check never looks identical to one computed under a satisfied check.

### The ontology is vendored, not fetched

`src/appfl_bio_suite/core/ga4gh/data/duo.json` is a snapshot of the DUO release named in
its `version_iri`, extracted from the ontology's own OWL. Two reasons it is not resolved
over the network: a partner's compute node has no outbound web access, and a consent
decision that changes because an ontology release moved under a running federation is a
governance incident rather than a feature. Every decision records the release it was made
under.

---

## DRS — data objects

### Registering

`simulate` does it automatically (step 5), writing `drs_registry.json` beside the bundles.
For data you already distributed:

```bash
appfl-bio-suite ga4gh drs register --data-root local/data/fine-mapping-run \
    --hostname drs.your-org.example
```

Ids are **content-addressed**: a blob's id *is* its sha-256, and a bundle's is a Merkle
hash over its sorted `(name, child id)` pairs. Consequences worth knowing:

- Two coordinators registering byte-identical data mint identical ids.
- A bundle whose id changed is a bundle whose content changed — including a file *added*,
  which per-file checksums alone would never notice.
- `reference_variants.tsv` is hard-linked into every bundle and must be identical for the
  allele harmonization to be sound. It therefore appears as **one object with three
  `file://` access methods**, which is the correct picture of that property.

### Wiring it up

```yaml
ga4gh:
  drs:
    hostname: drs.your-org.example
    registry: local/data/fine-mapping-run/drs_registry.json
    https_base: https://drs.your-org.example   # optional; adds an https access method
    globus_collection: <collection uuid>       # optional; adds a globus access method

experiments:
  fine-mapping:
    ga4gh:
      verify_bundles: metadata
    sites:
      - site: site-north
        drs_uri: drs://drs.your-org.example/5ee8bab1…
```

`verify_bundles` decides how much a site re-checksums before computing:

| | | |
| --- | --- | --- |
| `off` | trust the filesystem | |
| `metadata` | everything except the genotype `.bed` | seconds — **the default** |
| `full` | everything | minutes at chromosome scale |

`metadata` is the useful default because it catches the failures that actually happen — a
stale bundle, a half-finished transfer, two runs' files mixed — for a cost nobody notices.
`full` is for the run whose result gets published.

The worker resolves nothing over the network: the coordinator inlines the object's
identity and per-member checksums into the client config, with the `file://` access
methods stripped (they point at the *coordinator's* filesystem).

### Serving

```bash
appfl-bio-suite ga4gh drs serve --port 8080     # read-only DRS 1.5
appfl-bio-suite ga4gh drs resolve drs://drs.your-org.example/5ee8bab1… --verify
```

**The served registry has no authorization.** Anything it serves is readable by anyone who
can reach the port. Bind it to localhost or an internal interface; a deployment holding
real genotypes puts its own gateway in front. There are also no writes — objects are
created by `simulate`, which knows what it produced, never by an HTTP request.

---

## TRS — the tool

```bash
appfl-bio-suite ga4gh trs publish --out local/trs
```

writes a servable TRS 2.0.1 tree, a CWL `CommandLineTool`, a `Containerfile`, a
`.dockstore.yml`, and `tool_pin.json`. The pin is the output that matters:

```yaml
experiments:
  fine-mapping:
    ga4gh:
      tool:
        id: "#workflow/github.com/your-org/appfl-bio-suite/fine-mapping-site-stage"
        version: "0.1.0"
        descriptor_checksum: 4d8fc94df1656448…
        image: registry.your-org.example/appfl-bio-suite:0.1.0
```

`descriptor_checksum` is a Merkle hash over the tool's files — **including the source of
`dataset.py` and `trainer.py`, the two modules APPFL actually ships to workers**. That is
what makes the pin cover the computation rather than a wrapper around it: editing a
comment in the trainer changes the checksum, and `preflight` fails until the pin is
refreshed.

`appfl-bio-suite ga4gh trs verify` checks the pin against the installed package on demand.

**No image is invented.** A freshly published tool version registers no container image
until you have built and pushed one, because a tool version naming an image nobody pushed
fails at the executor, minutes into a run, rather than at the pin:

```bash
podman build -t registry.your-org.example/appfl-bio-suite:0.1.0 -f local/trs/descriptors/Containerfile .
podman push  registry.your-org.example/appfl-bio-suite:0.1.0
appfl-bio-suite ga4gh trs publish --out local/trs \
    --image registry.your-org.example/appfl-bio-suite:0.1.0 \
    --image-digest <the sha256 the push reported>
```

Pin the digest, not the tag. A tag can be repushed; a federation that agreed on a tag has
agreed on nothing durable.

### Dockstore

The generated `.dockstore.yml` registers the tool with [Dockstore](https://dockstore.org),
itself a TRS implementation, once its GitHub App is installed on the repository. That is
what makes the pin in your `federation.yaml` resolvable by someone who is not in your
federation — the only reason a tool registry beats a git tag.

---

## TES — execution

Globus Compute stays the default. TES is a second transport, for a partner who already
runs one.

```bash
appfl-bio-suite ga4gh tes task --experiment fine-mapping      # write the task documents
appfl-bio-suite run fine-mapping --driver tes                 # submit and aggregate
```

`ga4gh tes task` submits nothing. It writes the complete task document — image, command,
inputs by DRS URI, resource request, and the provenance tags — which is what to send a
partner who asks what will actually run on their cluster.

```yaml
ga4gh:
  tes:
    url: https://tes.partner.example
    poll_seconds: 15
    timeout_seconds: 7200
    cpu_cores: 4
    ram_gb: 16
    disk_gb: 64

experiments:
  fine-mapping:
    sites:
      - site: site-north
        tes_url: https://tes.north.example       # optional per-site override
        tes_outputs_url: file:///shared/fm/anl   # where its outputs land
```

Two requirements the Globus path does not have:

1. **A published container image**, named by the TRS pin. Without it the driver refuses to
   submit and says so.
2. **Outputs the driver can read.** TES stages a task's outputs to a URL the *service* can
   write; this driver then reads them back to aggregate. A `file://` URL must resolve on
   the driver's machine too — a shared filesystem, or an object store both ends reach.

Neither is a limitation of this implementation; they are what running someone else's
container on someone else's cluster costs.

The TES driver does not run APPFL's round loop, and that is a property of the experiment
rather than of TES: fine-mapping is a **single-round aggregate** exchange, so there is no
global model to distribute and no second round. It submits one task per site, waits, reads
the payloads back, and calls the same `FineMappingAggregator` the other two drivers call.

A site whose DUO terms refuse the study exits **77**, so a task log distinguishes "this
site declined" from "this task crashed" — completely different responses from a
coordinator.

---

## What a run records

Every GA4GH-configured run writes `ga4gh_provenance.json` beside its results:

```
dispatched   what this coordinator sent: the data use request, the tool pin, the DRS
             objects it expected each site to hold, and the version of all four specs
attested     what each site says it did: its DUO decision term by term, the DRS object
             it verified against, the tool pin it was given
outputs      a DRS object per results file, in drs_outputs.json
```

The two halves are compared on arrival. A site that computed over a **different DRS
object** than the run's provenance names is a hard error — every number it contributed
would be attributed to data it did not read, and no downstream reader could tell. A
`denied` decision that arrived anyway (enforcement waived) is logged loudly and recorded.

The coordinator can only record what it dispatched; the site records what it actually did.
That is why the attestation travels in the payload rather than being reconstructed here.

---

## Adopting this on an existing federation

In order, each step useful on its own:

1. **`ga4gh trs publish`**, pin the tool. Costs nothing, and from then on every result
   names the version that produced it.
2. **`ga4gh drs register`**, add each site's `drs_uri`. Sites start verifying their own
   bundles; a stale one now fails loudly.
3. **Ask each partner for their `DATA_USE.json`**, record the path in `data_use_profile`,
   and declare your `data_use_request`. Run `ga4gh duo check` before you launch.
4. **TES only if a partner asks for it.** It needs an image you publish and an outputs
   path you can read; the Globus Compute path needs neither.

Steps 1–3 change no partner's setup. Step 4 replaces it.

---

## What is deliberately not implemented

- **No DRS writes, and no authorization on the served registry.** DRS 1.5 is a read API
  and this is a read implementation.
- **No TRS server with a database.** A static tree at the spec's paths, servable by any
  file host, plus Dockstore for anything beyond that.
- **No TES *service*.** This suite is a TES client. Partners run Funnel, TESK, or a cloud
  service.
- **No compact-identifier DRS URIs** (`drs://prefix:accession`). Resolving them needs
  identifiers.org, and a partner's worker has no outbound network.
- **No GA4GH Passports/AAI.** Authorization here is Globus identity mapping, which the
  federation already uses and which partners already control. DUO decides what a study may
  do; it is not an authentication system, and wiring one in would mean asking every
  partner to trust a second identity provider.
