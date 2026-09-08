# The federation network map

A page that shows every institution in the federation on a map: where they are, which
experiments they run, what kind of data they hold, how many samples they declared, and —
if you ask for it — whether their endpoint is answering right now.

Built on **[hivewatch](https://github.com/APPFL/hivewatch)**, APPFL's monitoring toolkit.
hivewatch supplies the event schema, the map metadata format, the live server and the
viewer. This suite supplies the one thing hivewatch cannot know: who is in *your*
federation and what each of them is doing. That comes from `federation.yaml`, the same
file everything else here reads from.

```
federation.yaml ──► appfl-bio-suite watch ──► hivewatch artifacts ──► the map
                                                    │
  run --watch ────────────────────────────────────► │  (live rounds, same page)
```

---

## The map exists before any run does

This is the part worth stating plainly, because it is the opposite of how monitoring
usually works. A dashboard normally has nothing to show until something runs. A
federation, though, is a real thing well before its first launch: partners have been
recruited, endpoints are up, bundles are unpacked, and there is a genuine question —
*who is actually in this thing?* — that nobody has a page for.

So the network view is built from `federation.yaml` alone. The day the config names a
partner, that partner is on the map.

```bash
pip install 'appfl-bio-suite[watch]' -c constraints.txt

appfl-bio-suite watch build          # draw it from federation.yaml
appfl-bio-suite watch serve          # open http://localhost:7070
```

`watch build` writes the network view into `local/watch/runs/` as an ordinary hivewatch
run under a fixed id, which is what lets it sit beside your real runs and appear
alongside them in the viewer's run list. Rerun it whenever `federation.yaml` changes; it
replaces rather than accumulates.

---

## Giving everyone the same page

`watch serve` binds a port on the machine you run it on. On an HPC login node that means
you can see it and nobody else can — partners cannot reach inbound ports on your cluster,
and the cluster is right not to let them.

What you hand out is a static export:

```bash
appfl-bio-suite watch export --out site/
```

Three files — the viewer, the data, and an index — with no server of ours behind them.
Put the directory on GitHub Pages, an institutional web host, or any object store, and
every coordinator and partner in the federation is looking at the same page.

```bash
python -m http.server -d site/ 8000     # preview it first
```

A plain `file://` open will not work: the page fetches its data from alongside itself, and
browsers refuse that for local files. Any static server does.

### What ends up in the published file

By design, because it is the map:

- institution names, countries, cities and coordinates
- which experiments each site runs
- what kind of data each holds, and the sample count they declared
- endpoint **fingerprints** (`cccccccc…0001`) — enough to check a site against a config
  you already hold, not enough to be an identifier

Never:

- service accounts, `data_dir` and `output_dir` paths, or anything else about the inside
  of a partner's cluster
- your Globus identity
- full endpoint UUIDs, unless you pass `--include-endpoint-uuids`

That last one is a default rather than a rule. A partner's endpoint refuses tasks from
anyone their identity mapping does not name, so a leaked UUID grants nothing — but
"grants nothing" is a claim about *their* configuration, not yours, and this is a page
meant to be handed out. For an internal deployment, opt in.

Read `site/network.map.json` before you publish it. It is small and it is the whole
payload.

---

## Coordinates are declared, not detected

Add a `location` block to each site in `federation.yaml`:

```yaml
sites:
  - id: site-north
    name: Northern Institute of Technology
    country: Canada
    location:
      lat: 43.6532
      lng: -79.3832
      city: Toronto
```

The coordinator gets one too, under `coordinator:`. It becomes the hub every partner
marker connects back to.

hivewatch's own examples resolve each client's location by calling `ipinfo.io` from the
client. That is a good default for a laptop demo and wrong here in three separate ways: a
partner's compute node has no outbound web access, so the call fails; if it succeeded it
would return the institution's border router rather than the cluster; and it puts an
outbound call inside a partner's security boundary for the sole purpose of drawing a dot,
which is a conversation with their security office that no dot is worth.

You already know where your partners are. Writing it down once is cheaper and more
accurate.

**There is no country-centroid fallback.** A site with no `location` is reported as
unplaced, with the YAML to paste, rather than drawn at a guessed position — putting a
partner's name on a map at coordinates they never gave is worse than a gap. Nothing about
a run depends on this: an unplaced site trains exactly as it always did.

---

## Live runs on the same map

```bash
appfl-bio-suite run gwas --watch
```

Per-round, per-site updates stream into the same runs directory, so `watch serve` shows
the standing federation and what it is doing right now on one page. The directory is
watched, so a run launched while the server is up appears without restarting it.

Two properties are load-bearing:

**A site's geography never travels with its results.** The worker returns statistics.
Where that site is comes from `federation.yaml` on your machine, and the two are joined
on the coordinator. A partner's payload gains nothing it did not already carry.

**Watching cannot break a run.** `--watch` with hivewatch missing prints a note and
launches anyway. A monitoring error mid-run disables the map, logs once, and the run
continues. A federated run costs a scheduler queue wait at every participating
institution; losing one because a dashboard raised would be an absurd trade.

Only an allowlist of numeric metrics is forwarded — accuracy, loss, sample count,
gradient norm, training time, bytes sent. Everything else a trainer returns is a *result*,
not telemetry, and results do not belong on a page that might be published.

---

## Endpoint liveness

```bash
appfl-bio-suite watch build --probe
```

Asks Globus Compute about every endpoint and colours the map by the answer. Costs one API
call per endpoint. Without it the map states **membership**, not liveness — which is the
honest reading when nothing has been asked.

A site with two endpoints and one of them down is neither up nor down; it is drawn as
`idle`, distinct from both, because a site half-online is exactly the state worth being
able to see at a glance.

---

## Command reference

```
appfl-bio-suite watch build   [--experiment X] [--probe] [--include-endpoint-uuids]
appfl-bio-suite watch serve   [--port 7070] [--runs-dir DIR]
appfl-bio-suite watch export  --out DIR [--experiment X] [--probe] [--title "..."]
appfl-bio-suite run EXPERIMENT --watch
```

`--experiment` narrows every command to one experiment's sites. The default is every
enabled experiment, which is the point: one map, every experiment, one federation.

---

## Troubleshooting

**`hivewatch is not installed`** — it is a coordinator-only extra and never reaches a
partner: `pip install 'appfl-bio-suite[watch]' -c constraints.txt`.

**The map is empty and the sidebar says "Waiting for clients…"** — no participating site
has coordinates. `watch build` names them and prints the YAML.

**A site is in the sidebar but not on the map** — that one site has no `location`. The
viewer lists a client it cannot place; it does not invent a position.

**The exported page is blank** — it was opened as `file://`. Serve the directory.

**The map draws but there are no base tiles** — the viewer loads Leaflet and its tiles
from the public internet. That is a property of hivewatch's viewer, not of the export;
a browser with no web access renders the markers on an empty background.
