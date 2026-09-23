# The federation network map

A page that shows every institution in the federation on a map: where they are, which
experiments they run, what kind of data they hold, how many samples they declared, and —
if you ask for it — whether their endpoint is answering right now.

Built on **[hivewatch](https://github.com/APPFL/hivewatch)**, APPFL's monitoring toolkit.
hivewatch supplies the event schema, the map metadata format, the live server and the
base viewer. This suite adds a federation viewer with an interactive globe and experiment
filters, and supplies who is in *your* federation and what each of them is doing. That
comes from `federation.yaml`, the same file everything else here reads from.

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

### Exploring the federation

The viewer opens in the familiar **Flat map** view. Switch to **Globe** for a rotating globe:
drag to turn it, zoom to look closer, and select a site to inspect its details. You can
pause rotation. Arrow keys rotate the focused globe, `+` and `-` zoom, `Space` toggles
rotation, and `Home` resets the view. Rotation starts paused when your browser requests
reduced motion. Both views use the same sites and experiment selection.
The theme button also changes the globe: light mode uses pale oceans, mint land,
and a bright sky, while dark mode retains the night palette.

Two query parameters let an embedding page open the viewer already matching itself
rather than flashing the default and then being corrected: `?view=globe` starts on the
globe instead of the flat map, and `?theme=light` starts in the light palette. Neither is
remembered -- the viewer's own controls own the state from there on -- and any other
value keeps the default. The project site's embed uses both.

The header carries the suite logo and, on the right, a **BUILT ON** credit linking to
[hivewatch](https://github.com/APPFL/hivewatch) and [APPFL](https://github.com/APPFL/APPFL).
Upstream's inlined wordmark is dropped from the export: it was most of the viewer's
bytes, and republishing someone else's mark on every host an export lands on is a claim
about their trademark terms rather than a credit. A named link is the credit that
travels.

The **Experiments** tab starts with **All experiments** selected, showing every site in
the loaded federation. Select **Fine-mapping**, **GWAS**, or any combination of experiments
to show sites participating in at least one selection. A site participating in several
selected experiments still appears once. These controls filter the page locally; they
do not change the federation configuration or launch experiments. Options are drawn from
the loaded data: a federation with only fine-mapping configured shows only fine-mapping.
Clear every checkbox to hide all sites; **All experiments** restores them. Selecting a
different run resets the experiment filter to All; changing the map view preserves it.

Run monitoring remains available under **Runs**, including selecting recorded runs and
playback controls. The federation view takes priority, and round, global accuracy and
global loss no longer occupy the page header. `watch serve` and `watch export` use the
same viewer.

### Partners and experiment results

The **Results** tab sits between Experiments and Runs. Choose an experiment to see
its published figures, tables, and reports in the main workspace. Figures expand;
CSV/TSV tables support search and numeric or text sorting; HTML reports open in a
sandboxed frame. Every artifact has a download link. Experiments without results
show an empty state. Map filters remain independent and are preserved when returning
to Experiments; a single selected map experiment becomes the initial Results selection.

An optional JSON catalogue adds planning partners and explicitly selected result files:

```bash
mkdir -p local/watch
cp watch.catalog.json.example local/watch/catalog.json
# Edit the catalogue's partners and result paths, then:
appfl-bio-suite watch build --catalog local/watch/catalog.json
appfl-bio-suite watch serve
# Or create the same static viewer:
appfl-bio-suite watch export --catalog local/watch/catalog.json --out local/watch-site
```

Paths are resolved relative to the catalogue file, and absolute paths also work.
For example, an artifact with `"path": "results/summary.tsv"` in that catalogue
reads `local/watch/results/summary.tsv`. The catalogue is loaded only
when `--catalog` is supplied. Rebuild or re-export after changing it.

Each partner has a stable `id`, `name`, `country`, `projects` array and `stage`.
Optional fields are `contacts`, `notes`, `location`, `location_basis`, and `source_url`.
Use the existing federation site id to enrich that site without duplicating it;
the federation's sample counts, coordinates and endpoint status stay authoritative.
New partners have no declared sample count or compute endpoint. A planning stage
does not imply that an experiment is running. Use `location: null` for an unresolved
location, or document a representative institutional marker with its source.

Each results group names an `experiment`, `title`, optional `description`, and an
`artifacts` array of `{ "title": "...", "path": "..." }` objects. Artifact descriptions
are optional. Supported files are PNG, JPEG, WebP, CSV, TSV and HTML, up to 16 MiB each.
Tables preview the first 200 rows; search and sorting apply to that preview, while
downloads contain the complete file. HTML report scripts are disabled in the viewer.
For interactive report content, publish a static figure or table alongside the report.

The catalogue publishes the listed contacts, notes, and entire selected file contents.
Choose files intended for the viewer's audience. The catalogue's filesystem paths
are not included in the generated metadata, and no result directories are scanned.
Artifacts are embedded in `network.map.json`, so the export still contains three files.
Use `local/` for your deployment's roster and catalogue; the repository example contains
only a fictional partner and an empty results group ready for your selected files.

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

This repository does exactly that for a fictional federation. `scripts/build_site.sh`
exports `website/demo-federation.yaml` -- twenty invented institutions on six continents
-- and `.github/workflows/pages.yml` runs that script on every push to `main`, publishing
the result at <https://mfjoel01.github.io/appfl-bio-suite/network.html>. Run the same
script to preview locally.

The fictional configs are the only ones the build is ever given:
`scripts/check_site_preview.py` fails it if the published metadata was assembled from
anything else, or carries a contact block or an email address. A partner's address
published to a public URL is not undone by deleting the page afterwards.

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

Read `site/network.map.json` before you publish it. It is the whole data payload;
with a catalogue, it also contains embedded results and the partner details you selected.

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

`--experiment` on `watch build` and `watch export` limits the generated data to one
experiment's sites. The default includes every enabled experiment. The viewer's
checkboxes can filter the experiments present in that data; they cannot restore sites
excluded when the data was generated.

---

## Troubleshooting

**`hivewatch is not installed`** — it is a coordinator-only extra and never reaches a
partner: `pip install 'appfl-bio-suite[watch]' -c constraints.txt`.

**The map is empty and the sidebar says "Waiting for clients…"** — no participating site
has coordinates. `watch build` names them and prints the YAML.

**A site is in the sidebar but not on the map** — that one site has no `location`. The
viewer lists a client it cannot place; it does not invent a position.

**The exported page is blank** — it was opened as `file://`. Serve the directory.

**The map draws but there are no base tiles** — the flat map loads Leaflet and its tiles
from the public internet. The globe's geography and rendering code are embedded in the
viewer, but the base viewer still requires Leaflet to initialize. Serve the page in a
browser with web access.

**`The installed hivewatch viewer has an unsupported layout`** — this suite extends the
viewer shipped in the pinned `hivewatch==0.2.1`. Reinstall the `[watch]` extra using
`constraints.txt`. The viewer is assembled into a separate file; the installed hivewatch
package is never edited.
