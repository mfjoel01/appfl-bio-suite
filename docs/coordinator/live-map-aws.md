# The live map on AWS

The project site on GitHub Pages draws the federation from committed files. It can say
who the partners are; it cannot say whose endpoint is answering, because nothing behind a
static page is able to ask. This is the same site on one small EC2 instance, with a
prober beside it that asks Globus Compute every minute and a map that shows the answer.

```
python scripts/aws/live_map.py deploy      # create everything, print the URL
python scripts/aws/live_map.py status      # is it up, and what does it report
python scripts/aws/live_map.py update      # serve a newer commit, new targets or a new secret
python scripts/aws/live_map.py destroy --yes
```

Nothing is done in the AWS console. Every resource is created through the API and tagged
`Project=appfl-bio-suite-live-map`, which is also how `destroy` finds it again.

## What "live" means here

Partner clusters accept no inbound connections, so there is nothing to ping. What can be
asked, from anywhere, is whether an endpoint's daemon is connected to the Globus Compute
service. The map shows that, per institution:

| On the map | Meaning |
|---|---|
| Endpoint online | every endpoint the site runs is connected |
| *n* of *m* endpoints online | some are |
| Endpoint offline | none is |
| Endpoint status unknown | the prober could not ask, or its last answer is too old |

It is the daemon's status and not the round trip `appfl-bio-suite endpoint smoke` makes:
an endpoint can be online and still refuse or fail a task. `core/endpoint.py` lists what
can be wrong behind an `online`.

**Unknown is not offline.** Globus Compute does not let every identity read every
endpoint's status. Partners' multi-user endpoints answer an identity that does not own
them. A single-user endpoint owned by a site's service account answers only that account
-- MBZUAI's refused the coordinator's own login with a 403 while it was running. A site
the prober may not ask about is drawn as unknown.

To show such a site, give its endpoint a **route** in the probe list: the client id of
the account that owns it, and the name of the environment variable holding a secret for
that client. That account can submit tasks to the endpoint, not only read its status, so
mint a separate secret on the client for the prober. It is then revocable on its own, and
the coordinator's secret never leaves Polaris. The same holds if the prober's
credential expires -- the page says it does not know, rather than that every partner is
down -- and if the prober stops altogether, the viewer retires its last answer after
three intervals and the header reads "Status stale".

Sites with no endpoint yet show their partnership stage and no status line, as before.

## How it is put together

```
browser ── 80/443 ──> Caddy ── static files ── /opt/live-map/site      (scripts/build_site.sh)
                                                   └─ map/endpoints.status.json
                                                            ▲ rewritten every minute
                      live-map-prober  (appfl-bio-suite watch live) ──> Globus Compute API
                            ▲ at start: targets + Globus secret
                      SSM Parameter Store  (read with the instance role)
```

* **The site** is the output of `scripts/build_site.sh`, from `website/federation.yaml`
  and `website/partners.json` -- the same inputs, and the same `check_site_preview.py`
  gate, as GitHub Pages. `SITE_LIVE_STATUS=1` adds one field to the map's metadata telling
  the viewer to poll `endpoints.status.json`.
* **The probe targets** are the one thing the instance needs from your real federation:
  public site id to endpoint UUIDs. `local/federation.yaml`, with its contacts and paths,
  never leaves your machine.
* **The status file** holds a state and two counts per site. No UUID, no fingerprint;
  `tests/test_live_status.py` asserts that against the bytes.
* **The secrets** -- the targets and the prober's Globus service account -- go from your
  machine to Parameter Store (encrypted) to the prober's memory. They are not in
  user-data, which is readable through the EC2 API.
* **No SSH.** The security group opens 80 and 443. `update` runs through Systems Manager.
* **HTTPS** needs a name. With `--domain map.example.org`, point an A record at the
  address `deploy` prints. Without one the deployment uses `<address>.sslip.io`, a public
  DNS service that resolves such names to the address in them; `--no-tls` serves plain
  HTTP on the address instead.

## Before the first deploy

1. **AWS credentials** in the environment or `~/.aws` (`--profile` selects one). The
   identity needs EC2, IAM role and instance-profile, and SSM Parameter Store permissions.

2. **A Globus service account for the prober.** At app.globus.org: Settings, Developers,
   Register a service account. Save it, mode 0600, as `local/live-map/globus-client.json`:

   ```json
   {"client_id": "…", "client_secret": "…"}
   ```

   Give it nothing else. If a partner's endpoint shows as unknown under it, that
   endpoint does not answer this identity's status requests.

3. **The probe targets**, from the real federation:

   ```bash
   appfl-bio-suite watch targets --out local/live-map/targets.json
   ```

   Read the file. Site ids must be the ones the *published* map uses
   (`website/partners.json`); `--rename LOCAL=PUBLIC` maps one that differs, and `deploy`
   warns about any id the map does not draw. Placeholder UUIDs are skipped.

   `--endpoint SITE=UUID` adds an endpoint the federation file does not name, such as one
   recorded only at registration. `--site-credentials MANIFEST` copies the route for each
   probed endpoint your coordinator manifest routes to a service account -- the client id
   and the variable *name*, never the secret. Before `deploy`, export that variable with
   the prober's own secret for the client; it is stored encrypted beside the targets.

4. **Push the commit.** The instance clones the public repository at the commit you
   deploy, so it has to be on GitHub.

`deploy --dry-run` checks all of this and prints the first-boot script without calling AWS.

## Running it without AWS

The prober and the export are ordinary CLI commands and do not care where they run:

```bash
SITE_LIVE_STATUS=1 scripts/build_site.sh local/site-preview
appfl-bio-suite watch live --targets local/live-map/targets.json --out local/site-preview/map &
python -m http.server -d local/site-preview 8000
```

## Cost and teardown

One `t3.small`, a 20 GB volume and one public IPv4 address: under 20 USD a month at
on-demand prices in us-east-1 at the time of writing. `destroy --yes` removes the
instance, the address, the security group, the role and the stored parameters.
