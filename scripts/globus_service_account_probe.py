#!/usr/bin/env python
"""Exercise the proposed "shared service account" endpoint-sharing design, for real.

WHAT IS BEING PROPOSED
----------------------
A federation stops using multi-user endpoints (which need root on every partner's login
node) and instead does this:

1. APPFL registers a Globus Auth *service account* -- a `client_identity` client -- per
   federation.
2. Every partner exports that client's ID and secret and runs an ordinary single-user
   endpoint.
3. Because all endpoints are then *owned* by the one client identity, anything holding
   those credentials can submit to any of them.

Globus have confirmed the mechanism works and suggested minting a separate secret per
partner so any one partner can be revoked individually.

WHY A PROBE RATHER THAN A DISCUSSION
------------------------------------
Three of the load-bearing claims are checkable in about a minute each, and the answers
decide whether the design is acceptable:

* Does a per-partner secret actually scope anything?  (`endpoints --as`)
* Does deleting one secret revoke one partner without disturbing the others? (`revoke`)
* What can a partner holding the federation secret see and do to *other* partners'
  endpoints?  (`endpoints --as`, again -- this is the whole question)

`--shape` selects which variant to build, so the two can be measured side by side rather
than argued about:

  shared    one client identity for the whole federation, N secrets  (as proposed)
  per-site  one client identity per site, one secret each            (the alternative)

WHAT THIS CREATES IN GLOBUS
---------------------------
`provision` creates a real Auth project and real client identities under the running
user's Globus account, and mints real secrets. `teardown` deletes all of it. Nothing here
touches the production federation: project and client names carry the `--federation`
prefix you pass, and `teardown` refuses to delete anything it did not create.

Secrets are written to `local/` -- which is gitignored in its entirety -- at mode 0600,
and are never printed. Commands that need to prove which credential is in play print a
fingerprint instead.

USAGE
-----
    python scripts/globus_service_account_probe.py inspect
    python scripts/globus_service_account_probe.py provision --federation probe \
        --sites alpha,beta --shape shared
    python scripts/globus_service_account_probe.py whoami --as alpha
    python scripts/globus_service_account_probe.py endpoints --as alpha
    python scripts/globus_service_account_probe.py revoke --as beta
    python scripts/globus_service_account_probe.py teardown --federation probe

The endpoint half of the experiment is deliberately not automated; see `env --as` for the
two exports a partner would really set, and REGISTERING AN ENDPOINT below.

REGISTERING AN ENDPOINT UNDER A PROBE IDENTITY
----------------------------------------------
    eval "$(python scripts/globus_service_account_probe.py env --as alpha)"
    globus-compute-endpoint configure probe-alpha
    globus-compute-endpoint start probe-alpha

No browser login happens: with both variables set the SDK builds a `ClientApp` and gets
tokens by client-credentials grant. Tokens land in `~/.globus_compute/storage.db` under
the namespace `clientprofile/production/<client_id>`, *beside* rather than on top of any
personal login already there -- so this is safe to run on a machine with a live endpoint.
"""

from __future__ import annotations

import argparse
import json
import os
import stat
import sys
from pathlib import Path

# Where minted secrets land. `local/*` is gitignored, so this cannot be committed by
# accident; see .gitignore and local/README.md.
STATE_DIR = Path(__file__).resolve().parent.parent / "local" / "globus_service_account"


def _state_file(federation: str) -> Path:
    """One file per federation, so the two `--shape`s can exist side by side.

    The whole point of this probe is to compare them, which means both sets of
    credentials have to be live at once.
    """
    return STATE_DIR / f"{federation}.json"


def _auth_client():
    """An AuthClient on the *human's* login, which is what owns projects and clients.

    Deliberately not affected by GLOBUS_COMPUTE_CLIENT_ID/SECRET being exported in the
    caller's shell: a client identity cannot create or revoke its own credentials, so if
    those are set we would fail later with a confusing 403 rather than here with a
    sentence. The unsetenv is scoped to this process.
    """
    os.environ.pop("GLOBUS_COMPUTE_CLIENT_ID", None)
    os.environ.pop("GLOBUS_COMPUTE_CLIENT_SECRET", None)

    from globus_compute_sdk.sdk.auth.auth_client import ComputeAuthClient
    from globus_compute_sdk.sdk.auth.globus_app import get_globus_app

    app = get_globus_app()
    if app.login_required():
        sys.exit("Not logged in. Run `globus-compute-endpoint login` first.")
    return ComputeAuthClient(app=app)


def _explain_session_error(e: Exception) -> None:
    """Turn Auth's session-policy 403 into the sentence that tells you what to do.

    This is not a rare edge. Every project-administration call -- creating a client,
    minting a credential -- is gated on having authenticated *interactively* within the
    last 30 minutes, and a stored refresh token does not satisfy it. It is the single
    most important constraint on automating partner onboarding, so it gets a real message
    rather than a traceback.
    """
    from globus_sdk import AuthAPIError

    if not isinstance(e, AuthAPIError) or e.http_status != 403:
        return
    if "session" not in str(e).lower():
        return
    sys.exit(
        "\nGlobus Auth refused: project administration requires an identity that "
        "authenticated\ninteractively within the last 30 minutes. A stored refresh token "
        "does not count.\n\n"
        "  globus-compute-endpoint login --force\n\n"
        "then re-run this command. Note this for the design discussion: minting "
        "per-partner\nsecrets cannot be fully unattended on the project owner's stored "
        "credentials alone."
    )


def _fingerprint(secret: str) -> str:
    """Enough to tell two secrets apart in a transcript, not enough to use one."""
    return f"{secret[:4]}...{secret[-2:]} (len {len(secret)})"


def _load_state(federation: str) -> dict:
    path = _state_file(federation)
    if not path.exists():
        sys.exit(f"No probe state at {path}. Run `provision --federation {federation}` first.")
    return json.loads(path.read_text())


def _save_state(state: dict) -> None:
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    path = _state_file(state["federation"])
    path.write_text(json.dumps(state, indent=2))
    path.chmod(stat.S_IRUSR | stat.S_IWUSR)


def _site_creds(state: dict, site: str) -> dict:
    if site not in state["sites"]:
        sys.exit(f"Unknown site {site!r}. Provisioned: {', '.join(state['sites'])}")
    return state["sites"][site]


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------


def cmd_inspect(args: argparse.Namespace) -> int:
    """Read-only. What does this Globus identity already own?"""
    ac = _auth_client()
    me = ac.userinfo()
    print(f"identity : {me.get('preferred_username') or me.get('sub')}")
    print(f"id       : {me.get('sub')}")

    projects = list(ac.get_projects())
    print(f"\nprojects : {len(projects)}")
    for p in projects:
        print(f"  {p.get('display_name')!r}  id={p.get('id')}")

    clients = list(ac.get_clients())
    print(f"\nclients  : {len(clients)}")
    for c in clients:
        creds = list(ac.get_client_credentials(c["id"]))
        print(f"  {c.get('name')!r}  id={c['id']}  type={c.get('client_type')}")
        for cred in creds:
            print(f"      credential {cred.get('name')!r} id={cred.get('id')}")

    from globus_compute_sdk import Client

    eps = Client(do_version_check=False).get_endpoints()
    print(f"\nendpoints owned by this identity: {len(eps)}")
    for ep in eps:
        print(f"  {ep.get('name')}  {ep.get('uuid')}")
    return 0


def cmd_provision(args: argparse.Namespace) -> int:
    """Build the thing being proposed, in whichever shape, so it can be measured."""
    ac = _auth_client()
    me = ac.userinfo()
    sites = [s.strip() for s in args.sites.split(",") if s.strip()]
    if not sites:
        sys.exit("--sites must name at least one site")

    contact = args.contact_email or me.get("email")
    if not contact:
        sys.exit("Could not determine a contact email; pass --contact-email")

    project_name = f"{args.federation} (APPFL sharing probe)"
    try:
        existing = {p["display_name"]: p for p in ac.get_projects()}
    except Exception as e:  # noqa: BLE001
        _explain_session_error(e)
        raise
    if project_name in existing:
        project = existing[project_name]
        print(f"reusing project {project_name!r} id={project['id']}")
    else:
        # admin_ids is not optional: Auth rejects a project with no administrator. Worth
        # noticing rather than working around -- "who may mint and revoke federation
        # secrets" is exactly the governance question this design raises, and Auth makes
        # you answer it at creation time. In a real deployment pass admin_group_ids
        # instead, so the answer is a Globus group and not one person's account.
        project = ac.create_project(project_name, contact_email=contact, admin_ids=me["sub"])[
            "project"
        ]
        print(f"created project {project_name!r} id={project['id']}")

    state = {
        "federation": args.federation,
        "shape": args.shape,
        "project_id": project["id"],
        "created_by": me.get("sub"),
        "sites": {},
    }

    def _new_client(name: str) -> dict:
        try:
            res = ac.create_client(
                name=name,
                project=project["id"],
                client_type="client_identity",
                visibility="private",
            )["client"]
        except Exception as e:  # noqa: BLE001
            _explain_session_error(e)
            raise
        print(f"created client identity {name!r} id={res['id']}")
        return res

    if args.shape == "shared":
        # The proposal as written: ONE identity, N secrets. Every endpoint below will be
        # owned by this single client.
        shared = _new_client(f"{args.federation}-federation")
        for site in sites:
            cred = ac.create_client_credential(shared["id"], f"{site}")["credential"]
            state["sites"][site] = {
                "client_id": shared["id"],
                "client_name": shared["name"],
                "credential_id": cred["id"],
                "secret": cred["secret"],
            }
            print(f"  minted secret for {site}: {_fingerprint(cred['secret'])}")
    else:
        # The alternative: one identity per site. Same "no root" property, but a secret
        # unlocks exactly one site.
        for site in sites:
            c = _new_client(f"{args.federation}-{site}")
            cred = ac.create_client_credential(c["id"], f"{site}")["credential"]
            state["sites"][site] = {
                "client_id": c["id"],
                "client_name": c["name"],
                "credential_id": cred["id"],
                "secret": cred["secret"],
            }
            print(f"  minted secret for {site}: {_fingerprint(cred['secret'])}")

    _save_state(state)
    print(
        f"\nstate written to {_state_file(args.federation)} (0600). "
        "Secrets are in it; it is gitignored."
    )
    print(
        f"shape={args.shape}: "
        f"{len({s['client_id'] for s in state['sites'].values()})} client identity/ies "
        f"for {len(sites)} sites"
    )
    return 0


def cmd_env(args: argparse.Namespace) -> int:
    """Print the two exports a partner would really set. `eval` this."""
    creds = _site_creds(_load_state(args.federation), args.as_site)
    print(f"export GLOBUS_COMPUTE_CLIENT_ID={creds['client_id']}")
    print(f"export GLOBUS_COMPUTE_CLIENT_SECRET={creds['secret']}")
    return 0


def _forget_tokens(client_id: str) -> int:
    """Drop cached tokens for a client identity, so the next call must re-authenticate.

    Needed to test revocation honestly. Without it a "revoked" credential keeps working
    from the local cache and the test proves nothing -- which is itself the finding:
    deleting a secret stops new tokens being issued, it does not invalidate ones already
    held.

    Note the namespace is keyed on the CLIENT id, not the credential
    (`sdk/auth/token_storage.py`). Under the shared shape every site therefore shares one
    cache entry, and nothing on disk records which partner's secret obtained it.
    """
    import sqlite3

    db = Path.home() / ".globus_compute" / "storage.db"
    if not db.exists():
        return 0
    con = sqlite3.connect(db)
    with con:
        cur = con.execute(
            "DELETE FROM token_storage WHERE namespace LIKE ?", (f"clientprofile/%/{client_id}",)
        )
    n = cur.rowcount
    con.close()
    return n


def _client_as(federation: str, site: str, fresh: bool = False):
    """A Compute Client authenticated as `site`'s credential, in-process."""
    state = _load_state(federation)
    creds = _site_creds(state, site)
    os.environ["GLOBUS_COMPUTE_CLIENT_ID"] = creds["client_id"]
    os.environ["GLOBUS_COMPUTE_CLIENT_SECRET"] = creds["secret"]
    if fresh:
        _forget_tokens(creds["client_id"])

    from globus_compute_sdk import Client

    return Client(do_version_check=False), creds


def cmd_whoami(args: argparse.Namespace) -> int:
    """Which identity does this site's credential actually authenticate as?

    Under `shared` every site prints the SAME identity. That is not a quirk of the probe;
    it is the design, and it is why there is no per-partner attribution downstream.
    """
    site = args.as_site
    client, creds = _client_as(args.federation, site, fresh=args.fresh)
    from globus_compute_sdk.sdk.auth.auth_client import ComputeAuthClient
    from globus_compute_sdk.sdk.auth.globus_app import get_globus_app

    info = ComputeAuthClient(app=get_globus_app()).userinfo()
    print(f"site           : {site}")
    print(f"secret         : {_fingerprint(creds['secret'])}")
    print(f"credential id  : {creds['credential_id']}")
    print("authenticates as")
    print(f"  username     : {info.get('preferred_username')}")
    print(f"  identity id  : {info.get('sub')}")
    return 0


def cmd_endpoints(args: argparse.Namespace) -> int:
    """What endpoints can this site's credential see -- and therefore submit to?

    This is the blast-radius measurement. Register an endpoint under each site's
    credential first, then run this for each site and compare the lists.
    """
    site = args.as_site
    client, creds = _client_as(args.federation, site, fresh=args.fresh)
    eps = client.get_endpoints()
    print(f"as {site} ({_fingerprint(creds['secret'])}): {len(eps)} endpoint(s) visible")
    for ep in eps:
        print(f"  {ep.get('name')}  uuid={ep.get('uuid')}  owner={ep.get('owner')}")
    return 0


def _probe_payload() -> dict:
    """Runs on the worker. Deliberately reports things a data owner would mind leaking."""
    import getpass
    import os
    import socket

    return {
        "host": socket.gethostname(),
        "posix_user": getpass.getuser(),
        "home": os.path.expanduser("~"),
        "cwd": os.getcwd(),
    }


def cmd_submit(args: argparse.Namespace) -> int:
    """Execute a function on an endpoint using one site's credential.

    The point is to run this with `--as` naming a site that did NOT register `--to`. If it
    returns, then holding any federation secret is sufficient to run code on any partner's
    machine, and the returned `posix_user` names the account it ran as -- which under this
    design is the partner who started the endpoint, not a mapped service account.
    """
    site = args.as_site
    client, creds = _client_as(args.federation, site, fresh=args.fresh)

    from globus_compute_sdk import Executor
    from globus_compute_sdk.serialize import CombinedCode, ComputeSerializer

    ex = Executor(endpoint_id=args.to, client=client)
    ex.serializer = ComputeSerializer(strategy_code=CombinedCode())
    print(f"submitting as {site} ({_fingerprint(creds['secret'])}) -> {args.to}")
    try:
        result = ex.submit(_probe_payload).result(timeout=args.timeout)
    except Exception as e:  # noqa: BLE001 -- the failure IS the result we are measuring
        print(f"REFUSED: {type(e).__name__}: {e}")
        return 1
    finally:
        ex.shutdown(wait=False, cancel_futures=True)
    print("ACCEPTED. The task ran. It reported:")
    for k, v in result.items():
        print(f"  {k:12s} {v}")
    return 0


def cmd_revoke(args: argparse.Namespace) -> int:
    """Delete one site's secret -- the revocation story Globus suggested.

    Prints what to re-check afterwards, because the interesting part is not that the
    delete succeeds. It is (a) whether the *other* sites are undisturbed, and (b) how
    long the revoked site keeps working on tokens it already holds.
    """
    state = _load_state(args.federation)
    site = args.as_site
    creds = _site_creds(state, site)
    ac = _auth_client()
    ac.delete_client_credential(creds["client_id"], creds["credential_id"])
    print(f"deleted credential {creds['credential_id']} ({site})")

    state["sites"][site]["revoked"] = True
    _save_state(state)

    others = [s for s in state["sites"] if s != site and not state["sites"][s].get("revoked")]
    print("\nnow check, in this order:")
    print(f"  whoami --as {site}        # should fail to get a token")
    for o in others[:1]:
        print(f"  whoami --as {o}        # should be unaffected")
    print(f"  endpoints --as {site}     # NOTE: an already-issued access token stays")
    print("                             # valid until it expires. Revocation is not")
    print("                             # instantaneous; delete the endpoint too if")
    print("                             # you need it to stop now.")
    if state["shape"] == "shared":
        print(f"\n  ALSO: {site}'s endpoint is owned by the shared client identity, which")
        print("  still exists. Revoking the secret removes their ability to authenticate;")
        print("  it does not change who owns the endpoints they registered.")
    return 0


def cmd_teardown(args: argparse.Namespace) -> int:
    """Delete everything `provision` made. Refuses to touch anything else."""
    state = _load_state(args.federation)
    if state["federation"] != args.federation:
        sys.exit(f"State is for federation {state['federation']!r}, not {args.federation!r}")
    ac = _auth_client()

    for client_id in {s["client_id"] for s in state["sites"].values()}:
        try:
            ac.delete_client(client_id)
            print(f"deleted client {client_id}")
        except Exception as e:  # noqa: BLE001 -- report and continue; partial teardown is worse
            print(f"could not delete client {client_id}: {e}")
    try:
        ac.delete_project(state["project_id"])
        print(f"deleted project {state['project_id']}")
    except Exception as e:  # noqa: BLE001
        print(f"could not delete project {state['project_id']}: {e}")

    _state_file(args.federation).unlink(missing_ok=True)
    print(f"removed {_state_file(args.federation)}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("inspect", help="read-only: projects, clients, credentials, endpoints")

    p = sub.add_parser("provision", help="create the project, client identity/ies and secrets")
    p.add_argument("--federation", required=True, help="name prefix for everything created")
    p.add_argument("--sites", required=True, help="comma-separated site names")
    p.add_argument(
        "--shape",
        choices=("shared", "per-site"),
        default="shared",
        help="shared: one identity for the federation (as proposed). "
        "per-site: one identity per site (the alternative).",
    )
    p.add_argument("--contact-email", default=None, help="project contact; defaults to yours")

    for name, fn, helptext in (
        ("env", cmd_env, "print the partner's two exports; eval this"),
        ("whoami", cmd_whoami, "which identity does this site's secret authenticate as"),
        ("endpoints", cmd_endpoints, "what endpoints can this site's secret reach"),
        ("revoke", cmd_revoke, "delete this site's secret"),
    ):
        sp = sub.add_parser(name, help=helptext)
        sp.add_argument("--as", dest="as_site", required=True, metavar="SITE")
        sp.add_argument("--federation", default="probe", help="which provisioned set to use")
        sp.set_defaults(func=fn)

    sp = sub.add_parser("submit", help="run a function on an endpoint as one site")
    sp.add_argument("--as", dest="as_site", required=True, metavar="SITE")
    sp.add_argument("--federation", default="probe", help="which provisioned set to use")
    sp.add_argument("--to", required=True, metavar="UUID", help="endpoint to submit to")
    sp.add_argument("--timeout", type=float, default=120.0)
    sp.add_argument(
        "--fresh",
        action="store_true",
        help="discard cached tokens first, forcing a real re-authentication",
    )
    sp.set_defaults(func=cmd_submit)

    t = sub.add_parser("teardown", help="delete everything provision created")
    t.add_argument("--federation", required=True)

    args = parser.parse_args(argv)
    dispatch = {
        "inspect": cmd_inspect,
        "provision": cmd_provision,
        "teardown": cmd_teardown,
    }
    fn = dispatch.get(args.command) or args.func
    return fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
