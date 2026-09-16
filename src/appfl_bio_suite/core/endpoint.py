"""Talk to Globus Compute endpoints: status, and a real round-trip smoke test.

Consolidates three near-identical probe scripts that had grown independently across two
projects. They agreed on the mechanism -- submit a trivial function, see what comes back --
and differed only in which endpoint UUID was hardcoded and how much diagnostic prose was
attached. The prose was the valuable part, so it is all here.

WHY A ROUND TRIP AND NOT A STATUS CHECK
---------------------------------------
``get_endpoint_status`` tells you a daemon is running. It does not tell you that your
identity is accepted, that identity mapping resolves, that the scheduler accepts the job,
that the environment activates, or that the worker can talk back. Every one of those has
failed independently on this project while the endpoint reported ``online``.

The round trip tests the whole chain, and the returned payload names the POSIX account the
task actually ran as -- which is the only real proof that a partner's identity mapping
took effect. "The endpoint started fine" is not evidence: a multi-user endpoint started by
an unprivileged user starts successfully, logs a warning nobody reads, and silently
ignores its identity mapping entirely.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

__all__ = [
    "EndpointStatus",
    "SmokeResult",
    "get_status",
    "smoke_test",
    "where_am_i",
    "explain_submit_failure",
]

# Long enough that a queued scheduler job is not mistaken for a failure. A partner's
# block may genuinely sit in the queue for minutes; that is not an error, and reporting
# it as one sends everyone chasing the wrong thing.
DEFAULT_TIMEOUT = 300.0


def where_am_i() -> dict[str, str]:
    """Run on the worker and report where it landed.

    Imports live inside the function because this is serialized and shipped to a remote
    interpreter -- it must not close over anything from the sending process. That is
    necessary but not sufficient: see :func:`smoke_test` for why the body being
    self-contained does not by itself make the function shippable.

    ``user`` is the load-bearing field. On a partner's multi-user endpoint it should be
    the experiment's service account, not the account that started the endpoint. If it is
    the latter, identity mapping was skipped and the endpoint is accepting you for the
    wrong reason.
    """
    import getpass
    import os
    import platform
    import sys

    return {
        "node": platform.node(),
        "user": getpass.getuser(),
        "cwd": os.getcwd(),
        "python": sys.version.split()[0],
    }


@dataclass
class EndpointStatus:
    uuid: str
    label: str
    status: str
    error: str | None = None

    @property
    def online(self) -> bool:
        return self.status == "online"

    def render(self) -> str:
        if self.error:
            return f"  {self.label:24} {self.uuid}  ERROR: {self.error}"
        mark = "online " if self.online else self.status
        return f"  {self.label:24} {self.uuid}  {mark}"


@dataclass
class SmokeResult:
    uuid: str
    label: str
    ok: bool
    elapsed: float
    payload: dict[str, Any] | None = None
    error: str | None = None
    diagnosis: str | None = None

    def render(self) -> str:
        head = f"{self.label} ({self.uuid})"
        if self.ok:
            p = self.payload or {}
            return (
                f"PASS  {head}\n"
                f"      ran as user '{p.get('user')}' on node '{p.get('node')}'\n"
                f"      python {p.get('python')}   cwd {p.get('cwd')}\n"
                f"      round-trip {self.elapsed:.1f}s"
            )
        out = [f"FAIL  {head}", f"      {self.error}"]
        if self.diagnosis:
            out.append("")
            out.extend(f"      {line}" for line in self.diagnosis.splitlines())
        return "\n".join(out)


def _client():
    from globus_compute_sdk import Client

    return Client()


def get_status(uuid: str, label: str = "") -> EndpointStatus:
    """Ask the Globus Compute service whether an endpoint is online.

    The cloud status is authoritative. Local ``globus-compute-endpoint list`` reads a
    pidfile that, on a shared filesystem, may record a PID from a different login node --
    so it can report an endpoint as running when it is not, and vice versa.
    """
    try:
        raw = _client().get_endpoint_status(uuid)
        return EndpointStatus(uuid=uuid, label=label or uuid, status=raw.get("status", "?"))
    except Exception as exc:
        return EndpointStatus(
            uuid=uuid, label=label or uuid, status="error", error=f"{type(exc).__name__}: {exc}"
        )


def smoke_test(uuid: str, label: str = "", timeout: float = DEFAULT_TIMEOUT) -> SmokeResult:
    """Submit a trivial task and wait for it to come back.

    Runs no experiment code and touches no data, so it is safe to run against a partner's
    endpoint at any time -- including before they have staged anything.

    THE PROBE SHIPS ITS OWN SOURCE
    ------------------------------
    The SDK's default code strategy (``DillCode``) serializes a module-level function *by
    reference* -- the payload is the string ``appfl_bio_suite.core.endpoint.where_am_i``,
    85 bytes of it -- so the worker has to import ``appfl_bio_suite`` to resolve the name.
    Writing every import inside :func:`where_am_i` does nothing about that; the body is
    never reached, because unpickling fails first.

    That made this check strictly stricter than the run it is supposed to clear. An
    experiment's dataset and model are shipped to a worker as source text precisely so a
    partner need not install this package -- see
    ``experiments/flamby_heart_disease/dataset.py`` -- so an endpoint that would train
    perfectly well failed its smoke test with ``ModuleNotFoundError: No module named
    'appfl_bio_suite'``, and the message pointed at serialization rather than at the
    install it was really about.

    ``CombinedCode`` packs the by-reference form *and* the function's source text, and the
    worker uses whichever it can. Partners who installed the suite are unaffected; the one
    who did not now passes, which is the correct answer for them.
    """
    from globus_compute_sdk import Executor
    from globus_compute_sdk.serialize import CombinedCode, ComputeSerializer

    label = label or uuid
    started = time.time()
    executor = Executor(endpoint_id=uuid)
    executor.serializer = ComputeSerializer(strategy_code=CombinedCode())
    try:
        future = executor.submit(where_am_i)
        payload = future.result(timeout=timeout)
        return SmokeResult(
            uuid=uuid, label=label, ok=True, elapsed=time.time() - started, payload=payload
        )
    except Exception as exc:
        return SmokeResult(
            uuid=uuid,
            label=label,
            ok=False,
            elapsed=time.time() - started,
            error=f"{type(exc).__name__}: {exc}",
            diagnosis=explain_submit_failure(exc),
        )
    finally:
        # NOT the Executor's context manager. Its __exit__ calls shutdown(wait=True),
        # which blocks until every outstanding future has a result -- so on a timeout it
        # would sit there waiting for exactly the queued task the timeout gave up on, and
        # --timeout would never return. A queued scheduler job is the normal case this
        # command has to survive, so the shutdown must not wait for one.
        try:
            executor.shutdown(wait=False, cancel_futures=True)
        except Exception:  # noqa: BLE001 - cleanup must not replace the real result
            pass


def explain_submit_failure(exc: BaseException) -> str:
    """Turn a submit failure into the thing to actually go check.

    Every branch below is a failure that has really happened on this project, and the
    ordering is by observed frequency rather than by what seems most likely in the
    abstract. The 422 case in particular sends people to the regex, which is the one
    thing that is usually already correct.
    """
    text = f"{type(exc).__name__}: {exc}"

    if "422" in text or "Identity failed to map" in text:
        return (
            "422 -- the endpoint is privileged and its mapping IS being read; the\n"
            "expression just did not resolve your identity.\n"
            "\n"
            "Almost always: ^ or $ anchors in the `match` field. The mapper escapes\n"
            "them into literal characters and anchors the pattern itself, so an\n"
            "anchored pattern can never match. Have the partner run:\n"
            "\n"
            "    appfl-bio-suite identity validate <their mapping.json>\n"
            "\n"
            "Next most likely: they edited a different file from the one\n"
            "`identity_mapping_config_path` names -- after switching to root, that path\n"
            "may still point into a personal home directory.\n"
            "\n"
            "No restart is needed after they fix it; the endpoint re-reads the file\n"
            "within about five seconds."
        )

    if "403" in text or "ENDPOINT_ACCESS_FORBIDDEN" in text:
        return (
            "403 -- this is NOT a mapping-expression problem. Do not start editing the\n"
            "JSON.\n"
            "\n"
            "It almost always means the endpoint was started by an ordinary user rather\n"
            "than a privileged one. In that case identity mapping is skipped entirely\n"
            "and only the starting identity may submit -- and the endpoint still starts\n"
            "successfully, with no error, which is what makes this confusing.\n"
            "\n"
            "Have them confirm who owns the running process, restart it as the\n"
            "privileged user, then RE-READ THE UUID: endpoint IDs are per-account, so\n"
            "the one they sent you is a different, now-dead endpoint."
        )

    if "404" in text:
        return (
            "404 -- no such endpoint. The UUID is wrong or the endpoint was deleted.\n"
            "Note that `delete` + `configure` issues a NEW UUID; stop/start/restart\n"
            "preserve it. Ask for a fresh `globus-compute-endpoint list`, run as the\n"
            "same user that starts the endpoint."
        )

    if "TimeoutError" in text or "timeout" in text.lower():
        return (
            "Timed out waiting for a result -- but the submission was ACCEPTED, which\n"
            "means authentication and identity mapping both worked. The task is sitting\n"
            "in the partner's scheduler queue.\n"
            "\n"
            "This is usually fine: retry with a longer --timeout. If it never returns,\n"
            "check whether their block is starting at all (a busy partition), and\n"
            "whether workers register once it does -- 'workers failed to register' or\n"
            "'AF_UNIX path too long' means `export TMPDIR=/tmp` is missing from their\n"
            "worker_init."
        )

    if "SystemError" in text and "opcode" in text:
        return (
            "Bytecode could not be decoded on the worker. Driver and workers are on\n"
            "different Python MINOR versions -- this fails before any of your code runs.\n"
            "Standardize the whole federation on one minor version.\n"
            "\n"
            "Patch-level differences (3.12.11 vs 3.12.13) are fine and produce only a\n"
            "warning; do not make partners chase those."
        )

    return (
        "Unrecognized failure. Worth checking, in order: the endpoint's cloud status,\n"
        "that the partner's globus-compute-endpoint version matches the federation pin,\n"
        "and the endpoint daemon log at ~/.globus_compute/<name>/endpoint.log.\n"
        "\n"
        "See docs/coordinator/troubleshooting.md."
    )
