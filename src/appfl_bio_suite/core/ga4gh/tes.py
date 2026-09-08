"""TES -- Task Execution Service: dispatch the site stage to a site that speaks TES.

WHERE THIS FITS, HONESTLY
-------------------------
Globus Compute is this federation's transport and remains the default. TES is a second
one, and it earns its place for two reasons rather than for standards-compliance:

1. **A partner who already runs a TES service can join without an endpoint.** Funnel,
   TESK, and several cloud services expose one. That partner does not have to stand up a
   multi-user Globus Compute endpoint, arrange an identity mapping, or install this
   suite -- they accept a task that names a container image and a DRS input.

2. **It forced the site stage to have a command-line form.** Expressing the computation
   as a TES task means it must be a command over files, which is what
   ``experiments/fine_mapping/site_stage.py`` now is. That entry point is testable from a
   shell, runnable under CWL, and reusable by any workflow engine -- none of which was
   true of a computation that existed only as an APPFL trainer object.

WHAT A TES RUN NEEDS THAT THE GLOBUS PATH DOES NOT
--------------------------------------------------
* **A published container image.** The TRS tool version names one; ``ga4gh trs publish``
  generates the Containerfile, and somebody has to build and push it. Until that happens
  the TES path has nothing to run, and ``preflight`` says so rather than failing at
  submission.
* **A place the outputs land that the coordinator can read.** TES stages outputs to a URL
  the server can write -- a shared filesystem path, an S3 bucket, whatever that
  deployment uses. The driver reads the aggregate payloads back from there.

Neither is a limitation of this implementation; they are what running somebody else's
container on somebody else's cluster costs.

THE TASK DOCUMENT IS THE INTERFACE
----------------------------------
``build_site_task`` produces a complete TES task document, and ``ga4gh tes task`` writes
one out without submitting anything. That is deliberate: a partner evaluating whether to
accept this workload should be able to read exactly what will run -- image, command,
inputs, resource request, and the DUO decision that authorized it, which travels in the
task's tags.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field

__all__ = [
    "TES_VERSION",
    "TesExecutor",
    "TesInput",
    "TesOutput",
    "TesResources",
    "TesTask",
    "TesError",
    "TesClient",
    "TERMINAL_STATES",
    "build_site_task",
    "service_info",
    "write_task",
]

TES_VERSION = "1.1.0"

# TES state machine. A task in one of these will not change again, which is the only
# thing a poll loop needs to know.
TERMINAL_STATES = frozenset(
    {"COMPLETE", "EXECUTOR_ERROR", "SYSTEM_ERROR", "CANCELED", "CANCELLED", "PREEMPTED"}
)

TesState = Literal[
    "UNKNOWN",
    "QUEUED",
    "INITIALIZING",
    "RUNNING",
    "PAUSED",
    "COMPLETE",
    "EXECUTOR_ERROR",
    "SYSTEM_ERROR",
    "CANCELED",
    "CANCELING",
    "PREEMPTED",
]


class TesError(RuntimeError):
    """A TES service rejected a task, or could not be reached."""


class _Model(BaseModel):
    model_config = ConfigDict(extra="allow")


class TesExecutor(_Model):
    """One container invocation. Tasks may have several; this experiment uses one."""

    image: str
    command: list[str]
    workdir: str | None = None
    stdin: str | None = None
    stdout: str | None = None
    stderr: str | None = None
    env: dict[str, str] = Field(default_factory=dict)
    ignore_error: bool = False


class TesInput(_Model):
    """Something staged into the task's filesystem before the executor runs.

    ``url`` may be a ``drs://`` URI -- TES servers are expected to resolve those, and it
    is the join between this specification and DRS. ``content`` is the alternative for
    small inline files, which is how the run config reaches the container without needing
    a shared filesystem for it.
    """

    path: str
    name: str | None = None
    description: str | None = None
    url: str | None = None
    type: Literal["FILE", "DIRECTORY"] = "FILE"
    content: str | None = None
    streamable: bool | None = None


class TesOutput(_Model):
    path: str
    url: str
    name: str | None = None
    description: str | None = None
    type: Literal["FILE", "DIRECTORY"] = "FILE"
    path_prefix: str | None = None


class TesResources(_Model):
    cpu_cores: int | None = None
    preemptible: bool | None = None
    ram_gb: float | None = None
    disk_gb: float | None = None
    zones: list[str] = Field(default_factory=list)
    backend_parameters: dict[str, str] = Field(default_factory=dict)
    backend_parameters_strict: bool | None = None


class TesTask(_Model):
    """A TES 1.1 task. ``id``, ``state`` and ``logs`` are server-assigned."""

    executors: list[TesExecutor]
    id: str | None = None
    state: TesState | None = None
    name: str | None = None
    description: str | None = None
    inputs: list[TesInput] = Field(default_factory=list)
    outputs: list[TesOutput] = Field(default_factory=list)
    resources: TesResources | None = None
    volumes: list[str] = Field(default_factory=list)
    tags: dict[str, str] = Field(default_factory=dict)
    creation_time: str | None = None
    logs: list[dict[str, Any]] = Field(default_factory=list)

    def request_body(self) -> dict[str, Any]:
        """The document to POST: server-assigned fields stripped.

        A server will usually ignore ``id`` and ``state`` on a create, but "usually" is
        not a property to rely on across four implementations.
        """
        body = json.loads(self.model_dump_json(exclude_none=True))
        for key in ("id", "state", "logs", "creation_time"):
            body.pop(key, None)
        return body

    def render(self) -> str:
        lines = [f"{self.name or '(unnamed task)'}"]
        if self.state:
            lines.append(f"  state    {self.state}")
        for executor in self.executors:
            lines.append(f"  image    {executor.image}")
            lines.append(f"  command  {' '.join(executor.command)}")
        for entry in self.inputs:
            source = entry.url or "(inline content)"
            lines.append(f"  in       {entry.path} <- {source}")
        for entry in self.outputs:
            lines.append(f"  out      {entry.path} -> {entry.url}")
        if self.resources:
            resources = self.resources
            lines.append(
                f"  ask      {resources.cpu_cores or '?'} cpu, "
                f"{resources.ram_gb or '?'} GB ram, {resources.disk_gb or '?'} GB disk"
            )
        for key, value in sorted(self.tags.items()):
            lines.append(f"  tag      {key}={value}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# building this experiment's task
# ---------------------------------------------------------------------------

# Where the container sees things. Fixed rather than configurable: these are paths inside
# an image this project generates, and a knob here would only let a coordinator break the
# agreement between the task document and the CWL descriptor.
CONTAINER_BUNDLE = "/data/bundle"
CONTAINER_CONFIG = "/config/site.json"
CONTAINER_OUT = "/out"


def build_site_task(
    client_id: str,
    experiment: str,
    image: str,
    bundle_url: str,
    outputs_url: str,
    run_config: dict[str, Any],
    cpu_cores: int = 4,
    ram_gb: float = 16.0,
    disk_gb: float = 64.0,
    preemptible: bool = False,
    tags: dict[str, str] | None = None,
) -> TesTask:
    """One site's task: read its bundle, write its aggregates, touch nothing else.

    ``bundle_url`` is normally a ``drs://`` URI, which the TES server resolves and stages
    at :data:`CONTAINER_BUNDLE`. ``run_config`` is inlined rather than staged from a URL:
    it is a few hundred bytes, it differs per site, and putting it in an object store
    would create one more thing to keep in step with ``federation.yaml``.
    """
    return TesTask(
        name=f"{experiment} site stage: {client_id}",
        description=(
            f"Compute {client_id}'s aggregate statistics for the federated {experiment} "
            "experiment. Reads one site's own cohort; returns second moments only. No "
            "genotype or phenotype value is written to the outputs."
        ),
        inputs=[
            TesInput(
                name="bundle",
                description="This site's own data bundle, addressed by DRS.",
                url=bundle_url,
                path=CONTAINER_BUNDLE,
                type="DIRECTORY",
            ),
            TesInput(
                name="run-config",
                description="Locus shard, ancestry columns, and the run's data use request.",
                path=CONTAINER_CONFIG,
                content=json.dumps(run_config, indent=2),
            ),
        ],
        outputs=[
            TesOutput(
                name="aggregates",
                description="Site aggregates and the GA4GH provenance record.",
                path=CONTAINER_OUT,
                url=outputs_url,
                type="DIRECTORY",
            )
        ],
        resources=TesResources(
            cpu_cores=cpu_cores,
            ram_gb=ram_gb,
            disk_gb=disk_gb,
            preemptible=preemptible,
        ),
        executors=[
            TesExecutor(
                image=image,
                command=[
                    "appfl-bio-suite",
                    "site-stage",
                    experiment,
                    "--config",
                    CONTAINER_CONFIG,
                    "--data-dir",
                    CONTAINER_BUNDLE,
                    "--out",
                    CONTAINER_OUT,
                ],
                workdir=CONTAINER_OUT,
                stdout=f"{CONTAINER_OUT}/stdout.log",
                stderr=f"{CONTAINER_OUT}/stderr.log",
                # The same caps the Globus Compute workers get, for the same reason: a
                # container scheduled onto a busy node with an unbounded BLAS thread pool
                # dies inside numpy's import.
                env={
                    "OMP_NUM_THREADS": "1",
                    "OPENBLAS_NUM_THREADS": "1",
                    "MKL_NUM_THREADS": "1",
                    "PYTHONNOUSERSITE": "1",
                },
            )
        ],
        # Tags are the only part of a TES task a server is required to hand back
        # untouched, so the provenance that must survive the round trip goes here: which
        # tool version, which shard, and the data use decision that authorized the task.
        tags={
            "experiment": experiment,
            "client_id": client_id,
            **(tags or {}),
        },
    )


def service_info() -> dict[str, Any]:
    """service-info for this suite acting as a TES *client*, used in provenance."""
    from appfl_bio_suite import __version__

    return {
        "id": "org.ga4gh.tes.client.appfl-bio-suite",
        "name": "appfl-bio-suite TES client",
        "type": {"group": "org.ga4gh", "artifact": "tes", "version": TES_VERSION},
        "version": __version__,
    }


# ---------------------------------------------------------------------------
# the client
# ---------------------------------------------------------------------------


class TesClient:
    """A TES 1.1 client. Standard library only.

    Bearer token from the constructor or ``$TES_TOKEN``. Nothing here writes a token to
    disk or into a task document -- a task is often visible to the whole service.
    """

    def __init__(self, base_url: str, token: str | None = None, timeout: float = 60.0):
        import os

        self.base_url = base_url.rstrip("/")
        self.token = token or os.environ.get("TES_TOKEN")
        self.timeout = timeout

    # -- transport --------------------------------------------------------

    def _request(self, method: str, path: str, payload: dict[str, Any] | None = None) -> Any:
        import urllib.error
        import urllib.request

        url = f"{self.base_url}{path}"
        data = json.dumps(payload).encode("utf-8") if payload is not None else None
        request = urllib.request.Request(url, data=data, method=method)
        request.add_header("Accept", "application/json")
        if data is not None:
            request.add_header("Content-Type", "application/json")
        if self.token:
            request.add_header("Authorization", f"Bearer {self.token}")
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                body = response.read().decode("utf-8")
                return json.loads(body) if body else {}
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode("utf-8", "replace")[:600]
            raise TesError(
                f"{method} {url} returned {exc.code}.\n{detail}\n"
                + (
                    "401/403 usually means the token is missing or not scoped for this "
                    "service; set $TES_TOKEN.\n"
                    if exc.code in (401, 403)
                    else ""
                )
            ) from exc
        except urllib.error.URLError as exc:
            raise TesError(
                f"could not reach {url}: {exc.reason}\n"
                "Check the service URL in federation.yaml (`ga4gh.tes.url`) and that it "
                "is reachable from the machine running the driver."
            ) from exc

    # -- the API ----------------------------------------------------------

    def service_info(self) -> dict[str, Any]:
        return self._request("GET", "/ga4gh/tes/v1/service-info")

    def create_task(self, task: TesTask) -> str:
        """Submit a task. Returns the server-assigned id."""
        response = self._request("POST", "/ga4gh/tes/v1/tasks", task.request_body())
        task_id = response.get("id")
        if not task_id:
            raise TesError(f"the service accepted the task but returned no id: {response}")
        return str(task_id)

    def get_task(self, task_id: str, view: str = "MINIMAL") -> TesTask:
        return TesTask.model_validate(
            self._request("GET", f"/ga4gh/tes/v1/tasks/{task_id}?view={view}")
        )

    def cancel_task(self, task_id: str) -> None:
        self._request("POST", f"/ga4gh/tes/v1/tasks/{task_id}:cancel", {})

    def list_tasks(self, state: str | None = None, name_prefix: str | None = None) -> list[TesTask]:
        query = []
        if state:
            query.append(f"state={state}")
        if name_prefix:
            query.append(f"name_prefix={name_prefix}")
        suffix = ("?" + "&".join(query)) if query else ""
        payload = self._request("GET", f"/ga4gh/tes/v1/tasks{suffix}")
        return [TesTask.model_validate(entry) for entry in payload.get("tasks", [])]

    def wait(
        self,
        task_id: str,
        poll_seconds: float = 15.0,
        timeout_seconds: float | None = None,
        on_state: Any = None,
    ) -> TesTask:
        """Poll until the task reaches a terminal state, or the timeout elapses.

        On timeout the task is **cancelled** before raising. A federated run that gave up
        waiting while a site kept computing would leave the partner burning an allocation
        on a result nobody will read.
        """
        import time

        started = time.time()
        last: str | None = None
        while True:
            task = self.get_task(task_id, view="MINIMAL")
            state = str(task.state or "UNKNOWN")
            if state != last and on_state is not None:
                on_state(task_id, state)
            last = state
            if state in TERMINAL_STATES:
                return self.get_task(task_id, view="FULL")
            if timeout_seconds is not None and (time.time() - started) > timeout_seconds:
                try:
                    self.cancel_task(task_id)
                except TesError:
                    pass
                raise TesError(
                    f"task {task_id} was still {state} after {timeout_seconds:.0f}s; "
                    "cancelled it. Raise the timeout, or check the service's queue."
                )
            time.sleep(poll_seconds)


def write_task(task: TesTask, path: str | Path) -> Path:
    """Write a task document to disk, for review or for `curl -d @task.json`."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(task.request_body(), indent=2) + "\n", encoding="utf-8")
    return path
