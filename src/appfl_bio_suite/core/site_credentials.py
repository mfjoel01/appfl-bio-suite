"""Experimental: Coordinator-only credential routing for single-user Compute endpoints.

The manifest contains public client IDs and environment-variable references, never
secrets. Nothing from this module is added to the configuration sent to workers.
"""

from __future__ import annotations

import json
import os
import re
from contextvars import ContextVar
from dataclasses import dataclass
from pathlib import Path
from uuid import UUID

CREDENTIALS_ENV = "APPFL_BIO_SUITE_GLOBUS_CREDENTIALS"


@dataclass(frozen=True)
class SiteCredential:
    client_id: str
    client_secret_env: str
    expected_user: str

    def client(self):
        from globus_compute_sdk import Client
        from globus_compute_sdk.serialize import CombinedCode
        from globus_sdk import ClientApp, GlobusAppConfig
        from globus_sdk.token_storage import MemoryTokenStorage

        secret = os.environ.get(self.client_secret_env)
        if not secret:
            raise ValueError(f"Set coordinator environment variable {self.client_secret_env}.")
        app = ClientApp(
            client_id=self.client_id,
            client_secret=secret,
            config=GlobusAppConfig(token_storage=MemoryTokenStorage()),
        )
        return Client(app=app, code_serialization_strategy=CombinedCode())


def _unique_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("Duplicate key in site credential manifest.")
        result[key] = value
    return result


def load_site_credentials(path: str | Path, endpoint_ids) -> dict[str, SiteCredential]:
    """Validate all requested routes before creating clients or submitting any tasks."""
    try:
        raw = json.loads(Path(path).read_text(), object_pairs_hook=_unique_keys)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("Cannot read site credential manifest as JSON.") from exc
    if not isinstance(raw, dict) or set(raw) != {"version", "endpoints"} or raw["version"] != 1:
        raise ValueError("Site credential manifest requires version: 1 and endpoints.")
    if not isinstance(raw["endpoints"], dict):
        raise ValueError("Manifest endpoints must be an object keyed by endpoint UUID.")
    credentials = {}
    identities = set()
    for endpoint, entry in raw["endpoints"].items():
        try:
            if str(UUID(endpoint)) != endpoint:
                raise ValueError()
            if not isinstance(entry, dict) or set(entry) != {
                "client_id",
                "client_secret_env",
                "expected_user",
            }:
                raise ValueError()
            identity = str(UUID(entry["client_id"]))
            if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", entry["client_secret_env"]):
                raise ValueError()
            if not isinstance(entry["expected_user"], str) or not entry["expected_user"].strip():
                raise ValueError()
        except (ValueError, TypeError, AttributeError) as exc:
            raise ValueError(
                "Each endpoint route requires UUIDs, client_secret_env, and expected_user. "
                "Inline secrets and extra fields are not accepted."
            ) from exc
        if identity in identities:
            raise ValueError("Use a distinct service-account client ID for each endpoint route.")
        identities.add(identity)
        credentials[endpoint] = SiteCredential(
            identity, entry["client_secret_env"], entry["expected_user"]
        )
    selected = set(endpoint_ids)
    missing = selected - credentials.keys()
    if missing:
        raise ValueError("Missing site credentials for endpoint(s): " + ", ".join(sorted(missing)))
    for endpoint in selected:
        name = credentials[endpoint].client_secret_env
        if not os.environ.get(name):
            raise ValueError(f"Set coordinator environment variable {name}.")
    return {endpoint: credentials[endpoint] for endpoint in selected}


class SiteExecutorRouter:
    """APPFL's executor interface, with one authenticated executor per endpoint.

    APPFL 1.10 sets endpoint_id immediately before submit. The selection is local to
    the calling context; submission and result watching belong to the selected SDK
    executor. The central APPFL scheduler and aggregation loop remain unchanged.
    """

    def __init__(self, credentials: dict[str, SiteCredential]):
        from globus_compute_sdk import Executor

        self._selected = ContextVar("site_endpoint", default=None)
        self.executors = {}
        self._closed = False
        try:
            for endpoint, credential in credentials.items():
                executor = Executor(endpoint_id=endpoint, client=credential.client())
                self.executors[endpoint] = executor
                # SDK 4.9 inherits a process-global task group across client identities.
                executor.task_group_id = None
        except Exception:
            self.shutdown(wait=False, cancel_futures=True)
            raise

    @property
    def endpoint_id(self):
        return self._selected.get()

    @endpoint_id.setter
    def endpoint_id(self, endpoint):
        if endpoint not in self.executors:
            raise ValueError("No authenticated executor configured for this endpoint.")
        self._selected.set(endpoint)

    def submit(self, function, *args, **kwargs):
        if self._closed:
            raise RuntimeError("Site executor router is shut down.")
        endpoint = self.endpoint_id
        if endpoint not in self.executors:
            raise ValueError("Select a configured endpoint before submitting.")
        return self.executors[endpoint].submit(function, *args, **kwargs)

    def shutdown(self, **kwargs):
        self._closed = True
        errors = []
        for executor in self.executors.values():
            try:
                executor.shutdown(**kwargs)
            except Exception as exc:
                errors.append(exc)
        if errors:
            raise RuntimeError(f"Could not shut down {len(errors)} site executor(s).") from errors[
                0
            ]


def site_communicator(credentials, **kwargs):
    """Use APPFL's existing round/result logic with explicit authenticated routing."""
    from appfl.comm.globus_compute import GlobusComputeServerCommunicator

    router = SiteExecutorRouter(credentials)

    class RoutedCommunicator(GlobusComputeServerCommunicator):
        def _load_gce(self, **unused):
            self.gce = router

    try:
        return RoutedCommunicator(**kwargs)
    except Exception:
        router.shutdown(wait=False, cancel_futures=True)
        raise
