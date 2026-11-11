"""Versioned onboarding specifications and resumable site installation.

Bundles are configuration, not executable scripts. Local configuration has an explicit
allowlist; it cannot override the experiment, identity, or software versions.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Literal
from uuid import UUID

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator

from appfl_bio_suite import __version__


class Model(BaseModel):
    model_config = ConfigDict(extra="forbid")


# The specification carries the site's client config, which for a DRS-verified bundle
# includes a checksum per data file -- 2.3 MB for a fine-mapping cohort of 11,850
# phenotype files. The cap guards against a hostile archive, not a large honest one.
MAX_SPEC_BYTES = 64_000_000

# The service-account client an endpoint authenticates as, in single-user mode. These are
# the Compute SDK's own variable names; nothing else is read from a credentials file.
CREDENTIAL_VARS = ("GLOBUS_COMPUTE_CLIENT_ID", "GLOBUS_COMPUTE_CLIENT_SECRET")


class Bundle(Model):
    schema_version: Literal[1] = 1
    experiment: str
    client_id: str
    suite_version: str
    python_version: str
    install_spec: str
    coordinator_identity: str = Field(min_length=1)
    coordinator_identity_match: str = Field(min_length=1)
    service_account: str
    endpoint_name: str
    scheduler: Literal["slurm", "pbspro", "local"]
    endpoint_mode: Literal["multi-user", "single-user"] = "multi-user"
    versions: dict[str, str]
    defaults: dict[str, str]
    constraints: dict[str, str]
    experiment_digest: str
    client_config: dict
    extras: list[str] = Field(default_factory=list)
    wheel_file: str | None = None
    wheel_sha256: str | None = None

    @field_validator("wheel_file")
    @classmethod
    def wheel_basename(cls, value):
        if value is not None and not re.fullmatch(r"appfl_bio_suite-[A-Za-z0-9_.+-]+\.whl", value):
            raise ValueError("wheel_file must be an appfl_bio_suite wheel basename")
        return value

    @field_validator("extras")
    @classmethod
    def valid_extras(cls, values):
        if any(not re.fullmatch(r"[a-z0-9_-]+", value) for value in values):
            raise ValueError("invalid package extra")
        return values

    @field_validator("endpoint_name", "service_account")
    @classmethod
    def safe_name(cls, value):
        if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", value):
            raise ValueError("must be a simple account/endpoint name")
        return value

    @field_validator("versions", "constraints")
    @classmethod
    def safe_pins(cls, value):
        for name, version in value.items():
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", name) or not re.fullmatch(
                r"[A-Za-z0-9_.+!-]+", version
            ):
                raise ValueError("Package constraints must contain exact names and versions")
        return value

    @field_validator("install_spec")
    @classmethod
    def package_requirement(cls, value):
        if not value.startswith("appfl-bio-suite") or any(c in value for c in "\r\n"):
            raise ValueError("must be one appfl-bio-suite requirement")
        return value


class SiteConfig(Model):
    # Required for slurm and pbspro, absent for local; see check_scheduler_settings.
    account: str | None = None
    partition: str | None = None
    queue: str | None = None
    data_dir: Path
    work_dir: Path
    interface: str | None = None

    @field_validator("account", "partition", "queue", "interface")
    @classmethod
    def scheduler_value(cls, value):
        if value is not None and not re.fullmatch(r"[A-Za-z0-9_.:/@+-]+", value):
            raise ValueError("must be a nonempty scheduler identifier")
        return value

    @field_validator("data_dir", "work_dir")
    @classmethod
    def absolute_path(cls, value):
        if not value.is_absolute() or any(c in str(value) for c in "\r\n"):
            raise ValueError("must be an absolute path without newlines")
        return value


def scheduler_fields(scheduler: str) -> tuple[str, ...]:
    """The local settings a scheduler requires, in prompt order."""
    return {"slurm": ("account", "partition"), "pbspro": ("account", "queue")}.get(scheduler, ())


def without_blanks(values: dict) -> dict:
    """Drop unset values, so an empty bundle default reads as missing rather than invalid."""
    return {key: value for key, value in values.items() if value not in (None, "")}


def check_scheduler_settings(scheduler: str, site: SiteConfig) -> None:
    required = scheduler_fields(scheduler)
    missing = [field for field in required if not getattr(site, field)]
    extra = [
        field
        for field in ("account", "partition", "queue")
        if field not in required and getattr(site, field)
    ]
    if scheduler == "local" and extra:
        raise ValueError(
            f"A local site runs workers on this host and has no scheduler; "
            f"remove {', '.join(extra)}"
        )
    if missing or extra:
        message = f"{scheduler} needs {' and '.join(required)}"
        if missing:
            message += f" (missing: {', '.join(missing)})"
        if extra:
            message += f"; do not set {', '.join(extra)}"
        raise ValueError(message)


class Receipt(Model):
    schema_version: Literal[1] = 1
    bundle_digest: str
    experiment: str
    client_id: str
    endpoint_id: UUID
    site: SiteConfig
    runtime: dict
    state_dir: Path

    @field_validator("state_dir")
    @classmethod
    def absolute_state(cls, value):
        if not value.is_absolute():
            raise ValueError("state_dir must be absolute")
        return value

    status: Literal["locally_verified"] = "locally_verified"


def digest(bundle: Bundle) -> str:
    canonical = json.dumps(bundle.model_dump(mode="json"), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def constraints_path() -> Path:
    from appfl_bio_suite.core.experiments import repo_root

    source = repo_root() / "constraints.txt"
    return source if source.is_file() else Path(__file__).parents[1] / "resources/constraints.txt"


def write_specification(
    destination: Path, context: dict, federation, client_config: dict
) -> list[Path]:
    pins = {}
    for line in constraints_path().read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if "==" in line:
            name, version = line.split("==")
            pins[name] = version
    # Check every package involved in serialization and endpoint operation.
    required = (
        "globus-compute-sdk",
        "globus-compute-endpoint",
        "globus-compute-common",
        "globus-sdk",
        "globus-identity-mapping",
        "dill",
        "appfl",
        "torch",
        "numpy",
        "parsl",
        "click",
    )
    defaults = {}
    if context["scheduler"] != "local":
        defaults["account"] = context["account"]
        defaults[context["scheduler_slot"]] = context["scheduler_slot_value"]
    defaults["data_dir"] = context["data_dir"]
    defaults["work_dir"] = context["output_dir"]
    if context["interface"]:
        defaults["interface"] = context["interface"]
    bundle = Bundle(
        **{
            key: context[key]
            for key in (
                "experiment",
                "client_id",
                "python_version",
                "install_spec",
                "coordinator_identity",
                "coordinator_identity_match",
                "service_account",
                "endpoint_name",
                "scheduler",
                "endpoint_mode",
            )
        },
        suite_version=__version__,
        versions={name: pins[name] for name in required},
        defaults=defaults,
        extras=context["partner_extras"].split(",") if context["partner_extras"] else [],
        constraints=pins,
        client_config=client_config,
        experiment_digest=experiment_digest(federation, context["experiment"]),
    )
    spec_path = destination / "onboarding.json"
    write_json(spec_path, bundle.model_dump())
    local = destination / "site.yaml"
    local.write_text(yaml.safe_dump(defaults, sort_keys=False))
    constraints = destination / "runtime-constraints.txt"
    constraints.write_text(constraints_path().read_text())
    return [spec_path, local, constraints]


def experiment_digest(federation, experiment: str) -> str:
    exp = federation.experiment(experiment)
    settings = exp.model_dump(mode="json", exclude={"sites"})
    # Registration changes UUIDs and paths, but membership/data assignments are immutable.
    settings["sites"] = [
        entry.model_dump(mode="json", exclude={"endpoint_uuid", "data_dir", "output_dir"})
        for entry in exp.sites
    ]
    return hashlib.sha256(json.dumps(settings, sort_keys=True).encode()).hexdigest()


def read_bundle(path: Path) -> Bundle:
    if path.is_dir():
        path = path / "onboarding.json"
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as archive:
            candidates = [
                n
                for n in archive.namelist()
                if n == "onboarding.json" or n.endswith("/onboarding.json")
            ]
            if len(candidates) != 1:
                raise ValueError("ZIP must contain exactly one onboarding.json")
            info = archive.getinfo(candidates[0])
            if info.file_size > MAX_SPEC_BYTES:
                raise ValueError(
                    f"onboarding.json is {info.file_size:,} bytes; the limit is {MAX_SPEC_BYTES:,}"
                )
            raw = archive.read(info)
    else:
        raw = path.read_bytes()
    bundle = Bundle.model_validate_json(raw)
    if bundle.suite_version != __version__:
        raise ValueError(f"Bundle needs suite {bundle.suite_version}; installed {__version__}")
    if bool(bundle.wheel_file) != bool(bundle.wheel_sha256):
        raise ValueError("Bundle wheel requires both a filename and checksum")
    if any(bundle.constraints.get(k) != v for k, v in bundle.versions.items()):
        raise ValueError("Runtime versions disagree with installation constraints")
    return bundle


def resolve_site(bundle: Bundle, local: Path) -> SiteConfig:
    values = yaml.safe_load(local.read_text())
    if not isinstance(values, dict):
        raise ValueError("site.yaml must contain a mapping")
    site = SiteConfig.model_validate(without_blanks({**bundle.defaults, **values}))
    check_scheduler_settings(bundle.scheduler, site)
    return site


def runtime_report(data_dir: str, work_dir: str, packages: list[str], state_dir: str) -> dict:
    # Self-contained so Globus can ship it as source to the worker.
    import hashlib
    import importlib.metadata
    import json
    import os
    import platform
    import pwd
    import tempfile
    from pathlib import Path

    try:
        raw_bundle = json.loads((Path(state_dir) / "onboarding.json").read_text())
        canonical = json.dumps(raw_bundle, sort_keys=True, separators=(",", ":"))
        bundle_digest = hashlib.sha256(canonical.encode()).hexdigest()
    except (OSError, ValueError):
        bundle_digest = "missing"
    versions = {}
    for name in packages:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "missing"
    data = Path(data_dir)
    readable = data.is_dir() and os.access(data, os.R_OK | os.X_OK)
    writable = False
    try:
        with tempfile.TemporaryFile(dir=work_dir):
            writable = True
    except OSError:
        pass
    return {
        "bundle_digest": bundle_digest,
        "python": platform.python_version(),
        "versions": versions,
        "user": pwd.getpwuid(os.geteuid()).pw_name,
        "data_readable": readable,
        "work_writable": writable,
    }


def check_runtime(
    bundle: Bundle, report: dict, *, check_user: bool = True, check_data: bool = True
) -> None:
    errors = []
    if report.get("bundle_digest") != digest(bundle):
        errors.append("Worker setup does not match this bundle revision")
    if ".".join(report.get("python", "").split(".")[:2]) != bundle.python_version:
        errors.append(f"Python must be {bundle.python_version}")
    expected = {**bundle.versions, "appfl-bio-suite": bundle.suite_version}
    for name, version in expected.items():
        actual = report.get("versions", {}).get(name, "missing").split("+", 1)[0]
        if actual != version:
            errors.append(f"{name}: expected {version}, found {actual}")
    if check_user and report.get("user") != bundle.service_account:
        errors.append(f"worker must run as {bundle.service_account}")
    if check_data and not report.get("data_readable"):
        errors.append("data_dir must exist and be readable by the service account")
    if not report.get("work_writable"):
        errors.append("work_dir must be writable by the service account")
    if errors:
        raise ValueError("\n".join(errors))


def run(args: list[str], **kwargs):
    return subprocess.run(args, check=True, **kwargs)


def inspect_runtime(python: Path, bundle: Bundle, site: SiteConfig) -> dict:
    import inspect

    script = (
        inspect.getsource(runtime_report)
        + "\nimport json\nprint(json.dumps(runtime_report("
        + repr(str(site.data_dir))
        + ","
        + repr(str(site.work_dir))
        + ","
        + repr([*bundle.versions, "appfl-bio-suite"])
        + ","
        + repr(str(python.parent.parent.parent))
        + ")))"
    )
    return json.loads(run([str(python), "-c", script], capture_output=True, text=True).stdout)


def _worker_settings(bundle: Bundle, site: SiteConfig, environment: Path) -> dict:
    worker_init = "\n".join(
        [
            f"source {shlex.quote(str(environment / 'bin/activate'))}",
            "export TMPDIR=/tmp",
            "export OPENBLAS_NUM_THREADS=1",
            "export OMP_NUM_THREADS=1",
            "export MKL_NUM_THREADS=1",
            "export PYTHONUNBUFFERED=1",
        ]
    )
    if bundle.scheduler == "local":
        # Workers on this host: no account, queue, walltime or launcher to choose.
        provider = {
            "type": "LocalProvider",
            "init_blocks": 0,
            "min_blocks": 0,
            "max_blocks": 1,
            "worker_init": worker_init,
        }
    else:
        slurm = bundle.scheduler == "slurm"
        provider = {
            "type": "SlurmProvider" if slurm else "PBSProProvider",
            "account": site.account,
            "partition" if slurm else "queue": site.partition if slurm else site.queue,
            "nodes_per_block": 1,
            "init_blocks": 0,
            "min_blocks": 0,
            "max_blocks": 1,
            "walltime": "00:30:00",
            "launcher": {"type": "SrunLauncher" if slurm else "MpiExecLauncher"},
            "worker_init": worker_init,
        }
    engine = {
        "type": "GlobusComputeEngine",
        "max_workers_per_node": 1,
        "strategy": "simple",
        "provider": provider,
    }
    if site.interface:
        engine["address"] = {"type": "address_by_interface", "ifname": site.interface}
    return {"engine": engine, "idle_heartbeats_soft": 0, "idle_heartbeats_hard": 5760}


def worker_config(bundle: Bundle, site: SiteConfig, environment: Path) -> str:
    """The MEP's per-user template: what each mapped user's endpoint runs."""
    return yaml.safe_dump(_worker_settings(bundle, site, environment), sort_keys=False)


def endpoint_config(bundle: Bundle, site: SiteConfig, environment: Path) -> str:
    """A single-user endpoint's whole config.yaml: the same engine, no mapping."""
    settings = {
        "display_name": f"{bundle.endpoint_name} ({bundle.experiment} / {bundle.client_id})",
        **_worker_settings(bundle, site, environment),
    }
    return yaml.safe_dump(settings, sort_keys=False)


def prepared_files(bundle: Bundle) -> tuple[str, ...]:
    """What a finished `setup` leaves in the state directory for this endpoint mode."""
    if bundle.endpoint_mode == "single-user":
        return ("runtime/bin/python", "endpoint_config.yaml")
    return ("runtime/bin/python", "identity_mapping.json", "user_config_template.yaml.j2")


def setup_site(
    bundle: Bundle,
    site: SiteConfig,
    destination: Path,
    *,
    python: str = sys.executable,
    install_source: str | None = None,
    wheel: bytes | None = None,
) -> Path:
    """Own only a new state directory and its venv; never overwrite another install."""
    destination = destination.resolve()
    if any(c in str(destination) for c in "\r\n{}"):
        raise ValueError("State path cannot contain newlines or template braces")
    marker = destination / "state.json"
    if bundle.wheel_file and not install_source:
        if wheel is None or hashlib.sha256(wheel).hexdigest() != bundle.wheel_sha256:
            raise ValueError("Bundled wheel is missing or its checksum does not match")
    identity = {
        "bundle_digest": digest(bundle),
        "site": site.model_dump(mode="json"),
        "install_source": install_source or bundle.wheel_sha256 or bundle.install_spec,
    }
    if marker.exists():
        if json.loads(marker.read_text()) != identity:
            raise ValueError(
                "Setup directory belongs to different inputs; use a new --state directory"
            )
    elif destination.exists() and any(destination.iterdir()):
        raise ValueError("Setup directory is not empty; choose a new --state directory")
    version = run(
        [python, "-c", "import sys; print('%s.%s' % sys.version_info[:2])"],
        capture_output=True,
        text=True,
    ).stdout.strip()
    if version != bundle.python_version:
        raise ValueError(
            f"Use --python pointing to Python {bundle.python_version}; found {version}"
        )
    destination.mkdir(parents=True, exist_ok=True)
    write_json(marker, identity)
    write_json(destination / "onboarding.json", bundle.model_dump())
    environment = destination / "runtime"
    if not (environment / "bin/python").exists():
        run([python, "-m", "venv", str(environment)])
    requirement = install_source or bundle.install_spec
    if bundle.wheel_file and not install_source:
        wheel_path = destination / bundle.wheel_file
        wheel_path.write_bytes(wheel)
        requirement = str(wheel_path)
        if bundle.extras:
            requirement += "[" + ",".join(bundle.extras) + "]"
    constraints = destination / "constraints.txt"
    constraints.write_text(
        "".join(f"{name}=={version}\n" for name, version in bundle.constraints.items())
    )
    run(
        [
            str(environment / "bin/python"),
            "-m",
            "pip",
            "install",
            "-c",
            str(constraints),
            "--extra-index-url",
            "https://download.pytorch.org/whl/cpu",
            requirement,
        ]
    )
    run([str(environment / "bin/python"), "-m", "pip", "check"])
    site.work_dir.mkdir(parents=True, exist_ok=True)
    # Data may arrive after the software: it is checked at activation and verification,
    # which is when a worker would first read it.
    check_runtime(
        bundle,
        inspect_runtime(environment / "bin/python", bundle, site),
        check_user=False,
        check_data=False,
    )
    if bundle.endpoint_mode == "single-user":
        # The endpoint runs as the service account itself: there is no identity to map.
        (destination / "endpoint_config.yaml").write_text(
            endpoint_config(bundle, site, environment)
        )
        return destination
    from appfl_bio_suite.core.identity import build_mapping_document

    write_json(
        destination / "identity_mapping.json",
        build_mapping_document(
            identity_match=bundle.coordinator_identity_match, output_account=bundle.service_account
        ),
    )
    from appfl_bio_suite.core.identity import validate_mapping_file

    mapping_report = validate_mapping_file(
        destination / "identity_mapping.json",
        identity=bundle.coordinator_identity,
        expected_output=bundle.service_account,
    )
    if not mapping_report.passed:
        raise ValueError(mapping_report.render())
    (destination / "user_config_template.yaml.j2").write_text(
        worker_config(bundle, site, environment)
    )
    return destination


def load_state(state: Path):
    state = state.resolve()
    saved = json.loads((state / "state.json").read_text())
    bundle = read_bundle(state / "onboarding.json")
    if saved["bundle_digest"] != digest(bundle):
        raise ValueError("State bundle has changed; repeat setup in a new directory")
    return bundle, SiteConfig.model_validate(saved["site"])


def activate_site(state: Path) -> Path:
    """Configure/start a MEP as its administrator, preserving its generated mapping path."""
    import pwd

    state = state.resolve()
    bundle, site = load_state(state)
    if bundle.endpoint_mode != "multi-user":
        raise ValueError(
            "This bundle declares a single-user endpoint. As "
            f"{bundle.service_account}, run: site activate --state {state} "
            "--credentials-file FILE"
        )
    if os.geteuid() != 0:
        raise ValueError(
            "Multi-user endpoint activation requires root; ask the cluster administrator"
        )
    pwd.getpwnam(bundle.service_account)
    # Verify as the account that will execute jobs, not as root.
    python = state / "runtime/bin/python"
    run(
        [
            "runuser",
            "-u",
            bundle.service_account,
            "--",
            str(python),
            "-m",
            "appfl_bio_suite.cli",
            "site",
            "verify",
            "--state",
            str(state),
            "--runtime-only",
        ]
    )
    endpoint = Path.home() / ".globus_compute" / bundle.endpoint_name
    owner = endpoint / "suite-onboarding.json"
    if endpoint.exists():
        if not owner.exists() or json.loads(owner.read_text()) != {"state": str(state)}:
            raise ValueError(f"Refusing to modify unmanaged endpoint {endpoint}")
    else:
        command = str(state / "runtime/bin/globus-compute-endpoint")
        run([command, "configure", "--multi-user", "true", bundle.endpoint_name])
        write_json(owner, {"state": str(state)})
    config = yaml.safe_load((endpoint / "config.yaml").read_text())
    mapping = Path(config["identity_mapping_config_path"])
    if not mapping.is_absolute():
        mapping = endpoint / mapping
    if not mapping.resolve().is_relative_to(endpoint.resolve()):
        raise ValueError("Generated identity mapping must be inside the endpoint directory")
    mapping.write_text((state / "identity_mapping.json").read_text())
    (endpoint / "user_config_template.yaml.j2").write_text(
        (state / "user_config_template.yaml.j2").read_text()
    )
    run([str(state / "runtime/bin/globus-compute-endpoint"), "start", bundle.endpoint_name])
    return endpoint


def read_credentials_file(path: Path) -> dict[str, str]:
    """The service-account client ID and secret, from a file only its owner can read.

    ``KEY=VALUE`` lines (the format systemd's EnvironmentFile reads), for exactly the two
    Compute SDK variables. Nothing from the file is echoed in an error: a malformed line
    may be the secret itself.
    """
    info = path.stat()
    if info.st_uid != os.geteuid():
        raise ValueError(f"{path} must be owned by the account that runs the endpoint")
    if info.st_mode & 0o077:
        raise ValueError(f"{path} must be readable only by its owner (chmod 600 {path})")
    values: dict[str, str] = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        key, separator, value = line.removeprefix("export ").partition("=")
        key, value = key.strip(), value.strip()
        if not separator or key not in CREDENTIAL_VARS:
            raise ValueError(
                f"{path} line {number}: only {' and '.join(CREDENTIAL_VARS)} are accepted"
            )
        if len(value) >= 2 and value[0] == value[-1] and value[0] in "'\"":
            value = value[1:-1]
        values[key] = value
    return values


def systemd_unit(bundle: Bundle, state: Path, credentials_file: Path, user_dir: Path) -> str:
    """A service that keeps the endpoint running across logouts, crashes and reboots."""
    lines = [
        "[Unit]",
        f"Description=Globus Compute endpoint {bundle.endpoint_name} "
        f"({bundle.experiment} / {bundle.client_id})",
        "Wants=network-online.target",
        "After=network-online.target",
        "",
        "[Service]",
        "Type=simple",
        f"User={bundle.service_account}",
        f"EnvironmentFile={credentials_file}",
        f"Environment=GLOBUS_COMPUTE_USER_DIR={user_dir}",
        f"ExecStart={state / 'runtime/bin/globus-compute-endpoint'} start {bundle.endpoint_name}",
        "Restart=on-failure",
        "RestartSec=30",
        "",
        "[Install]",
        "WantedBy=multi-user.target",
    ]
    return "\n".join(lines) + "\n"


def activate_single_user(
    state: Path, credentials_file: Path | None = None, *, wait: float = 120.0
) -> dict:
    """Configure and start a single-user endpoint as the service account itself.

    No root and no identity mapping: tasks run as whoever starts the endpoint, so this
    refuses to run as anyone but the bundle's service account. The endpoint registers
    under the service-account client in GLOBUS_COMPUTE_CLIENT_ID/SECRET; without them the
    SDK would fall back to an interactive personal login, producing an endpoint only that
    person can submit to.
    """
    import pwd
    import time

    state = state.resolve()
    bundle, site = load_state(state)
    if bundle.endpoint_mode != "single-user":
        raise ValueError("This bundle declares a multi-user endpoint; activate it as root")
    if os.geteuid() == 0:
        raise ValueError(
            f"Run single-user activation as {bundle.service_account}, not root: the endpoint "
            "runs every task as the account that starts it"
        )
    user = pwd.getpwuid(os.geteuid()).pw_name
    if user != bundle.service_account:
        raise ValueError(
            f"Run as {bundle.service_account} (this shell is {user}): the endpoint runs "
            "every task as the account that starts it"
        )
    env = dict(os.environ)
    if credentials_file is not None:
        credentials_file = credentials_file.resolve()
        env.update(read_credentials_file(credentials_file))
    missing = [name for name in CREDENTIAL_VARS if not env.get(name)]
    if missing:
        raise ValueError(
            f"Set {' and '.join(missing)} for the site's Globus service-account client, or "
            "pass --credentials-file. Without them the endpoint would register under a "
            "personal Globus login, which the coordinator cannot submit to."
        )
    try:
        UUID(env["GLOBUS_COMPUTE_CLIENT_ID"])
    except ValueError as exc:
        raise ValueError("GLOBUS_COMPUTE_CLIENT_ID must be a client UUID") from exc

    # The data and work directories, as the account that will read and write them.
    verify_site(state)

    user_dir = Path(env.get("GLOBUS_COMPUTE_USER_DIR") or Path.home() / ".globus_compute")
    env["GLOBUS_COMPUTE_USER_DIR"] = str(user_dir)
    endpoint = user_dir / bundle.endpoint_name
    owner = endpoint / "suite-onboarding.json"
    command = str(state / "runtime/bin/globus-compute-endpoint")
    if endpoint.exists():
        if not owner.exists() or json.loads(owner.read_text()) != {"state": str(state)}:
            raise ValueError(f"Refusing to modify unmanaged endpoint {endpoint}")
        # Restarting an endpoint this state already manages; `stop` fails harmlessly if
        # it is not running.
        subprocess.run([command, "stop", bundle.endpoint_name], env=env, check=False)
    else:
        run([command, "configure", "--multi-user", "false", bundle.endpoint_name], env=env)
        write_json(owner, {"state": str(state)})
    (endpoint / "config.yaml").write_text((state / "endpoint_config.yaml").read_text())
    unit = None
    if credentials_file is not None:
        unit = state / f"{bundle.endpoint_name}.service"
        unit.write_text(systemd_unit(bundle, state, credentials_file, user_dir))
    run([command, "start", "--detach", bundle.endpoint_name], env=env)
    info = endpoint / "endpoint.json"
    deadline = time.monotonic() + wait
    while not info.exists() and time.monotonic() < deadline:
        time.sleep(1)
    if not info.exists():
        raise TimeoutError(
            f"The endpoint did not record its UUID in {info} within {wait:.0f}s; "
            f"see {endpoint / 'endpoint.log'}"
        )
    endpoint_id = str(UUID(json.loads(info.read_text())["endpoint_id"]))
    return {"endpoint": endpoint, "endpoint_id": endpoint_id, "unit": unit}


def verify_site(state: Path, endpoint_id: str | None = None) -> dict:
    state = state.resolve()
    bundle, site = load_state(state)
    report = inspect_runtime(state / "runtime/bin/python", bundle, site)
    check_runtime(bundle, report)
    if endpoint_id is None:
        return report
    receipt = Receipt(
        bundle_digest=digest(bundle),
        experiment=bundle.experiment,
        client_id=bundle.client_id,
        endpoint_id=UUID(endpoint_id),
        site=site,
        runtime=report,
        state_dir=state,
    )
    write_json(state / "registration.json", receipt.model_dump(mode="json"))
    return receipt.model_dump(mode="json")


def register_site(
    bundle: Bundle, receipt_path: Path, output: Path, timeout: float, credential=None
) -> dict:
    """Probe the site's real worker; ``credential`` is its --site-credentials route."""
    receipt = Receipt.model_validate_json(receipt_path.read_text())
    if (receipt.bundle_digest, receipt.experiment, receipt.client_id) != (
        digest(bundle),
        bundle.experiment,
        bundle.client_id,
    ):
        raise ValueError("Receipt does not match this site's bundle revision")
    check_runtime(bundle, receipt.runtime)
    if bundle.endpoint_mode == "single-user" and credential is None:
        raise ValueError(
            f"{bundle.client_id} uses a single-user endpoint, which accepts tasks only from "
            "its service-account client. Pass --site-credentials with a route for "
            f"endpoint {receipt.endpoint_id}."
        )
    if credential is not None and credential.expected_user != bundle.service_account:
        raise ValueError(
            f"The credential route expects user {credential.expected_user}, but this "
            f"site's worker account is {bundle.service_account}"
        )
    from globus_compute_sdk import Executor
    from globus_compute_sdk.serialize import CombinedCode, ComputeSerializer

    if credential is None:
        executor = Executor(endpoint_id=str(receipt.endpoint_id))
    else:
        executor = Executor(endpoint_id=str(receipt.endpoint_id), client=credential.client())
        # SDK 4.9 shares one task group across client identities; see site_credentials.
        executor.task_group_id = None
    executor.serializer = ComputeSerializer(strategy_code=CombinedCode())
    try:
        future = executor.submit(
            runtime_report,
            str(receipt.site.data_dir),
            str(receipt.site.work_dir),
            [*bundle.versions, "appfl-bio-suite"],
            str(receipt.state_dir),
        )
        report = future.result(timeout=timeout)
    finally:
        executor.shutdown(wait=False, cancel_futures=True)
    check_runtime(bundle, report)
    result = {**receipt.model_dump(mode="json"), "status": "ready", "runtime": report}
    write_json(output, result)
    return result


def updated_federation(bundle: Bundle, receipt: Receipt, source: Path) -> dict:
    """Prepare a validated coordinator config with the partner's actual UUID and paths."""
    from appfl_bio_suite.core.config import Federation, load_federation

    fed = load_federation(source)
    if experiment_digest(fed, bundle.experiment) != bundle.experiment_digest:
        raise ValueError("Experiment settings changed since packaging; generate new bundles")
    if fed.coordinator.identity != bundle.coordinator_identity:
        raise ValueError("Coordinator identity changed since packaging")
    fed.experiment_site(bundle.experiment, bundle.client_id)
    if fed.experiment(bundle.experiment).service_account != bundle.service_account:
        raise ValueError("Service account changed since packaging")
    document = yaml.safe_load(source.read_text())
    for item in document["experiments"][bundle.experiment]["sites"]:
        if item["client_id"] == bundle.client_id:
            item["endpoint_uuid"] = str(receipt.endpoint_id)
            item["data_dir"] = str(receipt.site.data_dir)
            item["output_dir"] = str(receipt.site.work_dir)
            break
    Federation.model_validate(document)
    return document


def attach_wheel(directory: Path, wheel: Path) -> None:
    """Bind a locally built release artifact to a site's specification and manifest."""
    from email.parser import BytesParser

    bundle = read_bundle(directory)
    with zipfile.ZipFile(wheel) as archive:
        names = [name for name in archive.namelist() if name.endswith(".dist-info/METADATA")]
        if len(names) != 1:
            raise ValueError("Wheel must contain exactly one distribution")
        metadata = BytesParser().parsebytes(archive.read(names[0]))
        if metadata["Name"] != "appfl-bio-suite" or metadata["Version"] != bundle.suite_version:
            raise ValueError("Wheel package/version does not match the bundle")
    raw = wheel.read_bytes()
    updated = Bundle.model_validate(
        {
            **bundle.model_dump(),
            "wheel_file": wheel.name,
            "wheel_sha256": hashlib.sha256(raw).hexdigest(),
        }
    )
    (directory / wheel.name).write_bytes(raw)
    write_json(directory / "onboarding.json", updated.model_dump())
    manifest_path = directory / "bundle_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["files"] = sorted(set([*manifest["files"], wheel.name]))
    write_json(manifest_path, manifest)


def read_wheel(path: Path, bundle: Bundle) -> bytes | None:
    if bundle.wheel_file is None:
        return None
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as archive:
            specs = [
                name
                for name in archive.namelist()
                if name == "onboarding.json" or name.endswith("/onboarding.json")
            ]
            prefix = specs[0][: -len("onboarding.json")]
            name = prefix + bundle.wheel_file
            if archive.namelist().count(name) != 1:
                raise ValueError("ZIP must contain exactly one matching wheel")
            if archive.getinfo(name).file_size > 100_000_000:
                raise ValueError("Suite wheel exceeds 100 MB")
            raw = archive.read(name)
    else:
        directory = path if path.is_dir() else path.parent
        raw = (directory / bundle.wheel_file).read_bytes()
    if hashlib.sha256(raw).hexdigest() != bundle.wheel_sha256:
        raise ValueError("Bundled wheel checksum does not match")
    return raw
