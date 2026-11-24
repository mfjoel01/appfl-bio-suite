"""Experimental: Prepare a small, synthetic two-site APPFL trial; no study data is required."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from uuid import uuid4

import yaml

from appfl_bio_suite.core.site_credentials import load_site_credentials

DATASET = """def get_dataset(seed=1):
    import torch
    from torch.utils.data import TensorDataset
    generator = torch.Generator().manual_seed(seed)
    x = torch.randn(16, 1, generator=generator)
    return TensorDataset(x, 2 * x + 1), None
"""

MODEL = """import torch

class TrialModel(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.linear = torch.nn.Linear(1, 1)

    def forward(self, x):
        return self.linear(x)
"""


METRIC = """def mean_squared_error(target, prediction):
    import numpy as np
    return float(np.mean((target - prediction) ** 2))
"""


def check_sites(manifest: Path, timeout: float = 300) -> int:
    """Check owner execution and both cross-site denials using harmless tasks."""
    import json

    from globus_compute_sdk import Executor

    from appfl_bio_suite.core.endpoint import smoke_test, where_am_i

    load_site_credentials(manifest, [])
    endpoints = list(json.loads(manifest.read_text())["endpoints"])
    if len(endpoints) != 2:
        raise ValueError("The trial requires exactly two endpoint routes.")
    credentials = load_site_credentials(manifest, endpoints)
    if any(credential.coordinator_login for credential in credentials.values()):
        raise ValueError(
            "The isolation trial compares two service-account clients; a coordinator_login "
            "route has no site identity to isolate."
        )
    clients = {endpoint: credentials[endpoint].client() for endpoint in endpoints}
    for endpoint in endpoints:
        result = smoke_test(endpoint, timeout=timeout, client=clients[endpoint])
        print(result.render(), flush=True)
        if not result.ok:
            return 1
        payload = result.payload or {}
        if payload.get("user") != credentials[endpoint].expected_user:
            print("FAIL: worker ran under an unexpected local account.", flush=True)
            return 1
        if not str(payload.get("python", "")).startswith("3.12."):
            print("FAIL: this trial requires Python 3.12 at both workers.", flush=True)
            return 1
    for caller, target in [endpoints, endpoints[::-1]]:
        executor = Executor(endpoint_id=target, client=clients[caller])
        executor.task_group_id = None
        try:
            executor.submit(where_am_i).result(timeout=timeout)
        except Exception as exc:
            if (
                getattr(exc, "http_status", None) == 403
                and getattr(exc, "code", None) == "ENDPOINT_ACCESS_FORBIDDEN"
            ):
                print(f"PASS: {caller} cannot submit to {target}.", flush=True)
            else:
                print(
                    f"INCONCLUSIVE: cross-site request returned {type(exc).__name__}; "
                    "expected ENDPOINT_ACCESS_FORBIDDEN.",
                    flush=True,
                )
                return 1
        else:
            print("FAIL: a site credential executed on the other site's endpoint.", flush=True)
            return 1
        finally:
            executor.shutdown(wait=False, cancel_futures=True)
    return 0


def prepare_trial(manifest: Path, out_dir: Path) -> tuple[Path, Path]:
    # Loading an empty selection validates the manifest without requiring secrets.
    # Read all endpoint IDs only after validating its structure.
    import json

    load_site_credentials(manifest, [])
    endpoints = list(json.loads(manifest.read_text())["endpoints"])
    if len(endpoints) != 2:
        raise ValueError("The synthetic trial requires exactly two endpoint routes.")
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=False)
    dataset = out_dir / "dataset.py"
    model = out_dir / "model.py"
    metric = out_dir / "metric.py"
    dataset.write_text(DATASET)
    model.write_text(MODEL)
    metric.write_text(METRIC)
    trial_id = uuid4().hex
    clients = []
    for index, endpoint in enumerate(endpoints, 1):
        clients.append(
            {
                "client_id": f"site-{index}",
                "endpoint_id": endpoint,
                "data_configs": {
                    "dataset_path": str(dataset),
                    "dataset_name": "get_dataset",
                    "dataset_kwargs": {"seed": index},
                },
                "train_configs": {
                    "logging_output_dirname": f"/tmp/appfl-site-trial-{trial_id}/site-{index}",
                    "logging_output_filename": "train",
                },
            }
        )
    server = {
        "client_configs": {
            "model_configs": {"model_path": str(model), "model_name": "TrialModel"},
            "train_configs": {
                "trainer": "VanillaTrainer",
                "mode": "step",
                "num_local_steps": 2,
                "optim": "SGD",
                "optim_args": {"lr": 0.01},
                "loss_fn": "MSELoss",
                "metric_path": str(metric),
                "metric_name": "mean_squared_error",
                "device": "cpu",
                "train_batch_size": 4,
                "train_data_shuffle": False,
                "do_validation": False,
                "do_pre_validation": False,
                "use_dp": False,
            },
            "comm_configs": {"compressor_configs": {"enable_compression": False}},
            "data_readiness_configs": {"generate_dr_report": False, "dr_metrics": {}},
        },
        "server_configs": {
            "num_clients": 2,
            "num_global_epochs": 2,
            "device": "cpu",
            "scheduler": "SyncScheduler",
            "scheduler_kwargs": {"same_init_model": True},
            "aggregator": "FedAvgAggregator",
            "aggregator_kwargs": {"client_weights_mode": "sample_size"},
            "logging_output_dirname": str(out_dir / "logs"),
            "logging_output_filename": "server",
        },
    }
    server_path, clients_path = out_dir / "server.yaml", out_dir / "clients.yaml"
    server_path.write_text(yaml.safe_dump(server, sort_keys=False))
    clients_path.write_text(yaml.safe_dump({"clients": clients}, sort_keys=False))
    return server_path, clients_path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--site-credentials",
        type=Path,
        default=os.environ.get("APPFL_BIO_SUITE_GLOBUS_CREDENTIALS"),
    )
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument(
        "--check-endpoints",
        action="store_true",
        help="Test owner execution, local accounts, and cross-site isolation.",
    )
    parser.add_argument("--timeout", type=float, default=300)
    args = parser.parse_args(argv)
    if not args.site_credentials:
        parser.error("Provide --site-credentials or APPFL_BIO_SUITE_GLOBUS_CREDENTIALS.")
    try:
        if args.check_endpoints:
            return check_sites(Path(args.site_credentials), args.timeout)
        if not args.out_dir:
            parser.error("Provide --out-dir to prepare a new trial.")
        server, clients = prepare_trial(Path(args.site_credentials), args.out_dir)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(f"Server config: {server}\nClient config: {clients}")
    print("Prepared synthetic data and model source for two sites and two FedAvg rounds.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
