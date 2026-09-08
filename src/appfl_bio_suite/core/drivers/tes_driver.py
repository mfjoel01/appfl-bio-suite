"""TES driver: dispatch the site stage as GA4GH tasks, then aggregate here.

HOW THIS DIFFERS FROM THE OTHER TWO DRIVERS
-------------------------------------------
The Globus Compute and serial drivers run APPFL's round loop: a ``ServerAgent`` hands out
a global model, clients train, the aggregator combines. This one does not, and the reason
is a property of the experiment rather than of TES.

Fine-mapping is a **single-round aggregate** exchange. There is no global model to
distribute, no convergence to reach, and no second round in which a site would need
anything back. So the honest shape for a TES run is: submit one task per site, wait, read
the payloads back, and call the same ``FineMappingAggregator`` the other drivers call. A
round loop wrapped around that would be ceremony -- and worse, it would imply this driver
could serve a multi-round experiment, which it cannot.

WHAT A TES RUN REQUIRES THAT THE OTHERS DO NOT
-----------------------------------------------
1. **A published container image.** The TRS pin names it. Until it exists there is
   nothing for an executor to run, and this driver says so before submitting anything.
2. **Outputs the driver can read.** TES stages a task's outputs to a URL the *service*
   can write. This driver then has to read them, so that URL must resolve here too --
   a shared filesystem path, or an object store both ends can reach. ``file://`` URLs are
   read directly; anything else must already be mounted or synced, and the driver says
   which path it could not find rather than hanging.

Both are properties of running somebody else's container on somebody else's cluster, and
neither can be papered over from this side.
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import appfl_bio_suite  # noqa: F401  -- applies the compat shim before APPFL is imported


def _local_path(url: str) -> Path | None:
    """The local path a ``file://`` (or bare path) output URL names, if it is one."""
    parsed = urlparse(url)
    if parsed.scheme in ("", "file"):
        return Path(parsed.path or url)
    return None


def _site_run_config(client: dict[str, Any], shared_train: dict[str, Any]) -> dict[str, Any]:
    """The JSON the container reads, built from the generated client config.

    Merged the same way APPFL's communicator merges before shipping to a worker: the
    server's shared ``train_configs`` are the base and the per-client block wins. Doing
    it differently here would mean a TES site computing over a different locus shard than
    a Globus Compute site in the same federation.
    """
    train = dict(shared_train)
    train.update(dict(client.get("train_configs", {}) or {}))
    # Worker-side paths from the Globus path mean nothing in a container.
    for key in ("logging_output_dirname", "trainer_output_dirname", "logging_output_filename"):
        train.pop(key, None)

    data_kwargs = dict((client.get("data_configs", {}) or {}).get("dataset_kwargs", {}) or {})
    data_kwargs.pop("data_dir", None)  # the task stages the bundle at a fixed path
    return {
        "client_id": client.get("client_id"),
        "train_configs": train,
        "dataset_kwargs": data_kwargs,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--server-config", required=True)
    parser.add_argument("--client-config", required=True)
    parser.add_argument("--watch", default=None)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write the task documents and stop. What a partner asks to see first.",
    )
    args = parser.parse_args(argv)

    from omegaconf import OmegaConf

    from appfl_bio_suite.core.ga4gh.tes import TesClient, build_site_task, write_task
    from appfl_bio_suite.core.watch import RunWatcher

    server_config = OmegaConf.load(args.server_config)
    client_configs = [
        OmegaConf.to_container(c, resolve=True)
        for c in OmegaConf.load(args.client_config)["clients"]
    ]
    server = OmegaConf.to_container(server_config, resolve=True)
    aggregator_kwargs = (server.get("server_configs") or {}).get("aggregator_kwargs") or {}
    shared_train = ((server.get("client_configs") or {}).get("train_configs") or {})
    ga4gh = aggregator_kwargs.get("ga4gh") or {}
    tes_settings = ga4gh.get("tes") or {}

    tool = ga4gh.get("tool") or {}
    image = tool.get("image")
    if not image:
        print(
            "no container image for the site stage.\n"
            "The TES path runs a container, and the image is named by the TRS tool pin:\n"
            "  experiments.<name>.ga4gh.tool.image in federation.yaml\n"
            "Generate the Containerfile with `appfl-bio-suite ga4gh trs publish`, build "
            "and push it, then pin it (ideally by digest).",
            file=sys.stderr,
        )
        return 1

    out_root = Path(aggregator_kwargs.get("output_dir", "local/output/fine-mapping")).resolve()
    tasks_dir = out_root / "tes"
    tasks_dir.mkdir(parents=True, exist_ok=True)

    # -- build one task per site ------------------------------------------
    plans = []
    for client in client_configs:
        client_id = str(client.get("client_id"))
        site = (tes_settings.get("sites") or {}).get(client_id, {})
        service_url = site.get("url") or tes_settings.get("url")
        bundle_url = site.get("drs_uri") or (ga4gh.get("drs") or {}).get("objects", {}).get(
            client_id
        )
        outputs_url = site.get("outputs_url")

        problems = []
        if not service_url:
            problems.append("no TES service URL (ga4gh.tes.url, or the site's tes_url)")
        if not bundle_url:
            problems.append("no DRS URI for its bundle (the site's drs_uri)")
        if not outputs_url:
            problems.append("no outputs URL (the site's tes_outputs_url)")
        if problems:
            print(f"{client_id}: cannot build a task -- " + "; ".join(problems), file=sys.stderr)
            return 1

        task = build_site_task(
            client_id=client_id,
            experiment="fine-mapping",
            image=image,
            bundle_url=bundle_url,
            outputs_url=outputs_url,
            run_config=_site_run_config(client, shared_train),
            cpu_cores=int(tes_settings.get("cpu_cores", 4)),
            ram_gb=float(tes_settings.get("ram_gb", 16.0)),
            disk_gb=float(tes_settings.get("disk_gb", 64.0)),
            preemptible=bool(tes_settings.get("preemptible", False)),
            tags={
                "trs_id": str(tool.get("id", "")),
                "trs_version": str(tool.get("version", "")),
                "drs_uri": str(bundle_url),
            },
        )
        path = write_task(task, tasks_dir / f"{client_id}.task.json")
        print(f"{client_id}: task -> {path}", file=sys.stderr)
        plans.append((client_id, service_url, outputs_url, task))

    if args.dry_run:
        print("\n--dry-run: task documents written, nothing submitted.", file=sys.stderr)
        return 0

    watcher = RunWatcher.load(args.watch)
    watcher.start(
        algorithm="single-round-aggregate",
        config={"rounds": 1, "sites": len(plans), "driver": "tes"},
    )
    watcher.round_start(1)

    # -- submit, and wait ---------------------------------------------------
    poll = float(tes_settings.get("poll_seconds", 15.0))
    timeout = float(tes_settings.get("timeout_seconds", 7200.0))
    submitted: dict[str, tuple[TesClient, str]] = {}

    def run_one(plan):
        client_id, service_url, outputs_url, task = plan
        client = TesClient(service_url)
        task_id = client.create_task(task)
        submitted[client_id] = (client, task_id)
        print(f"{client_id}: submitted {task_id} to {service_url}", file=sys.stderr)

        def on_state(_id, state):
            print(f"{client_id}: {state}", file=sys.stderr)

        finished = client.wait(
            task_id, poll_seconds=poll, timeout_seconds=timeout, on_state=on_state
        )
        return client_id, outputs_url, finished

    # Concurrently, because sites are independent and a serial wait would make the run as
    # long as the sum of the sites rather than the slowest one. Cancellation on failure
    # is handled below rather than here: a site that fails should not take down the
    # others' tasks silently.
    results = []
    with ThreadPoolExecutor(max_workers=max(1, len(plans))) as pool:
        for outcome in pool.map(run_one, plans):
            results.append(outcome)

    failed = [(cid, task) for cid, _, task in results if str(task.state) != "COMPLETE"]
    if failed:
        for client_id, task in failed:
            print(f"\n{client_id}: {task.state}", file=sys.stderr)
            for entry in task.logs or []:
                for log_entry in entry.get("logs", []) or []:
                    if log_entry.get("stderr"):
                        print(f"  stderr: {log_entry['stderr'][-2000:]}", file=sys.stderr)
                    if log_entry.get("exit_code") == 77:
                        print(
                            "  exit 77 is the site stage's DATA USE REFUSED code: this "
                            "site's DUO terms do not permit the declared study.",
                            file=sys.stderr,
                        )
        watcher.finish()
        return 1

    # -- read the payloads back --------------------------------------------
    from appfl_bio_suite.experiments.fine_mapping.site_stage import load_payload

    local_models: dict[str, Any] = {}
    for client_id, outputs_url, _task in results:
        base = _local_path(outputs_url)
        if base is None:
            print(
                f"{client_id}: outputs are at {outputs_url}, which this driver cannot "
                "read directly. Stage them to a path this machine can see, or fetch them "
                "and re-run with --driver serial over the downloaded payloads.",
                file=sys.stderr,
            )
            return 1
        payload_path = base / "aggregates.npz"
        if not payload_path.is_file():
            print(
                f"{client_id}: the task completed but {payload_path} is not there. The "
                "service staged its outputs somewhere this driver cannot see.",
                file=sys.stderr,
            )
            return 1
        local_models[client_id] = load_payload(payload_path)
        watcher.client_update(client_id, 1, num_samples=None)

    # -- aggregate, in this process ----------------------------------------
    import logging

    from appfl_bio_suite.experiments.fine_mapping.aggregator import FineMappingAggregator

    logging.basicConfig(
        level=logging.INFO, format="[%(asctime)s %(levelname)s] %(message)s", datefmt="%H:%M:%S"
    )
    aggregator = FineMappingAggregator(
        aggregator_configs=aggregator_kwargs, logger=logging.getLogger("fine-mapping")
    )
    aggregator.aggregate(local_models)

    watcher.round_end(1)
    watcher.finish()
    print(f"\nresults in {out_root}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
