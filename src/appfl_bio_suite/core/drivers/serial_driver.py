"""Loopback driver: a real federation, on one machine, with no external dependencies.

WHY THIS EXISTS
---------------
Everything else in this suite needs a partner. Standing up a federation means recruiting
institutions, getting service accounts created, waiting on identity mappings. Before any
of that, someone needs to be able to answer: *is my installation correct?*

This runs the genuine APPFL round loop with two or more simulated sites in a single
process. No Globus Compute, no scheduler, no network, no credentials. If it completes,
the install works: the configs are valid, the shipped modules import, the data loads, the
trainer runs, and the aggregator produces a global model.

It is also the only federated path CI can exercise, and it is what a stranger should run
first.

BUILT ON APPFL'S OWN RUNNER, NOT A MOCK
---------------------------------------
This drives ``ClientAgent`` and ``ServerAgent`` directly -- the same objects the Globus
Compute path drives remotely. A mock transport would prove only that the mock works. The
one thing it cannot exercise is the transport itself: serialization, identity mapping,
and scheduler behaviour are exactly what a loopback run does not test, which is why
``endpoint smoke`` exists separately.

WHAT DIFFERS FROM THE REMOTE PATH
---------------------------------
``dataset_path`` is read directly here rather than being inlined and shipped, so a
loopback run does *not* catch a shipped module that imports something a partner will not
have. ``tests/test_shipped_modules.py`` covers that instead.

Only synchronous schedulers can be simulated serially -- an asynchronous aggregator has
nothing to be asynchronous about when clients run one after another.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import appfl_bio_suite  # noqa: F401  -- applies the compat shim before APPFL is imported


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--server-config", required=True)
    parser.add_argument("--client-config", required=True)
    args = parser.parse_args(argv)

    from appfl.agent import ClientAgent, ServerAgent
    from omegaconf import OmegaConf

    server_config = OmegaConf.load(args.server_config)
    client_configs = list(OmegaConf.load(args.client_config)["clients"])

    if not client_configs:
        print("no clients in the client config", file=sys.stderr)
        return 1

    server_config.server_configs.num_clients = len(client_configs)

    scheduler = server_config.server_configs.get("scheduler", "SyncScheduler")
    if "Async" in str(scheduler):
        print(
            f"scheduler is {scheduler}, but a serial run executes clients one after\n"
            "another -- there is nothing for an asynchronous scheduler to overlap.\n"
            "Use a synchronous server config for loopback runs.",
            file=sys.stderr,
        )
        return 1

    server_agent = ServerAgent(server_agent_config=server_config)
    log = server_agent.logger

    rounds = server_config.server_configs.num_global_epochs
    log.info(f"[loopback] {len(client_configs)} simulated site(s), {rounds} round(s)")
    log.info("[loopback] no Globus Compute, no scheduler, no network")

    # Merged BEFORE the agents are built, mirroring what the Globus Compute communicator
    # does before it ships a config to a worker:
    #
    #     client_config=OmegaConf.merge(server_agent_config.client_configs, client_config)
    #
    # It has to happen first. `ClientAgent.__init__` calls `_load_trainer()` immediately,
    # and `trainer_path` lives in the server's shared block -- so constructing the agent
    # from the per-client config alone and merging afterwards fails on any experiment
    # whose trainer is loaded from a file rather than named among APPFL's built-ins
    # ("Invalid trainer name: SiteGWASTrainer"). Per-client values still win, as they do
    # remotely: the shared block is the base.
    shared = server_agent.get_client_configs()
    agents = []
    for config in client_configs:
        # Endpoint ids are meaningless here and their presence in a loopback config is
        # confusing, so they are dropped rather than ignored.
        config = OmegaConf.create({k: v for k, v in config.items() if k != "endpoint_id"})
        agents.append(ClientAgent(client_agent_config=OmegaConf.merge(shared, config)))

    initial = server_agent.get_parameters(serial_run=True)
    if isinstance(initial, tuple):
        initial = initial[0]
    for agent in agents:
        agent.load_parameters(initial)

    for agent in agents:
        try:
            size = agent.get_sample_size()
        except Exception:  # noqa: BLE001 - not every loader reports one
            continue
        server_agent.set_sample_size(client_id=agent.get_id(), sample_size=size)
        log.info(f"[loopback] {agent.get_id()}: n = {size}")

    round_no = 0
    while not server_agent.training_finished():
        round_no += 1
        log.info(f"[loopback] round {round_no}/{rounds}")
        futures = []
        for agent in agents:
            agent.train()
            local = agent.get_parameters()
            if isinstance(local, tuple):
                local, metadata = local
            else:
                metadata = {}
            futures.append(
                server_agent.global_update(
                    client_id=agent.get_id(), local_model=local, blocking=False, **metadata
                )
            )
        for agent, future in zip(agents, futures, strict=True):
            agent.load_parameters(future.result())

    log.info("Federated Learning Training Completed!")
    out = server_config.server_configs.get("logging_output_dirname")
    if out:
        log.info(f"[loopback] results in {Path(out).resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
