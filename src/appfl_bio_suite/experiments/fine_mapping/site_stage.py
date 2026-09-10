"""The site stage as a command, not only as an APPFL trainer object.

WHY THIS EXISTS
---------------
``trainer.py`` computes a site's aggregates inside APPFL's ``ClientAgent``: constructed by
a framework, handed a config object, driven through ``train()``. That is right for the
Globus Compute path and useless to everything else. A TES executor runs a *command*; a
CWL descriptor wraps a *command*; a partner who wants to see what will run on their
cluster before agreeing to it wants to run a *command*.

So the same computation is reachable as::

    appfl-bio-suite site-stage fine-mapping --config site.json --data-dir <bundle> --out <dir>

and it is genuinely the same computation, not a reimplementation: this module constructs
``SiteFineMappingDataset`` and ``SiteFineMappingTrainer`` -- the very classes APPFL ships
to a worker -- and writes what they produce. There is no second copy of the statistics
here, which is the property that keeps the TES path and the Globus Compute path from
drifting into two different experiments.

WHAT IT WRITES
--------------
    <out>/<name>.npz            the aggregate payload: every tensor the trainer produced
    <out>/<name>.manifest.json  the trainer's manifest, unpacked from the payload
    <out>/ga4gh_provenance.json the DUO decision, the DRS verification, the TRS pin

The npz is the wire format for this path. APPFL's communicators serialize tensors for the
Globus Compute path; a TES task has no communicator, so the payload is written as arrays
and read back by the driver. Keys are the trainer's own ``<kind>::<scope>::<pop>`` strings
and are preserved exactly, so the aggregator's decoder needs no special case: it receives
the same dict either way.

THE GA4GH CHECKS RUN, AND THEY RUN FIRST
-----------------------------------------
Constructing the dataset is what enforces them (see ``dataset.py``), so a task whose data
use request is not permitted by this site's terms exits non-zero having opened no
genotype file. That is true of this path exactly as it is of the APPFL path, because it
is the same constructor -- which is the reason the check lives in the loader rather than
in either driver.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any

__all__ = ["run_site_stage", "load_payload", "main"]

log = logging.getLogger(__name__)

MANIFEST_SUFFIX = ".manifest.json"
PROVENANCE_FILENAME = "ga4gh_provenance.json"


def _logger(name: str) -> logging.Logger:
    """A logger shaped like the one APPFL hands a trainer.

    The trainer calls ``.info`` and ``.warning`` and nothing else, so this is the whole
    of the interface -- and building it here rather than importing APPFL's keeps this
    module runnable in a container that has the partner extra and no server config.
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(
            logging.Formatter("[%(asctime)s %(levelname)s] %(message)s", "%H:%M:%S")
        )
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger


def run_site_stage(
    config: dict[str, Any],
    data_dir: str | Path,
    out_dir: str | Path,
    out_name: str = "aggregates",
) -> dict[str, Path]:
    """Run one site's computation and write its outputs. Returns the paths written.

    ``config`` is the same shape as an APPFL client config's two blocks, flattened::

        {"client_id": "anl",
         "train_configs": {...},        locus shard, pops, uplink dtype, ga4gh_tool
         "dataset_kwargs": {...}}       data_use_request, drs_object, verify_bundles

    Deliberately the same shape rather than a new one: a coordinator debugging a TES task
    should be able to lift the block straight out of the generated client config.
    """
    import numpy as np

    from appfl_bio_suite.experiments.fine_mapping.dataset import get_dataset
    from appfl_bio_suite.experiments.fine_mapping.trainer import (
        KEY_MANIFEST,
        SiteFineMappingTrainer,
    )

    client_id = str(config.get("client_id") or config.get("site_id") or "site")
    train_configs = dict(config.get("train_configs", {}) or {})
    dataset_kwargs = dict(config.get("dataset_kwargs", {}) or {})
    dataset_kwargs.pop("data_dir", None)  # the command line is authoritative here
    dataset_kwargs.pop("site_id", None)

    out_dir = Path(out_dir).resolve()
    logger = _logger(f"site-stage.{client_id}")

    # Constructing the dataset is what runs the DUO and DRS checks. A refusal raises out
    # of here and the process exits non-zero, having read nothing.
    #
    # Before the output directory is created, deliberately: a refused task should leave
    # nothing behind at all, so that an empty directory in an outputs tree is never
    # something a coordinator has to interpret.
    train_dataset, _ = get_dataset(data_dir=data_dir, site_id=client_id, **dataset_kwargs)
    out_dir.mkdir(parents=True, exist_ok=True)

    # The trainer writes its logs beside the outputs unless told otherwise; inside a
    # container there is nowhere else sensible for them to go.
    train_configs.setdefault("trainer_output_dirname", str(out_dir))

    trainer = SiteFineMappingTrainer(
        train_dataset=train_dataset,
        train_configs=_as_config(train_configs),
        logger=logger,
        client_id=client_id,
    )
    trainer.train()
    payload = trainer.get_parameters()

    written: dict[str, Path] = {}

    manifest_blob = payload.get(KEY_MANIFEST)
    manifest = (
        json.loads(bytes(manifest_blob.numpy().tobytes()).decode("utf-8"))
        if manifest_blob is not None
        else {}
    )

    arrays = {
        key: value.detach().cpu().numpy() for key, value in payload.items() if key != KEY_MANIFEST
    }
    payload_path = out_dir / f"{out_name}.npz"
    # Uncompressed: the Gram is dense float64 with no exploitable structure, so
    # compression spends CPU at both ends for a few percent -- the same reasoning that
    # keeps the APPFL compressor off for this experiment.
    np.savez(payload_path, **arrays)
    written["payload"] = payload_path

    manifest_path = out_dir / f"{out_name}{MANIFEST_SUFFIX}"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    written["manifest"] = manifest_path

    provenance_path = out_dir / PROVENANCE_FILENAME
    provenance_path.write_text(
        json.dumps(
            {
                "client_id": client_id,
                "data_dir": str(Path(data_dir).resolve()),
                **(manifest.get("ga4gh") or {}),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    written["provenance"] = provenance_path

    total_mb = payload_path.stat().st_size / (1 << 20)
    logger.info(f"{client_id}: wrote {len(arrays)} array(s), {total_mb:.1f} MB -> {payload_path}")
    return written


def _as_config(mapping: dict[str, Any]):
    """Wrap a plain dict so the trainer's ``.get`` calls work unchanged.

    APPFL hands the trainer an OmegaConf node. A dict satisfies every access the trainer
    makes, and using one here avoids requiring omegaconf inside a container that only
    needs to compute second moments.
    """
    return mapping


def load_payload(path: str | Path) -> dict[str, Any]:
    """Read an ``.npz`` payload back into the dict the aggregator decodes.

    Returns torch tensors, because that is what the aggregator's ``_to_numpy`` and the
    APPFL path both hand it -- the decoder must not have to know which driver produced
    its input.
    """
    import numpy as np
    import torch

    path = Path(path)
    with np.load(path, allow_pickle=False) as archive:
        payload = {key: torch.from_numpy(archive[key]) for key in archive.files}

    manifest_path = path.with_name(path.name.replace(".npz", MANIFEST_SUFFIX))
    if manifest_path.is_file():
        from appfl_bio_suite.experiments.fine_mapping.trainer import KEY_MANIFEST

        blob = manifest_path.read_bytes()
        # Re-encoded exactly as the trainer sent it, so the aggregator's decoder is
        # byte-for-byte the same code path on both transports.
        payload[KEY_MANIFEST] = torch.frombuffer(bytearray(blob), dtype=torch.uint8)
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="appfl-bio-suite site-stage",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("experiment", nargs="?", default="fine-mapping")
    parser.add_argument("--config", required=True, help="JSON run config for this site.")
    parser.add_argument("--data-dir", required=True, help="This site's bundle directory.")
    parser.add_argument("--out", default=".", help="Where to write the payload.")
    parser.add_argument("--out-name", default="aggregates")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO, format="[%(asctime)s %(levelname)s] %(message)s", datefmt="%H:%M:%S"
    )

    if args.experiment != "fine-mapping":
        print(
            f"site-stage is implemented for fine-mapping; got '{args.experiment}'.\n"
            "The other experiments have no command-line site stage: FLamby's is a "
            "multi-round training loop with no single-shot form, and the GWAS site "
            "stage has not been given one.",
            file=sys.stderr,
        )
        return 2

    try:
        config = json.loads(Path(args.config).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        print(f"could not read --config {args.config}: {exc}", file=sys.stderr)
        return 2

    try:
        run_site_stage(config, args.data_dir, args.out, args.out_name)
    except PermissionError as exc:
        # The DUO refusal. Given its own exit code so that a TES engine's task log
        # distinguishes "this site declined" from "this task crashed" -- they call for
        # completely different responses from a coordinator.
        print(f"\nDATA USE REFUSED\n\n{exc}", file=sys.stderr)
        return 77
    except Exception as exc:  # noqa: BLE001 - a CLI boundary; the traceback is in the log
        log.exception("site stage failed")
        print(f"\nsite stage failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
