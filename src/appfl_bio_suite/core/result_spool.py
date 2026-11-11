"""Return site results larger than Globus Compute's result cap, over Globus Compute.

A Globus Compute endpoint refuses to return a task result larger than 10 MiB
(``MaxResultSizeExceeded``), and the limit is not configurable. A fine-mapping site's
payload is O(M^2) per locus and ancestry: tens of MiB for one locus, gigabytes for a
sweep. Sharding by locus cannot get under the cap.

So the site writes its payload to its own work directory and returns a small *handle*
instead. The coordinator's driver sees the handle, reads the file back in chunks that
each fit under the cap -- through the same endpoint, with the same credentials, so no
new service, transfer collection or trust relationship is involved -- checks the
sha-256 the site reported, and hands the reassembled payload to the aggregator exactly
as if it had arrived inline. The site's copy is removed once it has been verified.

The handle travels as one uint8 tensor under :data:`SPOOL_KEY`, the same way the
trainer's manifest does, so it passes through APPFL untouched.
"""

from __future__ import annotations

import hashlib
import json
from collections import deque
from collections.abc import Callable
from pathlib import Path
from typing import Any

# Must equal trainer.KEY_SPOOL; the trainer is shipped as source and cannot import this.
SPOOL_KEY = "result_spool"
SPOOL_FORMAT = "torch-save-v1"

# 6 MiB of raw bytes serializes to about 8.5 MB through the SDK's base64 data
# serializer, comfortably under the endpoint's 10 MiB (10,485,760 byte) limit.
CHUNK_BYTES = 6 * 1024 * 1024

# A site spools when its payload's tensor bytes exceed this. Below it the payload goes
# inline as before, so small runs are byte-for-byte unchanged.
INLINE_LIMIT_BYTES = 6 * 1024 * 1024

_FIELDS = {"format", "path", "size", "sha256", "client_id"}


def spool_info(payload: Any) -> dict | None:
    """The handle in ``payload``, validated; ``None`` when the payload is inline."""
    if not isinstance(payload, dict) or SPOOL_KEY not in payload:
        return None
    blob = payload[SPOOL_KEY]
    raw = blob.numpy().tobytes() if hasattr(blob, "numpy") else bytes(blob)
    info = json.loads(raw.decode("utf-8"))
    if (
        not isinstance(info, dict)
        or set(info) != _FIELDS
        or info["format"] != SPOOL_FORMAT
        or not str(info["path"]).startswith("/")
        or not isinstance(info["size"], int)
        or info["size"] <= 0
        or len(str(info["sha256"])) != 64
    ):
        raise ValueError(f"malformed result-spool handle: {sorted(info)}")
    return info


def read_spool_chunk(path: str, offset: int, length: int) -> bytes:
    """Runs at the site. Self-contained, so Globus can ship it as source."""
    from pathlib import Path

    spool = Path(path)
    if spool.suffix != ".pt" or spool.parent.name != "spool":
        raise ValueError(f"{path} is not a result-spool file")
    with spool.open("rb") as handle:
        handle.seek(offset)
        return handle.read(length)


def remove_spool(path: str) -> bool:
    """Runs at the site, after the coordinator has verified its copy."""
    from pathlib import Path

    spool = Path(path)
    if spool.suffix != ".pt" or spool.parent.name != "spool":
        raise ValueError(f"{path} is not a result-spool file")
    spool.unlink(missing_ok=True)
    return True


def spool_executor(endpoint_id: str, credential=None):
    """An executor for one site's chunk reads, separate from APPFL's own.

    Separate so its function serialization (source, via CombinedCode) cannot change how
    APPFL ships its tasks, and so the reads work on workers without this package --
    the functions above import nothing from it.
    """
    from globus_compute_sdk import Executor
    from globus_compute_sdk.serialize import CombinedCode, ComputeSerializer

    if credential is None:
        executor = Executor(endpoint_id=endpoint_id)
    else:
        executor = Executor(endpoint_id=endpoint_id, client=credential.client())
        # SDK 4.9 shares one task group across client identities; see site_credentials.
        executor.task_group_id = None
    executor.serializer = ComputeSerializer(strategy_code=CombinedCode())
    return executor


def fetch_spooled(
    info: dict,
    submit: Callable,
    dest_dir: Path,
    *,
    chunk_bytes: int = CHUNK_BYTES,
    window: int = 8,
    retries: int = 3,
    timeout: float = 900.0,
    log: Callable[[str], None] = print,
) -> dict:
    """Pull a spooled payload back through ``submit`` and return it, verified.

    ``submit(fn, *args)`` returns a future for ``fn(*args)`` run at the site. Up to
    ``window`` reads are in flight; they are written in order, and a failed read is
    retried before the fetch gives up.
    """
    import torch

    size, path = info["size"], info["path"]
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{info['client_id']}-{info['sha256'][:16]}.pt"
    offsets = list(range(0, size, chunk_bytes))
    digest = hashlib.sha256()
    pending: deque = deque()
    queued = written = 0
    reported = -1
    try:
        with dest.open("wb") as out:
            while queued < len(offsets) or pending:
                while queued < len(offsets) and len(pending) < window:
                    offset = offsets[queued]
                    length = min(chunk_bytes, size - offset)
                    pending.append(
                        (offset, length, submit(read_spool_chunk, path, offset, length), 1)
                    )
                    queued += 1
                offset, length, future, attempt = pending.popleft()
                try:
                    data = future.result(timeout=timeout)
                except Exception as exc:  # noqa: BLE001 - a transient task failure is retried
                    if attempt >= retries:
                        raise RuntimeError(
                            f"{info['client_id']}: reading bytes {offset}-{offset + length} of "
                            f"its spooled result failed {attempt} time(s): "
                            f"{type(exc).__name__}: {exc}"
                        ) from exc
                    retry = submit(read_spool_chunk, path, offset, length)
                    pending.appendleft((offset, length, retry, attempt + 1))
                    continue
                if len(data) != length:
                    raise RuntimeError(
                        f"{info['client_id']}: expected {length} bytes at offset {offset}, "
                        f"got {len(data)}; the spool file changed or was truncated"
                    )
                out.write(data)
                digest.update(data)
                written += length
                percent = written * 100 // size
                if percent // 10 > reported:
                    reported = percent // 10
                    done, total = written / 2**20, size / 2**20
                    log(f"{info['client_id']}: fetched {done:.0f} of {total:.0f} MiB")
        if digest.hexdigest() != info["sha256"]:
            raise RuntimeError(
                f"{info['client_id']}: reassembled result does not match the sha-256 the site "
                "reported; refusing to aggregate it"
            )
        return torch.load(dest, map_location="cpu", weights_only=True)
    finally:
        dest.unlink(missing_ok=True)
