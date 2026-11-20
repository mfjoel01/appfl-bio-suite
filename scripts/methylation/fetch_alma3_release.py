#!/usr/bin/env python3
"""Extract the immutable v3.1.2 release from GHCR without Docker/Apptainer.

Downloads only the release layer, verifies its OCI digest, then each release
checksum. The destination must be new and outside the repository. Never executes
container contents. On a cluster with Apptainer, use the Apptainer definition.
"""

import argparse
import hashlib
import json
import shutil
import tarfile
import urllib.request
from pathlib import Path

# Release layer from ghcr.io/almagx/alma3:3.1.2, inspected 2026-11-11.
LAYER = "d9d0336ded9042c7e09c0bed53168680adac170526c9327c884e11e9f32ecb6d"


def checksum(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    out = args.destination.resolve()
    repository = Path(__file__).resolve().parents[2]
    if out.is_relative_to(repository) or out.exists():
        parser.error("destination must be new and outside the repository")
    out.parent.mkdir(parents=True, exist_ok=True)
    archive_path = out.parent / f"alma3-{LAYER}.tar.gz"
    if not archive_path.exists():
        token_url = "https://ghcr.io/token?scope=repository:almagx/alma3:pull"
        with urllib.request.urlopen(token_url, timeout=60) as response:
            token = json.load(response)["token"]
        request = urllib.request.Request(
            f"https://ghcr.io/v2/almagx/alma3/blobs/sha256:{LAYER}",
            headers={"Authorization": "Bearer " + token},
        )
        partial = archive_path.with_suffix(".partial")
        with urllib.request.urlopen(request, timeout=60) as response, partial.open("wb") as stream:
            shutil.copyfileobj(response, stream, length=1024 * 1024)
        if checksum(partial) != LAYER:
            raise ValueError("OCI layer checksum mismatch")
        partial.replace(archive_path)
    elif checksum(archive_path) != LAYER:
        raise ValueError("cached OCI layer checksum mismatch")
    out.mkdir()
    prefix = "opt/alma3/release/"
    with tarfile.open(archive_path, "r|gz") as archive:
        for entry in archive:
            if entry.isdir():
                continue
            if not entry.isfile() or not entry.name.startswith(prefix):
                raise ValueError(f"unexpected container entry: {entry.name}")
            target = out / entry.name[len(prefix) :]
            if not target.resolve().is_relative_to(out) or target.exists():
                raise ValueError("unsafe or duplicate release path")
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.extractfile(entry) as source, target.open("wb") as stream:
                shutil.copyfileobj(source, stream, length=1024 * 1024)
    hashes = json.loads((out / "SHA256SUMS.json").read_text())
    for name, expected in hashes.items():
        target = (out / name).resolve()
        if not target.is_relative_to(out) or checksum(target) != expected:
            raise ValueError(f"release checksum mismatch: {name}")
    print(f"Verified ALMA3 v3.1.2 release layer: {out}")
    print(f"Cached OCI layer (may be removed after verification): {archive_path}")


if __name__ == "__main__":
    main()
