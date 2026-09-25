#!/usr/bin/env python
"""Fail if the published network-map preview was built from a real federation.

THE POINT
---------
The project site embeds a live export of the federation viewer. That export is supposed
to come from one of the two committed, fictional configs -- ``website/demo-federation.yaml``
for the site, or ``federation.yaml.example``. The file it must never come from is
``local/federation.yaml``, which names real institutions and carries partner contact
names and email addresses.

The two differ by one ``--federation`` argument in .github/workflows/pages.yml, and the
export succeeds either way. Publishing a partner's email address to a public URL is not
undone by deleting the page afterwards, so the distinction is checked here rather than
remembered.

WHAT IS CHECKED
---------------
* the metadata records one of the fictional configs as its source;
* no ``contacts`` block survived into it; and
* nothing anywhere in it looks like an email address.

Usage::

    python scripts/check_site_preview.py _site/map/network.map.json
"""

from __future__ import annotations

import json
import string
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import Any

EXPECTED_SOURCES = ("website/demo-federation.yaml", "federation.yaml.example")

_LOCAL_CHARS = frozenset(string.ascii_letters + string.digits + "._%+-")
_DOMAIN_CHARS = frozenset(string.ascii_letters + string.digits + ".-")


def _looks_like_email(value: str) -> str | None:
    """The first address-shaped run in ``value``, or None.

    Scanned by hand rather than by regex on purpose. A map's results catalogue can embed
    a multi-megabyte ``data:`` URI, and the obvious ``[\\w.+-]+@...`` pattern backtracks
    quadratically over one of those -- a check that hangs the deploy instead of passing
    it. This walks each string once.
    """
    for index, char in enumerate(value):
        if char != "@":
            continue
        start = index
        while start > 0 and value[start - 1] in _LOCAL_CHARS:
            start -= 1
        end = index + 1
        while end < len(value) and value[end] in _DOMAIN_CHARS:
            end += 1
        domain = value[index + 1 : end].strip(".")
        if start < index and "." in domain:
            return value[start:end]
    return None


def _walk(node: Any, path: str = "$") -> Iterator[tuple[str, str]]:
    """Every (json path, string value) pair in the document."""
    if isinstance(node, dict):
        for key, value in node.items():
            yield from _walk(value, f"{path}.{key}")
    elif isinstance(node, list):
        for index, value in enumerate(node):
            yield from _walk(value, f"{path}[{index}]")
    elif isinstance(node, str):
        yield path, node


def check(path: Path) -> tuple[str, list[str]]:
    """The recorded source, and findings. No findings means the preview is publishable."""
    metadata = json.loads(path.read_text(encoding="utf-8"))
    findings: list[str] = []

    source = str(metadata.get("config", {}).get("source", ""))
    if not any(source.endswith(expected) for expected in EXPECTED_SOURCES):
        findings.append(
            f"built from {source or '(no source recorded)'}, which is not one of "
            f"{', '.join(EXPECTED_SOURCES)}. Fix the --federation argument in "
            f"scripts/build_site.sh."
        )

    for location, value in _walk(metadata):
        address = _looks_like_email(value)
        if address is not None:
            findings.append(f"{location} contains an email address: {address!r}")
        elif ".contacts" in location:
            findings.append(f"{location} is a partner contact block: {value!r}")

    return source, findings


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2

    path = Path(argv[1])
    if not path.is_file():
        print(f"no such file: {path}", file=sys.stderr)
        return 2

    source, findings = check(path)
    if findings:
        print(f"{path} is not safe to publish:\n", file=sys.stderr)
        for finding in findings:
            print(f"  - {finding}", file=sys.stderr)
        return 1

    print(f"{path}: built from {source}, no contacts, no addresses.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
