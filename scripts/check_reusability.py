#!/usr/bin/env python
"""Fail if coordinator-specific values leak into the tracked tree.

THE POINT
---------
This repository has to work for someone who is not its author. The most likely way that
quietly stops being true is a real identity, endpoint UUID, hostname or institution name
ending up somewhere that is actually loaded -- in `src/`, in a config, or in a
partner-facing template.

Everything works fine on the original author's machine either way, which is exactly why
this needs to be mechanical rather than remembered.

WHERE SUCH VALUES ARE LEGITIMATE
--------------------------------
* federation.yaml.example        -- fictional by construction
* docs/coordinator/reference-deployment.md
                                 -- one deployment cited as an example, clearly labelled
* this file                      -- it has to name the patterns to look for

Anywhere else is a finding.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

# Patterns that indicate a specific person, institution, machine, or allocation.
# Word-boundaried so ordinary English is not matched: "polarised" is not "polaris".
PATTERNS = {
    "coordinator username": r"\bmjoel\b",
    # Negative lookbehind for the public documentation host: citing a vendor user
    # guide in a bibliography is not an identity leak, whereas name@anl.gov is.
    "institutional domain": r"(?<!docs\.alcf\.)\banl\.gov\b",
    "partner institution": r"\bc-?dac\b|\bcovenant\s+university\b|\bmbzuai\b",
    "cluster hostname": r"\bpolaris\b|\bsophia\b",
    "allocation name": r"\bGeomicVar\b",
    "live endpoint UUID": r"\bc540b588\b|\b84a675b1\b|\bbcc264f2\b",
    "identity UUID": r"\b03c0f4c1\b",
    "shared drive link": r"drive\.google\.com/drive/folders",
}

ALLOWED = {
    "federation.yaml.example",
    "docs/coordinator/reference-deployment.md",
    "scripts/check_reusability.py",
}

# Patterns that are legitimate under a path prefix, for a stated reason. Narrower than
# allowlisting the file: a live endpoint UUID or a coordinator username appearing under
# one of these prefixes is still a finding.
#
# This exists because "coordinator-specific" and "matches the pattern" came apart when the
# third experiment landed, and widening PATTERNS or ALLOWED would both have been the wrong
# fix -- one stops catching real leaks everywhere, the other stops catching them here.
SCOPED_EXEMPTIONS = [
    (
        "docs/experiments/fine-mapping/",
        {"partner institution", "cluster hostname"},
        "the three cohort names are the published study design's, not this deployment's",
    ),
    (
        "src/appfl_bio_suite/experiments/fine_mapping/",
        {"partner institution"},
        "anl/covenant/mbzuai are simulated-cohort identifiers baked into the vendored "
        "sampler's SITE_ORDER; a stranger runs the experiment with those labels and no "
        "relationship to those institutions",
    ),
    (
        "src/appfl_bio_suite/experiments/fine_mapping/fedfm/",
        {"cluster hostname"},
        "vendored byte-for-byte from upstream; a docstring there names the cluster it was "
        "written for, and editing it would break the guarantee that the copy is verbatim",
    ),
    (
        "tests/test_fine_mapping_",
        {"partner institution"},
        "the ported upstream tests assert against those cohort identifiers",
    ),
    (
        "scripts/fine-mapping/",
        {"partner institution", "cluster hostname"},
        "PBS scripts for a named scheduler necessarily name it; the allocation and every "
        "path in them is a variable",
    ),
]


def _exempt(rel: str, label: str) -> bool:
    """Is this pattern legitimate at this path?"""
    return any(
        rel.startswith(prefix) and label in labels
        for prefix, labels, _reason in SCOPED_EXEMPTIONS
    )

# Never scanned: not shipped, not tracked, or not text.
SKIP_DIRS = {".git", "local", ".venv", "venv", "__pycache__", ".pytest_cache",
             ".ruff_cache", ".mypy_cache", "build", "dist", "node_modules"}
SKIP_SUFFIXES = {".png", ".jpg", ".gz", ".bed", ".bim", ".fam", ".pyc", ".whl", ".pdf"}


def tracked_files(root: Path) -> list[Path]:
    """Prefer git's view of the tree; fall back to a walk before the first commit."""
    try:
        result = subprocess.run(
            ["git", "-C", str(root), "ls-files"],
            capture_output=True, text=True, timeout=30,
        )
        if result.returncode == 0 and result.stdout.strip():
            return [root / line for line in result.stdout.splitlines() if line.strip()]
    except (OSError, subprocess.SubprocessError):
        pass

    out = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if any(part in SKIP_DIRS for part in path.relative_to(root).parts):
            continue
        out.append(path)
    return out


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    findings: list[tuple[str, int, str, str]] = []
    scanned = 0

    for path in tracked_files(root):
        if not path.is_file() or path.suffix.lower() in SKIP_SUFFIXES:
            continue
        rel = path.relative_to(root).as_posix()
        if any(part in SKIP_DIRS for part in path.relative_to(root).parts):
            continue
        if rel in ALLOWED:
            continue

        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        scanned += 1

        for lineno, line in enumerate(text.splitlines(), 1):
            for label, pattern in PATTERNS.items():
                if _exempt(rel, label):
                    continue
                if re.search(pattern, line, re.IGNORECASE):
                    findings.append((rel, lineno, label, line.strip()[:100]))

    print(f"Scanned {scanned} tracked text files.")
    print(f"Allowlisted: {', '.join(sorted(ALLOWED))}")
    for prefix, labels, reason in SCOPED_EXEMPTIONS:
        print(f"Scoped:      {prefix}  [{', '.join(sorted(labels))}] -- {reason}")
    print()

    if not findings:
        print("PASS: no coordinator-specific values outside the allowlist.")
        return 0

    print(f"FAIL: {len(findings)} occurrence(s) found.\n")
    for rel, lineno, label, snippet in findings:
        print(f"  {rel}:{lineno}")
        print(f"    [{label}] {snippet}")
    print()
    print("These make the repository depend on one particular person, institution, or")
    print("machine. Move the value into federation.yaml, or into one of the allowlisted")
    print("documents where it is clearly labelled as an example.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
