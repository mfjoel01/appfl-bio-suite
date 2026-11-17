"""ALMA3 study: validation, freeze, then one test evaluation."""

import argparse
import json
import os
import shlex
from pathlib import Path

from .config import Config


def main():
    # Avoid synchronous huge-page compaction for large model-update arrays on
    # fragmented HPC login/compute nodes. Explicit caller overrides are retained.
    os.environ.setdefault("NUMPY_MADVISE_HUGEPAGE", "0")
    from .study import freeze, run, test_frozen

    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    study = sub.add_parser("run")
    study.add_argument("--out", type=Path, required=True)
    study.add_argument("--config", type=Path, help="JSON Config fields")
    study.add_argument("--fixture", action="store_true", help="Tiny stand-in; NOT ALMA3 weights")
    study.add_argument("--release", type=Path)
    study.add_argument("--member")
    study.add_argument(
        "--worker-command", help="Isolated Python/Apptainer Python command, no shell"
    )
    study.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    study.add_argument("--validation-only", action="store_true", required=True)
    lock = sub.add_parser("freeze")
    lock.add_argument("directory", type=Path)
    test = sub.add_parser("evaluate")
    test.add_argument("frozen", type=Path)
    test.add_argument("--out", type=Path, required=True)
    test.add_argument("--worker-command")
    args = parser.parse_args()
    if args.action == "freeze":
        print(freeze(args.directory))
    elif args.action == "evaluate":
        print(
            test_frozen(
                args.frozen,
                args.out,
                command=shlex.split(args.worker_command) if args.worker_command else None,
            )
        )
    else:
        config = Config(**json.loads(args.config.read_text())) if args.config else Config()
        print(
            run(
                args.out,
                config,
                fixture=args.fixture,
                release=args.release,
                member=args.member,
                device=args.device,
                command=shlex.split(args.worker_command) if args.worker_command else None,
            )
        )


if __name__ == "__main__":
    main()
