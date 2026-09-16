#!/usr/bin/env python
"""Check a Globus Compute identity-mapping file, offline, in about a second.

Standalone on purpose: it needs only `globus-identity-mapping`, which ships with
globus-compute-endpoint, so it runs in the endpoint's own environment with nothing extra
installed.

    python validate_mapping.py <path> --identity <coordinator identity> --expect <account>

Pass the exact path that `identity_mapping_config_path` in your endpoint's config.yaml
names. Editing a different copy of the file changes nothing -- and after switching to a
privileged user, that path may still point into a personal home directory.

THE TRAP THIS CATCHES, which reading the file cannot:

  `^` and `$` anchors in `match`. That field is NOT a full regex. The library escapes
  anchors into literal characters and then wraps the expression in its own `^...$`, so an
  anchored pattern can never match and every submission fails with a 422. This is the
  cause in nearly every case.

If this prints PASS, no restart is needed -- the endpoint re-reads the file within about
five seconds.

A 403 rather than a 422 is a different problem entirely and this script cannot see it: it
means the endpoint was started by an unprivileged user, so mapping was ignored.
"""

import argparse
import json
import pathlib
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mapping_file")
    parser.add_argument("--identity", required=True, help="the coordinator's Globus identity")
    parser.add_argument("--expect", required=True, help="the local service account")
    parser.add_argument("--identity-id", default="00000000-0000-0000-0000-000000000000")
    args = parser.parse_args()

    try:
        from globus_identity_mapping.loader import load_mappers
    except ImportError:
        print("ERROR: globus-identity-mapping is not installed.")
        print("It ships with globus-compute-endpoint -- run this in the endpoint's environment.")
        return 2

    path = pathlib.Path(args.mapping_file)
    if not path.is_file():
        print(f"ERROR: no such file: {path}")
        print("The authoritative path is whatever `identity_mapping_config_path` in the")
        print("endpoint's config.yaml points at. Check there first.")
        return 2

    try:
        document = json.loads(path.read_bytes())
    except json.JSONDecodeError as exc:
        print(f"ERROR: {path} is not valid JSON: {exc}")
        return 2

    # `status` matters: the mapper silently skips any identity record that is not
    # "used" or "private", which looks exactly like a bad regex.
    record = [
        {
            "id": args.identity_id,
            "sub": args.identity_id,
            "username": args.identity,
            "status": "used",
        }
    ]

    mapped, anchored = [], False
    for entry in document if isinstance(document, list) else [document]:
        try:
            mappers = load_mappers([entry], None, None)
        except Exception as exc:
            print(f"ERROR: could not load the mapping document: {type(exc).__name__}: {exc}")
            print()
            print("If this mentions an invalid '\\' expression, the `match` field contains an")
            print(
                "escape the mapper does not accept. Only \\. \\? \\* \\| \\( \\) \\\\ are allowed --"
            )
            print("in particular \\- is rejected and makes the whole file invalid.")
            return 2
        for mapper in mappers:
            for mapping in getattr(mapper, "mappings", []):
                print(f"  source={mapping['source']!r}  match={mapping['match']!r}")
                compiled = mapping.get("compiled_match")
                if compiled is not None:
                    print(f"    -> compiles to: {compiled.pattern}")
                if any(a in mapping["match"] for a in ("^", "$")):
                    anchored = True
                    print("    !! contains ^ or $ -- these become LITERAL characters.")
                    print("       Remove them; matching is already anchored for you.")
                print(f"    output: {mapping['output']}")
            try:
                for result in mapper.map_identities(record):
                    for names in result.values():
                        mapped.extend(names)
            except Exception as exc:
                print(f"    !! mapping raised {type(exc).__name__}: {exc}")

    print()
    if args.expect in mapped:
        print(f"PASS: {args.identity} maps to '{args.expect}'.")
        print("No restart needed -- the endpoint re-reads this file within about 5 seconds.")
        return 0

    if mapped:
        print(f"FAIL: maps to {mapped}, but tasks must run as '{args.expect}'.")
        print("Check the `output` field.")
    elif anchored:
        print(f"FAIL: {args.identity} does not map to any local user.")
        print("The `match` field contains ^ or $ anchors. Remove them -- that is the cause.")
    else:
        print(f"FAIL: {args.identity} does not map to any local user.")
        print("Check for a typo, an unescaped dot (write \\.), or that this is really the")
        print("file `identity_mapping_config_path` names.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
