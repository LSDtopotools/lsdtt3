#!/usr/bin/env python3

"""Generate/check the explicit Phase 3 Python API documentation backlog."""

import argparse
import json
from pathlib import Path
import sys


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--update", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    manifest = json.loads((root / "api" / "public_api.json").read_text(encoding="utf-8"))
    output = root / "api" / "documentation_exclusions.json"
    data = {
        "schema_version": 1,
        "reason": "Explicit Phase 3 documentation backlog",
        "symbols": manifest["documentation_facade"]["missing"],
    }
    rendered = json.dumps(data, indent=2) + "\n"
    if args.update:
        output.write_text(rendered, encoding="utf-8")
        print(f"updated {output} with {len(data['symbols'])} symbols")
        return 0
    if output.read_text(encoding="utf-8") != rendered:
        raise RuntimeError("Python API documentation exclusions are stale")
    print(f"Python API documentation exclusions are current: {len(data['symbols'])}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
