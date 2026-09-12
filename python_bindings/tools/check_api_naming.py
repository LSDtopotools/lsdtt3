#!/usr/bin/env python3

"""Check public Python parameter names against the reviewed violation baseline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=root / "api" / "public_api.json")
    parser.add_argument("--baseline", type=Path, default=root / "api" / "naming_baseline.json")
    parser.add_argument("--update-baseline", action="store_true")
    return parser.parse_args()


def parameter_names(definition: str) -> list[str]:
    return re.findall(r"\b([A-Za-z_]\w*)(?=:)", definition)


def violation_rule(name: str) -> str | None:
    if re.search(r"_(m|px|deg)$", name):
        return "spell_out_unit_suffix"
    if re.search(r"_metres[23]$", name):
        return "use_m2_or_m3_suffix"
    if name.startswith("tol_") or name.endswith("_tol"):
        return "name_tolerance_semantics"
    if name == "threshold_pixels":
        return "qualify_contributing_pixel_threshold"
    if name == "relief_window_pixels":
        return "state_window_side_length"
    return None


def collect_violations(manifest: dict) -> list[dict[str, str]]:
    violations = set()
    for export in manifest["exports"]:
        callable_sets = [("", export["signatures"])]
        callable_sets.extend(
            (member["name"], member["signatures"]) for member in export["members"]
        )
        for member, signatures in callable_sets:
            for signature in signatures:
                for parameter in parameter_names(signature["definition"]):
                    if parameter in {"self", "cls"}:
                        continue
                    rule = violation_rule(parameter)
                    if rule is not None:
                        violations.add((export["name"], member, parameter, rule))
    return [
        {"symbol": symbol, "member": member, "parameter": parameter, "rule": rule}
        for symbol, member, parameter, rule in sorted(violations)
    ]


def main() -> int:
    args = parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    current = collect_violations(manifest)
    rendered = json.dumps({"schema_version": 1, "violations": current}, indent=2) + "\n"
    if args.update_baseline:
        args.baseline.write_text(rendered, encoding="utf-8")
        print(f"updated {args.baseline} with {len(current)} violations")
        return 0
    baseline = json.loads(args.baseline.read_text(encoding="utf-8"))
    expected = baseline["violations"]
    if current != expected:
        current_set = {tuple(item.values()) for item in current}
        expected_set = {tuple(item.values()) for item in expected}
        added = sorted(current_set - expected_set)
        stale = sorted(expected_set - current_set)
        raise RuntimeError(f"Python API naming baseline differs\nnew: {added}\nstale: {stale}")
    print(f"Python API naming baseline is current: {len(current)} violations")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
