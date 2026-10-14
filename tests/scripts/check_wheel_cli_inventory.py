#!/usr/bin/env python3
"""Check that the Python wheel ships every CLI the C++ build defines.

The wheel declares its command-line drivers in four places, all maintained by
hand except the last:

  python_bindings/CMakeLists.txt        lsdtt_python_cli(...)   builds the binary
  python_bindings/pyproject.toml        [project.scripts]       the pip entry point
  python_bindings/src/lsdtt3/_cli.py    _DRIVERS                dispatch allow-list
  python_bindings/api/public_api.json   console_scripts and declared_drivers

Each must list exactly the root CMakeLists.txt lsdtt_cli() targets. When they
drift, the wheel silently lacks a tool: four CLIs (advection-diffusion,
bankfull-width, distributed-hydrology, slope-area) were missing until
2026-10-14. python_bindings/tests/test_cli_packaging.py checks most of this
too, but it needs pytest and runs only in the bindings suite, so nothing in the
root CTest suite noticed.

Standard library only (tier_python_stdlib).

Usage:
    check_wheel_cli_inventory.py --repo-root <dir>
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys


def root_targets(root: Path) -> set[str]:
    text = (root / "CMakeLists.txt").read_text(encoding="utf-8")
    return set(re.findall(r"^\s*lsdtt_cli\((lsdtt-[\w-]+)", text, re.MULTILINE))


def wheel_cmake(root: Path) -> set[str]:
    text = (root / "python_bindings" / "CMakeLists.txt").read_text(encoding="utf-8")
    return set(re.findall(r"^\s*lsdtt_python_cli\((lsdtt-[\w-]+)", text, re.MULTILINE))


def console_scripts(root: Path) -> set[str]:
    text = (root / "python_bindings" / "pyproject.toml").read_text(encoding="utf-8")
    return set(re.findall(r'^(lsdtt-[\w-]+)\s*=\s*"lsdtt3\._cli:main"$', text, re.MULTILINE))


def dispatch_drivers(root: Path) -> set[str]:
    text = (root / "python_bindings" / "src" / "lsdtt3" / "_cli.py").read_text(encoding="utf-8")
    match = re.search(r"_DRIVERS\s*=\s*frozenset\(\s*\{(.*?)\}\s*\)", text, re.DOTALL)
    if match is None:
        raise SystemExit("ERROR: could not find _DRIVERS = frozenset({...}) in _cli.py")
    return set(re.findall(r'"(lsdtt-[\w-]+)"', match.group(1)))


def manifest_lists(root: Path) -> tuple[set[str], set[str]]:
    path = root / "python_bindings" / "api" / "public_api.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    return set(manifest["console_scripts"]), set(manifest["declared_drivers"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--repo-root", required=True, type=Path)
    root = parser.parse_args().repo_root.resolve()

    expected = root_targets(root)
    if not expected:
        print("ERROR: found no lsdtt_cli() targets in the root CMakeLists.txt")
        return 1

    manifest_scripts, manifest_drivers = manifest_lists(root)
    declared = {
        "python_bindings/CMakeLists.txt lsdtt_python_cli()": wheel_cmake(root),
        "python_bindings/pyproject.toml [project.scripts]": console_scripts(root),
        "python_bindings/src/lsdtt3/_cli.py _DRIVERS": dispatch_drivers(root),
        "python_bindings/api/public_api.json console_scripts": manifest_scripts,
        "python_bindings/api/public_api.json declared_drivers": manifest_drivers,
    }

    failed = False
    for where, names in declared.items():
        missing = sorted(expected - names)
        extra = sorted(names - expected)
        if missing or extra:
            failed = True
            print(f"MISMATCH: {where}")
            for name in missing:
                print(f"  missing: {name}")
            for name in extra:
                print(f"  not a root lsdtt_cli() target: {name}")

    if failed:
        print(
            "\nEvery root lsdtt_cli() target must be declared in all four wheel "
            "locations.\nSee 'Adding a new CLI' in "
            "docs/developer_notes/developer_onboarding.md. Regenerate "
            "public_api.json with\n  python python_bindings/tools/"
            "generate_public_api_manifest.py --update")
        return 1

    print(f"Wheel CLI inventory matches the {len(expected)} root lsdtt_cli() targets "
          f"in all {len(declared)} locations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
