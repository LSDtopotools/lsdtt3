#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Generate docs/developer_notes/test_inventory.md: what each part of lsdtt3 is
tested with, by evidence category (how a test knows its answer) and level.

Sources:
  * CTest: `ctest --show-only=json-v1` on a configured build directory, which
    gives every registered test's labels and command;
  * pytest: the markers of python_bindings/tests, read statically by
    tests/scripts/check_python_test_markers.py (no built binding needed).

The set of registered CTest tests depends on configure options, so the
committed inventory is generated from ONE stated configuration (CANONICAL
below). With any other configuration, --check skips (exit 77) and --update
refuses, rather than producing a different file.

Usage:
  generate_test_inventory.py --build-dir build --update
  generate_test_inventory.py --build-dir build --check      # CTest: unit_test_inventory_is_current
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUTPUT = REPO / "docs" / "developer_notes" / "test_inventory.md"
EXIT_SKIP = 77

# The configuration the committed inventory describes: everything a developer
# build registers except the benchmark and performance tiers, which are
# reported, not asserted, and need reference packages.
CANONICAL = {
    "LSDTT_ENABLE_PYTHON_TESTS": "ON",
    "LSDTT_ENABLE_PYTHON_GEO_TESTS": "ON",
    "LSDTT_ENABLE_NETWORK_TESTS": "ON",
    "LSDTT_ENABLE_BENCHMARKS": "OFF",
    "LSDTT_ENABLE_PERF_TESTS": "OFF",
}

# Labels that are neither domains, evidence nor level: tiers, requirements,
# cost and interface flags.
NON_DOMAIN = {
    "tier_native", "tier_python_stdlib", "tier_python_geo", "tier_pixi",
    "hermetic", "heavy", "needs_network", "needs_gdal_tools", "cli", "python",
    "fixture", "lem_component", "configure_contract", "benchmark",
}

sys.path.insert(0, str(REPO / "tests" / "scripts"))
from check_python_test_markers import (  # noqa: E402
    EVIDENCE_MARKERS, LEVEL_MARKERS, check_source)

EVIDENCE = list(EVIDENCE_MARKERS)
LEVELS = list(LEVEL_MARKERS)
SHORT = {"exact_solution": "exact", "property": "prop", "metamorphic": "meta",
         "differential": "diff", "characterisation": "char", "smoke": "smoke",
         "negative": "neg", "performance": "perf", "conformance": "conf"}
STRONG = {"exact_solution", "property", "metamorphic", "differential"}


def read_cache(build: Path) -> dict[str, str]:
    values = {}
    cache = build / "CMakeCache.txt"
    if not cache.exists():
        return values
    for line in cache.read_text().splitlines():
        m = re.match(r"^([A-Z0-9_]+):[A-Z]+=(.*)$", line)
        if m:
            values[m.group(1)] = m.group(2)
    return values


def configuration_mismatch(build: Path) -> list[str]:
    cache = read_cache(build)
    out = []
    for key, want in CANONICAL.items():
        have = cache.get(key, "OFF").upper()
        have = "ON" if have in ("ON", "TRUE", "1", "YES") else "OFF"
        if have != want:
            out.append(f"{key}={have} (inventory uses {want})")
    return out


def ctest_tests(build: Path, ctest: str) -> list[dict]:
    raw = subprocess.run([ctest, "--test-dir", str(build), "--show-only=json-v1"],
                         check=True, capture_output=True, text=True).stdout
    tests = []
    for t in json.loads(raw)["tests"]:
        labels = next((p["value"] for p in t.get("properties", [])
                       if p["name"] == "LABELS"), [])
        tests.append({"name": t["name"], "labels": labels,
                      "command": t.get("command", []) or []})
    return tests


def cli_names() -> list[str]:
    text = (REPO / "CMakeLists.txt").read_text()
    return sorted(set(re.findall(r"^\s*lsdtt_cli\((lsdtt-[A-Za-z0-9_-]+)", text, re.M)))


def clis_in_command(command: list[str], clis: set[str]) -> set[str]:
    found = set()
    for arg in command:
        for token in re.split(r"[=:;,\s]", arg):
            name = Path(token).name
            if name.endswith(".exe"):
                name = name[:-4]
            if name in clis:
                found.add(name)
    return found


def python_items() -> dict[str, dict[str, list[str]]]:
    """file -> {item id -> markers}"""
    out = {}
    for path in sorted((REPO / "python_bindings" / "tests").glob("test_*.py")):
        result = check_source(path.read_text(), path.name)
        out[path.name] = result.items
    return out


def table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |",
             "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(lines)


def counts_row(counter: Counter, keys: list[str]) -> list[str]:
    return [str(counter[k]) if counter[k] else "." for k in keys]


def render(tests: list[dict], pyfiles: dict[str, dict[str, list[str]]]) -> str:
    clis = cli_names()
    cli_set = set(clis)
    ev_keys = EVIDENCE
    out: list[str] = []
    out.append("# Test inventory\n")
    out.append("<!-- GENERATED by tools/generate_test_inventory.py -- do not edit by hand. -->\n")
    out.append(
        "What lsdtt3 is tested with, by **evidence** (how a test knows its answer) "
        "and **level** (unit or system). The categories are defined in "
        "[`lsdtt3_architecture.md`](lsdtt3_architecture.md#testing-categories) and "
        "[`testing_and_benchmarking.md`](testing_and_benchmarking.md). A test with "
        "several evidence labels is counted under each, so rows sum to more than "
        "the number of tests. `.` means none.\n")
    cfg = ", ".join(f"`{k}={v}`" for k, v in CANONICAL.items())
    out.append(f"CTest counts are for the configuration {cfg}; the benchmark and "
               "performance tiers are not included. Python counts are test functions "
               "and methods (a parametrised function counts once).\n")
    out.append("Regenerate with `python3 tools/generate_test_inventory.py --build-dir "
               "build --update` after adding or relabelling tests; the CTest test "
               "`unit_test_inventory_is_current` fails while this file is stale.\n")

    # ---- totals ---------------------------------------------------------
    ct_ev, ct_lv = Counter(), Counter()
    for t in tests:
        for lab in t["labels"]:
            if lab in EVIDENCE_MARKERS:
                ct_ev[lab] += 1
            elif lab in LEVEL_MARKERS:
                ct_lv[lab] += 1
    py_ev, py_lv = Counter(), Counter()
    n_py = 0
    for items in pyfiles.values():
        for marks in items.values():
            n_py += 1
            for m in set(marks):
                if m in EVIDENCE_MARKERS:
                    py_ev[m] += 1
                elif m in LEVEL_MARKERS:
                    py_lv[m] += 1
    out.append("## Totals\n")
    rows = [["CTest", str(len(tests))] + counts_row(ct_ev, ev_keys) + counts_row(ct_lv, LEVELS),
            ["pytest", str(n_py)] + counts_row(py_ev, ev_keys) + counts_row(py_lv, LEVELS)]
    out.append(table(["Suite", "Tests"] + [f"`{k}`" for k in ev_keys] + [f"`{k}`" for k in LEVELS], rows) + "\n")

    # ---- CTest by domain --------------------------------------------------
    by_domain: dict[str, Counter] = defaultdict(Counter)
    n_domain: Counter = Counter()
    for t in tests:
        labs = set(t["labels"])
        domains = sorted(labs - NON_DOMAIN - set(EVIDENCE_MARKERS) - set(LEVEL_MARKERS))
        for d in domains:
            n_domain[d] += 1
            for lab in labs & set(EVIDENCE_MARKERS):
                by_domain[d][lab] += 1
    out.append("## CTest by domain\n")
    out.append("A test in several domains (for example `lithology`, which cuts across "
               "compartments) is counted in each.\n")
    rows = [[f"`{d}`", str(n_domain[d])] + counts_row(by_domain[d], ev_keys)
            for d in sorted(n_domain)]
    out.append(table(["Domain", "Tests"] + [SHORT[k] for k in ev_keys], rows) + "\n")

    weak = [d for d in sorted(n_domain) if not (set(by_domain[d]) & STRONG)]
    out.append("### Domains with no exact-solution, property, metamorphic or differential test\n")
    out.append("These are checked only by smoke, negative, characterisation or "
               "conformance tests: they show that the code runs, refuses bad input "
               "or has not changed, but not that its answers are right.\n")
    out.append(("\n".join(f"- `{d}` ({n_domain[d]} test{'' if n_domain[d] == 1 else 's'})" for d in weak) if weak
                else "None.") + "\n")

    # ---- CLIs ------------------------------------------------------------
    by_cli: dict[str, Counter] = defaultdict(Counter)
    n_cli: Counter = Counter()
    for t in tests:
        for c in clis_in_command(t["command"], cli_set):
            n_cli[c] += 1
            for lab in set(t["labels"]) & set(EVIDENCE_MARKERS):
                by_cli[c][lab] += 1
    out.append("## Command-line programs\n")
    out.append("CTest tests whose command runs the program, directly or through the "
               "harness. Unit tests of the routines a CLI calls are not counted here; "
               "see the domain table.\n")
    rows = [[f"`{c}`", str(n_cli[c])] + counts_row(by_cli[c], ev_keys) for c in clis]
    out.append(table(["CLI", "Tests"] + [SHORT[k] for k in ev_keys], rows) + "\n")
    untested = [c for c in clis if n_cli[c] == 0]
    only_weak = [c for c in clis if n_cli[c] and not (set(by_cli[c]) & STRONG)]
    out.append("### CLIs with no system test\n")
    out.append(("\n".join(f"- `{c}`" for c in untested) if untested else "None.") + "\n")
    out.append("### CLIs whose system tests are only smoke, negative, characterisation or conformance\n")
    out.append(("\n".join(f"- `{c}`" for c in only_weak) if only_weak else "None.") + "\n")

    # ---- pytest by file --------------------------------------------------
    out.append("## Python binding tests by file\n")
    rows = []
    for fname, items in pyfiles.items():
        c: Counter = Counter()
        lv: Counter = Counter()
        for marks in items.values():
            for m in set(marks):
                if m in EVIDENCE_MARKERS:
                    c[m] += 1
                elif m in LEVEL_MARKERS:
                    lv[m] += 1
        rows.append([f"`{fname}`", str(len(items))] + counts_row(c, ev_keys)
                    + counts_row(lv, LEVELS))
    out.append(table(["File", "Tests"] + [SHORT[k] for k in ev_keys] + LEVELS, rows) + "\n")

    out.append("Abbreviations: " + ", ".join(f"{v} = `{k}`" for k, v in SHORT.items()) + ".\n")
    return "\n".join(out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-dir", required=True, type=Path)
    ap.add_argument("--ctest", default="ctest")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--update", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = ap.parse_args(argv)

    mismatch = configuration_mismatch(args.build_dir)
    if mismatch:
        msg = ("test inventory: build is not configured as the committed inventory "
               "describes: " + "; ".join(mismatch))
        if args.check:
            print("SKIP " + msg)
            return EXIT_SKIP
        print("ERROR " + msg + ". Reconfigure, or regenerate from a build that matches.")
        return 2

    text = render(ctest_tests(args.build_dir, args.ctest), python_items())
    if args.update:
        OUTPUT.write_text(text)
        print(f"wrote {OUTPUT.relative_to(REPO)}")
        return 0
    current = OUTPUT.read_text() if OUTPUT.exists() else ""
    if current != text:
        print(f"FAIL {OUTPUT.relative_to(REPO)} is stale. Regenerate with:\n"
              f"  python3 tools/generate_test_inventory.py --build-dir <build> --update")
        return 1
    print(f"PASS {OUTPUT.relative_to(REPO)} is current")
    return 0


if __name__ == "__main__":
    sys.exit(main())
