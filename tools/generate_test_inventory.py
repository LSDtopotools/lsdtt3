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
import tempfile
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


def git_describe() -> str:
    """`git describe --tags --dirty --always` of the source tree: the release tag,
    commits since it, the commit, and -dirty when the tree has uncommitted
    changes. "not a git checkout" when there is no git or no repository."""
    try:
        done = subprocess.run(["git", "-C", str(REPO), "describe", "--tags", "--dirty",
                               "--always"], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return "not a git checkout"
    return done.stdout.strip() if done.returncode == 0 and done.stdout.strip() \
        else "not a git checkout"


# The git part of the provenance line changes with every commit, so --check
# ignores it; the version number is still compared.
GIT_STAMP = re.compile(r"\(git `[^`]*`\)")


def without_git_stamp(text: str) -> str:
    return GIT_STAMP.sub("(git)", text)


def lsdtt3_version() -> str:
    """The version in project(lsdtt3 VERSION x.y.z), the one shared version."""
    text = (REPO / "CMakeLists.txt").read_text()
    m = re.search(r"^project\(lsdtt3\s+VERSION\s+(\S+)", text, re.M)
    return m.group(1) if m else "unknown"


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
    # Any ctest run in <build> rewrites <build>/Testing/Temporary/LastTest.log,
    # even --show-only, and these tools also run AS tests (the inventory check,
    # unit_summarise_tests_runs): run there, the suite's own log came out empty
    # (found 2026-11-13). So list the tests through a throwaway directory whose
    # CTestTestfile.cmake only points at the build: ctest writes its scratch
    # files there, and the listing is identical (names, commands, properties).
    with tempfile.TemporaryDirectory(prefix="lsdtt_ctest_list_") as proxy:
        Path(proxy, "CTestTestfile.cmake").write_text(
            'subdirs("' + str(build.resolve()).replace("\\", "/") + '")\n')
        raw = subprocess.run([ctest, "--test-dir", proxy, "--show-only=json-v1"],
                             check=True, capture_output=True, text=True).stdout
    configured = configured_commands(build)
    tests = []
    for t in json.loads(raw)["tests"]:
        labels = next((p["value"] for p in t.get("properties", [])
                       if p["name"] == "LABELS"), [])
        # ctest gives no command for a test whose executable has not been built,
        # so a configure-only directory would credit no test to the CLIs it runs
        # directly (lsdtt-fetch-raster -h, ...). The generated CTestTestfile.cmake
        # still holds the command as configured.
        command = t.get("command") or configured.get(t["name"], [])
        tests.append({"name": t["name"], "labels": labels, "command": command})
    return tests


ADD_TEST_LINE = re.compile(r'^add_test\((\S+)\s+(.*)\)\s*$', re.M)
CMAKE_ARG = re.compile(r'"((?:[^"\\]|\\.)*)"|\[(=*)\[(.*?)\]\2\]|(\S+)', re.S)


def configured_commands(build: Path) -> dict[str, list[str]]:
    """Test name -> command, from the CTestTestfile.cmake files of a build tree."""
    out: dict[str, list[str]] = {}
    for testfile in build.rglob("CTestTestfile.cmake"):
        for m in ADD_TEST_LINE.finditer(testfile.read_text(errors="replace")):
            name = m.group(1).strip("[]=")
            args = []
            for a in CMAKE_ARG.finditer(m.group(2)):
                args.append(a.group(1) if a.group(1) is not None
                            else a.group(3) if a.group(3) is not None else a.group(4))
            out.setdefault(name, args)
    return out


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
    out.append(f"**Generated from the lsdtt3 {lsdtt3_version()} source tree** "
               f"(git `{git_describe()}`). See "
               "[how current it is](#how-this-page-is-made-and-how-current-it-is).\n")
    out.append(
        "What lsdtt3 is tested with, by **evidence** (how a test knows its answer) "
        "and **level** (unit or system). Where the categories come from, and what "
        "each means, is in [`testing_and_benchmarking.md`]"
        "(testing_and_benchmarking.md#where-the-test-categories-come-from). A test "
        "with several evidence labels is counted under each, so rows sum to more "
        "than the number of tests. `.` means none.\n")
    cfg = ", ".join(f"`{k}={v}`" for k, v in CANONICAL.items())
    out.append(f"CTest counts are for the configuration {cfg}; the benchmark and "
               "performance tiers are not included. Python counts are test functions "
               "and methods (a parametrised function counts once).\n")
    out.append("## How this page is made and how current it is\n")
    out.append(
        "This page is generated by `tools/generate_test_inventory.py` from two sources: "
        "the CTest tests registered in a build directory configured as above (it "
        "only needs configuring, nothing is compiled), and the pytest markers, "
        "read from the files in `python_bindings/tests`. To regenerate it, from the "
        "repository root:\n")
    opts = " ".join(f"-D{k}={v}" for k, v in CANONICAL.items())
    out.append("```bash\n"
               f"cmake -S . -B build_tests {opts}\n"
               "python3 tools/generate_test_inventory.py --build-dir build_tests --update\n"
               "```\n")
    out.append(
        "**How current it is.** `tools/set_version.py`, which sets a new version "
        "number in every file before a release is tagged, regenerates this page as its "
        "last step, so each release's page is generated from that release's source "
        "and names its version. Between releases, the CTest test "
        "`unit_test_inventory_is_current` fails when this page disagrees with the "
        "tests, or with the version number, of the source tree it is in; it runs in a "
        "build configured as above, when the suite is run locally (continuous "
        "integration does not run it). So between releases the page matches its "
        "commit when that test passed before the commit.\n")
    out.append(
        "**The git stamp** at the top is `git describe --tags --dirty --always` when the "
        "page was generated: the latest release tag, the number of commits since it, "
        "and the commit (`v0.9.5-3-g4c4db5c` is three commits after v0.9.5, at commit "
        "4c4db5c). The page cannot name the commit that contains it, which does not "
        "exist yet when it is generated, so the stamp names the commit it was generated "
        "from. `-dirty` means the working tree had uncommitted changes at the time. That "
        "is the normal case: the page is regenerated after tests are added or "
        "relabelled and before those changes are committed, so it describes the stamped "
        "commit plus the changes committed with it. The currency check ignores the git "
        "stamp (it changes with every commit) but not the version number.\n")
    out.append(
        "**Why the release does not generate it.** The documentation site is built from "
        "the release tag, so a page generated by the release would have to be written "
        "after the source is tagged and then committed. That commit would come after "
        "the tag, so the tagged release, and the site built from it, would not contain "
        "it; moving the tag to include it would change the release after it was made. "
        "The page is therefore generated before the release is tagged, by "
        "`tools/set_version.py`, and committed with the new version number, like the "
        "rest of the documentation.\n")
    out.append(
        "**On the documentation site.** The site is built from release tags (`v*`), so "
        "this page there is the inventory as committed at that tag. If the version at "
        "the top differs from the site's version, the page describes the tests of that "
        "earlier version: the inventory was not regenerated before the release.\n")
    out.append(
        "For counts from any version or configuration, as Markdown, CSV or LaTeX, and a "
        "list of every test with what it does, use `tools/summarise_tests.py` (see "
        "[counting and listing the tests](testing_and_benchmarking.md"
        "#counting-and-listing-the-tests-for-a-paper-or-a-report)).\n")

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
    if without_git_stamp(current) != without_git_stamp(text):
        print(f"FAIL {OUTPUT.relative_to(REPO)} is stale. Regenerate with:\n"
              f"  python3 tools/generate_test_inventory.py --build-dir <build> --update")
        return 1
    print(f"PASS {OUTPUT.relative_to(REPO)} is current")
    return 0


if __name__ == "__main__":
    sys.exit(main())
