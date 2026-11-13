#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Summarise the lsdtt3 test suite: counts by evidence category, or every test
with what it does.

Two tables, from the tests registered in a configured build directory (CTest)
and the Python binding tests (pytest, read statically: no built binding needed):

  counts  how many tests carry each evidence category (how a test knows its
          answer) and each level (unit or system). For a paper or a report.
  list    every test: its suite, level, categories, domains, the programs it
          runs, what the CLI harness checks, and a description taken from the
          test's own comments (see "Descriptions" below).

Formats: md (Markdown, the default), csv, or tex (LaTeX, booktabs).

Unlike tools/generate_test_inventory.py, which writes the committed
docs/developer_notes/test_inventory.md for one fixed configuration, this reports
on whatever the build directory is configured with, and says so: the options
LSDTT_ENABLE_PYTHON_TESTS, _PYTHON_GEO_TESTS, _NETWORK_TESTS, _BENCHMARKS and
_PERF_TESTS change which CTest tests exist.

Descriptions, in order of preference:
  * CTest: the comment block directly above the test's registration in
    tests/CMakeLists.txt; else the header comment of its parameter file
    (tests/params/*.param.in); else the nearest comment above it, which often
    describes a group of tests ("section"); else none. Tests registered in a
    loop are found by the longest trailing part of their name in the file.
  * pytest: the test's docstring, else its class's docstring.
  * Either, failing those: the test's name, as words ("name").
The "source" column says which was used.

Usage:
  python3 tools/summarise_tests.py --build-dir build                     # counts, Markdown
  python3 tools/summarise_tests.py --build-dir build --format tex        # counts, LaTeX
  python3 tools/summarise_tests.py --build-dir build --table list --format csv -o tests.csv
"""

from __future__ import annotations

import argparse
import ast
import csv
import io
import re
import sys
from collections import Counter
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tools"))
sys.path.insert(0, str(REPO / "tests" / "scripts"))
from check_python_test_markers import EVIDENCE_MARKERS, LEVEL_MARKERS  # noqa: E402
from generate_test_inventory import (  # noqa: E402
    CANONICAL, NON_DOMAIN, cli_names, clis_in_command, ctest_tests, python_items,
    read_cache)

EVIDENCE = list(EVIDENCE_MARKERS)
LEVELS = list(LEVEL_MARKERS)
CMAKELISTS = REPO / "tests" / "CMakeLists.txt"
PARAMS_DIR = REPO / "tests" / "params"
PY_TESTS = REPO / "python_bindings" / "tests"

# A comment this many lines or fewer above a registration belongs to it.
DIRECT_GAP = 2
# How far up to look for a section comment.
SECTION_REACH = 120


# -----------------------------------------------------------------------------
# CTest descriptions
# -----------------------------------------------------------------------------
REGISTRATION = re.compile(r"^\s*(add_test\(\s*(?:NAME)?|add_lsdtt_\w+\(|lsdtt_add_\w+\()\s*([A-Za-z0-9_]+)?")


def registration_lines(lines: list[str]) -> dict[str, int]:
    """Test name -> index of the line that registers it (literal names only)."""
    out: dict[str, int] = {}
    for i, line in enumerate(lines):
        m = REGISTRATION.match(line)
        if not m:
            continue
        name = m.group(2)
        if not name and i + 1 < len(lines) and m.group(1).startswith("add_test"):
            nxt = re.match(r"\s*NAME\s+([A-Za-z0-9_]+)", lines[i + 1])
            name = nxt.group(1) if nxt else None
        elif not name and i + 1 < len(lines):
            nxt = re.match(r"\s*([A-Za-z0-9_]+)\s*$", lines[i + 1]) or \
                re.match(r"\s*([A-Za-z0-9_]+)\s", lines[i + 1])
            name = nxt.group(1) if nxt else None
        if name and name not in out:
            out[name] = i
    return out


def clean_comment(block: list[str]) -> str:
    text = []
    for raw in block:
        body = raw.lstrip()[1:].strip() if raw.lstrip().startswith("#") else raw.strip()
        if not body or set(body) <= set("=-*#"):
            continue
        text.append(body)
    return re.sub(r"\s+", " ", " ".join(text)).strip()


def comment_block_ending_at(lines: list[str], j: int) -> tuple[list[str], int]:
    """The contiguous comment lines ending at line j, and the index above them."""
    block = []
    while j >= 0 and lines[j].lstrip().startswith("#"):
        block.append(lines[j])
        j -= 1
    return list(reversed(block)), j


CONTROL = re.compile(r"^\s*(if|elseif|else|foreach|endforeach|endif)\s*\(")


def direct_comment(lines: list[str], i: int) -> str:
    """The comment just above line i: blank lines and if()/foreach() lines between
    them are allowed, other code is not."""
    j, gap = i - 1, 0
    while j >= 0 and (not lines[j].strip() or CONTROL.match(lines[j])):
        gap += 0 if CONTROL.match(lines[j]) else 1
        if gap > DIRECT_GAP:
            return ""
        j -= 1
    if j < 0 or not lines[j].lstrip().startswith("#"):
        return ""
    return clean_comment(comment_block_ending_at(lines, j)[0])


def section_comment(lines: list[str], i: int) -> str:
    j = i - 1
    while j >= 0 and i - j <= SECTION_REACH:
        if lines[j].lstrip().startswith("#"):
            return clean_comment(comment_block_ending_at(lines, j)[0])
        j -= 1
    return ""


def param_header(command: list[str]) -> str:
    for arg in command:
        m = re.search(r"([A-Za-z0-9_.-]+\.param)(\.in)?$", arg)
        if not m:
            continue
        for cand in (PARAMS_DIR / (m.group(1) + ".in"), PARAMS_DIR / m.group(1), Path(arg)):
            if cand.is_file():
                header = []
                for line in cand.read_text(errors="replace").splitlines():
                    if line.lstrip().startswith("#"):
                        header.append(line)
                    elif line.strip() or header:
                        break
                return clean_comment(header)
    return ""


# What the CLI harness (tests/scripts/run_cli_smoke.py) is asked to check.
CHECK_FLAGS = {
    "--expect": "output written",
    "--expect-absent": "output not written",
    "--expect-file-contains": "file contains",
    "--expect-file-not-contains": "file lacks",
    "--stdout-contains": "stdout contains",
    "--stderr-contains": "stderr contains",
    "--stderr-not-contains": "stderr lacks",
    "--expect-raster-checksum": "raster checksum",
    "--expect-sha256": "sha256",
    "--expect-raster-size": "raster size",
    "--expect-raster-stat": "raster statistic",
    "--expect-raster-values": "raster values",
    "--expect-raster-min": "raster minimum",
    "--expect-raster-epsg": "raster EPSG",
    "--expect-valid-percent": "valid-pixel percent",
    "--expect-valid-percent-100": "no nodata",
    "--expect-vector-feature-min": "vector feature count",
    "--compare-rasters": "rasters compared",
    "--compare-raster-values": "raster values compared",
    "--expect-exit-code": "exit code",
}


def harness_checks(command: list[str]) -> str:
    counts: Counter = Counter()
    failure = ""
    for k, arg in enumerate(command):
        flag = arg.split("=")[0]
        if flag == "--expect-failure":
            failure = "fails as expected"
        elif flag in CHECK_FLAGS:
            counts[CHECK_FLAGS[flag]] += 1
            if flag == "--stderr-contains" and failure and k + 1 < len(command):
                failure = f'fails with "{command[k + 1][:140]}"'
    parts = [failure] if failure else []
    parts += [f"{name} ({n})" if n > 1 else name for name, n in counts.items()]
    return "; ".join(parts)


def programs(command: list[str], clis: set[str]) -> str:
    found = set(clis_in_command(command, clis))
    for arg in command:
        name = Path(arg).name
        if name.startswith("lsdtt-test-"):
            found.add(name)
    return ", ".join(sorted(found))


def anchor_line(lines: list[str], name: str) -> tuple[int | None, bool]:
    """For a test with no literal registration: the first line naming it, or else
    the first naming the longest trailing part of its name (at least two words),
    as a test generated in a loop from a list is. The flag is True for a part."""
    words = name.split("_")
    for k in range(len(words) - 1):
        part = "_".join(words[k:])
        pattern = re.compile(r"(?<![A-Za-z0-9_])" + re.escape(part) + r"(?![A-Za-z0-9_])")
        for i, line in enumerate(lines):
            if not line.lstrip().startswith("#") and pattern.search(line):
                return i, k > 0
    return None, False


def ctest_rows(build: Path, ctest: str) -> list[dict]:
    lines = CMAKELISTS.read_text().splitlines()
    where = registration_lines(lines)
    clis = set(cli_names())
    rows = []
    for t in ctest_tests(build, ctest):
        labels = set(t["labels"])
        description, source = "", ""
        i = where.get(t["name"])
        if i is None:
            i, _ = anchor_line(lines, t["name"])
        if i is not None:
            description = direct_comment(lines, i)
            source = "comment" if description else ""
        if not description:
            description = param_header(t["command"])
            source = "parameter file" if description else ""
        if not description and i is not None:
            description = section_comment(lines, i)
            source = "section" if description else ""
        if not description:
            description, source = name_as_words(t["name"]), "name"
        rows.append({
            "suite": "ctest",
            "test": t["name"],
            "level": ";".join(sorted(labels & set(LEVELS))),
            "evidence": ";".join(e for e in EVIDENCE if e in labels),
            "domains": ";".join(sorted(labels - NON_DOMAIN - set(EVIDENCE) - set(LEVELS))),
            "programs": programs(t["command"], clis),
            "checks": harness_checks(t["command"]),
            "description": description,
            "source": source,
        })
    return rows


def name_as_words(name: str) -> str:
    """test_swath_profile_leaves_out_pixels -> "Swath profile leaves out pixels"."""
    words = name.split("_")
    if words and words[0] in ("test", "unit", "smoke", "failure", "rigorous", "regression"):
        words = words[1:]
    text = " ".join(w for w in words if w)
    return text[:1].upper() + text[1:] if text else ""


# -----------------------------------------------------------------------------
# pytest descriptions
# -----------------------------------------------------------------------------
def docstrings(path: Path) -> dict[str, str]:
    """pytest item id (as check_python_test_markers forms it) -> docstring."""
    tree = ast.parse(path.read_text(), filename=path.name)
    out = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            out[f"{path.name}::{node.name}"] = ast.get_docstring(node) or ""
        elif isinstance(node, ast.ClassDef):
            class_doc = ast.get_docstring(node) or ""
            for member in node.body:
                if isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    doc = ast.get_docstring(member) or ""
                    out[f"{path.name}::{node.name}::{member.name}"] = (
                        doc, "docstring") if doc else (class_doc, "class docstring" if class_doc else "")
    return out


def pytest_rows() -> list[dict]:
    rows = []
    for fname, items in python_items().items():
        docs = docstrings(PY_TESTS / fname)
        for item, marks in items.items():
            doc = docs.get(item, "")
            text, source = doc if isinstance(doc, tuple) else (doc, "docstring" if doc else "")
            if not text:
                text, source = name_as_words(item.split("::")[-1]), "name"
            rows.append({
                "suite": "pytest",
                "test": item,
                "level": ";".join(m for m in LEVELS if m in marks),
                "evidence": ";".join(e for e in EVIDENCE if e in marks),
                "domains": Path(fname).stem.removeprefix("test_"),
                "programs": "",
                "checks": "",
                "description": re.sub(r"\s+", " ", text).strip(),
                "source": source,
            })
    return rows


# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------
def configuration(build: Path) -> str:
    cache = read_cache(build)
    parts = []
    for key in CANONICAL:
        value = cache.get(key, "OFF").upper()
        parts.append(f"{key}={'ON' if value in ('ON', 'TRUE', '1', 'YES') else 'OFF'}")
    return ", ".join(parts)


def count_table(rows: list[dict]) -> list[list[str]]:
    by_suite = {s: [r for r in rows if r["suite"] == s] for s in ("ctest", "pytest")}
    out = []
    for key in EVIDENCE + LEVELS:
        field = "evidence" if key in EVIDENCE else "level"
        n = {s: sum(key in r[field].split(";") for r in rs) for s, rs in by_suite.items()}
        out.append([key, "evidence" if field == "evidence" else "level",
                    str(n["ctest"]), str(n["pytest"]), str(n["ctest"] + n["pytest"])])
    n_ct, n_py = len(by_suite["ctest"]), len(by_suite["pytest"])
    out.append(["all tests", "", str(n_ct), str(n_py), str(n_ct + n_py)])
    return out


def definition(key: str) -> str:
    return EVIDENCE_MARKERS.get(key) or LEVEL_MARKERS.get(key) or "Every registered test."


def tex_escape(s: str) -> str:
    for a, b in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("$", r"\$"),
                 ("#", r"\#"), ("_", r"\_"), ("{", r"\{"), ("}", r"\}"), ("~", r"\textasciitilde{}"),
                 ("^", r"\textasciicircum{}")):
        s = s.replace(a, b)
    return s


def short(text: str, limit: int = 300) -> str:
    if len(text) <= limit:
        return text
    cut = text[:limit]
    stop = cut.rfind(". ")
    return (cut[:stop + 1] if stop > limit // 2 else cut.rstrip() + " …")


def render(rows: list[dict], table: str, fmt: str, config: str) -> str:
    header_note = (f"lsdtt3 tests on {date.today().isoformat()}; CTest configuration: {config}. "
                   "A test with several evidence categories is counted under each, so the "
                   "category rows sum to more than the number of tests.")
    buf = io.StringIO()
    if table == "counts":
        counts = count_table(rows)
        if fmt == "csv":
            w = csv.writer(buf)
            w.writerow(["category", "kind", "ctest", "pytest", "total", "definition"])
            for r in counts:
                w.writerow(r + [definition(r[0])])
        elif fmt == "tex":
            buf.write("% " + header_note + "\n")
            buf.write("\\begin{tabular}{llrrr}\n\\toprule\n"
                      "Category & Kind & CTest & pytest & Total \\\\\n\\midrule\n")
            for r in counts:
                if r[0] == "all tests":
                    buf.write("\\midrule\n")
                buf.write(" & ".join(tex_escape(c) for c in r) + " \\\\\n")
            buf.write("\\bottomrule\n\\end{tabular}\n")
        else:
            buf.write("# lsdtt3 tests by category\n\n" + header_note + "\n\n")
            buf.write("| Category | Kind | CTest | pytest | Total | What it means |\n"
                      "|---|---|---|---|---|---|\n")
            for r in counts:
                name = f"**{r[0]}**" if r[0] == "all tests" else f"`{r[0]}`"
                buf.write(f"| {name} | {r[1]} | {r[2]} | {r[3]} | {r[4]} | {definition(r[0])} |\n")
        return buf.getvalue()

    fields = ["suite", "test", "level", "evidence", "domains", "programs", "checks",
              "description", "source"]
    if fmt == "csv":
        w = csv.DictWriter(buf, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    elif fmt == "tex":
        buf.write("% " + header_note + "\n")
        buf.write("\\begin{longtable}{p{0.22\\linewidth}p{0.18\\linewidth}p{0.52\\linewidth}}\n"
                  "\\toprule\nTest & Evidence & Description \\\\\n\\midrule\n\\endhead\n")
        for r in rows:
            desc = short(r["description"] or r["checks"])
            buf.write(f"{tex_escape(r['test'])} & {tex_escape(r['evidence'].replace(';', ', '))} "
                      f"& {tex_escape(desc)} \\\\\n")
        buf.write("\\bottomrule\n\\end{longtable}\n")
    else:
        sources = Counter(r["source"] or "none" for r in rows)
        buf.write("# lsdtt3 tests\n\n" + header_note + "\n\n")
        buf.write("Descriptions come from: " + ", ".join(
            f"{k} {v}" for k, v in sources.most_common()) + ". A *section* description "
            "is the nearest comment above the test and often describes a group of tests.\n\n")
        buf.write("| Test | Suite | Level | Evidence | Domains | Programs | Checks | Description |\n"
                  "|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            desc = short(r["description"]).replace("|", "\\|")
            if r["source"] == "section" and desc:
                desc = "*(section)* " + desc
            cells = [f"`{r['test']}`", r["suite"], r["level"], r["evidence"].replace(";", ", "),
                     r["domains"].replace(";", ", "), r["programs"],
                     r["checks"].replace("|", "\\|"), desc]
            buf.write("| " + " | ".join(cells) + " |\n")
    return buf.getvalue()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-dir", required=True, type=Path,
                    help="A configured lsdtt3 build directory (for the CTest tests).")
    ap.add_argument("--table", choices=("counts", "list"), default="counts")
    ap.add_argument("--format", choices=("md", "csv", "tex"), default="md")
    ap.add_argument("--suite", choices=("all", "ctest", "pytest"), default="all")
    ap.add_argument("-o", "--output", type=Path, help="Write here instead of to stdout.")
    ap.add_argument("--ctest", default="ctest")
    args = ap.parse_args(argv)

    if not (args.build_dir / "CTestTestfile.cmake").exists():
        print(f"ERROR {args.build_dir} is not a configured lsdtt3 build directory "
              "(no CTestTestfile.cmake). Run cmake there first.", file=sys.stderr)
        return 2
    rows = []
    if args.suite in ("all", "ctest"):
        rows += ctest_rows(args.build_dir, args.ctest)
    if args.suite in ("all", "pytest"):
        rows += pytest_rows()
    text = render(rows, args.table, args.format, configuration(args.build_dir))
    if args.output:
        args.output.write_text(text)
        print(f"wrote {args.output} ({len(rows)} tests)", file=sys.stderr)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
