#!/usr/bin/env python3
# check_python_test_markers.py
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Check the evidence and level pytest markers of python_bindings/tests.

Every test in the Python suite carries the same two kinds of label as the CTest
suite (tests/cmake/evidence_labels.cmake), written as pytest markers:

  EVIDENCE, one or more: how the test knows the right answer (its oracle).
  LEVEL, exactly one: unit or system.

This script reads each python_bindings/tests/test_*.py with `ast`, without
importing it, so it needs neither pytest nor a built lsdtt3 binding and runs in
every configure (CTest unit_python_test_markers, tier_python_stdlib). Each
`test*` function, and each `test*` method of a test class (a unittest.TestCase
subclass or a pytest `Test*` class), gets its markers the way pytest resolves
them, by walking the node chain:

  1. the module `pytestmark` (one mark or a list of marks);
  2. the class decorators, the class-level `pytestmark`, and those of any base
     class defined in the same file;
  3. the function or method decorators.

pytest's own marks (skip, skipif, xfail, parametrize, usefixtures,
filterwarnings, timeout) are ignored. Reported, and exit 1:

  - an item with no evidence marker;
  - an item without exactly one level marker;
  - a `pytest.mark.<x>` where x is neither a marker below nor a pytest mark
    (a misspelling such as `characterization`);
  - a `pytestmark` the checker cannot read as a mark or a list of marks.

python_bindings/tests/conftest.py imports EVIDENCE_MARKERS and LEVEL_MARKERS
from this file to register the markers, and repeats the first two checks at
collection time for parametrised and generated items. This file is the one
list; --cmake-labels confirms it equals the CMake lists.

Standard library only (tier_python_stdlib).

Usage:
    check_python_test_markers.py [--repo-root DIR] [--files F ...]
                                 [--counts] [--cmake-labels FILE]
    check_python_test_markers.py --self-test

Exit codes: 0 everything passes; 1 a violation (every offender is listed);
2 bad arguments or a file that cannot be read.
"""

from __future__ import annotations

import argparse
import ast
from dataclasses import dataclass, field
from pathlib import Path
import re
import sys

# The one list of marker names and definitions (plan section 1, and the header
# of tests/cmake/evidence_labels.cmake). Order matches LSDTT_EVIDENCE_LABELS and
# LSDTT_LEVEL_LABELS there.
EVIDENCE_MARKERS: dict[str, str] = {
    "exact_solution": "Output compared with an answer worked out analytically or "
                      "by hand from the fixture's construction (includes "
                      "hand-written expected values for parsers, formatters and "
                      "table lookups).",
    "property": "An assertion that holds for any correct output: no pits, "
                "conservation, weights sum to 1, monotonicity, range, "
                "idempotence, no nodata or NaN in an output, write/read round "
                "trip.",
    "metamorphic": "Two lsdtt3 runs on related inputs or settings whose outputs "
                   "must be related: rotation or roll, thread count, two methods "
                   "or code paths that must agree, two interfaces, scaling "
                   "relations, a model round trip.",
    "differential": "Comparison with an implementation lsdtt3 does not own, "
                    "including frozen outputs of those programs and a reference "
                    "written in the test that does not call the code under test.",
    "characterisation": "lsdtt3's own earlier output stored and compared: "
                        "checksums, sha256, frozen files, pinned measured "
                        "numbers, measured bands, or the replaced code kept "
                        "verbatim as the reference.",
    "smoke": "Runs to completion; in the Python suite, returns a well-formed "
             "result (shape, columns, not empty).",
    "negative": "Invalid input or configuration is refused with the expected "
                "message, exception or exit code.",
    "performance": "Timings reported, never asserted. Carries no evidence of "
                   "correctness.",
    "conformance": "The repository, help text, documentation, version strings, "
                   "naming conventions, warnings or test tooling conform to a "
                   "stated rule; no computed output is judged.",
}
LEVEL_MARKERS: dict[str, str] = {
    "unit": "Level: calls lsdtt3 routines in process, or inspects repository "
            "files or test tooling without running a shipped program.",
    "system": "Level: judges an output written by a shipped program run end to "
              "end (a CLI or console script, or a stored CLI output compared "
              "with a binding result).",
}
# pytest's own marks (and pytest-timeout's), which carry no evidence or level.
PYTEST_BUILTIN_MARKS = frozenset(
    {"skip", "skipif", "xfail", "parametrize", "usefixtures", "filterwarnings", "timeout"})

TEST_PREFIX = "test"   # pytest's python_functions default, and unittest's


# -----------------------------------------------------------------------------
# Reading marks from the AST
# -----------------------------------------------------------------------------
@dataclass
class Mark:
    name: str
    line: int


@dataclass
class FileResult:
    items: dict[str, list[str]] = field(default_factory=dict)   # id -> marker names
    problems: list[str] = field(default_factory=list)           # file-level problems


class _Unreadable(Exception):
    pass


def _dotted(node: ast.AST) -> str | None:
    parts: list[str] = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if isinstance(node, ast.Name):
        parts.append(node.id)
        return ".".join(reversed(parts))
    return None


def _mark_name(node: ast.AST, mark_prefixes: tuple[str, ...]) -> str | None:
    """The marker name of `pytest.mark.x` or `pytest.mark.x(...)`, else None."""
    if isinstance(node, ast.Call):
        node = node.func
    dotted = _dotted(node)
    if dotted is None:
        return None
    for prefix in mark_prefixes:
        if dotted.startswith(prefix) and dotted.count(".") == prefix.count("."):
            return dotted[len(prefix):]
    return None


def _marks_of_pytestmark(value: ast.AST, mark_prefixes: tuple[str, ...]) -> list[Mark]:
    elements = value.elts if isinstance(value, (ast.List, ast.Tuple)) else [value]
    marks = []
    for element in elements:
        name = _mark_name(element, mark_prefixes)
        if name is None:
            raise _Unreadable(element.lineno)
        marks.append(Mark(name, element.lineno))
    return marks


def _pytestmark_in(body: list[ast.stmt], mark_prefixes: tuple[str, ...],
                   where: str, problems: list[str], filename: str) -> list[Mark]:
    """Marks from `pytestmark = ...` (and `pytestmark += ...`) in a body."""
    marks: list[Mark] = []
    for stmt in body:
        value = None
        augmented = False
        if isinstance(stmt, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == "pytestmark" for t in stmt.targets):
            value = stmt.value
        elif (isinstance(stmt, ast.AnnAssign) and isinstance(stmt.target, ast.Name)
              and stmt.target.id == "pytestmark" and stmt.value is not None):
            value = stmt.value
        elif (isinstance(stmt, ast.AugAssign) and isinstance(stmt.target, ast.Name)
              and stmt.target.id == "pytestmark"):
            value, augmented = stmt.value, True
        if value is None:
            continue
        try:
            found = _marks_of_pytestmark(value, mark_prefixes)
        except _Unreadable as exc:
            problems.append(f"{filename}:{exc.args[0]}: {where} pytestmark is not a "
                            "pytest.mark or a list of them that this checker can read")
            continue
        marks = marks + found if augmented else found
    return marks


def _decorator_marks(node: ast.AST, mark_prefixes: tuple[str, ...]) -> list[Mark]:
    marks = []
    for decorator in getattr(node, "decorator_list", []):
        name = _mark_name(decorator, mark_prefixes)
        if name is not None:
            marks.append(Mark(name, decorator.lineno))
    return marks


def _mark_prefixes(tree: ast.Module) -> tuple[str, ...]:
    """`pytest.mark.`, plus `<alias>.` for `from pytest import mark [as alias]`
    and `<alias>.mark.` for `import pytest as <alias>`."""
    prefixes = {"pytest.mark."}
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module == "pytest":
            for alias in node.names:
                if alias.name == "mark":
                    prefixes.add((alias.asname or "mark") + ".")
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "pytest" and alias.asname:
                    prefixes.add(alias.asname + ".mark.")
    return tuple(sorted(prefixes))


def _is_testcase_base(base: ast.AST) -> bool:
    dotted = _dotted(base)
    return dotted in {"unittest.TestCase", "TestCase", "unittest.IsolatedAsyncioTestCase",
                      "IsolatedAsyncioTestCase"}


def check_source(source: str, filename: str) -> FileResult:
    """Every test item of one file with its resolved marker names."""
    result = FileResult()
    tree = ast.parse(source, filename=filename)
    prefixes = _mark_prefixes(tree)

    def note_unknown(marks: list[Mark]) -> None:
        for mark in marks:
            if (mark.name not in EVIDENCE_MARKERS and mark.name not in LEVEL_MARKERS
                    and mark.name not in PYTEST_BUILTIN_MARKS):
                message = (f"{filename}:{mark.line}: pytest.mark.{mark.name} is neither an "
                           "evidence or level marker nor a pytest mark")
                if message not in result.problems:
                    result.problems.append(message)

    module_marks = _pytestmark_in(tree.body, prefixes, "module", result.problems, filename)
    note_unknown(module_marks)

    classes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}
    own_class_marks: dict[str, list[Mark]] = {}
    for name, node in classes.items():
        marks = _decorator_marks(node, prefixes) + _pytestmark_in(
            node.body, prefixes, f"class {name}", result.problems, filename)
        note_unknown(marks)
        own_class_marks[name] = marks

    def class_chain(name: str, seen: frozenset[str] = frozenset()) -> tuple[list[Mark], bool]:
        """Marks of a class and its same-file bases; whether it is a TestCase."""
        node = classes[name]
        marks: list[Mark] = []
        is_testcase = False
        for base in node.bases:
            if _is_testcase_base(base):
                is_testcase = True
            elif isinstance(base, ast.Name) and base.id in classes and base.id not in seen:
                base_marks, base_is_testcase = class_chain(base.id, seen | {name})
                marks += base_marks
                is_testcase = is_testcase or base_is_testcase
        return marks + own_class_marks[name], is_testcase

    def add_item(item_id: str, marks: list[Mark]) -> None:
        result.items[item_id] = [mark.name for mark in marks]

    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name.startswith(TEST_PREFIX):
                marks = _decorator_marks(node, prefixes)
                note_unknown(marks)
                add_item(f"{filename}::{node.name}", module_marks + marks)
        elif isinstance(node, ast.ClassDef):
            marks, is_testcase = class_chain(node.name)
            if not is_testcase and not node.name.startswith("Test"):
                continue
            for member in node.body:
                if (isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef))
                        and member.name.startswith(TEST_PREFIX)):
                    method_marks = _decorator_marks(member, prefixes)
                    note_unknown(method_marks)
                    add_item(f"{filename}::{node.name}::{member.name}",
                             module_marks + marks + method_marks)
    return result


def item_violations(item_id: str, names: list[str]) -> list[str]:
    """What is wrong with one item's resolved markers (empty when nothing)."""
    out = []
    if not any(name in EVIDENCE_MARKERS for name in names):
        out.append(f"{item_id}: no evidence marker")
    levels = sorted({name for name in names if name in LEVEL_MARKERS})
    if len(levels) != 1:
        given = ", ".join(levels) if levels else "none"
        out.append(f"{item_id}: needs exactly one level marker (unit or system); has {given}")
    return out


# -----------------------------------------------------------------------------
# Cross-check against tests/cmake/evidence_labels.cmake
# -----------------------------------------------------------------------------
def cmake_list(text: str, variable: str) -> list[str] | None:
    text = re.sub(r"#[^\n]*", "", text)
    match = re.search(r"\bset\(\s*" + re.escape(variable) + r"\b([^)]*)\)", text)
    return match.group(1).split() if match else None


def cmake_problems(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    out = []
    for variable, ours, name in (("LSDTT_EVIDENCE_LABELS", EVIDENCE_MARKERS, "EVIDENCE_MARKERS"),
                                 ("LSDTT_LEVEL_LABELS", LEVEL_MARKERS, "LEVEL_MARKERS")):
        theirs = cmake_list(text, variable)
        if theirs is None:
            out.append(f"{path}: no set({variable} ...) found")
        elif theirs != list(ours):
            out.append(f"{path}: {variable} is [{' '.join(theirs)}] but {name} in "
                       f"{Path(__file__).name} is [{' '.join(ours)}]")
    return out


# -----------------------------------------------------------------------------
# Running over files
# -----------------------------------------------------------------------------
def check_files(paths: list[Path], root: Path) -> tuple[dict[str, list[str]], list[str]]:
    items: dict[str, list[str]] = {}
    problems: list[str] = []
    for path in paths:
        try:
            label = str(path.resolve().relative_to(root))
        except ValueError:
            label = str(path)
        result = check_source(path.read_text(encoding="utf-8"), label)
        items.update(result.items)
        problems += result.problems
    for item_id, names in items.items():
        problems += item_violations(item_id, names)
    return items, problems


def print_counts(items: dict[str, list[str]]) -> None:
    print(f"Functions and methods per marker ({len(items)} in total; a "
          "function with several evidence markers is counted under each):")
    for name in list(EVIDENCE_MARKERS) + list(LEVEL_MARKERS):
        count = sum(1 for names in items.values() if name in names)
        print(f"  {name:<17} {count}")


# -----------------------------------------------------------------------------
# Self-test: built-in snippets, each with the verdict it must produce
# -----------------------------------------------------------------------------
SELF_TEST_CASES: list[tuple[str, str, set[str]]] = [
    ("module mark reaches TestCase methods",
     "import unittest\nimport pytest\n"
     "pytestmark = [pytest.mark.conformance, pytest.mark.unit]\n"
     "class A(unittest.TestCase):\n"
     "    def test_one(self):\n        pass\n"
     "    def helper(self):\n        pass\n"
     "def test_two():\n    pass\n",
     set()),
    ("class mark plus method mark",
     "import unittest\nimport pytest\n"
     "@pytest.mark.unit\n"
     "class A(unittest.TestCase):\n"
     "    @pytest.mark.negative\n"
     "    def test_refused(self):\n        pass\n"
     "    def test_unmarked(self):\n        pass\n"
     "class B(A):\n"
     "    pytestmark = pytest.mark.property\n"
     "    def test_inherits_level(self):\n        pass\n",
     {"A::test_unmarked: no evidence marker"}),
    ("missing level",
     "import pytest\n"
     "@pytest.mark.exact_solution\n"
     "def test_x():\n    pass\n",
     {"test_x: needs exactly one level marker (unit or system); has none"}),
    ("two levels",
     "import pytest\n"
     "pytestmark = pytest.mark.unit\n"
     "@pytest.mark.metamorphic\n@pytest.mark.system\n"
     "def test_x():\n    pass\n",
     {"test_x: needs exactly one level marker (unit or system); has system, unit"}),
    ("misspelt marker",
     "import pytest\n"
     "@pytest.mark.characterization\n@pytest.mark.unit\n"
     "def test_x():\n    pass\n",
     {"line 2: pytest.mark.characterization is neither an evidence or level marker "
      "nor a pytest mark",
      "test_x: no evidence marker"}),
    ("skipif and parametrize are ignored",
     "import pytest\n"
     "pytestmark = [pytest.mark.skipif(True, reason='x'), pytest.mark.unit]\n"
     "@pytest.mark.parametrize('a', [1, 2])\n"
     "@pytest.mark.smoke\n"
     "def test_x(a):\n    pass\n"
     "@pytest.mark.skipif(False, reason='y')\n"
     "def test_y():\n    pass\n",
     {"test_y: no evidence marker"}),
    ("unreadable pytestmark",
     "import pytest\n"
     "MARKS = [pytest.mark.unit]\n"
     "pytestmark = MARKS\n"
     "@pytest.mark.smoke\n"
     "def test_x():\n    pass\n",
     {"line 3: module pytestmark is not a pytest.mark or a list of them that this "
      "checker can read",
      "test_x: needs exactly one level marker (unit or system); has none"}),
]


def run_self_test() -> int:
    failures = 0
    for title, source, expected in SELF_TEST_CASES:
        result = check_source(source, "snippet.py")
        got = set()
        for problem in result.problems:
            got.add(re.sub(r"^snippet\.py:(\d+): ", r"line \1: ", problem))
        for item_id, names in result.items.items():
            for violation in item_violations(item_id, names):
                got.add(violation.removeprefix("snippet.py::"))
        if got == expected:
            print(f"PASS  {title}")
        else:
            failures += 1
            print(f"FAIL  {title}")
            for line in sorted(expected - got):
                print(f"        expected, not reported: {line}")
            for line in sorted(got - expected):
                print(f"        reported, not expected: {line}")
    # The CMake cross-check on a canned file: one equal, one with a misspelling.
    good = ("set(LSDTT_EVIDENCE_LABELS\n    " + " ".join(EVIDENCE_MARKERS)
            + ")  # trailing comment\nset(LSDTT_LEVEL_LABELS " + " ".join(LEVEL_MARKERS) + ")\n")
    for title, text, want_problems in (
            ("cmake lists equal", good, 0),
            ("cmake list differs", good.replace("characterisation", "characterization"), 1)):
        got = sum(1 for variable, ours in (("LSDTT_EVIDENCE_LABELS", EVIDENCE_MARKERS),
                                           ("LSDTT_LEVEL_LABELS", LEVEL_MARKERS))
                  if cmake_list(text, variable) != list(ours))
        if got == want_problems:
            print(f"PASS  {title}")
        else:
            failures += 1
            print(f"FAIL  {title}: {got} lists differ, expected {want_problems}")
    total = len(SELF_TEST_CASES) + 2
    print(f"{total - failures} of {total} self-test cases passed")
    return 0 if failures == 0 else 1


def main(argv: list[str] | None = None) -> int:
    default_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Check the evidence and level pytest markers of python_bindings/tests.")
    parser.add_argument("--repo-root", type=Path, default=default_root,
                        help="repository root (default: two levels above this script)")
    parser.add_argument("--files", type=Path, nargs="+", metavar="FILE",
                        help="check only these test files (default: every "
                             "python_bindings/tests/test_*.py)")
    parser.add_argument("--counts", action="store_true",
                        help="print the number of functions and methods per marker")
    parser.add_argument("--cmake-labels", type=Path, metavar="FILE",
                        help="also check that the marker lists equal LSDTT_EVIDENCE_LABELS "
                             "and LSDTT_LEVEL_LABELS in this file "
                             "(tests/cmake/evidence_labels.cmake)")
    parser.add_argument("--self-test", action="store_true",
                        help="run the checker on built-in snippets and exit")
    args = parser.parse_args(argv)

    if args.self_test:
        if args.files or args.counts or args.cmake_labels:
            parser.error("--self-test takes no other option")
        return run_self_test()

    root = args.repo_root.resolve()
    if args.files:
        paths = list(args.files)
    else:
        tests_dir = root / "python_bindings" / "tests"
        if not tests_dir.is_dir():
            parser.error(f"no python_bindings/tests under {root}")
        paths = sorted(tests_dir.glob("test_*.py"))
    for path in paths:
        if not path.is_file():
            parser.error(f"not a file: {path}")
    if args.cmake_labels is not None and not args.cmake_labels.is_file():
        parser.error(f"not a file: {args.cmake_labels}")

    try:
        items, problems = check_files(paths, root)
    except SyntaxError as exc:
        print(f"ERROR: cannot parse {exc.filename}: {exc}", file=sys.stderr)
        return 2
    if args.cmake_labels is not None:
        problems = cmake_problems(args.cmake_labels) + problems

    if args.counts:
        print_counts(items)
    offenders = {p.split(": ", 1)[0] for p in problems}
    if problems:
        print(f"FAIL: {len(problems)} marker problems in {len(items)} test functions "
              f"and methods ({len(offenders)} offenders) across {len(paths)} files:")
        for problem in problems:
            print(f"  {problem}")
        print("Markers and their meanings: tests/scripts/check_python_test_markers.py "
              "and docs/developer_notes/testing_and_benchmarking.md.")
        return 1
    checked = " and the CMake label lists" if args.cmake_labels is not None else ""
    print(f"PASS: {len(items)} test functions and methods in {len(paths)} files have "
          f"evidence and exactly one level marker{checked}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
