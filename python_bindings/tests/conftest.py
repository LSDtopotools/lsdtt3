# conftest.py
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Evidence and level markers for python_bindings/tests.

Every test carries one or more EVIDENCE markers (how it knows the right answer)
and exactly one LEVEL marker (unit or system), the same names as the CTest
labels in tests/cmake/evidence_labels.cmake:

    python -m pytest -m metamorphic python_bindings/tests
    python -m pytest -m "negative and unit" python_bindings/tests

The names and definitions are imported from tests/scripts/check_python_test_markers.py,
the one list (CTest unit_python_test_markers checks it against the CMake lists
and checks every test file statically, without a built binding). There is no
inline copy to fall back on: every run of this suite (cibuildwheel's
`cd {project} && python -m pytest -q python_bindings/tests`, the Apptainer wheel
build, tests/scripts/check_python_bindings_smoke.sh) starts from a repository
tree, so a failed import means a broken tree and stops the run.

pytest_collection_modifyitems is the runtime guard: collection fails, naming
every offender, for an item with no evidence marker, not exactly one level
marker, or a marker that is neither of those nor a pytest mark. It also covers
parametrised and generated items, which the static check cannot see. There is
deliberately no switch to turn it off.
"""

from __future__ import annotations

from pathlib import Path
import sys

import pytest

_SCRIPTS = Path(__file__).resolve().parents[2] / "tests" / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
try:
    from check_python_test_markers import (  # noqa: E402
        EVIDENCE_MARKERS, LEVEL_MARKERS, PYTEST_BUILTIN_MARKS)
except ImportError as exc:
    raise pytest.UsageError(
        f"python_bindings/tests/conftest.py could not import the marker list from "
        f"{_SCRIPTS / 'check_python_test_markers.py'} ({exc}). Run the suite from an "
        "lsdtt3 repository tree.") from exc


def pytest_configure(config: pytest.Config) -> None:
    for name, definition in {**EVIDENCE_MARKERS, **LEVEL_MARKERS}.items():
        config.addinivalue_line("markers", f"{name}: {definition}")


def pytest_collection_modifyitems(session: pytest.Session, config: pytest.Config,
                                  items: list[pytest.Item]) -> None:
    offenders: list[str] = []
    for item in items:
        names = {mark.name for mark in item.iter_markers()}
        problems = []
        if not names & EVIDENCE_MARKERS.keys():
            problems.append("no evidence marker")
        levels = sorted(names & LEVEL_MARKERS.keys())
        if len(levels) != 1:
            problems.append("needs exactly one level marker (unit or system); has "
                            + (", ".join(levels) if levels else "none"))
        unknown = sorted(names - EVIDENCE_MARKERS.keys() - LEVEL_MARKERS.keys()
                         - PYTEST_BUILTIN_MARKS)
        if unknown:
            problems.append("unknown marker " + ", ".join(unknown))
        if problems:
            offenders.append(f"  {item.nodeid}: {'; '.join(problems)}")
    if offenders:
        raise pytest.UsageError(
            f"{len(offenders)} of {len(items)} collected tests lack the evidence and "
            "level markers (tests/scripts/check_python_test_markers.py lists them; "
            "docs/developer_notes/testing_and_benchmarking.md defines them):\n"
            + "\n".join(offenders))
