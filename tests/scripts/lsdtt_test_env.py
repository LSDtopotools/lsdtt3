"""Shared process launching for the lsdtt3 test scripts.

Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3

Copyright (C) 2026 Simon M. Mudd
SPDX-License-Identifier: AGPL-3.0-or-later

Every test script starts its child processes through this module rather than
calling subprocess directly, for one reason: a test must never wait on the
terminal.

CTest hands the invoking terminal to each test as its stdin, and a child
started with subprocess inherits it. A CLI that prompts when stdin is a
terminal (lsdtt-fetch-raster asks for an OpenTopography key) then blocks
forever with its prompt swallowed by capture_output -- in an interactive ctest
run only. Agent and CI runs have no terminal, so they pass, and the hang is
invisible to them. That is how failure_fetch_raster_ot_missing_key_strict came
to hang a full run for 40 minutes on 2026-09-18.

So every launch here:

  * gives the child /dev/null as stdin unless the caller supplies stdin or
    input, and
  * sets LSDTT_NONINTERACTIVE=1 in the child's environment, which lsdtt3's
    credential lookup treats as "never prompt" (see lsdtt_secret.hpp). CTest
    sets the same variable on every registered test, so this matters only when
    a script is run by hand.

run_logged() additionally streams the child's output to files while it runs,
so a hung or killed test leaves its command and partial output on disk.

ISOLATION FROM AN ACTIVE PIXI / CONDA ENVIRONMENT

Tests are often run inside Pixi (`pixi run -e test ctest ...`) so that the
Python side has its packages. Pixi's activation sets GDAL_DATA, PROJ_DATA and
GDAL_DRIVER_PATH for ITS GDAL (3.12 on 2026-09-18). An lsdtt3 CLI built against
the SYSTEM GDAL (3.4.1 / PROJ 8.2.1 here) inherits them and silently reads
another installation's PROJ database and GDAL data -- no error, just possibly
different CRS definitions. Older Pixi releases also set LD_LIBRARY_PATH, which
made such a CLI load Pixi's libgdal and fail with `CXXABI_1.3.15 not found`
(recorded in the benchmark drivers; not reproducible with Pixi 0.70).

So an lsdtt3 child (argv[0] named lsdtt-*) gets system_env() -- those
variables removed -- whenever the build's GDAL is not the active environment's.
CTest tells every test which GDAL the build used (LSDTT_TEST_GDAL_LIBRARY);
when that is unset (a script run by hand) nothing is stripped. Other children,
such as the environment's own gdalinfo, keep the environment they need.
"""

from __future__ import annotations

import importlib
import os
import pathlib
import re
import shlex
import subprocess
import sys
from typing import Any, Mapping, Sequence

NONINTERACTIVE_ENV = "LSDTT_NONINTERACTIVE"
GDAL_LIBRARY_ENV = "LSDTT_TEST_GDAL_LIBRARY"

# What an activated Pixi/conda environment sets that a binary built against a
# DIFFERENT GDAL must not see. PROJ_NETWORK=ON (set by pixi) would let the
# system PROJ fetch transformation grids over the network mid-test.
ENVIRONMENT_LIBRARY_VARS = (
    "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "GDAL_DATA", "GDAL_DRIVER_PATH",
    "PROJ_DATA", "PROJ_LIB", "PROJ_NETWORK", "CONDA_PREFIX",
)


def system_env(base: Mapping[str, str] | None = None) -> dict[str, str]:
    """A copy of `base` (default os.environ) without the Pixi/conda library and
    data variables, for running an lsdtt3 binary built against the system GDAL.

    The one copy of this logic; the benchmark drivers call it too."""
    env = dict(os.environ if base is None else base)
    for key in ENVIRONMENT_LIBRARY_VARS:
        env.pop(key, None)
    return env


def build_gdal_is_foreign_to_env(env: Mapping[str, str]) -> bool:
    """True when an environment is active (CONDA_PREFIX) and the lsdtt3 build's
    GDAL (LSDTT_TEST_GDAL_LIBRARY, set by CTest) is not inside it."""
    prefix = env.get("CONDA_PREFIX", "")
    gdal_library = env.get(GDAL_LIBRARY_ENV, "")
    if not prefix or not gdal_library:
        return False
    prefix = os.path.realpath(prefix) + os.sep
    return not os.path.realpath(gdal_library).startswith(prefix)


def _is_lsdtt_binary(cmd: Sequence[Any] | str) -> bool:
    first = cmd if isinstance(cmd, str) else (cmd[0] if cmd else "")
    return os.path.basename(str(first)).startswith("lsdtt-")


def _prepare(cmd: Sequence[Any] | str, kwargs: dict[str, Any]) -> dict[str, Any]:
    if "stdin" not in kwargs and "input" not in kwargs:
        kwargs["stdin"] = subprocess.DEVNULL
    base: Mapping[str, str] = kwargs.get("env") or os.environ
    env = dict(base)
    env.setdefault(NONINTERACTIVE_ENV, "1")
    if _is_lsdtt_binary(cmd) and build_gdal_is_foreign_to_env(env):
        env = system_env(env)
    kwargs["env"] = env
    return kwargs


# ---------------------------------------------------------------------------
# Diagnostics: say what is missing and how to get it, instead of a traceback.
# ---------------------------------------------------------------------------
GDAL_PROGRAMS = ("gdalinfo", "ogrinfo", "gdal_translate")
PIXI_TEST_HINT = ("pixi install -e test, configure with -DPython3_EXECUTABLE="
                  "$(pixi run -e test python -c 'import sys; print(sys.executable)'), "
                  "and run: pixi run -e test ctest --test-dir <build> ...")


def _program_of(cmd: Sequence[Any] | str) -> str:
    first = cmd if isinstance(cmd, str) else (cmd[0] if cmd else "")
    return os.path.basename(str(first).split()[0]) if first else ""


def _missing_program(cmd: Sequence[Any] | str, exc: FileNotFoundError) -> FileNotFoundError:
    prog = _program_of(cmd)
    if prog in GDAL_PROGRAMS:
        msg = (f"{prog} was not found on PATH. This test runs GDAL's command-line "
               f"programs (label needs_gdal_tools). Install them (Ubuntu: sudo apt "
               f"install gdal-bin), or use the Pixi test environment: {PIXI_TEST_HINT}")
    else:
        msg = f"{prog or cmd!r} could not be started: {exc}"
    return FileNotFoundError(msg)


_LOADER_FAILURE = re.compile(
    r"(CXXABI|GLIBCXX)_[0-9.]+' not found|error while loading shared libraries")


def explain_failure(stderr: str) -> str:
    """A remediation line for a child that failed to LOAD rather than to run:
    the Pixi-versus-system library mix-up. Empty when the failure is anything
    else."""
    if not _LOADER_FAILURE.search(stderr or ""):
        return ""
    return ("The program failed to load its shared libraries, not to run. This is "
            "usually a binary built against the system GDAL picking up a Pixi/conda "
            "environment's libraries (LD_LIBRARY_PATH), or the reverse. Configure the "
            "C++ build OUTSIDE Pixi and run the tests with: pixi run -e test ctest ... "
            "(the launcher then isolates lsdtt3 binaries), or rebuild inside the "
            "environment the tests run in.")


def require_modules(modules: Sequence[str], tier: str = "tier_python_geo") -> None:
    """Exit with a clear message if any module cannot be imported. Call before
    the script's own third-party imports."""
    missing = []
    for name in modules:
        try:
            importlib.import_module(name)
        except ImportError:
            missing.append(name)
    if missing:
        print(f"FAIL  this test is {tier} and needs Python packages that "
              f"{sys.executable} cannot import: {', '.join(missing)}. "
              f"Use the Pixi test environment: {PIXI_TEST_HINT}", file=sys.stderr)
        sys.exit(1)


def run(cmd: Sequence[Any] | str, **kwargs: Any) -> subprocess.CompletedProcess:
    """subprocess.run, non-interactive. Same arguments and return value."""
    try:
        return subprocess.run(cmd, **_prepare(cmd, kwargs))
    except FileNotFoundError as exc:
        raise _missing_program(cmd, exc) from None


def popen(cmd: Sequence[Any] | str, **kwargs: Any) -> subprocess.Popen:
    """subprocess.Popen, non-interactive. Same arguments and return value."""
    try:
        return subprocess.Popen(cmd, **_prepare(cmd, kwargs))
    except FileNotFoundError as exc:
        raise _missing_program(cmd, exc) from None


def check_output(cmd: Sequence[Any] | str, **kwargs: Any) -> Any:
    """subprocess.check_output, non-interactive. Same arguments and return value."""
    try:
        return subprocess.check_output(cmd, **_prepare(cmd, kwargs))
    except FileNotFoundError as exc:
        raise _missing_program(cmd, exc) from None


def run_logged(cmd: Sequence[Any],
               work_dir: pathlib.Path,
               cwd: pathlib.Path | None = None,
               **kwargs: Any) -> subprocess.CompletedProcess:
    """Run cmd non-interactively, streaming its output to files as it runs.

    Writes, in work_dir:
      command.txt  argv (shell-quoted) and cwd, BEFORE the child starts
      stdout.log   the child's stdout, written while it runs
      stderr.log   the child's stderr, written while it runs

    Returns a CompletedProcess whose stdout and stderr are the decoded logs,
    so callers can assert on them exactly as they would after
    subprocess.run(capture_output=True, text=True).
    """
    work_dir = pathlib.Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    argv = [str(a) for a in cmd]
    run_dir = pathlib.Path(cwd) if cwd is not None else pathlib.Path.cwd()

    (work_dir / "command.txt").write_text(
        f"cwd: {run_dir}\n{shlex.join(argv)}\n", encoding="utf-8")

    out_path = work_dir / "stdout.log"
    err_path = work_dir / "stderr.log"
    with open(out_path, "wb") as out, open(err_path, "wb") as err:
        result = run(argv, cwd=str(run_dir), stdout=out, stderr=err,
                     check=False, **kwargs)

    stdout = out_path.read_text(encoding="utf-8", errors="replace")
    stderr = err_path.read_text(encoding="utf-8", errors="replace")
    if result.returncode != 0:
        hint = explain_failure(stderr)
        if hint:
            print(f"NOTE  {hint}", file=sys.stderr)
    return subprocess.CompletedProcess(argv, result.returncode, stdout, stderr)
