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

import os
import pathlib
import shlex
import subprocess
from typing import Any, Mapping, Sequence

NONINTERACTIVE_ENV = "LSDTT_NONINTERACTIVE"
GDAL_LIBRARY_ENV = "LSDTT_TEST_GDAL_LIBRARY"

# What an activated Pixi/conda environment sets that a binary built against a
# DIFFERENT GDAL must not see.
ENVIRONMENT_LIBRARY_VARS = (
    "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "GDAL_DATA", "GDAL_DRIVER_PATH",
    "PROJ_DATA", "PROJ_LIB", "CONDA_PREFIX",
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


def run(cmd: Sequence[Any] | str, **kwargs: Any) -> subprocess.CompletedProcess:
    """subprocess.run, non-interactive. Same arguments and return value."""
    return subprocess.run(cmd, **_prepare(cmd, kwargs))


def popen(cmd: Sequence[Any] | str, **kwargs: Any) -> subprocess.Popen:
    """subprocess.Popen, non-interactive. Same arguments and return value."""
    return subprocess.Popen(cmd, **_prepare(cmd, kwargs))


def check_output(cmd: Sequence[Any] | str, **kwargs: Any) -> Any:
    """subprocess.check_output, non-interactive. Same arguments and return value."""
    return subprocess.check_output(cmd, **_prepare(cmd, kwargs))


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
    return subprocess.CompletedProcess(argv, result.returncode, stdout, stderr)
