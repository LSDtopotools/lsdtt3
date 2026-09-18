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
"""

from __future__ import annotations

import os
import pathlib
import shlex
import subprocess
from typing import Any, Mapping, Sequence

NONINTERACTIVE_ENV = "LSDTT_NONINTERACTIVE"


def _noninteractive(kwargs: dict[str, Any]) -> dict[str, Any]:
    if "stdin" not in kwargs and "input" not in kwargs:
        kwargs["stdin"] = subprocess.DEVNULL
    base: Mapping[str, str] = kwargs.get("env") or os.environ
    env = dict(base)
    env.setdefault(NONINTERACTIVE_ENV, "1")
    kwargs["env"] = env
    return kwargs


def run(cmd: Sequence[Any] | str, **kwargs: Any) -> subprocess.CompletedProcess:
    """subprocess.run, non-interactive. Same arguments and return value."""
    return subprocess.run(cmd, **_noninteractive(kwargs))


def popen(cmd: Sequence[Any] | str, **kwargs: Any) -> subprocess.Popen:
    """subprocess.Popen, non-interactive. Same arguments and return value."""
    return subprocess.Popen(cmd, **_noninteractive(kwargs))


def check_output(cmd: Sequence[Any] | str, **kwargs: Any) -> Any:
    """subprocess.check_output, non-interactive. Same arguments and return value."""
    return subprocess.check_output(cmd, **_noninteractive(kwargs))


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
