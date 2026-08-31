#!/usr/bin/env python3
"""Standalone cross-implementation benchmark suite for lsdtt3.

This is the entry point a user or a peer reviewer runs.  It is deliberately NOT
part of the test suite: `ctest` never invokes it, it asserts nothing, and it
cannot fail a build.  Its only job is to produce a defensible measurement of how
lsdtt3 compares against other topographic-analysis implementations on the
machine it is run on.

    pixi run bench-laptop                 # ~minutes, reduced fixtures
    pixi run bench-server                 # the peer-review run
    pixi run bench-laptop --dry-run       # show the plan, measure nothing

Both profiles write a single self-contained directory named after the host and
timestamp, containing the raw per-measurement CSV, a hardware manifest, a run
manifest recording every command, and a rendered summary.  Nothing outside that
directory is written, so a campaign can be tarred and sent as-is.


THE TWO PROFILES

  laptop  Reduced fixtures and thread counts.  Two purposes: it gives a laptop
          user a real baseline without a multi-hour run, and it is the smoke
          test -- run it before dispatching a server campaign to catch a broken
          path in minutes rather than after a queue wait.  Its numbers are
          honest but small-N; see OVERHEAD below for why the small fixtures are
          the ones you must not over-read.

  server  Full fixture ladder and the 1/2/4/8 thread sweep.  This is the
          profile whose output is quotable.

  Slurm is deliberately NOT handled here.  Cluster submission lives in
  tests/scripts/slurm_benchmark.sbatch, which does nothing but arrange for this
  script to be run; every machine without Slurm runs the identical command.


WHAT MAKES THIS DEFENSIBLE (and where it still is not)

  Threads are pinned for every engine, never left at the library default.
  Whitebox's `max_procs` defaults to -1 (all cores) and lsdtt3 will take every
  OpenMP thread available; an unpinned comparison measures core count.  Every
  row records the thread count it ran at.  RichDEM and TopoToolbox are serial
  for these routines and are only ever run at 1 thread -- their rows are
  labelled `serial_reference` so nobody reads a flat scaling curve as a finding.

  Repeats, not single shots.  Each measurement is a discarded warm-up followed
  by N timed repeats; the CSV carries best, median, stdev and the full list.
  The SUMMARY quotes the MEDIAN, not the best: best-of-N is the right statistic
  for "how fast can this go" and the wrong one for "what will I get".

  I/O is measured apart from compute.  The engines have very different I/O
  layers (GDAL, libtopotoolbox, Whitebox's own reader, lsdtt3's loader) and a
  bundled number compares file formats.  `seconds` is compute; `seconds_read`
  and `seconds_write` are separate columns.

  Memory is the child process peak RSS, measured with wait4(2) around the whole
  child, so it covers the C++/Rust allocation the process actually made rather
  than what a Python-level probe can see.

  Fixture identity is hashed.  Every row carries the fixture's SHA256 prefix.
  A campaign compared against another campaign on a different fixture build is
  comparing different rasters, and the hash makes that visible.

  *** WHERE IT IS NOT APPLES TO APPLES -- state this alongside any table ***

  The lsdtt3 numbers are whole-CLI-process wall times: process start, GDAL
  raster read, the routine, the write, and teardown.  The reference numbers are
  in-process timings around the routine alone.  lsdtt3 is therefore charged for
  startup that the references are not.  `seconds_process` and `seconds` are
  both recorded so the gap is visible, and the summary flags any row where
  process overhead is a large fraction of the total.  Reducing that gap would
  mean building a C++ timing harness per routine, which is a bigger change than
  this suite; until then, small-fixture lsdtt3 rows understate lsdtt3.


OVERHEAD -- WHY SMALL FIXTURES LIE, MEASURED NOT ASSUMED

  Two fixed costs dominate small grids and neither is an algorithm:

  1. lsdtt3 CLI process startup (fork, dynamic linking, GDAL driver
     registration).  Measured per campaign by `--calibrate` and written to the
     manifest as `lsdtt_startup_seconds`.

  2. Whitebox tool-handle resolution.  `wbe.hydrology.sink` is not an attribute
     read; it is a registry resolution through the Rust layer costing ~190 ms
     PER ACCESS on the development machine.  whitebox_reference.py binds the
     handle once, outside the timing loop, which took fill on 256^2 from
     0.207 s to 0.019 s -- a 10.9x artefact that had flattered lsdtt3.  If a
     Whitebox row ever comes back suspiciously flat at ~0.2 s regardless of
     grid size, that regression is back.

  Consequently any row whose median is within MIN_TRUSTWORTHY_FACTOR of the
  measured fixed overhead is marked `overhead_dominated` in the CSV and is
  excluded from the summary's ratio table.  On the laptop profile that will be
  most of the small-fixture rows.  This is the intended behaviour: the laptop
  profile proves the pipeline works, the server profile produces the numbers.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import pathlib
import platform
import shutil
import socket
import statistics
import subprocess
import sys
import time
from typing import Any

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "dem_cache"

SKIP_EXIT = 77

# A measurement is excluded from the quotable ratio table unless its median is
# at least this many times the measured fixed overhead.  3x is not a standard,
# it is a judgement: below it, a 30% change in startup moves the ratio by 10%.
MIN_TRUSTWORTHY_FACTOR = 3.0

# Reference engines that have no thread control for these routines.  They are
# run once at 1 thread and reused across the sweep rather than re-measured,
# because re-running a serial tool at --threads 8 measures nothing and triples
# the campaign length.
SERIAL_REFERENCES = ("richdem", "topotoolbox")

REFERENCE_SCRIPTS = {
    "richdem": "richdem_reference.py",
    "topotoolbox": "topotoolbox_reference.py",
    "whitebox": "whitebox_reference.py",
}


# -----------------------------------------------------------------------------
# Case registry
#
# One entry per comparable operation.  `engines` maps an engine name to the
# routine string that engine understands; an engine absent from the map simply
# produces no row, which is how "RichDEM has no slope tool" is expressed.
#
# `lsdtt` describes how to drive the lsdtt3 CLI: which executable, which
# parameters, and which output file to expect.  These mirror run_flow_benchmark
# .py's COMPARISONS table; they are restated rather than imported because that
# table is shaped for correctness assertions and this one for timing, and
# entangling them would mean a tolerance change silently altering a benchmark.
# -----------------------------------------------------------------------------
CASES: dict[str, dict[str, Any]] = {
    # -- suite: flow -------------------------------------------------------
    "fill": {
        "suite": "flow",
        "note": "depression fill, Barnes (2014) lineage, true flats",
        "engines": {"richdem": "fill", "topotoolbox": "fill",
                    "whitebox": "fill"},
        "lsdtt": {"exe": "lsdtt-raster-preprocessing",
                  "params": {"write_fill_raster": "true",
                             "fill_min_slope": "0"},
                  "output": "_fill.tif"},
    },
    "fill_epsilon": {
        "suite": "flow",
        "note": "fill with a gradient imposed across flats",
        "engines": {"richdem": "fill_epsilon", "whitebox": "fill_epsilon"},
        "lsdtt": {"exe": "lsdtt-raster-preprocessing",
                  "params": {"write_fill_raster": "true",
                             "fill_min_slope": "0.0001"},
                  "output": "_fill.tif"},
    },
    "breach": {
        "suite": "flow",
        "note": "least-cost breaching, Lindsay (2016)",
        "engines": {"richdem": "breach", "whitebox": "breach"},
        "lsdtt": {"exe": "lsdtt-raster-preprocessing",
                  "params": {"write_breached_raster": "true",
                             "breach_algorithm_choice": "lindsay2016"},
                  "output": "_breached.tif"},
    },
    "accum_d8": {
        "needs_filled": True,
        "suite": "flow",
        "note": "D8 accumulation, O'Callaghan & Mark (1984)",
        "engines": {"richdem": "accum_d8", "topotoolbox": "accum_d8",
                    "whitebox": "accum_d8"},
        "lsdtt": {"exe": "lsdtt-drainage-analysis",
                  "params": {"write_d8_drainage_area_raster": "true",
                             "raster_is_filled": "true"},
                  "output": "_d8_drainage_area.tif"},
    },

    # -- suite: mfd --------------------------------------------------------
    # Multi-direction routing. Every implementation uses a different
    # proportioning rule (Tarboton D-inf vs Quinn MFD vs Freeman FD8), so these
    # are TIMING rows. They are still worth measuring: the cost of spreading
    # flow across multiple receivers is the interesting quantity, not whether
    # the three agree on the split.
    "accum_dinf": {
        "needs_filled": True,
        "suite": "mfd",
        "note": "D-infinity, Tarboton (1997); TIMING ONLY, rules differ",
        "engines": {"richdem": "accum_dinf", "whitebox": "accum_dinf"},
        "lsdtt": {"exe": "lsdtt-drainage-analysis",
                  "params": {"write_dinf_drainage_area_raster": "true",
                             "raster_is_filled": "true"},
                  "output": "_dinf_drainage_area.tif"},
        "report_only": True,
    },
    "accum_mfd": {
        "needs_filled": True,
        "suite": "mfd",
        "note": "MFD: lsdtt3 Quinn (1991) vs Freeman FD8; TIMING ONLY",
        "engines": {"richdem": "accum_mfd", "whitebox": "accum_mfd"},
        "lsdtt": {"exe": "lsdtt-drainage-analysis",
                  "params": {"write_mfd_drainage_area_raster": "true",
                             "raster_is_filled": "true",
                             "mfd_exponent": "1.0"},
                  "output": "_mfd_drainage_area.tif"},
        "report_only": True,
    },

    # -- suite: surface ----------------------------------------------------
    # The three engines are architecturally different here and the case list is
    # shaped to make that visible rather than to hide it.
    #
    #   lsdtt3      fits a 6-term quadratic over a moving window of
    #               `polyfit_lambda_metres` and derives EVERY metric from the
    #               one fit.
    #   TopoToolbox runs scipy convolutions with the Evans 3x3 kernels, in
    #               Python, serial -- and recomputes all five convolutions per
    #               curvature type.
    #   Whitebox    runs 3x3 finite-difference operators in Rust, and likewise
    #               recomputes per curvature type.
    #
    # Two consequences, both stated in the summary rather than corrected for:
    #
    # 1. The window sizes differ. lsdtt3's default lambda of 70 m is ~2x the
    #    30 m grid spacing, a wider window than the references' fixed 3x3, so
    #    lsdtt3 is doing MORE work per cell. Its surface numbers are therefore
    #    conservative. Lambda is pinned below rather than left at the default
    #    so this is a stated parameter, not an inherited one.
    # 2. `polyfit_all` is the case where lsdtt3's one-fit-many-metrics design
    #    actually pays, and `slope_only` is the case where it does not. Both
    #    are measured; quoting only the first would be cherry-picking.
    #
    # Values are never compared: the curvature definitions and aspect
    # handedness genuinely differ between the three.
    "slope_only": {
        "suite": "surface",
        "note": ("slope alone. lsdtt3 pays for a full quadratic fit to produce "
                 "one metric, so this is its WORST surface case. TIMING ONLY"),
        "engines": {"topotoolbox": "evansslope", "whitebox": "slope"},
        "lsdtt": {"exe": "lsdtt-surface-metrics",
                  "params": {"polyfit_metrics": "slope",
                             "polyfit_method": "direct",
                             "polyfit_lambda_metres": "70.0"},
                  "output": "_direct_polyfit_slope.tif"},
        "report_only": True,
    },
    "polyfit_all": {
        "suite": "surface",
        "note": ("slope + aspect + profile + plan curvature. lsdtt3 derives "
                 "all four from ONE fit; the references recompute their "
                 "derivatives per metric. TIMING ONLY"),
        "engines": {"topotoolbox": "curvature_all4",
                    "whitebox": "curvature_all2"},
        "lsdtt": {"exe": "lsdtt-surface-metrics",
                  "params": {"polyfit_metrics":
                             "slope,aspect_degrees,profile_curvature,"
                             "plan_curvature",
                             "polyfit_method": "direct",
                             "polyfit_lambda_metres": "70.0"},
                  "output": "_direct_polyfit_slope.tif"},
        "report_only": True,
    },
    "hillshade": {
        "suite": "surface",
        "note": "hillshade rendering, not an algorithmic comparison. TIMING ONLY",
        "engines": {"topotoolbox": "hillshade_fused", "whitebox": "hillshade"},
        "lsdtt": {"exe": "lsdtt-surface-metrics",
                  "params": {"write_hillshade_raster": "true"},
                  "output": "_hillshade.tif"},
        "report_only": True,
    },

    # -- suite: depression -------------------------------------------------
    # lsdtt3's depression hierarchy builds a nested tree of depressions and
    # their spill relationships. Nothing else here does that, so these rows
    # compare lsdtt3 against the cheapest thing the other tools offer and the
    # ratio necessarily flatters them. Recorded because the ABSOLUTE cost of
    # the hierarchy is a number users ask for.
    "depression_hierarchy": {
        "suite": "depression",
        "note": ("lsdtt3 builds a nested depression tree; the references only "
                 "label or fill. TIMING ONLY, and NOT a fair ratio"),
        "engines": {"whitebox": "depression_sink"},
        "lsdtt": {"exe": "lsdtt-raster-preprocessing",
                  "params": {"write_depression_hierarchy_raster": "true"},
                  "output": "_depression_labels.tif"},
        "report_only": True,
    },
    "fill_wang_liu": {
        "suite": "depression",
        "note": "Wang & Liu (2006) fill lineage, a second fill implementation",
        "engines": {"whitebox": "depression_fill_wl"},
        "lsdtt": {"exe": "lsdtt-raster-preprocessing",
                  "params": {"write_fill_raster": "true",
                             "fill_min_slope": "0"},
                  "output": "_fill.tif"},
    },
}


# -----------------------------------------------------------------------------
# Profiles
# -----------------------------------------------------------------------------
PROFILES: dict[str, dict[str, Any]] = {
    "laptop": {
        "description": ("reduced fixtures and threads; a real baseline for a "
                        "laptop and the smoke test before a server campaign"),
        "fixtures": ["alpine_small", "lakeland_small"],
        "threads": [1, 2, 4],
        "repeats": 3,
        "suites": ["flow", "surface", "mfd", "depression"],
        # A laptop should not silently start a multi-gigabyte download.
        "allow_fetch": True,
        "max_fixture_px": 720,
    },
    "server": {
        "description": "full ladder and thread sweep; the quotable profile",
        "fixtures": ["alpine_small", "alpine_medium", "alpine_large",
                     "lakeland_small", "lakeland_medium"],
        "threads": [1, 2, 4, 8],
        "repeats": 5,
        "suites": ["flow", "surface", "mfd", "depression"],
        "allow_fetch": True,
        "max_fixture_px": 2400,
    },
}


def now_stamp() -> str:
    return dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def sha_prefix(path: pathlib.Path, n: int = 12) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:n]


def system_env(threads: int) -> dict[str, str]:
    """Environment for running lsdtt3 binaries from inside the pixi env.

    Two separate jobs, both load-bearing.

    1. Strip pixi's library and data paths.  The lsdtt3 CLIs are C++ binaries
       built against the SYSTEM GDAL.  Left alone they inherit pixi's
       LD_LIBRARY_PATH and load pixi's libgdal.so.38, which wants a newer
       libstdc++ than the system one they were linked against:

           libstdc++.so.6: version `CXXABI_1.3.15' not found

       This is lifted from run_flow_benchmark.py's system_env(), where it was
       found the hard way.

    2. Pin the thread count.  OMP_NUM_THREADS is what lsdtt3's OpenMP loops
       read.  OPENBLAS_NUM_THREADS and MKL_NUM_THREADS are set alongside it
       because a BLAS pulled in transitively will otherwise spin up its own
       pool and the "1 thread" row will not be one thread.
    """
    env = dict(os.environ)
    for key in ("LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "GDAL_DATA",
                "GDAL_DRIVER_PATH", "PROJ_DATA", "PROJ_LIB", "CONDA_PREFIX"):
        env.pop(key, None)
    env["OMP_NUM_THREADS"] = str(threads)
    env["OPENBLAS_NUM_THREADS"] = str(threads)
    env["MKL_NUM_THREADS"] = str(threads)
    return env


def run_with_rusage(cmd: list[str], env: dict[str, str] | None = None,
                    cwd: pathlib.Path | None = None) -> dict[str, Any]:
    """Run a command, returning output, wall time and the child's peak RSS.

    Uses wait4(2) so the memory figure is the kernel's accounting for the whole
    child process.  That is the only way to see a C++ or Rust allocation: a
    Python-level probe in the parent cannot, and getrusage(RUSAGE_CHILDREN) in
    the parent is a high-water mark across ALL children, so differencing it
    across calls attributes one child's peak to another.

    subprocess's own helpers reap the child, which consumes the rusage before
    it can be read, so the pipes are drained manually first.  stderr goes to a
    temp file rather than a second pipe: reading two pipes sequentially
    deadlocks the moment the unread one fills its buffer.

    Falls back to subprocess.run with peak_rss_mb=None where wait4 is
    unavailable (Windows).  A missing memory column is better than a wrong one.
    """
    import tempfile

    if not hasattr(os, "wait4"):
        t0 = time.perf_counter()
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              env=env, cwd=str(cwd) if cwd else None)
        return {"returncode": proc.returncode, "stdout": proc.stdout,
                "stderr": proc.stderr,
                "wall_seconds": time.perf_counter() - t0,
                "peak_rss_mb": None}

    with tempfile.TemporaryFile(mode="w+") as errf:
        t0 = time.perf_counter()
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errf,
                                text=True, env=env,
                                cwd=str(cwd) if cwd else None)
        stdout = proc.stdout.read() if proc.stdout else ""
        if proc.stdout:
            proc.stdout.close()
        pid, status, ru = os.wait4(proc.pid, 0)
        wall = time.perf_counter() - t0
        # Tell Popen the child is already reaped so its destructor does not
        # warn about an abandoned process.
        proc.returncode = os.waitstatus_to_exitcode(status)
        errf.seek(0)
        stderr = errf.read()

    return {
        "returncode": proc.returncode,
        "stdout": stdout,
        "stderr": stderr,
        "wall_seconds": wall,
        # ru_maxrss is kibibytes on Linux, bytes on macOS. Linux is where the
        # server campaign runs; on macOS this reads ~1024x high, which is
        # obvious rather than subtle.
        "peak_rss_mb": ru.ru_maxrss / 1024.0,
    }


# -----------------------------------------------------------------------------
# Fixtures
# -----------------------------------------------------------------------------
def ensure_fixture(key: str, fetch_exe: pathlib.Path | None,
                   allow_fetch: bool) -> pathlib.Path | None:
    """Return the fixture path, fetching it if allowed and missing.

    Fetching is delegated to run_flow_benchmark.ensure_fixture rather than
    reimplemented: it already encodes the bounding-box maths, the COP30 source
    and the reprojection settings, and a second copy of that would drift and
    silently benchmark a different raster.
    """
    path = FIXTURE_DIR / f"{key}.tif"
    if path.exists():
        return path
    if not allow_fetch or fetch_exe is None or not fetch_exe.exists():
        return None

    sys.path.insert(0, str(SCRIPT_DIR))
    try:
        import run_flow_benchmark as rfb
    except ImportError:
        return None
    try:
        return rfb.ensure_fixture(key, fetch_exe)
    except SystemExit:
        # ensure_fixture calls die_skip (SystemExit 77) when a source is
        # unreachable. A fixture we cannot fetch is a smaller campaign, not a
        # failed one.
        return None


# -----------------------------------------------------------------------------
# Measurement
# -----------------------------------------------------------------------------
def prepare_filled_fixture(fixture: pathlib.Path, work_dir: pathlib.Path,
                           build_dir: pathlib.Path) -> pathlib.Path | None:
    """Produce the filled DEM that every accumulation case is measured on.

    Flow accumulation cannot run on a raw DEM: lsdtt-drainage-analysis aborts
    with "could not assign a receiver ... unresolved sink or flat", and the
    other engines would each silently apply their OWN depression handling
    first.  That would make an accumulation comparison a comparison of fills.

    So the fill is done ONCE, and every engine accumulates over the identical
    raster.  lsdtt3's fill is used as that common input, which needs
    justifying since it is one of the tools under test: it is bit-exact against
    both RichDEM's Barnes 2014 fill and TopoToolbox's `fillsinks` (see
    docs/testing_and_benchmarking.md), so choosing it favours nobody.  It is
    NOT Whitebox's fill, which returns a different surface on the test fixture
    -- see the note in whitebox_reference.py.

    The fill itself is not timed here; the `fill` case measures that.
    """
    out_dir = work_dir / "filled"
    out_dir.mkdir(parents=True, exist_ok=True)
    filled = out_dir / f"{fixture.stem}_filled.tif"
    if filled.exists():
        return filled

    exe = build_dir / "lsdtt-raster-preprocessing"
    if not exe.exists():
        return None

    params = {
        "read_path": str(fixture.parent),
        "read_fname": fixture.name,
        "write_path": str(out_dir),
        "write_fname": fixture.stem,
        "enforce_projection": "false",
        "write_output_format": "gtiff",
        "write_fill_raster": "true",
        # NOT min_slope 0. A fill with true flats is unroutable by D8 --
        # lsdtt-drainage-analysis aborts with "could not assign a receiver ...
        # unresolved sink or flat" -- because a flat cell has no downslope
        # neighbour to pick. The epsilon gradient is what makes the filled
        # surface routable, and every engine accumulates over this same raster,
        # so it advantages none of them.
        "fill_min_slope": "0.0001",
    }
    param_path = out_dir / f"{fixture.stem}_fill.param"
    param_path.write_text("".join(f"{k}: {v}\n" for k, v in params.items()))

    r = run_with_rusage([str(exe), str(out_dir), param_path.name],
                        env=system_env(1))
    if r["returncode"] != 0:
        return None

    produced = out_dir / f"{fixture.stem}_fill.tif"
    if not produced.exists():
        return None
    produced.rename(filled)
    return filled


COMPUTE_TIME_RE = None


def parse_reported_times(stdout: str) -> dict[str, float | None]:
    """Pull lsdtt3's OWN reported compute and write times out of its stdout.

    Some lsdtt3 CLIs print, per output:

        Compute time: 0.00264887 s
        Write time: 0.0089856 s

    Where they do, that compute figure is directly comparable with the
    references' in-process timings, and it is far better than whole-process
    wall time for a routine that takes milliseconds against ~80 ms of process
    startup.  Times are SUMMED across the outputs a run produced, because a
    multi-metric run prints one pair per raster and the total is what that run
    cost.

    Returns None for a CLI that prints nothing -- currently
    lsdtt-raster-preprocessing and lsdtt-drainage-analysis, which is why
    `timing_scope` is recorded per row rather than assumed for all lsdtt3 rows.
    """
    import re

    global COMPUTE_TIME_RE
    if COMPUTE_TIME_RE is None:
        COMPUTE_TIME_RE = re.compile(
            r"^\s*(Compute|Write) time:\s*([0-9.eE+-]+)\s*s\s*$",
            re.MULTILINE)

    compute, write = [], []
    for kind, value in COMPUTE_TIME_RE.findall(stdout):
        try:
            v = float(value)
        except ValueError:
            continue
        (compute if kind == "Compute" else write).append(v)

    return {
        "compute": sum(compute) if compute else None,
        "write": sum(write) if write else None,
    }


def measure_lsdtt(case_name: str, case: dict, fixture: pathlib.Path,
                  work_dir: pathlib.Path, threads: int, repeats: int,
                  build_dir: pathlib.Path) -> dict[str, Any] | None:
    """Time the lsdtt3 CLI for one case.

    Whole-process wall time, repeated.  See the module docstring: this includes
    startup and I/O that the in-process reference timings exclude, and the
    campaign's measured startup cost is recorded so the gap can be reasoned
    about rather than ignored.
    """
    spec = case["lsdtt"]
    exe = build_dir / spec["exe"]
    if not exe.exists():
        return None

    out_dir = work_dir / f"lsdtt_{case_name}_{fixture.stem}_t{threads}"
    out_dir.mkdir(parents=True, exist_ok=True)

    params = {
        "read_path": str(fixture.parent),
        "read_fname": fixture.name,
        "write_path": str(out_dir),
        "write_fname": fixture.stem,
        # The fixture is already in metric UTM. Letting the CLI reproject would
        # score the engines on different rasters.
        "enforce_projection": "false",
        "write_output_format": "gtiff",
    }
    params.update(spec["params"])
    param_path = out_dir / "benchmark.param"
    param_path.write_text("".join(f"{k}: {v}\n" for k, v in params.items()))

    env = system_env(threads)
    # The CLI takes <directory> <filename>, not a single path.
    cmd = [str(exe), str(out_dir), param_path.name]

    r = run_with_rusage(cmd, env=env)          # warm-up, discarded
    if r["returncode"] != 0:
        return {"error": f"exit {r['returncode']}",
                "stderr": r["stderr"][-1500:], "stdout": r["stdout"][-1500:]}

    timings, rss, computes, writes = [], [], [], []
    for _ in range(max(1, repeats)):
        r = run_with_rusage(cmd, env=env)
        if r["returncode"] != 0:
            return {"error": f"exit {r['returncode']}",
                    "stderr": r["stderr"][-1500:]}
        timings.append(r["wall_seconds"])
        if r["peak_rss_mb"] is not None:
            rss.append(r["peak_rss_mb"])
        reported = parse_reported_times(r["stdout"])
        if reported["compute"] is not None:
            computes.append(reported["compute"])
        if reported["write"] is not None:
            writes.append(reported["write"])

    info = {
        "repeats": len(timings),
        "peak_rss_mb": max(rss) if rss else None,
        "seconds_process": min(timings),
        "seconds_process_median": statistics.median(timings),
    }

    if computes:
        # The CLI told us what the routine itself cost. Prefer it: it is the
        # figure directly comparable with the references' in-process timings,
        # and it excludes ~80 ms of process startup that would otherwise
        # dominate every small-fixture row.
        info.update({
            "seconds": min(computes),
            "seconds_median": statistics.median(computes),
            "seconds_all": computes,
            "seconds_stdev": (statistics.stdev(computes)
                              if len(computes) > 1 else 0.0),
            "seconds_write": min(writes) if writes else None,
            "timing_scope": "cli_reported_compute",
        })
    else:
        # This CLI prints no timing, so whole-process wall time is all there
        # is. Labelled as such so a reader is never left guessing which
        # convention a row uses, and so these rows are not silently compared
        # against in-process reference numbers as though they matched.
        info.update({
            "seconds": min(timings),
            "seconds_median": statistics.median(timings),
            "seconds_all": timings,
            "seconds_stdev": (statistics.stdev(timings)
                              if len(timings) > 1 else 0.0),
            "timing_scope": "whole_process",
        })
    return info


def measure_reference(engine: str, routine: str, fixture: pathlib.Path,
                      threads: int, repeats: int,
                      python_exe: str) -> dict[str, Any] | None:
    """Time one reference implementation via its standalone script.

    Each reference is a subprocess so that a missing optional package is a
    clean skip rather than an ImportError in this driver, and so that the
    child's peak RSS is separately attributable.
    """
    script = SCRIPT_DIR / REFERENCE_SCRIPTS[engine]
    if not script.exists():
        return None

    cmd = [python_exe, str(script), "--routine", routine,
           "--in", str(fixture), "--repeats", str(repeats)]
    # Only whitebox_reference.py takes --threads; the other two are serial and
    # rejecting an unknown flag is how argparse tells us we got that wrong.
    if engine == "whitebox":
        cmd += ["--threads", str(threads)]

    r = run_with_rusage(cmd, env=dict(os.environ))
    if r["returncode"] == SKIP_EXIT:
        return {"skipped": True, "reason": r["stderr"].strip()[-400:]}
    if r["returncode"] != 0:
        return {"error": f"exit {r['returncode']}",
                "stderr": r["stderr"][-1500:]}

    lines = [ln for ln in r["stdout"].strip().splitlines() if ln.strip()]
    if not lines:
        return {"error": "no JSON on stdout", "stderr": r["stderr"][-800:]}
    try:
        info = json.loads(lines[-1])
    except json.JSONDecodeError as exc:
        return {"error": f"unparseable JSON: {exc}",
                "stdout": r["stdout"][-800:]}

    info["peak_rss_mb"] = r["peak_rss_mb"]
    info["seconds_process"] = r["wall_seconds"]
    info["timing_scope"] = "in_process_routine"
    return info


def calibrate_overheads(build_dir: pathlib.Path, fixture: pathlib.Path,
                        python_exe: str) -> dict[str, Any]:
    """Measure the fixed costs that make small-fixture rows untrustworthy.

    Two numbers, both measured rather than assumed:

      lsdtt_startup_seconds  -- an lsdtt3 CLI invoked with no work to do. This
                                is process start, dynamic linking and GDAL
                                driver registration, which every lsdtt3 row
                                pays and no reference row does.

      whitebox_dispatch_seconds -- one Whitebox routine on the smallest
                                fixture. If whitebox_reference.py ever
                                regresses to resolving tool handles inside the
                                timing loop this jumps to ~0.2 s and every
                                Whitebox row becomes meaningless; recording it
                                makes that regression visible in the manifest
                                rather than silent in the data.
    """
    out: dict[str, Any] = {}

    exe = build_dir / "lsdtt-raster-preprocessing"
    if exe.exists():
        # No parameter file: the CLI prints usage and exits. That is exactly
        # the fixed cost we want, with none of the work.
        times = []
        for _ in range(5):
            r = run_with_rusage([str(exe)], env=system_env(1))
            times.append(r["wall_seconds"])
        out["lsdtt_startup_seconds"] = min(times)
        out["lsdtt_startup_median"] = statistics.median(times)

    info = measure_reference("whitebox", "depression_sink", fixture, 1, 3,
                             python_exe)
    if info and "seconds" in info:
        out["whitebox_dispatch_seconds"] = info["seconds"]
        out["whitebox_dispatch_note"] = (
            "REGRESSION TRIPWIRE, not a fixed cost: one real routine on the "
            "smallest fixture, so most of it is genuine work. A few tens of ms "
            "is correct; ~0.2 s regardless of grid size means "
            "whitebox_reference.py is resolving tool handles inside the "
            "timing loop again (see its _make_call docstring)")
    return out


# -----------------------------------------------------------------------------
# Hardware manifest
# -----------------------------------------------------------------------------
def run_capture(cmd: list[str]) -> dict[str, Any]:
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except (FileNotFoundError, PermissionError) as exc:
        return {"ok": False, "error": str(exc)}
    return {"ok": proc.returncode == 0, "stdout": proc.stdout,
            "stderr": proc.stderr}


def read_mem_total_kb() -> int | None:
    path = pathlib.Path("/proc/meminfo")
    if not path.exists():
        return None
    try:
        for line in path.read_text().splitlines():
            if line.startswith("MemTotal:"):
                return int(line.split()[1])
    except (OSError, ValueError, IndexError):
        return None
    return None


def detect_cpu_governors() -> dict[str, Any]:
    """CPU frequency governors, one of the few things that silently invalidates
    a timing campaign.

    A `powersave` governor or an active thermal cap makes a long campaign drift
    slower as it runs, which reads as poor scaling at high thread counts. This
    does not correct for it -- it records it, so a reviewer can see that the
    machine was or was not pinned to `performance`.
    """
    base = pathlib.Path("/sys/devices/system/cpu")
    governors: dict[str, int] = {}
    try:
        for child in base.glob("cpu[0-9]*/cpufreq/scaling_governor"):
            txt = child.read_text().strip()
            governors[txt] = governors.get(txt, 0) + 1
    except OSError:
        return {"available": False}
    if not governors:
        return {"available": False}
    return {"available": True, "counts": governors}


def package_versions(python_exe: str) -> dict[str, Any]:
    """Versions of the reference implementations, from the environment that
    will actually run them.

    A timing table is meaningless without these: "lsdtt3 is 3x faster than
    Whitebox" is a claim about a specific Whitebox version, and 2.x is
    explicitly described upstream as still moving.
    """
    code = (
        "import json,sys;out={'python':sys.version.split()[0]}\n"
        "for name,mod in (('richdem','richdem'),('topotoolbox','topotoolbox'),"
        "('whitebox_workflows','whitebox_workflows'),('numpy','numpy'),"
        "('scipy','scipy')):\n"
        "    try:\n"
        "        import importlib.metadata as md\n"
        "        out[name]=md.version(name)\n"
        "    except Exception as e:\n"
        "        out[name]=None\n"
        "print(json.dumps(out))\n"
    )
    r = run_capture([python_exe, "-c", code])
    if not r.get("ok"):
        return {"available": False, "error": r.get("stderr", "")[-400:]}
    try:
        return {"available": True, **json.loads(r["stdout"])}
    except json.JSONDecodeError:
        return {"available": False, "raw": r["stdout"][-400:]}


def git_context() -> dict[str, Any]:
    sha = run_capture(["git", "-C", str(REPO_ROOT), "rev-parse", "--short",
                       "HEAD"])
    dirty = run_capture(["git", "-C", str(REPO_ROOT), "status", "--porcelain"])
    return {
        "sha": sha["stdout"].strip() if sha.get("ok") else "unknown",
        # A dirty tree means the binaries may not match the recorded SHA. That
        # is fine for a local run and disqualifying for a published one.
        "dirty": bool(dirty["stdout"].strip()) if dirty.get("ok") else None,
    }


def read_cmake_cache(build_dir: pathlib.Path) -> dict[str, Any]:
    """Compiler, flags and build type.

    A benchmark of a Debug build is a benchmark of assertions. This is recorded
    so `CMAKE_BUILD_TYPE=Debug` in a published campaign is visible rather than
    inferred from suspiciously slow numbers.
    """
    cache = build_dir / "CMakeCache.txt"
    if not cache.exists():
        return {"available": False}
    # CMAKE_CXX_COMPILER_ID and _VERSION are NOT written to CMakeCache.txt --
    # they live in CMakeFiles/*/CMakeCXXCompiler.cmake -- so asking for them
    # here yields nothing and the summary renders "compiler ? ?". The compiler
    # path and the flags are in the cache and are what actually get recorded;
    # the version is recovered separately by compiler_version().
    wanted = {"CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS",
              "CMAKE_CXX_FLAGS_RELEASE"}
    out: dict[str, Any] = {"available": True}
    for line in cache.read_text(errors="ignore").splitlines():
        if "=" not in line or ":" not in line or line.startswith(("//", "#")):
            continue
        key_type, value = line.split("=", 1)
        key = key_type.split(":", 1)[0]
        if key in wanted:
            out[key] = value

    # Ask the compiler itself rather than guessing from the path: `c++` is a
    # symlink and its version is the thing that changes optimisation between
    # campaigns.
    exe = out.get("CMAKE_CXX_COMPILER")
    if exe:
        v = run_capture([exe, "--version"])
        if v.get("ok"):
            out["CXX_VERSION"] = v["stdout"].splitlines()[0].strip()
    return out


def capture_hardware_manifest(build_dir: pathlib.Path,
                              python_exe: str) -> dict[str, Any]:
    return {
        "captured_at": dt.datetime.now().isoformat(timespec="seconds"),
        "host": {
            "hostname": socket.gethostname(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "python": sys.version,
        },
        "cpu": {
            "logical_cores": os.cpu_count(),
            "lscpu": run_capture(["lscpu"]).get("stdout", ""),
            "governors": detect_cpu_governors(),
        },
        "memory": {"mem_total_kb": read_mem_total_kb()},
        "toolchain": {
            "cmake_cache": read_cmake_cache(build_dir),
            "git": git_context(),
            "packages": package_versions(python_exe),
        },
    }


# -----------------------------------------------------------------------------
# CSV
# -----------------------------------------------------------------------------
CSV_FIELDS = [
    "suite", "case", "engine", "fixture", "fixture_sha", "n_cells",
    "threads", "repeats",
    "seconds_median", "seconds_best", "seconds_stdev", "seconds_all",
    "seconds_read", "seconds_write", "seconds_process", "timing_scope",
    "peak_rss_mb", "mcells_per_second",
    "report_only", "serial_reference", "overhead_dominated",
    "status", "note",
]


def n_cells_of(path: pathlib.Path) -> int | None:
    try:
        from osgeo import gdal
        gdal.UseExceptions()
        ds = gdal.Open(str(path))
        return int(ds.RasterXSize) * int(ds.RasterYSize)
    except Exception:
        return None


def build_row(suite, case_name, engine, fixture_key, fixture_sha, n_cells,
              threads, info, case, overhead) -> dict[str, Any]:
    row = {f: "" for f in CSV_FIELDS}
    row.update({
        "suite": suite, "case": case_name, "engine": engine,
        "fixture": fixture_key, "fixture_sha": fixture_sha,
        "n_cells": n_cells if n_cells is not None else "",
        "threads": threads,
        "report_only": bool(case.get("report_only", False)),
        "serial_reference": engine in SERIAL_REFERENCES,
        "note": case.get("note", ""),
    })

    if info is None:
        row["status"] = "unavailable"
        return row
    if info.get("skipped"):
        row["status"] = "skipped"
        row["note"] = info.get("reason", "")[:200]
        return row
    if info.get("error"):
        row["status"] = "error"
        row["note"] = str(info.get("error"))[:200]
        return row

    median = info.get("seconds_median", info.get("seconds"))
    row.update({
        "status": "ok",
        "repeats": info.get("repeats", ""),
        "seconds_median": f"{median:.9g}" if median is not None else "",
        "seconds_best": f"{info['seconds']:.9g}" if "seconds" in info else "",
        "seconds_stdev": f"{info.get('seconds_stdev', 0.0):.9g}",
        "seconds_all": json.dumps([round(v, 9)
                                   for v in info.get("seconds_all", [])]),
        "seconds_read": (f"{info['seconds_read']:.9g}"
                         if info.get("seconds_read") is not None else ""),
        "seconds_write": (f"{info['seconds_write']:.9g}"
                          if info.get("seconds_write") is not None else ""),
        "seconds_process": (f"{info['seconds_process']:.9g}"
                            if info.get("seconds_process") is not None else ""),
        "timing_scope": info.get("timing_scope", ""),
        "peak_rss_mb": (f"{info['peak_rss_mb']:.2f}"
                        if info.get("peak_rss_mb") is not None else ""),
    })

    if median and n_cells:
        row["mcells_per_second"] = f"{(n_cells / median) / 1e6:.6g}"

    # Mark, do not drop. A row excluded from the ratio table is still evidence
    # that the pipeline ran; deleting it would hide a working measurement.
    # The overhead floor applies to exactly one thing: an lsdtt3 row whose
    # time is whole-process wall clock, and therefore genuinely contains ~80 ms
    # of process startup that no reference row pays.
    #
    # It deliberately does NOT apply to:
    #   - lsdtt3 rows timed by the CLI's own reported compute time, which
    #     exclude startup already;
    #   - whitebox rows. `whitebox_dispatch_seconds` is a REGRESSION TRIPWIRE
    #     for the tool-handle bug, not a fixed cost: it is measured by running
    #     a real routine on a real grid, so most of it is genuine work. Using
    #     it as a floor flagged every Whitebox row as overhead-dominated and
    #     emptied the ratio table.
    #   - richdem and topotoolbox rows, which are timed in-process with no
    #     comparable fixed floor measured for them.
    if engine == "lsdtt3" and info.get("timing_scope") == "whole_process":
        floor = overhead.get("lsdtt_startup_seconds")
    else:
        floor = None
    if floor and median:
        row["overhead_dominated"] = bool(median < MIN_TRUSTWORTHY_FACTOR * floor)
    else:
        row["overhead_dominated"] = False
    return row


# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
def write_summary(path: pathlib.Path, rows: list[dict], manifest: dict,
                  overhead: dict) -> None:
    """Render the campaign as Markdown.

    Quotes the MEDIAN, never the best-of-N: best-of-N answers "how fast can
    this go" and a user wants "what will I get". Rows flagged report_only or
    overhead_dominated are excluded from the ratio table and listed separately,
    so an excluded row is visibly excluded rather than quietly missing.
    """
    hw = manifest["hardware"]["host"]
    cpu = manifest["hardware"]["cpu"]
    pkgs = manifest["hardware"]["toolchain"].get("packages", {})
    cmake = manifest["hardware"]["toolchain"].get("cmake_cache", {})

    L: list[str] = []
    L.append(f"# lsdtt3 benchmark campaign - {manifest['profile']}")
    L.append("")
    L.append(f"- host: `{hw['hostname']}` ({hw['machine']}, "
             f"{cpu.get('logical_cores')} logical cores)")
    L.append(f"- platform: {hw['platform']}")
    L.append(f"- lsdtt3 git: `{manifest['hardware']['toolchain']['git']['sha']}`"
             f"{' **(dirty tree)**' if manifest['hardware']['toolchain']['git'].get('dirty') else ''}")
    L.append(f"- build type: `{cmake.get('CMAKE_BUILD_TYPE', 'unknown')}`, "
             f"compiler `{cmake.get('CXX_VERSION', cmake.get('CMAKE_CXX_COMPILER', 'unknown'))}`")
    if cmake.get("CMAKE_BUILD_TYPE", "").lower() not in ("release",
                                                         "relwithdebinfo"):
        L.append("")
        L.append(f"> **Build type is "
                 f"`{cmake.get('CMAKE_BUILD_TYPE') or 'unset'}`, not Release.** "
                 "These numbers measure a non-optimised build and must not be "
                 "quoted as lsdtt3's performance.")
    L.append(f"- references: richdem `{pkgs.get('richdem')}`, "
             f"topotoolbox `{pkgs.get('topotoolbox')}`, "
             f"whitebox_workflows `{pkgs.get('whitebox_workflows')}`")
    L.append(f"- started: {manifest['created_at']}")
    L.append("")

    gov = cpu.get("governors", {})
    if gov.get("available") and set(gov.get("counts", {})) - {"performance"}:
        L.append(f"> **CPU governor is not uniformly `performance`** "
                 f"({gov['counts']}). Timings may drift as the campaign runs.")
        L.append("")

    L.append("## Measured fixed overhead")
    L.append("")
    for k, v in sorted(overhead.items()):
        if isinstance(v, float):
            L.append(f"- `{k}`: {v * 1000:.1f} ms")
    L.append("")
    L.append("`lsdtt_startup_seconds` is process start, dynamic linking and "
             "GDAL driver registration. An lsdtt3 row whose time is "
             "whole-process wall clock and below "
             f"{MIN_TRUSTWORTHY_FACTOR:g}x that figure is flagged "
             "`overhead_dominated` and excluded from the ratio table. Rows "
             "timed by the CLI's own reported compute time already exclude "
             "startup and are never flagged.")
    L.append("")
    wb = overhead.get("whitebox_dispatch_seconds")
    if wb is not None and wb > 0.1:
        L.append(f"> **`whitebox_dispatch_seconds` is {wb * 1000:.0f} ms.** "
                 "That is the signature of `whitebox_reference.py` resolving "
                 "tool handles inside its timing loop, which inflates every "
                 "Whitebox measurement by a fixed ~190 ms. **Do not quote the "
                 "Whitebox rows in this campaign** until that is fixed; see "
                 "`_make_call` in that file.")
        L.append("")

    ok = [r for r in rows if r["status"] == "ok"]

    # Only rows measured on the SAME basis may share a ratio table. An lsdtt3
    # row timed as whole-process wall clock includes startup, read and write
    # that the in-process reference rows exclude; putting the two in one table
    # produces a ratio that is mostly process startup. Being comfortably above
    # the startup floor does not fix that -- it only makes startup a smaller
    # share of a number that is still measuring something else.
    COMPARABLE_SCOPES = {"cli_reported_compute", "in_process_routine"}

    L.append("## Throughput (median, Mcells/s, higher is better)")
    L.append("")
    L.append("Like-for-like only: every row here is an in-process timing of "
             "the routine itself.")
    L.append("")
    quotable = [r for r in ok
                if not r["report_only"] and not r["overhead_dominated"]
                and r["timing_scope"] in COMPARABLE_SCOPES]
    if not quotable:
        L.append("_No row was measured on a comparable basis. This is expected "
                 "on the `laptop` profile with small fixtures; run the "
                 "`server` profile for quotable numbers._")
    else:
        _ratio_table(L, quotable)

    whole = [r for r in ok if r["timing_scope"] == "whole_process"]
    if whole:
        L.append("")
        L.append(f"### Whole-process rows ({len(whole)}), not comparable above")
        L.append("")
        L.append("`lsdtt-raster-preprocessing` and `lsdtt-drainage-analysis` "
                 "do not print their own compute times, so the only figure "
                 "available for them is whole-process wall clock: process "
                 "start, GDAL read, the routine, the write and teardown. "
                 "Against a reference timed in-process that is not a "
                 "comparison, so these rows are reported alone. Subtract "
                 "roughly `lsdtt_startup_seconds` above to estimate the "
                 "routine, and treat the result as an estimate. The clean fix "
                 "is for those CLIs to print `Compute time:` the way "
                 "`lsdtt-surface-metrics` already does.")
        L.append("")
        _plain_table(L, whole)

    excluded = [r for r in ok if r["overhead_dominated"] and not r["report_only"]
                and r["timing_scope"] != "whole_process"]
    if excluded:
        L.append("")
        L.append(f"### Excluded: overhead-dominated ({len(excluded)} rows)")
        L.append("")
        L.append("Measured successfully, but too close to fixed overhead to "
                 "quote as an algorithm comparison.")
        L.append("")
        _plain_table(L, excluded)

    report_only = [r for r in ok if r["report_only"]]
    if report_only:
        L.append("")
        L.append(f"### Report-only ({len(report_only)} rows)")
        L.append("")
        L.append("The engines are not performing the same work here, so a "
                 "ratio would be a claim about the definitions rather than the "
                 "implementations. Absolute times only.")
        L.append("")
        _plain_table(L, report_only)

    L.append("")
    L.append("## Parallel scaling (lsdtt3 and whitebox)")
    L.append("")
    _scaling_table(L, ok)

    L.append("")
    L.append("## Peak memory (MiB, whole child process)")
    L.append("")
    _memory_table(L, ok)

    bad = [r for r in rows if r["status"] in ("error", "skipped", "unavailable")]
    if bad:
        L.append("")
        L.append(f"## Not measured ({len(bad)} rows)")
        L.append("")
        L.append("| suite | case | engine | fixture | status | reason |")
        L.append("|---|---|---|---|---|---|")
        seen = set()
        for r in bad:
            key = (r["case"], r["engine"], r["status"], r["note"][:60])
            if key in seen:
                continue
            seen.add(key)
            L.append(f"| {r['suite']} | {r['case']} | {r['engine']} | "
                     f"{r['fixture']} | {r['status']} | {r['note'][:80]} |")

    L.append("")
    L.append("---")
    L.append("")
    L.append("### Reading these numbers")
    L.append("")
    L.append("- The `timing_scope` column says how each row was measured, and "
             "rows are only ever compared with others of the same scope. "
             "`cli_reported_compute` is the compute time the lsdtt3 CLI printed "
             "itself, excluding startup and I/O. `in_process_routine` is a "
             "reference script timing the routine in-process. `whole_process` "
             "is wall clock for an entire lsdtt3 CLI run, used only where that "
             "CLI prints no timing of its own; those rows are tabled "
             "separately and never enter a ratio.")
    L.append("- I/O is excluded from `seconds` for every engine and reported "
             "separately in `seconds_read` / `seconds_write` in the CSV.")
    L.append("- RichDEM and TopoToolbox are serial for these routines; their "
             "rows exist only at 1 thread and are marked `serial_reference`. A "
             "flat scaling curve for them is not a finding.")
    L.append("- Values are not asserted anywhere in this campaign. Correctness "
             "lives in the CTest suite; see `docs/testing_and_benchmarking.md`.")
    path.write_text("\n".join(L) + "\n")


def _fixture_order(rows: list[dict]) -> list[str]:
    seen = {}
    for r in rows:
        try:
            seen[r["fixture"]] = int(r["n_cells"]) if r["n_cells"] else 0
        except ValueError:
            seen[r["fixture"]] = 0
    return [k for k, _ in sorted(seen.items(), key=lambda kv: kv[1])]


def _ratio_table(L: list[str], rows: list[dict]) -> None:
    """Throughput per engine at the highest common thread count."""
    engines = sorted({r["engine"] for r in rows})
    fixtures = _fixture_order(rows)
    L.append("| case | fixture | threads | " + " | ".join(engines) + " |")
    L.append("|---|---|---|" + "---|" * len(engines))
    for case_name in sorted({r["case"] for r in rows}):
        for fx in fixtures:
            for threads in sorted({int(r["threads"]) for r in rows
                                   if r["case"] == case_name
                                   and r["fixture"] == fx}):
                cells = []
                any_val = False
                for eng in engines:
                    match = [r for r in rows if r["case"] == case_name
                             and r["fixture"] == fx and r["engine"] == eng
                             and int(r["threads"]) == threads]
                    if match and match[0]["mcells_per_second"]:
                        cells.append(f"{float(match[0]['mcells_per_second']):.2f}")
                        any_val = True
                    else:
                        cells.append("-")
                if any_val:
                    L.append(f"| {case_name} | {fx} | {threads} | "
                             + " | ".join(cells) + " |")


def _plain_table(L: list[str], rows: list[dict]) -> None:
    L.append("| case | engine | fixture | threads | median (ms) | +/- (ms) |")
    L.append("|---|---|---|---|---|---|")
    for r in sorted(rows, key=lambda r: (r["case"], r["engine"], r["fixture"],
                                         int(r["threads"]))):
        med = float(r["seconds_median"]) * 1000
        sd = float(r["seconds_stdev"] or 0) * 1000
        L.append(f"| {r['case']} | {r['engine']} | {r['fixture']} | "
                 f"{r['threads']} | {med:.2f} | {sd:.2f} |")


def _scaling_table(L: list[str], rows: list[dict]) -> None:
    """Speedup relative to each engine's own 1-thread time.

    Relative to ITS OWN baseline, not to lsdtt3's: a scaling curve answers "how
    well does this implementation use more cores", which is a separate question
    from "which is faster", and mixing them makes both unreadable.
    """
    par = [r for r in rows if r["engine"] in ("lsdtt3", "whitebox")]
    if not par:
        L.append("_No parallel-capable rows measured._")
        return
    threads = sorted({int(r["threads"]) for r in par})
    L.append("| case | engine | fixture | "
             + " | ".join(f"{t}t" for t in threads) + " |")
    L.append("|---|---|---|" + "---|" * len(threads))
    for case_name in sorted({r["case"] for r in par}):
        for eng in sorted({r["engine"] for r in par}):
            for fx in _fixture_order(par):
                base = [r for r in par if r["case"] == case_name
                        and r["engine"] == eng and r["fixture"] == fx
                        and int(r["threads"]) == 1]
                if not base or not base[0]["seconds_median"]:
                    continue
                t1 = float(base[0]["seconds_median"])
                cells = []
                for t in threads:
                    m = [r for r in par if r["case"] == case_name
                         and r["engine"] == eng and r["fixture"] == fx
                         and int(r["threads"]) == t]
                    if m and m[0]["seconds_median"]:
                        cells.append(f"{t1 / float(m[0]['seconds_median']):.2f}x")
                    else:
                        cells.append("-")
                if any(c != "-" for c in cells[1:]):
                    L.append(f"| {case_name} | {eng} | {fx} | "
                             + " | ".join(cells) + " |")


def _memory_table(L: list[str], rows: list[dict]) -> None:
    mem = [r for r in rows if r["peak_rss_mb"]]
    if not mem:
        L.append("_Peak RSS unavailable on this platform (needs wait4)._")
        return
    engines = sorted({r["engine"] for r in mem})
    L.append("| case | fixture | " + " | ".join(engines) + " |")
    L.append("|---|---|" + "---|" * len(engines))
    for case_name in sorted({r["case"] for r in mem}):
        for fx in _fixture_order(mem):
            cells, any_val = [], False
            for eng in engines:
                # Peak across thread counts: memory is reported per engine and
                # fixture, and the worst case is the number that matters for
                # "will this fit".
                m = [float(r["peak_rss_mb"]) for r in mem
                     if r["case"] == case_name and r["fixture"] == fx
                     and r["engine"] == eng]
                if m:
                    cells.append(f"{max(m):.0f}")
                    any_val = True
                else:
                    cells.append("-")
            if any_val:
                L.append(f"| {case_name} | {fx} | " + " | ".join(cells) + " |")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profile", choices=sorted(PROFILES), default="laptop")
    ap.add_argument("--build-dir", default="build",
                    help="directory holding the built lsdtt-* executables")
    ap.add_argument("--out-dir", default="benchmark_campaigns",
                    help="campaign root; one subdirectory per run")
    ap.add_argument("--python-exe", default=sys.executable,
                    help="interpreter for the reference scripts; must be the "
                         "one with richdem/topotoolbox/whitebox installed")
    ap.add_argument("--suites", default="",
                    help="comma-separated subset of suites (default: profile)")
    ap.add_argument("--cases", default="",
                    help="comma-separated subset of case names")
    ap.add_argument("--fixtures", default="",
                    help="comma-separated fixture keys, overriding the profile")
    ap.add_argument("--threads", default="",
                    help="comma-separated thread counts, overriding the profile")
    ap.add_argument("--repeats", type=int, default=0,
                    help="timed repeats, overriding the profile")
    ap.add_argument("--no-fetch", action="store_true",
                    help="never download a missing fixture; skip it instead")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the measurement plan and exit")
    ap.add_argument("--list", action="store_true",
                    help="print the case registry as JSON and exit")
    args = ap.parse_args()

    if args.list:
        print(json.dumps({"profiles": PROFILES, "cases": CASES}, indent=2))
        return 0

    profile = dict(PROFILES[args.profile])
    if args.fixtures:
        profile["fixtures"] = [s.strip() for s in args.fixtures.split(",")
                               if s.strip()]
    if args.threads:
        profile["threads"] = [int(s) for s in args.threads.split(",")
                              if s.strip()]
    if args.repeats:
        profile["repeats"] = args.repeats
    if args.suites:
        profile["suites"] = [s.strip() for s in args.suites.split(",")
                             if s.strip()]
    if args.no_fetch:
        profile["allow_fetch"] = False

    selected_cases = {
        name: case for name, case in CASES.items()
        if case["suite"] in profile["suites"]
        and (not args.cases or name in
             {s.strip() for s in args.cases.split(",")})
    }
    if not selected_cases:
        print("ERROR: no cases selected", file=sys.stderr)
        return 2

    build_dir = pathlib.Path(args.build_dir)
    if not build_dir.is_absolute():
        build_dir = (REPO_ROOT / build_dir).resolve()

    # Warn early rather than producing a campaign of `unavailable` rows.
    missing = sorted({c["lsdtt"]["exe"] for c in selected_cases.values()
                      if not (build_dir / c["lsdtt"]["exe"]).exists()})
    if missing:
        print(f"WARNING: not found in {build_dir}: {', '.join(missing)}\n"
              f"         lsdtt3 rows for those cases will be 'unavailable'.\n"
              f"         Build them, or pass --build-dir.", file=sys.stderr)

    if args.dry_run:
        n = 0
        print(f"profile: {args.profile} - {profile['description']}")
        print(f"fixtures: {profile['fixtures']}")
        print(f"threads:  {profile['threads']}   repeats: {profile['repeats']}")
        print(f"build:    {build_dir}")
        print("")
        for name, case in sorted(selected_cases.items()):
            engines = ["lsdtt3"] + sorted(case["engines"])
            per_engine = []
            for e in engines:
                t = 1 if e in SERIAL_REFERENCES else len(profile["threads"])
                per_engine.append(f"{e}x{t}")
                n += t * len(profile["fixtures"])
            print(f"  {case['suite']:11s} {name:22s} {' '.join(per_engine)}")
        print(f"\n{n} measurements, each "
              f"{profile['repeats']} repeats + 1 warm-up")
        return 0

    stamp = now_stamp()
    run_dir = (pathlib.Path(args.out_dir) /
               f"{args.profile}_{socket.gethostname()}_{stamp}")
    if not run_dir.is_absolute():
        run_dir = REPO_ROOT / run_dir
    work_dir = run_dir / "work"
    run_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    print(f"campaign: {run_dir}")
    print(f"profile:  {args.profile} ({profile['description']})")

    hardware = capture_hardware_manifest(build_dir, args.python_exe)
    (run_dir / "hardware_manifest.json").write_text(
        json.dumps(hardware, indent=2))

    # Resolve fixtures up front so a download does not interrupt the timings.
    fetch_exe = build_dir / "lsdtt-fetch-raster"
    fixtures: dict[str, pathlib.Path] = {}
    for key in profile["fixtures"]:
        path = ensure_fixture(key, fetch_exe, profile["allow_fetch"])
        if path is None:
            print(f"  fixture {key}: unavailable, skipping")
            continue
        fixtures[key] = path
        print(f"  fixture {key}: {path}")
    if not fixtures:
        print("ERROR: no fixtures available. Build lsdtt-fetch-raster to let "
              "the suite download them, or place them in "
              f"{FIXTURE_DIR}.", file=sys.stderr)
        return 2

    fixture_meta = {k: {"sha256_12": sha_prefix(p), "n_cells": n_cells_of(p),
                        "path": str(p)}
                    for k, p in fixtures.items()}

    print("calibrating fixed overhead...")
    smallest = min(fixtures.values(),
                   key=lambda p: n_cells_of(p) or (1 << 62))
    overhead = calibrate_overheads(build_dir, smallest, args.python_exe)
    for k, v in overhead.items():
        if isinstance(v, float):
            print(f"  {k}: {v * 1000:.1f} ms")

    manifest = {
        "created_at": dt.datetime.now().isoformat(timespec="seconds"),
        "profile": args.profile,
        "profile_settings": profile,
        "build_dir": str(build_dir),
        "fixtures": fixture_meta,
        "overhead": overhead,
        "hardware": hardware,
        "min_trustworthy_factor": MIN_TRUSTWORTHY_FACTOR,
    }

    rows: list[dict[str, Any]] = []
    # Matches --dry-run exactly: lsdtt3 runs at every thread count, a serial
    # reference runs once, a parallel reference runs at every thread count.
    n_threads = len(profile["threads"])
    total = sum(
        len(fixtures) * (n_threads + sum(1 if e in SERIAL_REFERENCES
                                         else n_threads
                                         for e in c["engines"]))
        for c in selected_cases.values())
    done = 0

    for case_name, case in sorted(selected_cases.items(),
                                  key=lambda kv: (kv[1]["suite"], kv[0])):
        for fx_key, fx_path in fixtures.items():
            meta = fixture_meta[fx_key]

            # Accumulation cases run on a pre-filled DEM, shared by every
            # engine. See prepare_filled_fixture().
            in_path = fx_path
            if case.get("needs_filled"):
                filled = prepare_filled_fixture(fx_path, work_dir, build_dir)
                if filled is None:
                    rows.append(build_row(
                        case["suite"], case_name, "lsdtt3", fx_key,
                        meta["sha256_12"], meta["n_cells"],
                        profile["threads"][0],
                        {"error": "could not prepare filled DEM "
                                  "(lsdtt-raster-preprocessing missing "
                                  "or failed)"},
                        case, overhead))
                    continue
                in_path = filled

            for threads in profile["threads"]:
                # lsdtt3
                info = measure_lsdtt(case_name, case, in_path, work_dir,
                                     threads, profile["repeats"], build_dir)
                rows.append(build_row(case["suite"], case_name, "lsdtt3",
                                      fx_key, meta["sha256_12"],
                                      meta["n_cells"], threads, info, case,
                                      overhead))
                done += 1
                _progress(done, total, case_name, "lsdtt3", fx_key, threads,
                          rows[-1])

                for engine, routine in sorted(case["engines"].items()):
                    # A serial reference is measured once, at 1 thread. Running
                    # it again at 2/4/8 would record the same number under a
                    # different label and imply it had been tested there.
                    if engine in SERIAL_REFERENCES and threads != 1:
                        continue
                    info = measure_reference(engine, routine, in_path,
                                             threads, profile["repeats"],
                                             args.python_exe)
                    rows.append(build_row(case["suite"], case_name, engine,
                                          fx_key, meta["sha256_12"],
                                          meta["n_cells"], threads, info, case,
                                          overhead))
                    done += 1
                    _progress(done, total, case_name, engine, fx_key, threads,
                              rows[-1])

    csv_path = run_dir / "benchmark_raw.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    (run_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
    write_summary(run_dir / "SUMMARY.md", rows, manifest, overhead)

    # The per-case output rasters are large and are not evidence: nothing in
    # this campaign compares values. Keeping them would make a campaign
    # directory too big to send.
    shutil.rmtree(work_dir, ignore_errors=True)

    n_ok = sum(1 for r in rows if r["status"] == "ok")
    n_err = sum(1 for r in rows if r["status"] == "error")
    print(f"\n{n_ok}/{len(rows)} measurements ok, {n_err} errors")
    print(f"  {csv_path}")
    print(f"  {run_dir / 'SUMMARY.md'}")

    # Exit 0 even with errors: an unavailable engine is information, not a
    # failure, and a non-zero exit would make this look like a test.
    return 0


def _progress(done, total, case_name, engine, fixture, threads, row) -> None:
    status = row["status"]
    detail = ""
    if status == "ok":
        detail = f"{float(row['seconds_median']) * 1000:8.2f} ms"
        if row["overhead_dominated"]:
            detail += "  (overhead-dominated)"
    elif row["note"]:
        detail = row["note"][:60]
    print(f"  [{done:4d}/{total:4d}] {case_name:22s} {engine:12s} "
          f"{fixture:16s} {threads}t  {status:12s} {detail}")


if __name__ == "__main__":
    sys.exit(main())
