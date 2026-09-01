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

  A fixed per-call cost is measured for EVERY engine, on an 8x8 crop (64 cells)
  where the algorithm itself is microseconds.  Do not assume an engine has none
  because it is "timed in-process" -- that was assumed once here and it was
  wrong.  Measured on crops of alpine_small, median of 7, fill:

      grid     richdem   topotoolbox   whitebox(1t)
      8^2       9.83 ms      0.036 ms       0.038 ms
      256^2    13.73 ms      1.966 ms       1.440 ms

  TopoToolbox and Whitebox have no meaningful fixed cost.  **RichDEM has
  ~10 ms**, which on a 256^2 grid is most of its measured time: its own
  internal timer reports the 8x8 fill as 1.2e-05 s, so ~99.8% of what a caller
  measures there is the Python binding -- a per-call rdarray construction plus
  the citation block it formats and prints on every invocation.

  lsdtt3 CLI process startup (fork, dynamic linking, GDAL driver registration)
  is measured separately as `lsdtt3_startup_seconds`.  It applies only to rows
  with no CLI-reported compute time; see TIMING SCOPE below.

  Whitebox tool-handle resolution is a THIRD, historical trap.
  `wbe.hydrology.sink` is a registry resolution through the Rust layer costing
  ~190 ms PER ACCESS.  whitebox_reference.py binds the handle once, outside the
  timing loop, which took fill on 256^2 from 0.207 s to 0.019 s -- a 10.9x
  artefact that had flattered lsdtt3.  `whitebox_dispatch_seconds` is the
  tripwire that catches it coming back.

  Any row whose median is below MIN_TRUSTWORTHY_FACTOR times ITS OWN engine's
  floor is marked `overhead_dominated` and excluded from the summary's ratio
  table.


TIMING SCOPE -- CHECK THIS BEFORE QUOTING A RATIO

  Rows are only ever compared with others measured the same way:

    cli_reported_compute  the compute time an lsdtt3 CLI printed itself,
                          excluding process startup and I/O
    in_process_routine    a reference script timing the routine in-process
    whole_process         wall clock for an entire lsdtt3 CLI run

  lsdtt-surface-metrics, lsdtt-raster-preprocessing and lsdtt-drainage-analysis
  all print "Compute time:" lines, which this driver parses and sums, so the
  flow, surface and depression suites are like-for-like.  A CLI that prints
  none falls back to whole_process, and those rows are tabled separately and
  never enter a ratio -- on a small raster they are mostly process startup
  (lsdtt3's fill measured 7.6 ms of compute inside a 279 ms process).
"""

from __future__ import annotations

import argparse
import contextlib
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

# WHERE THINGS ARE WRITTEN
#
# Three separate locations, because they have three different lifetimes and
# three different size and speed requirements:
#
#   FIXTURE_DIR   the downloaded region sources and their crops. Large
#                 (~512 MB for the server profile), reused across every
#                 campaign, and worth keeping. Default is inside the repo,
#                 which is WRONG on a cluster where the repo sits on a quota'd
#                 home filesystem -- override with --fixture-dir.
#
#   out-dir       the campaign directories. Small, one per run, the thing you
#                 actually want to keep and send. A RELATIVE --out-dir resolves
#                 against the CURRENT DIRECTORY, so running from a data
#                 directory puts results in that data directory.
#
#   scratch       staged fixtures and intermediate rasters (see
#                 choose_scratch_root). Node-local, disposable, and must be
#                 fast -- measurements are taken here.
#
# FIXTURE_DIR is module-level because the fixture helpers all reach for it;
# main() rebinds it once from --fixture-dir before anything uses it.
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
    "fill_spill_merge": {
        "suite": "depression",
        "note": ("lsdtt3 builds a nested depression hierarchy (Barnes 2020) "
                 "then Fill-Spill-Merge (Barnes 2021); the references only "
                 "label or fill. TIMING ONLY, and NOT a fair ratio"),
        "engines": {"whitebox": "depression_sink"},
        # The depression hierarchy is not separately switchable: it is built as
        # the first half of the barnesetal2021 fill, which is the only CLI path
        # that reaches it. An earlier version of this table asked for
        # `write_depression_hierarchy_raster`, which does not exist -- the CLI
        # ignored it, did no work, and the case happily reported a timing for
        # a run that had computed nothing.
        "lsdtt": {"exe": "lsdtt-raster-preprocessing",
                  "params": {"write_fill_raster": "true",
                             "fill_algorithm_choice": "barnesetal2021",
                             "fill_min_slope": "0"},
                  "output": "_fill.tif"},
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
        "description": ("reduced ladder and threads; a real baseline for a "
                        "laptop and the smoke test before a server campaign"),
        "regions": ["alpine", "lakeland"],
        # Small rungs, so the one-off source download is a few MB rather than
        # a few hundred: the source is fetched at the LARGEST rung.
        "ladder": [500, 1000],
        "threads": [1, 2, 4],
        "repeats": 3,
        "suites": ["flow", "surface", "mfd", "depression"],
        "allow_fetch": True,
    },
    "server": {
        "description": "full ladder and thread sweep; the quotable profile",
        "regions": ["alpine", "lakeland"],
        # 8000 px at 30 m is a 240 km box, ~256 MB per region, fetched ONCE and
        # cropped into every rung. 64M cells is enough work that a 32-thread
        # run has something to do; the current 2400 px fixture was not.
        "ladder": [1000, 2000, 4000, 8000],
        "threads": [1, 2, 4, 8],
        "repeats": 5,
        "suites": ["flow", "surface", "mfd", "depression"],
        "allow_fetch": True,
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
                    cwd: pathlib.Path | None = None,
                    timeout: float | None = None) -> dict[str, Any]:
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

    `timeout` kills the child after that many seconds and returns
    timed_out=True.  A server campaign runs unattended for hours and the
    engines are not uniformly fast -- Whitebox's breach on low-relief terrain
    measured ~300x slower than lsdtt3's -- so without a per-measurement
    deadline one pathological cell can stall a whole campaign.  A killed
    measurement becomes a `timeout` row and the campaign continues.

    Falls back to subprocess.run with peak_rss_mb=None where wait4 is
    unavailable (Windows).  A missing memory column is better than a wrong one.
    """
    import tempfile
    import threading

    if not hasattr(os, "wait4"):
        t0 = time.perf_counter()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  env=env, cwd=str(cwd) if cwd else None,
                                  timeout=timeout)
        except subprocess.TimeoutExpired:
            return {"returncode": None, "stdout": "", "stderr": "",
                    "wall_seconds": time.perf_counter() - t0,
                    "peak_rss_mb": None, "timed_out": True}
        return {"returncode": proc.returncode, "stdout": proc.stdout,
                "stderr": proc.stderr,
                "wall_seconds": time.perf_counter() - t0,
                "peak_rss_mb": None, "timed_out": False}

    timed_out = False

    with tempfile.TemporaryFile(mode="w+") as errf:
        t0 = time.perf_counter()
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errf,
                                text=True, env=env,
                                cwd=str(cwd) if cwd else None)

        # A watchdog thread rather than Popen.communicate(timeout=...):
        # communicate() reaps the child, which consumes the rusage before
        # wait4 can read it, and the memory column is the whole reason this
        # function exists.
        timer = None
        if timeout is not None:
            def _kill():
                nonlocal timed_out
                timed_out = True
                try:
                    proc.kill()
                except OSError:
                    pass
            timer = threading.Timer(timeout, _kill)
            timer.start()

        try:
            stdout = proc.stdout.read() if proc.stdout else ""
            if proc.stdout:
                proc.stdout.close()
            pid, status, ru = os.wait4(proc.pid, 0)
        finally:
            if timer is not None:
                timer.cancel()

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
        "timed_out": timed_out,
        # ru_maxrss is kibibytes on Linux, bytes on macOS. Linux is where the
        # server campaign runs; on macOS this reads ~1024x high, which is
        # obvious rather than subtle.
        "peak_rss_mb": ru.ru_maxrss / 1024.0,
    }


# -----------------------------------------------------------------------------
# Fixtures: fetch one large raster per region, then crop the size ladder
#
# WHY CROP RATHER THAN FETCH EACH SIZE
#
# A size ladder built from independent downloads confounds size with terrain.
# At 30 m, 2400 px is a 72 km box and 8000 px is 240 km: scale the box up and
# an "alpine" fixture stops being alpine, so a scaling curve measured across
# them is partly measuring the landscape changing. Cropping every rung from ONE
# source keeps the terrain as constant as it can be and needs one download per
# region instead of one per rung.
#
# WHY NOT SYNTHETIC
#
# lsdtt-raster-synthetic-tests can generate a diamond-square surface of any
# size in ~1.3 s with no network, which is tempting. It is not usable here.
# Measured depression density (fraction of cells a no-epsilon fill raises,
# RichDEM, 256^2):
#
#     REAL alpine  (Italian Alps)       0.94 %
#     REAL lakeland (Finland)           1.29 %
#     SYNTH diamond-square, no noise   11.54 %
#     SYNTH roughness 0.70 + noise     34.24 %
#
# 10-30x more depressions than real terrain, at 1 m cell size with a few metres
# of relief. Since depression structure is the single largest driver of the
# fill and breach timings -- Whitebox measured 265x slower on low-relief real
# terrain than on high-relief -- a fractal surface would inflate every
# depression-handling measurement for every engine. Synthetic surfaces are for
# the analytic oracles in the CTest suite, not for this.
#
# THE CONFOUND THAT REMAINS, AND HOW IT IS MADE VISIBLE
#
# Centred crops of increasing size progressively reach further from the centre,
# so terrain character still drifts up the ladder -- a 240 km box centred in the
# Himalaya necessarily includes both the Gangetic plain and the Tibetan
# plateau. This cannot be engineered away with real terrain. It is instead
# MEASURED: `describe_fixture()` records relief and depression density for
# every rung into the run manifest, so a reader can see whether a scaling curve
# is confounded rather than having to assume it is not.
# -----------------------------------------------------------------------------

# Region centres and their UTM zones. Both are deliberate end-members.
REGIONS = {
    "alpine": {
        "centre": (86.9, 27.9),
        "epsg": 32645,                 # UTM 45N
        "description": "High Himalaya, Khumbu/Makalu region, Nepal",
    },
    "lakeland": {
        "centre": (-92.0864236880431, 50.57343374931578),
        "epsg": 32615,                 # UTM 15N
        "description": ("Canadian Shield lakeland, northwestern Ontario -- "
                        "lake-covered, and by far the harder depression test"),
    },
}

GRID_SPACING_M = 30.0

# Measured on a 256^2 probe of each region, RichDEM no-epsilon fill:
#
#     alpine   (Himalaya)  relief 1662 m   6.40 % depression cells
#     lakeland (Ontario)   relief   94 m  18.51 % depression cells
#
# For comparison the previous fixtures (Italian Alps / Finland) were 0.94 % and
# 1.29 %. These regions are 7x and 14x more depression-rich, which is the point
# -- but expect the depression-handling cases to be far slower here, and expect
# Whitebox's breach on lakeland to hit --timeout-seconds at the larger rungs.


def source_raster_path(region: str, source_px: int) -> pathlib.Path:
    return FIXTURE_DIR / f"{region}_src{source_px}.tif"


def ensure_source_raster(region: str, source_px: int,
                         fetch_exe: pathlib.Path | None,
                         allow_fetch: bool) -> pathlib.Path | None:
    """Fetch the one large source raster for a region, if not already cached.

    The bounding box is widened by 25% so the exact centre crop cannot run off
    the edge onto nodata after reprojection -- the same margin
    run_flow_benchmark.py uses, for the same reason.

    Fetching goes through the lsdtt3 CLI (`lsdtt-fetch-raster`), not a Python
    helper, so the suite depends on no extra download stack.
    """
    import math

    out = source_raster_path(region, source_px)
    if out.exists():
        return out
    if not allow_fetch or fetch_exe is None or not fetch_exe.exists():
        return None
    if region not in REGIONS:
        return None

    lon, lat = REGIONS[region]["centre"]
    epsg = REGIONS[region]["epsg"]
    half_m = 0.5 * source_px * GRID_SPACING_M * 1.25
    dlat = half_m / 111_320.0
    dlon = dlat / max(math.cos(math.radians(lat)), 1e-6)

    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    param_dir = FIXTURE_DIR / "_fetch_params"
    param_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{region}_src{source_px}"
    param_path = param_dir / f"{stem}.param"
    param_path.write_text(
        f"write_path: {FIXTURE_DIR}\n"
        f"write_fname: {stem}_raw\n"
        f"west: {lon - dlon:.10f}\n"
        f"south: {lat - dlat:.10f}\n"
        f"east: {lon + dlon:.10f}\n"
        f"north: {lat + dlat:.10f}\n"
        f"dem_source: cop30_aws\n"
        f"target_epsg: {epsg}\n"
        f"grid_spacing: {GRID_SPACING_M}\n"
        f"size_pixels: {source_px}\n"
    )

    km = source_px * GRID_SPACING_M / 1000.0
    print(f"  fetching {region} source {source_px}px (~{km:.0f} km box, "
          f"{REGIONS[region]['description']})")
    print(f"    this is a one-off download of roughly "
          f"{source_px * source_px * 4 / 1e6:.0f} MB; every ladder rung is "
          f"cropped from it")

    proc = subprocess.run([str(fetch_exe), str(param_dir), param_path.name],
                          capture_output=True, text=True, env=system_env(1))
    if proc.returncode != 0:
        print(f"    fetch failed:\n{proc.stdout[-1200:]}\n{proc.stderr[-1200:]}",
              file=sys.stderr)
        return None

    written = FIXTURE_DIR / f"{stem}_raw_DEM.tif"
    if not written.exists():
        print(f"    fetch reported success but {written} is missing",
              file=sys.stderr)
        return None
    written.rename(out)
    return out


def crop_ladder_rung(source: pathlib.Path, region: str,
                     size_px: int) -> pathlib.Path | None:
    """A centred square crop of `size_px` from the region's source raster.

    Centred so that every rung shares the same midpoint and the ladder grows
    outward symmetrically. Written uncompressed GeoTIFF so that read time is
    the filesystem's, not the decompressor's -- the engines are being compared
    on their algorithms, and a compressed fixture would tax whichever of them
    reads through a slower GDAL path.
    """
    out = FIXTURE_DIR / f"{region}_{size_px}.tif"
    if out.exists():
        return out
    try:
        from osgeo import gdal
        gdal.UseExceptions()
        ds = gdal.Open(str(source))
        w, h = ds.RasterXSize, ds.RasterYSize
        if size_px > w or size_px > h:
            return None
        x0 = (w - size_px) // 2
        y0 = (h - size_px) // 2
        # `ds`, not `str(ds)`: stringifying a Dataset yields its repr, which
        # GDAL then tries to open as a filename and fails with a bare
        # "No such file or directory" naming the repr.
        gdal.Translate(str(out), ds, srcWin=[x0, y0, size_px, size_px],
                       creationOptions=["TILED=NO", "COMPRESS=NONE"])
    except Exception as exc:
        print(f"    crop {region}_{size_px} failed: {exc}", file=sys.stderr)
        return None
    return out if out.exists() else None


@contextlib.contextmanager
def suppress_fds():
    """Silence a C/Rust library writing straight to the process's stdout.

    RichDEM prints multi-kilobyte ANSI progress bars and a citation block from
    the C++ layer, which writes past any Python-level capture. Redirect the
    file descriptors themselves.
    """
    saved = os.dup(1), os.dup(2)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        os.close(devnull)
        os.close(saved[0])
        os.close(saved[1])


def describe_fixture(path: pathlib.Path) -> dict[str, Any]:
    """Relief and TRUE depression density for one fixture.

    Recorded per rung so that terrain drift up the ladder is visible in the
    manifest. A scaling curve across rungs whose depression density doubles is
    not a clean scaling curve, and the reader should be able to see that
    without re-deriving it. Depression density is the right statistic to track
    because it is the single largest driver of the fill and breach timings.

    Uses RichDEM's priority-flood fill and counts cells the fill actually
    raised. An earlier version used a 5x5 `scipy.ndimage.grey_closing` as a
    cheap proxy and it was badly wrong -- it flags any local concavity rather
    than a closed basin, and reported 29-55% where the true figures are 6.4%
    (Himalaya) and 18.5% (Ontario). A statistic that is off by 3x is worse than
    no statistic, because it looks authoritative in the manifest.

    RichDEM is already a hard dependency of the benchmark environment. The cost
    is one fill per rung per campaign -- tens of seconds at the 8000px rung,
    against a campaign measured in hours. Returns what it can and never raises;
    a missing description is a cosmetic loss, not a failed campaign.
    """
    try:
        from osgeo import gdal
        import numpy as np
        gdal.UseExceptions()
        ds = gdal.Open(str(path))
        band = ds.GetRasterBand(1)
        arr = band.ReadAsArray().astype("float64")
        nodata = band.GetNoDataValue()
        mask = np.isfinite(arr)
        if nodata is not None:
            mask &= (arr != nodata)
        if not mask.any():
            return {"error": "no valid cells"}

        valid = arr[mask]
        out: dict[str, Any] = {
            "relief_m": float(valid.max() - valid.min()),
            "mean_elev_m": float(valid.mean()),
            "nodata_fraction": float(1.0 - mask.mean()),
        }

        try:
            import richdem as rd
            with suppress_fds():
                dem = rd.LoadGDAL(str(path))
                filled = rd.FillDepressions(dem, epsilon=False, in_place=False)
            a = np.asarray(dem, dtype="float64")
            b = np.asarray(filled, dtype="float64")
            raised = (b > a + 1e-6) & mask
            out["depression_fraction"] = float(raised.sum() / mask.sum())
            out["depression_volume_m3_per_cell"] = float(
                (b - a)[raised].sum() / mask.sum()) if raised.any() else 0.0
            out["depression_method"] = "richdem priority-flood, epsilon=False"
        except Exception as exc:
            out["depression_error"] = str(exc)[:200]

        return out
    except Exception as exc:
        return {"error": str(exc)[:200]}


def resolve_fixtures(profile: dict, fetch_exe: pathlib.Path | None,
                     explicit: list[str] | None) -> dict[str, pathlib.Path]:
    """Build the fixture set for a profile: fetch each source, crop each rung.

    The source is fetched at the profile's LARGEST rung, not at a fixed size,
    so the laptop profile needs a few MB rather than the server profile's few
    hundred. A rung cropped from a 1000px source is not byte-identical to the
    same rung cropped from an 8000px source -- different reprojection extent --
    which is why every fixture's SHA is recorded and why the manifest names the
    source each rung came from.
    """
    fixtures: dict[str, pathlib.Path] = {}

    if explicit:
        # An explicit --fixtures list bypasses the ladder entirely and is taken
        # as literal names in the cache. This is the escape hatch for the
        # pre-existing fixtures and for anything hand-placed.
        for key in explicit:
            path = FIXTURE_DIR / f"{key}.tif"
            if path.exists():
                fixtures[key] = path
            else:
                print(f"  fixture {key}: not in {FIXTURE_DIR}, skipping")
        return fixtures

    rungs = profile["ladder"]
    source_px = max(rungs)
    for region in profile["regions"]:
        source = ensure_source_raster(region, source_px, fetch_exe,
                                      profile["allow_fetch"])
        if source is None:
            print(f"  region {region}: source raster unavailable, skipping "
                  f"(build lsdtt-fetch-raster, or pass --no-fetch off)")
            continue
        for size_px in rungs:
            rung = crop_ladder_rung(source, region, size_px)
            if rung is None:
                print(f"  {region}_{size_px}: crop failed, skipping")
                continue
            fixtures[f"{region}_{size_px}"] = rung
            print(f"  fixture {region}_{size_px}: {rung}")
    return fixtures


# -----------------------------------------------------------------------------
# Measurement
# -----------------------------------------------------------------------------
# Filesystem types that are far slower than local disk for the many small
# reads and writes some engines perform inside a timed call. Measured on this
# repository's own WSL2 checkout (/mnt/c, 9p): Whitebox's fill on a 256^2
# raster took 16.4 ms with the fixture on 9p and 2.6 ms on local ext4 -- a 6x
# penalty attributable entirely to the filesystem, because Whitebox touches its
# working directory inside the tool call. RichDEM (11.0 vs 12.5 ms) and
# TopoToolbox (1.96 vs 2.15 ms) were unaffected.
#
# A campaign run straight off such a mount would report that lsdtt3 comfortably
# beats Whitebox when the real difference is the filesystem. Fixtures are
# therefore STAGED to local storage before anything is measured.
SLOW_FILESYSTEMS = {
    "9p", "v9fs",       # WSL2 /mnt/c
    "drvfs", "9pfs",    # WSL1/WSL2 Windows drive mounts
    "cifs", "smb2",     # network shares
    "nfs", "nfs4",      # network filesystems
    "fuseblk",          # ntfs-3g
    "overlayfs",        # some container layers
}


def filesystem_type(path: pathlib.Path) -> str:
    """Filesystem type for a path, or "unknown".

    `stat -f` is used rather than parsing /proc/mounts because it resolves the
    mount for the path directly, including bind mounts, and it works on macOS.
    """
    try:
        r = subprocess.run(["stat", "-f", "-c", "%T", str(path)],
                           capture_output=True, text=True, check=False)
        if r.returncode == 0 and r.stdout.strip():
            return r.stdout.strip()
    except (FileNotFoundError, PermissionError):
        pass
    return "unknown"


def choose_scratch_root(explicit: str | None) -> pathlib.Path:
    """Where to stage fixtures so measurements are not filesystem-bound.

    Order: an explicit --scratch-dir, then $TMPDIR / the system temp directory.
    The system temp directory is local disk on essentially every machine this
    will run on, including cluster nodes, where it is usually node-local and
    therefore the right place anyway.
    """
    import tempfile

    if explicit:
        root = pathlib.Path(explicit).expanduser()
    else:
        root = pathlib.Path(tempfile.gettempdir()) / "lsdtt3_benchmark_fixtures"
    root.mkdir(parents=True, exist_ok=True)
    return root


def stage_fixture(src: pathlib.Path, scratch_root: pathlib.Path) -> pathlib.Path:
    """Copy a fixture to local storage, unless it is already there.

    Cached by name and size so repeated campaigns do not re-copy a 20 MB
    raster. Content is not hashed for the cache check: the campaign records the
    SHA of the STAGED file it actually measured, so a stale copy would be
    visible as a hash change rather than passing silently.
    """
    dst = scratch_root / src.name
    if dst.exists() and dst.stat().st_size == src.stat().st_size:
        return dst
    shutil.copy2(src, dst)
    return dst


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
        # Two shapes, because the lsdtt3 CLIs print two.
        #
        #   "  Compute time: 0.0026 s"        the common form
        #   "  Write time: 0.0090 s"
        #   "  Computed hillshade... done (0.0035 s)"   surface-metrics only
        #   "  Running direct fit... done (0.0041 s)"
        #
        # The "... done (X s)" phrasing is NOT tidied away in the CLI, because
        # run_surface_benchmark.py's FIT_TIMING_RE parses that exact wording
        # and deliberately raises rather than falling back to wall time if it
        # stops matching. Changing the CLI to suit this parser broke that
        # driver; the parser accommodates the CLI instead.
        COMPUTE_TIME_RE = re.compile(
            r"^\s*(?:(Compute|Write) time:\s*([0-9.eE+-]+)\s*s"
            r"|(?:Running \w+ fit|Computed \w+)\.\.\. done "
            r"\(([0-9.eE+-]+) s\))\s*$",
            re.MULTILINE)

    compute, write = [], []
    for kind, value, done_value in COMPUTE_TIME_RE.findall(stdout):
        raw = value or done_value
        try:
            v = float(raw)
        except ValueError:
            continue
        # A "... done (X s)" line is always a compute time; only the explicit
        # "Write time:" form is a write.
        (write if kind == "Write" else compute).append(v)

    return {
        "compute": sum(compute) if compute else None,
        "write": sum(write) if write else None,
    }


def measure_lsdtt(case_name: str, case: dict, fixture: pathlib.Path,
                  work_dir: pathlib.Path, threads: int, repeats: int,
                  build_dir: pathlib.Path,
                  timeout: float | None = None) -> dict[str, Any] | None:
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

    r = run_with_rusage(cmd, env=env, timeout=timeout)   # warm-up, discarded
    if r.get("timed_out"):
        return {"timed_out": True, "timeout_seconds": timeout}
    if r["returncode"] != 0:
        return {"error": f"exit {r['returncode']}",
                "stderr": r["stderr"][-1500:], "stdout": r["stdout"][-1500:]}

    timings, rss, computes, writes = [], [], [], []
    for _ in range(max(1, repeats)):
        r = run_with_rusage(cmd, env=env, timeout=timeout)
        if r.get("timed_out"):
            return {"timed_out": True, "timeout_seconds": timeout}
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


def reference_env(threads: int) -> dict[str, str]:
    """Environment for a reference subprocess, with thread counts pinned.

    NOT system_env(): that strips LD_LIBRARY_PATH and CONDA_PREFIX so the C++
    binaries find the system GDAL, and the reference scripts are Python running
    inside the pixi environment -- stripping those would break their imports.
    The two need opposite treatment and must stay separate functions.

    What this shares with system_env() is the thread pinning, which is the
    point. Whitebox is pinned through `WbEnvironment.max_procs` via --threads,
    but RichDEM and TopoToolbox sit on numpy/scipy, and an OpenBLAS or MKL
    pulled in transitively will happily spin up a pool sized to the machine.
    Without these variables a "1 thread" reference row is not one thread, and
    the campaign silently compares a pinned lsdtt3 against an unpinned
    reference.
    """
    env = dict(os.environ)
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[key] = str(threads)
    return env


def measure_reference(engine: str, routine: str, fixture: pathlib.Path,
                      threads: int, repeats: int, python_exe: str,
                      timeout: float | None = None) -> dict[str, Any] | None:
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

    r = run_with_rusage(cmd, env=reference_env(threads), timeout=timeout)
    if r.get("timed_out"):
        return {"timed_out": True, "timeout_seconds": timeout}
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


def make_tiny_grid(fixture: pathlib.Path, work_dir: pathlib.Path,
                   size: int = 8) -> pathlib.Path | None:
    """An 8x8 crop of a real fixture, for measuring fixed per-call cost.

    8x8 is 64 cells: small enough that the algorithm itself is microseconds, so
    whatever a routine costs on it is essentially all fixed overhead.  A real
    crop rather than a synthetic surface, so the nodata mask, dtype and
    georeferencing are exactly what the engines see on the full fixtures.
    """
    out = work_dir / f"calib_{size}.tif"
    if out.exists():
        return out
    try:
        from osgeo import gdal
        gdal.UseExceptions()
        work_dir.mkdir(parents=True, exist_ok=True)
        gdal.Translate(str(out), str(fixture), srcWin=[0, 0, size, size])
    except Exception:
        return None
    return out if out.exists() else None


def calibrate_overheads(build_dir: pathlib.Path, fixture: pathlib.Path,
                        work_dir: pathlib.Path,
                        python_exe: str) -> dict[str, Any]:
    """Measure the fixed per-call cost of EVERY engine, not just lsdtt3.

    WHY THIS EXISTS IN THIS FORM

    An earlier version of this suite measured a startup cost for lsdtt3 only,
    on the assumption that the three reference scripts time the routine
    in-process and therefore carry no comparable fixed cost.  That assumption
    was checked and is FALSE for one of them.  Measured on crops of
    alpine_small, median of 7, fill:

        grid     richdem   topotoolbox   whitebox(1t)
        8^2       9.83 ms      0.036 ms       0.038 ms
        16^2     11.56 ms      0.029 ms       0.032 ms
        64^2     10.19 ms      0.111 ms       0.126 ms
        256^2    13.73 ms      1.966 ms       1.440 ms

    TopoToolbox and Whitebox extrapolate to ~0.03 ms at zero cells -- no
    meaningful fixed cost.  RichDEM extrapolates to ~10 ms.  Its own internal
    timer prints "Zhou2016 wall-time = 1.2e-05 s" for the 8x8 fill, so ~99.8%
    of what a caller measures on a small grid is RichDEM's Python binding
    overhead: the per-call rdarray construction and the citation block it
    formats and prints on every single invocation.  On alpine_small that is
    roughly half of RichDEM's reported fill time.

    So a floor is measured for every engine, from the same 8x8 grid, and the
    `overhead_dominated` flag uses the engine's own floor.  Doing this for
    lsdtt3 alone would have quietly held RichDEM to a standard it does not
    meet, and RichDEM is the reference lsdtt3 is most often compared against.

    Also returns `whitebox_dispatch_seconds` as a REGRESSION TRIPWIRE for the
    tool-handle bug; see whitebox_reference.py's `_make_call`.
    """
    out: dict[str, Any] = {}

    exe = build_dir / "lsdtt-raster-preprocessing"
    if exe.exists():
        # No parameter file: the CLI prints usage and exits. Process start,
        # dynamic linking and GDAL driver registration, with none of the work.
        # This still matters for any row that has no CLI-reported compute time.
        times = []
        for _ in range(5):
            r = run_with_rusage([str(exe)], env=system_env(1))
            times.append(r["wall_seconds"])
        out["lsdtt3_startup_seconds"] = min(times)

    tiny = make_tiny_grid(fixture, work_dir)
    if tiny is None:
        return out

    # One routine every engine implements, on 64 cells.
    for engine in ("richdem", "topotoolbox", "whitebox"):
        info = measure_reference(engine, "fill", tiny, 1, 7, python_exe)
        if info and info.get("seconds_median") is not None:
            out[f"{engine}_fixed_seconds"] = info["seconds_median"]

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
    if info.get("timed_out"):
        # A distinct status from "error": the engine was working, it was just
        # too slow for the deadline. Reading it as a crash would hide a real
        # and quotable performance fact.
        row["status"] = "timeout"
        row["note"] = (f"exceeded --timeout-seconds "
                       f"({info.get('timeout_seconds')})")
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
    # Each engine is held to ITS OWN measured fixed per-call cost.
    #
    # For lsdtt3 the floor applies only to a row timed as whole-process wall
    # clock, which genuinely contains process startup; a row timed by the CLI's
    # own reported compute time already excludes it.
    #
    # For the references the floor is the measured 8x8 cost. It is ~0 for
    # TopoToolbox and Whitebox and ~10 ms for RichDEM, whose Python binding
    # rebuilds an rdarray and reprints its citation block on every call. See
    # calibrate_overheads().
    if engine == "lsdtt3":
        floor = (None if info.get("timing_scope") == "cli_reported_compute"
                 else overhead.get("lsdtt3_startup_seconds"))
    else:
        floor = overhead.get(f"{engine}_fixed_seconds")
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

    fs = manifest.get("filesystems", {})
    if fs:
        L.append(f"- measured on: `{fs['measured_on']['type']}` at "
                 f"`{fs['measured_on']['path']}` "
                 f"({'staged from ' + fs['fixture_source']['type'] if fs.get('staged') else 'not staged'})")
        L.append("")
        if fs.get("slow"):
            L.append(f"> **Measurements were taken on a "
                     f"`{fs['measured_on']['type']}` filesystem.** Engines that "
                     "touch the filesystem inside a timed call are inflated "
                     "several-fold there -- Whitebox's 256² fill measured "
                     "16.4 ms on 9p against 2.6 ms on local ext4, a 6x penalty "
                     "that is the filesystem, not the algorithm. **Do not "
                     "quote this campaign.** Re-run with `--scratch-dir` "
                     "pointing at local storage.")
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
    L.append("`*_fixed_seconds` is each engine's per-call cost measured on an "
             "8x8 crop (64 cells), where the algorithm itself is microseconds. "
             "`lsdtt3_startup_seconds` is process start, dynamic linking and "
             "GDAL driver registration, and applies only to rows with no "
             "CLI-reported compute time. A row below "
             f"{MIN_TRUSTWORTHY_FACTOR:g}x **its own engine's** floor is "
             "flagged `overhead_dominated` and excluded from the ratio table.")
    L.append("")
    rd_fixed = overhead.get("richdem_fixed_seconds")
    if rd_fixed and rd_fixed > 0.002:
        L.append(f"> **RichDEM carries {rd_fixed * 1000:.1f} ms of fixed "
                 "per-call cost** -- a per-call `rdarray` construction plus the "
                 "citation block it prints on every invocation. Its own timer "
                 "reports the 8x8 fill as ~1.2e-05 s. On a small grid most of "
                 "a RichDEM measurement is its Python binding, not its "
                 "algorithm; subtract this figure before reading a RichDEM row "
                 "as an algorithm cost.")
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
        L.append("These rows come from a CLI that printed no `Compute time:` "
                 "line, so the only figure available is whole-process wall "
                 "clock: process start, GDAL read, the routine, the write and "
                 "teardown. Against a reference timed in-process that is not a "
                 "comparison, so they are reported alone. Subtract roughly "
                 "`lsdtt3_startup_seconds` above to estimate the routine, and "
                 "treat the result as an estimate.")
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

    bad = [r for r in rows
           if r["status"] in ("error", "skipped", "unavailable", "timeout")]
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
                    help="campaign root; one subdirectory per run. A relative "
                         "path is resolved against the CURRENT DIRECTORY. Note "
                         "pixi tasks run from the repo root, so pass an "
                         "absolute path to collect results elsewhere.")
    ap.add_argument("--fixture-dir", default="",
                    help="where the downloaded region sources and their crops "
                         "live. Defaults to tests/data/dem_cache inside the "
                         "repo, which is the wrong place on a cluster: the "
                         "server profile stores ~512 MB here and a repo on a "
                         "quota'd home filesystem will not take it. Point this "
                         "at project or scratch storage. Reused across "
                         "campaigns, so it is worth keeping.")
    ap.add_argument("--python-exe", default=sys.executable,
                    help="interpreter for the reference scripts; must be the "
                         "one with richdem/topotoolbox/whitebox installed")
    ap.add_argument("--suites", default="",
                    help="comma-separated subset of suites (default: profile)")
    ap.add_argument("--cases", default="",
                    help="comma-separated subset of case names")
    ap.add_argument("--regions", default="",
                    help=f"comma-separated regions, overriding the profile "
                         f"(available: {', '.join(sorted(REGIONS))})")
    ap.add_argument("--ladder", default="",
                    help="comma-separated crop sizes in pixels, overriding the "
                         "profile. The source raster is fetched at the largest "
                         "of these.")
    ap.add_argument("--fixtures", default="",
                    help="bypass the region/ladder system entirely and use "
                         "these raster names from the fixture cache verbatim. "
                         "The escape hatch for hand-placed or legacy rasters.")
    ap.add_argument("--threads", default="",
                    help="comma-separated thread counts, overriding the profile")
    ap.add_argument("--repeats", type=int, default=0,
                    help="timed repeats, overriding the profile")
    ap.add_argument("--timeout-seconds", type=float, default=1800.0,
                    help="kill any single measurement that exceeds this and "
                         "record it as a `timeout` row. Default 1800. The "
                         "engines are not uniformly fast -- Whitebox's breach "
                         "on low-relief terrain measured ~300x slower than "
                         "lsdtt3's -- so an unattended campaign needs a "
                         "deadline. Set 0 to disable.")
    ap.add_argument("--scratch-dir", default="",
                    help="local directory to stage fixtures and intermediate "
                         "rasters into. Defaults to a subdirectory of the "
                         "system temp directory. Set this if your temp "
                         "directory is not local disk.")
    ap.add_argument("--no-stage", action="store_true",
                    help="measure fixtures where they lie instead of staging "
                         "them to local storage. Only for checking what "
                         "staging is worth: on a slow mount it inflates some "
                         "engines several-fold (see SLOW_FILESYSTEMS).")
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

    # Line-buffer our own stdout.
    #
    # Python block-buffers stdout (~8 KB) whenever it is not a terminal, which
    # is exactly the case for every way this script is actually run for real:
    # `nohup ... &`, a Slurm job log, or any redirect to a file. A campaign
    # that takes hours would then show nothing at all for a long stretch and
    # look hung, when it is simply running. Progress output is the only way to
    # tell a slow measurement from a stalled one, so it has to arrive as it
    # happens.
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except (AttributeError, ValueError):
        # Older interpreters, or a stdout that cannot be reconfigured. Not
        # worth failing a campaign over -- the run is still correct, the log is
        # just chunky.
        pass

    # Rebind before ANY fixture helper runs -- they all reach for the module
    # global rather than taking it as an argument.
    global FIXTURE_DIR
    if args.fixture_dir:
        FIXTURE_DIR = pathlib.Path(args.fixture_dir).expanduser().resolve()

    profile = dict(PROFILES[args.profile])
    explicit_fixtures = ([f.strip() for f in args.fixtures.split(",")
                          if f.strip()] if args.fixtures else None)
    if args.regions:
        profile["regions"] = [r.strip() for r in args.regions.split(",")
                              if r.strip()]
    if args.ladder:
        profile["ladder"] = [int(v) for v in args.ladder.split(",")
                             if v.strip()]
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
        # Name the valid values. A bare "no cases selected" gives a user who
        # mistyped a case or a suite nothing to correct against, and the case
        # names are not guessable from the docs alone.
        print("ERROR: no cases selected.", file=sys.stderr)
        print(f"  suites: {sorted({c['suite'] for c in CASES.values()})}",
              file=sys.stderr)
        print(f"  cases:  {sorted(CASES)}", file=sys.stderr)
        print("  (--cases names must also lie within the selected suites)",
              file=sys.stderr)
        return 2

    # An "8 thread" row on a 4-core machine is not a slow result, it is a
    # mislabelled one, and it is indistinguishable from a genuine scaling
    # ceiling. The Slurm wrapper enforces this for its allocation; enforce it
    # here too, because most runs will not go through Slurm.
    cores = os.cpu_count() or 1
    if max(profile["threads"]) > cores:
        print(f"WARNING: --threads goes up to {max(profile['threads'])} but "
              f"this machine reports {cores} logical cores.\n"
              f"         Rows above {cores} threads will be oversubscribed and "
              f"their scaling figures are not meaningful.\n"
              f"         Pass --threads with a smaller list to avoid this.",
              file=sys.stderr)

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
        if explicit_fixtures:
            print(f"fixtures: {explicit_fixtures} (explicit, ladder bypassed)")
        else:
            src = max(profile["ladder"])
            print(f"regions:  {profile['regions']}")
            print(f"ladder:   {profile['ladder']} px "
                  f"(source fetched at {src} px, "
                  f"~{src * GRID_SPACING_M / 1000:.0f} km box, "
                  f"~{src * src * 4 / 1e6:.0f} MB per region)")
        print(f"threads:  {profile['threads']}   repeats: {profile['repeats']}")
        print(f"build:    {build_dir}")
        print("")
        # Counted, not fetched: --dry-run must never touch the network.
        n_fixtures = (len(explicit_fixtures) if explicit_fixtures
                      else len(profile["regions"]) * len(profile["ladder"]))
        for name, case in sorted(selected_cases.items()):
            engines = ["lsdtt3"] + sorted(case["engines"])
            per_engine = []
            for e in engines:
                t = 1 if e in SERIAL_REFERENCES else len(profile["threads"])
                per_engine.append(f"{e}x{t}")
                n += t * n_fixtures
            print(f"  {case['suite']:11s} {name:22s} {' '.join(per_engine)}")
        print(f"\n{n} measurements, each "
              f"{profile['repeats']} repeats + 1 warm-up")
        return 0

    stamp = now_stamp()
    # Relative to the CURRENT DIRECTORY, not the repo. Resolving a relative
    # path against anything other than cwd is surprising, and here it was
    # actively wrong: running from a data directory silently deposited the
    # campaign inside the git working tree instead.
    #
    # Note pixi tasks run with cwd set to the manifest root, so
    # `pixi run bench-server --out-dir results` still lands in the repo. To
    # collect results elsewhere pass an ABSOLUTE --out-dir, or invoke the
    # script directly from the directory you want.
    run_dir = (pathlib.Path(args.out_dir).expanduser().resolve() /
               f"{args.profile}_{socket.gethostname()}_{stamp}")
    # Created AFTER staging decides where work goes: staging relocates the work
    # directory to node-local scratch, and creating it here unconditionally
    # left an empty `work/` behind in every campaign directory.
    work_dir = run_dir / "work"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Both paths, up front and unconditionally.
    #
    # The fixture cache defaults INSIDE THE REPOSITORY, while the campaign
    # defaults to the current directory -- two different roots, which is
    # genuinely surprising and has already cost someone a confused hunt for
    # DEMs that were never going to be next to their results. The sbatch
    # wrapper overrides the fixture path to sit beside the campaign; a direct
    # run does not. Printing both makes the difference visible in the first
    # three lines of any log instead of being something you have to know.
    print(f"campaign: {run_dir}")
    print(f"fixtures: {FIXTURE_DIR}"
          + ("" if args.fixture_dir else "   (default: inside the repo -- "
                                         "use --fixture-dir to move it)"))
    print(f"profile:  {args.profile} ({profile['description']})")

    hardware = capture_hardware_manifest(build_dir, args.python_exe)
    (run_dir / "hardware_manifest.json").write_text(
        json.dumps(hardware, indent=2))

    # Resolve fixtures up front so a download does not interrupt the timings.
    fetch_exe = build_dir / "lsdtt-fetch-raster"
    fixtures = resolve_fixtures(profile, fetch_exe, explicit_fixtures)
    if not fixtures:
        print("ERROR: no fixtures available. Build lsdtt-fetch-raster so the "
              "suite can fetch the region sources, or pass --fixtures naming "
              f"rasters already in {FIXTURE_DIR}.", file=sys.stderr)
        return 2

    # Stage onto local disk BEFORE any measurement. See SLOW_FILESYSTEMS.
    source_fs = filesystem_type(FIXTURE_DIR)
    scratch_root = choose_scratch_root(args.scratch_dir or None)
    scratch_fs = filesystem_type(scratch_root)

    if args.no_stage:
        print(f"  staging disabled; measuring on {source_fs}")
        work_dir.mkdir(parents=True, exist_ok=True)
    else:
        print(f"  staging fixtures: {source_fs} -> {scratch_root} "
              f"({scratch_fs})")
        for key in list(fixtures):
            fixtures[key] = stage_fixture(fixtures[key], scratch_root)
        # The work directory holds the filled DEMs the accumulation cases read,
        # so it has to be local too or the staging is undone for those cases.
        work_dir = scratch_root / f"work_{stamp}"
        work_dir.mkdir(parents=True, exist_ok=True)

    filesystems = {
        "fixture_source": {"path": str(FIXTURE_DIR), "type": source_fs},
        "measured_on": {"path": str(fixtures[next(iter(fixtures))].parent),
                        "type": filesystem_type(
                            fixtures[next(iter(fixtures))].parent)},
        "staged": not args.no_stage,
        "slow": filesystem_type(
            fixtures[next(iter(fixtures))].parent) in SLOW_FILESYSTEMS,
    }
    if filesystems["slow"]:
        print(f"  WARNING: measuring on {filesystems['measured_on']['type']}, "
              f"which is slow for engines that touch the filesystem inside a "
              f"timed call. Pass --scratch-dir to point at local storage.")

    # `terrain` is what makes a ladder honest: centred crops of increasing size
    # reach into different landscape, so relief and depression density are
    # recorded per rung and a reader can see whether a scaling curve is
    # confounded instead of assuming it is not.
    fixture_meta = {k: {"sha256_12": sha_prefix(p), "n_cells": n_cells_of(p),
                        "path": str(p), "terrain": describe_fixture(p)}
                    for k, p in fixtures.items()}

    print("calibrating fixed overhead...")
    smallest = min(fixtures.values(),
                   key=lambda p: n_cells_of(p) or (1 << 62))
    overhead = calibrate_overheads(build_dir, smallest, work_dir,
                                   args.python_exe)
    for k, v in overhead.items():
        if isinstance(v, float):
            print(f"  {k}: {v * 1000:.1f} ms")

    manifest = {
        "created_at": dt.datetime.now().isoformat(timespec="seconds"),
        # The exact invocation. A campaign is only reproducible if a reader can
        # see what was actually asked for -- every profile setting can be
        # overridden from the command line, so `profile_settings` alone does
        # not say whether a value came from the profile or from a flag.
        "invocation": {
            "argv": sys.argv,
            "cwd": os.getcwd(),
            "python": sys.executable,
        },
        "profile": args.profile,
        "profile_settings": profile,
        "build_dir": str(build_dir),
        "fixtures": fixture_meta,
        "overhead": overhead,
        "hardware": hardware,
        "filesystems": filesystems,
        "min_trustworthy_factor": MIN_TRUSTWORTHY_FACTOR,
        "logical_cores": cores,
        "oversubscribed": max(profile["threads"]) > cores,
    }

    rows: list[dict[str, Any]] = []
    # Matches --dry-run exactly: lsdtt3 runs at every thread count, a serial
    # reference runs once, a parallel reference runs at every thread count.
    baseline_threads = min(profile["threads"])
    measurement_timeout = args.timeout_seconds if args.timeout_seconds > 0 else None
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
                    # Record the failure for EVERY engine, not just lsdtt3.
                    # All of them accumulate over this shared raster, so none
                    # of them can run; emitting one lsdtt3 row made the
                    # reference engines disappear from the case with no
                    # explanation anywhere in the CSV.
                    err = {"error": "could not prepare filled DEM "
                                    "(lsdtt-raster-preprocessing missing "
                                    "or failed)"}
                    for eng in ["lsdtt3"] + sorted(case["engines"]):
                        rows.append(build_row(
                            case["suite"], case_name, eng, fx_key,
                            meta["sha256_12"], meta["n_cells"],
                            baseline_threads, err, case, overhead))
                    continue
                in_path = filled

            for threads in profile["threads"]:
                # lsdtt3
                info = measure_lsdtt(case_name, case, in_path, work_dir,
                                     threads, profile["repeats"], build_dir,
                                     timeout=measurement_timeout)
                rows.append(build_row(case["suite"], case_name, "lsdtt3",
                                      fx_key, meta["sha256_12"],
                                      meta["n_cells"], threads, info, case,
                                      overhead))
                done += 1
                _progress(done, total, case_name, "lsdtt3", fx_key, threads,
                          rows[-1])

                for engine, routine in sorted(case["engines"].items()):
                    # A serial reference is measured ONCE, at the baseline
                    # thread count. Running it again at 2/4/8 would record the
                    # same number under a different label and imply it had been
                    # tested there.
                    #
                    # Against `baseline_threads`, not a hard-coded 1: with
                    # `--threads 2,4` the old test never fired, so RichDEM and
                    # TopoToolbox produced NO rows at all and vanished from the
                    # campaign without a word.
                    if engine in SERIAL_REFERENCES and threads != baseline_threads:
                        continue
                    info = measure_reference(engine, routine, in_path,
                                             threads, profile["repeats"],
                                             args.python_exe,
                                             timeout=measurement_timeout)
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
    # directory too big to send. The work directory may have been relocated to
    # local scratch by the staging step, so remove whichever one was used.
    #
    # The STAGED FIXTURES are deliberately left in place: they are the
    # expensive part to recreate, they are cached across campaigns by name and
    # size, and they live under the system temp directory, which the OS clears
    # on its own schedule.
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
