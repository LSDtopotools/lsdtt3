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
  than what a Python-level probe can see.  Every command is started through
  rusage_launcher.py, because a Linux child inherits its launcher's lifetime
  peak RSS; the resulting floor is measured and tripwired every campaign
  (`rss_floor` in the manifest, see calibrate_rss_floor()).

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
import lsdtt_test_env

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent

# WHERE THINGS ARE WRITTEN
#
# Three separate locations, because they have three different lifetimes and
# three different size and speed requirements:
#
#   FIXTURE_DIR   the downloaded region sources and their crops. Large
#                 (~300 MB for the server profile), reused across every
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

# Parallel-capable references whose thread sweep is CAPPED, engine -> highest
# thread count measured. lsdtt3 still sweeps the full profile list.
#
# Whitebox did not scale on any surface case on forth02 (2026-09-09): 0.96x to
# 1.04x at 2, 4 and 8 threads for slope, hillshade and curvature, every
# fixture. Separately, a laptop check found `max_procs=1` already uses ~1.3
# cores for d8_flow_accum and ~1.9 for slope, and `max_procs=-1` changes
# neither. Measuring it at 4 and 8 threads repeated the same number under
# larger labels for ~12 min per campaign. 2 threads is kept so the flat curve
# stays visible in the scaling table rather than being assumed.
REFERENCE_MAX_THREADS: dict[str, int] = {"whitebox": 2}


def _engine_thread_count(engine: str, case: dict, threads: list[int]) -> int:
    """How many thread counts `engine` is measured at for `case`.

    Mirrors the measurement loop in main(), which --dry-run and the progress
    counter must agree with."""
    if case.get("sequential") or engine in SERIAL_REFERENCES:
        return 1
    base = min(threads)
    cap = REFERENCE_MAX_THREADS.get(engine)
    if cap is None:
        return len(threads)
    return sum(1 for t in threads if t == base or t <= cap)


# REUSE_FROM_NOTE
#
# A case may carry `reuse_from: {engine: source_case}`. For those engines the
# row measured for `source_case` (same fixture, same thread count) is COPIED
# into this case instead of being measured again, with `reused_from` set in the
# CSV. Only for engines whose command is byte-for-byte identical between the
# two cases -- that is the whole justification, and each use states it.
#
# Source cases run first (see the case ordering in main()). If the source case
# is not selected in this run, or produced no row, the engine is measured
# normally, so `--cases` subsets never lose a row.

# Cases where lsdtt3's OWN implementation is inherently sequential too, so no
# engine in the case has anything to gain from more threads.
#
# Depression filling (priority-flood and grayscale reconstruction alike),
# least-cost breaching, and D8/D-inf/MFD graph construction are ordered
# traversals: each step depends on the result of the one before it. The first
# server campaign measured every one of them flat -- 0.92x to 1.09x across
# 1/2/4 threads, i.e. noise about 1.00x -- which is the CORRECT answer for a
# sequential algorithm, not a defect.
#
# Sweeping threads over them therefore measures nothing. It multiplies the
# campaign length for those cases and, worse, fills the scaling table with flat
# rows that invite a reader to conclude something is broken. They are measured
# ONCE, at the baseline thread count, and flagged `sequential_case` in the CSV.
#
# Do NOT add a case here merely because its scaling looks poor. Poor scaling in
# a genuinely parallel kernel is a finding worth keeping -- the surface metrics
# are memory-bandwidth-bound and stay in the sweep for exactly that reason.
# This list is only for algorithms with no parallelism to measure.
SEQUENTIAL_CASES_NOTE = (
    "marked per-case with `sequential: True`; see the comment above this line"
)

REFERENCE_SCRIPTS = {
    "richdem": "richdem_reference.py",
    "topotoolbox": "topotoolbox_reference.py",
    "whitebox": "whitebox_reference.py",
}


# -----------------------------------------------------------------------------
# WHY AN ENGINE IS ABSENT FROM A CASE.
#
# An engine missing from a case's `engines` map produces no row, and the summary
# used to print "-" for that cell -- the SAME mark it prints when a run was
# skipped or failed. For a peer-reviewed table those are very different claims:
# "not reachable from the interface we benchmark" is a result; "we did not
# measure it" is a gap.
#
# *** THE CLAIM IS ABOUT THE INTERFACE, NOT ABOUT THE PROJECT. ***
# This benchmark drives TopoToolbox through its PYTHON package, which is early
# (0.0.12) and exposes a fraction of what TopoToolbox 3 offers in MATLAB.
# TopoToolbox 3 DOES implement multiple-flow-direction and D-infinity routing
# -- see md/usersguide_3_mfd.md in TopoToolbox/topotoolbox3 -- they are simply
# not bound in the Python package yet. Writing "TopoToolbox does not implement
# MFD" would be false, and a reviewer who uses the MATLAB version would know it
# at a glance. Word every reason as a statement about the BINDING.
#
# The same care applies to the others: RichDEM and Whitebox are also driven
# through Python bindings, and lsdtt3 through its CLI.
#
# Every entry has been CHECKED against the installed package. Do not add one on
# the strength of a recollection -- the registry previously asserted "RichDEM
# has no slope tool", which is false (rd.TerrainAttribute offers slope_riserun,
# slope_percentage, slope_degrees, slope_radians, planform_curvature and
# profile_curvature). That claim silently kept a whole reference out of the
# surface comparison.
#
# Verified 2026-09-08 against richdem 2.3.0, topotoolbox 0.0.12,
# whitebox_workflows 2.0.6.
NOT_IMPLEMENTED: dict[str, dict[str, str]] = {
    "fill": {
        "whitebox": "no grayscale-reconstruction fill in whitebox_workflows "
                    "2.0.6 (hydrology offers fill_depressions [priority-flood, "
                    "measured in fill_priority_flood], "
                    "fill_depressions_wang_and_liu, "
                    "fill_depressions_planchon_and_darboux, fill_pits)",
    },
    "breach": {
        "whitebox": "no Lindsay (2016) complete breaching in whitebox_workflows "
                    "2.0.6 (hydrology offers breach_depressions_least_cost "
                    "[measured in breach_least_cost] and "
                    "breach_single_cell_pits)",
    },
    "breach_least_cost": {
        "richdem": "no least-cost breaching in richdem 2.3.0 "
                   "(rd.BreachDepressions is Lindsay 2016 complete breaching, "
                   "measured in breach)",
    },
    "accum_dinf": {
        "topotoolbox": "not bound in the Python package 0.0.12 "
                       "(FlowObject takes method='d8' only); TopoToolbox 3 "
                       "implements D-infinity in MATLAB",
    },
    "accum_mfd": {
        "topotoolbox": "not bound in the Python package 0.0.12 "
                       "(FlowObject takes method='d8' only); TopoToolbox 3 "
                       "implements MFD in MATLAB",
    },
    "fill_epsilon": {
        "topotoolbox": "not bound in the Python package 0.0.12 "
                       "(fillsinks(bc, hybrid) takes no flat-increment "
                       "argument)",
    },
    "hillshade": {
        "richdem": "no hillshade in richdem 2.3.0 (rd.TerrainAttribute's "
                   "attribute set is slope_*/aspect/curvature only)",
    },
}


# -----------------------------------------------------------------------------
# Case registry
#
# One entry per comparable operation.  `engines` maps an engine name to the
# routine string that engine understands; an engine absent from the map simply
# produces no row.
#
# *** THAT ABSENCE IS NOT SELF-DOCUMENTING, AND ONE OF THEM WAS WRONG. ***
# This comment used to give "RichDEM has no slope tool" as its example. RichDEM
# 2.3.0 HAS one: rd.TerrainAttribute offers slope_riserun, slope_percentage,
# slope_degrees, slope_radians, aspect, curvature, planform_curvature and
# profile_curvature. The surface comparison ran with three engines where four
# were available, for two campaigns, on the strength of a code comment nobody
# rechecked. RichDEM is now wired into `slope_only` and `polyfit_all`
# (2026-09-09) and the routines are exercised in richdem_reference.py.
#
# Where an engine genuinely cannot do a case, record it in NOT_IMPLEMENTED
# above WITH the reason, and check it against the installed package first --
# an unexplained absence reads as "we did not bother", which is the more
# damaging of the two readings for a paper.
#
# `lsdtt` describes how to drive the lsdtt3 CLI: which executable, which
# parameters, and which output file to expect.  These mirror run_flow_benchmark
# .py's COMPARISONS table; they are restated rather than imported because that
# table is shaped for correctness assertions and this one for timing, and
# entangling them would mean a tolerance change silently altering a benchmark.
# -----------------------------------------------------------------------------
CASES: dict[str, dict[str, Any]] = {
    # -- suite: flow -------------------------------------------------------
    # -- fill: TWO cases, because there are two algorithm families ---------
    #
    # Every tool here implements both a grayscale-reconstruction fill (raster
    # scans + FIFO queue, O(N), sequential access) and a priority-flood fill
    # (binary heap, O(N log N)). They produce the SAME surface and differ only
    # in cost, so a fill benchmark must pair like with like.
    #
    # An earlier version of this table did not. It ran lsdtt3's DEFAULT
    # (`barnes2014`, the heap) against TopoToolbox's DEFAULT (`hybrid=True`,
    # the reconstruction) and reported lsdtt3 2.9x slower. Measured on
    # alpine_1000 with fill_min_slope 0, bit-identical output:
    #
    #     lsdtt3 barnes2014    0.1282 s    7.80 Mcells/s
    #     lsdtt3 vincent1993   0.0446 s   22.42 Mcells/s   <- 2.87x
    #     topotoolbox hybrid                22.84 Mcells/s
    #
    # Matched, lsdtt3 is within 2% of TopoToolbox. The 2.9x was the algorithm
    # choice, not the implementation. See docs/developer_notes/testing_and_benchmarking.md.
    "fill": {
        "sequential": True,
        "suite": "flow",
        "note": ("RECONSTRUCTION family: grayscale reconstruction by erosion "
                 "(Vincent 1993 hybrid -- raster scans + FIFO queue, O(N)). "
                 "lsdtt3 vincent1993 vs TopoToolbox hybrid=True"),
        # Whitebox was paired here until 2026-09-13, but its fill_depressions
        # is a PRIORITY-FLOOD tool (its own tool description: "priority-flood
        # strategy"), so it compared a heap fill with a reconstruction fill --
        # the exact mismatch that once reported lsdtt3 2.9x slow. It is now in
        # `fill_priority_flood`. Whitebox has no reconstruction fill; see
        # NOT_IMPLEMENTED.
        "engines": {"topotoolbox": "fill"},
        "lsdtt": {"exe": "lsdtt-raster-preprocessing",
                  "params": {"write_fill_raster": "true",
                             "fill_min_slope": "0",
                             # NOT the CLI default. Reconstruction is what
                             # TopoToolbox's default runs, and it cannot impose
                             # a gradient, hence min_slope 0.
                             "fill_algorithm_choice": "vincent1993"},
                  "output": "_fill.tif"},
    },
    "fill_priority_flood": {
        "sequential": True,
        "suite": "flow",
        "note": ("PRIORITY-FLOOD family: heap-based flood, O(N log N). "
                 "lsdtt3 barnes2014 vs TopoToolbox hybrid=False vs RichDEM "
                 "(Zhou 2016) vs Whitebox fill_depressions. Same surface as "
                 "the reconstruction case, roughly 2-3x the cost in every "
                 "implementation"),
        # Whitebox `fill` = fill_depressions(fix_flats=False), a priority-flood
        # fill (moved here from `fill` 2026-09-13).
        "engines": {"richdem": "fill", "topotoolbox": "fill_priority_flood",
                    "whitebox": "fill"},
        "lsdtt": {"exe": "lsdtt-raster-preprocessing",
                  "params": {"write_fill_raster": "true",
                             "fill_min_slope": "0",
                             "fill_algorithm_choice": "barnes2014"},
                  "output": "_fill.tif"},
    },
    "fill_epsilon": {
        "sequential": True,
        "suite": "flow",
        "note": ("fill with a gradient imposed across flats. PRIORITY-FLOOD "
                 "ONLY in every tool: reconstruction fills to exactly the "
                 "spill level and has no mechanism for a gradient, so there "
                 "is no reconstruction row to compare against here. Whitebox "
                 "runs fill_depressions_wang_and_liu with flat_increment "
                 "matched to fill_min_slope * cell size"),
        # Whitebox: the Wang & Liu tool, not `fill_epsilon` (its Barnes-style
        # fill_depressions), which disagrees with three other implementations
        # by up to 9.78 m (benchmarking plan section 2.2) -- and until
        # 2026-09-13 it also ran at Whitebox's default flat_increment of 1e-4 m
        # per cell against lsdtt3's 1e-4 GRADIENT (3e-3 m per cell at 30 m),
        # a 30x different amount of flat work. The increment is now matched
        # per fixture: see `whitebox_flat_increment_from` and
        # reference_extra_args().
        "engines": {"richdem": "fill_epsilon",
                    "whitebox": "depression_fill_wl_epsilon"},
        "whitebox_flat_increment_from": "fill_min_slope",
        "lsdtt": {"exe": "lsdtt-raster-preprocessing",
                  "params": {"write_fill_raster": "true",
                             "fill_min_slope": "0.0001",
                             "fill_algorithm_choice": "barnes2014"},
                  "output": "_fill.tif"},
    },
    "breach": {
        "sequential": True,
        "suite": "flow",
        "note": "complete breaching, Lindsay (2016): lsdtt3 vs RichDEM",
        # Whitebox was paired here until 2026-09-13, but it runs
        # breach_depressions_least_cost -- Lindsay (2020), a different
        # algorithm with a bounded search. It is now compared in
        # `breach_least_cost` against lsdtt3's own least-cost breach.
        "engines": {"richdem": "breach"},
        "lsdtt": {"exe": "lsdtt-raster-preprocessing",
                  "params": {"write_breached_raster": "true",
                             "breach_algorithm_choice": "lindsay2016"},
                  "output": "_breached.tif"},
    },
    # Added 2026-09-13 so Whitebox's breach is paired with the SAME algorithm.
    # Matched explicitly: max search distance 100 cells (whitebox_reference.py
    # --breach-max-dist default), unresolved depressions filled afterwards
    # (Whitebox fill_deps=True, lsdtt3 fill_remainder_after_breach, default
    # true), no cost cap (lsdtt3 default; Whitebox's max_cost left unset).
    # NOT verified equal, because whitebox_workflows 2.0.6 does not report its
    # defaults: minimize_dist (lsdtt3 min_dist=true) and flat_increment
    # (lsdtt3 0). A timing row, not a value comparison.
    "breach_least_cost": {
        "sequential": True,
        "suite": "flow",
        "note": ("least-cost breaching, Lindsay (2020), max_dist 100, "
                 "unresolved depressions filled: lsdtt3 vs Whitebox"),
        "engines": {"whitebox": "breach"},
        "lsdtt": {"exe": "lsdtt-raster-preprocessing",
                  "params": {"write_breached_raster": "true",
                             "breach_algorithm_choice": "least_cost",
                             "breach_lc_max_dist": "100",
                             "fill_remainder_after_breach": "true"},
                  "output": "_breached.tif"},
    },
    "accum_d8": {
        "sequential": True,
        "needs_filled": True,
        "suite": "flow",
        "note": ("D8 accumulation, O'Callaghan & Mark (1984). NOT fully "
                 "like-for-like: lsdtt3 also builds an S-vector the references "
                 "do not build (~57% of its build time), retained by design "
                 "because downstream analyses depend on it -- see the note in "
                 "CASES and docs/developer_notes/testing_and_benchmarking.md"),
        "engines": {"richdem": "accum_d8", "topotoolbox": "accum_d8",
                    "whitebox": "accum_d8"},
        "lsdtt": {"exe": "lsdtt-drainage-analysis",
                  "params": {"write_d8_drainage_area_raster": "true",
                             "raster_is_filled": "true"},
                  "output": "_d8_drainage_area.tif"},
    },

    # The like-for-like partner to `accum_d8`. Same operation, same output --
    # pinned bit-identical by P11 in tests/cpp/test_flow_properties.cpp -- with
    # the s-vector switched off, which is the structure no reference builds.
    #
    # REPORT BOTH ROWS OR NEITHER. `accum_d8` is what a user of lsdtt3 actually
    # pays and is the honest headline; this row is what a fair comparison
    # against tools that build no equivalent structure costs. Quoting only this
    # one would flatter lsdtt3 by hiding work its users do get; quoting only the
    # other charges it for work the references never do.
    "accum_d8_no_svector": {
        "sequential": True,
        "needs_filled": True,
        "suite": "flow",
        "note": ("D8 accumulation with build_flow_s_vector: false. The "
                 "like-for-like row: same raster as accum_d8, without the "
                 "s-vector (~57% of the D8 build) that no reference builds"),
        "engines": {"richdem": "accum_d8", "topotoolbox": "accum_d8",
                    "whitebox": "accum_d8"},
        # The references run EXACTLY the call they run in `accum_d8` -- none
        # of them has an s-vector to switch off -- so measuring them again only
        # measured noise (RichDEM 10.50 vs 10.17 Mcells/s on forth02). Their
        # `accum_d8` rows are copied here instead. See REUSE_FROM_NOTE.
        "reuse_from": {"richdem": "accum_d8", "topotoolbox": "accum_d8",
                       "whitebox": "accum_d8"},
        "lsdtt": {"exe": "lsdtt-drainage-analysis",
                  "params": {"write_d8_drainage_area_raster": "true",
                             "raster_is_filled": "true",
                             "build_flow_s_vector": "false"},
                  "output": "_d8_drainage_area.tif"},
    },

    # -- suite: mfd --------------------------------------------------------
    # Multi-direction routing. Every implementation uses a different
    # proportioning rule (Tarboton D-inf vs Quinn MFD vs Freeman FD8), so these
    # are TIMING rows. They are still worth measuring: the cost of spreading
    # flow across multiple receivers is the interesting quantity, not whether
    # the three agree on the split.
    "accum_dinf": {
        "sequential": True,
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
        "sequential": True,
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
    # The four engines are architecturally different here and the case list is
    # shaped to make that visible rather than to hide it.
    #
    #   lsdtt3      fits a 6-term quadratic by least squares over a moving
    #               window of `polyfit_lambda_metres` and derives EVERY metric
    #               from the one fit.
    #   RichDEM     fits a 9-term quadratic analytically to the 3x3
    #               neighbourhood (Zevenbergen & Thorne 1987) and
    #               differentiates it, in C++, serial -- and refits per
    #               attribute, one attribute per call.
    #   TopoToolbox runs scipy convolutions with the Evans 3x3 kernels, in
    #               Python, serial -- and recomputes all five convolutions per
    #               curvature type.
    #   Whitebox    runs 3x3 finite-difference operators in Rust, and likewise
    #               recomputes per curvature type.
    #
    # RichDEM is the nearest of the three references to lsdtt3 in FORM -- both
    # fit a quadratic and differentiate it -- and is still not value-comparable
    # to it: RichDEM's is an exact interpolation of a fixed 3x3, lsdtt3's a
    # least-squares fit over a window the user chooses. Curvature sign and
    # aspect handedness conventions were NOT checked across the four.
    #
    # Measured decomposition of lsdtt3's cost on alpine_1000 -- the shared fit
    # is paid ONCE and the per-metric work is what scales:
    #
    #     shared `direct` fit      82.5 ms
    #     + slope                  28.0 ms
    #     + aspect_degrees         44.3 ms
    #     + profile_curvature      22.3 ms
    #     + plan_curvature         20.4 ms
    #
    # Three consequences, all stated in the summary rather than corrected for:
    #
    # 1. Window size does NOT tilt this comparison. An earlier version of this
    #    comment claimed lsdtt3's 70 m lambda was a wider window than the
    #    references' 3x3 and that its numbers were therefore "conservative".
    #    That was WRONG: `direct` is a separable rolling-moment fit whose cost
    #    is independent of window size -- measured flat from lambda 35 m to
    #    280 m (0.0191 to 0.0193 s), a 729x range of window cells. Lambda is
    #    pinned so it is a stated parameter, not because it matters to cost.
    # 2. `polyfit_all` is where lsdtt3's one-fit-many-metrics design starts to
    #    pay and `slope_only` is where it does not -- there the full quadratic
    #    fit produces a single number. Both are measured; quoting only the
    #    first would be cherry-picking.
    # 3. The metric SETS differ between tools, so even `polyfit_all` is not
    #    like-for-like in what it produces. Hence report_only throughout.
    #
    # Values are never compared: the curvature definitions and aspect
    # handedness genuinely differ between the three.
    "slope_only": {
        "suite": "surface",
        "note": ("slope alone. lsdtt3 pays the FULL 6-coefficient quadratic "
                 "fit (~82 ms of ~110 ms on 1000^2) to produce one metric, "
                 "where TopoToolbox runs 2 convolutions and RichDEM and "
                 "Whitebox run one 3x3 operator. This is lsdtt3's WORST "
                 "surface case by construction. TIMING ONLY"),
        "engines": {"richdem": "slope", "topotoolbox": "evansslope",
                    "whitebox": "slope"},
        "lsdtt": {"exe": "lsdtt-surface-metrics",
                  "params": {"polyfit_metrics": "slope",
                             "polyfit_method": "direct",
                             # Window size is FREE for the `direct` method --
                             # measured flat from lambda 35 m to 280 m (0.0191
                             # to 0.0193 s), because it is a separable
                             # rolling-moment fit. So this value does not tilt
                             # the comparison either way, contrary to what an
                             # earlier version of this file claimed.
                             "polyfit_lambda_metres": "70.0"},
                  "output": "_direct_polyfit_slope.tif"},
        "report_only": True,
    },
    "polyfit_all": {
        "suite": "surface",
        "note": ("slope + aspect + profile + plan curvature. lsdtt3 derives "
                 "all four from ONE shared fit; the references recompute per "
                 "metric. NOT like-for-like in two ways: the metric SETS "
                 "differ (lsdtt3 slope+aspect+2 curvatures, TopoToolbox 4 "
                 "curvatures, RichDEM and Whitebox 2 curvatures), and "
                 "lsdtt3's fixed fit cost has not yet amortised at 4 "
                 "metrics. TIMING ONLY"),
        # RichDEM and Whitebox contribute TWO curvatures where TopoToolbox
        # contributes four; the metric counts differ across the row and the
        # note says so. Comparing the three reference columns against each
        # other is not meaningful either.
        "engines": {"richdem": "curvature_all2",
                    "topotoolbox": "curvature_all4",
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
    # `polyfit_many` (eight metrics from one lsdtt3 fit, lsdtt3-only) was
    # REMOVED 2026-09-13. It had no reference to compare against, cost ~44 min
    # of a server campaign -- most of it writing eight rasters (33 s of a 47 s
    # process at 8000^2) -- and the amortisation it showed is already visible
    # between `slope_only` and `polyfit_all`.
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
        "sequential": True,
        "suite": "depression",
        "note": ("lsdtt3 builds a nested depression hierarchy (Barnes 2020) "
                 "then Fill-Spill-Merge (Barnes 2021). lsdtt3-only: the "
                 "absolute cost of the hierarchy. TIMING ONLY"),
        # lsdtt3-ONLY since 2026-09-13. Whitebox used to be paired here through
        # `sink`, which only labels cells lying in depressions -- a different,
        # far cheaper operation, so the row compared nothing and cost ~22 min
        # of a server campaign. whitebox_reference.py keeps the routine because
        # calibrate_overheads() uses it as the dispatch tripwire.
        "engines": {},
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
        "sequential": True,
        "suite": "depression",
        "note": "Wang & Liu (2006) fill lineage, a second fill implementation",
        "engines": {"whitebox": "depression_fill_wl"},
        # lsdtt3 has no Wang & Liu fill. Its side of this case is the CLI
        # default, barnes2014 at fill_min_slope 0 -- byte-for-byte the command
        # `fill_priority_flood` already runs (7.63 vs 7.65 Mcells/s on
        # forth02). That row is copied here instead of re-measured.
        "reuse_from": {"lsdtt3": "fill_priority_flood"},
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
        # 4000 px at 30 m is a 120 km box, fetched ONCE and cropped into every
        # rung.
        #
        # 8000 px was DROPPED 2026-09-13. On forth02 it was 4.2 h of a 5.5 h
        # campaign, and for the surface suite, D8 and D-inf every engine's
        # throughput at 8000^2 was within 0.95-1.04x of 4000^2 with the same
        # thread-scaling curve. What it did still show -- fill-family
        # throughput falling further (0.59-0.82x of 4000^2) -- is confounded
        # with terrain, because each rung is a wider crop. Pass --ladder to add
        # it back for a one-off size study. Rungs already cropped from an 8000
        # source are reused as they are (same SHA), so campaigns stay
        # comparable; see ensure_source_raster().
        "ladder": [1000, 2000, 4000],
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

       This was found the hard way in run_flow_benchmark.py; the stripping now
       lives in lsdtt_test_env.system_env(), shared by every driver.

    2. Pin the thread count.  OMP_NUM_THREADS is what lsdtt3's OpenMP loops
       read.  OPENBLAS_NUM_THREADS and MKL_NUM_THREADS are set alongside it
       because a BLAS pulled in transitively will otherwise spin up its own
       pool and the "1 thread" row will not be one thread.
    """
    env = lsdtt_test_env.system_env()
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

    *** THE COMMAND IS LAUNCHED THROUGH rusage_launcher.py, NOT DIRECTLY. ***
    A Linux child inherits its parent's lifetime peak RSS and wait4 reports it
    as the child's own, so launching straight from this driver made every
    measurement report at least this driver's peak. The 2026-09-09 forth02
    campaign lost its entire memory column to that: 441 of 504 rows read
    exactly 3556.98 MiB, the peak describe_fixture() had reached. The launcher
    is a fresh small process, so the floor under each measurement is ITS peak,
    which calibrate_overheads() measures as `rss_floor_mb`. See the launcher's
    docstring for the kernel mechanism and the measurement that confirmed it.

    The launcher does the wait4 and writes the result to a temp file; its wall
    time brackets only the measured command, so the Python start-up of the
    launcher itself is not charged to `wall_seconds`.

    stderr goes to a temp file rather than a second pipe: reading two pipes
    sequentially deadlocks the moment the unread one fills its buffer.

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
            proc = lsdtt_test_env.run(cmd, capture_output=True, text=True,
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

    import signal

    timed_out = False
    fd, result_name = tempfile.mkstemp(prefix="rusage_", suffix=".json")
    os.close(fd)
    result_path = pathlib.Path(result_name)
    launcher = [sys.executable, str(SCRIPT_DIR / "rusage_launcher.py"),
                str(result_path), "--"]

    try:
        with tempfile.TemporaryFile(mode="w+") as errf:
            t0 = time.perf_counter()
            # start_new_session: the launcher leads its own process group, so a
            # timeout can kill the launcher AND the command it is waiting on.
            # Killing the launcher alone would orphan the measured command and
            # leave it running for the rest of the campaign.
            proc = lsdtt_test_env.popen(launcher + list(cmd), stdout=subprocess.PIPE,
                                    stderr=errf, text=True, env=env,
                                    cwd=str(cwd) if cwd else None,
                                    start_new_session=True)

            timer = None
            if timeout is not None:
                def _kill():
                    nonlocal timed_out
                    timed_out = True
                    try:
                        os.killpg(proc.pid, signal.SIGKILL)
                    except OSError:
                        pass
                timer = threading.Timer(timeout, _kill)
                timer.start()

            try:
                stdout = proc.stdout.read() if proc.stdout else ""
                if proc.stdout:
                    proc.stdout.close()
                proc.wait()
            finally:
                if timer is not None:
                    timer.cancel()

            outer_wall = time.perf_counter() - t0
            errf.seek(0)
            stderr = errf.read()

        try:
            res = json.loads(result_path.read_text())
        except (OSError, ValueError):
            res = None
    finally:
        result_path.unlink(missing_ok=True)

    if res is None:
        # The launcher died before writing a result: killed on timeout, or it
        # could not start. No memory figure is better than a wrong one.
        return {"returncode": proc.returncode, "stdout": stdout,
                "stderr": stderr, "wall_seconds": outer_wall,
                "timed_out": timed_out, "peak_rss_mb": None,
                "cpu_seconds": None}

    return {
        "returncode": res["exit_code"],
        "stdout": stdout,
        "stderr": stderr,
        "wall_seconds": res["wall_seconds"],
        "timed_out": timed_out,
        # ru_maxrss is kibibytes on Linux, bytes on macOS. Linux is where the
        # server campaign runs; on macOS this reads ~1024x high, which is
        # obvious rather than subtle.
        "peak_rss_mb": res["peak_rss_kib"] / 1024.0,
        # user + system CPU of the measured command. cpu_seconds / wall well
        # above 1 means the command used more cores than its thread label says.
        "cpu_seconds": res["cpu_seconds"],
        "launcher_peak_rss_mb": res["launcher_peak_rss_kib"] / 1024.0,
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
# lsdtt-synthetic-surfaces can generate a diamond-square surface of any
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

    # A LARGER cached source is used in preference to a new download. Rungs are
    # cached by name (`alpine_2000.tif`), so on a machine that ran the old
    # 8000 px ladder they already exist, cropped from the 8000 source; fetching
    # a fresh 4000 source would change nothing but the download bill. Taking
    # the smallest adequate source keeps any NEW crop as close as possible to
    # what a fresh fetch would give.
    larger = []
    for p in FIXTURE_DIR.glob(f"{region}_src*.tif"):
        try:
            px = int(p.stem.rsplit("_src", 1)[1])
        except (IndexError, ValueError):
            continue
        if px >= source_px:
            larger.append((px, p))
    if larger:
        px, p = min(larger)
        print(f"  {region}: using cached {px}px source for a {source_px}px "
              f"ladder ({p.name})")
        return p
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

    proc = lsdtt_test_env.run([str(fetch_exe), str(param_dir), param_path.name],
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
        r = lsdtt_test_env.run(["stat", "-f", "-c", "%T", str(path)],
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
    docs/developer_notes/testing_and_benchmarking.md), so choosing it favours nobody.  It is
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
        # MUST be set explicitly -- OMP_NUM_THREADS alone is NOT enough.
        #
        # A CLI that declares n_threads passes it into the library, and
        # resolve_threads() (src/cpu/lsdtt_raster_surface_fitting_cpu.cpp)
        # honours any POSITIVE value verbatim, applying it as an explicit
        # `num_threads(n)` clause. That clause OVERRIDES the environment. Only
        # n_threads == 0 falls through to omp_get_max_threads(), which is what
        # actually reads OMP_NUM_THREADS.
        #
        # So for any CLI defaulting to n_threads = 1 -- lsdtt-surface-metrics
        # and lsdtt-spectral-analysis do -- setting OMP_NUM_THREADS in the
        # environment achieves NOTHING: every measurement runs single-threaded
        # and the scaling curve is flat by construction. That silently wasted
        # the whole thread sweep of the 2026-08-31 campaigns, which swept only
        # the `surface` suite, i.e. only lsdtt-surface-metrics.
        #
        # Writing it here is a no-op for CLIs that default to 0 (they resolve to
        # the same count via the environment) and is the fix for those that do
        # not. Set it for every run so no future CLI can reintroduce the trap.
        "n_threads": str(threads),
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

    timings, rss, computes, writes, cpus = [], [], [], [], []
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
        if r.get("cpu_seconds") is not None:
            cpus.append(r["cpu_seconds"])
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
        "cpu_seconds_process": statistics.median(cpus) if cpus else None,
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


def reference_extra_args(engine: str, case: dict,
                         fixture: pathlib.Path) -> list[str]:
    """Per-fixture arguments a reference needs to do the SAME work as lsdtt3.

    `whitebox_flat_increment_from: <lsdtt3 param>` converts lsdtt3's flat
    GRADIENT into Whitebox's flat_increment, which is metres PER CELL:
    increment = gradient * cell size. Handing Whitebox the raw gradient (or
    leaving its default) makes it impose a 30x smaller increment on a 30 m
    grid -- a different amount of flat-resolution work. Same conversion as
    run_flow_benchmark.py uses for its value comparison.
    """
    param = case.get("whitebox_flat_increment_from")
    if engine != "whitebox" or not param:
        return []
    gradient = float(case["lsdtt"]["params"][param])
    try:
        from osgeo import gdal
        gdal.UseExceptions()
        ds = gdal.Open(str(fixture))
        cell = abs(float(ds.GetGeoTransform()[1]))
    except Exception as exc:
        raise RuntimeError(f"cannot read cell size of {fixture} to match "
                           f"Whitebox flat_increment: {exc}") from exc
    return ["--flat-increment", repr(gradient * cell)]


def measure_reference(engine: str, routine: str, fixture: pathlib.Path,
                      threads: int, repeats: int, python_exe: str,
                      timeout: float | None = None,
                      extra_args: list[str] | None = None) -> dict[str, Any] | None:
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
    if extra_args:
        cmd += list(extra_args)

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
    info["cpu_seconds_process"] = r.get("cpu_seconds")
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


# Peak RSS a do-nothing command may report through rusage_launcher.py before
# the memory column is declared untrustworthy. The launcher is a bare Python
# interpreter (~10 MiB measured on WSL2); anything near this limit means the
# launcher has started importing things or the inheritance trap is back.
RSS_FLOOR_LIMIT_MB = 64.0

# A row whose peak is below this multiple of the floor may be reporting the
# floor rather than itself, so it is flagged `rss_near_floor`.
RSS_NEAR_FLOOR_FACTOR = 1.5


def calibrate_rss_floor() -> dict[str, Any]:
    """Measure the smallest peak RSS any measurement can report. TRIPWIRE.

    A Linux child inherits its parent's lifetime peak RSS (see
    rusage_launcher.py). So the memory column is only honest if the process
    each command is launched FROM has a small peak. This measures that floor
    directly, by launching `/bin/true` -- whose own footprint is ~1 MiB --
    through exactly the path every real measurement takes.

    It is called AFTER describe_fixture() has run, which is the driver's own
    largest allocation, so it measures the floor in the state the campaign is
    actually measured in rather than at start-up when it would look fine.

    Three numbers are recorded, because together they show the mechanism
    rather than just a verdict:

      floor_mb          /bin/true through the launcher. Should be ~10 MiB.
      direct_mb         /bin/true launched straight from this driver, i.e. the
                        old behaviour. Equals the driver's peak when the trap
                        is present on this kernel.
      driver_peak_mb    this driver's own lifetime peak.

    `ok` is False when floor_mb exceeds RSS_FLOOR_LIMIT_MB; write_summary()
    then says, above the memory table, that the column must not be quoted.
    """
    import resource

    true_exe = shutil.which("true") or "/bin/true"
    out: dict[str, Any] = {"command": true_exe,
                           "limit_mb": RSS_FLOOR_LIMIT_MB}

    if not hasattr(os, "wait4"):
        out.update({"available": False, "ok": False})
        return out

    floors = []
    for _ in range(3):
        r = run_with_rusage([true_exe])
        if r.get("peak_rss_mb") is not None:
            floors.append(r["peak_rss_mb"])

    # The old, broken route, kept only to show whether the trap is live here.
    proc = lsdtt_test_env.popen([true_exe])
    _, _, ru = os.wait4(proc.pid, 0)
    proc.returncode = 0

    out.update({
        "available": bool(floors),
        "floor_mb": max(floors) if floors else None,
        "direct_mb": ru.ru_maxrss / 1024.0,
        "driver_peak_mb":
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
    })
    out["ok"] = bool(floors) and out["floor_mb"] <= RSS_FLOOR_LIMIT_MB
    return out


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
    tool-handle bug; see whitebox_reference.py's `_make_call`, and `rss_floor`,
    the tripwire for the peak-RSS inheritance trap; see calibrate_rss_floor().
    """
    out: dict[str, Any] = {}
    # A dict, not a float, so the "fixed overhead" printers (which render
    # floats as milliseconds) leave it alone.
    out["rss_floor"] = calibrate_rss_floor()

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
        proc = lsdtt_test_env.run(cmd, capture_output=True, text=True, check=False)
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
    "peak_rss_mb", "rss_near_floor", "cpu_seconds_process",
    "mcells_per_second",
    "report_only", "serial_reference", "sequential_case",
    "overhead_dominated", "reused_from",
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
        "sequential_case": bool(case.get("sequential", False)),
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
        "cpu_seconds_process": (f"{info['cpu_seconds_process']:.9g}"
                                if info.get("cpu_seconds_process") is not None
                                else ""),
    })

    # A peak this close to the launcher floor may be the floor, not the
    # command. See calibrate_rss_floor().
    floor_mb = (overhead.get("rss_floor") or {}).get("floor_mb")
    if floor_mb and info.get("peak_rss_mb") is not None:
        row["rss_near_floor"] = bool(
            info["peak_rss_mb"] < RSS_NEAR_FLOOR_FACTOR * floor_mb)
    else:
        row["rss_near_floor"] = ""

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
    L.append("`n/a` means the operation is NOT REACHABLE THROUGH THE INTERFACE "
             "BENCHMARKED HERE -- a result, not a gap. `-` means no comparable "
             "row was produced (skipped, failed, or overhead-dominated). The "
             "two are not the same claim and are not interchangeable.")
    L.append("")
    L.append("**`n/a` is a statement about the binding, not about the "
             "project.** TopoToolbox is driven through its Python package, "
             "which is early (0.0.12); TopoToolbox 3 implements more in MATLAB "
             "than is bound there. Read each reason below literally.")
    L.append("")
    _ni = sorted({(c, e, why) for c, m in NOT_IMPLEMENTED.items()
                  for e, why in m.items()})
    if _ni:
        L.append("| case | engine | why `n/a` |")
        L.append("|---|---|---|")
        for c, e, why in _ni:
            L.append(f"| {c} | {e} | {why} |")
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
    seq_cases = sorted({r["case"] for r in rows
                        if str(r.get("sequential_case")) == "True"})
    if seq_cases:
        L.append("Depression filling, breaching and flow-graph construction "
                 "are inherently sequential -- ordered traversals where each "
                 "step depends on the previous one. They are measured once, "
                 "at one thread, and are **not** in this table: "
                 f"`{'`, `'.join(seq_cases)}`. A flat curve for them would be "
                 "the correct answer, not a finding, so the sweep is not run.")
        L.append("")
    _scaling_table(L, ok)

    L.append("")
    L.append("## Peak memory (MiB, whole child process)")
    L.append("")
    fl = overhead.get("rss_floor") or {}
    if fl.get("available"):
        L.append(f"Measurement floor: `/bin/true` launched the same way as every "
                 f"row reports **{fl['floor_mb']:.1f} MiB** "
                 f"(tripwire limit {fl['limit_mb']:.0f} MiB). Launched directly "
                 f"from the driver it reports {fl['direct_mb']:.1f} MiB, against a "
                 f"driver peak of {fl['driver_peak_mb']:.1f} MiB -- the gap is the "
                 "peak-RSS inheritance that `rusage_launcher.py` exists to remove. "
                 f"Cells marked `<` are below {RSS_NEAR_FLOOR_FACTOR:g}x the floor "
                 "and may be reporting the floor rather than the engine.")
        L.append("")
        if not fl.get("ok"):
            L.append(f"> **MEMORY TRIPWIRE FAILED: the floor is "
                     f"{fl['floor_mb']:.1f} MiB.** Every row reports at least "
                     "this, so the memory column below measures the launcher, "
                     "not the engines. **Do not quote it.** Check that "
                     "`rusage_launcher.py` has not grown imports and that "
                     "`run_with_rusage()` still launches through it.")
            L.append("")
    else:
        L.append("> **The memory floor could not be measured**, so nothing "
                 "below shows whether a row reports the engine or an inherited "
                 "peak. Treat the column as unverified.")
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
    L.append("### What this campaign does NOT cover")
    L.append("")
    L.append("Absent from every table above, and not because they failed -- "
             "they are measured by a different driver whose axes do not fit "
             "this one (nesting depth, and per-phase decomposition):")
    L.append("")
    L.append("| not measured here | where it lives |")
    L.append("|---|---|")
    L.append("| nested basin extraction, swept over nesting depth | "
             "`tests/scripts/run_svector_benchmark.py --cases basins` |")
    L.append("| chi, end-to-end and decomposed by phase | "
             "`... --cases chi` |")
    L.append("| chi VALUE agreement with TopoToolbox | "
             "`... --cases chi_value --check` |")
    L.append("| fill / breach / D8 VALUE agreement | "
             "`tests/scripts/run_flow_benchmark.py --check` |")
    L.append("")
    L.append("This campaign asserts no values at all. Correctness lives in "
             "those two drivers and in the CTest suite.")
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
             "lives in the CTest suite; see `docs/developer_notes/testing_and_benchmarking.md`.")
    path.write_text("\n".join(L) + "\n")


def _fixture_order(rows: list[dict]) -> list[str]:
    seen = {}
    for r in rows:
        try:
            seen[r["fixture"]] = int(r["n_cells"]) if r["n_cells"] else 0
        except ValueError:
            seen[r["fixture"]] = 0
    return [k for k, _ in sorted(seen.items(), key=lambda kv: kv[1])]


def _engine_columns(rows: list[dict]) -> list[str]:
    """Engine columns for a table: those with rows, PLUS any known to be unable.

    An engine that cannot perform a case produces no row, so deriving the column
    list from observed rows alone drops it from the table entirely -- and the
    `n/a` that says "this package does not implement this" never appears. The
    absence then reads as though the engine were never considered, which is the
    opposite of the finding. See NOT_IMPLEMENTED.
    """
    engines = {r["engine"] for r in rows}
    cases = {r["case"] for r in rows}
    for case_name in cases:
        engines.update(NOT_IMPLEMENTED.get(case_name, {}))
    return sorted(engines)


def _ratio_table(L: list[str], rows: list[dict]) -> None:
    """Throughput per engine at the highest common thread count."""
    engines = _engine_columns(rows)
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
                    elif eng in NOT_IMPLEMENTED.get(case_name, {}):
                        # The package cannot do this. That is a RESULT, and must
                        # not look like a measurement we failed to take.
                        cells.append("n/a")
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
    # Sequential cases are excluded: they were measured at one thread count by
    # design, so a row for them here would either be empty or -- worse -- look
    # like a flat scaling curve that had been tested and failed.
    par = [r for r in rows if r["engine"] in ("lsdtt3", "whitebox")
           and str(r.get("sequential_case")) != "True"]
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
    engines = _engine_columns(mem)
    L.append("| case | fixture | " + " | ".join(engines) + " |")
    L.append("|---|---|" + "---|" * len(engines))
    for case_name in sorted({r["case"] for r in mem}):
        for fx in _fixture_order(mem):
            cells, any_val = [], False
            for eng in engines:
                # Peak across thread counts: memory is reported per engine and
                # fixture, and the worst case is the number that matters for
                # "will this fit".
                match = [r for r in mem
                         if r["case"] == case_name and r["fixture"] == fx
                         and r["engine"] == eng]
                if match:
                    worst = max(match, key=lambda r: float(r["peak_rss_mb"]))
                    near = str(worst.get("rss_near_floor")) == "True"
                    cells.append(("<" if near else "")
                                 + f"{float(worst['peak_rss_mb']):.0f}")
                    any_val = True
                elif eng in NOT_IMPLEMENTED.get(case_name, {}):
                    cells.append("n/a")
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
                         "live. Defaults to <out-dir>/benchmark_fixtures. "
                         "Because --out-dir itself defaults to "
                         "'benchmark_campaigns', fixtures always land beside "
                         "the campaign and never inside tests/data/. "
                         "The server profile stores ~300 MB "
                         "here, so on a cluster keep it off a quota'd home: "
                         "point it at project or scratch storage. Reused "
                         "across campaigns, so it is worth keeping.")
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
    elif args.out_dir:
        # DEFAULT WITH --out-dir: keep the DEMs beside the campaign.
        #
        # It used to default into the repository (tests/data/dem_cache) even
        # when the caller had said where their data should go, while the sbatch
        # wrapper defaulted to <submit dir>/benchmark_fixtures. So the two
        # routes disagreed about where ~512 MB of DEMs belonged, and a direct
        # run quietly filled the checkout -- on a cluster, usually a quota'd
        # home filesystem.
        #
        # Still a CACHE, and still reused: --out-dir is the PARENT of the
        # per-run campaign directory, not the campaign directory itself, so
        # successive runs into the same --out-dir share one set of fixtures and
        # download nothing after the first.
        FIXTURE_DIR = (pathlib.Path(args.out_dir).expanduser().resolve()
                       / "benchmark_fixtures")

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
            seq = bool(case.get("sequential"))
            per_engine = []
            for e in engines:
                src = case.get("reuse_from", {}).get(e)
                if src and src in selected_cases:
                    per_engine.append(f"{e}<-{src}")
                    continue
                t = _engine_thread_count(e, case, profile["threads"])
                per_engine.append(f"{e}x{t}")
                n += t * n_fixtures
            tag = "  [sequential: no thread sweep]" if seq else ""
            print(f"  {case['suite']:11s} {name:22s} "
                  f"{' '.join(per_engine)}{tag}")
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
    # DEMs that were never going to be next to their results. Both routes now
    # put the cache beside the campaign when a destination is given, but the
    # paths are still printed together so the log shows where they went rather
    # than it being something you have to know.
    print(f"campaign: {run_dir}")
    print(f"fixtures: {FIXTURE_DIR}"
          + ("" if (args.fixture_dir or args.out_dir)
             else "   (default: inside the repo -- "
                  "pass --out-dir or --fixture-dir to move it)"))
    print(f"profile:  {args.profile} ({profile['description']})")

    # -------------------------------------------------------------------------
    # Report which references this interpreter can actually import, BEFORE
    # measuring anything.
    #
    # The suite treats a missing reference as a clean skip, so without this the
    # only symptom is absent rows in a campaign that exited successfully -- and
    # you find out after the run, not before it. The Slurm wrapper has always
    # printed this; a direct `python ...` invocation did not, which is exactly
    # the route someone takes when they have just run `pixi install -e
    # benchmark` and assumed that was enough. It is not: `pixi install` creates
    # the environment, `pixi run -e benchmark python ...` is what uses it.
    # -------------------------------------------------------------------------
    # Print the interpreter PATH, not a derived nickname: which python is
    # running is the whole question here, and the path answers it unambiguously.
    print(f"\ninterpreter: {args.python_exe}")
    print("reference implementations importable by it:")
    _available = []
    for _mod in ("richdem", "topotoolbox", "whitebox_workflows"):
        _probe = lsdtt_test_env.run(
            [args.python_exe, "-c", f"import {_mod}"],
            capture_output=True, text=True)
        _ok = _probe.returncode == 0
        _available.append(_ok)
        print(f"    {_mod:20s} {'yes' if _ok else 'NO'}")
    if not any(_available):
        print("\n*** NO REFERENCE IMPLEMENTATION IS IMPORTABLE. ***\n"
              "    This campaign will contain lsdtt3 rows ONLY, and will still\n"
              "    exit successfully -- there will be nothing to compare against.\n"
              "    Run it through pixi rather than a bare python:\n\n"
              "        pixi run -e benchmark python "
              "tests/scripts/run_benchmark_suite.py ...\n\n"
              "    `pixi install -e benchmark` only CREATES that environment.\n",
              file=sys.stderr)
    elif not all(_available):
        print("    (a missing reference is SKIPPED, not an error -- it produces\n"
              "     no rows, so check this list before trusting the tables)")
    print()

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
    _fl = overhead.get("rss_floor") or {}
    if _fl.get("available"):
        print(f"  rss floor: {_fl['floor_mb']:.1f} MiB through the launcher "
              f"(direct launch {_fl['direct_mb']:.1f} MiB, driver peak "
              f"{_fl['driver_peak_mb']:.1f} MiB)"
              + ("" if _fl.get("ok") else
                 "   *** TRIPWIRE: memory column will not be trustworthy ***"))

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
    def _case_measurements(c):
        return sum(0 if c.get("reuse_from", {}).get(e) in selected_cases
                   else _engine_thread_count(e, c, profile["threads"])
                   for e in ["lsdtt3"] + list(c["engines"]))

    total = sum(len(fixtures) * _case_measurements(c)
                for c in selected_cases.values())
    done = 0

    def _reused_row(case_name, case, engine, fx_key, threads):
        """Copy of the source case's row for this engine, or None to measure.

        See REUSE_FROM_NOTE."""
        src = case.get("reuse_from", {}).get(engine)
        if not src or src not in selected_cases:
            return None
        for r in rows:
            if (r["case"] == src and r["engine"] == engine
                    and r["fixture"] == fx_key
                    and int(r["threads"]) == int(threads)):
                copy = dict(r)
                copy.update({
                    "suite": case["suite"], "case": case_name,
                    "report_only": bool(case.get("report_only", False)),
                    "sequential_case": bool(case.get("sequential", False)),
                    "reused_from": src,
                    "note": (f"REUSED from `{src}`: identical command, not "
                             f"re-measured. {case.get('note', '')}")[:400],
                })
                return copy
        return None

    # Cases that reuse another case's rows run AFTER everything else, so their
    # sources exist by the time they are reached.
    for case_name, case in sorted(selected_cases.items(),
                                  key=lambda kv: (bool(kv[1].get("reuse_from")),
                                                  kv[1]["suite"], kv[0])):
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

            # A sequential case has no parallelism in ANY engine, so the
            # sweep is skipped entirely and it is measured once. See
            # SEQUENTIAL_CASES_NOTE.
            case_threads = ([baseline_threads] if case.get("sequential")
                            else profile["threads"])

            for threads in case_threads:
                # lsdtt3
                reused = _reused_row(case_name, case, "lsdtt3", fx_key, threads)
                if reused is not None:
                    rows.append(reused)
                else:
                    info = measure_lsdtt(case_name, case, in_path, work_dir,
                                         threads, profile["repeats"], build_dir,
                                         timeout=measurement_timeout)
                    rows.append(build_row(case["suite"], case_name, "lsdtt3",
                                          fx_key, meta["sha256_12"],
                                          meta["n_cells"], threads, info, case,
                                          overhead))
                    done += 1
                    _progress(done, total, case_name, "lsdtt3", fx_key, threads,
                              rows[-1],
                              report_only=bool(case.get("report_only")))

                for engine, routine in sorted(case["engines"].items()):
                    # A capped reference (see REFERENCE_MAX_THREADS) is not
                    # measured above its cap; the scaling table shows `-`. The
                    # baseline is always measured, or `--threads 4,8` would
                    # drop the engine from the campaign without a word.
                    if (threads != baseline_threads
                            and threads > REFERENCE_MAX_THREADS.get(engine,
                                                                    threads)):
                        continue
                    reused = _reused_row(case_name, case, engine, fx_key,
                                         threads)
                    if reused is not None:
                        rows.append(reused)
                        continue
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
                    info = measure_reference(
                        engine, routine, in_path, threads, profile["repeats"],
                        args.python_exe, timeout=measurement_timeout,
                        extra_args=reference_extra_args(engine, case, in_path))
                    rows.append(build_row(case["suite"], case_name, engine,
                                          fx_key, meta["sha256_12"],
                                          meta["n_cells"], threads, info, case,
                                          overhead))
                    done += 1
                    _progress(done, total, case_name, engine, fx_key, threads,
                              rows[-1], routine=routine,
                              report_only=bool(case.get("report_only")))

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


def _progress(done, total, case_name, engine, fixture, threads, row,
              routine: str = "", report_only: bool = False) -> None:
    """One line per measurement.

    NAME WHAT THE ENGINE ACTUALLY RAN, NOT JUST THE CASE.

    Case names are lsdtt3-centric, and for a report-only case the reference is
    often running something quite different. Printing

        [   2/ 488] fill_spill_merge  whitebox  alpine_1000  1t  ok  16.97 ms

    says Whitebox did fill-spill-merge in 17 ms. It did not: it ran `sink`,
    which only labels cells lying in depressions, while lsdtt3 built a nested
    depression hierarchy AND ran Fill-Spill-Merge AND wrote a filled DEM. Next
    to lsdtt3's 201 ms that reads as lsdtt3 being 12x slower at the same job,
    and a reader watching the log has no way to know otherwise -- the
    report-only flag only ever appeared in the rendered tables, long after.

    So the routine is shown whenever it differs from the case name, and
    report-only rows say so inline.
    """
    status = row["status"]
    detail = ""
    if status == "ok":
        detail = f"{float(row['seconds_median']) * 1000:8.2f} ms"
        if row["overhead_dominated"]:
            detail += "  (overhead-dominated)"
    elif row["note"]:
        detail = row["note"][:60]

    label = engine
    if routine and routine != case_name:
        label = f"{engine}:{routine}"
    if report_only:
        detail = (detail + "  [report-only: not the same operation]").strip()

    print(f"  [{done:4d}/{total:4d}] {case_name:22s} {label:26s} "
          f"{fixture:16s} {threads}t  {status:12s} {detail}")


if __name__ == "__main__":
    sys.exit(main())
