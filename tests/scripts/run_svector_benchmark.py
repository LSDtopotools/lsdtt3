#!/usr/bin/env python3
"""What the s-vector BUYS: nested basin extraction, and chi.

docs/testing_and_benchmarking.md is the durable record; the working plan is the
(uncommitted) benchmarking plan, section 3b.

WHY THIS IS A SEPARATE DRIVER

  `run_benchmark_suite.py` sweeps threads x DEM size x region. Neither case here
  fits that model: the basin case sweeps NESTING DEPTH, and the chi case reports
  PER-PHASE timings that the campaign's row schema does not carry. Forcing them
  in would distort the campaign for two cases; this follows the precedent of
  run_flow_benchmark.py and run_hydraulic_dagger_benchmark.py instead.

WHAT IT MEASURES, AND WHY THESE TWO CASES

  lsdtt3's D8 build produces an s-vector: each cell's upslope set as a
  CONTIGUOUS INTERVAL. It costs roughly 57% of the D8 build and ~24 of the ~43
  bytes per cell of the flowinfo object. `run_flow_benchmark.py` and the
  campaign measure that COST. A cost with no counterweight is half a result.

  These are the two places the structure pays back, and both have a reference
  implementation to be scored against.

  1. NESTED BASINS. A nested basin is a SUB-INTERVAL of its parent's interval
     (lsdtt_basin_cpu.cpp:196-220), so nesting is free. The references cannot
     express it: Whitebox's `unnest_basins` emits ONE RASTER PER NESTING LEVEL,
     and TopoToolbox needs one `dependencemap` pass per level.

     *** THE RESULT IS THE SLOPE AGAINST DEPTH, NOT THE ABSOLUTE TIME. ***
     Every engine here rebuilds its flow structure per measurement -- lsdtt3
     because the CLI does, the references because they are made to, so that all
     three include the same work. That fixed cost dominates a small fixture and
     is NOT the quantity of interest. What separates the designs is how the
     time GROWS with depth: flat for lsdtt3 (a sub-interval read per basin),
     rising for the references (a pass per level). Fit and report the gradient;
     an intercept comparison on one depth says nothing.

     *** DO NOT REPLACE THIS WITH "EXTRACT N DISJOINT BASINS". ***  All three
     tools already do N disjoint basins in one whole-raster pass
     (`drainagebasins`, `basins`, `watershed`), so that case measures O(cells)
     everywhere, shows nothing, and invites the reader to conclude the s-vector
     buys nothing. DEPTH is the axis the references are structurally forced to
     pay for repeatedly. Sweeping basin COUNT instead would mix in basin size
     and blur the one axis that separates the designs.

  2. CHI. Chi extraction hard-requires the s-vector
     (lsdtt_chi_extraction_cpu.cpp:64) and uses it twice: `upslope_nodes()` for
     the node set, and `s_vector_size(flat) * pixel_area` for DRAINAGE AREA AT
     EVERY NODE IN O(1) -- and chi is an integral of A^(-m/n), so it needs A
     everywhere. TopoToolbox's `chitransform` takes a PRECOMPUTED area raster,
     which the caller must produce with `flow_accumulation()`.

     So the comparison is END-TO-END, DEM to chi, and decomposed. Timing
     `chitransform` alone against lsdtt3's chi loop alone would flatter
     TopoToolbox by excluding the accumulation lsdtt3 has already paid for.

     Whitebox and RichDEM have NO chi tool. This case is two-way.

A_0 MUST BE 1 m^2 ON BOTH SIDES

  Not a convention, a definition: the gradient of a chi-elevation profile equals
  the channel steepness index k_sn ONLY when A_0 = 1; otherwise it is k_sn
  scaled by A_0^(m/n). lsdtt3 already defaults to 1.0 everywhere. TopoToolbox's
  `chitransform` defaults to a0 = 1e6 and is passed 1.0 explicitly here. At
  m/n = 0.45 the difference is a factor of ~631, and it is SILENT -- the profile
  is unchanged and only the axis numbers move.

  `correctcellsize` is pinned too: it multiplies the supplied area by
  cellsize^2, so area in m^2 requires correctcellsize=False.

*** TIMINGS FROM A THERMALLY LIMITED LAPTOP ARE NOT RESULTS. ***

  This driver exists so the cases can be BUILT and verified anywhere. Publishable
  numbers come from a server, on an idle machine. The output records the host so
  a campaign cannot be misattributed later.

Usage:
    pixi run -e benchmark python tests/scripts/run_svector_benchmark.py \\
        --build-dir build --fixture alpine_small --depths 1,2,3,4,5
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import pathlib
import platform
import re
import shutil
import subprocess
import sys
import time

SKIP_EXIT = 77
REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "dem_cache"

# Sentinel used by lsdtt3's junction network for "this junction is an outlet".
K_OUTLET_JUNCTION = -2


def die_skip(msg: str):
    print(f"SKIP: {msg}")
    sys.exit(SKIP_EXIT)


def run_cli(exe: pathlib.Path, work: pathlib.Path, param_name: str,
            params: dict) -> tuple[str, float]:
    """Run an lsdtt3 CLI from a generated parameter file.

    Returns (stdout, wall_seconds). Wall time is recorded but is NOT the
    measurement: the CLIs print their own "Compute time:" lines, which exclude
    process startup and GDAL I/O, and those are what the rows use. On a small
    raster the process is overwhelmingly startup -- lsdtt3's fill has measured
    7.6 ms of compute inside a 279 ms process -- so a whole-process number would
    compare interpreter and loader behaviour rather than algorithms.
    """
    path = work / param_name
    path.write_text("".join(f"{k}: {v}\n" for k, v in params.items()))
    t0 = time.perf_counter()
    proc = subprocess.run([str(exe), str(work), param_name],
                          capture_output=True, text=True)
    wall = time.perf_counter() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"{exe.name} failed:\n{proc.stdout[-2000:]}\n"
                           f"{proc.stderr[-2000:]}")
    return proc.stdout, wall


def compute_times(stdout: str) -> list[float]:
    """Every "Compute time: X s" the CLI printed, in order."""
    return [float(m) for m in
            re.findall(r"Compute time:\s*([0-9.eE+-]+)\s*s", stdout)]


# -----------------------------------------------------------------------------
# Setup: fill, junctions, and the nesting depth of every junction
# -----------------------------------------------------------------------------
def prepare_filled(build_dir: pathlib.Path, fixture: pathlib.Path,
                   work: pathlib.Path) -> pathlib.Path:
    out = work / "filled"
    out.mkdir(parents=True, exist_ok=True)
    filled = out / f"{fixture.stem}_fill.tif"
    if filled.exists():
        return filled
    run_cli(build_dir / "lsdtt-raster-preprocessing", out, "fill.param", {
        "read_path": fixture.parent, "read_fname": fixture.name,
        "write_path": out, "write_fname": fixture.stem,
        "write_fill_raster": "true", "fill_min_slope": "0.0001",
        "write_output_format": "gtiff",
    })
    if not filled.exists():
        raise RuntimeError(f"expected {filled}")
    return filled


def junction_table(build_dir: pathlib.Path, filled: pathlib.Path,
                   work: pathlib.Path, threshold: int) -> list[dict]:
    """Junctions with their receiver, from lsdtt-channel-extraction."""
    out = work / "junctions"
    out.mkdir(parents=True, exist_ok=True)
    csv_path = out / "jn_channel_network_junctions.csv"
    if not csv_path.exists():
        run_cli(build_dir / "lsdtt-channel-extraction", out, "jn.param", {
            "read_path": filled.parent, "read_fname": filled.name,
            "write_path": out, "write_fname": "jn",
            "raster_is_filled": "true",
            "threshold_contributing_pixels": threshold,
            "write_channel_network_junctions": "true",
            "channel_network_junctions_output_format": "csv",
        })
    if not csv_path.exists():
        raise RuntimeError(f"expected {csv_path}")
    with csv_path.open() as fh:
        return list(csv.DictReader(fh))


def nesting_depth(junctions: list[dict]) -> dict[int, int]:
    """Depth of each junction below its outlet: an outlet is depth 1.

    Depth is measured DOWNSTREAM-outward, so depth d is "every junction within
    d-1 links of an outlet". That makes the depth-d seed set a strict superset
    of depth d-1, which is what makes the sweep monotone and the reference
    engines' per-level cost visible.
    """
    receiver = {int(j["junction_number"]): int(j["receiver_junction"])
                for j in junctions}
    depth: dict[int, int] = {}

    def resolve(j: int, guard: int = 0) -> int:
        if j in depth:
            return depth[j]
        r = receiver.get(j, K_OUTLET_JUNCTION)
        if r == K_OUTLET_JUNCTION or r not in receiver or guard > len(receiver):
            depth[j] = 1
        else:
            depth[j] = resolve(r, guard + 1) + 1
        return depth[j]

    sys.setrecursionlimit(max(10000, 4 * len(receiver) + 1000))
    for j in receiver:
        resolve(j)
    return depth


def write_seed_csv(junctions: list[dict], depth: dict[int, int],
                   max_depth: int, path: pathlib.Path) -> int:
    """Seed points = every junction at depth <= max_depth. Returns the count."""
    rows = [j for j in junctions
            if depth.get(int(j["junction_number"]), 10**9) <= max_depth]
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["latitude", "longitude"])
        for j in rows:
            w.writerow([j["latitude"], j["longitude"]])
    return len(rows)


# -----------------------------------------------------------------------------
# Case 1: nested basins
# -----------------------------------------------------------------------------
def basins_lsdtt3(build_dir: pathlib.Path, filled: pathlib.Path,
                  seeds: pathlib.Path, work: pathlib.Path, depth: int,
                  threshold: int, repeats: int) -> dict:
    out = work / f"basins_lsdtt3_d{depth}"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    # lsdtt-basin-extraction reports three phases, in this order. They are
    # NOT equally comparable, which is the whole point of splitting them:
    #
    #   flow_routing      the references build a flow structure too -> COMPARABLE
    #   junction_network  NOTHING else here builds one              -> EXCLUDED
    #   basin_selection   dependencemap / watershed equivalent      -> COMPARABLE
    #
    # `seconds` is flow + basin_selection. The junction network is reported
    # separately rather than folded in, because charging lsdtt3 for a structure
    # no reference builds is the same error `accum_d8` made with the s-vector.
    PHASES = ["flow_routing", "junction_network", "basin_selection"]

    totals, comparable, phase_rows = [], [], []
    scope = "cli_reported_compute"
    for _ in range(max(1, repeats)):
        stdout, wall = run_cli(build_dir / "lsdtt-basin-extraction", out,
                               "b.param", {
            "read_path": filled.parent, "read_fname": filled.name,
            "write_path": out, "write_fname": f"b{depth}",
            "raster_is_filled": "true",
            "threshold_contributing_pixels": threshold,
            "basin_outlet_file": seeds,
            "select_basins_from_outlets": "true",
            "write_basin_index_raster": "true",
            "write_output_format": "gtiff",
        })
        ct = compute_times(stdout)
        if not ct:
            # Say so rather than mislabelling. An earlier version of this
            # function fell back to wall clock and still reported the scope as
            # `cli_reported_compute`, which put whole-process rows -- 99 % of
            # which were process startup and GDAL I/O -- into a table that
            # claimed to exclude them.
            scope = "whole_process"
            totals.append(wall)
            comparable.append(wall)
            continue
        phase_rows.append(ct)
        totals.append(sum(ct))
        # Positional: the phases print in the fixed order above. Guard the
        # length so a future CLI change that adds a phase degrades to the total
        # rather than silently attributing the wrong number.
        if len(ct) == len(PHASES):
            comparable.append(ct[0] + ct[2])
        else:
            comparable.append(sum(ct))
            scope = "cli_reported_compute_unattributed"

    best = phase_rows[0] if phase_rows else []
    return {"engine": "lsdtt3",
            "seconds": min(comparable),
            "seconds_all": comparable,
            "seconds_all_phases": min(totals),
            "phase_names": PHASES if len(best) == len(PHASES) else [],
            "phase_seconds": best,
            "timing_scope": scope}


def basins_reference(engine: str, filled: pathlib.Path, seeds: pathlib.Path,
                     depth: int, repeats: int, python_exe: str,
                     timeout: float) -> dict:
    script = {"topotoolbox": "topotoolbox_reference.py",
              "whitebox": "whitebox_reference.py"}[engine]
    proc = subprocess.run(
        [python_exe, str(SCRIPT_DIR / script), "--routine", "basins_from_seeds",
         "--in", str(filled), "--seeds", str(seeds),
         "--levels", str(depth), "--repeats", str(repeats)],
        capture_output=True, text=True, timeout=timeout)
    if proc.returncode == SKIP_EXIT:
        return {"engine": engine, "error": "not available"}
    if proc.returncode != 0:
        return {"engine": engine, "error": proc.stderr[-400:].strip()}
    info = json.loads(proc.stdout.strip().splitlines()[-1])
    info["engine"] = engine
    return info


# -----------------------------------------------------------------------------
# Case 2: chi, end-to-end and decomposed
# -----------------------------------------------------------------------------
def chi_lsdtt3(build_dir: pathlib.Path, fixture: pathlib.Path,
               work: pathlib.Path, threshold: int, m_over_n: float,
               a0: float, repeats: int) -> dict:
    """End-to-end DEM -> chi, with the CLI's own per-phase compute times.

    Run from the RAW fixture, not the filled one: the end-to-end number has to
    include depression handling, because TopoToolbox's does.
    """
    exe = build_dir / "lsdtt-chi-analysis"
    if not exe.exists():
        return {"engine": "lsdtt3", "error": "lsdtt-chi-analysis not built "
                                             "(needs Eigen3)"}
    out = work / "chi_lsdtt3"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    totals, phases = [], []
    for _ in range(max(1, repeats)):
        try:
            stdout, wall = run_cli(exe, out, "c.param", {
                "read_path": fixture.parent, "read_fname": fixture.name,
                "write_path": out, "write_fname": "chi",
                "threshold_contributing_pixels": threshold,
                "chi_A0": a0, "m_over_n": m_over_n,
                "write_output_format": "gtiff",
            })
        except RuntimeError as exc:
            return {"engine": "lsdtt3", "error": str(exc)[-400:]}
        ct = compute_times(stdout)
        totals.append(sum(ct) if ct else wall)
        phases.append(ct)
    return {"engine": "lsdtt3", "seconds": min(totals), "seconds_all": totals,
            "phase_seconds": phases[0],
            "timing_scope": "cli_reported_compute"}


def chi_reference(filled_or_raw: pathlib.Path, threshold: int,
                  m_over_n: float, a0: float, repeats: int,
                  python_exe: str, timeout: float) -> dict:
    proc = subprocess.run(
        [python_exe, str(SCRIPT_DIR / "topotoolbox_reference.py"),
         "--routine", "chi_network", "--in", str(filled_or_raw),
         "--stream-threshold", str(threshold), "--mn", str(m_over_n),
         "--a0", str(a0), "--repeats", str(repeats)],
        capture_output=True, text=True, timeout=timeout)
    if proc.returncode == SKIP_EXIT:
        return {"engine": "topotoolbox", "error": "not available"}
    if proc.returncode != 0:
        return {"engine": "topotoolbox", "error": proc.stderr[-400:].strip()}
    info = json.loads(proc.stdout.strip().splitlines()[-1])
    info["engine"] = "topotoolbox"
    return info


# -----------------------------------------------------------------------------
# The chi-value fixture: generated, not committed.
#
# The chi comparison needs a DEM containing at least one COMPLETE catchment --
# one whose entire upslope area is inside the valid data. lsdtt3 drops any basin
# that is not (see chi_value_compare), so on a cropped real DEM every chi output
# is empty and there is nothing to compare.
#
# A synthetic island supplies one, and `lsdtt-synthetic-networks` is
# DETERMINISTIC: the same `random_seed` produces a byte-identical DEM (verified
# -- two runs, same sha256). So the fixture is generated from the recipe below
# rather than committed, which keeps a 640 KB raster out of the tree while
# leaving the comparison exactly reproducible.
#
# THE RECIPE IS PART OF THE PIN. Change any of it and the tolerances below mean
# nothing, because they were measured against this raster and no other. That is
# why the sha256 is checked: applying a pinned tolerance to a different DEM is
# how a benchmark comes to assert nothing at all.
#
# It has to be a FAT, SMOOTH island. Measured while building this:
# island_fill_fraction 0.55 with the default harmonics gave 16 % land in a
# lobate shape whose basins all still reached the coast -- zero complete basins.
# 0.95 with 3 harmonics gives ~57 % land and 150 px of interior, and 9 basins
# survive.
# -----------------------------------------------------------------------------
ISLAND_RECIPE = {
    "n_rows": 400, "n_cols": 400, "cell_size": 30.0,
    "domain_type": "efd_island",
    "island_fill_fraction": 0.95,
    "efd_n_harmonics": 3,
    "efd_roughness": 3.0,
    "outlet_mode": "auto",
    "steady_state_mode": "relief",
    "target_relief": 900.0,
    "m_exponent": 0.45,
    "n_exponent": 1.0,
    "random_seed": 1,
    "write_output_format": "gtiff",
}

# sha256 (first 12) of the DEM the tolerances were pinned against.
EXPECTED_ISLAND_SHA = "2834b2b08c64"

# Chi-value tolerances, pinned just above OBSERVED divergence -- the same rule
# run_flow_benchmark.py uses. A tolerance set to a round number far above
# reality asserts nothing.
#
# These are TIGHT on purpose. Since lsdtt3 adopted the trapezoidal rule the two
# implementations compute the identical quantity, so the expected slope is 1
# exactly and the only spread left is float32 rounding in lsdtt3's CSV output.
# MEASURED 2026-09-08 from THIS DRIVER's own run (threshold 100, no explicit
# outlets), 5 basins / 783 nodes:
#
#     slope                0.999999666 - 1.000000146
#     R^2                  1.0000000000 (all basins)
#     max abs residual     5.428e-06 chi units
#
# The residual is float32 rounding in lsdtt3's CSV, not a difference in value.
#
# The basin and node counts are pinned as FLOORS because they were measured
# here, from this configuration. An earlier hand-run with explicit outlets and
# n_largest_basins gave 9 basins / 804 nodes; pinning those numbers made the
# check fail against the driver's own output. Pin what the thing being checked
# actually produces, not what a differently-configured run produced.
#
# If this ever fails on slope, the first thing to check is the QUADRATURE RULE
# on both sides, not the terrain: a left-endpoint rule against a trapezoidal one
# reads as slope ~1.004 here, and mismatched A0 reads as A0_ratio^(m/n) (~631
# for a0=1e6 vs 1.0 at m/n=0.45).
CHI_VALUE_TOLERANCES = {
    # Just above observed on both sides of 1.
    "slope_min":            0.999999,
    "slope_max":            1.000001,
    "r_squared_min":        0.999999999,
    "max_abs_residual_chi": 1.0e-5,
    # Floors, not equalities: a fixture or filter change that REDUCES the
    # comparison is the failure worth catching. More basins is not a problem.
    "min_basins":           5,
    "min_matched_nodes":    780,
}


def prepare_island(build_dir: pathlib.Path, work: pathlib.Path) -> pathlib.Path:
    """Generate the island fixture, and refuse to proceed on a different one."""
    out = work / "island"
    out.mkdir(parents=True, exist_ok=True)
    dem = out / "isl_ocn_dem.tif"
    if not dem.exists():
        exe = build_dir / "lsdtt-synthetic-networks"
        if not exe.exists():
            die_skip(f"{exe} not found; needed to generate the chi fixture")
        params = {"write_path": out, "write_fname": "isl"}
        params.update(ISLAND_RECIPE)
        run_cli(exe, out, "island.param", params)
    if not dem.exists():
        raise RuntimeError(f"expected {dem}")

    import hashlib
    digest = hashlib.sha256(dem.read_bytes()).hexdigest()[:12]
    if digest != EXPECTED_ISLAND_SHA:
        # A skip, not a failure: a different GDAL or compiler produces a
        # different raster, which is an environment difference rather than a
        # defect in lsdtt3. Failing here would turn someone else's toolchain
        # into an lsdtt3 bug report.
        die_skip(
            f"island fixture sha256 {digest} != pinned {EXPECTED_ISLAND_SHA}. "
            f"The chi-value tolerances were measured against the pinned raster "
            f"and mean nothing against another. Re-measure and re-pin both "
            f"together, or investigate why the generator changed.")
    return dem


# -----------------------------------------------------------------------------
# Case 3: chi VALUE comparison
# -----------------------------------------------------------------------------
def chi_value_compare(build_dir: pathlib.Path, filled: pathlib.Path,
                      work: pathlib.Path, threshold: int, m_over_n: float,
                      a0: float, python_exe: str, timeout: float) -> list[dict]:
    """Do lsdtt3 and TopoToolbox agree on chi, node for node?

    WHY A PER-BASIN REGRESSION AND NOT A DIFFERENCE.

    chi is defined relative to a BASELEVEL: it is 0 at the outlet the
    integration started from and accumulates upstream. The two tools do not
    necessarily start from the same node -- lsdtt3 seeds from its own basin
    outlets, TopoToolbox from the stream network's outlets -- so a raw
    difference would be dominated by a constant offset that says nothing about
    the integrand.

    What the comparison is actually testing is whether the two integrate the
    SAME QUANTITY: chi = integral of (A0/A)^(m/n) ds. So within each basin,
    regress lsdtt3's chi on TopoToolbox's:

        SLOPE     must be 1.0. It is the ratio of the integrands. A slope of
                  A0_ratio^(m/n) is the signature of mismatched A0 -- at
                  m/n = 0.45 an unnoticed a0=1e6 on one side shows up as ~631,
                  which is why this test is worth running at all.
        INTERCEPT is the baselevel offset. Free; not a defect.
        RESIDUAL  is the real disagreement, in chi units.

    Node join is on the flat index (row * n_cols + col), so no coordinate
    matching and no tolerance on position.

    *** THIS NEEDS A FIXTURE WITH A COMPLETE CATCHMENT. ***

    lsdtt3 drops any basin whose upslope area is `is_influenced_by_invalid`,
    and that test flags a cell if ANY of its 8 neighbours is nodata, or if it
    sits on an open boundary (lsdtt_flowinfo_cpu.cpp:690-728). On a cropped DEM
    every basin's headwaters reach an edge, so EVERY chi output is empty --
    raster, fit nodes, segments and ksn summary alike. That is a deliberate
    stance (a truncated basin has the wrong drainage area, so its chi is wrong)
    and TopoToolbox does not take it, which is itself worth reporting.

    A synthetic island works, but only a FAT one. Measured while building this:
    `island_fill_fraction: 0.55` with the default harmonics produced 16% land in
    a lobate shape whose basins all still reached the coast, and zero complete
    basins. The recipe that works:

        lsdtt-synthetic-networks
            n_rows/n_cols        400
            domain_type          efd_island
            island_fill_fraction 0.95      -> ~57% land, 150 px of interior
            efd_n_harmonics      3         -> smooth coast, few lobes
            efd_roughness        3.0
            steady_state_mode    relief
            target_relief        900

    then chi with `basin_outlet_snapping_method: none` -- the default
    `downstream` snapping drags supplied outlets back onto the trunk near the
    coast, which reintroduces the very problem the inland outlet avoids.

    That gave 9 complete basins and 804 chi nodes on a 400^2 island.
    """
    import math

    out = work / "chi_value"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    exe = build_dir / "lsdtt-chi-analysis"
    if not exe.exists():
        return [{"engine": "lsdtt3", "case": "chi_value",
                 "error": "lsdtt-chi-analysis not built (needs Eigen3)"}]

    # lsdtt3: the FIT NODE table, not the chi coordinate raster.
    #
    # The raster path runs through write_chi_data_maps, which drops any basin
    # whose upslope area touches invalid data -- and `is_influenced_by_invalid`
    # counts the raster edge. On a CROPPED fixture essentially every basin's
    # headwaters reach an edge, so that path yields nothing at all. See the
    # note recorded with this case: it is a deliberate correctness stance in
    # lsdtt3, not a defect, but it makes the raster unusable here.
    try:
        stdout, _ = run_cli(exe, out, "cv.param", {
            "read_path": filled.parent, "read_fname": filled.name,
            "write_path": out, "write_fname": "cv",
            "raster_is_filled": "true",
            "threshold_contributing_pixels": threshold,
            "chi_A0": a0, "m_over_n": m_over_n,
            "fit_chi_networks": "true", "write_fit_tables": "true",
            "write_output_format": "gtiff",
        })
    except RuntimeError as exc:
        return [{"engine": "lsdtt3", "case": "chi_value",
                 "error": str(exc)[-300:]}]

    nodes = sorted(out.glob("cv_chi_fit_nodes_mn_*.csv"))
    if not nodes:
        return [{"engine": "lsdtt3", "case": "chi_value",
                 "error": "no chi fit node table written"}]

    # TopoToolbox: chi as a full grid, same filled DEM, same a0 and m/n.
    ttb_tif = out / "ttb_chi.tif"
    proc = subprocess.run(
        [python_exe, str(SCRIPT_DIR / "topotoolbox_reference.py"),
         "--routine", "chi_grid", "--in", str(filled), "--out", str(ttb_tif),
         "--stream-threshold", str(threshold), "--mn", str(m_over_n),
         "--a0", str(a0)],
        capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        return [{"engine": "topotoolbox", "case": "chi_value",
                 "error": proc.stderr[-300:].strip()}]

    from osgeo import gdal
    gdal.UseExceptions()
    ds = gdal.Open(str(ttb_tif))
    ttb = ds.GetRasterBand(1).ReadAsArray().astype(float)
    n_cols = ds.RasterXSize
    ds = None

    per_basin: dict[int, list[tuple[float, float]]] = {}
    with nodes[0].open() as fh:
        for row in csv.DictReader(fh):
            flat = int(row["flat"])
            r, c = divmod(flat, n_cols)
            if not (0 <= r < ttb.shape[0] and 0 <= c < ttb.shape[1]):
                continue
            b = ttb[r, c]
            if not math.isfinite(b):
                continue
            per_basin.setdefault(int(row["basin_index"]), []).append(
                (b, float(row["chi"])))

    rows_out = []
    for basin, pairs in sorted(per_basin.items()):
        if len(pairs) < 3:
            continue
        xs = [p[0] for p in pairs]
        ys = [p[1] for p in pairs]
        n = len(xs)
        mx, my = sum(xs) / n, sum(ys) / n
        sxx = sum((x - mx) ** 2 for x in xs)
        sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        if sxx <= 0.0:
            continue
        slope = sxy / sxx
        intercept = my - slope * mx
        resid = [y - (slope * x + intercept) for x, y in zip(xs, ys)]
        sst = sum((y - my) ** 2 for y in ys)
        ssr = sum(r * r for r in resid)
        rows_out.append({
            "case": "chi_value", "engine": "lsdtt3_vs_topotoolbox",
            "basin_index": basin, "n_matched": n,
            "slope": slope, "intercept": intercept,
            "r_squared": (1.0 - ssr / sst) if sst > 0 else float("nan"),
            "max_abs_residual_chi": max(abs(r) for r in resid),
            "a0": a0, "m_over_n": m_over_n,
        })
    if not rows_out:
        rows_out = [{"case": "chi_value", "engine": "lsdtt3_vs_topotoolbox",
                     "error": "no basin had >=3 nodes present in both"}]
    return rows_out


# -----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-dir", default="build")
    ap.add_argument("--fixture", default="alpine_small")
    ap.add_argument("--depths", default="1,2,3,4,5",
                    help="nesting depths to sweep")
    ap.add_argument("--threshold", type=int, default=200,
                    help="contributing-pixel threshold defining the network")
    ap.add_argument("--m-over-n", type=float, default=0.45)
    ap.add_argument("--a0", type=float, default=1.0,
                    help="chi reference area (m2). MUST be 1.0 for the chi "
                         "gradient to be k_sn; see the module docstring")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--timeout-seconds", type=float, default=1800.0)
    ap.add_argument("--python-exe", default=sys.executable)
    ap.add_argument("--out-dir", default="")
    ap.add_argument("--cases", default="basins,chi,chi_value")
    ap.add_argument("--check", action="store_true",
                    help="enforce CHI_VALUE_TOLERANCES and exit non-zero on a "
                         "violation. Value metrics only -- nothing here "
                         "asserts a timing, which would be meaningless off a "
                         "dedicated machine.")
    args = ap.parse_args()

    build_dir = pathlib.Path(args.build_dir).resolve()
    fixture = FIXTURE_DIR / f"{args.fixture}.tif"
    if not fixture.exists():
        die_skip(f"fixture {fixture} not found; run the flow benchmark once "
                 f"to fetch it")
    for exe in ("lsdtt-raster-preprocessing", "lsdtt-channel-extraction",
                "lsdtt-basin-extraction"):
        if not (build_dir / exe).exists():
            die_skip(f"{exe} not found in {build_dir}")

    if args.a0 != 1.0:
        print(f"WARNING: --a0 {args.a0} != 1.0. The chi-elevation gradient is "
              f"then k_sn scaled by a0^(m/n), NOT k_sn. See --help.\n")

    stamp = dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    out_dir = pathlib.Path(args.out_dir).resolve() if args.out_dir else \
        (build_dir / "svector_benchmark" / f"{args.fixture}_{stamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    work = out_dir / "work"
    work.mkdir(exist_ok=True)

    print(f"fixture   : {fixture}")
    print(f"output    : {out_dir}")
    print(f"host      : {platform.node()}  ({platform.machine()})")
    print("NOTE: timings from a thermally limited machine are NOT results.\n")

    rows: list[dict] = []
    cases = {c.strip() for c in args.cases.split(",") if c.strip()}

    filled = prepare_filled(build_dir, fixture, work)

    # -- case 1: nested basins ---------------------------------------------
    if "basins" in cases:
        print("=== nested basins, swept over nesting depth ===")
        junctions = junction_table(build_dir, filled, work, args.threshold)
        depth = nesting_depth(junctions)
        max_available = max(depth.values()) if depth else 0
        print(f"  junction network: {len(junctions)} junctions, "
              f"max nesting depth {max_available}")

        for d in [int(x) for x in args.depths.split(",")]:
            if d > max_available:
                print(f"  depth {d}: beyond the network's {max_available}; "
                      f"skipped")
                continue
            seeds = work / f"seeds_d{d}.csv"
            n_seeds = write_seed_csv(junctions, depth, d, seeds)
            print(f"  depth {d}: {n_seeds} seed junctions")

            res = [basins_lsdtt3(build_dir, filled, seeds, work, d,
                                 args.threshold, args.repeats)]
            for eng in ("topotoolbox", "whitebox"):
                res.append(basins_reference(eng, filled, seeds, d,
                                            args.repeats, args.python_exe,
                                            args.timeout_seconds))
            for r in res:
                r.update({"case": "basins_nested", "depth": d,
                          "n_seeds": n_seeds, "fixture": args.fixture})
                # Recorded so the reader is reminded at the point of use that a
                # single depth's absolute time is not the measurement.
                r["result_is_slope_vs_depth"] = True
                rows.append(r)
                s = (f"{r['seconds']:.4f} s" if "seconds" in r
                     else f"ERROR {r.get('error','')[:60]}")
                print(f"    {r['engine']:14s} {s}")

    # -- case 2: chi --------------------------------------------------------
    if "chi" in cases:
        print("\n=== chi, end-to-end and decomposed ===")
        print(f"  A0 = {args.a0} m2, m/n = {args.m_over_n}  "
              f"(both sides; TopoToolbox's own default a0 is 1e6)")
        res = [chi_lsdtt3(build_dir, fixture, work, args.threshold,
                          args.m_over_n, args.a0, args.repeats),
               chi_reference(fixture, args.threshold, args.m_over_n,
                             args.a0, args.repeats, args.python_exe,
                             args.timeout_seconds)]
        for r in res:
            r.update({"case": "chi_network", "depth": "", "n_seeds": "",
                      "fixture": args.fixture})
            rows.append(r)
            s = (f"{r['seconds']:.4f} s" if "seconds" in r
                 else f"ERROR {r.get('error','')[:120]}")
            print(f"    {r['engine']:14s} {s}")
            for name, sec in zip(r.get("phase_names", []),
                                 r.get("phase_seconds", [])):
                print(f"        {name:24s} {sec:.4f} s")

    # -- case 3: chi values -------------------------------------------------
    if "chi_value" in cases:
        print("\n=== chi VALUE comparison (lsdtt3 vs TopoToolbox) ===")
        print("  per-basin regression of lsdtt3 chi on TopoToolbox chi:")
        print("  slope must be 1.0 (same integrand); intercept is the "
              "baselevel offset and is free")
        # Its OWN fixture, not the campaign's: a cropped DEM yields no complete
        # basin and therefore no chi at all. See prepare_island().
        island = prepare_island(build_dir, work)
        print(f"  fixture: generated island, sha256 {EXPECTED_ISLAND_SHA} "
              f"(deterministic, random_seed="
              f"{ISLAND_RECIPE['random_seed']})")
        vrows = chi_value_compare(build_dir, island, work, args.threshold,
                                  args.m_over_n, args.a0, args.python_exe,
                                  args.timeout_seconds)
        for r in vrows:
            r["fixture"] = args.fixture
            rows.append(r)
            if "error" in r:
                print(f"    ERROR {r['error'][:110]}")
            else:
                print(f"    basin {r['basin_index']:3d}  n={r['n_matched']:5d}  "
                      f"slope={r['slope']:.6f}  intercept={r['intercept']:+.4f}  "
                      f"R2={r['r_squared']:.6f}  "
                      f"max|resid|={r['max_abs_residual_chi']:.4f}")

    # -- tolerance enforcement (VALUE metrics only) --------------------------
    failures: list[str] = []
    if args.check:
        vals = [r for r in rows if r.get("case") == "chi_value"
                and "error" not in r]
        errs = [r for r in rows if r.get("case") == "chi_value"
                and "error" in r]
        for r in errs:
            failures.append(f"chi_value: {r['error']}")

        t = CHI_VALUE_TOLERANCES
        if len(vals) < t["min_basins"]:
            failures.append(
                f"chi_value: {len(vals)} basins compared, expected at least "
                f"{t['min_basins']}. Fewer basins means the fixture or the "
                f"completeness filter changed, and the per-basin tolerances "
                f"below were pinned on a different population.")
        matched = sum(r["n_matched"] for r in vals)
        if matched < t["min_matched_nodes"]:
            failures.append(
                f"chi_value: {matched} nodes matched, expected at least "
                f"{t['min_matched_nodes']}")
        for r in vals:
            b = r["basin_index"]
            if not (t["slope_min"] <= r["slope"] <= t["slope_max"]):
                failures.append(
                    f"chi_value basin {b}: slope {r['slope']:.6f} outside "
                    f"[{t['slope_min']}, {t['slope_max']}]. Check the "
                    f"QUADRATURE RULE and A0 on both sides before the terrain.")
            if r["r_squared"] < t["r_squared_min"]:
                failures.append(
                    f"chi_value basin {b}: R^2 {r['r_squared']:.8f} below "
                    f"{t['r_squared_min']}")
            if r["max_abs_residual_chi"] > t["max_abs_residual_chi"]:
                failures.append(
                    f"chi_value basin {b}: max abs residual "
                    f"{r['max_abs_residual_chi']:.3e} above "
                    f"{t['max_abs_residual_chi']:.1e} chi units")

    # -- output -------------------------------------------------------------
    fields = ["fixture", "case", "depth", "n_seeds", "engine", "seconds",
              "timing_scope", "result_is_slope_vs_depth", "phase_names",
              "phase_seconds", "basin_index", "n_matched", "slope",
              "intercept", "r_squared", "max_abs_residual_chi", "a0",
              "m_over_n", "error"]
    csv_path = out_dir / "svector_benchmark.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            r = dict(r)
            for k in ("phase_names", "phase_seconds"):
                if isinstance(r.get(k), list):
                    r[k] = ";".join(str(x) for x in r[k])
            w.writerow(r)

    (out_dir / "run_manifest.json").write_text(json.dumps({
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "host": platform.node(), "machine": platform.machine(),
        "platform": platform.platform(),
        "fixture": args.fixture, "depths": args.depths,
        "threshold_contributing_pixels": args.threshold,
        "m_over_n": args.m_over_n, "a0": args.a0, "repeats": args.repeats,
        "command": " ".join(sys.argv),
        "timings_are_results": False,
        "note": "Timings are only publishable from an idle server; see the "
                "module docstring.",
    }, indent=2))

    print(f"\nwrote {csv_path}")
    print("Timings here are NOT results unless this ran on an idle server.")

    if args.check:
        if failures:
            print("\nTOLERANCE FAILURES:")
            for f in failures:
                print(f"  - {f}")
            return 1
        print("\nchi value tolerances: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
