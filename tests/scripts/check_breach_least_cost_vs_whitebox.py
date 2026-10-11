#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Fuzz lsdtt3's least-cost breach against Whitebox 2.0.6, bit for bit.

WHAT IS COMPARED

  lsdtt-raster-preprocessing with breach_algorithm_choice least_cost,
  breach_lc_edge_policy whitebox and fill_remainder_after_breach false, against
  whitebox_workflows' breach_depressions_least_cost with fill_deps False, on
  the same float32 input. The requirement is float32(Whitebox's f64 output) ==
  lsdtt3's float32 output on every valid cell, with identical nodata masks.
  Any difference is a porting bug, not a tolerance to pin.

  The real-terrain rows in run_flow_benchmark.py (breach_least_cost) cover two
  DEMs under one parameter set. This covers the parameter space on many small
  DEMs, where tie-breaking, max_dist and max_cost edges, flats and nodata are
  dense:

    * N seeded random DEMs (default 200), 20x20 to 64x64 cells, cell size 1 or
      30 m. Values are integers or multiples of 0.125, so every elevation and
      every channel target pit_z - L * 0.125 is exact in float32 and the
      comparison cannot be blurred by rounding. Integer DEMs have many flats
      and equal-elevation pits, i.e. many heap ties.
    * About a third carry nodata blobs (interior holes and border notches),
      which act as outlets.
    * Every DEM runs every combination of minimize_dist {false, true},
      max_dist {3, 10, 100} and max_cost {unlimited, 2.5}: 12 runs per DEM.
      max_cost 2.5 is metres (minimize_dist false) or metres x cells (true)
      on both sides.

  A coverage guard fails the run if the fuzz exercised nothing: there must be
  runs that breached pits, runs with pits left unsolved, and DEMs with nodata.

  On any mismatch the DEM is kept in --work-dir with the parameters and the
  first differing cell, so it can be rerun by hand.

EXIT CODES

  0 pass, 1 mismatch (or coverage guard), 77 skip (whitebox_workflows, numpy
  or GDAL not importable, or --exe missing).

Run from the pixi benchmark environment:
  pixi run -e benchmark python -I tests/scripts/check_breach_least_cost_vs_whitebox.py \\
      --exe build/lsdtt-raster-preprocessing --work-dir <dir>
"""

from __future__ import annotations

import argparse
import concurrent.futures
import itertools
import json
import math
import os
import pathlib
import re
import shutil
import subprocess
import sys

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import lsdtt_test_env  # noqa: E402

SKIP_EXIT = 77
NODATA = -9999.0
FLAT_INCREMENT = 0.125
SMALL_MAX_COST = 2.5
MAX_DISTS = (3, 10, 100)
MAX_COSTS = (None, SMALL_MAX_COST)     # None = unlimited
MIN_DISTS = (False, True)

SUMMARY_RE = re.compile(
    r"(\d+) pits found, (\d+) skipped \(no longer pits\), (\d+) breached, "
    r"(\d+) unsolved")


def skip(msg: str) -> int:
    print(f"SKIP: {msg}", file=sys.stderr)
    return SKIP_EXIT


def make_dem(rng, np):
    """One random DEM: (values float32, cell size, has_nodata, kind)."""
    rows = int(rng.integers(20, 65))
    cols = int(rng.integers(20, 65))
    # Smoothed noise plus a tilt gives drainage with many local pits; the
    # roughness and relief vary so some DEMs are pitted everywhere and others
    # only here and there.
    noise = rng.normal(0.0, 1.0, (rows, cols))
    for _ in range(int(rng.integers(0, 4))):
        noise = (noise + np.roll(noise, 1, 0) + np.roll(noise, -1, 0)
                 + np.roll(noise, 1, 1) + np.roll(noise, -1, 1)) / 5.0
    relief = float(rng.choice([2.0, 5.0, 20.0, 80.0]))
    noise = noise / (np.abs(noise).max() or 1.0) * relief
    rr, cc = np.mgrid[0:rows, 0:cols]
    tilt = float(rng.uniform(-0.3, 0.3)) * rr + float(rng.uniform(-0.3, 0.3)) * cc
    # Base 200 m keeps every cell above min_elevation (0 m by default, which
    # would otherwise turn low cells into nodata on the lsdtt3 side only).
    z = 200.0 + noise + tilt
    kind = "integer" if rng.random() < 0.5 else "dyadic"
    if kind == "integer":
        z = np.round(z)
    else:
        z = np.round(z * 8.0) / 8.0
    z = z.astype(np.float32)

    has_nodata = bool(rng.random() < 0.35)
    if has_nodata:
        for _ in range(int(rng.integers(1, 4))):
            r0 = int(rng.integers(0, rows))
            c0 = int(rng.integers(0, cols))
            rad = float(rng.uniform(1.0, 4.0))
            z[(rr - r0) ** 2 + (cc - c0) ** 2 <= rad * rad] = NODATA
    cs = float(rng.choice([1.0, 30.0]))
    return z, cs, has_nodata, kind


def write_dem(path: pathlib.Path, z, cs: float) -> None:
    from osgeo import gdal, osr
    gdal.UseExceptions()
    rows, cols = z.shape
    # No DEFLATE predictor: Whitebox 2.0.6 reads a predictor-compressed
    # float raster as all NaN (see whitebox_reference.py).
    ds = gdal.GetDriverByName("GTiff").Create(str(path), cols, rows, 1,
                                              gdal.GDT_Float32)
    ds.SetGeoTransform((500000.0, cs, 0.0, 6000000.0, 0.0, -cs))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32630)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(NODATA)
    band.WriteArray(z)
    band.FlushCache()
    ds = None


def run_lsdtt(exe: pathlib.Path, dem: pathlib.Path, out_dir: pathlib.Path,
              md: bool, max_dist: int, max_cost):
    out_dir.mkdir(parents=True, exist_ok=True)
    params = {
        "read_path": str(dem.parent),
        "read_fname": dem.name,
        "write_path": str(out_dir),
        "write_fname": "lc",
        "enforce_projection": "false",
        "check_vertical_units_against_usgs": "false",
        "write_output_format": "gtiff",
        "write_breached_raster": "true",
        "breach_algorithm_choice": "least_cost",
        "breach_lc_edge_policy": "whitebox",
        "fill_remainder_after_breach": "false",
        "breach_lc_flat_increment": repr(FLAT_INCREMENT),
        "breach_lc_max_dist": str(max_dist),
        "breach_lc_minimize_dist": "true" if md else "false",
    }
    if max_cost is not None:
        params["breach_lc_max_cost"] = repr(max_cost)
    (out_dir / "lc.param").write_text("".join(f"{k}: {v}\n" for k, v in params.items()))
    env = lsdtt_test_env.system_env()
    env["OMP_NUM_THREADS"] = "1"
    proc = lsdtt_test_env.run([str(exe), str(out_dir), "lc.param"],
                              capture_output=True, text=True, env=env,
                              stdin=subprocess.DEVNULL)
    if proc.returncode != 0:
        raise RuntimeError(f"lsdtt3 failed (exit {proc.returncode}) on {dem}:\n"
                           f"{proc.stdout[-1500:]}\n{proc.stderr[-1500:]}")
    m = SUMMARY_RE.search(proc.stdout)
    counts = tuple(int(g) for g in m.groups()) if m else None
    return out_dir / "lc_breached.tif", counts


def read_float32(path: pathlib.Path):
    from osgeo import gdal
    gdal.UseExceptions()
    ds = gdal.Open(str(path))
    band = ds.GetRasterBand(1)
    a = band.ReadAsArray().astype("float32")
    nd = band.GetNoDataValue()
    ds = None
    return a, nd


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, help="path to lsdtt-raster-preprocessing")
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--n-dems", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261011)
    ap.add_argument("--jobs", type=int, default=4,
                    help="concurrent lsdtt3 CLI runs (each single-threaded)")
    args = ap.parse_args()

    try:
        import numpy as np
        from osgeo import gdal  # noqa: F401
        import whitebox_workflows as wb
    except ImportError as exc:
        return skip(f"{exc}. Use `pixi run -e benchmark ...`")

    exe = pathlib.Path(args.exe).resolve()
    if not exe.exists():
        return skip(f"{exe} not found; build lsdtt3 first")

    work = pathlib.Path(args.work_dir).resolve() / "breach_least_cost_fuzz"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    wbe = wb.WbEnvironment()
    wbe.verbose = False
    wbe.max_procs = 1

    rng = np.random.default_rng(args.seed)
    combos = list(itertools.product(MIN_DISTS, MAX_DISTS, MAX_COSTS))
    n_runs = n_fail = 0
    n_breached_runs = n_unsolved_runs = n_nodata_dems = n_changed_runs = 0
    failures = []

    for idx in range(args.n_dems):
        z, cs, has_nodata, kind = make_dem(rng, np)
        n_nodata_dems += has_nodata
        dem_dir = work / f"dem{idx:03d}"
        dem_dir.mkdir()
        dem_path = dem_dir / "dem.tif"
        write_dem(dem_path, z, cs)
        valid_in = z != np.float32(NODATA)

        # lsdtt3: all combinations for this DEM, concurrently.
        with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            futs = {
                combo: pool.submit(run_lsdtt, exe, dem_path,
                                   dem_dir / f"md{int(combo[0])}_d{combo[1]}_"
                                             f"c{'inf' if combo[2] is None else combo[2]}",
                                   *combo)
                for combo in combos
            }
            lsdtt_results = {c: f.result() for c, f in futs.items()}

        wbe.working_directory = str(dem_dir)
        wdem = wbe.read_raster("dem.tif")
        dem_failed = False
        for combo in combos:
            md, max_dist, max_cost = combo
            kw = dict(dem=wdem, max_dist=max_dist, minimize_dist=md,
                      flat_increment=FLAT_INCREMENT, fill_deps=False)
            if max_cost is not None:
                kw["max_cost"] = max_cost
            wout = np.asarray(wbe.hydrology.breach_depressions_least_cost(**kw).to_numpy(),
                              dtype=np.float64)
            wb32 = wout.astype(np.float32)
            wb_valid = wout != NODATA

            out_path, counts = lsdtt_results[combo]
            ls32, nd = read_float32(out_path)
            ls_valid = ls32 != np.float32(nd if nd is not None else NODATA)

            n_runs += 1
            mask_ok = np.array_equal(ls_valid, wb_valid) and np.array_equal(ls_valid, valid_in)
            both = ls_valid & wb_valid
            differ = both & (ls32.view(np.uint32) != wb32.view(np.uint32))
            if counts is not None:
                n_breached_runs += counts[2] > 0
                n_unsolved_runs += counts[3] > 0
            n_changed_runs += bool((both & (ls32 != z)).any())

            if not mask_ok or differ.any():
                n_fail += 1
                dem_failed = True
                first = ""
                if differ.any():
                    r, c = np.argwhere(differ)[0]
                    first = (f"({r},{c}) lsdtt3 {float(ls32[r, c])!r} "
                             f"whitebox {float(wb32[r, c])!r}")
                failures.append({
                    "dem": str(dem_path), "kind": kind, "cell_size": cs,
                    "shape": list(z.shape), "minimize_dist": md,
                    "max_dist": max_dist,
                    "max_cost": "unlimited" if max_cost is None else max_cost,
                    "mask_identical": bool(mask_ok),
                    "cells_differing": int(differ.sum()),
                    "first_differing_cell": first,
                    "lsdtt3_counts": counts,
                })
        if not dem_failed:
            shutil.rmtree(dem_dir)

    print(f"breach_least_cost fuzz vs Whitebox 2.0.6: {args.n_dems} DEMs "
          f"(seed {args.seed}), {n_runs} runs, {n_fail} mismatching")
    print(f"  coverage: {n_breached_runs} runs breached pits, "
          f"{n_unsolved_runs} runs left pits unsolved, "
          f"{n_changed_runs} runs changed the DEM, "
          f"{n_nodata_dems} DEMs with nodata")

    rc = 0
    if failures:
        (work / "failures.json").write_text(json.dumps(failures, indent=1))
        print(f"FAIL: {len(failures)} mismatching runs; details in "
              f"{work / 'failures.json'}")
        for f in failures[:10]:
            print(f"  - {json.dumps(f)}")
        rc = 1
    if n_breached_runs == 0 or n_unsolved_runs == 0 or (args.n_dems >= 20 and n_nodata_dems == 0):
        print("FAIL: coverage guard -- the fuzz did not exercise breached pits, "
              "unsolved pits and nodata")
        rc = 1
    print("RESULT: " + ("PASS" if rc == 0 else "FAIL"))
    return rc


if __name__ == "__main__":
    sys.exit(main())
