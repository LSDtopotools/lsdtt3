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
    times = []
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
        times.append(sum(ct) if ct else wall)
    return {"engine": "lsdtt3", "seconds": min(times),
            "seconds_all": times,
            "timing_scope": "cli_reported_compute" if times else "whole_process"}


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
    ap.add_argument("--cases", default="basins,chi")
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

    # -- output -------------------------------------------------------------
    fields = ["fixture", "case", "depth", "n_seeds", "engine", "seconds",
              "timing_scope", "result_is_slope_vs_depth", "phase_names",
              "phase_seconds", "error"]
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
    return 0


if __name__ == "__main__":
    sys.exit(main())
