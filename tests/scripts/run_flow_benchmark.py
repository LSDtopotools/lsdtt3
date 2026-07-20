#!/usr/bin/env python3
"""Flow-routing benchmark driver: lsdtt3 vs RichDEM.

Phase 2 of docs/flow_routing_benchmark_plan.md.

Runs an lsdtt3 CLI and the corresponding RichDEM reference over the same fetched
fixture, computes the comparison metrics of plan Sec4, and writes one row per
comparison to benchmark_results.csv.

Phase 2 scope is the harness: metrics are computed and REPORTED.  Tolerances are
deliberately not asserted yet -- plan Sec8 Q2 requires pinning them just above
observed divergence once real numbers exist, otherwise a tolerance asserts
nothing.  Pass --check to enforce whatever tolerances have been pinned in
TOLERANCES so far.

Skips cleanly (exit 77, CTest's skip code) when RichDEM or the fixtures are
unavailable, so an optional dependency never turns CI red.

Usage:
    pixi run -e benchmark python tests/scripts/run_flow_benchmark.py \
        --exe build/lsdtt-raster-preprocessing
    python tests/scripts/run_flow_benchmark.py --list
"""

import argparse
import csv
import hashlib
import json
import pathlib
import shutil
import subprocess
import sys
import time

SKIP_EXIT = 77

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "dem_cache"
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent

# -----------------------------------------------------------------------------
# Fixture registry
#
# Fixtures are fetched by the lsdtt3 CLI itself (`lsdtt-fetch-raster`), not by a
# Python helper -- the Python fetcher was retired once the C++ path was shown to
# produce byte-identical output (plan Sec4a). Two terrain regimes, because they
# diverge very differently: see plan Sec2.3.
#
# Each fixture is a centre point plus a pixel count. The bounding box is widened
# by a margin so the exact centre crop cannot run off the edge onto nodata after
# reprojection.
# -----------------------------------------------------------------------------
FIXTURES = {
    "alpine":   {"centre": (10.40, 46.50), "epsg": 32632,
                 "description": "high relief, Italian Alps (Val Venosta / Ortler)"},
    "lakeland": {"centre": (25.40, 61.50), "epsg": 32635,
                 "description": "low relief, Finnish Lakeland (many closed basins)"},
}

SIZES = {"small": 256, "medium": 720, "large": 2400}

GRID_SPACING_M = 30.0


def fixture_bbox(name: str, size_px: int):
    """WGS84 bbox for a fixture, wide enough to survive reprojection."""
    import math
    lon, lat = FIXTURES[name]["centre"]
    half_m = 0.5 * size_px * GRID_SPACING_M * 1.25
    dlat = half_m / 111_320.0
    dlon = dlat / max(math.cos(math.radians(lat)), 1e-6)
    return (lon - dlon, lat - dlat, lon + dlon, lat + dlat)


def ensure_fixture(fixture_key: str, fetch_exe: pathlib.Path) -> pathlib.Path:
    """Fetch a fixture with lsdtt-fetch-raster if it is not already cached."""
    out_path = FIXTURE_DIR / f"{fixture_key}.tif"
    if out_path.exists():
        return out_path

    if "_" not in fixture_key:
        die_skip(f"fixture name {fixture_key!r} should be <name>_<size>")
    name, size_name = fixture_key.rsplit("_", 1)

    if name not in FIXTURES or size_name not in SIZES:
        die_skip(f"unknown fixture {fixture_key!r}; "
                 f"names={sorted(FIXTURES)} sizes={sorted(SIZES)}")

    if not fetch_exe or not pathlib.Path(fetch_exe).exists():
        die_skip(f"fixture {out_path} is missing and lsdtt-fetch-raster was not "
                 f"found (looked for {fetch_exe}). Build it, or pass --fetch-exe.")

    size_px = SIZES[size_name]
    west, south, east, north = fixture_bbox(name, size_px)

    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    param_dir = FIXTURE_DIR / "_fetch_params"
    param_dir.mkdir(parents=True, exist_ok=True)
    param_path = param_dir / f"{fixture_key}.param"

    param_path.write_text(
        f"write_path: {FIXTURE_DIR}\n"
        f"write_fname: {fixture_key}_raw\n"
        f"west: {west:.10f}\n"
        f"south: {south:.10f}\n"
        f"east: {east:.10f}\n"
        f"north: {north:.10f}\n"
        f"dem_source: cop30_aws\n"
        f"target_epsg: {FIXTURES[name]['epsg']}\n"
        f"grid_spacing: {GRID_SPACING_M}\n"
        f"size_pixels: {size_px}\n"
    )

    print(f"fetching fixture {fixture_key} ({FIXTURES[name]['description']})")
    proc = subprocess.run([str(fetch_exe), str(param_dir), param_path.name],
                          capture_output=True, text=True)

    if proc.returncode != 0:
        # An unreachable source is a missing dependency, not a test failure.
        die_skip(f"could not fetch fixture {fixture_key}:\n{proc.stdout[-1500:]}\n"
                 f"{proc.stderr[-1500:]}")

    written = FIXTURE_DIR / f"{fixture_key}_raw_DEM.tif"
    if not written.exists():
        die_skip(f"lsdtt-fetch-raster reported success but {written} is missing")
    written.rename(out_path)
    return out_path

# Comparison definitions.  `richdem_routine` names a routine in
# richdem_reference.py; `lsdtt_suffix` is the output file lsdtt3 writes.
#
# NOTE on fill (plan Sec4.2): lsdtt3's fill_min_slope defaults to 1e-4, which
# imposes an epsilon gradient across flats, while RichDEM's epsilon=False gives
# true flats.  These SHOULD disagree.  The tight comparison is min_slope=0 vs
# epsilon=False; the epsilon pair is a loose shape-agreement check only.
COMPARISONS = {
    "fill": {
        "richdem_routine": "fill",
        "lsdtt_suffix": "_fill.tif",
        "params": {"write_fill_raster": "true", "fill_min_slope": "0"},
        "kind": "elevation",
        "note": "Barnes 2014, no epsilon -- the tight comparison",
    },
    "fill_epsilon": {
        "richdem_routine": "fill_epsilon",
        "lsdtt_suffix": "_fill.tif",
        "params": {"write_fill_raster": "true", "fill_min_slope": "0.0001"},
        "kind": "elevation",
        "note": "epsilon gradients differ in detail; shape agreement only",
    },
    "breach": {
        "richdem_routine": "breach",
        "lsdtt_suffix": "_breached.tif",
        "params": {"write_breached_raster": "true",
                   "breach_algorithm_choice": "lindsay2016"},
        "kind": "elevation",
        "note": "both implementations pit-flatten before carving (plan Sec7a)",
    },
}

# Pinned tolerances (plan Sec4).  None = not yet pinned, reported only.
TOLERANCES = {
    "fill":         {"max_abs": None, "frac_differing": None},
    "fill_epsilon": {"max_abs": None, "frac_differing": None},
    "breach":       {"max_abs": None, "frac_differing": None},
}


def die_skip(msg: str):
    print(f"SKIP: {msg}", file=sys.stderr)
    sys.exit(SKIP_EXIT)


def load_raster(path: pathlib.Path):
    """Return (array, nodata_mask) as float64 / bool."""
    import numpy as np
    from osgeo import gdal
    gdal.UseExceptions()

    ds = gdal.Open(str(path))
    if ds is None:
        raise RuntimeError(f"could not open {path}")
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype("float64")
    nodata = band.GetNoDataValue()
    valid = np.isfinite(arr)
    if nodata is not None:
        valid &= (arr != nodata)
    return arr, valid


def elevation_metrics(a_path, b_path, tol=0.01):
    """Plan Sec4.1 metrics for two elevation rasters."""
    import numpy as np

    a, va = load_raster(a_path)
    b, vb = load_raster(b_path)

    if a.shape != b.shape:
        raise RuntimeError(f"shape mismatch {a.shape} vs {b.shape}")

    # The nodata masks must match exactly.  A mismatch is a bug, never
    # something to absorb into a tolerance (plan Sec4.1).
    mask_identical = bool((va == vb).all())
    both = va & vb
    n = int(both.sum())
    if n == 0:
        raise RuntimeError("no overlapping valid cells")

    d = a[both] - b[both]
    ad = np.abs(d)

    # Fill volume ratio: integral check that catches systematic offsets a
    # cell-wise max would miss.
    sum_a, sum_b = float(a[both].sum()), float(b[both].sum())
    volume_ratio = (sum_a / sum_b) if sum_b != 0 else float("nan")

    return {
        "mask_identical": mask_identical,
        "n_valid": n,
        "max_abs": float(ad.max()),
        "mean_abs": float(ad.mean()),
        "frac_differing": float((ad > tol).sum() / n),
        "signed_min": float(d.min()),
        "signed_max": float(d.max()),
        "volume_ratio": volume_ratio,
    }


def run_lsdtt(exe: pathlib.Path, fixture: pathlib.Path, out_dir: pathlib.Path,
              params: dict) -> float:
    """Run the lsdtt3 CLI on the fixture. Returns wall time in seconds."""
    out_dir.mkdir(parents=True, exist_ok=True)

    base = {
        "read_path": str(fixture.parent),
        "read_fname": fixture.name,
        "write_path": str(out_dir),
        "write_fname": fixture.stem,
        # The fixture is already in metric UTM (fetched that way). Letting
        # the CLI reproject would score the two tools on different rasters.
        "enforce_projection": "false",
        "write_output_format": "gtiff",
    }
    base.update(params)

    param_path = out_dir / "benchmark.param"
    param_path.write_text("".join(f"{k}: {v}\n" for k, v in base.items()))

    # The CLI takes the parameter file as <directory> <filename>, not a single
    # path (see tests/scripts/run_cli_smoke.py).
    t0 = time.perf_counter()
    proc = subprocess.run([str(exe), str(out_dir), param_path.name],
                          capture_output=True, text=True)
    elapsed = time.perf_counter() - t0

    if proc.returncode != 0:
        raise RuntimeError(
            f"{exe.name} failed (exit {proc.returncode})\n"
            f"--- stdout ---\n{proc.stdout[-2000:]}\n"
            f"--- stderr ---\n{proc.stderr[-2000:]}")
    return elapsed


def run_richdem(routine: str, fixture: pathlib.Path, out_path: pathlib.Path,
                cache: bool = True):
    """Run the RichDEM reference. Returns (seconds, cached_bool)."""
    if cache and out_path.exists():
        meta_path = out_path.with_suffix(".json")
        if meta_path.exists():
            return json.loads(meta_path.read_text())["seconds"], True

    out_path.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT_DIR / "richdem_reference.py"),
         "--routine", routine, "--in", str(fixture), "--out", str(out_path)],
        capture_output=True, text=True)

    if proc.returncode == SKIP_EXIT:
        die_skip("richdem not available. Use `pixi run -e benchmark ...`")
    if proc.returncode != 0:
        raise RuntimeError(f"richdem_reference.py failed:\n{proc.stderr[-2000:]}")

    info = json.loads(proc.stdout.strip().splitlines()[-1])
    out_path.with_suffix(".json").write_text(json.dumps(info))
    return info["seconds"], False


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", help="path to lsdtt-raster-preprocessing")
    ap.add_argument("--fetch-exe",
                    help="path to lsdtt-fetch-raster (used to fetch missing "
                         "fixtures; defaults to build/lsdtt-fetch-raster)")
    ap.add_argument("--fixture", default="alpine_small",
                    help="fixture basename in tests/data/dem_cache (default alpine_small)")
    ap.add_argument("--comparison", choices=list(COMPARISONS) + ["all"], default="all")
    ap.add_argument("--work-dir", default=None)
    ap.add_argument("--check", action="store_true",
                    help="fail on pinned tolerance violations (Phase 3)")
    ap.add_argument("--list", action="store_true", help="list comparisons and exit")
    args = ap.parse_args()

    if args.list:
        for name, spec in COMPARISONS.items():
            print(f"{name:14s} {spec['note']}")
        return 0

    try:
        import numpy  # noqa: F401
        from osgeo import gdal  # noqa: F401
    except ImportError:
        die_skip("numpy/GDAL not importable. Use `pixi run -e benchmark ...`")

    fetch_exe = pathlib.Path(args.fetch_exe).resolve() if args.fetch_exe else \
        (REPO_ROOT / "build" / "lsdtt-fetch-raster")
    fixture = ensure_fixture(args.fixture, fetch_exe)

    # Short content hash of the fixture, used in the RichDEM cache key below.
    fixture_digest = hashlib.sha256(fixture.read_bytes()).hexdigest()[:12]

    if not args.exe:
        die_skip("--exe not given")
    exe = pathlib.Path(args.exe).resolve()
    if not exe.exists():
        die_skip(f"executable {exe} not found; build lsdtt3 first")

    work_dir = pathlib.Path(args.work_dir).resolve() if args.work_dir else \
        REPO_ROOT / "build" / "flow_benchmark"
    work_dir.mkdir(parents=True, exist_ok=True)

    names = list(COMPARISONS) if args.comparison == "all" else [args.comparison]

    rows, failures = [], []
    for name in names:
        spec = COMPARISONS[name]
        print(f"\n=== {name} on {args.fixture} ===")
        print(f"    {spec['note']}")

        case_dir = work_dir / f"{args.fixture}_{name}"
        if case_dir.exists():
            shutil.rmtree(case_dir)

        lsdtt_s = run_lsdtt(exe, fixture, case_dir, spec["params"])
        lsdtt_out = case_dir / f"{fixture.stem}{spec['lsdtt_suffix']}"
        if not lsdtt_out.exists():
            raise RuntimeError(f"expected lsdtt3 output {lsdtt_out} not written")

        # Key the cache on the fixture's CONTENT, not just its name. Keying on
        # the name alone silently reuses a reference computed from a different
        # raster whenever the fixture is regenerated -- which happened when
        # fixture generation moved from pixi's GDAL to the C++ tool's, and
        # produced a bogus "99.9 % of cells differ" instead of the true 0 %.
        ref_out = (work_dir / "richdem_cache" /
                   f"{args.fixture}_{fixture_digest}_{spec['richdem_routine']}.tif")
        rd_s, cached = run_richdem(spec["richdem_routine"], fixture, ref_out)

        m = elevation_metrics(lsdtt_out, ref_out)

        print(f"    lsdtt3   {lsdtt_s:8.3f} s")
        print(f"    richdem  {rd_s:8.3f} s{'  (cached)' if cached else ''}")
        print(f"    mask identical:  {m['mask_identical']}")
        print(f"    max |diff|:      {m['max_abs']:.4f} m")
        print(f"    mean |diff|:     {m['mean_abs']:.6f} m")
        print(f"    cells >1cm:      {100.0 * m['frac_differing']:.4f} %")
        print(f"    signed range:    [{m['signed_min']:.4f}, {m['signed_max']:.4f}] m")
        print(f"    volume ratio:    {m['volume_ratio']:.9f}")

        # Record which fixture build produced the row: benchmark_results.csv is
        # append-only, and rows from different fixtures are otherwise
        # indistinguishable and easy to compare by mistake.
        row = {"fixture": args.fixture, "fixture_sha": fixture_digest,
               "comparison": name,
               "lsdtt_seconds": lsdtt_s, "richdem_seconds": rd_s, **m}
        rows.append(row)

        if not m["mask_identical"]:
            failures.append(f"{name}: nodata masks differ")

        if args.check:
            tol = TOLERANCES.get(name, {})
            for key, limit in tol.items():
                if limit is not None and m[key] > limit:
                    failures.append(
                        f"{name}: {key} = {m[key]:.6g} exceeds pinned {limit:g}")

    csv_path = work_dir / "benchmark_results.csv"
    write_header = not csv_path.exists()
    with open(csv_path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        if write_header:
            w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {csv_path}")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        return 1

    if not args.check:
        print("\n(reporting only -- tolerances not pinned yet, plan Sec8 Q2)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
