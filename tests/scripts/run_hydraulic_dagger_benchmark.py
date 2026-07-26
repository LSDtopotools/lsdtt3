#!/usr/bin/env python3
"""Hydraulic benchmark: lsdtt3 vs DAGGER GraphFlood reference.

Plan document: docs/plan_docs/hydraulic_routing_dagger_validation_plan.md

This Phase 3 harness runs matched lsdtt3 and DAGGER scenarios, compares key
output rasters, and writes one metrics row to hydraulic_benchmark_results.csv.

Default mode is report-only: metrics are recorded, no numeric tolerances are
asserted. Use --check to enforce hard comparison preconditions and pinned
Phase-4 parity tolerances.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import pathlib
import subprocess
import sys
import time
import hashlib

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import run_flow_benchmark as flow  # noqa: E402

try:
    import numpy as np
    from osgeo import gdal
except Exception:  # pragma: no cover - optional benchmark dependency
    np = None
    gdal = None


SKIP_EXIT = 77
REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent


SCENARIOS = {
    "stationary_fixed": {
        "hydraulic_run_mode": "stationary",
        "dynamic_routing_enabled": "false",
        "run_until_converged": "false",
        "n_steps": "30",
        "reroute_every": "1",
        "hydraulic_dt": "1e-3",
        "mannings_n": "0.033",
        "hydraulic_min_slope": "1e-8",
        "uniform_rainfall_rate": "1e-6",
        "flow_routing_method": "mfd",
        "mfd_exponent": "1.0",
    },
    "transient_dynamic": {
        "hydraulic_run_mode": "transient",
        "dynamic_routing_enabled": "true",
        "run_until_converged": "false",
        "n_steps": "30",
        "reroute_every": "1",
        "hydraulic_dt": "1e-3",
        "mannings_n": "0.033",
        "hydraulic_min_slope": "1e-8",
        "uniform_rainfall_rate": "1e-6",
        "flow_routing_method": "mfd",
        "mfd_exponent": "1.0",
    },
}


# Phase-4 pinned tolerances (2026-07-26) on tiny_dem.asc baselines.
PARITY_TOLERANCES = {
    "tiny_dem.asc": {
        "stationary_fixed": {
            "water_depth_max_abs": 2.5e-5,
            "water_depth_rmse": 1.1e-5,
            "water_surface_rmse": 1.1e-5,
            "qout_rmse": 5.0e-7,
            "storage_abs_diff": 1.5e-4,
        },
        "transient_dynamic": {
            "water_depth_max_abs": 1.0e-12,
            "water_depth_rmse": 1.0e-12,
            "water_surface_rmse": 1.0e-12,
            "qout_rmse": 1.0e-12,
            "storage_abs_diff": 1.0e-12,
        },
    },
    "alpine_small": {
        "stationary_fixed": {
            "water_depth_max_abs": 0.07,
            "water_depth_rmse": 0.0035,
            "water_surface_rmse": 2.1,
            "qout_rmse": 0.2,
            "storage_abs_diff": 1.5e4,
        },
        "transient_dynamic": {
            "water_depth_max_abs": 0.002,
            "water_depth_rmse": 8.0e-5,
            "water_surface_rmse": 2.1,
            "qout_rmse": 0.2,
            "storage_abs_diff": 400.0,
        },
    },
}


EXPECTED_FIXTURE_SHA = {
    "alpine_small": "899cc313e27c",
}


class SkipDependency(RuntimeError):
    pass


def die_skip(msg: str) -> int:
    print(f"SKIP: {msg}", file=sys.stderr)
    return SKIP_EXIT


def system_env() -> dict[str, str]:
    env = dict(__import__("os").environ)
    for key in (
        "LD_LIBRARY_PATH",
        "DYLD_LIBRARY_PATH",
        "GDAL_DATA",
        "GDAL_DRIVER_PATH",
        "PROJ_DATA",
        "PROJ_LIB",
        "CONDA_PREFIX",
    ):
        env.pop(key, None)
    return env


def fixture_sha12(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:12]


def resolve_fixture(fixture_arg: str, fetch_exe_arg: str | None) -> tuple[pathlib.Path, str, str]:
    direct = pathlib.Path(fixture_arg)
    if direct.exists():
        path = direct.resolve()
        return path, path.name, fixture_sha12(path)

    tests_data = REPO_ROOT / "tests" / "data"
    in_data = tests_data / fixture_arg
    if in_data.exists():
        path = in_data.resolve()
        return path, fixture_arg, fixture_sha12(path)

    if fixture_arg.endswith(".tif") and (tests_data / "dem_cache" / fixture_arg).exists():
        path = (tests_data / "dem_cache" / fixture_arg).resolve()
        return path, fixture_arg[:-4], fixture_sha12(path)

    in_cache = tests_data / "dem_cache" / f"{fixture_arg}.tif"
    if in_cache.exists():
        path = in_cache.resolve()
        return path, fixture_arg, fixture_sha12(path)

    fetch_exe = pathlib.Path(fetch_exe_arg).resolve() if fetch_exe_arg else (REPO_ROOT / "build" / "lsdtt-fetch-raster")
    try:
        fetched = flow.ensure_fixture(fixture_arg, fetch_exe)
    except SystemExit as exc:
        code = int(exc.code) if isinstance(exc.code, int) else 1
        raise SkipDependency(f"fixture unavailable ({fixture_arg}), fetch exit {code}") from exc

    path = pathlib.Path(fetched).resolve()
    return path, fixture_arg, fixture_sha12(path)


def short_git_sha(path: pathlib.Path) -> str:
    try:
        proc = subprocess.run(
            ["git", "-C", str(path), "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        return "unknown"
    if proc.returncode != 0:
        return "unknown"
    sha = proc.stdout.strip()
    return sha if sha else "unknown"


def write_csv_row(
    out_csv: pathlib.Path,
    row: dict,
) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = list(row.keys())
    write_header = not out_csv.exists()
    with out_csv.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        if write_header:
            w.writeheader()
        w.writerow(row)


def load_raster(path: pathlib.Path) -> tuple[np.ndarray, np.ndarray, tuple, str, float]:
    gdal.UseExceptions()
    ds = gdal.Open(str(path))
    if ds is None:
        raise RuntimeError(f"could not open raster: {path}")

    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype("float64")
    ndv = band.GetNoDataValue()

    valid = np.isfinite(arr)
    if ndv is not None:
        valid &= arr != ndv

    gt = ds.GetGeoTransform()
    proj = ds.GetProjectionRef() or ""
    cell_area = abs(gt[1] * gt[5])

    signature = (ds.RasterXSize, ds.RasterYSize, tuple(gt))
    return arr, valid, signature, proj, cell_area


def rankdata_simple(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty_like(order, dtype="float64")
    ranks[order] = np.arange(order.size, dtype="float64")
    return ranks


def corrcoef_safe(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 2:
        return float("nan")
    sa = float(np.std(a))
    sb = float(np.std(b))
    if sa == 0.0 or sb == 0.0:
        return 1.0 if np.allclose(a, b) else float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def compare_rasters(lsdtt_path: pathlib.Path, dagger_path: pathlib.Path) -> dict:
    a, va, sig_a, proj_a, _ = load_raster(lsdtt_path)
    b, vb, sig_b, proj_b, _ = load_raster(dagger_path)

    if sig_a != sig_b:
        raise RuntimeError(f"grid mismatch: {lsdtt_path.name} vs {dagger_path.name}")
    if proj_a != proj_b:
        raise RuntimeError(f"projection mismatch: {lsdtt_path.name} vs {dagger_path.name}")

    mask_identical = bool(np.array_equal(va, vb))
    both = va & vb
    n_valid = int(np.sum(both))
    if n_valid == 0:
        raise RuntimeError(f"no overlapping valid cells: {lsdtt_path.name}")

    d = a[both] - b[both]
    ad = np.abs(d)
    rmse = float(np.sqrt(np.mean(d * d)))
    pearson = corrcoef_safe(a[both], b[both])
    spearman = corrcoef_safe(rankdata_simple(a[both]), rankdata_simple(b[both]))

    return {
        "mask_identical": mask_identical,
        "n_valid": n_valid,
        "max_abs": float(np.max(ad)),
        "rmse": rmse,
        "pearson": pearson,
        "spearman": spearman,
    }


def water_volume(depth_path: pathlib.Path) -> float:
    depth, valid, _, _, cell_area = load_raster(depth_path)
    return float(np.sum(depth[valid] * cell_area))


def run_lsdtt(
    exe: pathlib.Path,
    fixture_path: pathlib.Path,
    scenario_key: str,
    out_dir: pathlib.Path,
) -> tuple[dict[str, pathlib.Path], float]:
    scenario = SCENARIOS[scenario_key]
    out_dir.mkdir(parents=True, exist_ok=True)

    params = {
        "read_path": str(fixture_path.parent),
        "read_fname": fixture_path.name,
        "write_path": str(out_dir),
        "write_fname": "lsdtt",
        "write_output_format": "gtiff",
        "n_threads": "1",
        "raster_is_filled": "false",
        "fill_min_slope": "0.0001",
        "boundary_conditions": "oooo",
        "write_water_depth_raster": "true",
        "write_water_surface_raster": "true",
        "write_qin_raster": "false",
        "write_qout_raster": "true",
        "write_dhw_raster": "false",
        "write_bed_raster": "false",
        "write_hydraulic_slope_raster": "false",
        "write_specific_discharge_raster": "false",
        "write_total_discharge_raster": "false",
        "write_specific_drainage_area_raster": "false",
        "write_effective_width_raster": "false",
    }
    params.update(scenario)

    param_name = "lsdtt_hydraulic_benchmark.param"
    param_path = out_dir / param_name
    param_path.write_text("".join(f"{k}: {v}\n" for k, v in params.items()), encoding="utf-8")

    cmd = [str(exe), str(out_dir), param_name]
    t0 = time.perf_counter()
    proc = subprocess.run(
        cmd,
        cwd=str(out_dir),
        capture_output=True,
        text=True,
        check=False,
        env=system_env(),
    )
    elapsed = time.perf_counter() - t0

    (out_dir / "lsdtt_stdout.log").write_text(proc.stdout, encoding="utf-8")
    (out_dir / "lsdtt_stderr.log").write_text(proc.stderr, encoding="utf-8")

    if proc.returncode != 0:
        raise RuntimeError(
            f"lsdtt-hydraulic-routing failed (exit {proc.returncode})\n"
            f"--- stdout ---\n{proc.stdout[-2000:]}\n"
            f"--- stderr ---\n{proc.stderr[-2000:]}"
        )

    outputs = {
        "water_depth": out_dir / "lsdtt_water_depth.tif",
        "water_surface": out_dir / "lsdtt_water_surface.tif",
        "qout": out_dir / "lsdtt_qout.tif",
    }
    for field, path in outputs.items():
        if not path.exists():
            raise RuntimeError(f"lsdtt output missing: {field} -> {path}")

    return outputs, elapsed


def run_dagger_reference(
    dagger_ref: pathlib.Path,
    dagger_root: pathlib.Path,
    fixture_path: pathlib.Path,
    scenario_key: str,
    out_dir: pathlib.Path,
    timeout_seconds: int,
) -> tuple[dict[str, pathlib.Path], float]:
    scenario = SCENARIOS[scenario_key]
    out_dir.mkdir(parents=True, exist_ok=True)

    if dagger_ref.suffix.lower() == ".py":
        launcher = [sys.executable, str(dagger_ref)]
    elif dagger_ref.is_file() and dagger_ref.stat().st_mode & 0o111:
        launcher = [str(dagger_ref)]
    else:
        raise SkipDependency(
            f"DAGGER reference is not executable/script: {dagger_ref}. "
            "Pass a .py script or executable implementing the benchmark interface."
        )

    cmd = launcher + [
        "--dagger-root", str(dagger_root),
        "--dem", str(fixture_path),
        "--scenario", scenario_key,
        "--out-dir", str(out_dir),
        "--out-prefix", "dagger",
        "--n-steps", scenario["n_steps"],
        "--dt", scenario["hydraulic_dt"],
        "--mannings-n", scenario["mannings_n"],
        "--rainfall-rate", scenario["uniform_rainfall_rate"],
        "--flow-method", scenario["flow_routing_method"],
        "--mfd-exponent", scenario["mfd_exponent"],
        "--dynamic-routing", scenario["dynamic_routing_enabled"],
        "--reroute-every", scenario["reroute_every"],
        "--run-mode", scenario["hydraulic_run_mode"],
    ]

    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(out_dir),
            capture_output=True,
            text=True,
            check=False,
            env=dict(__import__("os").environ),
            timeout=timeout_seconds,
        )
    except OSError as exc:
        raise SkipDependency(f"could not launch DAGGER reference: {exc}") from exc
    except subprocess.TimeoutExpired as exc:
        raise SkipDependency(
            f"DAGGER reference timed out after {timeout_seconds}s on {fixture_path.name}"
        ) from exc
    elapsed = time.perf_counter() - t0

    (out_dir / "dagger_stdout.log").write_text(proc.stdout, encoding="utf-8")
    (out_dir / "dagger_stderr.log").write_text(proc.stderr, encoding="utf-8")

    if proc.returncode == SKIP_EXIT:
        raise SkipDependency("DAGGER runner reported unavailable dependency (exit 77).")
    if proc.returncode != 0:
        raise RuntimeError(
            f"DAGGER reference failed (exit {proc.returncode})\n"
            f"--- stdout ---\n{proc.stdout[-2000:]}\n"
            f"--- stderr ---\n{proc.stderr[-2000:]}"
        )

    outputs = {
        "water_depth": out_dir / "dagger_water_depth.tif",
        "water_surface": out_dir / "dagger_water_surface.tif",
        "qout": out_dir / "dagger_qout.tif",
    }
    for field, path in outputs.items():
        if not path.exists():
            raise RuntimeError(f"DAGGER output missing: {field} -> {path}")

    return outputs, elapsed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exe", required=True,
                    help="path to lsdtt-hydraulic-routing executable")
    ap.add_argument("--dagger-root", required=True,
                    help="path to DAGGER repository root")
    ap.add_argument("--dagger-ref", required=True,
                    help="path to DAGGER reference driver (script or executable)")
    ap.add_argument("--fixture", default="tiny_dem.asc",
                    help="fixture key or raster path (tiny_dem.asc, alpine_small)")
    ap.add_argument("--fetch-exe", default=None,
                    help="path to lsdtt-fetch-raster for on-demand dem_cache fixture fetch")
    ap.add_argument("--scenario", default="stationary_fixed",
                    choices=sorted(SCENARIOS.keys()),
                    help="matched hydraulic scenario to run")
    ap.add_argument("--work-dir", default="build/hydraulic_benchmark",
                    help="output working directory")
    ap.add_argument("--write-json", action="store_true",
                    help="also write run_manifest.json in work-dir")
    ap.add_argument("--check", action="store_true",
                    help="enforce hard parity preconditions (grid/mask match)")
    ap.add_argument("--dagger-timeout-seconds", type=int, default=300,
                    help="timeout for DAGGER reference subprocess")
    args = ap.parse_args()

    exe = pathlib.Path(args.exe).resolve()
    dagger_root = pathlib.Path(args.dagger_root).resolve()
    dagger_ref = pathlib.Path(args.dagger_ref).resolve()
    work_dir = pathlib.Path(args.work_dir).resolve()

    if not exe.exists():
        return die_skip(f"lsdtt executable not found: {exe}")
    if not dagger_root.exists():
        return die_skip(f"DAGGER root not found: {dagger_root}")
    if not dagger_ref.exists():
        return die_skip(f"DAGGER reference driver not found: {dagger_ref}")

    if np is None or gdal is None:
        return die_skip("numpy/GDAL not importable in this Python environment")

    try:
        fixture_path, fixture_key, fixture_sha = resolve_fixture(args.fixture, args.fetch_exe)
    except SkipDependency as exc:
        return die_skip(str(exc))
    except Exception as exc:
        return die_skip(f"fixture resolution failed for {args.fixture}: {exc}")

    lsdtt_sha = short_git_sha(REPO_ROOT)
    dagger_sha = short_git_sha(dagger_root)

    work_dir.mkdir(parents=True, exist_ok=True)
    out_csv = work_dir / "hydraulic_benchmark_results.csv"

    lsdtt_dir = work_dir / "lsdtt"
    dagger_dir = work_dir / "dagger"

    try:
        lsdtt_outputs, lsdtt_seconds = run_lsdtt(exe, fixture_path, args.scenario, lsdtt_dir)
    except Exception as exc:
        print(f"FAIL: lsdtt run failed: {exc}", file=sys.stderr)
        return 1

    try:
        dagger_outputs, dagger_seconds = run_dagger_reference(
            dagger_ref, dagger_root, fixture_path, args.scenario, dagger_dir,
            args.dagger_timeout_seconds,
        )
    except SkipDependency as exc:
        return die_skip(f"DAGGER run unavailable: {exc}")
    except Exception as exc:
        print(f"FAIL: DAGGER run failed: {exc}", file=sys.stderr)
        return 1

    depth = compare_rasters(lsdtt_outputs["water_depth"], dagger_outputs["water_depth"])
    surface = compare_rasters(lsdtt_outputs["water_surface"], dagger_outputs["water_surface"])
    qout = compare_rasters(lsdtt_outputs["qout"], dagger_outputs["qout"])

    storage_delta_lsdtt = water_volume(lsdtt_outputs["water_depth"])
    storage_delta_dagger = water_volume(dagger_outputs["water_depth"])
    storage_delta_diff = storage_delta_lsdtt - storage_delta_dagger

    run_status = "phase3_report_only"
    if args.check:
        if not (depth["mask_identical"] and surface["mask_identical"] and qout["mask_identical"]):
            print("FAIL: nodata masks differ between lsdtt3 and DAGGER outputs", file=sys.stderr)
            return 1

        expected_sha = EXPECTED_FIXTURE_SHA.get(fixture_key)
        if expected_sha is not None and fixture_sha != expected_sha:
            return die_skip(
                f"fixture {fixture_key} sha mismatch: got {fixture_sha}, expected {expected_sha}; "
                "regenerate/pin tolerances for this fixture content"
            )

        tol = PARITY_TOLERANCES.get(fixture_key, {}).get(args.scenario)
        if tol is None:
            print(
                f"FAIL: no pinned parity tolerances for fixture '{fixture_key}' scenario '{args.scenario}'",
                file=sys.stderr,
            )
            return 1

        checks = [
            ("water_depth_max_abs", depth["max_abs"], tol["water_depth_max_abs"]),
            ("water_depth_rmse", depth["rmse"], tol["water_depth_rmse"]),
            ("water_surface_rmse", surface["rmse"], tol["water_surface_rmse"]),
            ("qout_rmse", qout["rmse"], tol["qout_rmse"]),
            ("storage_abs_diff", abs(storage_delta_diff), tol["storage_abs_diff"]),
        ]
        failures = [
            f"{name}={actual:.10g} > {limit:.10g}"
            for name, actual, limit in checks
            if not np.isfinite(actual) or actual > limit
        ]
        if failures:
            print("FAIL: hydraulic parity tolerances exceeded:", file=sys.stderr)
            for line in failures:
                print(f"  - {line}", file=sys.stderr)
            return 1

        run_status = "phase4_checked"

    row = {
        "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
        "fixture": fixture_key,
        "fixture_sha": fixture_sha,
        "scenario": args.scenario,
        "lsdtt_sha": lsdtt_sha,
        "dagger_sha": dagger_sha,
        "n_steps": int(SCENARIOS[args.scenario]["n_steps"]),
        "dt": float(SCENARIOS[args.scenario]["hydraulic_dt"]),
        "water_depth_mask_identical": depth["mask_identical"],
        "water_surface_mask_identical": surface["mask_identical"],
        "qout_mask_identical": qout["mask_identical"],
        "water_depth_max_abs": depth["max_abs"],
        "water_depth_rmse": depth["rmse"],
        "water_surface_rmse": surface["rmse"],
        "qout_rmse": qout["rmse"],
        "water_depth_pearson": depth["pearson"],
        "water_depth_spearman": depth["spearman"],
        "water_surface_pearson": surface["pearson"],
        "water_surface_spearman": surface["spearman"],
        "qout_pearson": qout["pearson"],
        "qout_spearman": qout["spearman"],
        "storage_delta_lsdtt": storage_delta_lsdtt,
        "storage_delta_dagger": storage_delta_dagger,
        "storage_delta_diff": storage_delta_diff,
        "lsdtt_seconds": lsdtt_seconds,
        "dagger_seconds": dagger_seconds,
        "dagger_root": str(dagger_root),
        "dagger_ref": str(dagger_ref),
        "notes": run_status,
    }
    write_csv_row(out_csv, row)

    if args.write_json:
        manifest = {
            "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
            "exe": str(exe),
            "dagger_root": str(dagger_root),
            "dagger_ref": str(dagger_ref),
            "fixture_path": str(fixture_path),
            "fixture": fixture_key,
            "fixture_sha": fixture_sha,
            "scenario": args.scenario,
            "metrics": {
                "water_depth": depth,
                "water_surface": surface,
                "qout": qout,
                "storage_delta_lsdtt": storage_delta_lsdtt,
                "storage_delta_dagger": storage_delta_dagger,
                "storage_delta_diff": storage_delta_diff,
            },
            "timing_seconds": {
                "lsdtt": lsdtt_seconds,
                "dagger": dagger_seconds,
            },
            "status": run_status,
        }
        (work_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"wrote {out_csv}")
    print(f"status: {run_status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
