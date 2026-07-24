#!/usr/bin/env python3
"""Server performance campaign runner (flow + surface).

Creates one timestamped campaign directory with:
  - hardware_manifest.json
  - run_manifest.json
  - flow_performance_raw.csv (if flow run enabled)
  - surface_benchmark_raw.csv (if surface run enabled)
  - scaling_summary.csv (derived strong-scaling metrics where possible)

This is intentionally report-only orchestration. It does not assert wall-clock
thresholds.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import pathlib
import platform
import shutil
import socket
import subprocess
import sys
import time
from typing import Any


REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent


def now_stamp() -> str:
    return dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def parse_threads(spec: str) -> list[int]:
    vals: list[int] = []
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        try:
            v = int(token)
        except ValueError as exc:
            raise ValueError(f"invalid thread count '{token}' in --threads") from exc
        if v < 1:
            raise ValueError("thread counts must be >= 1")
        vals.append(v)
    if not vals:
        raise ValueError("--threads produced an empty thread set")
    return sorted(set(vals))


def run_capture(cmd: list[str]) -> dict[str, Any]:
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        return {"ok": False, "error": "not found", "command": cmd}
    return {
        "ok": proc.returncode == 0,
        "returncode": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "command": cmd,
    }


def maybe_json_from_cmd(cmd: list[str]) -> Any:
    r = run_capture(cmd)
    if not r.get("ok"):
        return {"available": False, "error": r.get("error", "command failed")}
    try:
        return {"available": True, "data": json.loads(r["stdout"])}
    except Exception:
        return {"available": True, "data": r["stdout"].strip()}


def read_mem_total_kb() -> int | None:
    path = pathlib.Path("/proc/meminfo")
    if not path.exists():
        return None
    try:
        for line in path.read_text().splitlines():
            if line.startswith("MemTotal:"):
                parts = line.split()
                return int(parts[1])
    except Exception:
        return None
    return None


def detect_cpu_governors() -> dict[str, Any]:
    base = pathlib.Path("/sys/devices/system/cpu")
    governors: dict[str, int] = {}
    try:
        for child in base.glob("cpu[0-9]*/cpufreq/scaling_governor"):
            txt = child.read_text().strip()
            governors[txt] = governors.get(txt, 0) + 1
    except Exception:
        return {"available": False}
    if not governors:
        return {"available": False}
    return {"available": True, "counts": governors}


def git_context(repo_root: pathlib.Path) -> dict[str, Any]:
    sha = run_capture(["git", "-C", str(repo_root), "rev-parse", "--short", "HEAD"])
    dirty = run_capture(["git", "-C", str(repo_root), "status", "--porcelain"])
    out: dict[str, Any] = {"available": False, "sha": "unknown", "dirty": None}
    if sha.get("ok"):
        out["available"] = True
        out["sha"] = sha["stdout"].strip() or "unknown"
    if dirty.get("ok"):
        out["dirty"] = bool(dirty["stdout"].strip())
    return out


def read_cmake_cache(path: pathlib.Path) -> dict[str, Any]:
    cache = path / "CMakeCache.txt"
    if not cache.exists():
        return {"available": False}
    wanted = {
        "CMAKE_BUILD_TYPE",
        "CMAKE_CXX_COMPILER",
        "CMAKE_CXX_COMPILER_ID",
        "CMAKE_CXX_COMPILER_VERSION",
        "LSDTT_ENABLE_PERF_TESTS",
        "LSDTT_ENABLE_BENCHMARKS",
    }
    out: dict[str, Any] = {"available": True}
    for line in cache.read_text(errors="ignore").splitlines():
        if not line or line.startswith("//") or line.startswith("#"):
            continue
        if "=" not in line or ":" not in line:
            continue
        key_type, value = line.split("=", 1)
        key, _typ = key_type.split(":", 1)
        if key in wanted:
            out[key] = value
    return out


def capture_hardware_manifest(build_dir: pathlib.Path) -> dict[str, Any]:
    return {
        "captured_at": dt.datetime.now().isoformat(timespec="seconds"),
        "host": {
            "hostname": socket.gethostname(),
            "platform": platform.platform(),
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "python": sys.version,
        },
        "cpu": {
            "lscpu": maybe_json_from_cmd(["lscpu", "--json"]),
            "nproc": run_capture(["nproc"]).get("stdout", "").strip() or "unknown",
            "governors": detect_cpu_governors(),
        },
        "memory": {
            "mem_total_kb": read_mem_total_kb(),
        },
        "numa": {
            "numactl_hardware": maybe_json_from_cmd(["numactl", "--hardware"]),
        },
        "toolchain": {
            "cmake_cache": read_cmake_cache(build_dir),
            "git": git_context(REPO_ROOT),
        },
    }


def run_logged(
    cmd: list[str],
    cwd: pathlib.Path,
    log_dir: pathlib.Path,
    stem: str,
    env: dict[str, str] | None = None,
) -> dict[str, Any]:
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True, env=env)
    elapsed = time.perf_counter() - t0

    stdout_path = log_dir / f"{stem}.stdout.log"
    stderr_path = log_dir / f"{stem}.stderr.log"
    stdout_path.write_text(proc.stdout)
    stderr_path.write_text(proc.stderr)

    return {
        "command": cmd,
        "cwd": str(cwd),
        "returncode": proc.returncode,
        "elapsed_seconds": elapsed,
        "stdout_log": str(stdout_path),
        "stderr_log": str(stderr_path),
    }


def copy_if_exists(src: pathlib.Path, dst: pathlib.Path) -> bool:
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def derive_flow_scaling(flow_csv: pathlib.Path) -> list[dict[str, Any]]:
    if not flow_csv.exists():
        return []
    rows: list[dict[str, Any]] = []
    with flow_csv.open(newline="") as fh:
        reader = csv.DictReader(fh)
        data = list(reader)

    baseline: dict[tuple[str, str], float] = {}
    for r in data:
        if r.get("n_threads") == "1":
            baseline[(r.get("routine", ""), r.get("n_cells", ""))] = float(r["seconds"])

    for r in data:
        routine = r.get("routine", "")
        n_cells = r.get("n_cells", "")
        n_threads = int(r.get("n_threads", "0") or "0")
        seconds = float(r.get("seconds", "nan"))
        base = baseline.get((routine, n_cells))
        speedup = ""
        efficiency = ""
        if base is not None and n_threads > 0 and seconds > 0:
            s = base / seconds
            speedup = f"{s:.6f}"
            efficiency = f"{(s / n_threads):.6f}"
        rows.append(
            {
                "suite": "flow",
                "label": routine,
                "n_cells": n_cells,
                "n_threads": n_threads,
                "seconds": f"{seconds:.9g}",
                "speedup": speedup,
                "efficiency": efficiency,
                "throughput_mcells_per_s": r.get("mcells_per_second", ""),
                "source": str(flow_csv),
            }
        )
    return rows


def derive_surface_scaling(surface_csv: pathlib.Path) -> list[dict[str, Any]]:
    if not surface_csv.exists():
        return []
    rows: list[dict[str, Any]] = []
    with surface_csv.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            if r.get("reference") != "topotoolbox":
                continue
            mt = r.get("lsdtt_seconds_mt", "")
            t1 = r.get("lsdtt_seconds_1t", "")
            nth = r.get("lsdtt_threads_mt", "")
            if not mt or not t1 or not nth:
                continue
            try:
                mt_v = float(mt)
                t1_v = float(t1)
                n_threads = int(nth)
            except ValueError:
                continue
            if mt_v <= 0 or t1_v <= 0 or n_threads <= 0:
                continue
            speedup = t1_v / mt_v
            rows.append(
                {
                    "suite": "surface",
                    "label": r.get("comparison", ""),
                    "n_cells": r.get("n_cells", ""),
                    "n_threads": n_threads,
                    "seconds": f"{mt_v:.9g}",
                    "speedup": f"{speedup:.6f}",
                    "efficiency": f"{(speedup / n_threads):.6f}",
                    "throughput_mcells_per_s": "",
                    "source": str(surface_csv),
                }
            )
    return rows


def write_scaling_summary(path: pathlib.Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "suite",
        "label",
        "n_cells",
        "n_threads",
        "seconds",
        "speedup",
        "efficiency",
        "throughput_mcells_per_s",
        "source",
    ]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--flow-exe", default="build/tests/lsdtt-test-flow-performance",
                    help="path to flow performance executable")
    ap.add_argument("--surface-exe", default="build/lsdtt-surface-metrics",
                    help="path to lsdtt-surface-metrics")
    ap.add_argument("--fetch-exe", default="build/lsdtt-fetch-raster",
                    help="path to lsdtt-fetch-raster")
    ap.add_argument("--python-exe", default=sys.executable,
                    help="python executable for child benchmark scripts")
    ap.add_argument("--threads", default="1,2,4,8,16",
                    help="comma-separated thread list for scaling runs")
    ap.add_argument("--repeats", type=int, default=5,
                    help="timed repeats for surface benchmark")
    ap.add_argument("--surface-fixture", default="alpine_medium",
                    help="fixture key for run_surface_benchmark.py")
    ap.add_argument("--out-dir", default="build/perf_campaign",
                    help="campaign root directory")
    ap.add_argument("--flow-mode", choices=["quick", "full"], default="full",
                    help="flow harness mode")
    ap.add_argument("--flow-raster", default="",
                    help="optional raster path for flow harness --raster")
    ap.add_argument("--skip-flow", action="store_true")
    ap.add_argument("--skip-surface", action="store_true")
    ap.add_argument("--dry-run", action="store_true",
                    help="write manifests but do not execute benchmarks")
    args = ap.parse_args()

    try:
        thread_values = parse_threads(args.threads)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    campaign_root = pathlib.Path(args.out_dir)
    run_dir = campaign_root / now_stamp()
    logs_dir = run_dir / "logs"
    flow_dir = run_dir / "flow"
    surface_dir = run_dir / "surface"
    for d in (run_dir, logs_dir, flow_dir, surface_dir):
        d.mkdir(parents=True, exist_ok=True)

    build_dir = REPO_ROOT / "build"
    hw = capture_hardware_manifest(build_dir)
    (run_dir / "hardware_manifest.json").write_text(json.dumps(hw, indent=2))

    run_manifest: dict[str, Any] = {
        "created_at": dt.datetime.now().isoformat(timespec="seconds"),
        "repo_root": str(REPO_ROOT),
        "run_dir": str(run_dir),
        "threads_requested": thread_values,
        "repeats": args.repeats,
        "surface_fixture": args.surface_fixture,
        "flow_mode": args.flow_mode,
        "commands": [],
        "notes": [
            "Scaffold runner: flow thread-list support depends on lsdtt-test-flow-performance CLI support.",
            "Surface scaling is gathered by repeated run_surface_benchmark.py invocations with different --threads.",
        ],
    }

    flow_raw_csv = run_dir / "flow_performance_raw.csv"
    surface_raw_csv = run_dir / "surface_benchmark_raw.csv"

    if args.dry_run:
        print(f"dry-run: wrote manifests to {run_dir}")
        (run_dir / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2))
        return 0

    failures: list[str] = []

    if not args.skip_flow:
        flow_exe = (REPO_ROOT / args.flow_exe).resolve() if not pathlib.Path(args.flow_exe).is_absolute() else pathlib.Path(args.flow_exe)
        if not flow_exe.exists():
            failures.append(f"flow exe not found: {flow_exe}")
        else:
            cmd = [str(flow_exe), str(flow_dir), f"--{args.flow_mode}"]
            cmd.extend(["--threads", args.threads])
            if args.flow_raster:
                cmd.extend(["--raster", args.flow_raster])
            rec = run_logged(cmd, REPO_ROOT, logs_dir, "flow")
            rec["kind"] = "flow"
            run_manifest["commands"].append(rec)
            if rec["returncode"] != 0:
                failures.append(f"flow run failed with exit {rec['returncode']}")
            src = flow_dir / "flow_performance.csv"
            if not copy_if_exists(src, flow_raw_csv):
                failures.append(f"expected flow CSV not found: {src}")

    if not args.skip_surface:
        surface_exe = (REPO_ROOT / args.surface_exe).resolve() if not pathlib.Path(args.surface_exe).is_absolute() else pathlib.Path(args.surface_exe)
        fetch_exe = (REPO_ROOT / args.fetch_exe).resolve() if not pathlib.Path(args.fetch_exe).is_absolute() else pathlib.Path(args.fetch_exe)
        run_surface = SCRIPT_DIR / "run_surface_benchmark.py"
        if not surface_exe.exists():
            failures.append(f"surface exe not found: {surface_exe}")
        elif not run_surface.exists():
            failures.append(f"missing helper script: {run_surface}")
        else:
            # One run per thread count so each row set has an explicit scaling point.
            # run_surface_benchmark also records 1-thread timing in every run.
            for t in thread_values:
                work_dir_t = surface_dir / f"threads_{t}"
                cmd = [
                    args.python_exe,
                    str(run_surface),
                    "--exe", str(surface_exe),
                    "--fetch-exe", str(fetch_exe),
                    "--fixture", args.surface_fixture,
                    "--threads", str(t),
                    "--repeats", str(args.repeats),
                    "--work-dir", str(work_dir_t),
                ]
                rec = run_logged(cmd, REPO_ROOT, logs_dir, f"surface_t{t}")
                rec["kind"] = "surface"
                rec["threads"] = t
                run_manifest["commands"].append(rec)
                if rec["returncode"] != 0:
                    failures.append(f"surface run (threads={t}) failed with exit {rec['returncode']}")

            # Consolidate all generated surface CSV rows into one raw CSV.
            merged_rows: list[dict[str, Any]] = []
            merged_fields: list[str] = []
            for t in thread_values:
                csv_t = surface_dir / f"threads_{t}" / "surface_benchmark.csv"
                if not csv_t.exists():
                    continue
                with csv_t.open(newline="") as fh:
                    reader = csv.DictReader(fh)
                    if not merged_fields:
                        merged_fields = list(reader.fieldnames or [])
                    for row in reader:
                        row["campaign_threads_request"] = str(t)
                        merged_rows.append(row)

            if merged_rows:
                if "campaign_threads_request" not in merged_fields:
                    merged_fields.append("campaign_threads_request")
                with surface_raw_csv.open("w", newline="") as fh:
                    writer = csv.DictWriter(fh, fieldnames=merged_fields)
                    writer.writeheader()
                    writer.writerows(merged_rows)
            else:
                failures.append("no surface benchmark CSV rows were produced")

    scaling_rows: list[dict[str, Any]] = []
    scaling_rows.extend(derive_flow_scaling(flow_raw_csv))
    scaling_rows.extend(derive_surface_scaling(surface_raw_csv))
    write_scaling_summary(run_dir / "scaling_summary.csv", scaling_rows)

    (run_dir / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2))

    print(f"campaign directory: {run_dir}")
    print(f"hardware manifest : {run_dir / 'hardware_manifest.json'}")
    print(f"run manifest      : {run_dir / 'run_manifest.json'}")
    print(f"scaling summary   : {run_dir / 'scaling_summary.csv'}")

    if failures:
        print("\nFailures:", file=sys.stderr)
        for f in failures:
            print(f"- {f}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
