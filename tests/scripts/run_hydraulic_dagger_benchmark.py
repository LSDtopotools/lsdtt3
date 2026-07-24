#!/usr/bin/env python3
"""Scaffold: lsdtt hydraulic routing vs DAGGER reference benchmark.

This script is intentionally minimal for Phase 0 wiring:
- validates required paths/options
- exits 77 (skip) when DAGGER inputs are unavailable
- writes a CSV skeleton with run metadata

It is a benchmark harness scaffold, not a full parity driver yet.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import pathlib
import subprocess
import sys


SKIP_EXIT = 77
REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent


def die_skip(msg: str) -> int:
    print(f"SKIP: {msg}", file=sys.stderr)
    return SKIP_EXIT


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


def write_scaffold_csv(
    out_csv: pathlib.Path,
    *,
    fixture: str,
    scenario: str,
    lsdtt_sha: str,
    dagger_sha: str,
    dagger_root: pathlib.Path,
    dagger_ref: pathlib.Path,
) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "timestamp",
        "fixture",
        "scenario",
        "lsdtt_sha",
        "dagger_sha",
        "status",
        "note",
        "dagger_root",
        "dagger_ref",
    ]
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerow(
            {
                "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
                "fixture": fixture,
                "scenario": scenario,
                "lsdtt_sha": lsdtt_sha,
                "dagger_sha": dagger_sha,
                "status": "scaffold",
                "note": "Harness scaffold only; DAGGER parity execution not implemented yet.",
                "dagger_root": str(dagger_root),
                "dagger_ref": str(dagger_ref),
            }
        )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exe", required=True,
                    help="path to lsdtt-hydraulic-routing executable")
    ap.add_argument("--dagger-root", required=True,
                    help="path to DAGGER repository root")
    ap.add_argument("--dagger-ref", required=True,
                    help="path to DAGGER reference driver (script or executable)")
    ap.add_argument("--fixture", default="tiny_dem.asc",
                    help="fixture identifier for result metadata")
    ap.add_argument("--scenario", default="stationary_fixed",
                    help="scenario identifier for result metadata")
    ap.add_argument("--work-dir", default="build/hydraulic_benchmark",
                    help="output working directory")
    ap.add_argument("--write-json", action="store_true",
                    help="also write run_manifest.json in work-dir")
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

    lsdtt_sha = short_git_sha(REPO_ROOT)
    dagger_sha = short_git_sha(dagger_root)

    work_dir.mkdir(parents=True, exist_ok=True)
    out_csv = work_dir / "hydraulic_benchmark_results.csv"
    write_scaffold_csv(
        out_csv,
        fixture=args.fixture,
        scenario=args.scenario,
        lsdtt_sha=lsdtt_sha,
        dagger_sha=dagger_sha,
        dagger_root=dagger_root,
        dagger_ref=dagger_ref,
    )

    if args.write_json:
        manifest = {
            "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
            "exe": str(exe),
            "dagger_root": str(dagger_root),
            "dagger_ref": str(dagger_ref),
            "fixture": args.fixture,
            "scenario": args.scenario,
            "status": "scaffold",
        }
        (work_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"wrote {out_csv}")
    print("status: scaffold (parity execution not implemented yet)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
