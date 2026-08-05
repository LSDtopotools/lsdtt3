#!/usr/bin/env python3
"""Rigorous SA2 check: sampled accumulated ratios are ratios of basin means.

Runs one multi-isotope step-change + accumulation invocation of
lsdtt-cosmo-prediction, reads the sampled CSV, and asserts:

  1) Be10_CONC / Al26_CONC / C14_CONC and R2610_ACC / R1410_ACC exist.
  2) R2610_ACC == Al26_CONC / Be10_CONC and R1410_ACC == C14_CONC / Be10_CONC
     at the sampled point (within float32/output precision tolerance).

Exit codes: 0 pass, 1 fail, 77 skip.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

EXIT_OK, EXIT_FAIL, EXIT_SKIP = 0, 1, 77

# Output is float32-backed rasters sampled then written as text.
REL_TOL = 5.0e-6
ABS_TOL = 5.0e-7


def fail(msg: str) -> int:
    print(f"FAIL: {msg}")
    return EXIT_FAIL


def run_cli(exe: pathlib.Path, workdir: pathlib.Path, param_file: pathlib.Path) -> None:
    proc = subprocess.run(
        [str(exe), str(workdir), param_file.name],
        cwd=str(workdir),
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"lsdtt-cosmo-prediction failed (exit {proc.returncode})\n"
            f"--- stdout ---\n{proc.stdout[-4000:]}\n"
            f"--- stderr ---\n{proc.stderr[-4000:]}"
        )

    needed = (
        "Using in-memory multi-isotope concentration rasters",
        "Accumulating multi-isotope concentration downstream",
    )
    for token in needed:
        if token not in proc.stdout:
            raise RuntimeError(
                "did not hit expected SA2 path; missing stdout token: " + token
            )


def close(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=REL_TOL, abs_tol=ABS_TOL)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred-exe", default=os.environ.get("LSDTT_COSMO_PREDICTION"))
    ap.add_argument(
        "--data-dir",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parents[1] / "data",
    )
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()

    if not args.pred_exe:
        print("SKIP: set --pred-exe (or LSDTT_COSMO_PREDICTION)")
        return EXIT_SKIP

    pred_exe = pathlib.Path(args.pred_exe)
    if not pred_exe.is_file():
        print(f"SKIP: executable not found: {pred_exe}")
        return EXIT_SKIP

    for fname in ("tiny_dem.asc", "cosmo_tiny_points.csv"):
        f = args.data_dir / fname
        if not f.is_file():
            return fail(f"fixture missing: {f}")

    workdir = pathlib.Path(tempfile.mkdtemp(prefix="cosmo_sa2_ratio_"))
    try:
        lines = [
            f"read_path: {args.data_dir}",
            "read_fname: tiny_dem.asc",
            f"write_path: {workdir}",
            "write_fname: sa2rig",
            "calculate_multi_isotope_step_change: true",
            "calculate_accumulated_crn_concentration: true",
            "sample_points_filename: cosmo_tiny_points.csv",
            "threshold_contributing_pixels: 2",
            "snap_search_radius_pixels: 3",
            "step_change_eff_erate_old_g_cm2_yr: 0.002",
            "step_change_eff_erate_new_g_cm2_yr: 0.02",
            "time_since_change_yr: 8000",
            "quartz_fertility: 0.7",
        ]
        param = workdir / "sa2rig.param"
        param.write_text("\n".join(lines) + "\n", encoding="utf-8")

        run_cli(pred_exe, workdir, param)

        sampled = workdir / "sa2rig_sampled_crn_conc.csv"
        if not sampled.is_file():
            return fail(f"missing sampled CSV: {sampled}")

        rows = list(csv.DictReader(sampled.open(newline="", encoding="utf-8")))
        if len(rows) != 1:
            return fail(f"expected 1 sampled row, got {len(rows)}")
        row = rows[0]

        for col in ("Be10_CONC", "Al26_CONC", "C14_CONC", "R2610_ACC", "R1410_ACC"):
            if col not in row:
                return fail(f"sampled CSV missing expected column: {col}")

        be = float(row["Be10_CONC"])
        al = float(row["Al26_CONC"])
        c14 = float(row["C14_CONC"])
        r2610 = float(row["R2610_ACC"])
        r1410 = float(row["R1410_ACC"])

        if not (be > 0.0):
            return fail(f"Be10_CONC must be > 0, got {be}")
        expect_r2610 = al / be
        expect_r1410 = c14 / be

        if not close(r2610, expect_r2610):
            return fail(
                f"R2610_ACC mismatch: got {r2610:.10g}, expected Al26_CONC/Be10_CONC "
                f"{expect_r2610:.10g}"
            )
        if not close(r1410, expect_r1410):
            return fail(
                f"R1410_ACC mismatch: got {r1410:.10g}, expected C14_CONC/Be10_CONC "
                f"{expect_r1410:.10g}"
            )

        print(
            "PASS: sampled accumulated ratios equal ratios of accumulated means "
            f"(R2610 {r2610:.8g}, R1410 {r1410:.8g})"
        )
        return EXIT_OK

    except Exception as exc:  # noqa: BLE001
        return fail(str(exc))
    finally:
        if args.keep:
            print(f"kept working directory: {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
