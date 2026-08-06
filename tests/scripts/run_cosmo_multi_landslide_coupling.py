#!/usr/bin/env python3
"""End-to-end smoke/regression for multi-isotope step-change + landslide coupling.

This path used to be rejected by a CLI guard. The test proves that one run can
now:
  1) generate landslide rasters with lsdtt-cosmo-landslides,
  2) consume them in lsdtt-cosmo-prediction with
     calculate_multi_isotope_step_change + landslide_coupling,
  3) write all three concentration rasters and both ratio rasters,
  4) write sampled accumulated outputs whose R2610_ACC/R1410_ACC columns equal
     Al26_CONC/Be10_CONC and C14_CONC/Be10_CONC.
"""

import argparse
import csv
import math
import pathlib
import shutil
import subprocess
import sys


LANDSLIDE_PARAMS = """read_path: .
read_fname: tiny_dem.asc
write_path: .
write_fname: lsmi
landslide_mode: single_event
threshold_contributing_pixels: 3
area_volume_preset: larsen_soil_scar
event_time_yr: 1000
event_volume_m3: 200
n_min_pixels: 1
write_time_of_failure_raster: true
write_excavation_depth_raster: true
"""


PREDICTION_PARAMS = """read_path: .
read_fname: tiny_dem.asc
write_path: .
write_fname: mls

calculate_multi_isotope_step_change: true
calculate_accumulated_crn_concentration: true
sample_points_filename: cosmo_tiny_points.csv
threshold_contributing_pixels: 2
snap_search_radius_pixels: 3
step_change_eff_erate_old_g_cm2_yr: 0.001
step_change_eff_erate_new_g_cm2_yr: 0.005
time_since_change_yr: 20000

landslide_coupling: both
landslide_time_raster_prefix: lsmi_LS_TFAIL
landslide_excavation_raster_prefix: lsmi_LS_EXCAVATION
landslide_report_time_yr: 5000
flush_time_threshold_yr: 250
scar_density_kg_m3: 1500
"""


REL_TOL = 5.0e-6
ABS_TOL = 5.0e-7


def close(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=REL_TOL, abs_tol=ABS_TOL)


def run(exe: pathlib.Path, work: pathlib.Path, text: str, stem: str) -> str:
    param = work / f"{stem}.param"
    param.write_text(text)
    res = subprocess.run([str(exe), param.name], cwd=str(work),
                         capture_output=True, text=True, check=False)
    (work / f"{stem}.stdout").write_text(res.stdout)
    (work / f"{stem}.stderr").write_text(res.stderr)
    if res.returncode != 0:
        sys.stderr.write(res.stdout)
        sys.stderr.write(res.stderr)
        raise SystemExit(f"{stem} failed with exit {res.returncode}")
    return res.stdout


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--landslides-exe", required=True)
    ap.add_argument("--prediction-exe", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--work-dir", required=True)
    args = ap.parse_args()

    ls_exe = pathlib.Path(args.landslides_exe).resolve()
    pr_exe = pathlib.Path(args.prediction_exe).resolve()
    data_dir = pathlib.Path(args.data_dir).resolve()
    work = pathlib.Path(args.work_dir).resolve()

    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    for stem in ("tiny_dem.asc", "tiny_dem.prj", "cosmo_tiny_points.csv"):
        shutil.copy(data_dir / stem, work / stem)

    run(ls_exe, work, LANDSLIDE_PARAMS, "lsmi")
    for need in ("lsmi_LS_TFAIL.tif", "lsmi_LS_EXCAVATION.tif"):
        if not (work / need).exists():
            raise SystemExit(f"lsdtt-cosmo-landslides did not write {need}")

    out = run(pr_exe, work, PREDICTION_PARAMS, "mls")
    if "Landslide coupling: both" not in out:
        raise SystemExit("prediction run did not report landslide coupling")
    if "failed pixels" not in out:
        raise SystemExit("prediction run did not report coupling summary")

    for need in (
        "mls_Be10_step_crn_conc.tif",
        "mls_Al26_step_crn_conc.tif",
        "mls_C14_step_crn_conc.tif",
        "mls_R2610_step.tif",
        "mls_R1410_step.tif",
    ):
        if not (work / need).exists():
            raise SystemExit(f"prediction run did not write {need}")

    sampled = work / "mls_sampled_crn_conc.csv"
    if not sampled.exists():
        raise SystemExit("prediction run did not write sampled accumulated CSV")

    rows = list(csv.DictReader(sampled.open(newline="", encoding="utf-8")))
    if len(rows) != 1:
        raise SystemExit(f"expected one sampled row, got {len(rows)}")
    row = rows[0]
    for col in ("Be10_CONC", "Al26_CONC", "C14_CONC", "R2610_ACC", "R1410_ACC"):
        if col not in row:
            raise SystemExit(f"sampled CSV missing expected column {col}")

    be = float(row["Be10_CONC"])
    al = float(row["Al26_CONC"])
    c14 = float(row["C14_CONC"])
    r2610 = float(row["R2610_ACC"])
    r1410 = float(row["R1410_ACC"])
    if not (be > 0.0):
        raise SystemExit(f"Be10_CONC must be > 0, got {be}")
    exp2610 = al / be
    exp1410 = c14 / be
    if not close(r2610, exp2610):
        raise SystemExit(
            f"R2610_ACC mismatch: got {r2610:.10g}, expected {exp2610:.10g}")
    if not close(r1410, exp1410):
        raise SystemExit(
            f"R1410_ACC mismatch: got {r1410:.10g}, expected {exp1410:.10g}")

    print("OK: SC multi+landslides path wrote 3 nuclides + 2 ratios and sampled ACC ratios are ratio-of-means")
    return 0


if __name__ == "__main__":
    sys.exit(main())
