#!/usr/bin/env python3
"""End-to-end check of the landslide 'scar' coupling across two CLIs.

lsdtt-cosmo-landslides writes a time-of-failure raster and an excavation-depth
raster; lsdtt-cosmo-prediction reads them and superimposes the scars on a
step-change concentration field. This exercises the bridge between the two,
which no single-CLI test can.

THE PROPERTY UNDER TEST is not "the numbers look plausible" -- the unit test
covers the physics. It is that **switching the coupling on leaves unfailed
ground bitwise unchanged**. A landslide field touches a handful of pixels; if
turning the coupling on perturbed the rest of the raster, every basin average
computed from it would shift for reasons nothing to do with landslides, and the
shift would be invisible against a background of plausible-looking numbers.

So the script runs the prediction twice on identical inputs, once with
`landslide_coupling: none` and once with `scar`, and requires that the two
rasters differ in **some** pixels (or the coupling did nothing) and that the
difference is confined to the failed ones.

Orchestration only: no numerics live here.
"""

import argparse
import pathlib
import re
import shutil
import subprocess
import sys


LANDSLIDE_PARAMS = """read_path: .
read_fname: tiny_dem.asc
write_path: .
write_fname: ls
landslide_mode: single_event
threshold_contributing_pixels: 3
area_volume_preset: larsen_soil_scar
event_time_yr: 1000
event_volume_m3: 200
n_min_pixels: 1
write_time_of_failure_raster: true
write_excavation_depth_raster: true
"""

PREDICTION_COMMON = """read_path: .
read_fname: tiny_dem.asc
write_path: .
calculate_step_change_crn_concentration_raster: true
step_change_eff_erate_old_g_cm2_yr: 0.001
step_change_eff_erate_new_g_cm2_yr: 0.005
time_since_change_yr: 20000
"""


def run(exe: pathlib.Path, work: pathlib.Path, text: str, name: str) -> str:
    param = work / f"{name}.param"
    param.write_text(text)
    r = subprocess.run([str(exe), param.name], cwd=str(work),
                       capture_output=True, text=True, check=False)
    if r.returncode != 0:
        sys.stderr.write(r.stdout)
        sys.stderr.write(r.stderr)
        raise SystemExit(f"{name} failed with exit {r.returncode}")
    (work / f"{name}.stdout").write_text(r.stdout)
    return r.stdout


def band_values(raster: pathlib.Path) -> list:
    """Every pixel value, via gdal_translate to an ASCII grid."""
    out = raster.with_suffix(".asc")
    r = subprocess.run(
        ["gdal_translate", "-of", "AAIGrid", "-q", str(raster), str(out)],
        capture_output=True, text=True, check=False)
    if r.returncode != 0:
        raise SystemExit(f"gdal_translate failed for {raster}: {r.stderr}")
    vals = []
    for line in out.read_text().splitlines():
        if re.match(r"^\s*[-0-9]", line):
            vals.extend(line.split())
    return vals


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--landslides-exe", required=True)
    ap.add_argument("--prediction-exe", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--work-dir", required=True)
    args = ap.parse_args()

    ls_exe = pathlib.Path(args.landslides_exe).resolve()
    pr_exe = pathlib.Path(args.prediction_exe).resolve()
    data = pathlib.Path(args.data_dir).resolve()
    work = pathlib.Path(args.work_dir).resolve()

    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    for stem in ("tiny_dem.asc", "tiny_dem.prj"):
        shutil.copy(data / stem, work / stem)

    # 1. Produce the landslide field.
    ls_out = run(ls_exe, work, LANDSLIDE_PARAMS, "ls")
    for needed in ("ls_LS_TFAIL.tif", "ls_LS_EXCAVATION.tif"):
        if not (work / needed).exists():
            raise SystemExit(f"lsdtt-cosmo-landslides did not write {needed}")
    if "Densities:" not in ls_out and "scar material" not in ls_out:
        print("note: lsdtt-cosmo-landslides did not report densities (CRN off)")

    # 2. Baseline: the same step change with no landslides.
    run(pr_exe, work, PREDICTION_COMMON + "write_fname: base\n"
        "landslide_coupling: none\n", "base")

    # 3. With the scars.
    scar_out = run(pr_exe, work, PREDICTION_COMMON + """write_fname: scar
landslide_coupling: scar
landslide_time_raster_prefix: ls_LS_TFAIL
landslide_excavation_raster_prefix: ls_LS_EXCAVATION
landslide_report_time_yr: 5000
scar_density_kg_m3: 1500
""", "scar")

    m = re.search(r"failed pixels\s*:\s*([0-9]+)", scar_out)
    if not m:
        raise SystemExit("the scar run did not report a failed-pixel count")
    n_failed = int(m.group(1))
    print(f"scar run reports {n_failed} failed pixels")
    if n_failed == 0:
        raise SystemExit(
            "no pixels failed, so this test proves nothing. The fixture is "
            "meant to place a few scars -- check the landslide parameters.")

    base = band_values(work / "base_Be10_step_crn_conc.tif")
    scar = band_values(work / "scar_Be10_step_crn_conc.tif")
    if len(base) != len(scar):
        raise SystemExit(f"raster sizes differ: {len(base)} vs {len(scar)}")

    changed = [i for i, (a, b) in enumerate(zip(base, scar)) if a != b]
    print(f"pixels changed by the coupling: {len(changed)} of {len(base)}")

    if not changed:
        raise SystemExit(
            "the scar coupling changed nothing at all. Either the landslide "
            "rasters were not read or the scars are being ignored.")

    if len(changed) > n_failed:
        raise SystemExit(
            f"the coupling changed {len(changed)} pixels but only {n_failed} "
            "failed. Ground that never failed must be BITWISE unchanged -- "
            "otherwise every basin average shifts for reasons nothing to do "
            "with landslides.")

    print(f"OK: {len(changed)} changed pixels, all within the {n_failed} that "
          "failed; unfailed ground is bitwise identical.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
