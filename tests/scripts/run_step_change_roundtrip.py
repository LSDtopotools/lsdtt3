#!/usr/bin/env python3
"""Disk round-trip regression for the generated step-change rasters.

The generated time-since-change field is handed to the concentration kernel in
memory. This asserts that doing so is equivalent to the older workflow of
writing the field to disk and reading it back in a second run: the two
step-change concentration rasters must be identical.

Why it is worth a test. The in-memory path and the on-disk path differ in one
respect that is easy to get wrong -- the on-disk field has been through a
float32 GeoTIFF and a nodata mask. If the generator ever started producing a
field the writer could not represent (a value at the nodata sentinel, a mask
that does not survive the round trip), concentrations would silently diverge
between two workflows that are supposed to agree.

Both runs compute production from the INPUT DEM, not the generated one, so that
the only difference between them is how the time-since-change field arrived.

WHY THE COMPARISON IS NOT BIT-EXACT, AND WHY THAT IS NOT A WEAKENING.
The two runs cannot use bit-identical effective erosion rates, and no parameter
file can make them. Run A forms the rate as a product of doubles promoted from
floats -- (double)0.0001f * (double)2650.0f * 0.1 = 0.026499999330553692 --
while run B parses the literal 0.0265 into a float and promotes it, giving
0.026499999687075615. They differ by 1.3e-8 relative, which is the float32
representation of a decimal in a text file and nothing to do with the raster.
So the band CHECKSUM is asserted exactly (it is insensitive to that last bit)
and the statistics to 1e-6 relative. A genuine round-trip failure -- a lost
mask, a value colliding with the nodata sentinel, a geometry mismatch -- moves
concentrations by orders of magnitude more than that and cannot hide inside the
tolerance.

Orchestration only: no numerics live here. Comparison is by GDAL band checksum
plus min/max/mean, which is what the CLI smoke harness already uses.
"""

import argparse
import pathlib
import re
import shutil
import subprocess
import sys


def run_cli(exe: pathlib.Path, work: pathlib.Path, param_text: str,
            name: str) -> None:
    param_path = work / f"{name}.param"
    param_path.write_text(param_text)
    result = subprocess.run(
        [str(exe), str(param_path.name)],
        cwd=str(work),
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise SystemExit(f"{name} run failed with exit {result.returncode}")
    (work / f"{name}.stdout").write_text(result.stdout)


def raster_signature(raster: pathlib.Path) -> dict:
    """Band checksum plus statistics, as a comparable dict."""
    result = subprocess.run(
        ["gdalinfo", "-checksum", "-stats", str(raster)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise SystemExit(f"gdalinfo failed for {raster}: {result.stderr}")

    sig = {}
    m = re.search(r"Checksum=([0-9]+)", result.stdout)
    if not m:
        raise SystemExit(f"no checksum for {raster}")
    sig["checksum"] = m.group(1)
    for key in ("STATISTICS_MINIMUM", "STATISTICS_MAXIMUM", "STATISTICS_MEAN"):
        m = re.search(rf"{key}=([-0-9.eE+]+)", result.stdout)
        if not m:
            raise SystemExit(f"no {key} for {raster}")
        sig[key] = m.group(1)
    return sig


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exe", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--work-dir", required=True)
    args = ap.parse_args()

    exe = pathlib.Path(args.exe).resolve()
    data = pathlib.Path(args.data_dir).resolve()
    work = pathlib.Path(args.work_dir).resolve()

    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    for stem in ("tiny_dem.asc", "tiny_dem.prj"):
        shutil.copy(data / stem, work / stem)

    common = """read_path: .
read_fname: tiny_dem.asc
write_path: .
step_change_channel_pixels: 2
calculate_step_change_crn_concentration_raster: true
"""

    # Run A: generate, write the fields, and consume the generated one in
    # memory. Production from the INPUT DEM so the two runs match.
    run_cli(exe, work, common + """write_fname: rtA
make_step_change_rasters: true
step_change_uplift_old_m_yr: 0.0001
step_change_uplift_new_m_yr: 0.0005
step_change_target_relief_metres: 400
step_change_age_yr: 1000000
step_change_n: 1.0
write_step_change_rasters: true
step_change_production_from_generated_dem: false
""", "rtA")

    steptime = work / "rtA_SC_STEPTIME.tif"
    if not steptime.exists():
        raise SystemExit(f"run A did not write {steptime}")

    # Run B: load the field run A wrote, with the same effective erosion rates.
    # 0.0001 m/yr * 2650 kg/m3 * 0.1 = 0.0265 g/cm2/yr, and 5x that.
    run_cli(exe, work, common + """write_fname: rtB
make_step_change_rasters: false
time_since_change_raster_prefix: rtA_SC_STEPTIME
step_change_eff_erate_old_g_cm2_yr: 0.0265
step_change_eff_erate_new_g_cm2_yr: 0.1325
""", "rtB")

    a = work / "rtA_Be10_step_crn_conc.tif"
    b = work / "rtB_Be10_step_crn_conc.tif"
    for p in (a, b):
        if not p.exists():
            raise SystemExit(f"missing concentration raster {p}")

    sig_a = raster_signature(a)
    sig_b = raster_signature(b)

    print("in-memory :", sig_a)
    print("round-trip:", sig_b)

    if sig_a["checksum"] != sig_b["checksum"]:
        raise SystemExit(
            "step-change concentration CHECKSUMS differ between the in-memory "
            "and on-disk routes:\n"
            f"  in-memory  {sig_a}\n"
            f"  round-trip {sig_b}\n"
            "The generated time-since-change field is not surviving a GeoTIFF "
            "round trip unchanged. This is NOT the float32 parameter effect "
            "documented at the top of this script -- that does not move the "
            "checksum."
        )

    tol = 1.0e-6
    worst = 0.0
    for key in ("STATISTICS_MINIMUM", "STATISTICS_MAXIMUM", "STATISTICS_MEAN"):
        va, vb = float(sig_a[key]), float(sig_b[key])
        denom = abs(vb) if vb != 0.0 else 1.0
        rel = abs(va - vb) / denom
        worst = max(worst, rel)
        if rel > tol:
            raise SystemExit(
                f"{key} differs by {rel:.3e} (> {tol:.0e}) between the "
                f"in-memory and on-disk routes: {va} vs {vb}.\n"
                "Too large to be the float32 parameter representation; the "
                "time-since-change field is not round-tripping."
            )

    print(f"OK: identical checksum, statistics agree to {worst:.3e} relative "
          f"(tolerance {tol:.0e}; expected residual ~1e-7 from the float32 "
          f"parameter representation, not the raster).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
