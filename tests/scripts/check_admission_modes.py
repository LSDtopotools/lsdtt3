#!/usr/bin/env python3
"""Assert that relief and gradient admission select genuinely different ground.

WHY THIS TEST EXISTS. `admission_algorithm_choice` exists because local relief
and a gradient bound are not interchangeable: relief charges a surface for its
own tilt, so a smooth sloping tread fails it while a level but rough surface
passes, and a gradient bound does the reverse. If a refactor ever collapsed the
two tests into one - or wired `relief_or_gradient` to the relief mask alone -
every count below would still be plausible and nothing else in the suite would
notice.

So the DEM is built to make the two tests disagree by construction:

  A  a SMOOTH RAMP        - 8 percent gradient, no roughness. Passes a 0.268
                            (15 degree) gradient bound. Over a 5-cell window it accumulates 8 m of
                            relief, so it FAILS a 4 m relief threshold.
  B  a ROUGH FLAT         - zero mean gradient with +/- 3 m of noise. Passes
                            nothing on gradient once the noise is fitted, and
                            FAILS the relief threshold too. Present so that
                            "admits everything flat-ish" does not pass.
  C  a SMOOTH FLAT        - flat and quiet. Passes BOTH.
  D  a STEEP ROUGH SLOPE  - gradient 1 (45 degrees). Passes NEITHER.

The assertions are on the SET RELATIONS, not on cell counts, because counts move
with the morphological post-processing and the relations do not:

    and  <=  min(relief, gradient)  <=  max(relief, gradient)  <=  or
    or   >   relief                     (gradient must ADD ground)
    or   >   gradient                   (relief must ADD ground)
    and  <   min(relief, gradient)      (the two must genuinely disagree)

It also runs the pre-2026-10-26 keys (admission_mode with `slope`,
admission_slope_cap_degrees, relief_window_side_pixels,
adaptive_relief_threshold, relief_threshold_metres) and requires the output to
be byte-identical to the canonical keys: the renames must carry every value.

Uses ReadRaster + struct rather than gdal.Array: the GDAL numpy bindings have an
ABI break in this environment, which is also why the other scripts here avoid
them.
"""

from __future__ import annotations

import argparse
import math
import struct
import subprocess
import sys
from pathlib import Path
import lsdtt_test_env


def build_dem(path: Path, cell: float = 10.0) -> None:
    """Write the four-quadrant synthetic DEM described in the module docstring.

    Args:
        path: GeoTIFF to create.
        cell: Cell size in metres.
    """
    from osgeo import gdal, osr

    n = 120
    half = n // 2
    rng_state = 12345

    def rnd() -> float:
        """A tiny reproducible LCG, so the test does not depend on numpy."""
        nonlocal rng_state
        rng_state = (1103515245 * rng_state + 12345) % (1 << 31)
        return rng_state / float(1 << 31) - 0.5

    rows = []
    for r in range(n):
        row = []
        for c in range(n):
            if r < half and c < half:            # A: smooth 8% ramp
                z = 100.0 + 0.08 * c * cell
            elif r < half:                        # B: rough flat
                z = 300.0 + 6.0 * rnd()
            elif c < half:                        # C: smooth flat
                z = 500.0
            else:                                 # D: steep rough slope
                z = 700.0 + 1.0 * r * cell + 2.0 * rnd()
            row.append(z)
        rows.append(row)

    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(str(path), n, n, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((0.0, cell, 0.0, n * cell, 0.0, -cell))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32644)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(-9999.0)
    for r, row in enumerate(rows):
        band.WriteRaster(0, r, n, 1, struct.pack(f"<{n}f", *row),
                         buf_type=gdal.GDT_Float32)
    band.FlushCache()
    ds = None


def count_set(path: Path) -> int:
    """Count cells above 0.5 in a mask raster.

    Args:
        path: Mask GeoTIFF.

    Returns:
        Number of set cells.
    """
    from osgeo import gdal

    ds = gdal.Open(str(path))
    if ds is None:
        raise RuntimeError(f"cannot open {path}")
    w, h = ds.RasterXSize, ds.RasterYSize
    raw = ds.GetRasterBand(1).ReadRaster(0, 0, w, h, buf_type=gdal.GDT_Float32)
    vals = struct.unpack(f"<{w * h}f", raw)
    ds = None
    return sum(1 for v in vals if v > 0.5)


CANONICAL_KEYS = (
    "relief_window_radius_metres: 20\n"
    "relief_threshold_algorithm_choice: fixed\n"
    "max_relief_metres: 4.0\n"
    "admission_algorithm_choice: {mode}\n"
    "max_admission_gradient: {gradient}\n")

# The same settings under the names used before 2026-10-26, with the angle the
# canonical gradient is the tangent of.
LEGACY_KEYS = (
    "relief_window_side_pixels: 5\n"
    "adaptive_relief_threshold: false\n"
    "relief_threshold_metres: 4.0\n"
    "admission_mode: {mode}\n"
    "admission_slope_cap_degrees: {degrees}\n")


def run_mode(exe: str, work: Path, dem: Path, mode: str, gradient: float,
             stem: str | None = None,
             window_setting: str =
             "admission_gradient_window_radius_metres: 25.0",
             keys: str | None = None) -> tuple[int, str, str, bytes]:
    """Run the finder in one admission mode and count the combined flat mask.

    Args:
        exe: Path to lsdtt-flat-surface-finder.
        work: Working directory.
        dem: Input DEM.
        mode: admission_algorithm_choice value.
        gradient: max_admission_gradient (m/m).
        stem: Output stem; defaults to one per mode.
        window_setting: The gradient-window line.
        keys: The admission lines; defaults to CANONICAL_KEYS.

    Returns:
        Set-cell count, stdout, stderr, and output bytes.
    """
    stem = stem or f"adm_{mode}"
    param = work / f"{stem}.param"
    if keys is None:
        keys = CANONICAL_KEYS.format(mode=mode, gradient=f"{gradient:.9g}")
    param.write_text(
        f"read_path: {work}/\n"
        f"read_fname: {dem.name}\n"
        f"write_path: {work}/\n"
        f"write_fname: {stem}\n"
        "n_threads: 1\n"
        "write_output_format: gtiff\n"
        "min_component_pixels: 1\n"
        "flat_closing_radius_metres: 0\n"
        "write_combined_flat_mask_raster: true\n"
        "write_healed_flat_mask: false\n"
        f"{keys}"
        f"{window_setting}\n",
        encoding="utf-8")
    proc = lsdtt_test_env.run([exe, param.name], capture_output=True,
                          text=True, cwd=str(work))
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout + proc.stderr)
        raise RuntimeError(f"admission {mode} exited {proc.returncode}")
    output = work / f"{stem}_combined_flat_mask.tif"
    return count_set(output), proc.stdout, proc.stderr, output.read_bytes()


def main(argv: list[str] | None = None) -> int:
    """Build the DEM, run all four modes, assert the set relations.

    Args:
        argv: Optional argument list.

    Returns:
        0 on success, 1 on a failed assertion.
    """
    lsdtt_test_env.require_modules(["osgeo.gdal"])
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--cap-degrees", type=float, default=15.0,
                    help="Gradient bound as an angle; the canonical runs use "
                         "its tangent, the legacy run the angle itself.")
    args = ap.parse_args(argv)

    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    dem = work / "admission_dem.tif"
    build_dem(dem)

    # The float the CLI holds, so the canonical text reads back to the same
    # value the legacy angle converts to.
    gradient = struct.unpack(
        "<f", struct.pack("<f", math.tan(math.radians(args.cap_degrees))))[0]
    n = {m: run_mode(args.exe, work, dem, m, gradient)[0]
          for m in ("relief", "gradient", "relief_and_gradient",
                    "relief_or_gradient")}
    for m, v in n.items():
        print(f"  {m:18s} {v:8d}")

    lo = min(n["relief"], n["gradient"])
    hi = max(n["relief"], n["gradient"])
    canonical = run_mode(args.exe, work, dem, "gradient", gradient,
                         stem="adm_window_equivalence")
    legacy = run_mode(
        args.exe, work, dem, "gradient", gradient,
        stem="adm_window_equivalence",
        window_setting="admission_slope_window_metres: 50.0")
    # Every pre-2026-10-26 admission key at once, including the old mode
    # spelling and the angle, against the canonical run of the same mode.
    old_names = run_mode(
        args.exe, work, dem, "relief_or_slope", gradient,
        stem="adm_old_names",
        window_setting="admission_slope_window_radius_metres: 25.0",
        keys=LEGACY_KEYS.format(mode="relief_or_slope",
                                degrees=f"{args.cap_degrees:g}"))
    new_names = run_mode(args.exe, work, dem, "relief_or_gradient", gradient,
                         stem="adm_new_names")

    checks = [
        ("and <= min(relief, gradient)", n["relief_and_gradient"] <= lo),
        ("max(relief, gradient) <= or",  hi <= n["relief_or_gradient"]),
        ("or > relief (gradient adds)",  n["relief_or_gradient"] > n["relief"]),
        ("or > gradient (relief adds)",  n["relief_or_gradient"] > n["gradient"]),
        ("and < min (they disagree)",    n["relief_and_gradient"] < lo),
        ("relief admits something",      n["relief"] > 0),
        ("gradient admits something",    n["gradient"] > 0),
        ("diameter alias keeps 50 m compute window",
         "polyfit window 50 m" in canonical[1] and
         "polyfit window 50 m" in legacy[1]),
        ("diameter alias output is byte-identical", canonical[3] == legacy[3]),
        ("diameter alias warning pins transformation",
         "'admission_slope_window_metres' is deprecated; use "
         "'admission_gradient_window_radius_metres'" in legacy[2] and
         "new = old * 0.5 + 0; zero remains zero" in legacy[2]),
        ("old admission keys give byte-identical output",
         old_names[3] == new_names[3]),
        ("old angle converted to the same gradient",
         f"max_admission_gradient {gradient:.9g}" in old_names[1]),
        ("old keys warn, naming the new keys",
         all(w in old_names[2] for w in (
             "'relief_window_side_pixels' is deprecated; use "
             "'relief_window_radius_metres'",
             "'adaptive_relief_threshold' is deprecated; use "
             "'relief_threshold_algorithm_choice'",
             "'relief_threshold_metres' is deprecated; use 'max_relief_metres'",
             "'admission_mode' is deprecated; use 'admission_algorithm_choice'",
             "'admission_slope_cap_degrees' is deprecated; use "
             "'max_admission_gradient'",
             "admission_algorithm_choice 'relief_or_slope' is deprecated"))),
    ]
    bad = [name for name, ok in checks if not ok]
    for name, ok in checks:
        print(f"  [{'ok' if ok else 'FAIL'}] {name}")
    if bad:
        print(f"FAILED: {len(bad)} check(s)")
        return 1
    print("all admission-mode checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
