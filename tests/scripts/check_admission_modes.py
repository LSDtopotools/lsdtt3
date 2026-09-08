#!/usr/bin/env python3
"""Assert that relief and slope admission select genuinely different ground.

WHY THIS TEST EXISTS. `admission_mode` was added because local relief and a
gradient cap are not interchangeable: relief charges a surface for its own tilt,
so a smooth sloping tread fails it while a level but rough surface passes, and a
gradient cap does the reverse. If a refactor ever collapsed the two tests into
one - or wired `relief_or_slope` to the relief mask alone - every count below
would still be plausible and nothing else in the suite would notice.

So the DEM is built to make the two tests disagree by construction:

  A  a SMOOTH RAMP        - 8 percent gradient, no roughness. Passes a 15 degree
                            cap. Over a 5-cell window it accumulates 8 m of
                            relief, so it FAILS a 4 m relief threshold.
  B  a ROUGH FLAT         - zero mean gradient with +/- 3 m of noise. Passes
                            nothing on gradient once the noise is fitted, and
                            FAILS the relief threshold too. Present so that
                            "admits everything flat-ish" does not pass.
  C  a SMOOTH FLAT        - flat and quiet. Passes BOTH.
  D  a STEEP ROUGH SLOPE  - 45 degrees. Passes NEITHER.

The assertions are on the SET RELATIONS, not on cell counts, because counts move
with the morphological post-processing and the relations do not:

    and  <=  min(relief, slope)  <=  max(relief, slope)  <=  or
    or   >   relief                  (slope must ADD ground)
    or   >   slope                   (relief must ADD ground)
    and  <   min(relief, slope)      (the two must genuinely disagree)

Uses ReadRaster + struct rather than gdal.Array: the GDAL numpy bindings have an
ABI break in this environment, which is also why the other scripts here avoid
them.
"""

from __future__ import annotations

import argparse
import struct
import subprocess
import sys
from pathlib import Path


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


def run_mode(exe: str, work: Path, dem: Path, mode: str, cap: float) -> int:
    """Run the finder in one admission mode and count the combined flat mask.

    Args:
        exe: Path to lsdtt-flat-surface-finder.
        work: Working directory.
        dem: Input DEM.
        mode: admission_mode value.
        cap: admission_slope_cap_degrees.

    Returns:
        Set-cell count of the combined flat mask.
    """
    stem = f"adm_{mode}"
    param = work / f"{stem}.param"
    param.write_text(
        f"read_path: {work}/\n"
        f"read_fname: {dem.name}\n"
        f"write_path: {work}/\n"
        f"write_fname: {stem}\n"
        "n_threads: 1\n"
        "write_output_format: gtiff\n"
        "relief_window_pixels: 5\n"
        "adaptive_relief_threshold: false\n"
        "relief_threshold_metres: 4.0\n"
        "min_component_pixels: 1\n"
        "flat_closing_radius_metres: 0\n"
        "write_combined_flat_mask_raster: true\n"
        "write_slope_break_centreline: false\n"
        "write_slope_break_centreline_connected: false\n"
        "write_slope_break_bridge_mask: false\n"
        "write_healed_flat_mask: false\n"
        f"admission_mode: {mode}\n"
        f"admission_slope_cap_degrees: {cap}\n"
        "admission_slope_window_metres: 50.0\n",
        encoding="utf-8")
    proc = subprocess.run([exe, param.name], capture_output=True,
                          text=True, cwd=str(work))
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout + proc.stderr)
        raise RuntimeError(f"admission_mode {mode} exited {proc.returncode}")
    return count_set(work / f"{stem}_combined_flat_mask.tif")


def main(argv: list[str] | None = None) -> int:
    """Build the DEM, run all four modes, assert the set relations.

    Args:
        argv: Optional argument list.

    Returns:
        0 on success, 1 on a failed assertion.
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--cap", type=float, default=15.0)
    args = ap.parse_args(argv)

    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    dem = work / "admission_dem.tif"
    build_dem(dem)

    n = {m: run_mode(args.exe, work, dem, m, args.cap)
         for m in ("relief", "slope", "relief_and_slope", "relief_or_slope")}
    for m, v in n.items():
        print(f"  {m:18s} {v:8d}")

    lo = min(n["relief"], n["slope"])
    hi = max(n["relief"], n["slope"])
    checks = [
        ("and <= min(relief, slope)",   n["relief_and_slope"] <= lo),
        ("max(relief, slope) <= or",    hi <= n["relief_or_slope"]),
        ("or > relief (slope adds)",    n["relief_or_slope"] > n["relief"]),
        ("or > slope (relief adds)",    n["relief_or_slope"] > n["slope"]),
        ("and < min (they disagree)",   n["relief_and_slope"] < lo),
        ("relief admits something",     n["relief"] > 0),
        ("slope admits something",      n["slope"] > 0),
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
