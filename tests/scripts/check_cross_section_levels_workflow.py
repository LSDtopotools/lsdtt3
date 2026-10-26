#!/usr/bin/env python3
"""Run the three-stage cross-section levels workflow on the staircase and check it.

Fluvial-surface methods plan, phase 5. The levels of a terrace flight can come
from cross sections instead of from the ladder:

  1. lsdtt-fluvial-surfaces       -> a valley mask
  2. lsdtt-valley-centreline      -> every section's width-step breaks, pooled
                                     per channel into <stem>_cross_section_levels.csv
  3. lsdtt-fluvial-surfaces       -> the admitted ground banded on those levels
     (fluvial_surface_algorithm_choice cross_section_levels)

The staircase fixture's treads sit 8, 16, 24 and 32 m above the channel and its
tread map is known, so both ends can be checked against the truth:

  * stage 2 finds four levels within 1 m of 8, 16, 24 and 32 m below 40 m, and
    reports the steps it dropped for lack of support (sections crossing the
    riser-gap ramps);
  * stage 3 runs every ladder window on the supplied levels, none below the
    resolvability floor, and each band 0..3 is at least 95 % one tread - the
    tread it should be.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

import lsdtt_test_env

TREADS_M = (8.0, 16.0, 24.0, 32.0)


def run_stage(exe: str, template: pathlib.Path, read_path: pathlib.Path,
              out_dir: pathlib.Path, work_dir: pathlib.Path, tag: str) -> str:
    """Render a template, run the CLI on it, and return its stdout."""
    text = (template.read_text(encoding="utf-8")
            .replace("__READ_PATH__", str(read_path))
            .replace("__WRITE_PATH__", str(out_dir)))
    param = work_dir / f"{tag}.param"
    param.write_text(text, encoding="utf-8")
    proc = lsdtt_test_env.run([exe, str(work_dir), param.name],
                              cwd=str(work_dir), capture_output=True,
                              text=True, check=False)
    (work_dir / f"stdout_{tag}.log").write_text(proc.stdout, encoding="utf-8")
    (work_dir / f"stderr_{tag}.log").write_text(proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout[-3000:] + proc.stderr[-3000:])
        raise SystemExit(f"stage {tag} exited {proc.returncode}")
    return proc.stdout


def read_float_raster(path: pathlib.Path):
    """Band 1 as a flat list of floats and its nodata value.

    Read with ReadRaster rather than ReadAsArray: gdal_array is tied to the
    NumPy it was built against, and that is not guaranteed to be the NumPy in
    the test interpreter (the same choice as check_admission_modes.py).
    """
    import struct
    from osgeo import gdal
    gdal.UseExceptions()
    ds = gdal.Open(str(path))
    w, h = ds.RasterXSize, ds.RasterYSize
    band = ds.GetRasterBand(1)
    raw = band.ReadRaster(0, 0, w, h, buf_type=gdal.GDT_Float32)
    return list(struct.unpack(f"<{w * h}f", raw)), band.GetNoDataValue()


def band_purity(band_path: pathlib.Path, tread_path: pathlib.Path):
    """Per band: (cells, majority tread id, purity)."""
    from collections import Counter, defaultdict
    band, nodata = read_float_raster(band_path)
    tread, _ = read_float_raster(tread_path)
    per_band = defaultdict(Counter)
    for b, t in zip(band, tread):
        if b < 0 or (nodata is not None and b == nodata):
            continue
        per_band[int(b)][int(t)] += 1
    out = {}
    for b, counts in per_band.items():
        tread_id, n_major = counts.most_common(1)[0]
        n = sum(counts.values())
        out[b] = (n, tread_id, n_major / n)
    return out


def main() -> int:
    lsdtt_test_env.require_modules(["osgeo.gdal"])
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--surfaces-exe", required=True)
    ap.add_argument("--centreline-exe", required=True)
    ap.add_argument("--param-dir", required=True)
    ap.add_argument("--read-path", required=True)
    ap.add_argument("--work-dir", required=True)
    args = ap.parse_args()

    params = pathlib.Path(args.param_dir)
    read_path = pathlib.Path(args.read_path)
    work = pathlib.Path(args.work_dir)
    out = work / "out"
    out.mkdir(parents=True, exist_ok=True)

    run_stage(args.surfaces_exe, params / "fluvial_surfaces_xs_stage1.param.in",
              read_path, out, work, "stage1")
    log2 = run_stage(args.centreline_exe,
                     params / "valley_centreline_xs_levels.param.in",
                     read_path, out, work, "stage2")
    log3 = run_stage(args.surfaces_exe,
                     params / "fluvial_surfaces_xs_levels.param.in",
                     read_path, out, work, "stage3")

    failures = []

    def check(ok: bool, what: str) -> None:
        print(("  ok    " if ok else "  FAIL  ") + what)
        if not ok:
            failures.append(what)

    with open(out / "xs2_cross_section_levels.csv", newline="") as f:
        levels = [float(r["level_hand_m"]) for r in csv.DictReader(f)]
    low = sorted(v for v in levels if v < 40.0)
    print(f"  stage 2 levels: {levels}")
    check(len(low) == len(TREADS_M) and
          all(abs(a - b) < 1.0 for a, b in zip(low, TREADS_M)),
          "stage 2 finds the four treads within 1 m (below 40 m)")
    check("dropped" in log2,
          "stage 2 reports the steps it dropped for lack of support")

    with open(out / "xs3_ladder_windows.csv", newline="") as f:
        windows = list(csv.DictReader(f))
    check(bool(windows) and all(w["status"] == "supplied" for w in windows),
          "stage 3 runs every ladder window on the supplied levels")
    check("(0 below the resolvability floor)" in log3,
          "no window's supplied levels fall below the resolvability floor")

    purity = band_purity(out / "xs3_band_index.tif",
                         read_path / "staircase_tread_id.tif")
    for b, (n, tread, p) in sorted(purity.items()):
        print(f"  band {b}: {n} cells, majority tread {tread}, purity {p:.4f}")
    check(all(b in purity for b in range(4)), "stage 3 produces bands 0..3")
    check(all(purity[b][1] == b + 1 and purity[b][2] >= 0.95
              for b in range(4) if b in purity),
          "each band 0..3 is at least 95 % the tread it should be")

    if failures:
        print(f"FAILED: {len(failures)} check(s)")
        return 1
    print("all cross-section levels workflow checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
