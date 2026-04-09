#!/usr/bin/env python3
"""
tif_to_png.py  --  Convert GeoTIFF rasters to PNG images for visual inspection.

Uses gdal_translate and gdalinfo from the GDAL command-line tools.
No Python GDAL bindings or numpy required.

Output images are 1500 pixels wide by default (5 inches at 300 DPI).
Height is calculated to preserve the raster's aspect ratio.
Pixel values are auto-stretched to the full data range with -scale.

Usage examples:
  python tif_to_png.py /path/to/outputs
  python tif_to_png.py /path/to/outputs --file my_dem_mountain_mask.tif
  python tif_to_png.py . --width 5 --dpi 150
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

DEFAULT_WIDTH_INCHES = 5.0
DEFAULT_DPI          = 300


def get_raster_size(tif_path):
    """Return (n_cols, n_rows) by calling gdalinfo -json."""
    result = subprocess.run(
        ["gdalinfo", "-json", str(tif_path)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"gdalinfo failed: {result.stderr.strip()}")

    info = json.loads(result.stdout)
    size = info.get("size")
    if size is None or len(size) < 2:
        raise RuntimeError(f"gdalinfo did not return a size field for {tif_path.name}")
    return int(size[0]), int(size[1])   # (cols, rows)


def tif_to_png(tif_path, width_inches=DEFAULT_WIDTH_INCHES, dpi=DEFAULT_DPI):
    """
    Convert a single GeoTIFF to PNG using gdal_translate.

    Output pixel width  = round(width_inches * dpi)
    Output pixel height = preserves aspect ratio of the source raster.

    Returns the Path of the written PNG.
    """
    tif_path = Path(tif_path)
    png_path = tif_path.with_suffix(".png")

    out_width_px = round(width_inches * dpi)

    n_cols, n_rows = get_raster_size(tif_path)
    if n_cols == 0 or n_rows == 0:
        raise RuntimeError(f"{tif_path.name}: raster has zero-size dimension")

    out_height_px = round(out_width_px * n_rows / n_cols)

    cmd = [
        "gdal_translate",
        "-of", "PNG",
        "-outsize", str(out_width_px), str(out_height_px),
        "-scale",
        "-ot", "Byte",
        str(tif_path),
        str(png_path),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"gdal_translate failed: {result.stderr.strip()}")

    return png_path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert GeoTIFF rasters to PNG images using gdal_translate.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "directory",
        nargs="?",
        default=".",
        help="Directory to search for .tif files (default: current directory).",
    )
    parser.add_argument(
        "--file", "-f",
        default=None,
        metavar="FILENAME",
        help="Process only this single .tif file. "
             "Can be a filename inside the directory argument "
             "or an absolute/relative path.",
    )
    parser.add_argument(
        "--width",
        type=float,
        default=DEFAULT_WIDTH_INCHES,
        metavar="INCHES",
        help=f"Output image width in inches (default: {DEFAULT_WIDTH_INCHES}). "
             "Height is calculated to preserve the raster aspect ratio.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=DEFAULT_DPI,
        metavar="DPI",
        help=f"Pixels per inch used to compute output pixel width "
             f"(default: {DEFAULT_DPI}). "
             "Output pixel width = round(width * dpi).",
    )
    return parser.parse_args()


def main():
    args    = parse_args()
    workdir = Path(args.directory).resolve()

    if not workdir.is_dir():
        print(f"Error: directory not found: {workdir}", file=sys.stderr)
        sys.exit(1)

    if args.file is not None:
        candidate = Path(args.file)
        if not candidate.is_absolute():
            candidate = workdir / candidate
        if not candidate.exists():
            print(f"Error: file not found: {candidate}", file=sys.stderr)
            sys.exit(1)
        tif_files = [candidate]
    else:
        tif_files = sorted(workdir.glob("*.tif"))
        if not tif_files:
            print(f"No .tif files found in {workdir}")
            sys.exit(0)

    out_px = round(args.width * args.dpi)
    print(f"Converting {len(tif_files)} file(s)  "
          f"({args.width} in × {args.dpi} DPI = {out_px} px wide)\n")

    n_ok  = 0
    n_err = 0
    for tif in tif_files:
        try:
            png = tif_to_png(tif, width_inches=args.width, dpi=args.dpi)
            print(f"  OK  {tif.name}  ->  {png.name}")
            n_ok += 1
        except Exception as exc:
            print(f"  ERR {tif.name}  :  {exc}", file=sys.stderr)
            n_err += 1

    print(f"\nDone: {n_ok} converted, {n_err} failed.")
    if n_err:
        sys.exit(1)


if __name__ == "__main__":
    main()
