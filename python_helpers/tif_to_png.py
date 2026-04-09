#!/usr/bin/env python3
"""
tif_to_png.py  --  Render GeoTIFF rasters as PNG images for visual inspection.

Uses GDAL for raster I/O, matplotlib for rendering.

Default behaviour: converts every .tif file found in a directory to a PNG
image saved alongside the original.  A single file can be processed instead
with --file.

Rendering:
  - Single-band rasters are displayed with a matplotlib colormap and a
    2 %/98 % percentile linear stretch.  Nodata pixels are shown white.
  - Three-band rasters are treated as RGB and displayed directly.
  - All output images are 5 inches wide at 300 DPI (1 500 px wide).
    Height is set to preserve the raster's aspect ratio.

Usage examples:
  python tif_to_png.py /path/to/outputs
  python tif_to_png.py /path/to/outputs --file my_dem_mountain_mask.tif
  python tif_to_png.py . --colormap grey --dpi 150
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")          # non-interactive backend — no display needed
import matplotlib.pyplot as plt
from osgeo import gdal

gdal.UseExceptions()           # raise Python exceptions instead of silent errors

# ─── defaults ────────────────────────────────────────────────────────────────
DEFAULT_WIDTH_INCHES = 5.0
DEFAULT_DPI          = 300
DEFAULT_COLORMAP     = "terrain"
PERCENTILE_LOW       = 2.0
PERCENTILE_HIGH      = 98.0


# ─── helpers ─────────────────────────────────────────────────────────────────

def _read_band(ds, band_index=1):
    """Return (array, nodata) for one GDAL band (1-based index)."""
    band   = ds.GetRasterBand(band_index)
    nodata = band.GetNoDataValue()
    arr    = band.ReadAsArray().astype(np.float64)
    return arr, nodata


def _mask_nodata(arr, nodata):
    """Return a masked array with nodata values masked out."""
    if nodata is not None:
        return np.ma.masked_where(arr == nodata, arr)
    return np.ma.masked_invalid(arr)


def _percentile_stretch(arr_masked):
    """Linear 2 %/98 % stretch → [0, 1] float array."""
    valid = arr_masked.compressed()
    if valid.size == 0:
        return arr_masked
    vmin = np.percentile(valid, PERCENTILE_LOW)
    vmax = np.percentile(valid, PERCENTILE_HIGH)
    if vmax == vmin:
        return np.zeros_like(arr_masked)
    stretched = (arr_masked - vmin) / (vmax - vmin)
    return np.clip(stretched, 0.0, 1.0)


def _build_rgba_single_band(arr_masked, colormap):
    """Convert a single-band masked array to an RGBA uint8 image."""
    stretched = _percentile_stretch(arr_masked)
    cmap      = plt.get_cmap(colormap)
    # cmap returns (R, G, B, A) floats in [0, 1]; mask maps to (1,1,1,1)=white
    rgba_f = cmap(stretched.filled(1.0))    # masked pixels → white in cmap
    # Make masked pixels fully transparent instead
    alpha        = (~arr_masked.mask).astype(np.float64) if np.ma.is_masked(stretched) \
                   else np.ones(stretched.shape)
    rgba_f[..., 3] = alpha
    return (rgba_f * 255).astype(np.uint8)


def _build_rgba_rgb(ds):
    """Build an RGBA image from a 3-band RGB GeoTIFF."""
    bands = []
    for i in range(1, 4):
        arr, nodata = _read_band(ds, i)
        bands.append(_mask_nodata(arr, nodata))

    # Stack and stretch each band independently
    rgb = np.stack([_percentile_stretch(b) for b in bands], axis=-1)
    # Alpha: valid where all three bands are valid
    mask  = np.any([np.ma.getmaskarray(b) for b in bands], axis=0)
    alpha = (~mask).astype(np.float64)

    rgba_f = np.dstack([rgb.filled(0.0), alpha])
    return (rgba_f * 255).astype(np.uint8)


# ─── main conversion function ────────────────────────────────────────────────

def tif_to_png(tif_path, colormap=DEFAULT_COLORMAP,
               width_inches=DEFAULT_WIDTH_INCHES, dpi=DEFAULT_DPI):
    """
    Convert a single GeoTIFF to a PNG image.

    Parameters
    ----------
    tif_path : Path or str
    colormap : str
        Matplotlib colormap name used for single-band rasters.
    width_inches : float
        Output image width in inches.
    dpi : int
        Output image resolution in dots per inch.

    Returns
    -------
    Path
        Path of the written PNG file.
    """
    tif_path = Path(tif_path)
    png_path = tif_path.with_suffix(".png")

    ds = gdal.Open(str(tif_path))
    if ds is None:
        raise RuntimeError(f"GDAL could not open: {tif_path}")

    n_bands = ds.RasterCount
    n_rows  = ds.RasterYSize
    n_cols  = ds.RasterXSize

    if n_cols == 0 or n_rows == 0:
        raise RuntimeError(f"{tif_path.name}: raster has zero-size dimension")

    aspect        = n_rows / n_cols
    height_inches = width_inches * aspect

    # Build RGBA image
    if n_bands >= 3:
        rgba = _build_rgba_rgb(ds)
    else:
        arr, nodata = _read_band(ds, 1)
        arr_masked  = _mask_nodata(arr, nodata)
        rgba        = _build_rgba_single_band(arr_masked, colormap)

    ds = None   # close GDAL dataset

    # Render with matplotlib
    fig, ax = plt.subplots(figsize=(width_inches, height_inches), dpi=dpi)
    ax.imshow(rgba, origin="upper", interpolation="nearest")
    ax.axis("off")
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)

    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", pad_inches=0)
    plt.close(fig)

    return png_path


# ─── CLI ─────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert GeoTIFF rasters to PNG images.",
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
        "--colormap", "-c",
        default=DEFAULT_COLORMAP,
        metavar="NAME",
        help=f"Matplotlib colormap for single-band rasters "
             f"(default: {DEFAULT_COLORMAP}). "
             "Examples: terrain, grey, viridis, plasma, RdBu.",
    )
    parser.add_argument(
        "--width",
        type=float,
        default=DEFAULT_WIDTH_INCHES,
        metavar="INCHES",
        help=f"Output image width in inches (default: {DEFAULT_WIDTH_INCHES}).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=DEFAULT_DPI,
        metavar="DPI",
        help=f"Output image resolution in DPI (default: {DEFAULT_DPI}).",
    )
    return parser.parse_args()


def main():
    args    = parse_args()
    workdir = Path(args.directory).resolve()

    if not workdir.is_dir():
        print(f"Error: directory not found: {workdir}", file=sys.stderr)
        sys.exit(1)

    # Resolve the list of files to process
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

    print(f"Converting {len(tif_files)} file(s) "
          f"({args.width} in wide, {args.dpi} DPI, colormap: {args.colormap})\n")

    n_ok  = 0
    n_err = 0
    for tif in tif_files:
        try:
            png = tif_to_png(tif,
                             colormap=args.colormap,
                             width_inches=args.width,
                             dpi=args.dpi)
            print(f"  OK  {tif.name}  →  {png.name}")
            n_ok += 1
        except Exception as exc:
            print(f"  ERR {tif.name}  :  {exc}", file=sys.stderr)
            n_err += 1

    print(f"\nDone: {n_ok} converted, {n_err} failed.")
    if n_err:
        sys.exit(1)


if __name__ == "__main__":
    main()
