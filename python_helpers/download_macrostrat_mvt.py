#!/usr/bin/env python3

import argparse
import math
import os
import sys
import urllib.error
import urllib.request

# Macrostrat public vector tile endpoint:
#   https://tiles.macrostrat.org/carto/{z}/{x}/{y}.mvt
#
# Corner input format:
#   "lat, lon"
#
# This is intended to match copy-paste from Google Maps, for example:
#   "55.965831207350035, -3.251400209050298"
#
# You should pass each corner in quotes because the shell will otherwise split
# on the space after the comma.
#
# Bounding box order:
#   1. lower-left  corner  (southwest)
#   2. upper-right corner  (northeast)
#
# Example:
#   python3 download_macrostrat_mvt.py \
#     "35.95753793852623, -121.61500241066403" \
#     "36.48163776466183, -121.21006900077077" \
#     10 macrostrat_tiles
#
# Zoom guidance:
#   0-5   very large area, low detail
#   6-8   regional scale
#   9-10  medium detail
#   11+   finer detail, many more tiles
#
# If unsure, start with 8, 9, or 10.

DEFAULT_TEMPLATE = "https://tiles.macrostrat.org/carto/{z}/{x}/{y}.mvt"


def clamp(value, lo, hi):
    return max(lo, min(hi, value))


def parse_corner(text, label):
    """
    Parse a corner in a Google Maps style format, e.g.

        "55.965831207350035, -3.251400209050298"

    Also tolerates surrounding parentheses:
        "(55.965831207350035, -3.251400209050298)"
    """
    s = text.strip()

    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1].strip()

    parts = s.split(",")
    if len(parts) != 2:
        raise ValueError(
            f'{label} must look like "lat, lon"; got: {text!r}'
        )

    lat_text = parts[0].strip()
    lon_text = parts[1].strip()

    try:
        lat = float(lat_text)
    except ValueError:
        raise ValueError(f"{label}: could not parse latitude from {lat_text!r}")

    try:
        lon = float(lon_text)
    except ValueError:
        raise ValueError(f"{label}: could not parse longitude from {lon_text!r}")

    if not (-90 <= lat <= 90):
        raise ValueError(f"{label}: latitude must be between -90 and 90")
    if not (-180 <= lon <= 180):
        raise ValueError(f"{label}: longitude must be between -180 and 180")

    return lat, lon


def lon_to_xtile(lon, z):
    n = 2 ** z
    lon = clamp(lon, -180.0, 180.0 - 1e-12)
    return int((lon + 180.0) / 360.0 * n)


def lat_to_ytile(lat, z):
    n = 2 ** z
    # Web Mercator is only defined to about +/-85.0511 degrees.
    lat = clamp(lat, -85.05112878, 85.05112878)
    lat_rad = math.radians(lat)
    y = (1.0 - math.asinh(math.tan(lat_rad)) / math.pi) / 2.0 * n
    return int(y)


def tile_range(lower_left_lat, lower_left_lon, upper_right_lat, upper_right_lon, z):
    if lower_left_lat >= upper_right_lat:
        raise ValueError("lower-left latitude must be less than upper-right latitude")
    if lower_left_lon >= upper_right_lon:
        raise ValueError("lower-left longitude must be less than upper-right longitude")

    # For slippy-map tiles:
    #   x increases west -> east
    #   y increases north -> south
    x_min = lon_to_xtile(lower_left_lon, z)
    x_max = lon_to_xtile(upper_right_lon - 1e-12, z)

    y_min = lat_to_ytile(upper_right_lat, z)
    y_max = lat_to_ytile(lower_left_lat + 1e-12, z)

    return x_min, x_max, y_min, y_max


def download_file(url, path, overwrite=False, timeout=60):
    if os.path.exists(path) and not overwrite:
        print(f"SKIP  {path}")
        return "skip"

    os.makedirs(os.path.dirname(path), exist_ok=True)

    req = urllib.request.Request(
        url,
        headers={"User-Agent": "macrostrat-mvt-downloader/1.0"}
    )

    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = resp.read()

        with open(path, "wb") as f:
            f.write(data)

        print(f"SAVED {path}")
        return "saved"

    except urllib.error.HTTPError as e:
        print(f"FAIL  {url}  HTTP {e.code}", file=sys.stderr)
        return "failed"
    except urllib.error.URLError as e:
        print(f"FAIL  {url}  {e}", file=sys.stderr)
        return "failed"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Download Macrostrat .mvt vector tiles for a bounding box.\n\n"
            "Corner input format:\n"
            '  "lat, lon"\n\n'
            "This is intended to match copy-paste from Google Maps.\n"
            "Example Google Maps corner text:\n"
            '  "55.965831207350035, -3.251400209050298"\n\n'
            "Pass the two corners in this order:\n"
            "  1. lower-left corner\n"
            "  2. upper-right corner\n\n"
            "Important:\n"
            "  Put each corner in quotes, otherwise the shell may split on the\n"
            "  space after the comma.\n\n"
            "Zoom guidance:\n"
            "  0-5   very large area, low detail\n"
            "  6-8   regional scale\n"
            "  9-10  medium detail\n"
            "  11+   much more detail, many more tiles\n\n"
            "If unsure, start with zoom 8, 9, or 10.\n\n"
            "Default tile template:\n"
            "  https://tiles.macrostrat.org/carto/{z}/{x}/{y}.mvt"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "lower_left",
        help='Lower-left corner as "lat, lon"'
    )
    parser.add_argument(
        "upper_right",
        help='Upper-right corner as "lat, lon"'
    )
    parser.add_argument(
        "zoom",
        type=int,
        help="Zoom level; start with 8, 9, or 10 if unsure"
    )
    parser.add_argument(
        "output_dir",
        help="Directory where downloaded .mvt files will be saved"
    )
    parser.add_argument(
        "--template",
        default=DEFAULT_TEMPLATE,
        help="Tile URL template (default: https://tiles.macrostrat.org/carto/{z}/{x}/{y}.mvt)"
    )
    parser.add_argument(
        "--flat",
        action="store_true",
        help="Save files directly in output_dir as z_x_y.mvt instead of z/x/y.mvt"
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing files"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print URLs and output paths without downloading"
    )

    args = parser.parse_args()

    if args.zoom < 0:
        parser.error("zoom must be >= 0")

    try:
        lower_left_lat, lower_left_lon = parse_corner(args.lower_left, "lower_left")
        upper_right_lat, upper_right_lon = parse_corner(args.upper_right, "upper_right")
        x_min, x_max, y_min, y_max = tile_range(
            lower_left_lat,
            lower_left_lon,
            upper_right_lat,
            upper_right_lon,
            args.zoom
        )
    except ValueError as e:
        parser.error(str(e))

    total = (x_max - x_min + 1) * (y_max - y_min + 1)

    print("Bounding box:")
    print(f"  lower-left : ({lower_left_lat}, {lower_left_lon})")
    print(f"  upper-right: ({upper_right_lat}, {upper_right_lon})")
    print(f"Zoom: {args.zoom}")
    print(f"Template: {args.template}")
    print("Tile range:")
    print(f"  x = {x_min}..{x_max}")
    print(f"  y = {y_min}..{y_max}")
    print(f"  total tiles = {total}")

    saved = 0
    skipped = 0
    failed = 0

    for x in range(x_min, x_max + 1):
        for y in range(y_min, y_max + 1):
            url = args.template.format(z=args.zoom, x=x, y=y)

            if args.flat:
                path = os.path.join(args.output_dir, f"{args.zoom}_{x}_{y}.mvt")
            else:
                path = os.path.join(args.output_dir, str(args.zoom), str(x), f"{y}.mvt")

            if args.dry_run:
                print(f"URL   {url}")
                print(f"PATH  {path}")
                continue

            result = download_file(url, path, overwrite=args.overwrite)

            if result == "saved":
                saved += 1
            elif result == "skip":
                skipped += 1
            else:
                failed += 1

    print()
    print(f"Done. Saved: {saved}, Skipped: {skipped}, Failed: {failed}")


if __name__ == "__main__":
    main()
