#!/usr/bin/env python3

from pathlib import Path
from urllib.parse import urlparse
import argparse
import requests
import planetary_computer
from pystac_client import Client

CATALOG_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"


def parse_latlon(text):
    """
    Parse a string like:
        "55.95244789111162, -3.2193971081922665"
    Returns (lat, lon)
    """
    parts = text.split(",")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            f"Invalid coordinate '{text}'. Expected format: 'lat, lon'"
        )

    try:
        lat = float(parts[0].strip())
        lon = float(parts[1].strip())
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid coordinate '{text}'. Could not parse numbers."
        )

    if not (-90 <= lat <= 90):
        raise argparse.ArgumentTypeError(f"Latitude out of range in '{text}'")
    if not (-180 <= lon <= 180):
        raise argparse.ArgumentTypeError(f"Longitude out of range in '{text}'")

    return lat, lon


def to_bbox(c1, c2):
    """
    Convert two opposite corners from (lat, lon) to STAC bbox:
    [west, south, east, north]
    """
    lat1, lon1 = c1
    lat2, lon2 = c2

    west = min(lon1, lon2)
    south = min(lat1, lat2)
    east = max(lon1, lon2)
    north = max(lat1, lat2)

    return [west, south, east, north]


def point_geojson(latlon):
    """
    Convert (lat, lon) to GeoJSON Point with coordinates [lon, lat]
    """
    lat, lon = latlon
    return {
        "type": "Point",
        "coordinates": [lon, lat],
    }


def infer_suffix_from_url(url):
    path = urlparse(url).path.lower()
    if path.endswith(".tif"):
        return ".tif"
    if path.endswith(".tiff"):
        return ".tiff"
    return ".tif"


def download_file(url, output_path):
    output_path = Path(output_path)

    with requests.get(url, stream=True, timeout=300) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        downloaded = 0

        with open(output_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if not chunk:
                    continue
                f.write(chunk)
                downloaded += len(chunk)

                if total:
                    percent = downloaded * 100 / total
                    print(f"\r  {output_path.name}: {percent:5.1f}% ", end="", flush=True)

    if total:
        print()
    print(f"  saved: {output_path.resolve()}")


def build_output_name(prefix, index, total, suffix):
    """
    If only one file is being downloaded, save as:
        prefix.tif
    If multiple files are being downloaded, save as:
        prefix_001.tif, prefix_002.tif, ...
    """
    if total == 1:
        return f"{prefix}{suffix}"
    return f"{prefix}_{index:03d}{suffix}"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Download DEM tiles from Microsoft Planetary Computer.\n\n"
            "Pass either:\n"
            "  - one coordinate: download the full tile(s) containing that point\n"
            "  - two coordinates: download all full tiles intersecting the bbox "
            "defined by those two corners\n\n"
            "Coordinates must be quoted and given as 'lat, lon'."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "coords",
        nargs="+",
        type=parse_latlon,
        help=(
            "One or two coordinates as 'lat, lon'\n"
            "Examples:\n"
            "  \"55.95, -3.21\"\n"
            "  \"55.95, -3.21\" \"55.96, -3.20\""
        ),
    )
    parser.add_argument(
        "-o", "--output-prefix",
        default="my_DEM",
        help="Prefix for output files (default: my_DEM)",
    )
    parser.add_argument(
        "-d", "--output-dir",
        default=".",
        help="Directory to save files into (default: current directory)",
    )
    parser.add_argument(
        "-c", "--collection",
        default="cop-dem-glo-30",
        help="Planetary Computer collection (default: cop-dem-glo-30)",
    )
    parser.add_argument(
        "-a", "--asset",
        default="data",
        help="Asset key to download from each STAC item (default: data)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=200,
        help="Maximum number of matching items to fetch (default: 200)",
    )

    args = parser.parse_args()

    if len(args.coords) not in (1, 2):
        parser.error("You must provide either 1 coordinate or 2 coordinates.")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    catalog = Client.open(CATALOG_URL)

    print(f"Collection: {args.collection}")
    print(f"Asset key:  {args.asset}")

    if len(args.coords) == 1:
        lat, lon = args.coords[0]
        print(f"Mode:       point search")
        print(f"Point:      lat={lat}, lon={lon}")
        print()

        search = catalog.search(
            collections=[args.collection],
            intersects=point_geojson(args.coords[0]),
            limit=args.limit,
        )

    else:
        bbox = to_bbox(args.coords[0], args.coords[1])
        print(f"Mode:       bbox search")
        print(f"BBOX:       {bbox}")
        print()

        search = catalog.search(
            collections=[args.collection],
            bbox=bbox,
            limit=args.limit,
        )

    items = list(search.items())
    if not items:
        raise SystemExit("No matching tile(s) found.")

    items.sort(key=lambda item: item.id)

    print(f"Found {len(items)} tile(s):")
    for i, item in enumerate(items, start=1):
        print(f"  {i:03d}: {item.id}")
    print()

    total = len(items)

    for i, raw_item in enumerate(items, start=1):
        item = planetary_computer.sign(raw_item)

        if args.asset not in item.assets:
            print(f"Skipping {item.id}: missing asset '{args.asset}'")
            print(f"Available assets: {list(item.assets.keys())}")
            continue

        url = item.assets[args.asset].href
        suffix = infer_suffix_from_url(url)
        filename = build_output_name(args.output_prefix, i, total, suffix)
        output_path = out_dir / filename

        print(f"Downloading {item.id} -> {filename}")
        download_file(url, output_path)

    print()
    print("Done.")


if __name__ == "__main__":
    main()
