
#!/usr/bin/env python3
"""
Download OSM waterways in a bounding box and save to GeoPackages.

Outputs:
- waterways_lines.gpkg      : river/stream line features
- water_polygons.gpkg       : lakes/ponds/reservoirs/riverbanks/etc. polygons

Input coordinates are expected in Google Maps style: "lat, lon"

Example:
    python osm_water_download.py \
        --lower-left "55.9421, -3.3347" \
        --upper-right "56.0234, -3.1022" \
        --lines-output edinburgh_waterways.gpkg \
        --polygons-output edinburgh_water_polygons.gpkg
"""

import argparse
import json
import sys
from pathlib import Path

import geopandas as gpd
import osmnx as ox
from shapely.geometry import box

# ---------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------
DEFAULT_LOWER_LEFT = "49.86, -8.65"
DEFAULT_UPPER_RIGHT = "60.86, 1.77"

DEFAULT_LINES_OUTPUT = "waterways_lines.gpkg"
DEFAULT_POLYGONS_OUTPUT = "water_polygons.gpkg"


def parse_latlon_string(s, label):
    """Parse 'lat, lon' into floats."""
    try:
        cleaned = s.strip().replace("(", "").replace(")", "")
        parts = [p.strip() for p in cleaned.split(",")]
        if len(parts) != 2:
            raise ValueError
        lat = float(parts[0])
        lon = float(parts[1])
    except Exception:
        raise ValueError(
            f"Invalid {label}: {s!r}. Expected format like '55.9533, -3.1883'"
        )

    if not (-90 <= lat <= 90):
        raise ValueError(f"{label} latitude out of range: {lat}")
    if not (-180 <= lon <= 180):
        raise ValueError(f"{label} longitude out of range: {lon}")

    return lat, lon


def build_bbox_polygon(lower_left_str, upper_right_str):
    """
    Build a Shapely polygon from:
      lower-left  = southwest corner, given as 'lat, lon'
      upper-right = northeast corner, given as 'lat, lon'
    """
    ll_lat, ll_lon = parse_latlon_string(lower_left_str, "lower-left")
    ur_lat, ur_lon = parse_latlon_string(upper_right_str, "upper-right")

    minx, miny = ll_lon, ll_lat
    maxx, maxy = ur_lon, ur_lat

    if minx > maxx:
        raise ValueError(
            f"Invalid bbox: lower-left lon {minx} > upper-right lon {maxx}"
        )
    if miny > maxy:
        raise ValueError(
            f"Invalid bbox: lower-left lat {miny} > upper-right lat {maxy}"
        )

    return box(minx, miny, maxx, maxy)


def filter_geom_types(gdf, allowed_types):
    """Keep only selected geometry types."""
    if gdf.empty:
        return gdf
    return gdf[gdf.geometry.geom_type.isin(allowed_types)].copy()


def make_gpkg_safe(gdf):
    """
    Convert complex columns like lists/dicts to strings so they can be written
    safely to a GeoPackage.
    """
    if gdf.empty:
        return gdf

    gdf = gdf.reset_index()

    for col in gdf.columns:
        if col == "geometry":
            continue

        gdf[col] = gdf[col].apply(
            lambda v: json.dumps(v) if isinstance(v, (list, dict, tuple, set)) else v
        )

    return gdf


def save_gpkg(gdf, output_path, layer_name):
    """Write GeoDataFrame to a GeoPackage."""
    if gdf.empty:
        print(f"No features to save for {output_path}")
        return

    gdf = make_gpkg_safe(gdf)

    try:
        gdf.to_file(output_path, layer=layer_name, driver="GPKG", engine="pyogrio")
    except Exception:
        # Fallback if pyogrio is not available
        gdf.to_file(output_path, layer=layer_name, driver="GPKG")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download OSM waterways and water polygons in a bounding box."
    )
    parser.add_argument(
        "--lower-left",
        default=DEFAULT_LOWER_LEFT,
        help=f'Lower-left / SW corner as "lat, lon". Default: {DEFAULT_LOWER_LEFT}',
    )
    parser.add_argument(
        "--upper-right",
        default=DEFAULT_UPPER_RIGHT,
        help=f'Upper-right / NE corner as "lat, lon". Default: {DEFAULT_UPPER_RIGHT}',
    )
    parser.add_argument(
        "--lines-output",
        default=DEFAULT_LINES_OUTPUT,
        help=f"Output GeoPackage for rivers/streams. Default: {DEFAULT_LINES_OUTPUT}",
    )
    parser.add_argument(
        "--polygons-output",
        default=DEFAULT_POLYGONS_OUTPUT,
        help=f"Output GeoPackage for lakes/ponds/etc. Default: {DEFAULT_POLYGONS_OUTPUT}",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite output files if they already exist.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    try:
        bbox_polygon = build_bbox_polygon(args.lower_left, args.upper_right)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)

    lines_output = Path(args.lines_output)
    polygons_output = Path(args.polygons_output)

    if not args.overwrite:
        for path in [lines_output, polygons_output]:
            if path.exists():
                print(
                    f"Error: {path} already exists. Use --overwrite to replace it.",
                    file=sys.stderr
                )
                sys.exit(2)
    else:
        for path in [lines_output, polygons_output]:
            if path.exists():
                path.unlink()

    # OSMnx settings
    ox.settings.use_cache = True
    ox.settings.log_console = True
    ox.settings.timeout = 180

    print("Bounding box accepted.")
    print(f"Lower-left:  {args.lower_left}")
    print(f"Upper-right: {args.upper_right}")
    print()

    # -----------------------------------------------------------------
    # 1. Rivers and streams as lines
    # -----------------------------------------------------------------
    line_tags = {
        "waterway": ["river", "stream"]
    }

    print("Downloading river/stream lines from OSM...")
    try:
        waterways = ox.features_from_polygon(bbox_polygon, tags=line_tags)
    except AttributeError:
        # Older OSMnx fallback
        waterways = ox.geometries_from_polygon(bbox_polygon, tags=line_tags)

    waterways = filter_geom_types(waterways, {"LineString", "MultiLineString"})

    print(f"Found {len(waterways):,} line features.")

    # -----------------------------------------------------------------
    # 2. Water polygons: lakes, ponds, reservoirs, riverbanks, etc.
    # -----------------------------------------------------------------
    polygon_tags = {
        "natural": ["water"],
        "water": ["lake", "pond", "reservoir", "lagoon", "basin", "oxbow"],
        "waterway": ["riverbank"],
    }

    print("Downloading water polygons from OSM...")
    try:
        water_polygons = ox.features_from_polygon(bbox_polygon, tags=polygon_tags)
    except AttributeError:
        water_polygons = ox.geometries_from_polygon(bbox_polygon, tags=polygon_tags)

    water_polygons = filter_geom_types(water_polygons, {"Polygon", "MultiPolygon"})

    print(f"Found {len(water_polygons):,} polygon features.")
    print()

    # -----------------------------------------------------------------
    # Save outputs
    # -----------------------------------------------------------------
    if waterways.empty:
        print("No river/stream line features found.")
    else:
        save_gpkg(waterways, lines_output, layer_name="waterways")
        print(f"Saved lines: {lines_output}")

    if water_polygons.empty:
        print("No water polygon features found.")
    else:
        save_gpkg(water_polygons, polygons_output, layer_name="water_polygons")
        print(f"Saved polygons: {polygons_output}")

    print("\nDone.")


if __name__ == "__main__":
    main()
