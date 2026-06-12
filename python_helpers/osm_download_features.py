#!/usr/bin/env python3
"""
Download OSM water features, roads, or both, clipped to a raster's extent.

Unlike osm_water_download.py (which takes a lat/lon bounding box), this script
reads the extent and CRS directly from a reference raster (e.g. your DEM), so the
downloaded features line up with the raster you are analysing. The raster bounding
box is reprojected to WGS84 for the OSM query, then the returned features are
reprojected into the raster's CRS and clipped exactly to its extent.

Outputs (FlatGeobuf, the lsdtt3-preferred vector format; one geometry type each):
  <prefix>_osm_roads.fgb            road / track centrelines (LineString)
  <prefix>_osm_waterways.fgb        river / stream lines     (LineString)
  <prefix>_osm_water_polygons.fgb   lakes / ponds / etc.     (Polygon)

Examples:
    # roads only, clipped to a DEM, output next to the DEM
    python osm_download_features.py --raster ../data/strathduie/dem_crop.tif \
        --features roads

    # both water and roads, custom output prefix
    python osm_download_features.py --raster mydem.tif --features both \
        --output-prefix mysite

Requires: rasterio, geopandas, osmnx, shapely, pyproj.
"""

import argparse
import json
import sys
from pathlib import Path

import geopandas as gpd
import osmnx as ox
import rasterio
from shapely.geometry import box

# ---------------------------------------------------------------------
# Tag sets
# ---------------------------------------------------------------------
# Road/track network. "track" is included deliberately — forest/agricultural
# tracks are often the features of interest. Footways/paths/cycleways are left
# out; add them here if you want them.
ROAD_HIGHWAY_VALUES = [
    "motorway", "trunk", "primary", "secondary", "tertiary",
    "unclassified", "residential", "service", "track", "road",
    "living_street",
]

WATER_LINE_TAGS = {
    "waterway": ["river", "stream", "canal", "drain", "ditch"],
}

WATER_POLYGON_TAGS = {
    "natural": ["water"],
    "water": ["lake", "pond", "reservoir", "lagoon", "basin", "oxbow"],
    "waterway": ["riverbank"],
}

WGS84 = "EPSG:4326"


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def get_raster_extent(raster_path):
    """
    Return (bbox_raster_crs, raster_crs, bbox_wgs84) for the raster.

    bbox_raster_crs : shapely box in the raster's own CRS (for final clipping)
    raster_crs      : rasterio/pyproj CRS of the raster
    bbox_wgs84      : shapely box in EPSG:4326 (for the OSM query)
    """
    with rasterio.open(raster_path) as src:
        if src.crs is None:
            raise ValueError(
                f"Raster {raster_path} has no CRS; cannot georeference the query."
            )
        b = src.bounds  # left, bottom, right, top in raster CRS
        raster_crs = src.crs
        print(f"Raster      : {raster_path}")
        print(f"Size        : {src.width} x {src.height} px, "
              f"cell {src.res[0]:g} x {src.res[1]:g}")
        print(f"CRS         : {raster_crs.to_string()}")
        print(f"Extent      : ({b.left:.1f}, {b.bottom:.1f}) - "
              f"({b.right:.1f}, {b.top:.1f})")

    bbox_raster_crs = box(b.left, b.bottom, b.right, b.top)

    # Reproject the bbox corners to WGS84 for the OSM query.
    bbox_gs = gpd.GeoSeries([bbox_raster_crs], crs=raster_crs)
    bbox_wgs84 = bbox_gs.to_crs(WGS84).iloc[0]
    wb = bbox_wgs84.bounds
    print(f"Query bbox  : lon [{wb[0]:.5f}, {wb[2]:.5f}], "
          f"lat [{wb[1]:.5f}, {wb[3]:.5f}] (WGS84)")
    print()
    return bbox_raster_crs, raster_crs, bbox_wgs84


def osm_features_from_polygon(polygon_wgs84, tags):
    """Query OSM features within a WGS84 polygon, across osmnx versions."""
    try:
        return ox.features_from_polygon(polygon_wgs84, tags=tags)
    except AttributeError:
        # Older osmnx
        return ox.geometries_from_polygon(polygon_wgs84, tags=tags)


def filter_geom_types(gdf, allowed_types):
    """Keep only the selected geometry types."""
    if gdf is None or gdf.empty:
        return gdf
    return gdf[gdf.geometry.geom_type.isin(allowed_types)].copy()


def make_vector_safe(gdf):
    """Stringify list/dict/tuple columns so the file driver can write them."""
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


def clip_and_save(gdf, allowed_types, bbox_raster_crs, raster_crs,
                  out_path, layer_name, overwrite):
    """Reproject to the raster CRS, clip to its extent, and write FlatGeobuf."""
    gdf = filter_geom_types(gdf, allowed_types)
    if gdf is None or gdf.empty:
        print(f"  No features for {out_path.name}; skipping.")
        return

    gdf = gdf.to_crs(raster_crs)
    clip_frame = gpd.GeoDataFrame(geometry=[bbox_raster_crs], crs=raster_crs)
    clipped = gpd.clip(gdf, clip_frame)
    clipped = clipped[~clipped.geometry.is_empty & clipped.geometry.notna()]

    if clipped.empty:
        print(f"  Nothing left inside the raster for {out_path.name}; skipping.")
        return

    if out_path.exists():
        if not overwrite:
            print(f"  {out_path} exists (use --overwrite to replace); skipping.")
            return
        out_path.unlink()

    clipped = make_vector_safe(clipped)
    clipped.to_file(out_path, driver="FlatGeobuf", layer=layer_name)
    print(f"  Saved {len(clipped):,} features -> {out_path}")


# ---------------------------------------------------------------------
# Download wrappers
# ---------------------------------------------------------------------
def download_roads(polygon_wgs84, bbox_raster_crs, raster_crs, prefix, overwrite):
    print("Downloading OSM roads/tracks...")
    tags = {"highway": ROAD_HIGHWAY_VALUES}
    feats = osm_features_from_polygon(polygon_wgs84, tags)
    clip_and_save(
        feats, {"LineString", "MultiLineString"},
        bbox_raster_crs, raster_crs,
        prefix.with_name(prefix.name + "_osm_roads.fgb"),
        "osm_roads", overwrite,
    )


def download_water(polygon_wgs84, bbox_raster_crs, raster_crs, prefix, overwrite):
    print("Downloading OSM water lines...")
    lines = osm_features_from_polygon(polygon_wgs84, WATER_LINE_TAGS)
    clip_and_save(
        lines, {"LineString", "MultiLineString"},
        bbox_raster_crs, raster_crs,
        prefix.with_name(prefix.name + "_osm_waterways.fgb"),
        "osm_waterways", overwrite,
    )

    print("Downloading OSM water polygons...")
    polys = osm_features_from_polygon(polygon_wgs84, WATER_POLYGON_TAGS)
    clip_and_save(
        polys, {"Polygon", "MultiPolygon"},
        bbox_raster_crs, raster_crs,
        prefix.with_name(prefix.name + "_osm_water_polygons.fgb"),
        "osm_water_polygons", overwrite,
    )


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Download OSM water/roads clipped to a reference raster.")
    p.add_argument("--raster", required=True,
                   help="Reference raster whose extent and CRS define the clip.")
    p.add_argument("--features", choices=["water", "roads", "both"],
                   default="both", help="Which features to download (default: both).")
    p.add_argument("--output-prefix", default=None,
                   help="Output path prefix. Default: the raster path without its "
                        "extension (outputs land beside the raster).")
    p.add_argument("--overwrite", action="store_true",
                   help="Overwrite existing output files.")
    return p.parse_args()


def main():
    args = parse_args()

    raster_path = Path(args.raster)
    if not raster_path.exists():
        print(f"Error: raster {raster_path} not found.", file=sys.stderr)
        sys.exit(2)

    prefix = Path(args.output_prefix) if args.output_prefix \
        else raster_path.with_suffix("")

    try:
        bbox_raster_crs, raster_crs, bbox_wgs84 = get_raster_extent(raster_path)
    except Exception as e:
        print(f"Error reading raster: {e}", file=sys.stderr)
        sys.exit(2)

    ox.settings.use_cache = True
    ox.settings.log_console = False
    ox.settings.timeout = 180

    if args.features in ("roads", "both"):
        download_roads(bbox_wgs84, bbox_raster_crs, raster_crs, prefix,
                       args.overwrite)
    if args.features in ("water", "both"):
        download_water(bbox_wgs84, bbox_raster_crs, raster_crs, prefix,
                       args.overwrite)

    print("\nDone.")


if __name__ == "__main__":
    main()
