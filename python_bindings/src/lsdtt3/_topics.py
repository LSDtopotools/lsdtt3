"""Which topic each public export belongs to.

One rule, shared by the generated API reference pages
(tools/generate_api_reference.py, one page per domain) and the LLM guide
(``helper_for_llm``), so that both group the exports the same way. Every export
gets exactly one domain; ``domain`` falls back to ``"raster"``.

This module imports nothing from the package: the reference generator loads it
by file path in CI, where the compiled extension is not installed.

Authors
-------
S. M. Mudd
"""

DOMAINS = {
    "raster": "Raster And Terrain",
    "flow": "Flow, Channels, And Basins",
    "chi": "Chi And Ksn",
    "fetch": "Fetching And Geospatial Data",
    "cosmogenic": "Cosmogenic Nuclides",
    "graphflood": "GraphFlood And Hydraulics",
    "landscape": "Landscape And Hillslopes",
    "drainage": "Drainage And Depressions",
    "sediment": "Sediment, Landslides, And Landforms",
    "synthetic": "Synthetic Workflows",
    "fluvial": "Fluvial, Valley, Marsh, And Vector Hydrology",
    "statistics": "Statistics And Utilities",
}


# GDAL/PROJ runtime diagnostics. Matched by exact name rather than by substring:
# a "proj" substring rule would also capture reproject / reproject_to_grid, which
# belong with the raster operations.
GEOSPATIAL_RUNTIME = {
    "gdal_data_description",
    "gdal_version",
    "proj_data_description",
    "proj_version",
}


def domain(name: str) -> str:
    """Return the domain key (a key of ``DOMAINS``) of a public export name."""
    lower = name.lower()
    if name in GEOSPATIAL_RUNTIME:
        return "fetch"
    if "synthetic" in lower:
        return "synthetic"
    if "graphflood" in lower or "hydraulic" in lower:
        return "graphflood"
    if lower.startswith("crn") or "cosmo" in lower:
        return "cosmogenic"
    if "chi" in lower:
        return "chi"
    if any(word in lower for word in ("fetch", "bounds", "cop30", "utm_epsg", "http_available", "rotat", "dem_sources")):
        return "fetch"
    if any(word in lower for word in ("landscape", "hillslope")):
        return "landscape"
    if any(word in lower for word in ("internaldrainage", "internal_drainage", "depression", "soil_topographic", "sti", "drainresult", "find_drains")):
        return "drainage"
    if any(word in lower for word in ("sediment", "landslide", "road", "mountainfront", "mountain_front", "ridgeline", "thinlines", "thin_lines")):
        return "sediment"
    if any(word in lower for word in ("valley", "tidal", "marsh", "water", "streamburn", "burn_stream", "repair_water",
                                      "floodplain", "surface_class", "surface_levels", "riser", "relief_above",
                                      "select_threshold", "perimeter")):
        return "fluvial"
    if any(word in lower for word in ("flow", "channel", "basin", "divide", "connectivity")):
        return "flow"
    if any(word in lower for word in ("compare", "linear_stats", "relief_ratio", "helper_for_llm")):
        return "statistics"
    return "raster"


# -----------------------------------------------------------------------------
# LLM guide topics (helper_for_llm). Finer than DOMAINS: each export is listed
# in exactly one topic's catalogue. A topic's hand-written guidance may mention
# exports of other topics.

#: The catalogue topics, in the order the guide lists them.
GUIDE_TOPICS = ("terrain", "fetch", "flow", "channel", "basin", "chi", "vector", "water",
                "hydraulic", "landscape", "sediment", "cosmo", "drainage", "synthetic", "stats")

#: The topic of an export, by its domain, unless a rule below places it.
GUIDE_TOPIC_OF_DOMAIN = {
    "raster": "terrain",
    "flow": "flow",
    "chi": "chi",
    "fetch": "fetch",
    "cosmogenic": "cosmo",
    "graphflood": "hydraulic",
    "landscape": "landscape",
    "drainage": "drainage",
    "sediment": "sediment",
    "synthetic": "synthetic",
    "fluvial": "vector",
    "statistics": "stats",
}

#: Exports placed in a topic other than their domain's, where the hand-written
#: guidance teaches them (frozen from the cards on 2026-11-19).
GUIDE_TOPIC_PLACEMENTS = {
    # basin
    "basin_ridgeline_extrema": "basin",
    # channel
    "conditioned_pixel_mask": "channel",
    "slope_area_bins": "channel",
    # drainage
    "ConnectivityResult": "drainage",
    "fill_interior_nodata": "drainage",
    "fill_pits_boundary_aware": "drainage",
    "fill_pits_with_outlets": "drainage",
    "fill_pits_with_sinks": "drainage",
    "FlowGraph": "drainage",
    "make_outlet_mask": "drainage",
    "sediment_connectivity_index": "drainage",
    "SoilTopographicIndexResult": "drainage",
    "topographic_wetness_index": "drainage",
    # fetch
    "opentopography_datasets": "fetch",
    # flow
    "make_boundary_conditions": "flow",
    # stats
    "absolute": "stats",
    "add": "stats",
    "divide": "stats",
    "invert_binary_mask": "stats",
    "mask_by_binary": "stats",
    "multiply": "stats",
    "scale": "stats",
    "subtract": "stats",
    "threshold_above": "stats",
    "threshold_below": "stats",
    # synthetic
    "add_uniform_noise": "synthetic",
    # terrain
    "breach_depressions": "terrain",
    # vector
    "clip_to_polygons": "vector",
    "mask_from_points": "vector",
    "rasterize_polygons": "vector",
    "sample_points": "vector",
    "swath_profile": "vector",
    "SwathProfile": "vector",
    # water
    "burn_to_corridor": "water",
    "BurnCorridorResult": "water",
    "fetch_river_lines": "water",
    "RiverLinesResult": "water",
}


def guide_topic(name: str) -> str:
    """Return the LLM guide topic (one of ``GUIDE_TOPICS``) of a public export name.

    ``GUIDE_TOPIC_PLACEMENTS`` first; then, within the flow domain, names with
    ``basin`` or ``divide`` go to ``basin`` and names with ``channel`` to
    ``channel``; within the fluvial domain, names with ``water`` or ``burn`` go
    to ``water``; otherwise ``GUIDE_TOPIC_OF_DOMAIN`` of the export's domain.
    """
    if name in GUIDE_TOPIC_PLACEMENTS:
        return GUIDE_TOPIC_PLACEMENTS[name]
    key = domain(name)
    if key == "flow":
        lower = name.lower()
        if "basin" in lower or "divide" in lower:
            return "basin"
        if "channel" in lower:
            return "channel"
    if key == "fluvial":
        lower = name.lower()
        if "water" in lower or "burn" in lower:
            return "water"
    return GUIDE_TOPIC_OF_DOMAIN[key]
