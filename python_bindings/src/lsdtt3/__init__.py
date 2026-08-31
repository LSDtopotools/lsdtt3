"""Python bindings for LSDTopoTools3.

Authors
-------
B. Gailleton
"""

from functools import wraps
import os
from pathlib import Path


def _configure_bundled_geospatial_data() -> None:
    """Configure bundled GDAL and PROJ data for a standalone wheel.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Authors
    -------
    B. Gailleton
    """
    package_dir = Path(__file__).resolve().parent
    gdal_data = package_dir / "gdal_data"
    proj_data = package_dir / "proj_data"

    if gdal_data.is_dir():
        os.environ.setdefault("GDAL_DATA", str(gdal_data))
    if proj_data.is_dir():
        os.environ.setdefault("PROJ_DATA", str(proj_data))
        # Older GDAL/PROJ combinations use the legacy spelling.
        os.environ.setdefault("PROJ_LIB", str(proj_data))


_configure_bundled_geospatial_data()

from ._core import (
    Raster,
    RasterMeta,
    VectorBounds,
    OutletFetchResult,
    FlowInfo,
    FlowGraph,
    HydraulicSurfaceFill,
    HydraulicResult,
    ChannelHeads,
    ChannelScope,
    ChannelNetwork,
    BasinCollection,
    JunctionBasinCollection,
    SideTributaryBasins,
    ChiCollection,
    ChiFit,
    ChiInversion,
    ConnectivityResult,
    FlatSurfaceClassification,
    FlatSurfaceResult,
    MarshPlatformResult,
    MountainFrontResult,
    SpectralResult,
    STILnkdTable,
    STISprTable,
    STICategoricalRaster,
    SoilTopographicIndexResult,
    SwathProfile,
    ChannelCandidatesResult,
    DrainResult,
    DivideMigrationResult,
    FrangiResult,
    RoadFlowResult,
    RoadResult,
    RoadSmoothnessResult,
    Ridgeline,
    ThinLinesResult,
    ValleyCentrelineResult,
    ValleyCrossSections,
    ValleyResult,
    WaterCentrelineResult,
    WaterLineRepairResult,
    StreamBurnResult,
    InternalDrainageResult,
    HillslopeSteadyStateResult,
    HillslopeTransportResult,
    LandscapeModel,
    LandscapeStepResult,
    SedimentRoutingResult,
    CRNColumn,
    CRNColumnSample,
    CRNBurialResult,
    CRNStepPairResult,
    add,
    add_uniform_noise,
    aspect,
    breach_depressions,
    breach_least_cost,
    carve_depth,
    compare_bearings,
    classify_surface_shape,
    classify_flat_surfaces,
    classify_flat_surfaces_with_scope,
    clip_to_maximum_rectangle,
    clip_to_polygons,
    detect_directional_periodicity,
    divide,
    divide_migration,
    extract_roads_smooth,
    extract_thin_lines,
    extract_valley_centrelines,
    extract_valley_cross_sections,
    valley_mask_from_cross_sections,
    extract_channel_heads,
    fill_depth,
    fill_hydraulic_surface,
    fill_pits,
    fill_pits_boundary_aware,
    fill_pits_with_outlets,
    fill_pits_with_sinks,
    fundamental_form_principal_curvatures,
    fill_interior_nodata,
    fetch_dem,
    fetch_dem_for_outlets,
    junction_basins,
    side_tributary_basins,
    line_crossing_basins,
    basin_channel_lengths,
    basin_surface_metrics,
    find_channel_candidates,
    find_drains,
    find_flat_surfaces,
    find_marsh_platform,
    find_road_flow_discrimination,
    find_roads,
    find_valleys,
    bounds_from_vector,
    fit_valley_baseline,
    frangi_filter,
    hillshade,
    run_graphflood,
    highpass_filter,
    invert_binary_mask,
    least_cost_carve_depth,
    linear_stats,
    lowpass_filter,
    make_boundary_conditions,
    make_outlet_mask,
    mask_by_binary,
    mask_from_points,
    mask_to_nodata,
    masked_surface_recovery,
    identify_mountain_front,
    mountain_front_outlets,
    multidirectional_hillshade,
    mosaic,
    multiply,
    normalize_unit_interval,
    plan_curvature,
    polyfit_metrics,
    rasterize_polygons,
    reproject,
    reproject_to_grid,
    relief_ratio,
    scale,
    sediment_connectivity_index,
    route_sediment,
    spectral_analysis,
    spectral_available,
    soil_topographic_index,
    swath_distances,
    swath_profile,
    set_nodata_rectangle,
    sample_points,
    slope,
    subtract,
    threshold_above,
    threshold_below,
    topographic_wetness_index,
    total_curvature,
    bandpass_filter,
    wiener_filter,
    load_sti_lnkd_table,
    load_sti_spr_table,
    rasterize_sti_categories,
    http_available,
    utm_epsg_for_point,
    cop30_tiles_for_bbox,
    water_centrelines,
    repair_water_lines,
    burn_streams,
    find_internal_drainage,
    hillslope_steady_state,
    hillslope_transport_explicit,
    hillslope_transport_implicit,
    hillslope_transport_adi,
    opentopography_datasets,
    crn_apparent_erosion_rate,
    crn_burial_pair,
    crn_step_change_pair,
    crn_production_scaling,
    crn_pressure,
    crn_uniform_raster,
    crn_steady_state,
    crn_step_change,
    crn_accumulated_concentration,
    crn_outlet_concentration,
    crn_isotope_ratio,
    crn_invert_basin_erosion,
    crn_pair_diagnostics,
)

__all__ = [
    "Raster",
    "RasterMeta",
    "VectorBounds",
    "OutletFetchResult",
    "FlowInfo",
    "FlowGraph",
    "HydraulicSurfaceFill",
    "HydraulicResult",
    "ChannelHeads",
    "ChannelScope",
    "ChannelNetwork",
    "BasinCollection",
    "JunctionBasinCollection",
    "SideTributaryBasins",
    "ChiCollection",
    "ChiFit",
    "ChiInversion",
    "ConnectivityResult",
    "FlatSurfaceClassification",
    "FlatSurfaceResult",
    "MarshPlatformResult",
    "MountainFrontResult",
    "SpectralResult",
    "STILnkdTable",
    "STISprTable",
    "STICategoricalRaster",
    "SoilTopographicIndexResult",
    "SwathProfile",
    "ChannelCandidatesResult",
    "DrainResult",
    "DivideMigrationResult",
    "FrangiResult",
    "RoadFlowResult",
    "RoadResult",
    "RoadSmoothnessResult",
    "Ridgeline",
    "ThinLinesResult",
    "ValleyCentrelineResult",
    "ValleyCrossSections",
    "ValleyResult",
    "WaterCentrelineResult",
    "WaterLineRepairResult",
    "StreamBurnResult",
    "InternalDrainageResult",
    "HillslopeSteadyStateResult",
    "HillslopeTransportResult",
    "CRNColumn",
    "CRNColumnSample",
    "CRNBurialResult",
    "CRNStepPairResult",
    "add",
    "add_uniform_noise",
    "aspect",
    "breach_depressions",
    "breach_least_cost",
    "carve_depth",
    "compare_bearings",
    "classify_surface_shape",
    "classify_flat_surfaces",
    "classify_flat_surfaces_with_scope",
    "clip_to_maximum_rectangle",
    "clip_to_polygons",
    "detect_directional_periodicity",
    "divide",
    "divide_migration",
    "extract_roads_smooth",
    "extract_thin_lines",
    "extract_valley_centrelines",
    "extract_valley_cross_sections",
    "valley_mask_from_cross_sections",
    "extract_channel_heads",
    "fill_depth",
    "fill_hydraulic_surface",
    "fill_pits",
    "fill_pits_boundary_aware",
    "fill_pits_with_outlets",
    "fill_pits_with_sinks",
    "fundamental_form_principal_curvatures",
    "fill_interior_nodata",
    "fetch_dem",
    "fetch_dem_for_outlets",
    "junction_basins",
    "side_tributary_basins",
    "line_crossing_basins",
    "basin_channel_lengths",
    "basin_surface_metrics",
    "find_channel_candidates",
    "find_drains",
    "find_flat_surfaces",
    "find_marsh_platform",
    "find_road_flow_discrimination",
    "find_roads",
    "find_valleys",
    "bounds_from_vector",
    "fit_valley_baseline",
    "frangi_filter",
    "hillshade",
    "run_graphflood",
    "highpass_filter",
    "invert_binary_mask",
    "least_cost_carve_depth",
    "linear_stats",
    "lowpass_filter",
    "make_boundary_conditions",
    "make_outlet_mask",
    "mask_by_binary",
    "mask_from_points",
    "mask_to_nodata",
    "masked_surface_recovery",
    "identify_mountain_front",
    "mountain_front_outlets",
    "multidirectional_hillshade",
    "mosaic",
    "multiply",
    "normalize_unit_interval",
    "plan_curvature",
    "polyfit_metrics",
    "rasterize_polygons",
    "reproject",
    "reproject_to_grid",
    "relief_ratio",
    "scale",
    "sediment_connectivity_index",
    "spectral_analysis",
    "spectral_available",
    "soil_topographic_index",
    "swath_distances",
    "swath_profile",
    "set_nodata_rectangle",
    "sample_points",
    "slope",
    "subtract",
    "threshold_above",
    "threshold_below",
    "topographic_wetness_index",
    "total_curvature",
    "bandpass_filter",
    "wiener_filter",
    "load_sti_lnkd_table",
    "load_sti_spr_table",
    "rasterize_sti_categories",
    "http_available",
    "utm_epsg_for_point",
    "cop30_tiles_for_bbox",
    "water_centrelines",
    "repair_water_lines",
    "burn_streams",
    "find_internal_drainage",
    "hillslope_steady_state",
    "hillslope_transport_explicit",
    "hillslope_transport_implicit",
    "hillslope_transport_adi",
    "opentopography_datasets",
    "crn_apparent_erosion_rate",
    "crn_burial_pair",
    "crn_step_change_pair",
    "crn_production_scaling",
    "crn_pressure",
    "crn_uniform_raster",
    "crn_steady_state",
    "crn_step_change",
    "crn_accumulated_concentration",
    "crn_outlet_concentration",
    "crn_isotope_ratio",
    "crn_invert_basin_erosion",
    "crn_pair_diagnostics",
]


def _numpy_style_docstring(function):
    """Give an exported binding a consistent public Python docstring.

    The low-level nanobind signature remains available through
    ``__wrapped__``.  The original binding documentation is retained in the
    Notes section so this normalisation does not discard operation-specific
    guidance while the C++ bindings are progressively documented in place.

    Parameters
    ----------
    function : callable
        Low-level function exported by :mod:`lsdtt3._core`.

    Returns
    -------
    callable
        A forwarding public function with a NumPy-style docstring.

    Authors
    -------
    B. Gailleton
    """
    original_doc = function.__doc__ or "No additional binding notes are available."
    summary = original_doc.strip().split("\n", maxsplit=1)[0].strip()
    if not summary:
        summary = f"Run :func:`lsdtt3.{function.__name__}`."

    @wraps(function)
    def documented(*args, **kwargs):
        return function(*args, **kwargs)

    documented.__doc__ = f"""{summary}

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    {original_doc.strip()}

    Authors
    -------
    B. Gailleton
    """
    return documented


for _public_name in __all__:
    _public_value = globals()[_public_name]
    if callable(_public_value) and not isinstance(_public_value, type):
        globals()[_public_name] = _numpy_style_docstring(_public_value)

del _public_name, _public_value
