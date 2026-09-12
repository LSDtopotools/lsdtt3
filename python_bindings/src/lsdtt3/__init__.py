"""Python bindings for LSDTopoTools3.

Authors
-------
B. Gailleton
"""

from functools import wraps
from importlib.metadata import PackageNotFoundError, distribution
import inspect
import os
import re
from pathlib import Path

from ._llm import helper_for_llm as _render_llm_api
from ._compat import MISSING as _MISSING, resolve_deprecated_keyword as _resolve_deprecated_keyword


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

try:
    _distribution = distribution("lsdtt3")
    _installed_package = Path(_distribution.locate_file("lsdtt3")).resolve()
    if _installed_package != Path(__file__).resolve().parent:
        raise PackageNotFoundError
    __version__ = _distribution.version
except PackageNotFoundError:
    # Source-tree and direct CMake installs have no wheel metadata.
    __version__ = "0.0.0"

from ._core import (
    Raster,
    RasterHierarchy,
    RasterMeta,
    VectorBounds,
    OutletFetchResult,
    FlowInfo,
    FlowGraph,
    HydraulicSurfaceFill,
    HydraulicResult,
    GraphFloodOptimisedPerformance,
    GraphFloodOptimisedResult,
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
    DepressionHierarchyResult,
    HillslopeSteadyStateResult,
    HillslopeTransportResult,
    HillslopeCellularResult,
    LandscapeModel,
    LandscapeStepResult,
    SedimentRoutingResult,
    LandslideResult,
    LandslideSeriesResult,
    SyntheticNetworkResult,
    TidalChannelResult,
    TidalSkeletonMetricsResult,
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
    fill_pits_priority_flood,
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
    basin_ridgeline_extrema,
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
    run_graphflood_optimised,
    run_graphflood_optimised_v2,
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
    resample_raster,
    double_raster,
    halve_raster,
    pad_to_nearest_power_of_two,
    interpolate_to_nearest_power_of_two,
    relief_ratio,
    scale,
    sediment_connectivity_index,
    route_sediment,
    trigger_landslides,
    run_landslide_series,
    landslide_crn_concentration,
    landslide_crn_concentrations,
    synthetic_noise,
    synthetic_mountain_range,
    synthetic_hillslope_strip,
    synthetic_network,
    find_tidal_channels,
    analyse_tidal_skeleton,
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
    depression_hierarchy,
    hillslope_steady_state,
    hillslope_transport_explicit,
    hillslope_transport_implicit,
    hillslope_transport_adi,
    hillslope_transport_cellular,
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
    crn_transient_history,
)

from ._graphflood_accelerators import (
    GraphFloodInducedResult,
    GraphFloodMultigridResult,
    run_graphflood_induced,
    run_graphflood_multigrid,
)

__all__ = [
    "Raster",
    "RasterHierarchy",
    "RasterMeta",
    "VectorBounds",
    "OutletFetchResult",
    "FlowInfo",
    "FlowGraph",
    "HydraulicSurfaceFill",
    "HydraulicResult",
    "GraphFloodOptimisedPerformance",
    "GraphFloodOptimisedResult",
    "GraphFloodMultigridResult",
    "GraphFloodInducedResult",
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
    "DepressionHierarchyResult",
    "HillslopeSteadyStateResult",
    "HillslopeTransportResult",
    "HillslopeCellularResult",
    "LandscapeModel",
    "LandscapeStepResult",
    "SedimentRoutingResult",
    "LandslideResult",
    "LandslideSeriesResult",
    "SyntheticNetworkResult",
    "TidalChannelResult",
    "TidalSkeletonMetricsResult",
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
    "fill_pits_priority_flood",
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
    "basin_ridgeline_extrema",
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
    "run_graphflood_optimised",
    "run_graphflood_optimised_v2",
    "run_graphflood_multigrid",
    "run_graphflood_induced",
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
    "resample_raster",
    "double_raster",
    "halve_raster",
    "pad_to_nearest_power_of_two",
    "interpolate_to_nearest_power_of_two",
    "relief_ratio",
    "scale",
    "sediment_connectivity_index",
    "route_sediment",
    "trigger_landslides",
    "run_landslide_series",
    "landslide_crn_concentration",
    "landslide_crn_concentrations",
    "synthetic_noise",
    "synthetic_mountain_range",
    "synthetic_hillslope_strip",
    "synthetic_network",
    "find_tidal_channels",
    "analyse_tidal_skeleton",
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
    "depression_hierarchy",
    "hillslope_steady_state",
    "hillslope_transport_explicit",
    "hillslope_transport_implicit",
    "hillslope_transport_adi",
    "hillslope_transport_cellular",
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
    "crn_transient_history",
    "helper_for_llm",
]


def helper_for_llm(topic="index"):
    """Return an ultra-compact, LLM-oriented API card.

    Parameters
    ----------
    topic : str, default="index"
        Topic name; use ``"symbol:<name>"`` for one public symbol,
        ``"search:<words>"`` for a compact name search, or ``"all"``.

    Returns
    -------
    str
        Machine-oriented API context with exact public identifiers.

    Authors
    -------
    B. Gailleton
    """
    return _render_llm_api(topic, _namespace=globals())


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


def _with_deprecated_keywords(function, *pairs):
    """Accept legacy keywords while preserving the native signature."""
    @wraps(function)
    def compatible(*args, **kwargs):
        for canonical_name, legacy_name in pairs:
            canonical = kwargs.pop(canonical_name, _MISSING)
            legacy = kwargs.pop(legacy_name, _MISSING)
            value = _resolve_deprecated_keyword(
                canonical_name, canonical, legacy_name, legacy, stacklevel=4
            )
            if value is not _MISSING:
                kwargs[canonical_name] = value
        return function(*args, **kwargs)
    native = inspect.unwrap(function)
    definition, _doc, defaults = native.__nb_signature__[0]
    parameters = []
    for match in re.finditer(
        r"(?:\(|, )([A-Za-z_]\w*):([^,)]*)", definition
    ):
        name, declaration = match.groups()
        default_match = re.search(r" = \\(\d+)", declaration)
        default = inspect.Parameter.empty
        if default_match is not None:
            default = defaults[int(default_match.group(1))]
        parameters.append(inspect.Parameter(
            name, inspect.Parameter.POSITIONAL_OR_KEYWORD, default=default
        ))
    compatible.__signature__ = inspect.Signature(parameters)
    return compatible


def _numpy_style_class_docstring(cls):
    """Give an exported native-backed class a consistent public docstring.

    Parameters
    ----------
    cls : type
        Native-backed class exported by :mod:`lsdtt3._core`.

    Returns
    -------
    type
        The same class, with its original documentation retained under
        ``Notes`` and a NumPy-style class-documentation envelope.

    Authors
    -------
    B. Gailleton
    """
    original_doc = cls.__doc__ or "No additional binding notes are available."
    summary = original_doc.strip().split("\n", maxsplit=1)[0].strip()
    if not summary:
        summary = f"Native-backed :mod:`lsdtt3` class ``{cls.__name__}``."

    cls.__doc__ = f"""{summary}

    Parameters
    ----------
    *args
        Constructor arguments described by the generated binding signature.
    **kwargs
        Constructor keyword arguments described by the generated binding
        signature.

    Attributes
    ----------
    See the generated binding signature for available properties and methods.

    Notes
    -----
    {original_doc.strip()}

    Authors
    -------
    B. Gailleton
    """
    return cls


for _graphflood_name in (
    "run_graphflood", "run_graphflood_optimised", "run_graphflood_optimised_v2"
):
    globals()[_graphflood_name] = _with_deprecated_keywords(
        globals()[_graphflood_name],
        ("median_dhw_convergence_tolerance_metres", "tol_median_abs_dhw"),
        ("median_dhw_relative_convergence_tolerance", "tol_median_rel_dhw"),
    )

find_flat_surfaces = _with_deprecated_keywords(
    find_flat_surfaces, ("relief_window_side_pixels", "relief_window_pixels")
)

for _public_name in __all__:
    _public_value = globals()[_public_name]
    if isinstance(_public_value, type):
        globals()[_public_name] = _numpy_style_class_docstring(_public_value)
    elif callable(_public_value):
        globals()[_public_name] = _numpy_style_docstring(_public_value)

del _graphflood_name, _public_name, _public_value
