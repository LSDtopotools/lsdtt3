"""LSDTopoTools3 (lsdtt3): topographic analysis of DEMs from Python.

lsdtt3 is the Python interface to the LSDTopoTools3 C++ library. It fetches
and conditions DEMs, routes flow, extracts channel networks and drainage
basins, and measures landscapes: surface metrics, chi and channel steepness
(ksn), swath profiles, cosmogenic nuclides, flooding (GraphFlood), hillslopes
and landscape evolution. help(lsdtt3.<name>) documents every function and
class; the wheel also installs the lsdtt-* command-line programs.

What it does, and the functions to start with
---------------------------------------------
Rasters, terrain      Raster.read, Raster.from_numpy, Raster.write, hillshade,
                      multidirectional_hillshade, slope, polyfit_metrics
                      (gradient, curvatures), classify_surface_shape,
                      mask_to_nodata, reproject, mosaic, spectral_analysis
Depressions           fill_pits, breach_depressions, breach_least_cost,
                      depression_hierarchy, find_internal_drainage
Fetching data         fetch_dem, fetch_dem_for_outlets, bounds_from_vector,
                      utm_epsg_for_point
Oblique frames        RotationFrame (a grid aligned with a fault, river or
                      range), fetch_dem(frame=...), rotate_to_frame,
                      rotate_back
Flow, channels,       FlowInfo (D8), FlowGraph (MFD, D-infinity),
basins, outlets       FlowInfo.channel_scope, FlowInfo.channel_network,
                      FlowInfo.snap_points (lat/lon or x/y to (row, col)),
                      FlowInfo.basins, junction_basins, extract_channel_heads
Chi and ksn           ChannelNetwork.chi -> ChiCollection.fit -> ChiFit
                      (nodes, segments, summary), ChiFit.invert_ksn,
                      ChannelNetwork.chi_sweep
Swath profiles        swath_profile, swath_distances
Cosmogenic nuclides   crn_production_scaling, crn_steady_state,
                      crn_invert_basin_erosion, CRNColumn
GraphFlood            run_graphflood, run_graphflood_optimised_v2,
                      fill_hydraulic_surface
Hillslopes, LEM       hillslope_transport_implicit, hillslope_steady_state,
                      LandscapeModel, route_sediment, trigger_landslides
Synthetic landscapes  synthetic_noise, synthetic_mountain_range,
                      synthetic_network
Statistics            BasinCollection.statistics, basin_surface_metrics,
                      sample_points, linear_stats, compare_bearings

Conventions
-----------
* Grids are projected, with coordinates and cell sizes in metres. fetch_dem
  returns UTM by default; reproject geographic rasters before analysis.
* Slopes are gradients in m/m (rise over run), never degrees.
* A = R.to_numpy() is a copy indexed A[row, col], and row 0 is the NORTH edge.
  With meta = R.metadata, the centre of pixel (row, col) is
      x = meta.x_min + (col + 0.5) * meta.cell_size
      y = meta.y_max - (row + 0.5) * meta.cell_size
  meta.row_col_to_xy(row, col) and meta.xy_to_row_col(x, y) apply these
  formulas to arrays. R.mask_numpy() is True on valid pixels. Pixels are
  given as (row, col): turn latitude/longitude or x/y points (outlets,
  samples) into (row, col) with FlowInfo.snap_points.
* Functions return new in-memory objects and never change their inputs.
  Nothing is written to disk unless you call R.write(path) (or save a
  GeoDataFrame yourself).
* Tables are pandas DataFrames; vector results are GeoDataFrames in the
  raster's CRS (R.metadata.crs_wkt), and vector inputs must use that CRS.
  Basin IDs in tables are zero-based.

Core workflow
-------------
Fetch a DEM, remove the sea, route flow, extract channels, snap outlets given
as (latitude, longitude) onto the channels, extract basins and fit chi
profiles for ksn::

    import lsdtt3

    # 1. Download a 30 m Copernicus DEM for a WGS84 box (west, south, east, north).
    #    It comes back projected to the UTM zone of the box centre, in metres.
    dem = lsdtt3.fetch_dem(15.63, 37.88, 16.21, 38.31, source="cop30_aws", grid_spacing=30.0)
    dem = lsdtt3.mask_to_nodata(dem, 0.0)  # Copernicus stores the sea as 0 m
    filled = lsdtt3.fill_pits(dem)

    # 2. D8 flow routing and the channel network.
    threshold = 500  # contributing pixels needed to start a channel
    F = lsdtt3.FlowInfo(filled)
    N = F.channel_network(F.channel_scope(threshold_contributing_pixels=threshold))

    # 3. Snap (latitude, longitude) outlets onto the channel network: each point
    #    follows the flow downstream to the first pixel that drains at least
    #    `threshold` pixels, as basin_outlet_snapping_method: downstream does in
    #    the CLIs. The result is an N x 2 array of (row, col).
    outlet_pixels = F.snap_points(
        latitude=[37.94234, 38.06143],
        longitude=[15.79008, 16.06060],
        method="downstream",
        threshold_contributing_pixels=threshold,
    )

    # 4. Basins: one per outlet, basin_id 0, 1, ... in outlet order.
    B = F.basins(outlet_pixels)
    print(B.table[["basin_id", "outlet_row", "outlet_col", "area_metres2", "complete"]])
    polygons = B.polygons  # GeoDataFrame in the DEM's CRS

    # 5. Chi and normalised channel steepness (ksn) by segment fitting.
    fit = N.chi(outlets=outlet_pixels, m_over_n=0.45).fit(filled, critical_divergence_metres=10.0)
    print(fit.segments[["basin_id", "segment_id", "ksn"]])
    print(fit.summary[["basin_id", "median_ksn"]])

    # 6. Nothing is written to disk unless you ask.
    filled.write("filled.tif")
    lsdtt3.hillshade(filled).write("hillshade.tif")
    polygons.to_file("basins.gpkg")

Help for AI assistants
----------------------
For a compact guide with exact signatures, this workflow and the common
mistakes, run print(lsdtt3.helper_for_llm("all")). Single topics:
helper_for_llm("workflow"), and likewise "gotchas", "terrain", "fetch",
"flow", "channel", "basin", "chi", "vector", "hydraulic", "landscape",
"cosmo", "drainage", "synthetic" and "stats"; helper_for_llm("symbol:ChiFit")
gives one symbol's documentation and helper_for_llm("search:chi ksn")
searches names and docstrings.

Authors
-------
B. Gailleton
"""

from functools import wraps
from importlib.metadata import PackageNotFoundError, distribution
import inspect
import math
import os
import re
import struct
import sys
import textwrap
import warnings
import types
import typing
from pathlib import Path

from ._authors import authors_for as _authors_for
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

    # The authoritative mechanism. lsdtt::gdal_init() probes <root>/share/proj
    # and <root>/share/gdal and pins what it finds, in preference to the
    # environment. Left to itself it probes relative to the RUNNING EXECUTABLE,
    # which in a Python process is the interpreter -- so it would pin the conda
    # or system environment's PROJ data rather than the wheel's. Pointing the
    # probe at this package directory is what stops that.
    #
    # This must happen before anything touches GDAL, which is why it is here and
    # not somewhere more convenient. Importing _core does not itself initialise
    # GDAL; the first call that needs it does.
    #
    # GDAL_DATA, PROJ_DATA and PROJ_LIB are deliberately NOT set. lsdtt3 does
    # not need them (the pinning above uses OSRSetPROJSearchPaths and
    # CPLSetConfigOption inside this process), and the variables are read by
    # every OTHER geospatial library in the process: a conda pyproj or GMT,
    # which has no data directory of its own, would then open this wheel's
    # proj.db, which may be older than its PROJ accepts, and fail.
    try:
        from . import _core

        _core.set_bundled_data_root(str(package_dir))
    except Exception:
        # A source tree or a partially built package: lsdtt3 then uses the
        # interpreter-relative probe or the environment, as before.
        pass

    # The Linux wheel bundles a libcurl built on AlmaLinux, which looks for CA
    # certificates only where AlmaLinux keeps them. On Debian/Ubuntu (Colab
    # included) and other distributions that path does not exist, and every
    # https:// read through /vsicurl fails certificate verification. GDAL reads
    # CURL_CA_BUNDLE, so point it at whichever system bundle exists. setdefault
    # and the SSL_CERT_FILE check leave a deliberate choice alone.
    if (
        sys.platform.startswith("linux")
        and "CURL_CA_BUNDLE" not in os.environ
        and "SSL_CERT_FILE" not in os.environ
    ):
        for bundle in (
            "/etc/ssl/certs/ca-certificates.crt",  # Debian, Ubuntu, Arch, Alpine
            "/etc/pki/tls/certs/ca-bundle.crt",  # RHEL, Fedora, Rocky, Alma
            "/etc/ssl/ca-bundle.pem",  # openSUSE
            "/etc/ssl/cert.pem",
        ):
            if os.path.isfile(bundle):
                os.environ["CURL_CA_BUNDLE"] = bundle
                break


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
    FloodplainThresholdResult,
    RiserResult,
    PerimeterSources,
    ValleyWidthNetwork,
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
    RotationFrame,
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
    BurnCorridorResult,
    RiverLinesResult,
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
    absolute,
    perona_malik_smoothing,
    coherence_diffusion,
    stencil_difference_magnitude,
    vertical_units_to_metres,
    coregister_dem,
    normalize_unit_interval,
    plan_curvature,
    polyfit_metrics,
    rasterize_polygons,
    reproject,
    rotate_back,
    rotate_to_frame,
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
    burn_to_corridor,
    fetch_river_lines,
    find_internal_drainage,
    depression_hierarchy,
    hillslope_steady_state,
    hillslope_transport_explicit,
    hillslope_transport_implicit,
    hillslope_transport_adi,
    hillslope_transport_cellular,
    opentopography_datasets,
    dem_sources,
    select_threshold,
    classify_floodplain_terraces_by_thresholds,
    cluster_surface_levels,
    compare_surface_classes,
    surface_class_agreement,
    score_surface_classes,
    relief_above_nearest_baseline,
    find_risers,
    score_risers_against_mapped,
    find_perimeter_sources,
    build_valley_width_network,
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
    # GDAL/PROJ runtime diagnostics. Public because they are what a user pastes
    # into a bug report about coordinates that look subtly wrong -- a mismatched
    # proj.db returns plausible values rather than failing.
    #
    # set_bundled_data_root is deliberately NOT re-exported: it only has an
    # effect before the first GDAL call, which for a user is always in the past
    # by the time they could reach for it. It stays available as
    # lsdtt3._core.set_bundled_data_root for the import-time wiring above.
    gdal_data_description,
    gdal_version,
    proj_data_description,
    proj_version,
)

from ._graphflood_accelerators import (
    GraphFloodInducedResult,
    GraphFloodMultigridResult,
    run_graphflood_induced,
    run_graphflood_multigrid,
)

__all__ = [
    "gdal_data_description",
    "gdal_version",
    "proj_data_description",
    "proj_version",
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
    "FloodplainThresholdResult",
    "RiserResult",
    "PerimeterSources",
    "ValleyWidthNetwork",
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
    "RotationFrame",
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
    "BurnCorridorResult",
    "RiverLinesResult",
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
    "absolute",
    "perona_malik_smoothing",
    "coherence_diffusion",
    "stencil_difference_magnitude",
    "vertical_units_to_metres",
    "coregister_dem",
    "normalize_unit_interval",
    "plan_curvature",
    "polyfit_metrics",
    "rasterize_polygons",
    "reproject",
    "rotate_back",
    "rotate_to_frame",
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
    "burn_to_corridor",
    "fetch_river_lines",
    "find_internal_drainage",
    "depression_hierarchy",
    "hillslope_steady_state",
    "hillslope_transport_explicit",
    "hillslope_transport_implicit",
    "hillslope_transport_adi",
    "hillslope_transport_cellular",
    "opentopography_datasets",
    "dem_sources",
    "select_threshold",
    "classify_floodplain_terraces_by_thresholds",
    "cluster_surface_levels",
    "compare_surface_classes",
    "surface_class_agreement",
    "score_surface_classes",
    "relief_above_nearest_baseline",
    "find_risers",
    "score_risers_against_mapped",
    "find_perimeter_sources",
    "build_valley_width_network",
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
    """Return a compact API guide for AI assistants (and for people in a hurry).

    Parameters
    ----------
    topic : str, default="index"
        ``"index"`` (the default) lists the topics and the notation;
        ``"workflow"`` is a runnable core analysis; ``"gotchas"`` lists common
        mistakes; ``"all"`` is every topic, sized for pasting into an
        assistant. ``"symbol:<name>"`` gives the documentation of one public
        symbol and ``"search:<words>"`` ranks public symbols by how many of
        the words their names and docstrings contain.

    Returns
    -------
    str
        Plain text with exact public identifiers. It is a ``str`` subclass
        that notebooks display as text rather than as an escaped literal, so
        ``lsdtt3.helper_for_llm("all")`` as the last line of a cell reads
        well; ``print()`` works too.

    Authors
    -------
    B. Gailleton
    """
    return _render_llm_api(topic, _namespace=globals())


# ---------------------------------------------------------------------------
# Generated public docstrings and signatures
#
# Every export is wrapped (functions) or re-documented (classes) below so that
# help(), inspect.signature() and the API manifest all describe the real
# parameters. The text is assembled at import time from the nanobind
# signatures (__nb_signature__) and the binding's own docstring, so it is the
# same on every Python version: binding docstrings are normalised with
# inspect.cleandoc (see _notes_block) before anything is added to them.
# ---------------------------------------------------------------------------

_EMPTY = inspect.Parameter.empty

# NumPy docstring section headings recognised in binding docstrings.
_SECTION_NAMES = (
    "Parameters", "Other Parameters", "Returns", "Yields", "Raises", "Warns",
    "Warnings", "See Also", "Notes", "References", "Examples", "Attributes",
    "Methods",
)

_PLACEHOLDER = "No additional binding notes are available."
_UNDESCRIBED = "Not yet described."


class _Annotation:
    """A type annotation already rendered as text; inspect prints its repr()."""

    __slots__ = ("text",)

    def __init__(self, text):
        self.text = text

    def __repr__(self):
        return self.text


class _DisplayFloat(float):
    """A float default shown as the shortest text that gives the same value.

    nanobind reports C++ ``float`` defaults as doubles, so ``1e-4f`` reads as
    ``9.999999747378752e-05``. The value is unchanged; only its repr is
    shortened, so signatures and docstrings show ``0.0001``.
    """

    def __repr__(self):
        return _format_float(float(self))


def _as_float32(value):
    return struct.unpack("f", struct.pack("f", value))[0]


def _format_float(value):
    """Shortest repr of ``value``, treating a single-precision value as exact."""
    try:
        single = _as_float32(value)
    except (OverflowError, struct.error):
        return repr(value)
    if not math.isfinite(value) or single != value:
        return repr(value)
    for digits in range(1, 10):
        text = f"{value:.{digits}g}"
        try:
            if _as_float32(float(text)) == value:
                return repr(float(text))
        except OverflowError:  # rounded beyond the single-precision range
            continue
    return repr(value)


def _display_default(value):
    if type(value) is float and repr(value) != _format_float(value):
        return _DisplayFloat(value)
    return value


def _short_type(text):
    """Render a type for readers: ``lsdtt3._core.Raster`` becomes ``Raster``."""
    text = re.sub(r"\blsdtt3\.(?:_core\.|_graphflood_accelerators\.)?", "", text)
    return text.replace("collections.abc.", "")


def _format_annotation(annotation):
    """Render an annotation identically on every supported Python.

    A union is always written ``A | B``, as the API manifest generator does
    (tools/generate_public_api_manifest.py): before 3.14 a union containing a
    typing generic prints as ``Union[...]``.
    """
    if isinstance(annotation, str):
        return annotation
    if typing.get_origin(annotation) in (typing.Union, types.UnionType):
        return " | ".join(_format_annotation(arg) for arg in typing.get_args(annotation))
    if annotation is type(None):
        return "None"
    return inspect.formatannotation(annotation)


class _Param(typing.NamedTuple):
    """One parameter of a call form; ``annotation`` is text or None."""

    name: str
    kind: inspect._ParameterKind
    annotation: typing.Optional[str]
    default: object


class _Form(typing.NamedTuple):
    """One way of calling an export: a nanobind overload or a Python signature."""

    name: str
    params: tuple
    returned: typing.Optional[str]
    doc: str


def _split_top_level(text):
    """Split ``text`` at commas that are not inside brackets."""
    parts, depth, start = [], 0, 0
    for index, character in enumerate(text):
        if character in "([{":
            depth += 1
        elif character in ")]}":
            depth -= 1
        elif character == "," and depth == 0:
            parts.append(text[start:index].strip())
            start = index + 1
    tail = text[start:].strip()
    if tail:
        parts.append(tail)
    return parts


def _parse_nb_definition(definition, defaults):
    """Parse one nanobind ``__nb_signature__`` definition.

    Parameters
    ----------
    definition : str
        Text such as ``def f(a: int, b: float = \\0, *, c: str) -> None``;
        ``\\N`` refers to ``defaults[N]``.
    defaults : tuple
        Default values of the definition.

    Returns
    -------
    tuple
        ``(name, params, returned)``, or None if the text is not understood.

    Authors
    -------
    S. M. Mudd
    """
    text = definition.strip()
    if text.startswith("def "):
        text = text[4:]
    open_at = text.find("(")
    if open_at < 0:
        return None
    depth, close_at = 0, -1
    for index in range(open_at, len(text)):
        if text[index] in "([{":
            depth += 1
        elif text[index] in ")]}":
            depth -= 1
            if depth == 0:
                close_at = index
                break
    if close_at < 0:
        return None
    rest = text[close_at + 1:].strip()
    returned = rest[2:].strip() if rest.startswith("->") else None

    params = []
    kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
    for entry in _split_top_level(text[open_at + 1:close_at]):
        if entry == "/":
            params = [param._replace(kind=inspect.Parameter.POSITIONAL_ONLY) for param in params]
            continue
        if entry == "*":
            kind = inspect.Parameter.KEYWORD_ONLY
            continue
        match = re.fullmatch(
            r"(\*{0,2})([A-Za-z_]\w*)\s*(?::\s*(.*?))?\s*(?:=\s*(.*))?", entry, re.DOTALL
        )
        if match is None:
            return None
        stars, name, annotation, default_text = match.groups()
        default = _EMPTY
        if default_text is not None:
            reference = re.fullmatch(r"\\(\d+)", default_text.strip())
            if reference is not None and int(reference.group(1)) < len(defaults):
                default = defaults[int(reference.group(1))]
            else:
                default = _Annotation(default_text.strip())
        if stars == "*":
            params.append(_Param(name, inspect.Parameter.VAR_POSITIONAL, annotation, _EMPTY))
            kind = inspect.Parameter.KEYWORD_ONLY
        elif stars == "**":
            params.append(_Param(name, inspect.Parameter.VAR_KEYWORD, annotation, _EMPTY))
        else:
            params.append(_Param(name, kind, annotation, default))
    return text[:open_at].strip(), tuple(params), returned


def _call_forms(value):
    """Return the call forms of a callable, from nanobind or from inspect."""
    native = inspect.unwrap(value)
    nb_signatures = getattr(native, "__nb_signature__", None)
    if nb_signatures is not None:
        forms = []
        for definition, doc, defaults in nb_signatures:
            parsed = _parse_nb_definition(definition, tuple(defaults or ()))
            if parsed is not None:
                forms.append(_Form(*parsed, inspect.cleandoc(doc) if doc else ""))
        return forms
    try:
        signature = inspect.signature(native)
    except (TypeError, ValueError):
        return []
    params = tuple(
        _Param(
            parameter.name,
            parameter.kind,
            None if parameter.annotation is _EMPTY else _format_annotation(parameter.annotation),
            parameter.default,
        )
        for parameter in signature.parameters.values()
    )
    returned = signature.return_annotation
    returned = None if returned is _EMPTY else _format_annotation(returned)
    return [_Form(getattr(native, "__name__", ""), params, returned, "")]


def _signature(params, returned=None, drop_self=True):
    """An inspect.Signature with readable annotations and defaults."""
    parameters = [
        inspect.Parameter(
            param.name,
            param.kind,
            default=_display_default(param.default),
            annotation=_EMPTY if param.annotation is None else _Annotation(_short_type(param.annotation)),
        )
        for param in params
        if not (drop_self and param.name == "self" and param.annotation is None)
    ]
    annotation = _EMPTY if returned is None else _Annotation(_short_type(returned))
    try:
        return inspect.Signature(parameters, return_annotation=annotation)
    except ValueError:
        # An order Python itself would reject (a default before a required
        # parameter); show it as nanobind declares it.
        return inspect.Signature(
            parameters, return_annotation=annotation, __validate_parameters__=False
        )


def _native_signature(function, drop_self=False):
    """inspect.Signature of the first call form of a native export, or None."""
    native = inspect.unwrap(function)
    if getattr(native, "__nb_signature__", None) is None:
        return None
    forms = _call_forms(native)
    if not forms:
        return None
    return _signature(forms[0].params, forms[0].returned, drop_self=drop_self)


def _notes_block(doc):
    """Normalise a binding docstring so it can be merged into a generated one.

    Python 3.13+ strips the indentation of source docstrings at compile time
    and earlier versions keep it, so a docstring pasted in as it stands would
    give text that differs between versions. Normalise it with
    inspect.cleandoc, so every Python version produces the same text, and
    drop its own Authors section (the generated docstring carries the authors).
    """
    return _without_authors(inspect.cleandoc(doc or ""))


def _without_authors(doc):
    """Drop a trailing ``Authors`` section from a binding's own docstring.

    The generated envelope carries the authors (see :mod:`lsdtt3._authors`),
    so a binding docstring that ends with its own ``Authors`` section would
    otherwise list them twice, and could disagree with the envelope.
    """
    return re.sub(r"\n*Authors\n-{3,}\n(?:[^\n]+(?:\n|$))+\s*$", "", doc).rstrip()


def _authors_block(name):
    """The author lines for export ``name``, indented for the docstring template."""
    return "\n    ".join(_authors_for(name))


def _as_docstring(text, name):
    """Finish a generated docstring: add the authors and indent like source.

    Lines after the first carry four spaces, as a docstring written in a
    module does on Python 3.12, so inspect.getdoc() gives the same text on
    every version.
    """
    first, _newline, rest = text.strip().partition("\n")
    body = f"{rest}\n\nAuthors\n-------\n" if rest else "\nAuthors\n-------\n"
    return f"{first}\n{textwrap.indent(body, '    ')}    {_authors_block(name)}\n    "


def _split_sections(doc):
    """Split a cleaned NumPy docstring into its preamble and its sections."""
    lines = doc.splitlines()
    preamble, sections, current = [], [], None
    index = 0
    while index < len(lines):
        line = lines[index]
        if (
            line.rstrip() in _SECTION_NAMES
            and index + 1 < len(lines)
            and re.fullmatch(r"-{3,}\s*", lines[index + 1])
        ):
            current = (line.rstrip(), [])
            sections.append(current)
            index += 2
            continue
        (current[1] if current is not None else preamble).append(line)
        index += 1
    return (
        "\n".join(preamble).strip("\n"),
        [(heading, "\n".join(body).strip("\n")) for heading, body in sections],
    )


def _first_paragraph(text):
    """Split prose into its first paragraph and the rest."""
    first, _blank, rest = text.strip().partition("\n\n")
    return first.strip(), rest.strip("\n")


def _one_line(text):
    """The first paragraph of ``text`` as a single line."""
    return " ".join(_first_paragraph(text)[0].split())


def _strip_signature_lines(doc, name):
    """Drop leading lines that repeat the signature, as nanobind's __doc__ has."""
    lines = doc.splitlines()
    while lines and name and lines[0].startswith(f"{name}("):
        lines.pop(0)
    return "\n".join(lines).strip("\n")


def _merged_params(forms, drop_self=True):
    """Parameters of all call forms, in order of first appearance."""
    merged = {}
    for form in forms:
        for param in form.params:
            if drop_self and param.name == "self" and param.annotation is None:
                continue
            entry = merged.setdefault(param.name, [param, []])
            if param.annotation is not None and param.annotation not in entry[1]:
                entry[1].append(param.annotation)
    return [(param, annotations) for param, annotations in merged.values()]


def _param_label(param):
    if param.kind is inspect.Parameter.VAR_POSITIONAL:
        return f"*{param.name}"
    if param.kind is inspect.Parameter.VAR_KEYWORD:
        return f"**{param.name}"
    return param.name


def _param_entry(param, annotations):
    """One generated ``Parameters`` entry, typed from the signature."""
    head = _param_label(param)
    details = []
    if annotations:
        details.append(" or ".join(_short_type(annotation) for annotation in annotations))
    if param.default is not _EMPTY and not isinstance(param.default, _Annotation):
        details.append(f"default={_display_default(param.default)!r}")
    elif isinstance(param.default, _Annotation):
        details.append(f"default={param.default!r}")
    if param.kind is inspect.Parameter.KEYWORD_ONLY:
        details.append("keyword-only")
    if details:
        head += " : " + ", ".join(details)
    return f"{head}\n    {_UNDESCRIBED}"


def _documented_names(section):
    """Names a binding's own ``Parameters`` section describes."""
    names = set()
    for line in section.splitlines():
        if not line.strip() or line.startswith("    "):
            continue
        head = line.split(" : ", 1)[0]
        names.update(name.strip().lstrip("*") for name in head.split(","))
    return names


def _parameters_section(existing, forms, drop_self=True):
    """The ``Parameters`` section: the binding's own, completed, or generated."""
    params = _merged_params(forms, drop_self=drop_self)
    if existing is not None:
        documented = _documented_names(existing)
        missing = [
            _param_entry(param, annotations)
            for param, annotations in params
            if param.name not in documented
        ]
        return "\n".join([existing] + missing) if existing.strip() else "\n".join(missing)
    if not params:
        return "None"
    return "\n".join(_param_entry(param, annotations) for param, annotations in params)


def _returns_section(forms):
    """A generated ``Returns`` section, typed from the signature(s)."""
    returned = []
    for form in forms:
        text = _short_type(form.returned or "object")
        if text not in returned:
            returned.append(text)
    text = " or ".join(returned)
    if text == "None":
        return "None"
    if text == "Raster":
        return "Raster\n    New in-memory raster; save it with ``.write(path)``."
    return f"{text}\n    {_UNDESCRIBED}"


def _section(heading, body):
    return f"{heading}\n{'-' * len(heading)}\n{body}"


def _section_entries(section):
    """Split a ``Parameters``-style section into ``(head, body lines)`` entries."""
    entries = []
    for line in section.splitlines():
        if line.strip() and not line.startswith(" "):
            entries.append((line.rstrip(), []))
        elif entries:
            entries[-1][1].append(line.rstrip())
    return [(head, "\n".join(body).strip("\n")) for head, body in entries]


def _entry_names(head):
    return tuple(name.strip().lstrip("*") for name in head.split(" : ", 1)[0].split(","))


def _entry_type(head):
    return head.split(" : ", 1)[1].strip() if " : " in head else ""


def _merged_form_parameters(form_sections):
    """One ``Parameters`` section from the sections of several call forms.

    Entries are kept in order of first appearance. A parameter documented by
    more than one form keeps one entry: identical text is shown once, and a
    parameter whose type or description differs between forms gets the types
    joined with "or" and each description labelled with its call form.
    """
    merged = {}
    for number, section in form_sections:
        for head, body in _section_entries(section):
            names = _entry_names(head)
            key = next((existing for existing in merged if set(names) <= set(existing)), None)
            if key is None:
                merged[names] = [names, [_entry_type(head)], [(number, body)]]
                continue
            entry = merged[key]
            if names != key:
                continue  # already documented as part of a grouped entry
            kind = _entry_type(head)
            if kind and kind not in entry[1]:
                entry[1].append(kind)
            if body and all(body != known for _form, known in entry[2]):
                entry[2].append((number, body))
    lines = []
    for names, types, bodies in merged.values():
        types = [kind for kind in types if kind]
        head = ", ".join(names) + (" : " + " or ".join(types) if types else "")
        lines.append(head)
        if len(bodies) == 1:
            lines.append(bodies[0][1])
        else:
            for number, body in bodies:
                text = textwrap.dedent(body).strip()
                lines.append(textwrap.indent(f"Call form {number}: {text}", "    "))
    return "\n".join(line for line in lines if line)


def _merged_form_returns(form_sections):
    """One ``Returns`` section from the sections of several call forms, or None."""
    distinct = []
    for number, section in form_sections:
        if section and all(section != known for _form, known in distinct):
            distinct.append((number, section))
    if not distinct:
        return None
    if len(distinct) == 1:
        return distinct[0][1]
    lines = []
    for number, section in distinct:
        for head, body in _section_entries(section):
            lines.append(head)
            text = textwrap.dedent(body).strip() or "Result."
            lines.append(textwrap.indent(f"Call form {number}: {text}", "    "))
    return "\n".join(lines)


def _function_text(name, function):
    """Assemble the docstring text (without authors) of an exported function."""
    forms = _call_forms(function)
    native = inspect.unwrap(function)
    if not forms:
        return _notes_block(native.__doc__) or _PLACEHOLDER

    signature_lines = [f"{name}{_signature(form.params, form.returned)}" for form in forms]
    if getattr(native, "__nb_signature__", None) is None:
        source_doc = _notes_block(native.__doc__)
    elif len(forms) == 1:
        source_doc = _without_authors(forms[0].doc)
    else:
        source_doc = ""
    source_doc = _strip_signature_lines(source_doc, name)
    preamble, sections = _split_sections(source_doc)

    if len(forms) > 1:
        # Each form keeps its own prose and other sections; the Parameters
        # and Returns of all forms are merged into the function's own
        # sections, so every parameter is described once.
        summary = _one_line(_without_authors(forms[0].doc)) or _PLACEHOLDER
        items = [f"This function has {len(forms)} call forms:"]
        form_parameters, form_returns = [], []
        for number, (form, line) in enumerate(zip(forms, signature_lines), 1):
            items.append(f"{number}. ``{line}``")
            form_preamble, form_sections = _split_sections(
                _strip_signature_lines(_without_authors(form.doc), name)
            )
            kept = [form_preamble] if form_preamble else []
            for heading, body in form_sections:
                if heading == "Parameters":
                    form_parameters.append((number, body))
                elif heading == "Returns":
                    form_returns.append((number, body))
                else:
                    kept.append(_section(heading, body))
            if kept:
                items.append(textwrap.indent("\n\n".join(kept), "   "))
        rest = "\n\n".join(items)
        if form_parameters:
            sections.append(("Parameters", _merged_form_parameters(form_parameters)))
        merged_returns = _merged_form_returns(form_returns)
        if merged_returns:
            sections.append(("Returns", merged_returns))
    else:
        summary, rest = _first_paragraph(preamble)
        summary = summary or _PLACEHOLDER

    # A single call form is shown by help(), inspect.signature() and notebook
    # tooltips from __signature__, so the docstring does not repeat it; the
    # call forms of an overloaded function are listed in ``rest``.
    blocks = [summary]
    if rest:
        blocks.append(rest)
    found = dict(sections)
    blocks.append(_section("Parameters", _parameters_section(found.get("Parameters"), forms)))
    blocks.append(_section("Returns", found.get("Returns") or _returns_section(forms)))
    blocks.extend(
        _section(heading, body) for heading, body in sections
        if heading not in ("Parameters", "Returns")
    )
    return "\n\n".join(blocks)


def _numpy_style_docstring(function, name=None):
    """Give an exported function a public docstring and a real signature.

    The docstring has a prose summary, the binding's description, a
    ``Parameters`` section naming every real parameter with its type and
    default (the binding's own section, completed from the signature if it
    misses any), a ``Returns`` section, the binding's other sections, and the
    authors. ``__signature__`` is set from the nanobind signature, so
    ``help()``, ``inspect.signature`` and notebook tooltips show the typed
    signature; the docstring lists it only for an overloaded function, which
    has several. The low-level function stays available through
    ``__wrapped__``.

    Parameters
    ----------
    function : callable
        Low-level function exported by :mod:`lsdtt3._core`, or a Python
        function wrapping one.
    name : str, optional
        Public name, used to look up its authors. Defaults to
        ``function.__name__``.

    Returns
    -------
    callable
        A forwarding public function with a NumPy-style docstring.

    Authors
    -------
    B. Gailleton
    S. M. Mudd
    """
    name = name or function.__name__

    @wraps(function)
    def documented(*args, **kwargs):
        return function(*args, **kwargs)

    signature = _native_signature(function)
    if signature is not None:
        documented.__signature__ = signature
    documented.__doc__ = _as_docstring(_function_text(name, function), name)
    return documented


def _with_deprecated_keywords(function, *pairs, is_method=False):
    """Accept legacy keywords while preserving the native signature."""
    @wraps(function)
    def compatible(*args, **kwargs):
        for canonical_name, legacy_name in pairs:
            canonical = kwargs.pop(canonical_name, _MISSING)
            legacy = kwargs.pop(legacy_name, _MISSING)
            value = _resolve_deprecated_keyword(
                canonical_name, canonical, legacy_name, legacy,
                stacklevel=3 if is_method else 4,
            )
            if value is not _MISSING:
                kwargs[canonical_name] = value
        return function(*args, **kwargs)
    # The native definition of a method names ``self``, so the signature
    # keeps it and a bound method drops it as usual.
    signature = _native_signature(function)
    if signature is not None:
        compatible.__signature__ = signature
    return compatible


def _producers(namespace):
    """Map each exported class name to the exports that return it.

    Values are ``(factories, others)``: ``factories`` are the class's own
    static methods that return it (``Raster.read()``); ``others`` are labels
    ``function()``, ``Class.method()`` or ``Class.property``, exported
    functions first, then members of other classes, then the class's own
    instance members.
    """
    classes = {name for name in namespace["__all__"] if isinstance(namespace[name], type)}
    found = {}

    def note(returned, label, rank):
        for class_name in re.findall(r"\b([A-Z]\w*)\b", _short_type(returned or "")):
            if class_name in classes:
                labels = found.setdefault(class_name, [])
                if (rank, label) not in labels:
                    labels.append((rank, label))

    for name in namespace["__all__"]:
        value = namespace[name]
        if isinstance(value, type):
            for member, raw, forms in _class_methods(value):
                for form in forms:
                    own = name in re.findall(r"\b([A-Z]\w*)\b", _short_type(form.returned or ""))
                    rank = (0 if _is_static(raw) else 3) if own else 2
                    note(form.returned, f"{name}.{member}()", rank)
            for member, returned, _doc in _class_properties(value):
                note(returned, f"{name}.{member}", 2)
        elif callable(value):
            for form in _call_forms(value):
                note(form.returned, f"{name}()", 1)
    producers = {}
    for class_name, labels in found.items():
        ordered = sorted(labels, key=lambda item: item[0])
        producers[class_name] = (
            [label for rank, label in ordered if rank == 0 and label.startswith(f"{class_name}.")],
            [label for rank, label in ordered if not (rank == 0 and label.startswith(f"{class_name}."))],
        )
    return producers


def _class_properties(cls):
    """``(name, return type text, doc)`` for every public property of ``cls``."""
    properties = []
    for member in sorted(name for name in dir(cls) if not name.startswith("_")):
        try:
            raw = inspect.getattr_static(cls, member)
        except AttributeError:
            continue
        if not isinstance(raw, property):
            continue
        returned, doc = None, raw.__doc__ or ""
        nb_signatures = getattr(raw.fget, "__nb_signature__", None)
        if nb_signatures:
            parsed = _parse_nb_definition(nb_signatures[0][0], ())
            returned = parsed[2] if parsed else None
            doc = nb_signatures[0][1] or ""
        elif raw.fget is not None:
            try:
                annotation = inspect.signature(raw.fget).return_annotation
            except (TypeError, ValueError):
                annotation = _EMPTY
            returned = None if annotation is _EMPTY else _format_annotation(annotation)
        if doc.startswith("(self)"):
            doc = ""
        properties.append((member, returned, _notes_block(doc)))
    return properties


def _class_methods(cls):
    """``(name, raw attribute, call forms)`` for every public method of ``cls``."""
    methods = []
    for member in sorted(name for name in dir(cls) if not name.startswith("_")):
        try:
            raw = inspect.getattr_static(cls, member)
        except AttributeError:
            continue
        if isinstance(raw, property) or not callable(getattr(cls, member, None)):
            continue
        forms = _call_forms(getattr(cls, member))
        if forms:
            methods.append((member, raw, forms))
    return methods


def _is_static(raw):
    return isinstance(raw, staticmethod) or type(raw).__name__ == "nb_func"


def _constructor_forms(cls):
    """Call forms of the constructor, or [] if Python cannot construct ``cls``."""
    init = inspect.getattr_static(cls, "__init__", None)
    if getattr(init, "__nb_signature__", None) is not None:
        return _call_forms(init)
    if inspect.isfunction(init):
        forms = _call_forms(init)
        return [form._replace(returned=None) for form in forms]
    return []


def _class_text(name, cls, producers):
    """Assemble the docstring text (without authors) of an exported class."""
    own_doc = _notes_block(cls.__doc__)
    preamble, sections = _split_sections(own_doc)
    found = dict(sections)
    constructor = _constructor_forms(cls)
    factories, others = producers if producers else ([], [])
    made_by = ", ".join(f"``{label}``" for label in others[:4])
    if len(others) > 4:
        made_by += f" and {len(others) - 4} more"

    if constructor:
        lines = [
            f"{name}{_signature(form.params, None)}" for form in constructor
        ]
        origin = "Construct it as\n\n" + "\n".join(f"    {line}" for line in lines)
        if made_by:
            origin += f"\n\nIt is also returned by {made_by}."
    elif factories:
        origin = (
            "Not constructible directly; create one with "
            + " or ".join(f"``{label}``" for label in factories)
        )
        origin += f". It is also returned by {made_by}." if made_by else "."
    else:
        origin = f"Returned by {made_by or 'other lsdtt3 calls'}; not constructible directly."

    summary, rest = _first_paragraph(preamble)
    if summary:
        blocks = [summary, rest, origin] if rest else [summary, origin]
    else:
        blocks = [origin, _PLACEHOLDER]
    if constructor and "Parameters" not in found:
        blocks.append(_section("Parameters", _parameters_section(None, constructor)))
    blocks.extend(
        _section(heading, body) for heading, body in sections
        if heading not in ("Attributes", "Methods")
    )

    if "Attributes" in found:
        blocks.append(_section("Attributes", found["Attributes"]))
    else:
        properties = _class_properties(cls)
        if properties:
            entries = []
            for member, returned, doc in properties:
                head = member if returned is None else f"{member} : {_short_type(returned)}"
                line = _one_line(doc)
                entries.append(f"{head}\n    {line}" if line else head)
            blocks.append(_section("Attributes", "\n".join(entries)))

    if "Methods" in found:
        blocks.append(_section("Methods", found["Methods"]))
    else:
        entries = []
        for member, raw, forms in _class_methods(cls):
            static = _is_static(raw)
            prefix = f"{name}.{member}" if static else member
            for form in forms:
                line = f"{prefix}{_signature(form.params, form.returned, drop_self=not static)}"
                summary_line = _one_line(_strip_signature_lines(_without_authors(form.doc), member))
                entries.append(f"{line}\n    {summary_line}" if summary_line else line)
        if entries:
            blocks.append(_section("Methods", "\n".join(entries)))
    return "\n\n".join(blocks)


def _numpy_style_class_docstring(cls, name=None, producers=None):
    """Give an exported class a public docstring that lists what it offers.

    The docstring says how to get an instance (its constructor signature and
    parameters, or which exports return it when Python cannot construct it),
    keeps the binding's own documentation, and lists every public property
    and method with its signature and the first line of its documentation.
    A constructible native class also gets ``__signature__``.

    Parameters
    ----------
    cls : type
        Class exported by :mod:`lsdtt3`.
    name : str, optional
        Public name, used to look up its authors. Defaults to
        ``cls.__name__``.
    producers : tuple of (list of str, list of str), optional
        Labels of the class's own factory methods and of the other exports
        that return it, as made by ``_producers``.

    Returns
    -------
    type
        The same class, re-documented.

    Authors
    -------
    B. Gailleton
    S. M. Mudd
    """
    name = name or cls.__name__
    text = _class_text(name, cls, producers)
    constructor = _constructor_forms(cls)
    if constructor and getattr(inspect.getattr_static(cls, "__init__"), "__nb_signature__", None):
        try:
            cls.__signature__ = _signature(constructor[0].params, None)
        except (AttributeError, TypeError):
            pass
    cls.__doc__ = _as_docstring(text, name)
    return cls


for _graphflood_name in (
    "run_graphflood", "run_graphflood_optimised", "run_graphflood_optimised_v2"
):
    globals()[_graphflood_name] = _with_deprecated_keywords(
        globals()[_graphflood_name],
        ("median_dhw_convergence_tolerance_metres", "tol_median_abs_dhw"),
        ("median_dhw_relative_convergence_tolerance", "tol_median_rel_dhw"),
    )

def _with_relief_window_side(function):
    """Accept the deprecated full-side relief window, converted to a radius.

    ``relief_window_side_pixels`` (and the older ``relief_window_pixels``)
    gave the window as a full side in cells; ``relief_window_radius_metres``
    gives it as a radius in metres. A side of s cells is a radius of
    (s - 1) / 2 cells, converted with the DEM cell size exactly as
    lsdtt-flat-surface-finder converts the same old parameter.
    """
    @wraps(function)
    def compatible(*args, **kwargs):
        side = kwargs.pop("relief_window_side_pixels", _MISSING)
        legacy = kwargs.pop("relief_window_pixels", _MISSING)
        side = _resolve_deprecated_keyword(
            "relief_window_side_pixels", side, "relief_window_pixels", legacy,
            stacklevel=3)
        if side is not _MISSING:
            if "relief_window_radius_metres" in kwargs:
                raise TypeError(
                    "'relief_window_radius_metres' and deprecated "
                    "'relief_window_side_pixels' cannot both be supplied")
            warnings.warn(
                "'relief_window_side_pixels' is deprecated; use "
                "'relief_window_radius_metres'", DeprecationWarning, stacklevel=2)
            side = int(side)
            if side < 3 or side % 2 == 0:
                raise ValueError(
                    f"relief_window_side_pixels must be odd and >= 3; got {side}")
            dem = args[0] if args else kwargs["dem"]
            kwargs["relief_window_radius_metres"] = (
                (side - 1) // 2) * dem.metadata.cell_size
        return function(*args, **kwargs)
    signature = _native_signature(function)
    if signature is not None:
        compatible.__signature__ = signature
    return compatible


find_flat_surfaces = _with_relief_window_side(_with_deprecated_keywords(
    find_flat_surfaces, ("max_relief_metres", "relief_threshold_metres")
))

for _classify_name in ("classify_flat_surfaces", "classify_flat_surfaces_with_scope"):
    globals()[_classify_name] = _with_deprecated_keywords(
        globals()[_classify_name],
        ("max_floodplain_hand_metres", "floodplain_max_height_metres"),
        ("min_terrace_hand_metres", "terrace_min_height_metres"),
        ("valley_corridor_half_width_metres", "terrace_max_distance_to_channel_metres"),
    )

wiener_filter = _with_deprecated_keywords(
    wiener_filter,
    ("signal_wavelength_long_metres", "signal_wavelength_long_m"),
    ("signal_wavelength_short_metres", "signal_wavelength_short_m"),
)

route_sediment = _with_deprecated_keywords(
    route_sediment,
    ("active_layer_depth_metres", "active_layer_depth_m"),
    ("d50_metres", "d50_m"),
    ("d84_metres", "d84_m"),
    ("d90_metres", "d90_m"),
)

find_internal_drainage = _with_deprecated_keywords(
    find_internal_drainage,
    ("min_area_m2", "min_area_metres2"),
    ("min_volume_m3", "min_volume_metres3"),
)

for _flow_threshold_method in (
    "channel_mask", "channel_scope", "distance_to_channel", "sources"
):
    setattr(
        FlowInfo,
        _flow_threshold_method,
        _with_deprecated_keywords(
            getattr(FlowInfo, _flow_threshold_method),
            ("threshold_contributing_pixels", "threshold_pixels"),
            is_method=True,
        ),
    )

for _crn_depth_method, _canonical_name, _legacy_name in (
    ("set_erosion_history", "removal_metres", "removal_m"),
    ("transient", "sampling_depths_metres", "sampling_depths_m"),
):
    setattr(
        CRNColumn,
        _crn_depth_method,
        _with_deprecated_keywords(
            getattr(CRNColumn, _crn_depth_method),
            (_canonical_name, _legacy_name),
            is_method=True,
        ),
    )

extract_channel_heads = _with_deprecated_keywords(
    extract_channel_heads,
    ("threshold_contributing_pixels", "threshold_pixels"),
    ("surface_fitting_radius_metres", "surface_fitting_radius_m"),
)

crn_transient_history = _with_deprecated_keywords(
    crn_transient_history,
    ("sampling_depth_metres", "sampling_depth_m"),
)

fetch_dem = _with_deprecated_keywords(fetch_dem, ("size_pixels", "size_px"))

synthetic_network = _with_deprecated_keywords(
    synthetic_network, ("relief_metres", "relief_m")
)

run_landslide_series = _with_deprecated_keywords(
    run_landslide_series,
    ("event_budget_relative_tolerance", "budget_tol"),
)

_class_producers = _producers(globals())
for _public_name in __all__:
    _public_value = globals()[_public_name]
    if isinstance(_public_value, type):
        globals()[_public_name] = _numpy_style_class_docstring(
            _public_value, _public_name, _class_producers.get(_public_name)
        )
    elif callable(_public_value):
        globals()[_public_name] = _numpy_style_docstring(_public_value, _public_name)

del _canonical_name, _class_producers, _crn_depth_method, _flow_threshold_method
del _graphflood_name, _legacy_name, _public_name, _public_value
