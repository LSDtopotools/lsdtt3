"""Documentation signatures for the public :mod:`lsdtt3` Python API.

This module is deliberately not imported by :mod:`lsdtt3`.  It is a
source-only façade consumed by MkDocs so that native nanobind objects can have
clear, stable Python documentation without changing their implementation or
runtime behaviour.

Authors
-------
B. Gailleton
"""

from __future__ import annotations

from typing import Any, Sequence


class RasterMeta:
    """Metadata for a north-up, square-cell raster grid.

    Attributes
    ----------
    n_rows, n_cols : int
        Grid dimensions.
    x_min, y_min, cell_size : float
        Lower-left grid origin and square cell size in CRS units.
    nodata_value : float
        Stored value representing invalid cells.
    crs_wkt : str
        Coordinate reference system as WKT, when known.
    epsg_code : int
        EPSG identifier, or zero when it is not known.

    Authors
    -------
    B. Gailleton
    """

    @classmethod
    def from_shape(
        cls, n_rows: int, n_cols: int, x_min: float = 0.0,
        y_min: float = 0.0, cell_size: float = 1.0,
        nodata_value: float = -9999.0, crs_wkt: str = "", epsg_code: int = 0,
    ) -> RasterMeta:
        """Create metadata from dimensions and a lower-left origin.

        Parameters
        ----------
        n_rows, n_cols : int
            Positive raster dimensions.
        x_min, y_min : float, default=0
            Lower-left grid origin.
        cell_size : float, default=1
            Square cell size in CRS units.
        nodata_value : float, default=-9999
            Value stored for invalid cells.
        crs_wkt : str, default=""
            Grid coordinate reference system in WKT.
        epsg_code : int, default=0
            Grid EPSG identifier when known.

        Returns
        -------
        RasterMeta
            Metadata compatible with :meth:`Raster.from_numpy`.

        Authors
        -------
        B. Gailleton
        """
        ...


class Raster:
    """An in-memory float raster with an explicit validity mask.

    Raster operations return new rasters and do not modify their inputs.
    Invalid cells are represented by :attr:`RasterMeta.nodata_value` in data
    arrays and by ``False`` in validity masks.

    Attributes
    ----------
    metadata : RasterMeta
        Immutable grid metadata.
    shape : tuple[int, int]
        Array shape as ``(rows, columns)``.
    valid_pixel_count : int
        Number of cells with valid data.

    Authors
    -------
    B. Gailleton
    """

    @classmethod
    def from_numpy(cls, data: Any, mask: Any, metadata: RasterMeta) -> Raster:
        """Create a raster by copying two-dimensional arrays.

        Parameters
        ----------
        data : array-like of float
            Values with shape ``(metadata.n_rows, metadata.n_cols)``.
        mask : array-like of bool
            Validity mask with the same shape; ``True`` marks valid cells.
        metadata : RasterMeta
            Metadata matching both array dimensions.

        Returns
        -------
        Raster
            A new raster containing copied data and mask values.

        Authors
        -------
        B. Gailleton
        """
        ...

    @classmethod
    def read(cls, path: str) -> Raster:
        """Read one band of a GDAL-compatible raster.

        Parameters
        ----------
        path : path-like
            Raster file to load.

        Returns
        -------
        Raster
            Values, grid metadata, and nodata mask.

        Authors
        -------
        B. Gailleton
        """
        ...

    def to_numpy(self) -> Any:
        """Return copied float32 data; invalid cells contain the nodata value.

        Returns
        -------
        numpy.ndarray
            Two-dimensional values array.

        Authors
        -------
        B. Gailleton
        """
        ...

    def mask_numpy(self) -> Any:
        """Return a copied Boolean validity mask.

        Returns
        -------
        numpy.ndarray
            ``True`` at cells with valid values.

        Authors
        -------
        B. Gailleton
        """
        ...

    def write(self, path: str, format: str = "cog") -> None:
        """Write the raster through GDAL.

        Parameters
        ----------
        path : path-like
            Destination file.
        format : str, default="cog"
            LSDTopoTools GDAL output format, for example ``"cog"`` or
            ``"gtiff"``.

        Authors
        -------
        B. Gailleton
        """
        ...


def resample_raster(source: Raster, target_metadata: RasterMeta,
                    method: str = "linear") -> Raster:
    """Resample a raster onto an explicit grid without GDAL.

    Parameters
    ----------
    source : Raster
        Source values and validity mask.
    target_metadata : RasterMeta
        Destination dimensions, extent, resolution, and CRS.
    method : str, default="linear"
        ``nearest``, ``linear``, ``cubic``, ``mean``, ``sum``, ``min``,
        ``max``, ``mode``, or ``conservative_mean``.

    Returns
    -------
    Raster
        Resampled values on the requested grid.

    Notes
    -----
    ``mean`` weights intersecting cells equally. ``conservative_mean`` weights
    them by overlap area. ``sum`` distributes each source-cell value according
    to the fraction of its area covered by the destination cell.

    Authors
    -------
    B. Gailleton
    """
    ...


def double_raster(source: Raster, method: str = "linear") -> Raster:
    """Double both dimensions using the optimized same-extent path.

    Parameters
    ----------
    source : Raster
        Source raster.
    method : str, default="linear"
        Resampling method accepted by :func:`resample_raster`.

    Returns
    -------
    Raster
        Raster with twice as many rows and columns and half-sized cells.

    Authors
    -------
    B. Gailleton
    """
    ...


def halve_raster(source: Raster, method: str = "conservative_mean") -> Raster:
    """Halve both even dimensions using the optimized same-extent path.

    Parameters
    ----------
    source : Raster
        Source raster with even row and column counts.
    method : str, default="conservative_mean"
        Resampling method accepted by :func:`resample_raster`.

    Returns
    -------
    Raster
        Raster with half as many rows and columns and double-sized cells.

    Authors
    -------
    B. Gailleton
    """
    ...


def pad_to_nearest_power_of_two(source: Raster) -> Raster:
    """Pad north and east with nodata to power-of-two dimensions.

    Parameters
    ----------
    source : Raster
        Source raster.

    Returns
    -------
    Raster
        Padded raster retaining the source cell size and south-west origin.

    Authors
    -------
    B. Gailleton
    """
    ...


def interpolate_to_nearest_power_of_two(
    source: Raster, direction: str = "nearest", method: str = "linear",
) -> Raster:
    """Resize to compatible power-of-two dimensions over the same extent.

    Parameters
    ----------
    source : Raster
        Source raster.
    direction : {"lower", "higher", "nearest"}, default="nearest"
        Power-of-two selection for both dimensions.
    method : str, default="linear"
        Resampling method accepted by :func:`resample_raster`.

    Returns
    -------
    Raster
        Resampled power-of-two raster.

    Notes
    -----
    The row and column scale factors must match because raster cells are square.

    Authors
    -------
    B. Gailleton
    """
    ...


class RasterHierarchy:
    """Manage same-extent resolutions for one float raster.

    Level zero is the input, finest grid. The hierarchy retains only the
    buffers selected by ``storage``; restriction and prolongation transfer
    from the active level.

    Authors
    -------
    B. Gailleton
    """

    def __init__(self, finest: Raster, storage: str = "retain_allocated",
                 keep_finest: bool = True, levels: int = 0) -> None:
        """Create a hierarchy.

        Parameters
        ----------
        finest : Raster
            Finest raster. Its dimensions must support each requested exact
            factor-two level. Add arbitrary integer factors later with
            :meth:`add_level_for_factor`.
        storage : {"active_only", "retain_allocated", "preallocated"}
            Buffer retention policy.
        keep_finest : bool, default=True
            Keep level zero when using ``active_only``.
        levels : int, default=0
            Number of levels including level zero; zero uses all exact levels.
        """
        ...

    def restrict_to(self, level: int,
                    method: str = "conservative_mean") -> Raster:
        """Transfer the active grid to a coarser level and activate it."""
        ...

    def prolong_to(self, level: int, method: str = "linear") -> Raster:
        """Transfer the active grid to a finer level and activate it."""
        ...

    def level_factor(self, level: int) -> int:
        """Return the integer coarsening factor relative to level zero."""
        ...

    def level_shape(self, level: int) -> tuple[int, int]:
        """Return ``(n_rows, n_cols)`` for a hierarchy level."""
        ...


def fill_pits(raster: Raster, min_slope: float = 1e-4, n_threads: int = 0,
              priority_flood_backend: str = "binary") -> Raster:
    """Fill depressions, optionally imposing a minimum drainage slope.

    Parameters
    ----------
    raster : Raster
        Input elevation raster.
    min_slope : float, default=1e-4
        Minimum imposed drainage gradient.
    n_threads : int, default=0
        OpenMP thread count; zero uses the runtime default.
    priority_flood_backend : {"binary", "radix"}, default="binary"
        Queue used when ``min_slope > 0`` dispatches to Priority-Flood.
        ``min_slope <= 0`` uses reconstruction instead.

    Returns
    -------
    Raster
        Conditioned elevation raster.

    Authors
    -------
    B. Gailleton
    """
    ...


def fill_pits_priority_flood(raster: Raster, min_slope: float = 1e-4,
                             n_threads: int = 0,
                             use_radix_heap: bool = False) -> Raster:
    """Explicitly run the Improved Priority-Flood depression fill.

    Parameters
    ----------
    raster : Raster
        Elevation raster.
    min_slope : float, default=1e-4
        Minimum slope imposed across filled flats. Zero keeps true flats.
    n_threads : int, default=0
        Reserved for API compatibility; Priority-Flood is serial.
    use_radix_heap : bool, default=False
        Replace the vanilla binary heap with the radix-hybrid queue. The
        filling algorithm and flat FIFO are otherwise unchanged.

    Returns
    -------
    Raster
        Filled raster.
    """
    ...


def breach_depressions(raster: Raster, n_threads: int = 0) -> Raster:
    """Carve drainage paths through depressions.

    Parameters
    ----------
    raster : Raster
        Input elevation raster.
    n_threads : int, default=0
        OpenMP thread count; zero uses the runtime default.

    Returns
    -------
    Raster
        Breached elevation raster.

    Authors
    -------
    B. Gailleton
    """
    ...


def multidirectional_hillshade(raster: Raster, n_directions: int = 8,
                               altitude: float = 45.0, z_factor: float = 1.0,
                               start_azimuth: float = 0.0) -> Raster:
    """Compute an averaged multi-illumination hillshade in the range 0–255.

    Parameters
    ----------
    raster : Raster
        Elevation raster.
    n_directions : int, default=8
        Number of evenly spaced illumination azimuths.
    altitude : float, default=45
        Sun altitude in degrees.
    z_factor : float, default=1
        Vertical exaggeration applied before shading.
    start_azimuth : float, default=0
        First illumination azimuth in degrees.

    Returns
    -------
    Raster
        Hillshade raster.

    Authors
    -------
    B. Gailleton
    """
    ...


def slope(raster: Raster) -> Raster:
    """Compute 3-by-3 gradient-magnitude slope in rise/run units.

    Parameters
    ----------
    raster : Raster
        Elevation raster.

    Returns
    -------
    Raster
        Slope raster.

    Authors
    -------
    B. Gailleton
    """
    ...


def polyfit_metrics(raster: Raster, metrics: Sequence[str], lambda_metres: float = 0.0,
                    z_factor: float = 1.0, n_threads: int = 0,
                    flat_slope_tolerance: float = 1e-12) -> dict[str, Raster]:
    """Compute requested quadratic-fit surface metrics in one shared fit.

    Parameters
    ----------
    raster : Raster
        Input elevation raster.
    metrics : sequence of str
        Requested metric names, such as ``"slope"``, ``"aspect"``,
        ``"curvature"``, or principal curvatures.
    lambda_metres : float, default=0
        Fit scale in metres; zero uses three raster cells.
    z_factor : float, default=1
        Elevation scale factor applied before fitting.
    n_threads : int, default=0
        OpenMP thread count; zero uses the runtime default.
    flat_slope_tolerance : float, default=1e-12
        Directional metrics are nodata below this fitted slope.

    Returns
    -------
    dict[str, Raster]
        Newly allocated metric rasters keyed by requested names.

    Authors
    -------
    B. Gailleton
    """
    ...


def classify_surface_shape(*args: Any, **kwargs: Any) -> dict[str, Raster]:
    """Classify surface shape from an elevation raster or principal curvatures.

    Parameters
    ----------
    raster : Raster
        Elevation raster for the first overload.
    maximum_curvature, minimum_curvature : Raster
        Existing principal-curvature rasters for the second overload.
    lambda_metres, z_factor, n_threads : optional
        Fit settings for the elevation-raster overload.
    kt : float, default=0
        Curvature threshold used for classification.

    Returns
    -------
    dict[str, Raster]
        Surface-shape categories and, when required, computed curvatures.

    Authors
    -------
    B. Gailleton
    """
    ...


def spectral_analysis(raster: Raster, log_bin_width: float = 0.05,
                      n_background_iterations: int = 20, seed: int = -1,
                      window: str = "hann", clip_to_valid_rectangle: bool = False,
                      n_threads: int = 0) -> SpectralResult:
    """Run FFT-based spectral analysis.

    Parameters
    ----------
    raster : Raster
        Elevation raster.
    log_bin_width : float, default=0.05
        Logarithmic spectral-bin width.
    n_background_iterations : int, default=20
        Number of background-spectrum iterations.
    seed : int, default=-1
        Random seed; a negative value uses the native default.
    window : str, default="hann"
        Spectral window name.
    clip_to_valid_rectangle : bool, default=False
        Restrict calculation to the largest valid rectangle.
    n_threads : int, default=0
        OpenMP thread count; zero uses the runtime default.

    Returns
    -------
    SpectralResult
        Spectral products and fitted background information.

    Authors
    -------
    B. Gailleton
    """
    ...


class SpectralResult:
    """Products returned by :func:`spectral_analysis`.

    Attributes
    ----------
    power_spectrum, radial_spectrum, background_spectrum : array-like
        Spectral outputs and background model.

    Authors
    -------
    B. Gailleton
    """


def frangi_filter(raster: Raster, sigma_min_metres: float = 1.0,
                  sigma_max_metres: float = 5.0, n_scales: int = 4,
                  beta: float = 0.5, c: float = 15.0, target: str = "ridges",
                  n_threads: int = 0) -> FrangiResult:
    """Run a multiscale Frangi ridge or valley filter.

    Parameters
    ----------
    raster : Raster
        Raster to filter.
    sigma_min_metres, sigma_max_metres : float
        Smallest and largest Gaussian scales in metres.
    n_scales : int, default=4
        Number of scales.
    beta, c : float
        Frangi vesselness-shape and contrast parameters.
    target : {"ridges", "valleys"}, default="ridges"
        Feature polarity to enhance.
    n_threads : int, default=0
        OpenMP thread count; zero uses the runtime default.

    Returns
    -------
    FrangiResult
        Filter response and scale products.

    Authors
    -------
    B. Gailleton
    """
    ...


class FrangiResult:
    """Multiscale Frangi-filter output rasters.

    Authors
    -------
    B. Gailleton
    """


def swath_profile(reference: Raster, baseline: Any, values: Raster | None = None,
                  half_width_metres: float = 0.0, bin_width_metres: float = 1000.0,
                  n_threads: int = 0) -> SwathProfile:
    """Summarise raster values in a corridor along one GeoPandas LineString.

    Parameters
    ----------
    reference : Raster
        Defines the swath grid and valid-pixel domain.
    baseline : geopandas.GeoDataFrame
        Exactly one LineString in the same projected metre CRS.
    values : Raster, optional
        Co-registered raster to summarise; defaults to ``reference``.
    half_width_metres : float, default=0
        Corridor half-width; zero includes all valid cells.
    bin_width_metres : float, default=1000
        Along-axis bin width.
    n_threads : int, default=0
        OpenMP thread count; zero uses the runtime default.

    Returns
    -------
    SwathProfile
        Along- and across-axis rasters plus binned statistics.

    Authors
    -------
    B. Gailleton
    """
    ...


class SwathProfile:
    """Distance and summary products for a swath profile.

    Attributes
    ----------
    table : pandas.DataFrame
        One row per along-axis bin, including mean and percentile values.
    along_axis, across_axis : Raster
        Per-cell distances relative to the baseline.

    Authors
    -------
    B. Gailleton
    """


class FlowInfo:
    """Immutable D8 drainage topology derived from a conditioned DEM.

    Construct with ``FlowInfo(dem, boundary_conditions="...")``.  It owns
    routing topology rather than elevations; methods explicitly receive a
    raster whenever elevation or another surface is required.

    Key methods
    -----------
    drainage_area()
        Return D8 drainage area.
    channel_scope(...), scope_from_heads(...)
        Define a channel domain from a threshold or detected heads.
    channel_network(scope)
        Build a network from that selected domain.
    basins(...)
        Return raster-labelled drainage basins.

    Authors
    -------
    B. Gailleton
    """


class FlowGraph:
    """Fractional D-infinity or MFD flow-routing graph.

    Construct with ``FlowGraph(dem, method="dinfinity" or "mfd")`` after
    conditioning depressions. Use it for distributed-flow drainage area;
    use :class:`FlowInfo` for D8 channel trees.

    Authors
    -------
    B. Gailleton
    """


def topographic_wetness_index(drainage_area: Raster, slope: Raster) -> Raster:
    """Compute ``log((A / cell_area) / slope)``.

    Parameters
    ----------
    drainage_area : Raster
        D8, MFD, or D-infinity drainage area.
    slope : Raster
        Co-registered slope in rise/run units.

    Returns
    -------
    Raster
        Topographic wetness index.

    Authors
    -------
    B. Gailleton
    """
    ...


def sediment_connectivity_index(flow: FlowInfo, roughness: Raster, slope: Raster,
                                drainage_area: Raster,
                                area_threshold_m2: float) -> ConnectivityResult:
    """Compute a D8 sediment-connectivity index.

    Parameters
    ----------
    flow : FlowInfo
        D8 routing topology.
    roughness, slope : Raster
        Co-registered resistance and slope rasters.
    drainage_area : Raster
        Co-registered drainage area; it may be D8, MFD, or D-infinity.
    area_threshold_m2 : float
        Minimum drainage area defining the channel network.

    Returns
    -------
    ConnectivityResult
        Connectivity index with upslope and downslope components.

    Authors
    -------
    B. Gailleton
    """
    ...


class ConnectivityResult:
    """Sediment-connectivity index and intermediate raster products.

    Authors
    -------
    B. Gailleton
    """


def extract_channel_heads(dem: Raster, flow: FlowInfo, method: str = "threshold",
                          threshold_pixels: int = 1000,
                          pruning_drainage_area_m2: float = 1000.0,
                          surface_fitting_radius_m: float = 6.0,
                          connected_components_threshold_pixels: int = 100,
                          diagnostics: bool = False) -> ChannelHeads:
    """Extract D8 channel heads by drainage threshold or Wiener curvature.

    Parameters
    ----------
    dem : Raster
        Projected elevation raster.
    flow : FlowInfo
        D8 flow topology on the same grid.
    method : {"threshold", "wiener"}, default="threshold"
        Channel-head detection method. Wiener detection requires FFTW.
    threshold_pixels : int, default=1000
        Drainage-pixel threshold for the threshold method.
    pruning_drainage_area_m2 : float, default=1000
        Minimum retained contributing area for Wiener candidates.
    surface_fitting_radius_m : float, default=6
        Curvature fitting scale for the Wiener method.
    connected_components_threshold_pixels : int, default=100
        Minimum connected candidate size for Wiener extraction.
    diagnostics : bool, default=False
        Retain filtering, curvature, mask, skeleton, and QQ products.

    Returns
    -------
    ChannelHeads
        Channel-head pixels and tabular/geospatial representations.

    Authors
    -------
    B. Gailleton
    """
    ...


class ChannelHeads:
    """Channel-head locations resolved on a D8 topology.

    Attributes
    ----------
    pixels : numpy.ndarray
        Head locations as ``(row, col)`` pairs.
    table : pandas.DataFrame
        Head attributes.
    channel_mask : numpy.ndarray
        Extracted candidate channel domain.
    filtered_dem, tangential_curvature, skeleton, qq_table : optional
        Diagnostics retained when requested.

    Authors
    -------
    B. Gailleton
    """


class ChannelScope:
    """Selected connected D8 channel domain.

    A scope is constructed from :meth:`FlowInfo.channel_scope` or
    :meth:`FlowInfo.scope_from_heads`, then passed to
    :meth:`FlowInfo.channel_network`.

    Attributes
    ----------
    mask : Raster
        Channel-domain mask on the flow grid.
    sources, outlets : numpy.ndarray
        Source and outlet pixel coordinates as ``(row, col)`` pairs.
    table : pandas.DataFrame
        Per-channel-node identifiers and attributes.

    Authors
    -------
    B. Gailleton
    """


class ChannelNetwork:
    """Topology and vector products for a selected D8 channel scope.

    Key methods
    -----------
    chi(...)
        Create a :class:`ChiCollection`.
    points_geodataframe(), lines_geodataframe(), junctions_geodataframe()
        Export node, link, and junction vector layers.
    link_profiles(values)
        Sample a raster along links.
    junction_angles(...), link_bearings(...)
        Calculate network geometry.

    Authors
    -------
    B. Gailleton
    """


def junction_basins(network: ChannelNetwork, selection: str = "all",
                    outlets: Any | None = None) -> JunctionBasinCollection:
    """Select junction-aware basins from a channel network.

    Parameters
    ----------
    network : ChannelNetwork
        Network defining junction and stream-order topology.
    selection : {"all", "outlets", "nested"}, default="all"
        Basin selection mode. The latter two require ``outlets``.
    outlets : array-like of int, optional
        N-by-2 ``(row, col)`` outlet pixels.

    Returns
    -------
    JunctionBasinCollection
        Selected basins with mutable chainable filters.

    Authors
    -------
    B. Gailleton
    """
    ...


class JunctionBasinCollection:
    """Junction-aware basins selected from a :class:`ChannelNetwork`.

    Methods ``filter_incomplete``, ``filter_stream_order``, ``filter_area``,
    ``filter_largest``, and nested-basin filters mutate the selection and
    return it for chaining.

    Attributes
    ----------
    labels : Raster
        Raster of selected basin identifiers.
    table : pandas.DataFrame
        Per-basin statistics and topology fields.

    Authors
    -------
    B. Gailleton
    """


def basin_surface_metrics(basins: Any, dem: Raster, scales: Sequence[float],
                          metrics: Sequence[str], include_outline: bool = True,
                          n_threads: int = 0) -> dict[str, Any]:
    """Calculate requested surface metrics for each basin and scale.

    Parameters
    ----------
    basins : BasinCollection or JunctionBasinCollection
        Basin selection on the DEM grid.
    dem : Raster
        Elevation raster.
    scales : sequence of float
        Surface-fitting scales in metres.
    metrics : sequence of str
        Requested basin metrics.
    include_outline : bool, default=True
        Include perimeter/outline measurements.
    n_threads : int, default=0
        OpenMP thread count; zero uses the runtime default.

    Returns
    -------
    dict[str, object]
        Metric tables and raster/vector products keyed by metric name.

    Authors
    -------
    B. Gailleton
    """
    ...


class ChiCollection:
    """Chi-coordinate channel products for one or more network basins.

    Create through :meth:`ChannelNetwork.chi`. Its :meth:`fit` method performs
    piecewise-linear chi-elevation segmentation.

    Attributes
    ----------
    raster : Raster
        Chi coordinate on the flow grid.
    nodes : pandas.DataFrame
        Channel nodes with basin, chi, and network attributes.

    Authors
    -------
    B. Gailleton
    """


class ChiFit:
    """Piecewise-linear chi-elevation segmentation result.

    Create with ``chi.fit(dem, critical_divergence_metres=...)``.

    Attributes
    ----------
    nodes : pandas.DataFrame
        Channel-node chi, elevation, segment, and fitted values.
    segments : pandas.DataFrame
        One row per fitted segment, including ksn-like slope.
    summary : pandas.DataFrame
        Basin-level fit summary.

    Authors
    -------
    B. Gailleton
    """


class ChiInversion:
    """Regularised ksn/u-star inversion across candidate m/n values.

    Create with :meth:`ChiFit.invert_ksn`.

    Attributes
    ----------
    trials : pandas.DataFrame
        Candidate m/n values and inversion diagnostics.
    summary, nodes : pandas.DataFrame
        Selected inversion results and channel-node values.

    Authors
    -------
    B. Gailleton
    """


def fetch_dem(west: float, south: float, east: float, north: float,
              source: str = "cop30_aws", target_epsg: int = 0,
              grid_spacing: float = 30.0, resampling: str = "cubic",
              nodata_value: float = -9999.0, size_px: int = 0,
              cache_dir: str = "", ot_dem_type: str = "COP30",
              pc_collection: str = "cop-dem-glo-30", pc_asset: str = "data",
              http_timeout_seconds: int = 60,
              http_connect_timeout_seconds: int = 15,
              http_max_retry: int = 3) -> Raster:
    """Fetch a WGS84 bounding box as an in-memory projected DEM.

    Parameters
    ----------
    west, south, east, north : float
        WGS84 bounds in decimal degrees. ``west`` must be less than ``east``.
    source : {"cop30_aws", "opentopography", "planetary_computer"}, default="cop30_aws"
        DEM provider. OpenTopography reads credentials from ``OT_API_KEY`` or
        ``OT_API_KEY_FILE``.
    target_epsg : int, default=0
        Output EPSG code; zero selects an appropriate UTM zone.
    grid_spacing : float, default=30
        Output pixel size in metres.
    resampling : str, default="cubic"
        GDAL resampling algorithm.
    nodata_value : float, default=-9999
        Output nodata value.
    size_px : int, default=0
        Optional maximum output dimension; zero leaves spacing unchanged.
    cache_dir : path-like, default=""
        Download cache directory; empty uses no explicit cache.
    ot_dem_type, pc_collection, pc_asset : str
        Provider-specific dataset selectors.
    http_timeout_seconds, http_connect_timeout_seconds, http_max_retry : int
        HTTP request controls.

    Returns
    -------
    Raster
        Downloaded and reprojected DEM. No file is written.

    Authors
    -------
    B. Gailleton
    """
    ...


def fetch_dem_for_outlets(outlets: Any, source: str = "cop30_aws",
                          target_epsg: int = 0, grid_spacing: float = 30.0,
                          initial_radius_km: float = 100.0,
                          auto_split_disconnected_basins: bool = True,
                          split_cluster_gap_km: float = 25.0,
                          max_final_dem_area_km2: float = 2500.0,
                          split_overlap_km: float = 2.0,
                          max_basins_per_dem: int = 0, cache_dir: str = "",
                          resampling: str = "cubic", nodata_value: float = -9999.0,
                          snapping_method: str = "downstream",
                          contributing_pixels_threshold: int = 1000,
                          edge_buffer_pixels: int = 3, nodata_buffer_pixels: int = 1,
                          tight_bbox_buffer_pixels: int = 3,
                          verify_expand_factor: float = 1.25,
                          verify_max_iters: int = 3, tight_bbox_refetch: bool = True,
                          min_search_radius_pixels: int = 0,
                          max_search_radius_pixels: int = 25,
                          mask_outside_basins: bool = False) -> OutletFetchResult:
    """Fetch and verify DEM groups containing requested outlet basins.

    Parameters
    ----------
    outlets : geopandas.GeoDataFrame
        WGS84 outlet points, optionally with grouping attributes.
    source, target_epsg, grid_spacing, cache_dir, resampling, nodata_value
        DEM-provider and output-grid settings; see :func:`fetch_dem`.
    initial_radius_km : float, default=100
        First fetch-envelope radius around each outlet group.
    auto_split_disconnected_basins : bool, default=True
        Split disconnected groups into suitable DEM requests.
    split_cluster_gap_km, max_final_dem_area_km2, split_overlap_km,
    max_basins_per_dem
        Controls automatic group splitting.
    snapping_method : str, default="downstream"
        Outlet-to-channel snapping strategy.
    contributing_pixels_threshold : int, default=1000
        Minimum contributing pixels defining a channel for snapping.
    edge_buffer_pixels, nodata_buffer_pixels, tight_bbox_buffer_pixels
        Safety buffers used during basin verification and trimming.
    verify_expand_factor : float, default=1.25
        Envelope-growth factor after an unsafe verification.
    verify_max_iters : int, default=3
        Maximum fetch/verify attempts.
    tight_bbox_refetch : bool, default=True
        Re-fetch a tight verified envelope.
    min_search_radius_pixels, max_search_radius_pixels : int
        Channel-search limits for snapping.
    mask_outside_basins : bool, default=False
        Set cells outside selected basins to nodata in returned DEMs.

    Returns
    -------
    OutletFetchResult
        Rasters, groups, original outlets, snapped outlets, and basin labels.

    Authors
    -------
    B. Gailleton
    """
    ...


class OutletFetchResult:
    """In-memory DEM groups fetched around outlet points.

    Attributes
    ----------
    rasters : list of Raster
        One verified DEM per output group.
    groups, outlets, snapped_outlets : geopandas.GeoDataFrame
        Grouping, initial inputs, and verified snapped outlet locations.
    basin_labels : list of Raster
        Basin identifiers aligned with ``rasters``.

    Authors
    -------
    B. Gailleton
    """


def bounds_from_vector(path: str, layer: str = "") -> VectorBounds:
    """Derive WGS84 bounds from polygons or exactly two points in a vector layer.

    Parameters
    ----------
    path : path-like
        OGR-readable vector dataset.
    layer : str, default=""
        Optional layer name.

    Returns
    -------
    VectorBounds
        WGS84 bounds and input metadata. Its ``bbox`` can be unpacked into
        :func:`fetch_dem`.

    Authors
    -------
    B. Gailleton
    """
    ...


class VectorBounds:
    """WGS84 bounds and metadata derived from a vector layer.

    Attributes
    ----------
    bbox : tuple[float, float, float, float]
        ``(west, south, east, north)`` WGS84 bounds.
    ok : bool
        Whether bounds were successfully derived.
    message, warning : str
        Diagnostic text.

    Authors
    -------
    B. Gailleton
    """


def opentopography_datasets() -> list[dict[str, Any]]:
    """Return OpenTopography DEM catalogue records.

    Parameters
    ----------
    None

    Returns
    -------
    list of dict
        Provider dataset metadata.

    Authors
    -------
    B. Gailleton
    """
    ...


def utm_epsg_for_point(longitude: float, latitude: float) -> int:
    """Return the UTM EPSG code containing a WGS84 point.

    Parameters
    ----------
    longitude, latitude : float
        WGS84 coordinates in decimal degrees.

    Returns
    -------
    int
        UTM EPSG code.

    Authors
    -------
    B. Gailleton
    """
    ...


def crn_invert_basin_erosion(flow: FlowInfo, samples: Any, production_scaling: Raster,
                             concentration_column: str = "concentration",
                             nuclide: str = "Be10",
                             muon_scheme: str = "BraucherBorchers",
                             channel_threshold_pixels: int = 1000,
                             snap_radius_pixels: int = 25,
                             topographic_shielding: Raster | None = None,
                             snow_shielding: Raster | None = None,
                             self_shielding: Raster | None = None,
                             rock_density_kg_m3: float = 2650.0) -> Any:
    """Invert basin-mean erosion rates from a point GeoDataFrame.

    Parameters
    ----------
    flow : FlowInfo
        D8 topology on the production-scaling grid.
    samples : geopandas.GeoDataFrame
        Sample points and measured concentration column.
    production_scaling : Raster
        Per-pixel cosmogenic production scaling.
    concentration_column : str, default="concentration"
        Measured concentration column name.
    nuclide : str, default="Be10"
        Cosmogenic nuclide.
    muon_scheme : str, default="BraucherBorchers"
        Muon production parameterisation.
    channel_threshold_pixels, snap_radius_pixels : int
        Channel definition and sample-snapping controls.
    topographic_shielding, snow_shielding, self_shielding : Raster, optional
        Multiplicative shielding rasters.
    rock_density_kg_m3 : float, default=2650
        Rock density.

    Returns
    -------
    geopandas.GeoDataFrame
        Input rows plus snapped outlets, inverted erosion, and status fields.

    Authors
    -------
    B. Gailleton
    """
    ...


def crn_pair_diagnostics(samples: Any, mode: str = "burial",
                         nuclide_fast: str = "Al26", nuclide_slow: str = "Be10",
                         concentration_fast_column: str = "Al26_CONC",
                         concentration_slow_column: str = "Be10_CONC",
                         total_scaling: float = 1.0,
                         muon_scheme: str = "BraucherBorchers",
                         effective_erosion_rate_old_g_cm2_yr: float = 0.01) -> Any:
    """Apply two-nuclide burial or step-change diagnostics to sample rows.

    Parameters
    ----------
    samples : pandas.DataFrame or geopandas.GeoDataFrame
        Sample rows containing fast and slow nuclide concentrations.
    mode : {"burial", "step_change"}, default="burial"
        Diagnostic model.
    nuclide_fast, nuclide_slow : str
        Fast- and slow-decaying nuclide identifiers.
    concentration_fast_column, concentration_slow_column : str
        Concentration-column names.
    total_scaling : float, default=1
        Combined production scaling.
    muon_scheme : str, default="BraucherBorchers"
        Muon production parameterisation.
    effective_erosion_rate_old_g_cm2_yr : float, default=0.01
        Pre-event effective erosion rate.

    Returns
    -------
    pandas.DataFrame or geopandas.GeoDataFrame
        Input rows with diagnostic results and status fields.

    Authors
    -------
    B. Gailleton
    """
    ...


def crn_production_scaling(dem: Raster, crn_parameters_path: str = "",
                           sample_age_yr: float = 0.0, atmosphere: str = "auto",
                           scaling: str = "stone", nuclide: str = "Be10") -> Raster:
    """Calculate a per-pixel cosmogenic production-scaling raster.

    Parameters
    ----------
    dem : Raster
        DEM with metre coordinates and elevation.
    crn_parameters_path : path-like, default=""
        Optional CRN parameter file.
    sample_age_yr : float, default=0
        Sample age before present.
    atmosphere : str, default="auto"
        Atmospheric model selection.
    scaling : str, default="stone"
        Spallation scaling scheme.
    nuclide : str, default="Be10"
        Cosmogenic nuclide.

    Returns
    -------
    Raster
        Dimensionless production-scaling factor.

    Authors
    -------
    B. Gailleton
    """
    ...


def crn_steady_state(erosion_rate: Raster, production_scaling: Raster,
                     nuclide: str = "Be10", muon_scheme: str = "BraucherBorchers",
                     topographic_shielding: Raster | None = None,
                     snow_shielding: Raster | None = None,
                     self_shielding: Raster | None = None) -> Raster:
    """Calculate whole-raster steady-state CRN concentration in atoms/g.

    Parameters
    ----------
    erosion_rate, production_scaling : Raster
        Co-registered erosion-rate and production-scaling rasters.
    nuclide : str, default="Be10"
        Cosmogenic nuclide.
    muon_scheme : str, default="BraucherBorchers"
        Muon production parameterisation.
    topographic_shielding, snow_shielding, self_shielding : Raster, optional
        Multiplicative shielding rasters.

    Returns
    -------
    Raster
        Steady-state concentration in atoms/g.

    Authors
    -------
    B. Gailleton
    """
    ...


def crn_outlet_concentration(concentration: Raster, erosion_rate: Raster,
                             flow: FlowInfo, row: int, col: int) -> float:
    """Calculate erosion-weighted upstream concentration at one outlet pixel.

    Parameters
    ----------
    concentration, erosion_rate : Raster
        Co-registered concentration and erosion-rate rasters.
    flow : FlowInfo
        D8 routing topology.
    row, col : int
        Outlet pixel position.

    Returns
    -------
    float
        Upslope erosion-weighted concentration in atoms/g.

    Authors
    -------
    B. Gailleton
    """
    ...


class CRNColumn:
    """One-dimensional cosmogenic-nuclide column model.

    Set a present-to-past erosion history using ``set_erosion_history`` before
    calling ``transient``. Point ``steady_state`` and ``step_change``
    calculations do not require a history.

    Authors
    -------
    B. Gailleton
    """


class CRNColumnSample:
    """Modeled cosmogenic sample within a :class:`CRNColumn`.

    Attributes
    ----------
    concentration_be10, concentration_al26, concentration_c14 : float
        Concentrations in atoms/g.
    depth_m, effective_depth_g_cm2 : float
        Physical and mass depth.

    Authors
    -------
    B. Gailleton
    """


def fill_hydraulic_surface(surface: Raster, node_state: Raster,
                            boundary_conditions: str = "nnnn",
                            slope_epsilon: float = 1e-8,
                            check_valid_outlet: bool = True,
                            priority_flood_backend: str = "binary") -> HydraulicSurfaceFill:
    """Priority-flood a hydraulic ``z + h`` surface towards selected outlets.

    Parameters
    ----------
    surface : Raster
        Hydraulic surface to project; normally bed elevation plus water depth.
    node_state : Raster
        Node-state code raster: ``0`` excluded/wall, ``1`` active, and ``3``
        selectable outlet.
    boundary_conditions : str, default="nnnn"
        Four-character north/east/south/west boundary code. ``"o"`` opens an
        edge; ``"n"`` closes it; periodic codes are supported by FlowInfo.
        These codes define neighbour connectivity, not outlet locations.
    slope_epsilon : float, default=1e-8
        Positive drainage slope imposed across filled flats.
    check_valid_outlet : bool, default=True
        Require every connected active component to contain a state-3 sink.
    priority_flood_backend : {"binary", "radix"}, default="binary"
        Queue used to project the hydraulic surface.

    Returns
    -------
    HydraulicSurfaceFill
        Filled surface, added depth, and non-conservative added volume.

    Authors
    -------
    B. Gailleton
    """
    ...


class HydraulicSurfaceFill:
    """Priority-flood projection of a hydraulic surface.

    Attributes
    ----------
    filled_surface : Raster
        Outlet-connected projected hydraulic surface.
    added_depth : Raster
        Water depth added by projection.
    added_volume : float
        Non-conservative water volume added during projection.

    Authors
    -------
    B. Gailleton
    """


def run_graphflood(bed: Raster, node_state: Raster, rainfall_rate: Raster,
                   friction_coefficient: Raster, method: str = "mfd",
                   boundary_conditions: str = "oooo", mfd_exponent: float = 1.0,
                   condition_bed: str = "none", fill_slope_epsilon: float = 1e-8,
                   check_valid_outlet: bool = True,
                   friction_depth_exponent: float = 5. / 3.,
                   friction_slope_exponent: float = 0.5, dt: float = 0.001,
                   min_slope: float = 1e-4,
                   point_inflows: Any | None = None, max_steps: int = 10000,
                   run_until_converged: bool = True,
                   convergence_check_every: int = 100,
                   tol_median_abs_dhw: float = 1e-7,
                   tol_median_rel_dhw: float = 1e-5,
                   n_threads: int = 0,
                   initial_water_depth: Raster | None = None) -> HydraulicResult:
    """Run stationary GraphFlood hydraulic routing on a selected domain.

    Parameters
    ----------
    bed : Raster
        Bed-elevation raster.
    node_state : Raster
        Node states: ``0`` excluded/wall, ``1`` active, ``3`` outlet.
    method : {"mfd", "dinfinity"}, default="mfd"
        Flow-partitioning method on the evolving hydraulic surface.
    boundary_conditions : str, default="oooo"
        NESW neighbour-connectivity codes: ``o`` allows flow to leave through
        that edge, ``n`` closes it, and paired ``p`` codes are periodic. These
        codes do not define drainage cells: ``node_state`` does.
    condition_bed : {"none", "fill"}, default="none"
        Optional one-time bed conditioning before routing.
    fill_slope_epsilon, check_valid_outlet
        Hydraulic-surface conditioning controls.
    mfd_exponent : float, default=1
        MFD partitioning exponent.
    rainfall_rate : Raster
        Co-registered rainfall/input-rate field in m/s.
    friction_coefficient : Raster
        Co-registered friction-coefficient field in
        ``q / width = coefficient * h**a * S**b``.
    friction_depth_exponent, friction_slope_exponent : float
        Exponents of the configurable friction law.
    dt, min_slope : float
        Stationary solver step and stability controls.
    point_inflows : array-like, optional
        Point inflow specifications.
    max_steps : int, default=10000
        Maximum solver iterations.
    run_until_converged : bool, default=True
        Stop when median depth-change criteria are met.
    convergence_check_every : int, default=100
        Iteration interval for convergence evaluation.
    tol_median_abs_dhw, tol_median_rel_dhw : float
        Absolute and relative median depth-change tolerances.
    n_threads : int, default=0
        OpenMP thread count; zero uses the runtime default.

    Returns
    -------
    HydraulicResult
        Water-depth, discharge, hydraulic-slope, and convergence products.

    Notes
    -----
    This Python interface is for stationary routing. Hydraulic priority-flood
    projection explicitly adds water to eliminate local minima in ``z + h``;
    its added volume is reported rather than mass-conserved.

    Authors
    -------
    B. Gailleton
    """
    ...


def run_graphflood_optimised(
        bed: Raster, node_state: Raster, rainfall_rate: Raster,
        friction_coefficient: Raster, boundary_conditions: str = "oooo",
        mfd_exponent: float = 1.0, condition_bed: str = "none",
        fill_slope_epsilon: float = 1e-8,
        check_valid_outlet: bool = True,
        friction_depth_exponent: float = 5. / 3.,
        friction_slope_exponent: float = 0.5, dt: float = 1e-3,
        min_slope: float = 1e-4, point_inflows: Any | None = None,
        max_steps: int = 10000, run_until_converged: bool = True,
        convergence_check_every: int = 20,
        tol_median_abs_dhw: float = 1e-7,
        tol_median_rel_dhw: float = 1e-5, reuse_topology: bool = True,
        initial_water_depth: Raster | None = None) -> GraphFloodOptimisedResult:
    """Run the separate CPU-optimised stationary MFD GraphFlood solver.

    This function does not replace :func:`run_graphflood`. It uses a compact
    raster traversal and validates topology before reuse.

    ``performance`` on the returned result contains cumulative stage timings,
    topology rebuild/reuse counts, processed-cell counts, and persistent
    working bytes.

    Authors
    -------
    B. Gailleton
    """
    ...


def run_graphflood_optimised_v2(
        bed: Raster, node_state: Raster, rainfall_rate: Raster,
        friction_coefficient: Raster, boundary_conditions: str = "oooo",
        mfd_exponent: float = 1.0, condition_bed: str = "none",
        fill_slope_epsilon: float = 1e-8,
        check_valid_outlet: bool = True,
        friction_depth_exponent: float = 5. / 3.,
        friction_slope_exponent: float = 0.5, dt: float = 1e-3,
        min_slope: float = 1e-4, point_inflows: Any | None = None,
        max_steps: int = 10000, run_until_converged: bool = True,
        convergence_check_every: int = 20,
        tol_median_abs_dhw: float = 1e-7,
        tol_median_rel_dhw: float = 1e-5, reuse_topology: bool = True,
        dirty_tile_size: int = 64, full_validation_every: int = 128,
        initial_water_depth: Raster | None = None,
        prescribed_qin: Raster | None = None) -> GraphFloodOptimisedResult:
    """Run the incremental v2 CPU GraphFlood solver.

    V2 adds a monotone radix-heap Priority-Flood, a bounds-free interior MFD
    kernel, specialised Manning exponents, and dirty-tile topology validation.
    A failed local validation always falls back to a complete Priority-Flood.

    ``node_state == 2`` marks a flux-only ghost cell. It routes
    ``prescribed_qin`` rather than flow received from its neighbours, while
    its depth and capacity-limited discharge evolve from ``initial_water_depth``.
    ``prescribed_qin=None`` supplies zero.

    Authors
    -------
    B. Gailleton
    """
    ...


def run_graphflood_multigrid(
        bed: Raster, node_state: Raster, rainfall_rate: Raster,
        friction_coefficient: Raster, *, backend: Any = "optimised_v2",
        factors: Any = (8, 4, 2, 1), steps_per_level: Any = 100,
        point_inflows: Any | None = None,
        initial_water_depth: Raster | None = None,
        **backend_options: Any) -> Any:
    """Warm-start GraphFlood over an integer coarse-to-fine grid hierarchy.

    ``backend`` accepts ``"reference"``, ``"optimised"``,
    ``"optimised_v2"``, or a compatible callable. Factors must divide the
    original shape and form a descending integer refinement sequence ending
    at one. The built-in v2 path performs field transfer and depth hand-off
    natively; other backends use the transparent ``RasterHierarchy`` driver.
    All hydraulic fields must be co-registered and valid on active cells.

    Authors
    -------
    B. Gailleton
    """
    ...


def run_graphflood_induced(
        bed: Raster, node_state: Raster, rainfall_rate: Raster,
        friction_coefficient: Raster, qin_threshold: float, *,
        backend: Any = "optimised_v2", halo_cells: int = 3,
        edge_depth_threshold: float = 1e-4, max_expansions: int = 2,
        point_inflows: Any | None = None, **backend_options: Any) -> Any:
    """Route on a discharge-selected channel and floodplain domain.

    A one-step full-domain prepass selects cells by ``qin_threshold``. Flux
    crossing from excluded hillslopes is accumulated as point discharge on
    the retained domain, and wet domain edges can expand by ``halo_cells``.

    Authors
    -------
    B. Gailleton
    """
    ...


class GraphFloodOptimisedPerformance:
    """Cumulative performance counters from the optimised solvers.

    Attributes
    ----------
    topology_check_seconds, projection_seconds, routing_seconds : float
        Cumulative routing-stage wall times.
    friction_seconds, depth_update_seconds, convergence_seconds : float
        Cumulative hydraulic-stage wall times.
    topology_rebuilds, topology_reuses, dirty_tiles_checked : int
        Topology work counters.
    radix_redistributions, interior_routed_cells, boundary_routed_cells : int
        V2 kernel counters; zero for v1 where they do not apply.
    persistent_working_bytes : int
        Capacity-based estimate of retained solver storage.
    """


class GraphFloodOptimisedResult:
    """Raster outputs and diagnostics from the optimised MFD solver."""


class GraphFloodMultigridResult:
    """Final GraphFlood result and results from every grid factor."""


class GraphFloodInducedResult:
    """Reduced-domain result, full prepass, mask, and expansion history."""


class HydraulicResult:
    """Products of stationary GraphFlood routing.

    Attributes
    ----------
    water_depth, water_surface, qin, qout, total_discharge : Raster
        Water and routed-discharge fields.
    hydraulic_slope, specific_discharge, dhw : Raster
        Hydraulic diagnostics.
    converged : bool
        Whether the requested convergence criterion was reached.
    step_count, current_time : int, float
        Solver progress.
    projection_volume, total_water_volume : float
        Hydraulic projection and final water volumes.

    Authors
    -------
    B. Gailleton
    """
