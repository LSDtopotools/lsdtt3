# Fluvial, Valley, Marsh, And Vector Hydrology API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

Valley identification is upstream of centreline tracing. extract_valley_centrelines requires an explicit, co-registered valley_mask.

## `MarshPlatformResult`

Results and diagnostics from marsh-platform extraction.

Returned by ``find_marsh_platform()``; not constructible directly.

Attributes
----------
cutoff_elevation_metres : float
cutoff_histogram : dict
dem_relative_relief : Raster
platform_labels : Raster
platform_mask : Raster
platform_provisional : Raster
platform_seed_mask : Raster
scarp_mask : Raster
scarp_order : Raster
scarp_seed_mask : Raster
search_crossover : Raster
search_histogram : dict
search_space : Raster
search_space_threshold : float
slope : Raster
slope_relative_relief : Raster

Authors
-------
B. Gailleton

### Members

#### `cutoff_elevation_metres`

(self) -> float

#### `cutoff_histogram`

(self) -> dict

#### `dem_relative_relief`

(self) -> lsdtt3.Raster

#### `platform_labels`

(self) -> lsdtt3.Raster

#### `platform_mask`

(self) -> lsdtt3.Raster

#### `platform_provisional`

(self) -> lsdtt3.Raster

#### `platform_seed_mask`

(self) -> lsdtt3.Raster

#### `scarp_mask`

(self) -> lsdtt3.Raster

#### `scarp_order`

(self) -> lsdtt3.Raster

#### `scarp_seed_mask`

(self) -> lsdtt3.Raster

#### `search_crossover`

(self) -> lsdtt3.Raster

#### `search_histogram`

(self) -> dict

#### `search_space`

(self) -> lsdtt3.Raster

#### `search_space_threshold`

(self) -> float

#### `slope`

(self) -> lsdtt3.Raster

#### `slope_relative_relief`

(self) -> lsdtt3.Raster

---

## `StreamBurnResult`

DEM conditioned along mapped stream profiles.

Returned by ``burn_streams()``; not constructible directly.

Attributes
----------
profiles : object
    Per-profile monotone burn diagnostics.
raster : Raster
    Copied DEM with channel cells lowered where required.

Authors
-------
B. Gailleton

### Members

#### `profiles`

Per-profile monotone burn diagnostics.

#### `raster`

Copied DEM with channel cells lowered where required.

---

## `TidalChannelResult`

In-memory tidal-channel extraction products.

Returned by ``find_tidal_channels()``; not constructible directly.

Attributes
----------
channel_mask : Raster
connected_skeleton : object
median_rasters : list
residual_rasters : list
scales_metres : object
skeleton : object
skeleton_bridge_mask : object
threshold_masks : list
thresholds_metres : object
upl : object

Authors
-------
B. Gailleton

### Members

#### `channel_mask`

(self) -> lsdtt3.Raster

#### `connected_skeleton`

(self) -> object

#### `median_rasters`

(self) -> list

#### `residual_rasters`

(self) -> list

#### `scales_metres`

(self) -> object

#### `skeleton`

(self) -> object

#### `skeleton_bridge_mask`

(self) -> object

#### `threshold_masks`

(self) -> list

#### `thresholds_metres`

(self) -> object

#### `upl`

(self) -> object

---

## `TidalSkeletonMetricsResult`

Tidal skeleton topology and station metrics.

Returned by ``analyse_tidal_skeleton()``; not constructible directly.

Attributes
----------
downstream_distance : object
endpoints : Raster
hack_order : Raster
junctions : Raster
n_nodes : int
n_sections : int
n_stations : int
nearest_station_ids : object
node_ids : Raster
nodes : object
section_ids : Raster
sections : object
stations : object
strahler_order : Raster

Authors
-------
B. Gailleton

### Members

#### `downstream_distance`

(self) -> object

#### `endpoints`

(self) -> lsdtt3.Raster

#### `hack_order`

(self) -> lsdtt3.Raster

#### `junctions`

(self) -> lsdtt3.Raster

#### `n_nodes`

(self) -> int

#### `n_sections`

(self) -> int

#### `n_stations`

(self) -> int

#### `nearest_station_ids`

(self) -> object

#### `node_ids`

(self) -> lsdtt3.Raster

#### `nodes`

(self) -> object

#### `section_ids`

(self) -> lsdtt3.Raster

#### `sections`

(self) -> object

#### `stations`

(self) -> object

#### `strahler_order`

(self) -> lsdtt3.Raster

---

## `ValleyCentrelineResult`

Valley masks and centreline diagnostics for supplied channel paths.

Made by :func:`extract_valley_centrelines`; pass it to
:func:`extract_valley_cross_sections`.

Returned by ``extract_valley_centrelines()``; not constructible directly.

Attributes
----------
table : pandas.DataFrame
    One row per channel: source_id, valley_hand_break_metres, n_valley_pixels, routing_succeeded, trough_retries, n_centreline_nodes.
valley_mask : Raster
    0/1 valley-floor Raster used for routing.

Methods
-------
centreline_pixels(source_id: int) -> object
    Return the valley centreline of one channel, upstream to downstream.

Authors
-------
B. Gailleton

### Members

#### `centreline_pixels`

centreline_pixels(self, source_id: int) -> object

Return the valley centreline of one channel, upstream to downstream.

Parameters
----------
source_id : int
    Zero-based channel index (``table.source_id``).

Returns
-------
numpy.ndarray
    N x 2 int64 array of ``(row, col)``.

```python
def centreline_pixels(self, source_id: int) -> object
```

#### `table`

One row per channel: source_id, valley_hand_break_metres, n_valley_pixels, routing_succeeded, trough_retries, n_centreline_nodes.

#### `valley_mask`

0/1 valley-floor Raster used for routing.

---

## `ValleyCrossSections`

Valley cross sections made by :func:`extract_valley_cross_sections`.

All tables are pandas DataFrames keyed by ``source_id`` (channel) and
``station_id`` (section); distances are metres.

Returned by ``extract_valley_cross_sections()``; not constructible directly.

Attributes
----------
features : pandas.DataFrame
    Return detected channel crossings, curvature peaks, and hillshade peaks.
hydraulic_breaks : pandas.DataFrame
    Return hydraulic break stages and local heights above the thalweg.
hydraulic_curve : pandas.DataFrame
    Return optional BankfullMapper-style hydraulic geometry curves.
hydraulic_edges : pandas.DataFrame
    Return hydraulic step-boundary picks for every cross section.
n_sections : int
    Number of cross sections.
samples : pandas.DataFrame
    One row per transect sample: source_id, station_id, offset_metres, x, y, elevation, elevation_valid, in_valley, curvature, curvature_valid, hillshade, hillshade_valid, is_channel.
stations : pandas.DataFrame
    One row per section: source_id, station_id, along_distance_metres, valley_width_metres, centre_x, centre_y, left_offset_metres, right_offset_metres, left_found, right_found.

Authors
-------
B. Gailleton

### Members

#### `features`

Return detected channel crossings, curvature peaks, and hillshade peaks.

#### `hydraulic_breaks`

Return hydraulic break stages and local heights above the thalweg.

#### `hydraulic_curve`

Return optional BankfullMapper-style hydraulic geometry curves.

#### `hydraulic_edges`

Return hydraulic step-boundary picks for every cross section.

#### `n_sections`

Number of cross sections.

#### `samples`

One row per transect sample: source_id, station_id, offset_metres, x, y, elevation, elevation_valid, in_valley, curvature, curvature_valid, hillshade, hillshade_valid, is_channel.

#### `stations`

One row per section: source_id, station_id, along_distance_metres, valley_width_metres, centre_x, centre_y, left_offset_metres, right_offset_metres, left_found, right_found.

---

## `ValleyResult`

Broad-valley detection output from :func:`find_valleys`.

``valley_mask`` is 1 on the selected concave-up surface-shape classes;
``shape_class`` uses the codes of :func:`classify_surface_shape`; the
curvatures (1/m) and ``fitted_surface`` (m) come from the smoothing polyfit.

Returned by ``find_valleys()``; not constructible directly.

Attributes
----------
fitted_surface : Raster
gaussian_curvature : Raster
maximum_curvature : Raster
mean_curvature : Raster
minimum_curvature : Raster
shape_class : Raster
valley_mask : Raster

Authors
-------
B. Gailleton

### Members

#### `fitted_surface`

(self) -> lsdtt3.Raster

#### `gaussian_curvature`

(self) -> lsdtt3.Raster

#### `maximum_curvature`

(self) -> lsdtt3.Raster

#### `mean_curvature`

(self) -> lsdtt3.Raster

#### `minimum_curvature`

(self) -> lsdtt3.Raster

#### `shape_class`

(self) -> lsdtt3.Raster

#### `valley_mask`

(self) -> lsdtt3.Raster

---

## `WaterCentrelineResult`

Vector water-area centreline output.

Returned by ``water_centrelines()``; not constructible directly.

Attributes
----------
lines : object
    Centreline, connector, and tributary GeoDataFrame.
summary : dict

Authors
-------
B. Gailleton

### Members

#### `lines`

Centreline, connector, and tributary GeoDataFrame.

#### `summary`

(self) -> dict

---

## `WaterLineRepairResult`

Repaired water-line output and diagnostics.

Returned by ``repair_water_lines()``; not constructible directly.

Attributes
----------
bridges : object
lines : object
summary : dict
tips : object

Authors
-------
B. Gailleton

### Members

#### `bridges`

(self) -> object

#### `lines`

(self) -> object

#### `summary`

(self) -> dict

#### `tips`

(self) -> object

---

## `analyse_tidal_skeleton`

Analyse an existing tidal-channel skeleton into sections, orders, and station metrics.

Parameters
----------
dem : Raster
    Elevation raster (m).
channel_mask : Raster
    0/1 channel mask on the same grid (``TidalChannelResult.channel_mask``).
skeleton : Raster
    0/1 skeleton on the same grid (``TidalChannelResult.skeleton``).
outlet_mask : Raster or None, default=None
    0/1 outlet mask, needed for ``outlet_mode="outlet_mask"``.
outlet_mode : {"edge_nodes", "endpoint_nodes", "outlet_mask", "lowest_endpoint"}, default="edge_nodes"
    Outlets for downstream distances: skeleton end points touching the
    raster edge, all end points, the masked nodes, or the single lowest
    end point.
station_spacing_metres : float, default=5.0
    Spacing (m) of measurement stations along each section (positive).
compute_cross_sections : bool, default=True
    Cast cross sections at each station to measure width.
tangent_support_points : int, default=2
    Neighbouring stations used to estimate the local channel direction.
cross_section_half_length_metres : float, default=25.0
    Half-length (m) of each cross section.
cross_section_sample_spacing_metres : float, default=1.0
    Sample spacing (m) along a cross section.
max_width_ratio : float, default=3.0
    Reject cross sections implausibly long relative to local width; <= 0
    disables.
compute_station_strip_metrics : bool, default=True
    Assign every pixel to its nearest station and compute strip metrics.
strip_values : dict of str to Raster or None, default=None
    Extra co-registered rasters integrated over the station strips.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
TidalSkeletonMetricsResult

Authors
-------
B. Gailleton

```python
def analyse_tidal_skeleton(dem: lsdtt3.Raster, channel_mask: lsdtt3.Raster, skeleton: lsdtt3.Raster, outlet_mask: object | None = None, outlet_mode: str = 'edge_nodes', station_spacing_metres: float = 5.0, compute_cross_sections: bool = True, tangent_support_points: int = 2, cross_section_half_length_metres: float = 25.0, cross_section_sample_spacing_metres: float = 1.0, max_width_ratio: float = 3.0, compute_station_strip_metrics: bool = True, strip_values: object | None = None, n_threads: int = 0) -> lsdtt3.TidalSkeletonMetricsResult
```

---

## `burn_streams`

Condition a copied DEM along mapped stream lines with monotone profiles.

Parameters
----------
dem : Raster
    Elevation raster.
lines : geopandas.GeoDataFrame
    Stream LineStrings in the DEM CRS.
n_threads : int, default=0
    OpenMP thread count; zero uses the runtime default.

Returns
-------
StreamBurnResult
    Conditioned raster and profile diagnostics.

Notes
-----
Only cells whose monotone stream elevation is lower are changed.

Authors
-------
B. Gailleton

```python
def burn_streams(dem: lsdtt3.Raster, lines: object, n_threads: int = 0) -> lsdtt3.StreamBurnResult
```

---

## `extract_valley_centrelines`

Route centrelines through a supplied valley mask using a resolved ``ChannelScope``.

This function has 2 call forms:

1. ``extract_valley_centrelines(dem: Raster, scope: ChannelScope, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: Raster) -> ValleyCentrelineResult``

   Route centrelines through a supplied valley mask using a resolved ``ChannelScope``.

2. ``extract_valley_centrelines(dem: Raster, channels: object, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: Raster) -> ValleyCentrelineResult``

   Route centrelines through a supplied valley mask from ordered channel paths.

   ``channels`` is a sequence of upstream-to-downstream ``(row, col)`` arrays.
   A synthetic trough is carved into the valley floor, deepest along the
   valley axis, and flow is routed down it, so the centreline follows the
   middle of the valley rather than the channel. This function has two call
   forms: pass a ``ChannelScope`` as ``scope`` (its traced channels are used) or
   the ordered channel paths as ``channels``.

Parameters
----------
dem : Raster
    Elevation raster (m).
scope : ChannelScope
    First call form: channels resolved by ``FlowInfo.channel_scope``.
channels : sequence of numpy.ndarray
    Second call form: N x 2 ``(row, col)`` arrays, upstream to downstream.
corridor_width_metres : float, default=5000.0
    Half-width (m) of the corridor around each channel searched for valley
    pixels.
fallback_valley_hand_metres : float or None, default=None
    Deprecated and ignored when ``valley_mask`` is supplied (a
    ``DeprecationWarning`` is emitted if set).
valley_dilation_pixels : int, default=4
    Pixels by which valley patches are dilated (towards the channel only) to
    fill holes.
trough_scaling_factor : float, default=0.1
    Gradient (m/m) of the synthetic trough carved towards the valley axis:
    the trough is lowered by this times the distance to the valley wall
    (Clubb et al. 2022 method); doubled on retries.
centreline_loops : int, default=5
    Number of trough-carving and routing iterations.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
valley_mask : Raster, keyword-only
    0/1 valley-floor raster on the DEM grid, for example from
    ``find_valleys`` or ``valley_mask_from_cross_sections``.

Returns
-------
ValleyCentrelineResult

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `fallback_valley_hand_metres` -> `valley_mask`; earliest removal: following minor.

```python
def extract_valley_centrelines(dem: lsdtt3.Raster, scope: lsdtt3.ChannelScope, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: lsdtt3.Raster) -> lsdtt3.ValleyCentrelineResult
def extract_valley_centrelines(dem: lsdtt3.Raster, channels: object, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: lsdtt3.Raster) -> lsdtt3.ValleyCentrelineResult
```

---

## `extract_valley_cross_sections`

Extract straight or boundary-fitted valley cross sections.

Optional curvature, hillshade, and channel-mask rasters enable feature picks.
Hydraulic break picks are always available; set ``compute_hydraulic_depth`` to
retain the full per-stage hydraulic-geometry curve.

Stations are placed along each centreline and a transect is cast across the
valley at each; all transects of a channel have the same length (a multiple
of its widest valley width).

Parameters
----------
dem : Raster
    Elevation raster (m).
centrelines : ValleyCentrelineResult
    From :func:`extract_valley_centrelines`.
spacing_metres : float, default=500.0
    Spacing (m) of stations along the centreline (positive).
half_length_multiplier : float, default=1.5
    Transect total length as a multiple of the channel's widest bounded
    valley width (positive).
min_half_length_metres : float, default=0.0
    Floor (m) on the transect half-length; 0 means none.
max_half_length_metres : float, default=5000.0
    Cap (m) on the outward search for a valley side and on the half-length.
tangent_window_nodes : int, default=8
    Centreline nodes on each side used to smooth the local direction.
sample_step_metres : float, default=0.0
    Sample spacing (m) along a transect; 0 uses the cell size.
curved_sections : bool, default=False
    Trace boundary-fitted sections that curve through bends instead of
    straight transects.
curvature, hillshade : Raster or None, default=None
    Optional rasters sampled along transects for feature picks.
channel_mask : Raster or None, default=None
    Optional 0/1 channel raster for channel-crossing picks.
side_gap_tolerance_steps : int, default=1
    Out-of-valley steps tolerated before a valley side is declared.
feature_peak_min_prominence_fraction : float, default=0.15
    Minimum peak prominence, as a fraction of the profile range.
feature_peak_edge_margin_samples : int, default=2
    Ignore peaks this close (samples) to the transect ends.
clip_overlapping_sections : bool, default=True
    Clip transects where they would cross neighbouring sections.
neighbour_clip_window : int, default=2
    Number of neighbouring stations on each side considered for clipping.
min_section_reach_fraction : float, default=0.35
    Drop a station whose clipped length is below this fraction of the full
    transect length.
curved_wall_clearance_fraction : float, default=0.10
    Curved sections only: a clipped side still extends at least this
    fraction beyond the valley wall.
curved_streamwise_smoothing_passes : int, default=40
    Curved sections only: smoothing passes of the along-valley coordinate.
curved_coarse_spacing_metres : float, default=700.0
    Curved sections only: spacing (m) of the first, coarse set of stations
    before midpoints are added.
hydraulic_stage_step_metres : float, default=0.1
    Flood-stage increment (m) used to build width against stage.
hydraulic_min_step_height_metres : float, default=2.0
    Minimum height (m) between accepted hydraulic step boundaries.
hydraulic_min_step_width_pixels : int, default=3
    Minimum width change (pixels) between accepted step boundaries.
hydraulic_midpoint_search_samples : int, default=81
    Search resolution for placing each step midpoint.
hydraulic_min_crossings_preferred : int, default=2
    Preferred number of profile crossings at a step midpoint (2 = one bank
    pair).
compute_hydraulic_depth : bool, default=False
    Store the full stage / width / area / hydraulic-depth curve of each
    section (``hydraulic_curve``; Delchiaro et al. 2025 BankfullMapper).

Returns
-------
ValleyCrossSections

Authors
-------
B. Gailleton

```python
def extract_valley_cross_sections(dem: lsdtt3.Raster, centrelines: lsdtt3.ValleyCentrelineResult, spacing_metres: float = 500.0, half_length_multiplier: float = 1.5, min_half_length_metres: float = 0.0, max_half_length_metres: float = 5000.0, tangent_window_nodes: int = 8, sample_step_metres: float = 0.0, curved_sections: bool = False, curvature: object | None = None, hillshade: object | None = None, channel_mask: object | None = None, side_gap_tolerance_steps: int = 1, feature_peak_min_prominence_fraction: float = 0.15, feature_peak_edge_margin_samples: int = 2, clip_overlapping_sections: bool = True, neighbour_clip_window: int = 2, min_section_reach_fraction: float = 0.35, curved_wall_clearance_fraction: float = 0.1, curved_streamwise_smoothing_passes: int = 40, curved_coarse_spacing_metres: float = 700.0, hydraulic_stage_step_metres: float = 0.1, hydraulic_min_step_height_metres: float = 2.0, hydraulic_min_step_width_pixels: int = 3, hydraulic_midpoint_search_samples: int = 81, hydraulic_min_crossings_preferred: int = 2, compute_hydraulic_depth: bool = False) -> lsdtt3.ValleyCrossSections
```

---

## `find_marsh_platform`

Extract marsh platforms and scarps from a projected DEM.

Scarps (the steep edges between marsh platform and tidal flat) are seeded and
grown from slope and relative relief, then platforms are grown from seeds
above them and cut at an elevation chosen from the elevation histogram.

Parameters
----------
dem : Raster
    Elevation raster (m) of the marsh.
slope : Raster or None, default=None
    Gradient (m/m) on the same grid; None computes it.
slope_method : {"gradient_3x3", "polyfit"}, default="gradient_3x3"
    How the slope is computed when ``slope`` is None.
slope_polyfit_lambda_metres : float, default=5.0
    Polyfit length scale (m) for ``slope_method="polyfit"``.
histogram_n_bins : int, default=100
    Number of histogram bins used to choose thresholds.
search_space_histogram_derivative_threshold : float, default=-2.0
    Histogram-derivative value that sets the search-space threshold (the
    first bin where the derivative rises above it).
max_scarp_growth_iterations : int, default=50
    Iterations of scarp growth after seeding.
max_platform_growth_iterations : int, default=100
    Iterations of platform growth.
min_platform_component_pixels : int, default=1
    Smallest platform component kept (pixels).
preserve_scarps_in_final_platform : bool, default=True
    Include the cleaned scarps in the final platform mask.
apply_scarp_closing : bool, default=False
    Close small gaps in the scarp mask before platform seeding.
scarp_closing_radius_metres : float, default=0.0
    Radius (m) of that closing.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
MarshPlatformResult

Authors
-------
B. Gailleton

```python
def find_marsh_platform(dem: lsdtt3.Raster, slope: object | None = None, slope_method: str = 'gradient_3x3', slope_polyfit_lambda_metres: float = 5.0, histogram_n_bins: int = 100, search_space_histogram_derivative_threshold: float = -2.0, max_scarp_growth_iterations: int = 50, max_platform_growth_iterations: int = 100, min_platform_component_pixels: int = 1, preserve_scarps_in_final_platform: bool = True, apply_scarp_closing: bool = False, scarp_closing_radius_metres: float = 0.0, n_threads: int = 0) -> lsdtt3.MarshPlatformResult
```

---

## `find_tidal_channels`

Extract tidal channels as multiscale local topographic lows or a flood mask.

At each scale the residual ``moving-window median(dem) - dem`` is computed;
pixels lying at least the residual threshold below their local median are
channel (a TidalGeoPro-style channelness signal). Optional cleanup,
skeletonisation and gap bridging follow.

Parameters
----------
dem : Raster
    Elevation raster (m) of a marsh or tidal flat.
use_flood_mask : bool, default=False
    Skip the median workflow: pixels with elevation <= the flood elevation
    are channel.
flood_elevation_metres : float, default=0.0
    Flood elevation (m) for ``use_flood_mask``.
median_window_radii_metres : list of float, default=[5.0, 10.0, 20.0]
    Radii (m) of the median windows, one per scale.
residual_thresholds_metres : list of float, default=[0.1]
    Residual (m) a pixel must reach below the median; one value for all
    scales or one per scale.
scale_mode : {"any", "all"}, default="any"
    A pixel is channel if it passes at any scale or at every scale.
apply_binary_closing : bool, default=False
    Close the channel mask with a disk.
binary_closing_radius_metres : float, default=0.0
    Radius (m) of that disk.
fill_holes : bool, default=False
    Fill enclosed holes in the channel mask.
min_channel_component_pixels : int, default=1
    Smallest channel component kept (pixels).
extract_skeleton : bool, default=True
    Thin the channel mask to a one-pixel skeleton.
skeleton_spur_max_pixels : int, default=10
    Remove skeleton spurs up to this length (pixels).
skeleton_terminal_loop_max_pixels : int, default=10
    Remove terminal skeleton loops up to this size (pixels).
skeleton_min_component_pixels : int, default=1
    Smallest skeleton component kept (pixels).
connect_skeleton_segments : bool, default=False
    Bridge gaps between facing skeleton end points.
connect_max_gap_pixels : float, default=20.0
    Largest gap (pixels) bridged.
connect_max_facing_angle_degrees : float, default=35.0
    Largest angle (degrees) between an end point's heading and the gap.
connect_max_heading_misalignment_degrees : float, default=45.0
    Largest misalignment (degrees) between the two end headings.
connect_max_candidate_score : float, default=1.0
    Largest combined (normalised) score of an accepted bridge.
connect_route_corridor_half_width_pixels : float, default=12.0
    Half-width (pixels) of the corridor in which a bridge is routed.
compute_upl : bool, default=True
    Compute ``upl``: the Euclidean distance (m) from each pixel to the
    nearest channel pixel.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
TidalChannelResult

Authors
-------
B. Gailleton

```python
def find_tidal_channels(dem: lsdtt3.Raster, use_flood_mask: bool = False, flood_elevation_metres: float = 0.0, median_window_radii_metres: collections.abc.Sequence[float] = [5.0, 10.0, 20.0], residual_thresholds_metres: collections.abc.Sequence[float] = [0.1], scale_mode: str = 'any', apply_binary_closing: bool = False, binary_closing_radius_metres: float = 0.0, fill_holes: bool = False, min_channel_component_pixels: int = 1, extract_skeleton: bool = True, skeleton_spur_max_pixels: int = 10, skeleton_terminal_loop_max_pixels: int = 10, skeleton_min_component_pixels: int = 1, connect_skeleton_segments: bool = False, connect_max_gap_pixels: float = 20.0, connect_max_facing_angle_degrees: float = 35.0, connect_max_heading_misalignment_degrees: float = 45.0, connect_max_candidate_score: float = 1.0, connect_route_corridor_half_width_pixels: float = 12.0, compute_upl: bool = True, n_threads: int = 0) -> lsdtt3.TidalChannelResult
```

---

## `find_valleys`

Find broad concave valley forms from a smoothed surface-shape classification.

The DEM is smoothed with a polyfit at ``smoothing_window_metres``, its
principal curvatures are classified with :func:`classify_surface_shape`,
and the selected concave-up classes form the valley mask. Larger windows
keep only broad valleys.

Parameters
----------
raster : Raster
    Elevation raster (m).
smoothing_window_metres : float, default=50.0
    Polyfit length scale (m) used for smoothing.
curvature_threshold : float, default=0.0
    Curvature threshold kt (1/m) of the classification.
include_basin, include_synformal_saddle, include_synform : bool, default=True
    Which concave-up classes (-3 basin, -2 synformal saddle, -1 synform)
    count as valley.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
ValleyResult

Authors
-------
B. Gailleton

```python
def find_valleys(raster: lsdtt3.Raster, smoothing_window_metres: float = 50.0, curvature_threshold: float = 0.0, include_basin: bool = True, include_synformal_saddle: bool = True, include_synform: bool = True, n_threads: int = 0) -> lsdtt3.ValleyResult
```

---

## `fit_valley_baseline`

Interpolate a non-increasing valley baseline through finite elevation pins.

Use ``numpy.nan`` for unpinned centreline nodes.

Parameters
----------
along_distance : array-like of float
    Distance (m) of each centreline node along the centreline.
pin_elevation : array-like of float
    Elevation (m) pinned at each node, NaN where unpinned.

Returns
-------
numpy.ndarray
    Baseline elevation (m) at every node, never increasing downstream.

Authors
-------
B. Gailleton

```python
def fit_valley_baseline(along_distance: object, pin_elevation: object) -> list[float]
```

---

## `repair_water_lines`

Repair short broken water-line gaps, preferentially across roads.

Parameters
----------
water_lines, roads : geopandas.GeoDataFrame
    Projected line layers in the same CRS.
snap_tolerance_metres : float, default=0.5
    Distance (m) within which line ends are treated as connected.
max_gap_length_metres : float, default=25
    Maximum road-associated bridge length (m).
max_nonroad_gap_metres : float, default=30
    Maximum length (m) of a bridge that does not cross a road.
road_intersection_tolerance_metres : float, default=1
    Distance (m) within which a bridge counts as crossing a road.
require_mutual_best : bool, default=True
    Only bridge two tips that are each other's best candidate.

Returns
-------
WaterLineRepairResult
    Original linework plus accepted bridge segments and diagnostics.

Authors
-------
B. Gailleton

```python
def repair_water_lines(water_lines: object, roads: object, snap_tolerance_metres: float = 0.5, max_gap_length_metres: float = 25.0, max_nonroad_gap_metres: float = 30.0, road_intersection_tolerance_metres: float = 1.0, require_mutual_best: bool = True) -> lsdtt3.WaterLineRepairResult
```

---

## `valley_mask_from_cross_sections`

Build a valley-floor mask from the median lowest hydraulic break per channel.

This complements the corridor/HAND mask returned by ``extract_valley_centrelines``.

Parameters
----------
dem : Raster
    Elevation raster (m).
centrelines : ValleyCentrelineResult
    Centrelines the sections were cut from.
sections : ValleyCrossSections
    From :func:`extract_valley_cross_sections`.

Returns
-------
Raster
    0/1 valley-floor mask.

Authors
-------
B. Gailleton

```python
def valley_mask_from_cross_sections(dem: lsdtt3.Raster, centrelines: lsdtt3.ValleyCentrelineResult, sections: lsdtt3.ValleyCrossSections) -> lsdtt3.Raster
```

---

## `water_centrelines`

Extract central flow paths from water polygons and connect tributary lines.

Parameters
----------
areas, tributaries : geopandas.GeoDataFrame
    Projected water polygons and linework in the same CRS.
epsg : int, default=0
    Output CRS identifier; zero infers it from areas.
cell_size_metres : float, default=1
    Internal skeletonisation-grid cell size (m).
padding_metres : float, default=5
    Margin (m) added around the polygons on the internal grid.
tributary_snap_tolerance_metres : float, default=2
    Distance (m) within which tributary ends snap to a water polygon.
centreline_snap_max_metres : float, default=100
    Largest distance (m) a tributary is extended to reach the centreline.
simplify_tolerance_metres : float, default=1
    Line simplification tolerance (m).
min_component_length_metres : float, default=10
    Centreline pieces shorter than this (m) are dropped.

Returns
-------
WaterCentrelineResult
    Connected linework and count diagnostics.

Authors
-------
B. Gailleton

```python
def water_centrelines(areas: object, tributaries: object, epsg: int = 0, cell_size_metres: float = 1.0, padding_metres: float = 5.0, tributary_snap_tolerance_metres: float = 2.0, centreline_snap_max_metres: float = 100.0, simplify_tolerance_metres: float = 1.0, min_component_length_metres: float = 10.0) -> lsdtt3.WaterCentrelineResult
```
