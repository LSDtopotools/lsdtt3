# Fluvial, Valley, Marsh, And Vector Hydrology API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

Valley identification is upstream of centreline tracing. extract_valley_centrelines requires an explicit, co-registered valley_mask.

## `FloodplainThresholdResult`

Floodplain and terrace masks from :func:`classify_floodplain_terraces_by_thresholds`.

``band_index`` is 0 on floodplain, 1 on terrace and -1 elsewhere, the
encoding of the ladder's band index. The counts report what the patch
filters removed.

Returned by ``classify_floodplain_terraces_by_thresholds()``; not constructible directly.

Attributes
----------
band_index : Raster
floodplain_mask : Raster
n_floodplain_candidates : int
n_floodplain_components_kept : int
n_floodplain_dropped_no_channel : int
n_floodplain_dropped_small : int
n_floodplain_pixels : int
n_terrace_components : int
n_terrace_pixels : int
terrace_mask : Raster

Authors
-------
B. Gailleton

### Members

#### `band_index`

(self) -> lsdtt3.Raster

#### `floodplain_mask`

(self) -> lsdtt3.Raster

#### `n_floodplain_candidates`

(self) -> int

#### `n_floodplain_components_kept`

(self) -> int

#### `n_floodplain_dropped_no_channel`

(self) -> int

#### `n_floodplain_dropped_small`

(self) -> int

#### `n_floodplain_pixels`

(self) -> int

#### `n_terrace_components`

(self) -> int

#### `n_terrace_pixels`

(self) -> int

#### `terrace_mask`

(self) -> lsdtt3.Raster

---

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

## `PerimeterSources`

Valley-floor channel sources on the DEM perimeter, from :func:`find_perimeter_sources`.

``pixels`` are the source cells as ``(row, col)``, ready for
``FlowInfo.channel_scope(..., sources=result.pixels)``. ``sources``,
``crossings`` (every gentle run of the perimeter trace, with its verdict),
``outlets`` (the outer outlet runs, by summed drainage area) and
``profile`` (the perimeter trace) are the tables lsdtt-channel-extraction
writes; ``report`` holds the resolved automatic parameters and counts.

Returned by ``find_perimeter_sources()``; not constructible directly.

Attributes
----------
crossings : pandas.DataFrame
    One row per gentle perimeter run: run_id, n_visits, width_metres, min_elevation, relief_metres, low_spot, verdict (too_narrow, not_low_spot, outlet_crossing, failed_inward, drains_elsewhere, too_close_to_outlet, source), and the inward-test results.
outlets : pandas.DataFrame
    One row per outer outlet run, in rank order: outlet_run_id, rank, summed_drainage_area_m2, n_cells, width_metres, on_crossing, run_id, selected, row, col, x, y, elevation (of its largest-area cell).
pixels : numpy.ndarray
    N x 2 int64 array of source (row, col).
profile : pandas.DataFrame
    The perimeter trace, one row per visit: distance_metres, row, col, x, y, elevation, gradient (m/m), gentle, run_id, outlet_cell, outlet_run_id, repeated_visit.
report : dict
    Resolved automatic parameters and counts, as a dict.
sources : pandas.DataFrame
    One row per source: source_id, run_id, outlet_run_id, row, col, x, y (cell centre), elevation, crossing_width_metres, crossing_relief_metres, path_length_to_outlet_metres, drop_to_outlet_metres.

Authors
-------
B. Gailleton

### Members

#### `crossings`

One row per gentle perimeter run: run_id, n_visits, width_metres, min_elevation, relief_metres, low_spot, verdict (too_narrow, not_low_spot, outlet_crossing, failed_inward, drains_elsewhere, too_close_to_outlet, source), and the inward-test results.

#### `outlets`

One row per outer outlet run, in rank order: outlet_run_id, rank, summed_drainage_area_m2, n_cells, width_metres, on_crossing, run_id, selected, row, col, x, y, elevation (of its largest-area cell).

#### `pixels`

N x 2 int64 array of source (row, col).

#### `profile`

The perimeter trace, one row per visit: distance_metres, row, col, x, y, elevation, gradient (m/m), gentle, run_id, outlet_cell, outlet_run_id, repeated_visit.

#### `report`

Resolved automatic parameters and counts, as a dict.

#### `sources`

One row per source: source_id, run_id, outlet_run_id, row, col, x, y (cell centre), elevation, crossing_width_metres, crossing_relief_metres, path_length_to_outlet_metres, drop_to_outlet_metres.

---

## `RiserResult`

Riser (break-of-slope step) detection output from :func:`find_risers`.

``centreline`` is 1 on the kept riser segments; ``segments`` has one row
and line per segment. ``edge_field``, ``adaptive_binary`` and ``skeleton``
are the break-of-slope products the segments were drawn from, and
``gap_closing_mask`` the cells added by gap closing.

Returned by ``find_risers()``; not constructible directly.

Attributes
----------
adaptive_binary : Raster
centreline : Raster
edge_field : Raster
gap_closing_mask : Raster
n_dropped_low_step : int
n_dropped_short : int
n_segments_found : int
segments : geopandas.GeoDataFrame
    One row per kept segment: segment_id, n_pixels, length_metres, step_height_metres (median step across the segment), face_gradient (m/m), mean_elevation_metres, is_loop, and the line geometry.
skeleton : Raster

Authors
-------
B. Gailleton

### Members

#### `adaptive_binary`

(self) -> lsdtt3.Raster

#### `centreline`

(self) -> lsdtt3.Raster

#### `edge_field`

(self) -> lsdtt3.Raster

#### `gap_closing_mask`

(self) -> lsdtt3.Raster

#### `n_dropped_low_step`

(self) -> int

#### `n_dropped_short`

(self) -> int

#### `n_segments_found`

(self) -> int

#### `segments`

One row per kept segment: segment_id, n_pixels, length_metres, step_height_metres (median step across the segment), face_gradient (m/m), mean_elevation_metres, is_loop, and the line geometry.

#### `skeleton`

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

Returned by ``extract_valley_centrelines()``, ``ValleyWidthNetwork.centrelines``; not constructible directly.

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
relative_elevation(dem: Raster, relative_elevation_domain_algorithm_choice: str = 'valley_mask', domain_mask: object | None = None, relative_elevation_max_distance_metres: float = 0.0) -> tuple
    Relative elevation above the valley centrelines' baselines (REM).
terrace_profile(dem: Raster, surface_patch_labels: Raster, relative_elevation_max_distance_metres: float = 0.0, terrace_profile_bin_metres: float = 100.0) -> tuple
    Surface patches along the valley centrelines: the terrace long profile.

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

#### `relative_elevation`

relative_elevation(self, dem: lsdtt3.Raster, relative_elevation_domain_algorithm_choice: str = 'valley_mask', domain_mask: object | None = None, relative_elevation_max_distance_metres: float = 0.0) -> tuple

Relative elevation above the valley centrelines' baselines (REM).

lsdtt-valley-centreline ``write_relative_elevation_raster``. On the valley
floor each cell is measured against the centreline node of the largest
channel whose floor it lies on; elsewhere against the node reached by the
shortest floor-width-scaled ground distance. Needs centrelines from a
``ChannelScope`` or :func:`build_valley_width_network` (the channels are
ranked by drainage area).

Parameters
----------
dem : Raster
    Elevation (m), the DEM the centrelines were made on.
relative_elevation_domain_algorithm_choice : {"valley_mask", "supplied", "whole_dem"}, default="valley_mask"
    Where the field is computed: the valley floor, ``domain_mask``, or
    every valid cell. The floor is always included.
domain_mask : Raster or None, default=None
    Cells > 0.5, for ``supplied``.
relative_elevation_max_distance_metres : float, default=0.0
    Largest ground distance (m) from a cell to its node; 0 means none.

Returns
-------
tuple of Raster
    ``(relative_elevation, owner)``: elevation above the owning node's
    baseline (m) and the owning centreline's ``source_id``; nodata outside
    the domain and the floor.

```python
def relative_elevation(self, dem: lsdtt3.Raster, relative_elevation_domain_algorithm_choice: str = 'valley_mask', domain_mask: object | None = None, relative_elevation_max_distance_metres: float = 0.0) -> tuple
```

#### `table`

One row per channel: source_id, valley_hand_break_metres, n_valley_pixels, routing_succeeded, trough_retries, n_centreline_nodes.

#### `terrace_profile`

terrace_profile(self, dem: lsdtt3.Raster, surface_patch_labels: lsdtt3.Raster, relative_elevation_max_distance_metres: float = 0.0, terrace_profile_bin_metres: float = 100.0) -> tuple

Surface patches along the valley centrelines: the terrace long profile.

lsdtt-valley-centreline ``write_terrace_profile_csv``. Every patch pixel is
tied to a centreline node (ownership computed over the whole DEM) and
binned by distance along that centreline.

Parameters
----------
dem : Raster
    Elevation (m).
surface_patch_labels : Raster
    Patch label per cell (0 = none), for example
    ``FlatSurfaceClassification.banded_labels``.
relative_elevation_max_distance_metres : float, default=0.0
    Largest ground distance (m) from a cell to its node; 0 means none.
terrace_profile_bin_metres : float, default=100.0
    Bin width (m) along the centreline.

Returns
-------
tuple of pandas.DataFrame
    ``(bins, patches)``. ``bins``: centreline_id, patch_label,
    bin_centre_metres, n_pixels, and elevation_* and
    relative_elevation_* five-number summaries (min, p25, median, p75,
    max). ``patches``: patch_label, centreline_id, n_pixels,
    n_pixels_reassigned, along_min_metres, along_max_metres,
    elevation_median, relative_elevation_median.

```python
def terrace_profile(self, dem: lsdtt3.Raster, surface_patch_labels: lsdtt3.Raster, relative_elevation_max_distance_metres: float = 0.0, terrace_profile_bin_metres: float = 100.0) -> tuple
```

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

Methods
-------
levels(ladder_break_algorithm_choice: str = 'antimode', ladder_max_level_count: int = 12, vertical_accuracy_metres: float = 0.0, ladder_resolvability_kappa: float = 3.0, cross_section_level_min_section_fraction: float = 0.5, include_dropped: bool = False) -> object
    Surface levels per channel from every section's hydraulic breaks.

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

#### `levels`

levels(self, ladder_break_algorithm_choice: str = 'antimode', ladder_max_level_count: int = 12, vertical_accuracy_metres: float = 0.0, ladder_resolvability_kappa: float = 3.0, cross_section_level_min_section_fraction: float = 0.5, include_dropped: bool = False) -> object

Surface levels per channel from every section's hydraulic breaks.

lsdtt-valley-centreline ``write_cross_section_levels_csv``: the breaks of
all of a channel's sections are pooled and split by the ladder's rule. The
result can be passed, grouped by ``source_id``, as ``supplied_levels`` to
``classify_flat_surfaces_with_scope``. The breaks are those in
``hydraulic_breaks``.

Parameters
----------
ladder_break_algorithm_choice : {"antimode", "jenks", "gap"}, default="antimode"
    How the pooled breaks are split.
ladder_max_level_count : int, default=12
    Most levels tried.
vertical_accuracy_metres : float, default=0.0
    DEM vertical accuracy E_v (m).
ladder_resolvability_kappa : float, default=3.0
    Levels must be at least kappa x E_v apart.
cross_section_level_min_section_fraction : float, default=0.5
    A level found on fewer than this fraction of the channel's sections is
    dropped.
include_dropped : bool, default=False
    Also return the dropped levels (``kept`` False).

Returns
-------
pandas.DataFrame
    One row per level: source_id, window_centre_metres, level_index,
    level_hand_metres, n_breaks, n_sections, spread_metres, status,
    required_gap_metres, min_level_gap_metres, kept.

```python
def levels(self, ladder_break_algorithm_choice: str = 'antimode', ladder_max_level_count: int = 12, vertical_accuracy_metres: float = 0.0, ladder_resolvability_kappa: float = 3.0, cross_section_level_min_section_fraction: float = 0.5, include_dropped: bool = False) -> object
```

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

## `ValleyWidthNetwork`

A valley centreline for every channel segment of a network, from
:func:`build_valley_width_network`.

``centrelines`` feeds :func:`extract_valley_cross_sections` as usual;
``segments`` describes the tree (each tributary stops at its parent's
channel); ``section_exclusions`` applies the network's rules to the cross
sections.

Returned by ``build_valley_width_network()``; not constructible directly.

Attributes
----------
centrelines : ValleyCentrelineResult
    ValleyCentrelineResult with one channel per segment (source_id = segment_id).
segments : pandas.DataFrame
    One row per segment: segment_id, parent_segment_id (-1 for a trunk), nesting_depth, n_channel_pixels, drains_to_edge, catchment_touches_edge, valley_width_estimate_metres, medial_path_failed, n_samples_snapped_to_mask_edge.

Methods
-------
section_exclusions(sections: ValleyCrossSections, run_tributary_stop_at_parent: bool = True, tributary_stop_width_fraction: float = 1.0, tributary_parent_width_window_multiple: float = 2.0, run_tributary_crossing_removal: bool = True) -> object
    Which cross sections the network rules exclude, and why.

Authors
-------
B. Gailleton

### Members

#### `centrelines`

ValleyCentrelineResult with one channel per segment (source_id = segment_id).

#### `section_exclusions`

section_exclusions(self, sections: lsdtt3.ValleyCrossSections, run_tributary_stop_at_parent: bool = True, tributary_stop_width_fraction: float = 1.0, tributary_parent_width_window_multiple: float = 2.0, run_tributary_crossing_removal: bool = True) -> object

Which cross sections the network rules exclude, and why.

As lsdtt-valley-centreline with ``run_valley_width_network``.

Parameters
----------
sections : ValleyCrossSections
    Sections made from ``centrelines``.
run_tributary_stop_at_parent : bool, default=True
    Exclude a tributary's sections within its parent's valley.
tributary_stop_width_fraction : float, default=1.0
    Multiple of the parent's half-width within which a tributary section
    is excluded.
tributary_parent_width_window_multiple : float, default=2.0
    Window (multiples of the parent's width) the parent's width is taken
    over.
run_tributary_crossing_removal : bool, default=True
    Exclude sections that cross a tributary's channel.

Returns
-------
pandas.DataFrame
    One row per section: source_id, station_id, exclusion (``none``,
    ``not_in_valley``, ``side_not_found``, ``near_parent_channel`` or
    ``tributary_crossing``).

```python
def section_exclusions(self, sections: lsdtt3.ValleyCrossSections, run_tributary_stop_at_parent: bool = True, tributary_stop_width_fraction: float = 1.0, tributary_parent_width_window_multiple: float = 2.0, run_tributary_crossing_removal: bool = True) -> object
```

#### `segments`

One row per segment: segment_id, parent_segment_id (-1 for a trunk), nesting_depth, n_channel_pixels, drains_to_edge, catchment_touches_edge, valley_width_estimate_metres, medial_path_failed, n_samples_snapped_to_mask_edge.

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

## `build_valley_width_network`

Valley centrelines for every channel of a network, for valley widths.

lsdtt-valley-centreline ``run_valley_width_network``: the channel network
is split into segments at confluences (each tributary's trace stops at its
parent), and each segment gets a centreline along the valley's medial axis
or a smoothed channel.

Parameters
----------
dem : Raster
    The routing DEM; its elevations pin the baselines.
scope : ChannelScope
    The channels (every channel of the scope is used).
valley_mask : Raster
    Valley floor (cells > 0.5) on the DEM grid.
network_centreline_algorithm_choice : {"medial_path", "smoothed_channel"}, default="medial_path"
    Medial axis of the valley mask, or the channel smoothed over a valley
    width.
network_centreline_window_margin_multiple : float, default=3.0
    ``medial_path``: walk window padding, as a multiple of the largest
    distance-to-wall at the segment's channel pixels.
network_centreline_smoothing_width_multiple : float, default=1.0
    ``smoothed_channel``: Gaussian sigma as a multiple of the valley width.
min_centreline_nodes : int, default=10
    Segments with fewer channel pixels get no centreline.
centreline_baseline_algorithm_choice : {"projected_monotonic", "channel_crossings"}, default="projected_monotonic"
    How each baseline is pinned to its channel.

Returns
-------
ValleyWidthNetwork

Authors
-------
B. Gailleton

```python
def build_valley_width_network(dem: lsdtt3.Raster, scope: lsdtt3.ChannelScope, valley_mask: lsdtt3.Raster, network_centreline_algorithm_choice: str = 'medial_path', network_centreline_window_margin_multiple: float = 3.0, network_centreline_smoothing_width_multiple: float = 1.0, min_centreline_nodes: int = 10, centreline_baseline_algorithm_choice: str = 'projected_monotonic') -> lsdtt3.ValleyWidthNetwork
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

## `classify_floodplain_terraces_by_thresholds`

Floodplains and terraces from a HAND bound and a gradient bound (Clubb et al. 2017).

The classifier of lsdtt-fluvial-surfaces' ``clubb2017`` method and of
lsdtt-hillslope-channel-coupling's floodplain rule. Floodplain: admitted
cells with HAND below ``max_floodplain_hand_metres`` and gradient below
``max_floodplain_gradient`` (when > 0), plus the channel cells when
``seed_channel_pixels``; 8-connected components smaller than
``min_band_component_pixels`` are dropped, and with
``require_channel_contact`` only components containing a channel cell are
kept. Terraces (``find_terraces``): admitted, non-floodplain cells of
gentle gradient whose ``terrace_field`` lies more than
``min_terrace_height_above_floodplain_metres`` above the floodplain bound
and at most ``max_terrace_height_metres`` (when > 0). Pick the bounds from
the data with :func:`select_threshold`.

Parameters
----------
admitted : Raster
    Cells > 0.5 may be floodplain or terrace; its mask is the output domain.
hand : Raster
    Height above the channel (m), co-registered.
channel_mask : Raster
    Channel cells (> 0.5), for seeding and the contact test.
gradient : Raster or None, default=None
    Gradient (m/m); required when ``max_floodplain_gradient`` > 0.
terrace_field : Raster or None, default=None
    The terrace reference (HAND, or :func:`relief_above_nearest_baseline`);
    required with ``find_terraces``.
max_floodplain_hand_metres : float, default=0.0
    Floodplain HAND bound (m).
max_floodplain_gradient : float, default=0.0
    Floodplain gradient bound (m/m); 0 applies none.
seed_channel_pixels, require_channel_contact : bool, default=True
    Add the channel cells to the floodplain; keep only components that
    contain one.
min_band_component_pixels : int, default=0
    Smallest component kept (pixels); 0 or 1 keeps all.
fill_floodplain_holes : bool, default=False
    Fill holes in the floodplain components.
find_terraces : bool, default=False
    Also classify terraces.
min_terrace_height_above_floodplain_metres : float, default=0.0
    Terraces start this far (m, >= 0) above the floodplain bound.
max_terrace_height_metres : float, default=0.0
    Terrace ceiling (m) of ``terrace_field``; 0 means none.

Returns
-------
FloodplainThresholdResult

Authors
-------
B. Gailleton

```python
def classify_floodplain_terraces_by_thresholds(admitted: lsdtt3.Raster, hand: lsdtt3.Raster, channel_mask: lsdtt3.Raster, gradient: object | None = None, terrace_field: object | None = None, max_floodplain_hand_metres: float = 0.0, max_floodplain_gradient: float = 0.0, seed_channel_pixels: bool = True, require_channel_contact: bool = True, min_band_component_pixels: int = 0, fill_floodplain_holes: bool = False, find_terraces: bool = False, min_terrace_height_above_floodplain_metres: float = 0.0, max_terrace_height_metres: float = 0.0) -> lsdtt3.FloodplainThresholdResult
```

---

## `cluster_surface_levels`

Group surface heights into levels and test whether they cluster.

The test behind lsdtt-fluvial-surfaces' ``write_patch_level_clusters_csv``:
sorted heights are split wherever consecutive values are at least
``required_gap_metres`` apart. The statistic is the empty fraction, the
share of the height range taken by those gaps, compared with ``n_null``
samples of as many heights drawn uniformly over the same range:
``p = (1 + #{null >= observed}) / (1 + n_null)``. A small p says the
heights sit on levels more than random heights would.

Parameters
----------
values : Raster or array-like of float
    Heights (m), for example the ``median_field`` of
    ``FlatSurfaceClassification.patches``.
required_gap_metres : float
    Gap (m, > 0) that separates two levels, typically kappa x E_v.
n_null : int, default=999
    Number of uniform null samples.
seed : int, default=1
    Random seed of the null samples.

Returns
-------
dict
    ``clusters`` (DataFrame: lower_metres, upper_metres, centre_metres (the
    median), n_values, ascending), ``n_values``, ``empty_fraction``,
    ``null_mean_empty_fraction``, ``p_value`` (1 with fewer than 3 values)
    and ``n_null``.

Authors
-------
B. Gailleton

```python
def cluster_surface_levels(values: object, required_gap_metres: float, n_null: int = 999, seed: int = 1) -> dict
```

---

## `compare_surface_classes`

Agreement between the surface classes of several methods, pair by pair.

The comparison lsdtt-fluvial-surfaces writes when it runs several methods.

Parameters
----------
surface_classes : dict of str to Raster
    Method name to a surface-class raster (0 none, 1 floodplain, 2
    terrace; ``FlatSurfaceClassification.surface_class``), all on one grid.
    Cells invalid in any raster are not counted.

Returns
-------
pandas.DataFrame
    One row per pair (method_a before method_b, in dict order):
    ``floodplain_iou``, ``terrace_iou`` (NaN when neither method has the
    class), ``classified_agreement`` (of the cells either method classed,
    the fraction both gave the same class) and ``confusion`` (3 x 3 counts,
    ``[class_a][class_b]``).

Authors
-------
B. Gailleton

```python
def compare_surface_classes(surface_classes: dict) -> object
```

---

## `extract_valley_centrelines`

Route centrelines through a supplied valley mask using a resolved ``ChannelScope``.

This function has 2 call forms:

1. ``extract_valley_centrelines(dem: Raster, scope: ChannelScope, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: Raster, centreline_baseline_algorithm_choice: str = 'projected_monotonic') -> ValleyCentrelineResult``

   Route centrelines through a supplied valley mask using a resolved ``ChannelScope``.

2. ``extract_valley_centrelines(dem: Raster, channels: object, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: Raster, centreline_baseline_algorithm_choice: str = 'projected_monotonic') -> ValleyCentrelineResult``

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
centreline_baseline_algorithm_choice : {"projected_monotonic", "channel_crossings"}, default="projected_monotonic"
    How each centreline's baseline elevation is pinned to its channel, as
    in lsdtt-valley-centreline: channel elevations projected onto the
    centreline and made non-increasing, or pinned only where the centreline
    crosses the channel.

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
def extract_valley_centrelines(dem: lsdtt3.Raster, scope: lsdtt3.ChannelScope, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: lsdtt3.Raster, centreline_baseline_algorithm_choice: str = 'projected_monotonic') -> lsdtt3.ValleyCentrelineResult
def extract_valley_centrelines(dem: lsdtt3.Raster, channels: object, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: lsdtt3.Raster, centreline_baseline_algorithm_choice: str = 'projected_monotonic') -> lsdtt3.ValleyCentrelineResult
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

## `find_perimeter_sources`

Find channel sources where valley floors cross the DEM perimeter.

lsdtt-channel-extraction ``run_perimeter_sources``; each keyword is that
CLI's ``perimeter_*`` parameter without the prefix (``n_outlets`` is
``n_perimeter_outlets``). For a floodplain that crosses the edge of the
DEM, the channel enters from outside, so no head exists inside the DEM.
The perimeter is traced, gentle wide low runs are found, and each is kept
as a source when the path inward from it is gentle and drains to one of
the ``n_outlets`` largest outlets.

Parameters
----------
dem : Raster
    The DEM after nodata handling: its mask is the valid region, and its
    elevations give the profile and gradients.
routed_dem : Raster
    The filled / flat-resolved DEM ``flow`` was built on.
flow : FlowInfo
    D8 routing on ``routed_dem`` with open boundaries on every edge
    (``boundary_conditions="oooo"``).
n_outlets : int, default=1
    Number of outer outlets (largest summed drainage area) sources must
    drain to.
min_valley_width_metres : float, default=300.0
    Narrowest gentle run (m) that counts as a valley-floor crossing.
max_gradient : float, default=0.05
    Largest perimeter gradient (m/m) that counts as gentle.
gradient_radius_metres : float, default=-1.0
    Plane-fit radius (m) of the perimeter gradient; negative is
    max(5 cells, 10 m).
gap_closing_metres : float, default=-1.0
    Steep gaps up to this length (m) inside a gentle run are absorbed;
    negative is max(2 cells, 20 m), 0 turns it off.
min_relief_metres : float, default=5.0
    Rise (m) needed on both sides for a run to be a low spot.
inward_check_distance_metres : float, default=500.0
    Length (m) of the path inward that is checked.
inward_max_gradient : float, default=0.10
    Steepest windowed gradient (m/m) allowed on that path.
inward_gradient_window_metres : float, default=-1.0
    Window (m) of that gradient; negative is max(5 cells, 25 m).
min_source_path_metres : float, default=-1.0
    Shortest path to the outlet (m) a source may have; negative is twice
    ``min_valley_width_metres``.
trim_pixels : int, default=0
    Erode the valid region by this many cells before tracing, moving the
    trace off interpolated edge values; trimmed cells still route flow.

Returns
-------
PerimeterSources

Authors
-------
B. Gailleton

```python
def find_perimeter_sources(dem: lsdtt3.Raster, routed_dem: lsdtt3.Raster, flow: lsdtt3.FlowInfo, n_outlets: int = 1, min_valley_width_metres: float = 300.0, max_gradient: float = 0.05, gradient_radius_metres: float = -1.0, gap_closing_metres: float = -1.0, min_relief_metres: float = 5.0, inward_check_distance_metres: float = 500.0, inward_max_gradient: float = 0.1, inward_gradient_window_metres: float = -1.0, min_source_path_metres: float = -1.0, trim_pixels: int = 0) -> lsdtt3.PerimeterSources
```

---

## `find_risers`

Find risers: short steep steps between gentler surfaces, as open line segments.

lsdtt-feature-finder ``feature_type: find_risers``; each keyword is that
CLI's ``riser_*`` parameter without the prefix. A break-of-slope edge field
is adaptively thresholded and thinned, optionally gap-closed, split at
junctions into segments, and each segment gets a step height
``|2 dz(d) - dz(2d)|`` across it (0 on a planar slope, H across a step of
height H between treads of equal tilt). Segments below
``min_step_height_metres`` or ``min_segment_pixels`` are dropped. Risers
rarely close, so the product is segments, not polygons.

Parameters
----------
dem : Raster
    Elevation (m).
search_mask : Raster or None, default=None
    Cells (> 0.5) where risers are looked for, for example a valley mask;
    None searches every valid cell.
edge_metric : {"vertical_excess_curvature", "curvedness", "maximum_curvature", "multidirectional_hillshade"}, default="vertical_excess_curvature"
    Edge field the lines are drawn on. The first three come from the
    polyfit; ``multidirectional_hillshade`` is the negated mean hillshade
    over ``multidirectional_hillshade_n_directions`` azimuths.
multidirectional_hillshade_n_directions : int, default=8
    Azimuths for ``multidirectional_hillshade``.
polyfit_lambda_metres : float, default=100.0
    Polyfit window (m) of the curvature metrics; also sets the default step
    window.
vec_adaptive_method : {"mean_std", "mean", "gaussian"}, default="mean_std"
    Adaptive threshold of the edge field.
vec_adaptive_std_factor : float, default=0.5
    ``mean_std``: local standard deviations above the local mean.
vec_adaptive_block_size_pixels : int, default=51
    Adaptive window side (odd, >= 3): wider than a riser, narrower than the
    spacing between risers.
vec_adaptive_C : float, default=0.0
    ``mean`` / ``gaussian``: offset below the local mean.
vec_adaptive_gaussian_sigma : float, default=0.0
    ``mean`` / ``gaussian``: pre-blur sigma; 0 means none.
closing_radius_metres : float, default=0.0
    Closing radius (m) of the binary mask before thinning; 0 means none.
binary_min_component_pixels : int, default=1
    Smallest binary component kept before thinning.
max_spur_pixels, max_terminal_loop_pixels : int, default=10
    Longest spur and terminal loop pruned from the skeleton.
min_component_pixels : int, default=1
    Smallest skeleton component kept.
gap_closing_algorithm_choice : {"none", "astar"}, default="none"
    ``astar`` bridges skeleton endpoints within ``connect_search_radius_metres``
    whose elevations differ by at most ``connect_elevation_tolerance_metres``.
connect_search_radius_metres : float, default=200.0
    Bridge search radius (m).
connect_elevation_tolerance_metres : float, default=2.0
    Largest elevation difference (m) a bridge may join; keep it below the
    smallest riser height.
connect_corridor_half_width_pixels : float, default=12.0
    Half-width (pixels) of the corridor a bridge is routed in.
connect_same_component : bool, default=False
    Also bridge endpoints of the same skeleton component.
step_window_radius_metres : float, default=0.0
    Step window d (m); 0 is half the polyfit window, at least three cells.
min_step_height_metres : float, default=1.0
    Smallest step height (m) kept.
min_segment_pixels : int, default=5
    Shortest segment kept (pixels).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
RiserResult

Authors
-------
B. Gailleton

```python
def find_risers(dem: lsdtt3.Raster, search_mask: object | None = None, edge_metric: str = 'vertical_excess_curvature', multidirectional_hillshade_n_directions: int = 8, polyfit_lambda_metres: float = 100.0, vec_adaptive_method: str = 'mean_std', vec_adaptive_std_factor: float = 0.5, vec_adaptive_block_size_pixels: int = 51, vec_adaptive_C: float = 0.0, vec_adaptive_gaussian_sigma: float = 0.0, closing_radius_metres: float = 0.0, binary_min_component_pixels: int = 1, max_spur_pixels: int = 10, max_terminal_loop_pixels: int = 10, min_component_pixels: int = 1, gap_closing_algorithm_choice: str = 'none', connect_search_radius_metres: float = 200.0, connect_elevation_tolerance_metres: float = 2.0, connect_corridor_half_width_pixels: float = 12.0, connect_same_component: bool = False, step_window_radius_metres: float = 0.0, min_step_height_metres: float = 1.0, min_segment_pixels: int = 5, n_threads: int = 0) -> lsdtt3.RiserResult
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

## `relief_above_nearest_baseline`

Height of each cell above its nearest baseline cell (swath baseline relief).

lsdtt-fluvial-surfaces' ``terrace_reference_field_algorithm_choice
swath_baseline`` with the channel as the baseline.

Parameters
----------
dem : Raster
    Elevation (m).
baseline_mask : Raster
    Baseline cells (> 0.5), for example the channel mask.
swath_half_width_metres : float, default=0.0
    Largest distance (m) from the baseline at which relief is computed; 0
    means no limit.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    Elevation minus the elevation of the nearest baseline cell (m).

Authors
-------
B. Gailleton

```python
def relief_above_nearest_baseline(dem: lsdtt3.Raster, baseline_mask: lsdtt3.Raster, swath_half_width_metres: float = 0.0, n_threads: int = 0) -> lsdtt3.Raster
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

## `score_risers_against_mapped`

Score detected risers against mapped terrace edges, cell by cell.

lsdtt-feature-finder's ``mapped_risers_*`` scoring. A detected cell is
matched when a mapped cell lies within the tolerance of it, and a mapped
cell when a detected cell does; only cells inside the surveyed extent count.

Parameters
----------
riser_centreline : Raster
    Detected riser lines (1/0), for example ``RiserResult.centreline``.
mapped_risers_fname : str
    Vector file of mapped terrace-edge lines (burned with all-touched cells).
mapped_risers_extent_fname : str, default=""
    Polygons of the surveyed extent (majority coverage); empty counts the
    whole raster.
mapped_risers_tolerance_metres : float, default=20.0
    Matching distance (m, centre to centre).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
dict
    ``tolerance_metres``, ``n_scored_cells``, ``n_detected``,
    ``n_detected_matched``, ``n_mapped``, ``n_mapped_matched``,
    ``precision``, ``recall`` and ``f1`` (NaN on an empty denominator).

Authors
-------
B. Gailleton

```python
def score_risers_against_mapped(riser_centreline: lsdtt3.Raster, mapped_risers_fname: str, mapped_risers_extent_fname: str = '', mapped_risers_tolerance_metres: float = 20.0, n_threads: int = 0) -> dict
```

---

## `score_surface_classes`

Score a surface classification against mapped floodplains and terraces.

As lsdtt-fluvial-surfaces' ``mapped_surfaces_*`` parameters: the polygons
are burned by cell centre, and only cells inside the surveyed extent (and
valid in ``surface_class``) count.

Parameters
----------
surface_class : Raster
    0 none, 1 floodplain, 2 terrace.
mapped_surfaces_fname : str
    Vector file of mapped polygons classed ``floodplain``, ``terrace`` or
    ``terraceN``.
mapped_surfaces_layer_name : str, default=""
    Its layer; empty reads the first.
mapped_surfaces_class_field_name : str, default="class"
    Its class field.
mapped_surfaces_extent_fname : str, default=""
    Polygons of the surveyed extent (majority coverage); empty counts the
    whole raster as surveyed.

Returns
-------
dict
    ``n_scored``; ``floodplain_precision``, ``floodplain_recall``,
    ``floodplain_iou`` and the same for terrace (NaN on an empty
    denominator); ``confusion`` (3 x 3 counts, ``[mapped][method]``).

Authors
-------
B. Gailleton

```python
def score_surface_classes(surface_class: lsdtt3.Raster, mapped_surfaces_fname: str, mapped_surfaces_layer_name: str = '', mapped_surfaces_class_field_name: str = 'class', mapped_surfaces_extent_fname: str = '') -> dict
```

---

## `select_threshold`

Pick a threshold from a sample, as lsdtt-fluvial-surfaces picks its floodplain bounds.

The routine behind ``floodplain_threshold_algorithm_choice`` (and the
``qq_upper_departure`` floodplain second pass).

Parameters
----------
values : Raster or array-like of float
    The sample; for a Raster, its valid cells. Non-finite values are dropped.
threshold_algorithm_choice : {"fixed", "qq_gaussian", "cdf_divergence", "qq_upper_departure"}, default="qq_gaussian"
    ``qq_gaussian`` is Clubb et al. (2017): the lowest value at which the
    quantile-quantile plot rejoins a Gaussian drawn through two percentiles
    and stays on it. ``cdf_divergence`` fits a family of distributions and
    takes the value where the empirical CDF most exceeds the best fit.
    ``qq_upper_departure`` mirrors ``qq_gaussian`` for the upper tail.
    ``fixed`` returns ``fixed_value``.
fixed_value : float, default=0.0
    The threshold for ``fixed``.
qq_lower_percentile, qq_upper_percentile : float, default=25.0, 75.0
    Percentiles the Gaussian reference line is drawn through.
qq_threshold_fraction : float, default=0.01
    Fraction of the data range within which a quantile counts as on the line.

Returns
-------
dict
    ``threshold``, ``found`` (False when no threshold could be chosen, and
    ``threshold`` is then meaningless), ``n_values``, and for
    ``cdf_divergence`` ``best_family``, ``threshold_percentile`` and
    ``max_cdf_divergence``.

Authors
-------
B. Gailleton

```python
def select_threshold(values: object, threshold_algorithm_choice: str = 'qq_gaussian', fixed_value: float = 0.0, qq_lower_percentile: float = 25.0, qq_upper_percentile: float = 75.0, qq_threshold_fraction: float = 0.01) -> dict
```

---

## `surface_class_agreement`

Which methods gave each cell a class, as a bit field.

Parameters
----------
surface_classes : dict of str to Raster
    As in :func:`compare_surface_classes` (at most 8 methods).
surface_class : {"floodplain", "terrace"}, default="floodplain"
    The class to count.

Returns
-------
Raster
    Bit k set where method k (dict order) gave ``surface_class``; 0 elsewhere.

Authors
-------
B. Gailleton

```python
def surface_class_agreement(surface_classes: dict, surface_class: str = 'floodplain') -> lsdtt3.Raster
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
