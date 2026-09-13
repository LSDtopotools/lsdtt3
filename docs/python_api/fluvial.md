# Fluvial, Valley, Marsh, And Vector Hydrology API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

Valley identification is upstream of centreline tracing. extract_valley_centrelines requires an explicit, co-registered valley_mask.

## `MarshPlatformResult`

Results and diagnostics from marsh-platform extraction.

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
Results and diagnostics from marsh-platform extraction.

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
    DEM conditioned along mapped stream profiles.

Authors
-------
B. Gailleton

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
In-memory tidal-channel extraction products.

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
Tidal skeleton topology and station metrics.

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
Valley masks and centreline diagnostics for supplied channel paths.

Authors
-------
B. Gailleton

### Members

#### `centreline_pixels`

centreline_pixels(self, source_id: int) -> object

```python
def centreline_pixels(self, source_id: int) -> object
```

#### `table`

(self) -> object

#### `valley_mask`

(self) -> lsdtt3.Raster

---

## `ValleyCrossSections`

No additional binding notes are available.

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
No additional binding notes are available.

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

(self) -> int

#### `samples`

(self) -> object

#### `stations`

(self) -> object

---

## `ValleyResult`

No additional binding notes are available.

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
No additional binding notes are available.

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
    Vector water-area centreline output.

Authors
-------
B. Gailleton

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
    Repaired water-line output and diagnostics.

Authors
-------
B. Gailleton

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

analyse_tidal_skeleton(dem: lsdtt3.Raster, channel_mask: lsdtt3.Raster, skeleton: lsdtt3.Raster, outlet_mask: object | None = None, outlet_mode: str = 'edge_nodes', station_spacing_metres: float = 5.0, compute_cross_sections: bool = True, tangent_support_points: int = 2, cross_section_half_length_metres: float = 25.0, cross_section_sample_spacing_metres: float = 1.0, max_width_ratio: float = 3.0, compute_station_strip_metrics: bool = True, strip_values: object | None = None, n_threads: int = 0) -> lsdtt3.TidalSkeletonMetricsResult

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
    analyse_tidal_skeleton(dem: lsdtt3.Raster, channel_mask: lsdtt3.Raster, skeleton: lsdtt3.Raster, outlet_mask: object | None = None, outlet_mode: str = 'edge_nodes', station_spacing_metres: float = 5.0, compute_cross_sections: bool = True, tangent_support_points: int = 2, cross_section_half_length_metres: float = 25.0, cross_section_sample_spacing_metres: float = 1.0, max_width_ratio: float = 3.0, compute_station_strip_metrics: bool = True, strip_values: object | None = None, n_threads: int = 0) -> lsdtt3.TidalSkeletonMetricsResult

Analyse an existing tidal-channel skeleton into sections, orders, and station metrics.

    Authors
    -------
    B. Gailleton


```python
def analyse_tidal_skeleton(dem: lsdtt3.Raster, channel_mask: lsdtt3.Raster, skeleton: lsdtt3.Raster, outlet_mask: object | None = None, outlet_mode: str = 'edge_nodes', station_spacing_metres: float = 5.0, compute_cross_sections: bool = True, tangent_support_points: int = 2, cross_section_half_length_metres: float = 25.0, cross_section_sample_spacing_metres: float = 1.0, max_width_ratio: float = 3.0, compute_station_strip_metrics: bool = True, strip_values: object | None = None, n_threads: int = 0) -> lsdtt3.TidalSkeletonMetricsResult
```

---

## `burn_streams`

burn_streams(dem: lsdtt3.Raster, lines: object, n_threads: int = 0) -> lsdtt3.StreamBurnResult

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
    burn_streams(dem: lsdtt3.Raster, lines: object, n_threads: int = 0) -> lsdtt3.StreamBurnResult

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

    Authors
    -------
    B. Gailleton


```python
def burn_streams(dem: lsdtt3.Raster, lines: object, n_threads: int = 0) -> lsdtt3.StreamBurnResult
```

---

## `extract_valley_centrelines`

extract_valley_centrelines(dem: lsdtt3.Raster, scope: lsdtt3.ChannelScope, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: lsdtt3.Raster) -> lsdtt3.ValleyCentrelineResult

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
    extract_valley_centrelines(dem: lsdtt3.Raster, scope: lsdtt3.ChannelScope, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: lsdtt3.Raster) -> lsdtt3.ValleyCentrelineResult
extract_valley_centrelines(dem: lsdtt3.Raster, channels: object, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: lsdtt3.Raster) -> lsdtt3.ValleyCentrelineResult

Overloaded function.

1. ``extract_valley_centrelines(dem: lsdtt3.Raster, scope: lsdtt3.ChannelScope, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: lsdtt3.Raster) -> lsdtt3.ValleyCentrelineResult``

Route centrelines through a supplied valley mask using a resolved ``ChannelScope``.


2. ``extract_valley_centrelines(dem: lsdtt3.Raster, channels: object, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: lsdtt3.Raster) -> lsdtt3.ValleyCentrelineResult``

Route centrelines through a supplied valley mask from ordered channel paths.

``channels`` is a sequence of upstream-to-downstream ``(row, col)`` arrays.

    Authors
    -------
    B. Gailleton


**Deprecated compatibility names:**
- `fallback_valley_hand_metres` -> `valley_mask`; earliest removal: following minor.

```python
def extract_valley_centrelines(dem: lsdtt3.Raster, scope: lsdtt3.ChannelScope, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: lsdtt3.Raster) -> lsdtt3.ValleyCentrelineResult
def extract_valley_centrelines(dem: lsdtt3.Raster, channels: object, corridor_width_metres: float = 5000.0, fallback_valley_hand_metres: object | None = None, valley_dilation_pixels: int = 4, trough_scaling_factor: float = 0.1, centreline_loops: int = 5, n_threads: int = 0, *, valley_mask: lsdtt3.Raster) -> lsdtt3.ValleyCentrelineResult
```

---

## `extract_valley_cross_sections`

extract_valley_cross_sections(dem: lsdtt3.Raster, centrelines: lsdtt3.ValleyCentrelineResult, spacing_metres: float = 500.0, half_length_multiplier: float = 1.5, min_half_length_metres: float = 0.0, max_half_length_metres: float = 5000.0, tangent_window_nodes: int = 8, sample_step_metres: float = 0.0, curved_sections: bool = False, curvature: object | None = None, hillshade: object | None = None, channel_mask: object | None = None, side_gap_tolerance_steps: int = 1, feature_peak_min_prominence_fraction: float = 0.15, feature_peak_edge_margin_samples: int = 2, clip_overlapping_sections: bool = True, neighbour_clip_window: int = 2, min_section_reach_fraction: float = 0.35, curved_wall_clearance_fraction: float = 0.1, curved_streamwise_smoothing_passes: int = 40, curved_coarse_spacing_metres: float = 700.0, hydraulic_stage_step_metres: float = 0.1, hydraulic_min_step_height_metres: float = 2.0, hydraulic_min_step_width_pixels: int = 3, hydraulic_midpoint_search_samples: int = 81, hydraulic_min_crossings_preferred: int = 2, compute_hydraulic_depth: bool = False) -> lsdtt3.ValleyCrossSections

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
    extract_valley_cross_sections(dem: lsdtt3.Raster, centrelines: lsdtt3.ValleyCentrelineResult, spacing_metres: float = 500.0, half_length_multiplier: float = 1.5, min_half_length_metres: float = 0.0, max_half_length_metres: float = 5000.0, tangent_window_nodes: int = 8, sample_step_metres: float = 0.0, curved_sections: bool = False, curvature: object | None = None, hillshade: object | None = None, channel_mask: object | None = None, side_gap_tolerance_steps: int = 1, feature_peak_min_prominence_fraction: float = 0.15, feature_peak_edge_margin_samples: int = 2, clip_overlapping_sections: bool = True, neighbour_clip_window: int = 2, min_section_reach_fraction: float = 0.35, curved_wall_clearance_fraction: float = 0.1, curved_streamwise_smoothing_passes: int = 40, curved_coarse_spacing_metres: float = 700.0, hydraulic_stage_step_metres: float = 0.1, hydraulic_min_step_height_metres: float = 2.0, hydraulic_min_step_width_pixels: int = 3, hydraulic_midpoint_search_samples: int = 81, hydraulic_min_crossings_preferred: int = 2, compute_hydraulic_depth: bool = False) -> lsdtt3.ValleyCrossSections

Extract straight or boundary-fitted valley cross sections.

Optional curvature, hillshade, and channel-mask rasters enable feature picks.
Hydraulic break picks are always available; set ``compute_hydraulic_depth`` to
retain the full per-stage hydraulic-geometry curve.

    Authors
    -------
    B. Gailleton


```python
def extract_valley_cross_sections(dem: lsdtt3.Raster, centrelines: lsdtt3.ValleyCentrelineResult, spacing_metres: float = 500.0, half_length_multiplier: float = 1.5, min_half_length_metres: float = 0.0, max_half_length_metres: float = 5000.0, tangent_window_nodes: int = 8, sample_step_metres: float = 0.0, curved_sections: bool = False, curvature: object | None = None, hillshade: object | None = None, channel_mask: object | None = None, side_gap_tolerance_steps: int = 1, feature_peak_min_prominence_fraction: float = 0.15, feature_peak_edge_margin_samples: int = 2, clip_overlapping_sections: bool = True, neighbour_clip_window: int = 2, min_section_reach_fraction: float = 0.35, curved_wall_clearance_fraction: float = 0.1, curved_streamwise_smoothing_passes: int = 40, curved_coarse_spacing_metres: float = 700.0, hydraulic_stage_step_metres: float = 0.1, hydraulic_min_step_height_metres: float = 2.0, hydraulic_min_step_width_pixels: int = 3, hydraulic_midpoint_search_samples: int = 81, hydraulic_min_crossings_preferred: int = 2, compute_hydraulic_depth: bool = False) -> lsdtt3.ValleyCrossSections
```

---

## `find_marsh_platform`

find_marsh_platform(dem: lsdtt3.Raster, slope: object | None = None, slope_method: str = 'gradient_3x3', slope_polyfit_lambda_metres: float = 5.0, histogram_n_bins: int = 100, search_space_histogram_derivative_threshold: float = -2.0, max_scarp_growth_iterations: int = 50, max_platform_growth_iterations: int = 100, min_platform_component_pixels: int = 1, preserve_scarps_in_final_platform: bool = True, apply_scarp_closing: bool = False, scarp_closing_radius_metres: float = 0.0, n_threads: int = 0) -> lsdtt3.MarshPlatformResult

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
    find_marsh_platform(dem: lsdtt3.Raster, slope: object | None = None, slope_method: str = 'gradient_3x3', slope_polyfit_lambda_metres: float = 5.0, histogram_n_bins: int = 100, search_space_histogram_derivative_threshold: float = -2.0, max_scarp_growth_iterations: int = 50, max_platform_growth_iterations: int = 100, min_platform_component_pixels: int = 1, preserve_scarps_in_final_platform: bool = True, apply_scarp_closing: bool = False, scarp_closing_radius_metres: float = 0.0, n_threads: int = 0) -> lsdtt3.MarshPlatformResult

Extract marsh platforms and scarps from a projected DEM.

    Authors
    -------
    B. Gailleton


```python
def find_marsh_platform(dem: lsdtt3.Raster, slope: object | None = None, slope_method: str = 'gradient_3x3', slope_polyfit_lambda_metres: float = 5.0, histogram_n_bins: int = 100, search_space_histogram_derivative_threshold: float = -2.0, max_scarp_growth_iterations: int = 50, max_platform_growth_iterations: int = 100, min_platform_component_pixels: int = 1, preserve_scarps_in_final_platform: bool = True, apply_scarp_closing: bool = False, scarp_closing_radius_metres: float = 0.0, n_threads: int = 0) -> lsdtt3.MarshPlatformResult
```

---

## `find_tidal_channels`

find_tidal_channels(dem: lsdtt3.Raster, use_flood_mask: bool = False, flood_elevation_metres: float = 0.0, median_window_radii_metres: collections.abc.Sequence[float] = [5.0, 10.0, 20.0], residual_thresholds_metres: collections.abc.Sequence[float] = [0.1], scale_mode: str = 'any', apply_binary_closing: bool = False, binary_closing_radius_metres: float = 0.0, fill_holes: bool = False, min_channel_component_pixels: int = 1, extract_skeleton: bool = True, skeleton_spur_max_pixels: int = 10, skeleton_terminal_loop_max_pixels: int = 10, skeleton_min_component_pixels: int = 1, connect_skeleton_segments: bool = False, connect_max_gap_pixels: float = 20.0, connect_max_facing_angle_degrees: float = 35.0, connect_max_heading_misalignment_degrees: float = 45.0, connect_max_candidate_score: float = 1.0, connect_route_corridor_half_width_pixels: float = 12.0, compute_upl: bool = True, n_threads: int = 0) -> lsdtt3.TidalChannelResult

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
    find_tidal_channels(dem: lsdtt3.Raster, use_flood_mask: bool = False, flood_elevation_metres: float = 0.0, median_window_radii_metres: collections.abc.Sequence[float] = [5.0, 10.0, 20.0], residual_thresholds_metres: collections.abc.Sequence[float] = [0.1], scale_mode: str = 'any', apply_binary_closing: bool = False, binary_closing_radius_metres: float = 0.0, fill_holes: bool = False, min_channel_component_pixels: int = 1, extract_skeleton: bool = True, skeleton_spur_max_pixels: int = 10, skeleton_terminal_loop_max_pixels: int = 10, skeleton_min_component_pixels: int = 1, connect_skeleton_segments: bool = False, connect_max_gap_pixels: float = 20.0, connect_max_facing_angle_degrees: float = 35.0, connect_max_heading_misalignment_degrees: float = 45.0, connect_max_candidate_score: float = 1.0, connect_route_corridor_half_width_pixels: float = 12.0, compute_upl: bool = True, n_threads: int = 0) -> lsdtt3.TidalChannelResult

Extract tidal channels as multiscale local topographic lows or a flood mask.

    Authors
    -------
    B. Gailleton


```python
def find_tidal_channels(dem: lsdtt3.Raster, use_flood_mask: bool = False, flood_elevation_metres: float = 0.0, median_window_radii_metres: collections.abc.Sequence[float] = [5.0, 10.0, 20.0], residual_thresholds_metres: collections.abc.Sequence[float] = [0.1], scale_mode: str = 'any', apply_binary_closing: bool = False, binary_closing_radius_metres: float = 0.0, fill_holes: bool = False, min_channel_component_pixels: int = 1, extract_skeleton: bool = True, skeleton_spur_max_pixels: int = 10, skeleton_terminal_loop_max_pixels: int = 10, skeleton_min_component_pixels: int = 1, connect_skeleton_segments: bool = False, connect_max_gap_pixels: float = 20.0, connect_max_facing_angle_degrees: float = 35.0, connect_max_heading_misalignment_degrees: float = 45.0, connect_max_candidate_score: float = 1.0, connect_route_corridor_half_width_pixels: float = 12.0, compute_upl: bool = True, n_threads: int = 0) -> lsdtt3.TidalChannelResult
```

---

## `find_valleys`

find_valleys(raster: lsdtt3.Raster, smoothing_window_metres: float = 50.0, curvature_threshold: float = 0.0, include_basin: bool = True, include_synformal_saddle: bool = True, include_synform: bool = True, n_threads: int = 0) -> lsdtt3.ValleyResult

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
    find_valleys(raster: lsdtt3.Raster, smoothing_window_metres: float = 50.0, curvature_threshold: float = 0.0, include_basin: bool = True, include_synformal_saddle: bool = True, include_synform: bool = True, n_threads: int = 0) -> lsdtt3.ValleyResult

Find broad concave valley forms from a smoothed surface-shape classification.

    Authors
    -------
    B. Gailleton


```python
def find_valleys(raster: lsdtt3.Raster, smoothing_window_metres: float = 50.0, curvature_threshold: float = 0.0, include_basin: bool = True, include_synformal_saddle: bool = True, include_synform: bool = True, n_threads: int = 0) -> lsdtt3.ValleyResult
```

---

## `fit_valley_baseline`

fit_valley_baseline(along_distance: object, pin_elevation: object) -> list[float]

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
    fit_valley_baseline(along_distance: object, pin_elevation: object) -> list[float]

Interpolate a non-increasing valley baseline through finite elevation pins.

Use ``numpy.nan`` for unpinned centreline nodes.

    Authors
    -------
    B. Gailleton


```python
def fit_valley_baseline(along_distance: object, pin_elevation: object) -> list[float]
```

---

## `repair_water_lines`

repair_water_lines(water_lines: object, roads: object, snap_tolerance_metres: float = 0.5, max_gap_length_metres: float = 25.0, max_nonroad_gap_metres: float = 30.0, road_intersection_tolerance_metres: float = 1.0, require_mutual_best: bool = True) -> lsdtt3.WaterLineRepairResult

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
    repair_water_lines(water_lines: object, roads: object, snap_tolerance_metres: float = 0.5, max_gap_length_metres: float = 25.0, max_nonroad_gap_metres: float = 30.0, road_intersection_tolerance_metres: float = 1.0, require_mutual_best: bool = True) -> lsdtt3.WaterLineRepairResult

Repair short broken water-line gaps, preferentially across roads.

Parameters
----------
water_lines, roads : geopandas.GeoDataFrame
    Projected line layers in the same CRS.
max_gap_length_metres : float, default=25
    Maximum road-associated bridge length.

Returns
-------
WaterLineRepairResult
    Original linework plus accepted bridge segments and diagnostics.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


```python
def repair_water_lines(water_lines: object, roads: object, snap_tolerance_metres: float = 0.5, max_gap_length_metres: float = 25.0, max_nonroad_gap_metres: float = 30.0, road_intersection_tolerance_metres: float = 1.0, require_mutual_best: bool = True) -> lsdtt3.WaterLineRepairResult
```

---

## `valley_mask_from_cross_sections`

valley_mask_from_cross_sections(dem: lsdtt3.Raster, centrelines: lsdtt3.ValleyCentrelineResult, sections: lsdtt3.ValleyCrossSections) -> lsdtt3.Raster

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
    valley_mask_from_cross_sections(dem: lsdtt3.Raster, centrelines: lsdtt3.ValleyCentrelineResult, sections: lsdtt3.ValleyCrossSections) -> lsdtt3.Raster

Build a valley-floor mask from the median lowest hydraulic break per channel.

This complements the corridor/HAND mask returned by ``extract_valley_centrelines``.

    Authors
    -------
    B. Gailleton


```python
def valley_mask_from_cross_sections(dem: lsdtt3.Raster, centrelines: lsdtt3.ValleyCentrelineResult, sections: lsdtt3.ValleyCrossSections) -> lsdtt3.Raster
```

---

## `water_centrelines`

water_centrelines(areas: object, tributaries: object, epsg: int = 0, cell_size_metres: float = 1.0, padding_metres: float = 5.0, tributary_snap_tolerance_metres: float = 2.0, centreline_snap_max_metres: float = 100.0, simplify_tolerance_metres: float = 1.0, min_component_length_metres: float = 10.0) -> lsdtt3.WaterCentrelineResult

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
    water_centrelines(areas: object, tributaries: object, epsg: int = 0, cell_size_metres: float = 1.0, padding_metres: float = 5.0, tributary_snap_tolerance_metres: float = 2.0, centreline_snap_max_metres: float = 100.0, simplify_tolerance_metres: float = 1.0, min_component_length_metres: float = 10.0) -> lsdtt3.WaterCentrelineResult

Extract central flow paths from water polygons and connect tributary lines.

Parameters
----------
areas, tributaries : geopandas.GeoDataFrame
    Projected water polygons and linework in the same CRS.
epsg : int, default=0
    Output CRS identifier; zero infers it from areas.
cell_size_metres : float, default=1
    Internal skeletonisation-grid cell size.

Returns
-------
WaterCentrelineResult
    Connected linework and count diagnostics.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


```python
def water_centrelines(areas: object, tributaries: object, epsg: int = 0, cell_size_metres: float = 1.0, padding_metres: float = 5.0, tributary_snap_tolerance_metres: float = 2.0, centreline_snap_max_metres: float = 100.0, simplify_tolerance_metres: float = 1.0, min_component_length_metres: float = 10.0) -> lsdtt3.WaterCentrelineResult
```
