# Sediment, Landslides, And Landforms API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `LandslideResult`

In-memory stochastic landslide-event outputs.

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
    In-memory stochastic landslide-event outputs.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


### Members

#### `achieved_volume_m3`

(self) -> float

#### `elevation`

(self) -> lsdtt3.Raster

#### `excavation_depth`

(self) -> lsdtt3.Raster

#### `failure_time`

(self) -> lsdtt3.Raster

#### `n_failed`

(self) -> int

#### `n_landslides`

(self) -> int

#### `requested_volume_m3`

(self) -> float

#### `summary`

(self) -> str

---

## `LandslideSeriesResult`

Outputs from a stochastic sequence of landslide events.

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
    Outputs from a stochastic sequence of landslide events.

The elevation raster accumulates every excavation. ``failure_time`` and
``excavation_depth`` describe the most recent scar at each pixel, which is the
state required by the cosmogenic response functions.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


### Members

#### `active_area_m2`

(self) -> float

#### `duration_yr`

(self) -> float

#### `elevation`

(self) -> lsdtt3.Raster

#### `event_times_yr`

(self) -> list

#### `excavation_depth`

(self) -> lsdtt3.Raster

#### `failure_time`

(self) -> lsdtt3.Raster

#### `fraction_area_failed`

(self) -> float

#### `mean_landslides_per_event`

(self) -> float

#### `n_events`

(self) -> int

#### `n_events_degenerate`

(self) -> int

#### `n_events_short`

(self) -> int

#### `n_failed`

(self) -> int

#### `n_landslides`

(self) -> int

#### `realised_erosion_rate_mm_kyr`

(self) -> float

#### `realised_pixels`

(self) -> list

#### `requested_erosion_rate_mm_kyr`

(self) -> float

#### `requested_volume_m3`

(self) -> float

#### `summary`

(self) -> str

#### `total_volume_m3`

(self) -> float

---

## `MountainFrontResult`

Results from multi-scale mountain-front identification.

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
Results from multi-scale mountain-front identification.

Authors
-------
B. Gailleton

### Members

#### `fitted_surface`

(self) -> lsdtt3.Raster

#### `low_ground_mask`

(self) -> lsdtt3.Raster

#### `mountain_mask`

(self) -> lsdtt3.Raster

#### `relief_rasters`

(self) -> list[lsdtt3.Raster]

#### `residual`

(self) -> lsdtt3.Raster

---

## `Ridgeline`

An ordered, four-connected raster ridgeline with profile-analysis methods.

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
An ordered, four-connected raster ridgeline with profile-analysis methods.

Authors
-------
B. Gailleton

```python
def __init__(self, grid: lsdtt3.Raster, pixels: object) -> None
```

### Members

#### `__init__`

__init__(self, grid: lsdtt3.Raster, pixels: object) -> None

```python
def __init__(self, grid: lsdtt3.Raster, pixels: object) -> None
```

#### `detrend`

detrend(self, smoothing_window_pixels: int) -> None

```python
def detrend(self, smoothing_window_pixels: int) -> None
```

#### `distance_metres`

(self) -> list[float]

#### `elevation`

(self) -> list[float]

#### `fill_elevation`

fill_elevation(self, dem: lsdtt3.Raster) -> None

```python
def fill_elevation(self, dem: lsdtt3.Raster) -> None
```

#### `fit_profile`

fit_profile(self, use_polynomial: bool = False, use_constant: bool = False) -> None

```python
def fit_profile(self, use_polynomial: bool = False, use_constant: bool = False) -> None
```

#### `is_four_connected`

is_four_connected(self, periodic_ns: bool = False, periodic_ew: bool = False) -> bool

```python
def is_four_connected(self, periodic_ns: bool = False, periodic_ew: bool = False) -> bool
```

#### `length_metres`

(self) -> float

#### `pixels`

(self) -> object

#### `positions`

(self) -> object

#### `r_squared_fit`

(self) -> float

#### `reverse`

reverse(self) -> None

```python
def reverse(self) -> None
```

#### `sample`

sample(self, raster: lsdtt3.Raster) -> list[float]

```python
def sample(self, raster: lsdtt3.Raster) -> list[float]
```

#### `segment_snout`

segment_snout(self, smoothing_window_fraction: float) -> None

```python
def segment_snout(self, smoothing_window_fraction: float) -> None
```

#### `slope_fit`

(self) -> float

#### `snout_index`

(self) -> int

#### `steepness`

steepness(self, basin_drainage_area: float, concavity: float) -> float

```python
def steepness(self, basin_drainage_area: float, concavity: float) -> float
```

---

## `RoadFlowResult`

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

#### `channel_mask`

(self) -> lsdtt3.Raster

#### `contributing_pixels`

(self) -> lsdtt3.Raster

#### `flow_path_bearing`

(self) -> lsdtt3.Raster

#### `flow_path_divergence`

(self) -> lsdtt3.Raster

#### `flow_vs_slope_angle`

(self) -> lsdtt3.Raster

#### `regional_slope`

(self) -> lsdtt3.Raster

#### `regional_slope_bearing`

(self) -> lsdtt3.Raster

#### `roughness`

(self) -> lsdtt3.Raster

#### `valley_axis_buffer_mask`

(self) -> lsdtt3.Raster

---

## `RoadResult`

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

#### `contour_alignment`

(self) -> lsdtt3.Raster

#### `cross_curvature`

(self) -> lsdtt3.Raster

#### `feature_orientation_degrees`

(self) -> lsdtt3.Raster

#### `hessian_elongation`

(self) -> lsdtt3.Raster

#### `mask`

(self) -> lsdtt3.Raster

#### `score`

(self) -> lsdtt3.Raster

#### `slope`

(self) -> lsdtt3.Raster

---

## `RoadSmoothnessResult`

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

#### `path_opened`

(self) -> lsdtt3.Raster

#### `road_mask`

(self) -> lsdtt3.Raster

#### `roughness`

(self) -> lsdtt3.Raster

#### `smoothness`

(self) -> lsdtt3.Raster

---

## `SedimentRoutingResult`

Reach-scale bedload routing outputs and channel-raster diagnostics.

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
    Reach-scale bedload routing outputs and channel-raster diagnostics.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


### Members

#### `active_layer_d50`

(self) -> lsdtt3.Raster

#### `deposition`

(self) -> lsdtt3.Raster

#### `n_reaches`

(self) -> int

#### `n_timesteps`

(self) -> int

#### `reaches`

(self) -> object

#### `total_flux_out_m3`

(self) -> float

---

## `ThinLinesResult`

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

#### `binary`

(self) -> lsdtt3.Raster

#### `blob_free`

(self) -> lsdtt3.Raster

#### `mask`

(self) -> lsdtt3.Raster

---

## `basin_ridgeline_extrema`

basin_ridgeline_extrema(basins: lsdtt3.BasinCollection, network: lsdtt3.ChannelNetwork, dem: lsdtt3.Raster, split_method: str = 'longest_channel', snout_smoothing_fraction: float = 0.10000000149011612, detrend_smoothing_pixels: int = 5) -> object

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
    basin_ridgeline_extrema(basins: lsdtt3.BasinCollection, network: lsdtt3.ChannelNetwork, dem: lsdtt3.Raster, split_method: str = 'longest_channel', snout_smoothing_fraction: float = 0.10000000149011612, detrend_smoothing_pixels: int = 5) -> object

Locate topology-aware maxima and minima along selected basin-outline ridgelines.

Parameters
----------
basins : BasinCollection
    Basin outlines to split into ridgelines.
network : ChannelNetwork
    Same-flow topology used to identify triple/quadruple divide junctions.
dem : Raster
    Elevation raster on the shared flow grid.
split_method : {"longest_channel", "maximum_elevation", "equal_length"}, optional
    Rule used to divide each basin outline into ridgelines.

Returns
-------
geopandas.GeoDataFrame
    Extrema points with ridge fitting and topological-distance diagnostics.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


```python
def basin_ridgeline_extrema(basins: lsdtt3.BasinCollection, network: lsdtt3.ChannelNetwork, dem: lsdtt3.Raster, split_method: str = 'longest_channel', snout_smoothing_fraction: float = 0.10000000149011612, detrend_smoothing_pixels: int = 5) -> object
```

---

## `extract_roads_smooth`

extract_roads_smooth(dem: lsdtt3.Raster, roughness_radius_metres: float = 1.5, smoothness_threshold: float = 0.5, path_min_length_pixels: int = 40, min_component_pixels: int = 40, n_threads: int = 0) -> lsdtt3.RoadSmoothnessResult

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
    extract_roads_smooth(dem: lsdtt3.Raster, roughness_radius_metres: float = 1.5, smoothness_threshold: float = 0.5, path_min_length_pixels: int = 40, min_component_pixels: int = 40, n_threads: int = 0) -> lsdtt3.RoadSmoothnessResult

Extract long smooth road-like ribbons from a DEM.

    Authors
    -------
    B. Gailleton


```python
def extract_roads_smooth(dem: lsdtt3.Raster, roughness_radius_metres: float = 1.5, smoothness_threshold: float = 0.5, path_min_length_pixels: int = 40, min_component_pixels: int = 40, n_threads: int = 0) -> lsdtt3.RoadSmoothnessResult
```

---

## `extract_thin_lines`

extract_thin_lines(response: lsdtt3.Raster, threshold: float = 20.0, path_length_pixels: int = 40, min_component_pixels: int = 40, min_elongation: float = 1.0, n_threads: int = 0) -> lsdtt3.ThinLinesResult

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
    extract_thin_lines(response: lsdtt3.Raster, threshold: float = 20.0, path_length_pixels: int = 40, min_component_pixels: int = 40, min_elongation: float = 1.0, n_threads: int = 0) -> lsdtt3.ThinLinesResult

Keep thin elongated lines from a response raster.

    Authors
    -------
    B. Gailleton


```python
def extract_thin_lines(response: lsdtt3.Raster, threshold: float = 20.0, path_length_pixels: int = 40, min_component_pixels: int = 40, min_elongation: float = 1.0, n_threads: int = 0) -> lsdtt3.ThinLinesResult
```

---

## `find_road_flow_discrimination`

find_road_flow_discrimination(dem: lsdtt3.Raster, regional_smoothing_metres: float = 100.0, threshold_contributing_pixels: int = 1000, flow_path_half_window: int = 5, intermediate_window_metres: float = 15.0, roughness_radius_metres: float = 2.5, raster_is_filled: bool = False, boundary_conditions: str = 'oooo', n_threads: int = 0) -> lsdtt3.RoadFlowResult

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
    find_road_flow_discrimination(dem: lsdtt3.Raster, regional_smoothing_metres: float = 100.0, threshold_contributing_pixels: int = 1000, flow_path_half_window: int = 5, intermediate_window_metres: float = 15.0, roughness_radius_metres: float = 2.5, raster_is_filled: bool = False, boundary_conditions: str = 'oooo', n_threads: int = 0) -> lsdtt3.RoadFlowResult

Build flow-direction and terrain diagnostics for road/channel discrimination.

    Authors
    -------
    B. Gailleton


```python
def find_road_flow_discrimination(dem: lsdtt3.Raster, regional_smoothing_metres: float = 100.0, threshold_contributing_pixels: int = 1000, flow_path_half_window: int = 5, intermediate_window_metres: float = 15.0, roughness_radius_metres: float = 2.5, raster_is_filled: bool = False, boundary_conditions: str = 'oooo', n_threads: int = 0) -> lsdtt3.RoadFlowResult
```

---

## `find_roads`

find_roads(dem: lsdtt3.Raster, smoothing_window_metres: float = 8.0, score_threshold: float = 0.5, use_bench_test: bool = True, bench_half_width_metres: float = 6.0, n_threads: int = 0) -> lsdtt3.RoadResult

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
    find_roads(dem: lsdtt3.Raster, smoothing_window_metres: float = 8.0, score_threshold: float = 0.5, use_bench_test: bool = True, bench_half_width_metres: float = 6.0, n_threads: int = 0) -> lsdtt3.RoadResult

Detect narrow, along-contour road benches from a DEM.

    Authors
    -------
    B. Gailleton


```python
def find_roads(dem: lsdtt3.Raster, smoothing_window_metres: float = 8.0, score_threshold: float = 0.5, use_bench_test: bool = True, bench_half_width_metres: float = 6.0, n_threads: int = 0) -> lsdtt3.RoadResult
```

---

## `identify_mountain_front`

identify_mountain_front(dem: lsdtt3.Raster, relief_scales: collections.abc.Sequence[collections.abc.Sequence[float]] = [], require_all_scales_pass: bool = False, relief_window_shape: str = 'square', use_elevation_filter: bool = False, elevation_upper_percentile: float = 0.5, fit_method: str = 'plane_ransac', ransac_inlier_tolerance: float = 20.0, ransac_max_iterations: int = 1000, ransac_seed: int = 1, spectral_lambda_metres: float = 5000.0, residual_percentile: float = 0.8, min_component_pixels: int = 0, infill_valley_radius_metres: float = 0.0, return_diagnostics: bool = False, n_threads: int = 0) -> lsdtt3.MountainFrontResult

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
    identify_mountain_front(dem: lsdtt3.Raster, relief_scales: collections.abc.Sequence[collections.abc.Sequence[float]] = [], require_all_scales_pass: bool = False, relief_window_shape: str = 'square', use_elevation_filter: bool = False, elevation_upper_percentile: float = 0.5, fit_method: str = 'plane_ransac', ransac_inlier_tolerance: float = 20.0, ransac_max_iterations: int = 1000, ransac_seed: int = 1, spectral_lambda_metres: float = 5000.0, residual_percentile: float = 0.8, min_component_pixels: int = 0, infill_valley_radius_metres: float = 0.0, return_diagnostics: bool = False, n_threads: int = 0) -> lsdtt3.MountainFrontResult

Identify elevated terrain above a fitted low-ground piedmont surface.

Parameters
----------
dem : Raster
relief_scales : sequence of [radius_metres, max_relief_metres], optional
    Multi-scale low-ground criteria. The native defaults are used when empty.
fit_method : {"plane_ols", "plane_ransac", "poly2_ols", "spectral_lowpass"}
    Method used to fit the low-ground surface.
residual_percentile : float, default=0.80
    Residual percentile above which pixels are classified as mountain.

Returns
-------
MountainFrontResult
    Low-ground mask, fitted surface, residual, mountain mask, and optional
    local-relief diagnostics.

    Authors
    -------
    B. Gailleton


```python
def identify_mountain_front(dem: lsdtt3.Raster, relief_scales: collections.abc.Sequence[collections.abc.Sequence[float]] = [], require_all_scales_pass: bool = False, relief_window_shape: str = 'square', use_elevation_filter: bool = False, elevation_upper_percentile: float = 0.5, fit_method: str = 'plane_ransac', ransac_inlier_tolerance: float = 20.0, ransac_max_iterations: int = 1000, ransac_seed: int = 1, spectral_lambda_metres: float = 5000.0, residual_percentile: float = 0.8, min_component_pixels: int = 0, infill_valley_radius_metres: float = 0.0, return_diagnostics: bool = False, n_threads: int = 0) -> lsdtt3.MountainFrontResult
```

---

## `landslide_crn_concentration`

landslide_crn_concentration(excavation_depth: lsdtt3.Raster, failure_time: lsdtt3.Raster, background_erosion_rate_g_cm2_yr: float = 0.001, production_scaling: float = 1.0, total_shielding: float = 1.0, nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', report_time_yr: float = 0.0) -> lsdtt3.Raster

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
    landslide_crn_concentration(excavation_depth: lsdtt3.Raster, failure_time: lsdtt3.Raster, background_erosion_rate_g_cm2_yr: float = 0.001, production_scaling: float = 1.0, total_shielding: float = 1.0, nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', report_time_yr: float = 0.0) -> lsdtt3.Raster

Compute one-nuclide surface concentration after landslide excavation.

Parameters
----------
excavation_depth, failure_time : Raster
    Same-grid event rasters from ``trigger_landslides``.
background_erosion_rate_g_cm2_yr, production_scaling, total_shielding : float
    Background CRN conditions.
nuclide, muon_scheme, report_time_yr
    Nuclide, production scheme, and absolute reporting time.

Returns
-------
Raster
    Surface concentration in atoms g-1.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


```python
def landslide_crn_concentration(excavation_depth: lsdtt3.Raster, failure_time: lsdtt3.Raster, background_erosion_rate_g_cm2_yr: float = 0.001, production_scaling: float = 1.0, total_shielding: float = 1.0, nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', report_time_yr: float = 0.0) -> lsdtt3.Raster
```

---

## `landslide_crn_concentrations`

landslide_crn_concentrations(excavation_depth: lsdtt3.Raster, failure_time: lsdtt3.Raster, background_erosion_rate_g_cm2_yr: float = 0.001, production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers', report_time_yr: float = 0.0) -> dict

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
    landslide_crn_concentrations(excavation_depth: lsdtt3.Raster, failure_time: lsdtt3.Raster, background_erosion_rate_g_cm2_yr: float = 0.001, production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers', report_time_yr: float = 0.0) -> dict

Compute Be10, Al26, and C14 concentrations after landslide excavation.

Parameters
----------
excavation_depth, failure_time : Raster
    Same-grid event rasters.
background_erosion_rate_g_cm2_yr, production_scaling, total_shielding, muon_scheme, report_time_yr
    Background CRN conditions.

Returns
-------
dict
    Concentration rasters keyed by nuclide.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


```python
def landslide_crn_concentrations(excavation_depth: lsdtt3.Raster, failure_time: lsdtt3.Raster, background_erosion_rate_g_cm2_yr: float = 0.001, production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers', report_time_yr: float = 0.0) -> dict
```

---

## `mountain_front_outlets`

mountain_front_outlets(mountain_mask: lsdtt3.Raster, contributing_pixels: lsdtt3.Raster, min_contributing_pixels: int) -> lsdtt3.Raster

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
    mountain_front_outlets(mountain_mask: lsdtt3.Raster, contributing_pixels: lsdtt3.Raster, min_contributing_pixels: int) -> lsdtt3.Raster

Return high-contributing-area pixels on a mountain-mask boundary.

``min_contributing_pixels`` is the minimum upstream pixel count required for
an edge pixel. This established Python keyword is independent of the CLI's
``basin_minimum_contributing_pixels`` parameter name.

    Authors
    -------
    B. Gailleton


```python
def mountain_front_outlets(mountain_mask: lsdtt3.Raster, contributing_pixels: lsdtt3.Raster, min_contributing_pixels: int) -> lsdtt3.Raster
```

---

## `route_sediment`

route_sediment(network: lsdtt3.ChannelNetwork, dem: lsdtt3.Raster, formula: str = 'wilcock_crowe', n_timesteps: int = 365, timestep_seconds: float = 86400.0, porosity: float = 0.4, active_layer_depth_metres: float = -1.0, width_coefficient: float = 2.76, width_exponent: float = 0.4, discharge_coefficient: float = 1e-05, discharge_exponent: float = 0.8, d50_metres: float = 0.016, d84_metres: float = 0.045, d90_metres: float = 0.055, mannings_n: float = 0.04, grain_sizes_mm: object | None = None, grain_fractions: object | None = None, discharge: object | None = None, n_threads: int = 0) -> lsdtt3.SedimentRoutingResult

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
    route_sediment(network: lsdtt3.ChannelNetwork, dem: lsdtt3.Raster, formula: str = 'wilcock_crowe', n_timesteps: int = 365, timestep_seconds: float = 86400.0, porosity: float = 0.4, active_layer_depth_metres: float = -1.0, width_coefficient: float = 2.76, width_exponent: float = 0.4, discharge_coefficient: float = 1e-05, discharge_exponent: float = 0.8, d50_metres: float = 0.016, d84_metres: float = 0.045, d90_metres: float = 0.055, mannings_n: float = 0.04, grain_sizes_mm: object | None = None, grain_fractions: object | None = None, discharge: object | None = None, n_threads: int = 0) -> lsdtt3.SedimentRoutingResult

Route bedload through an existing D8 channel network.

Each channel reach begins with an active sediment layer. Incoming material and
locally mobilised material are limited by the selected transport capacity,
producing reach-scale deposition, mobilisation, grain-size, and outlet-flux
diagnostics. This operation does not modify the DEM.

Parameters
----------
network : ChannelNetwork
    D8 channel network constructed from the DEM's FlowInfo topology.
dem : Raster
    Same-grid elevation raster used to derive reach geometry.
formula : {"wong_parker", "wilcock_crowe", "engelund_hansen"}, optional
    Bedload transport-capacity relation.
n_timesteps, timestep_seconds : int, float, optional
    Number and duration of routing timesteps.
 porosity, active_layer_depth_metres : float, optional
    Sediment porosity and initial active-layer depth; a non-positive depth uses
    twice D90.
width_coefficient, width_exponent : float, optional
    Width relation ``W = coefficient * drainage_area ** exponent``.
discharge_coefficient, discharge_exponent : float, optional
    Default discharge relation when ``discharge`` is omitted.
 d50_metres, d84_metres, d90_metres, mannings_n : float, optional
    Uniform reach sediment and hydraulic properties.
grain_sizes_mm, grain_fractions : array-like, optional
    Equal-length grain-size classes and mass fractions summing to one.
discharge : array-like, optional
    Explicit non-negative discharge array with shape
    ``(n_timesteps, network.n_junctions)`` in m³ s⁻¹.
n_threads : int, optional
    OpenMP thread count; zero uses the runtime default.

Returns
-------
SedimentRoutingResult
    Channel-raster diagnostics, reach table, and total outlet sediment flux.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


**Deprecated compatibility names:**
- `active_layer_depth_m` -> `active_layer_depth_metres`; earliest removal: following minor.
- `d50_m` -> `d50_metres`; earliest removal: following minor.
- `d84_m` -> `d84_metres`; earliest removal: following minor.
- `d90_m` -> `d90_metres`; earliest removal: following minor.

```python
def route_sediment(network: lsdtt3.ChannelNetwork, dem: lsdtt3.Raster, formula: str = 'wilcock_crowe', n_timesteps: int = 365, timestep_seconds: float = 86400.0, porosity: float = 0.4, active_layer_depth_metres: float = -1.0, width_coefficient: float = 2.76, width_exponent: float = 0.4, discharge_coefficient: float = 1e-05, discharge_exponent: float = 0.8, d50_metres: float = 0.016, d84_metres: float = 0.045, d90_metres: float = 0.055, mannings_n: float = 0.04, grain_sizes_mm: object | None = None, grain_fractions: object | None = None, discharge: object | None = None, n_threads: int = 0) -> lsdtt3.SedimentRoutingResult
```

---

## `run_landslide_series`

run_landslide_series(dem: lsdtt3.Raster, flow: lsdtt3.FlowInfo, start_time_yr: float = 0.0, end_time_yr: float = 100000.0, mean_recurrence_yr: float = 1000.0, poisson_timing: bool = True, target_volume_per_event_m3: float = 0.0, long_term_erosion_rate_mm_kyr: float = 0.0, channel_threshold_pixels: int = 1000, size_model: str = 'inverse_gamma', area_volume_preset: str = 'larsen_soil_scar', min_pixels: int = 1, max_area_m2: float = 0.0, anneal_time_yr: float = 0.0, event_budget_relative_tolerance: float = 0.02, budget_max_iter: int = 50, seed: int = 1, active_mask: lsdtt3.Raster | None = None) -> lsdtt3.LandslideSeriesResult

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
    run_landslide_series(dem: lsdtt3.Raster, flow: lsdtt3.FlowInfo, start_time_yr: float = 0.0, end_time_yr: float = 100000.0, mean_recurrence_yr: float = 1000.0, poisson_timing: bool = True, target_volume_per_event_m3: float = 0.0, long_term_erosion_rate_mm_kyr: float = 0.0, channel_threshold_pixels: int = 1000, size_model: str = 'inverse_gamma', area_volume_preset: str = 'larsen_soil_scar', min_pixels: int = 1, max_area_m2: float = 0.0, anneal_time_yr: float = 0.0, event_budget_relative_tolerance: float = 0.02, budget_max_iter: int = 50, seed: int = 1, active_mask: lsdtt3.Raster | None = None) -> lsdtt3.LandslideSeriesResult

Run a reproducible, annealing-aware landslide time series.

Provide exactly one budget: a fixed ``target_volume_per_event_m3`` or a
``long_term_erosion_rate_mm_kyr`` over the eligible host area. The final DEM
accumulates scars; the accompanying failure/depth rasters retain the latest
event at every pixel and can feed the cosmogenic functions directly.

Parameters
----------
dem : Raster
    Depression-conditioned elevation raster.
flow : FlowInfo
    Same-grid D8 routing topology.
start_time_yr, end_time_yr : float
    Inclusive simulation time interval in years.
mean_recurrence_yr : float
    Mean separation of events in years.
poisson_timing : bool
    Draw exponential waiting times when True; otherwise use fixed recurrence.
target_volume_per_event_m3, long_term_erosion_rate_mm_kyr : float
    Provide exactly one event-budget mode.
channel_threshold_pixels, size_model, area_volume_preset :
    Native scar-placement and size-distribution controls.
min_pixels, max_area_m2, anneal_time_yr :
    Scar bounds and time before failed ground becomes eligible again.
 event_budget_relative_tolerance, budget_max_iter, seed, active_mask :
    Event-budget solver controls, reproducibility, and optional host mask.

Returns
-------
LandslideSeriesResult
    Final terrain and scar rasters, timing, volume/rate closure, and event
    diagnostics.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


**Deprecated compatibility names:**
- `budget_tol` -> `event_budget_relative_tolerance`; earliest removal: following minor.

```python
def run_landslide_series(dem: lsdtt3.Raster, flow: lsdtt3.FlowInfo, start_time_yr: float = 0.0, end_time_yr: float = 100000.0, mean_recurrence_yr: float = 1000.0, poisson_timing: bool = True, target_volume_per_event_m3: float = 0.0, long_term_erosion_rate_mm_kyr: float = 0.0, channel_threshold_pixels: int = 1000, size_model: str = 'inverse_gamma', area_volume_preset: str = 'larsen_soil_scar', min_pixels: int = 1, max_area_m2: float = 0.0, anneal_time_yr: float = 0.0, event_budget_relative_tolerance: float = 0.02, budget_max_iter: int = 50, seed: int = 1, active_mask: lsdtt3.Raster | None = None) -> lsdtt3.LandslideSeriesResult
```

---

## `sediment_connectivity_index`

sediment_connectivity_index(flow: lsdtt3.FlowInfo, roughness: lsdtt3.Raster, slope: lsdtt3.Raster, drainage_area: lsdtt3.Raster, area_threshold_m2: float) -> lsdtt3.ConnectivityResult

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
    sediment_connectivity_index(flow: lsdtt3.FlowInfo, roughness: lsdtt3.Raster, slope: lsdtt3.Raster, drainage_area: lsdtt3.Raster, area_threshold_m2: float) -> lsdtt3.ConnectivityResult

Compute sediment connectivity using D8 topology; drainage area may be D8, MFD, or D-infinity.

    Authors
    -------
    B. Gailleton


```python
def sediment_connectivity_index(flow: lsdtt3.FlowInfo, roughness: lsdtt3.Raster, slope: lsdtt3.Raster, drainage_area: lsdtt3.Raster, area_threshold_m2: float) -> lsdtt3.ConnectivityResult
```

---

## `trigger_landslides`

trigger_landslides(dem: lsdtt3.Raster, flow: lsdtt3.FlowInfo, target_volume_m3: float, event_time_yr: float = 0.0, channel_threshold_pixels: int = 1000, size_model: str = 'inverse_gamma', area_volume_preset: str = 'larsen_soil_scar', min_pixels: int = 1, max_area_m2: float = 0.0, anneal_time_yr: float = 0.0, seed: int = 1, active_mask: lsdtt3.Raster | None = None) -> lsdtt3.LandslideResult

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
    trigger_landslides(dem: lsdtt3.Raster, flow: lsdtt3.FlowInfo, target_volume_m3: float, event_time_yr: float = 0.0, channel_threshold_pixels: int = 1000, size_model: str = 'inverse_gamma', area_volume_preset: str = 'larsen_soil_scar', min_pixels: int = 1, max_area_m2: float = 0.0, anneal_time_yr: float = 0.0, seed: int = 1, active_mask: lsdtt3.Raster | None = None) -> lsdtt3.LandslideResult

Trigger a volume-budgeted stochastic landslide event.

Parameters
----------
dem : Raster
    Depression-conditioned elevation raster.
flow : FlowInfo
    Same-grid D8 routing topology.
target_volume_m3 : float
    Requested total scar volume.
event_time_yr, channel_threshold_pixels, size_model, area_volume_preset
    Event timing and native placement/size controls.
min_pixels, max_area_m2, anneal_time_yr, seed, active_mask
    Scar bounds, recovery, reproducibility, and optional host restriction.

Returns
-------
LandslideResult
    Excavated surface, failure-time/depth rasters, and budget diagnostics.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


```python
def trigger_landslides(dem: lsdtt3.Raster, flow: lsdtt3.FlowInfo, target_volume_m3: float, event_time_yr: float = 0.0, channel_threshold_pixels: int = 1000, size_model: str = 'inverse_gamma', area_volume_preset: str = 'larsen_soil_scar', min_pixels: int = 1, max_area_m2: float = 0.0, anneal_time_yr: float = 0.0, seed: int = 1, active_mask: lsdtt3.Raster | None = None) -> lsdtt3.LandslideResult
```
