# Sediment, Landslides, And Landforms API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `LandslideResult`

In-memory stochastic landslide-event outputs.

Returned by ``trigger_landslides()``; not constructible directly.

Attributes
----------
achieved_volume_m3 : float
elevation : Raster
excavation_depth : Raster
failure_time : Raster
n_failed : int
n_landslides : int
requested_volume_m3 : float
summary : str

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

The elevation raster accumulates every excavation. ``failure_time`` and
``excavation_depth`` describe the most recent scar at each pixel, which is the
state required by the cosmogenic response functions.

Returned by ``run_landslide_series()``; not constructible directly.

Attributes
----------
active_area_m2 : float
duration_yr : float
elevation : Raster
event_times_yr : list
excavation_depth : Raster
failure_time : Raster
fraction_area_failed : float
mean_landslides_per_event : float
n_events : int
n_events_degenerate : int
n_events_short : int
n_failed : int
n_landslides : int
realised_erosion_rate_mm_kyr : float
realised_pixels : list
requested_erosion_rate_mm_kyr : float
requested_volume_m3 : float
summary : str
total_volume_m3 : float

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

Returned by ``identify_mountain_front()``; not constructible directly.

Attributes
----------
fitted_surface : Raster
low_ground_mask : Raster
mountain_mask : Raster
relief_rasters : list[Raster]
residual : Raster

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

Construct it as

    Ridgeline(grid: Raster, pixels: object)

Parameters
----------
grid : Raster
    Raster whose grid the pixels refer to.
pixels : array-like of int
    N x 2 ``(row, col)`` pixels along the ridge, in order (at least two).

Attributes
----------
distance_metres : list[float]
elevation : list[float]
length_metres : float
pixels : numpy.ndarray
    N x 2 int64 array of (row, col).
positions : numpy.ndarray
    N x 2 float64 array of (x, y) map coordinates.
r_squared_fit : float
slope_fit : float
snout_index : int

Methods
-------
detrend(smoothing_window_pixels: int) -> None
fill_elevation(dem: Raster) -> None
fit_profile(use_polynomial: bool = False, use_constant: bool = False) -> None
is_four_connected(periodic_ns: bool = False, periodic_ew: bool = False) -> bool
reverse() -> None
sample(raster: Raster) -> list[float]
segment_snout(smoothing_window_fraction: float) -> None
steepness(basin_drainage_area: float, concavity: float) -> float

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

N x 2 int64 array of (row, col).

#### `positions`

N x 2 float64 array of (x, y) map coordinates.

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

Flow-direction diagnostics from :func:`find_road_flow_discrimination`.

Compares the local D8 flow-path bearing with the regional slope bearing:
roads divert flow along the contour, so ``flow_vs_slope_angle`` is large on
roads and small on natural slopes and channels. Bearings are compass
degrees; slopes are m/m.

Returned by ``find_road_flow_discrimination()``; not constructible directly.

Attributes
----------
channel_mask : Raster
contributing_pixels : Raster
flow_path_bearing : Raster
flow_path_divergence : Raster
flow_vs_slope_angle : Raster
regional_slope : Raster
regional_slope_bearing : Raster
roughness : Raster
valley_axis_buffer_mask : Raster

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

Hillside road / bench detection output from :func:`find_roads`.

``score`` is elongation x contour alignment x edge strength, ``mask`` its
gated threshold; ``contour_alignment`` and ``hessian_elongation`` are in
[0, 1], ``cross_curvature`` in 1/m, ``feature_orientation_degrees`` the
along-feature azimuth (0-180, from north) and ``slope`` the gradient (m/m).

Returned by ``find_roads()``; not constructible directly.

Attributes
----------
contour_alignment : Raster
cross_curvature : Raster
feature_orientation_degrees : Raster
hessian_elongation : Raster
mask : Raster
score : Raster
slope : Raster

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

Smooth-ribbon road extraction output from :func:`extract_roads_smooth`.

``roughness`` is the local elevation range (m), ``smoothness`` =
1 / (1 + roughness / h), ``path_opened`` the smoothness after a path opening
that keeps long ribbons, and ``road_mask`` the final 0/1 mask.

Returned by ``extract_roads_smooth()``; not constructible directly.

Attributes
----------
path_opened : Raster
road_mask : Raster
roughness : Raster
smoothness : Raster

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

Returned by ``route_sediment()``; not constructible directly.

Attributes
----------
active_layer_d50 : Raster
deposition : Raster
n_reaches : int
n_timesteps : int
reaches : object
total_flux_out_m3 : float

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

Thin-line extraction output from :func:`extract_thin_lines`.

``binary`` is the thresholded response, ``blob_free`` the result of the path
opening, and ``mask`` the final 0/1 line mask.

Returned by ``extract_thin_lines()``; not constructible directly.

Attributes
----------
binary : Raster
blob_free : Raster
mask : Raster

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
snout_smoothing_fraction : float, optional
    Smoothing window, as a fraction of ridge length, used to find the ridge
    snout from its curvature.
detrend_smoothing_pixels : int, optional
    Smoothing window (pixels) used to detrend each ridge profile.

Returns
-------
geopandas.GeoDataFrame
    Extrema points with ridge fitting and topological-distance diagnostics.

Authors
-------
B. Gailleton

```python
def basin_ridgeline_extrema(basins: lsdtt3.BasinCollection, network: lsdtt3.ChannelNetwork, dem: lsdtt3.Raster, split_method: str = 'longest_channel', snout_smoothing_fraction: float = 0.10000000149011612, detrend_smoothing_pixels: int = 5) -> object
```

---

## `extract_roads_smooth`

Extract long smooth road-like ribbons from a DEM.

Road surfaces are smoother than the surrounding terrain: smoothness is
computed from local roughness, a path opening keeps only long smooth
ribbons, and small components are removed.

Parameters
----------
dem : Raster
    Elevation raster (m), usually a high-resolution lidar DEM.
roughness_radius_metres : float, default=1.5
    Radius (m) of the local elevation-range window.
smoothness_threshold : float, default=0.5
    Path-opened smoothness (0-1) at or above which a pixel is road.
path_min_length_pixels : int, default=40
    Minimum path length (pixels) the path opening keeps.
min_component_pixels : int, default=40
    Smallest connected component kept (pixels).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
RoadSmoothnessResult

Authors
-------
B. Gailleton

```python
def extract_roads_smooth(dem: lsdtt3.Raster, roughness_radius_metres: float = 1.5, smoothness_threshold: float = 0.5, path_min_length_pixels: int = 40, min_component_pixels: int = 40, n_threads: int = 0) -> lsdtt3.RoadSmoothnessResult
```

---

## `extract_thin_lines`

Keep thin elongated lines from a response raster.

Parameters
----------
response : Raster
    Line-response raster, for example a Frangi response or a path-opening
    output.
threshold : float, default=20.0
    Response at or above which a pixel is a line candidate.
path_length_pixels : int, default=40
    Minimum path length (pixels) of the binary path opening; 0 skips it.
min_component_pixels : int, default=40
    Smallest connected component kept (pixels).
min_elongation : float, default=1.0
    Smallest component elongation kept; 1 disables the test.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
ThinLinesResult

Authors
-------
B. Gailleton

```python
def extract_thin_lines(response: lsdtt3.Raster, threshold: float = 20.0, path_length_pixels: int = 40, min_component_pixels: int = 40, min_elongation: float = 1.0, n_threads: int = 0) -> lsdtt3.ThinLinesResult
```

---

## `find_road_flow_discrimination`

Build flow-direction and terrain diagnostics for road/channel discrimination.

Parameters
----------
dem : Raster
    Elevation raster (m).
regional_smoothing_metres : float, default=100.0
    Polyfit length scale (m) of the regional slope field.
threshold_contributing_pixels : int, default=1000
    Channel threshold in contributing pixels.
flow_path_half_window : int, default=5
    Pixels up and down each flow path used to measure its bearing.
intermediate_window_metres : float, default=15.0
    Polyfit length scale (m) of the intermediate-scale surface.
roughness_radius_metres : float, default=2.5
    Radius (m) of the local roughness window.
raster_is_filled : bool, default=False
    Skip the internal pit filling if the DEM is already conditioned.
boundary_conditions : str, default='oooo'
    NESW boundary codes for the internal flow routing.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
RoadFlowResult

Authors
-------
B. Gailleton
S. M. Mudd

```python
def find_road_flow_discrimination(dem: lsdtt3.Raster, regional_smoothing_metres: float = 100.0, threshold_contributing_pixels: int = 1000, flow_path_half_window: int = 5, intermediate_window_metres: float = 15.0, roughness_radius_metres: float = 2.5, raster_is_filled: bool = False, boundary_conditions: str = 'oooo', n_threads: int = 0) -> lsdtt3.RoadFlowResult
```

---

## `find_roads`

Detect narrow, along-contour road benches from a DEM.

Roads cut into hillslopes are elongated, low-curvature benches whose long axis
follows the contour, unlike channels which run downslope. The score
combines Hessian elongation, alignment with the contour and cut/fill edge
strength from a polyfit at ``smoothing_window_metres``.

Parameters
----------
dem : Raster
    Elevation raster (m), usually a high-resolution lidar DEM.
smoothing_window_metres : float, default=8.0
    Polyfit length scale (m).
score_threshold : float, default=0.5
    Score at or above which a pixel enters ``mask``.
use_bench_test : bool, default=True
    Require both flanks (sampled across the feature) to be steeper than the
    tread, rejecting one-sided edges such as valley-floor margins.
bench_half_width_metres : float, default=6.0
    Distance (m) to each side at which the flank slopes are sampled.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
RoadResult

Authors
-------
B. Gailleton

```python
def find_roads(dem: lsdtt3.Raster, smoothing_window_metres: float = 8.0, score_threshold: float = 0.5, use_bench_test: bool = True, bench_half_width_metres: float = 6.0, n_threads: int = 0) -> lsdtt3.RoadResult
```

---

## `identify_mountain_front`

Identify elevated terrain above a fitted low-ground piedmont surface.

Low ground (the piedmont) is found where local relief is small, a surface
is fitted to it, and pixels standing high above that surface are mountain.

Parameters
----------
dem : Raster
    Elevation raster (m).
relief_scales : sequence of [radius_metres, max_relief_metres], optional
    Multi-scale low-ground criteria: a pixel is low ground when its local
    relief within ``radius_metres`` is below ``max_relief_metres``. Empty
    uses ``[[500, 50], [2000, 200]]``.
require_all_scales_pass : bool, default=False
    Require every scale to pass (True) or only the largest (False).
relief_window_shape : {"square", "circular"}, default="square"
    Shape of the relief windows.
use_elevation_filter : bool, default=False
    Also exclude pixels above an elevation percentile from low ground.
elevation_upper_percentile : float, default=0.5
    That percentile, as a fraction (0-1).
fit_method : {"plane_ols", "plane_ransac", "poly2_ols", "spectral_lowpass"}
    Method used to fit the low-ground surface.
ransac_inlier_tolerance : float, default=20.0
    ``plane_ransac``: largest residual (m) of an inlier.
ransac_max_iterations : int, default=1000
    ``plane_ransac``: number of iterations.
ransac_seed : int, default=1
    ``plane_ransac``: random seed.
spectral_lambda_metres : float, default=5000.0
    ``spectral_lowpass``: cut-off wavelength (m); larger than the piedmont
    width.
residual_percentile : float, default=0.80
    Residual percentile above which pixels are classified as mountain
    (0.8 = the top 20% of residuals).
min_component_pixels : int, default=0
    Remove mountain components smaller than this (pixels); 0 skips.
infill_valley_radius_metres : float, default=0.0
    Radius (m) of a closing that seals valley notches in the mountain mask;
    0 skips.
return_diagnostics : bool, default=False
    Keep the per-scale local-relief rasters (``relief_rasters``).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

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

```python
def landslide_crn_concentration(excavation_depth: lsdtt3.Raster, failure_time: lsdtt3.Raster, background_erosion_rate_g_cm2_yr: float = 0.001, production_scaling: float = 1.0, total_shielding: float = 1.0, nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', report_time_yr: float = 0.0) -> lsdtt3.Raster
```

---

## `landslide_crn_concentrations`

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

```python
def landslide_crn_concentrations(excavation_depth: lsdtt3.Raster, failure_time: lsdtt3.Raster, background_erosion_rate_g_cm2_yr: float = 0.001, production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers', report_time_yr: float = 0.0) -> dict
```

---

## `mountain_front_outlets`

Return high-contributing-area pixels on a mountain-mask boundary.

``min_contributing_pixels`` is the minimum upstream pixel count required for
an edge pixel. This established Python keyword is independent of the CLI's
``basin_minimum_contributing_pixels`` parameter name.

Parameters
----------
mountain_mask : Raster
    0/1 mountain mask (``MountainFrontResult.mountain_mask``).
contributing_pixels : Raster
    Contributing-pixel counts on the same grid
    (``FlowInfo.contributing_pixels()``).
min_contributing_pixels : int
    Smallest contributing-pixel count of a returned pixel.

Returns
-------
Raster
    0/1 raster marking the outlet pixels; use
    ``np.argwhere(r.to_numpy() > 0.5)`` for ``(row, col)`` outlets.

Authors
-------
B. Gailleton
S. M. Mudd

```python
def mountain_front_outlets(mountain_mask: lsdtt3.Raster, contributing_pixels: lsdtt3.Raster, min_contributing_pixels: int) -> lsdtt3.Raster
```

---

## `route_sediment`

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
S. M. Mudd

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
channel_threshold_pixels : int
    Channel threshold (contributing pixels) used by scar placement.
size_model : {"inverse_gamma", "truncated_pareto"}
    Landslide area distribution.
area_volume_preset : str
    Area-volume scaling preset: ``"larsen_global"``, ``"larsen_soil_scar"``,
    ``"larsen_bedrock_scar"``, ``"larsen_soil_all"``,
    ``"larsen_bedrock_all"``, ``"hovius_slab"`` or ``"hovius_parabolic"``.
min_pixels : int
    Smallest scar (pixels).
max_area_m2 : float
    Largest scar area (m^2); 0 means no cap.
anneal_time_yr : float
    Years before failed ground can fail again.
event_budget_relative_tolerance : float
    Relative tolerance of each event's volume against its budget.
budget_max_iter : int
    Largest number of budget-solver iterations per event.
seed : int
    Random seed.
active_mask : Raster, optional
    0/1 raster of pixels allowed to fail; None allows every pixel.

Returns
-------
LandslideSeriesResult
    Final terrain and scar rasters, timing, volume/rate closure, and event
    diagnostics.

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `budget_tol` -> `event_budget_relative_tolerance`; earliest removal: following minor.

```python
def run_landslide_series(dem: lsdtt3.Raster, flow: lsdtt3.FlowInfo, start_time_yr: float = 0.0, end_time_yr: float = 100000.0, mean_recurrence_yr: float = 1000.0, poisson_timing: bool = True, target_volume_per_event_m3: float = 0.0, long_term_erosion_rate_mm_kyr: float = 0.0, channel_threshold_pixels: int = 1000, size_model: str = 'inverse_gamma', area_volume_preset: str = 'larsen_soil_scar', min_pixels: int = 1, max_area_m2: float = 0.0, anneal_time_yr: float = 0.0, event_budget_relative_tolerance: float = 0.02, budget_max_iter: int = 50, seed: int = 1, active_mask: lsdtt3.Raster | None = None) -> lsdtt3.LandslideSeriesResult
```

---

## `sediment_connectivity_index`

Compute sediment connectivity using D8 topology; drainage area may be D8, MFD, or D-infinity.

Cavalli et al. (2013) index ``IC = log10(D_up / D_down)`` with
``D_up = mean(W) * mean(S) * sqrt(A)`` over the upslope area and
``D_down`` the sum of ``d / (W S)`` along the D8 path to the nearest sink;
``W = max(1 - R / R_max, 0.001)`` from the roughness R.

Parameters
----------
flow : FlowInfo
    D8 flow topology.
roughness : Raster
    Surface roughness on the flow grid (for example a local standard
    deviation of residual elevation, m).
slope : Raster
    Gradient (m/m); values below 0.001 are clipped in ``D_down``.
drainage_area : Raster
    Drainage area (m^2) from any routing scheme.
area_threshold_m2 : float
    Pixels with at least this drainage area (m^2) are sinks (channels).

Returns
-------
ConnectivityResult
    ``index``, ``weight`` (W), ``upstream`` (D_up) and ``downstream``
    (D_down) Rasters.

Authors
-------
B. Gailleton

```python
def sediment_connectivity_index(flow: lsdtt3.FlowInfo, roughness: lsdtt3.Raster, slope: lsdtt3.Raster, drainage_area: lsdtt3.Raster, area_threshold_m2: float) -> lsdtt3.ConnectivityResult
```

---

## `trigger_landslides`

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

```python
def trigger_landslides(dem: lsdtt3.Raster, flow: lsdtt3.FlowInfo, target_volume_m3: float, event_time_yr: float = 0.0, channel_threshold_pixels: int = 1000, size_model: str = 'inverse_gamma', area_volume_preset: str = 'larsen_soil_scar', min_pixels: int = 1, max_area_m2: float = 0.0, anneal_time_yr: float = 0.0, seed: int = 1, active_mask: lsdtt3.Raster | None = None) -> lsdtt3.LandslideResult
```
