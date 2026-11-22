# Sediment, Landslides, And Landforms API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `GroundwaterDupuit`

An unconfined aquifer on an impermeable base: the Dupuit - Boussinesq model.

lsdtt-slope-stability ``water_table_method: dupuit_steady`` and
``dupuit_transient`` (after Landlab's GroundwaterDupuitPercolator, Litwin et
al. 2020). The base lies ``regolith_thickness_metres`` below the surface;
water flows down the water-table gradient; ``fixed_head`` cells (the CLI's
channels) are seepage outlets with the water table at the surface; the
raster edge and nodata are no-flux. Where the regolith saturates, water
seeps out. The keyword names are the CLI's less the ``groundwater_`` prefix.

Members: ``run_to_steady``, ``run``, ``set_thickness``,
``water_table_depth``, ``saturated_thickness``, ``seepage``, ``storage``,
``budget``, ``steady_rates``.

Construct it as

    GroundwaterDupuit(surface: Raster, fixed_head: Raster, *, regolith_thickness_metres: float = 2.0, regolith_thickness: Raster | None = None, hydraulic_conductivity_m_per_s: float = 0.001, drainable_porosity: float = 0.2, transmissivity_profile: str = 'uniform', conductivity_decay_depth_metres: float = 1.0, n_threads: int = 0)

Parameters
----------
surface : Raster
    Ground elevation (m), normally the input DEM.
fixed_head : Raster
    Nonzero = fixed-head (channel) cell, e.g. ``RupturePlanes.channel``.
regolith_thickness_metres : float, default=2.0, keyword-only
    Aquifer thickness d (m), used where ``regolith_thickness`` is None.
regolith_thickness : Raster or None, default=None, keyword-only
    Thickness per cell (m); nodata or <= 0 cells are left out.
hydraulic_conductivity_m_per_s : float, default=0.001, keyword-only
    K, or K0 at the surface for the exponential profile.
drainable_porosity : float, default=0.2, keyword-only
    n.
transmissivity_profile : {"uniform", "exponential"}, default="uniform", keyword-only
    exponential: K decays with depth over conductivity_decay_depth_metres
    (TOPMODEL's assumption).
conductivity_decay_depth_metres : float, default=1.0, keyword-only
    z_m of the exponential profile.
n_threads : int, default=0, keyword-only
    OpenMP threads; the result does not depend on it.

Attributes
----------
storage : float
    Water stored, n h summed over the active cells (m^3).

Methods
-------
budget() -> dict
    The cumulative water budget of the explicit steady march and the runs (m^3).
run(duration_s: float, recharge_m_per_s: float | Raster, *, transient_method: str = 'explicit', implicit_max_dt_s: float = 3600.0, implicit_dt_tolerance_metres: float = 0.0) -> dict
    Advance the water table by ``duration_s`` under a constant recharge.
run_to_steady(recharge_m_per_s: float | Raster, *, steady_method: str = 'implicit', steady_tolerance: float = 1e-06, max_substeps: int = 10000000, max_newton_iterations: int = 200) -> dict
    Find the steady water table under a constant recharge.
saturated_thickness() -> Raster
    Saturated thickness h (m).
seepage() -> Raster
    Seepage rate (m/s) from the last call.
set_thickness(thickness: Raster) -> None
    Set the saturated thickness h of the active cells (clipped to [0, d]).
steady_rates() -> dict
    The budget of the last implicit steady solve, as rates (m^3/s).
water_table_depth() -> Raster
    Depth of the water table below the surface (m): the input to ``FlowInfo.rupture_planes``.

Authors
-------
B. Gailleton

```python
def __init__(self, surface: Raster, fixed_head: Raster, *, regolith_thickness_metres: float = 2.0, regolith_thickness: Raster | None = None, hydraulic_conductivity_m_per_s: float = 0.001, drainable_porosity: float = 0.2, transmissivity_profile: str = 'uniform', conductivity_decay_depth_metres: float = 1.0, n_threads: int = 0) -> None
```

### Members

#### `__init__`

__init__(self, surface: Raster, fixed_head: Raster, *, regolith_thickness_metres: float = 2.0, regolith_thickness: Raster | None = None, hydraulic_conductivity_m_per_s: float = 0.001, drainable_porosity: float = 0.2, transmissivity_profile: str = 'uniform', conductivity_decay_depth_metres: float = 1.0, n_threads: int = 0) -> None

```python
def __init__(self, surface: Raster, fixed_head: Raster, *, regolith_thickness_metres: float = 2.0, regolith_thickness: Raster | None = None, hydraulic_conductivity_m_per_s: float = 0.001, drainable_porosity: float = 0.2, transmissivity_profile: str = 'uniform', conductivity_decay_depth_metres: float = 1.0, n_threads: int = 0) -> None
```

#### `budget`

budget(self) -> dict

The cumulative water budget of the explicit steady march and the runs (m^3).

recharge_m3 + clipped_added_m3 = storage_change_m3 + outflow_fixed_head_m3 +
seepage_m3 to round-off. An implicit steady solve reports ``steady_rates()``
instead and does not add to this.

Parameters
----------
self : GroundwaterDupuit
    The aquifer.

Returns
-------
dict
    ``recharge_m3, storage_change_m3, outflow_fixed_head_m3, seepage_m3,
    clipped_added_m3``.

```python
def budget(self) -> dict
```

#### `run`

run(self, duration_s: float, recharge_m_per_s: float | Raster, *, transient_method: str = 'explicit', implicit_max_dt_s: float = 3600.0, implicit_dt_tolerance_metres: float = 0.0) -> dict

Advance the water table by ``duration_s`` under a constant recharge.

lsdtt-slope-stability ``dupuit_transient`` runs its recharge series as one
call per constant segment. ``budget()`` accumulates over calls.

Parameters
----------
duration_s : float
    Time to advance (s), > 0.
recharge_m_per_s : float or Raster
    Recharge f (m/s); nodata is 0.
transient_method : {"explicit", "implicit"}, default="explicit", keyword-only
    explicit: substeps at the stability limit (set by stability, not
    accuracy: on coarse grids they can be hours long). implicit: backward
    Euler, each step solved by Newton's method; first order in the step and
    damps fast transients.
implicit_max_dt_s : float, default=3600.0, keyword-only
    implicit: longest step (s).
implicit_dt_tolerance_metres : float, default=0.0, keyword-only
    implicit: step-doubling error tolerance on the water table (m); 0 = off.

Returns
-------
dict
    ``transient_method``, ``steps``; implicit also ``newton_iterations``,
    ``linear_iterations``, ``splits``, ``rejected_by_error``,
    ``max_error_estimate_metres``.

```python
def run(self, duration_s: float, recharge_m_per_s: float | Raster, *, transient_method: str = 'explicit', implicit_max_dt_s: float = 3600.0, implicit_dt_tolerance_metres: float = 0.0) -> dict
```

#### `run_to_steady`

run_to_steady(self, recharge_m_per_s: float | Raster, *, steady_method: str = 'implicit', steady_tolerance: float = 1e-06, max_substeps: int = 10000000, max_newton_iterations: int = 200) -> dict

Find the steady water table under a constant recharge.

lsdtt-slope-stability ``dupuit_steady`` (and the antecedent state of
``dupuit_transient``). Raises RuntimeError if it does not converge.

Parameters
----------
recharge_m_per_s : float or Raster
    Recharge f (m/s), > 0 somewhere; nodata is 0.
steady_method : {"implicit", "explicit"}, default="implicit", keyword-only
    implicit: solved directly by Newton's method (pseudo-transient
    continuation with an active-set step for saturated and dry cells); tens
    of iterations whatever the grid. explicit: the time-dependent model
    marched until it stops changing; its steps are limited by the Courant
    and von Neumann conditions. Both solve the same discrete equations.
steady_tolerance : float, default=1e-06, keyword-only
    Converged when the residual is at most this times the largest recharge.
max_substeps : int, default=10000000, keyword-only
    explicit: most substeps.
max_newton_iterations : int, default=200, keyword-only
    implicit: most Newton (pseudo-time) steps.

Returns
-------
dict
    ``steady_method``, ``steps``, ``residual``, ``pseudo_time_s``; implicit
    also ``linear_iterations``, ``used_direct_solver``, ``n_saturated``,
    ``n_dry``. The implicit solve's budget is in ``steady_rates()``.

```python
def run_to_steady(self, recharge_m_per_s: float | Raster, *, steady_method: str = 'implicit', steady_tolerance: float = 1e-06, max_substeps: int = 10000000, max_newton_iterations: int = 200) -> dict
```

#### `saturated_thickness`

saturated_thickness(self) -> Raster

Saturated thickness h (m).

Parameters
----------
self : GroundwaterDupuit
    The aquifer.

Returns
-------
Raster
    h (m); nodata outside the aquifer.

```python
def saturated_thickness(self) -> Raster
```

#### `seepage`

seepage(self) -> Raster

Seepage rate (m/s) from the last call.

The implicit steady state's seepage, or the mean over the last ``run``.

Parameters
----------
self : GroundwaterDupuit
    The aquifer.

Returns
-------
Raster
    Seepage (m/s); nodata outside the aquifer.

```python
def seepage(self) -> Raster
```

#### `set_thickness`

set_thickness(self, thickness: Raster) -> None

Set the saturated thickness h of the active cells (clipped to [0, d]).

Parameters
----------
thickness : Raster
    Saturated thickness (m); nodata is 0.

Returns
-------
None

```python
def set_thickness(self, thickness: Raster) -> None
```

#### `steady_rates`

steady_rates(self) -> dict

The budget of the last implicit steady solve, as rates (m^3/s).

Parameters
----------
self : GroundwaterDupuit
    The aquifer.

Returns
-------
dict
    ``recharge_m3_per_s, outflow_fixed_head_m3_per_s, seepage_m3_per_s,
    imbalance_m3_per_s``.

```python
def steady_rates(self) -> dict
```

#### `storage`

Water stored, n h summed over the active cells (m^3).

#### `water_table_depth`

water_table_depth(self) -> Raster

Depth of the water table below the surface (m): the input to ``FlowInfo.rupture_planes``.

Parameters
----------
self : GroundwaterDupuit
    The aquifer.

Returns
-------
Raster
    Water-table depth (m); 0 where saturated; nodata outside the aquifer.

```python
def water_table_depth(self) -> Raster
```

---

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
moving_average_fits : dict
    Window (pixels) to (fitted, detrended) elevations per node from ``fit_profile``'s moving averages; nodata before the snout.
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
fit_profile(use_polynomial: bool = False, use_constant: bool = False, detrend_moving_average_window_pixels: Sequence[int] = []) -> None
    Fit the ridgeline profile beyond the snout (slope and r^2 from the linear fit).
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

fit_profile(self, use_polynomial: bool = False, use_constant: bool = False, detrend_moving_average_window_pixels: collections.abc.Sequence[int] = []) -> None

Fit the ridgeline profile beyond the snout (slope and r^2 from the linear fit).

Parameters
----------
use_polynomial, use_constant : bool, default=False
    Detrend against a quadratic, or a constant, instead of the line.
detrend_moving_average_window_pixels : list of int, default=[]
    Also detrend against moving averages of these widths (positive odd
    pixel counts, each once), as lsdtt-ridgeline-analysis
    ``ridgeline_detrend_moving_average_window_pixels``; read them from
    ``moving_average_fits``.

```python
def fit_profile(self, use_polynomial: bool = False, use_constant: bool = False, detrend_moving_average_window_pixels: collections.abc.Sequence[int] = []) -> None
```

#### `is_four_connected`

is_four_connected(self, periodic_ns: bool = False, periodic_ew: bool = False) -> bool

```python
def is_four_connected(self, periodic_ns: bool = False, periodic_ew: bool = False) -> bool
```

#### `length_metres`

(self) -> float

#### `moving_average_fits`

Window (pixels) to (fitted, detrended) elevations per node from ``fit_profile``'s moving averages; nodata before the snout.

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
degrees from grid north (true north with
``bearings_relative_to_true_north=True``; the angle between them is the
same either way); slopes are m/m.

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
along-feature azimuth (0-180, from grid north, or true north when
``find_roads`` was called with ``bearings_relative_to_true_north=True``) and
``slope`` the gradient (m/m).

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

## `RupturePlanes`

Planar rupture surfaces (scaLr) on a DEM, from ``FlowInfo.rupture_planes``.

Holds the inputs of the run, so its probability of failure, landslide
catalogues and slip surfaces are computed from the same flow, geometry,
channels, regions and water table.

Members: ``safety_factor``, ``depth``, ``rupture_slope``,
``topographic_slope``, ``length``, ``unstable``, ``channel``, ``daylight``,
``n_sources``, ``n_unstable``, ``n_channel``, ``failure_probability``,
``landslide_catalogues``, ``slip_surfaces``.

Returned by ``FlowInfo.rupture_planes()``; not constructible directly.

Attributes
----------
channel : Raster
    The channel mask used (1 channel, 0 hillslope).
daylight : object
    Flat index of each cell's daylight point (numpy int64, -1 where none).
depth : Raster
    Optimal rupture depth h* below the source cell (m).
length : Raster
    Horizontal distance dx from the source to its daylight point (m).
n_channel : int
    Channel cells (not sources).
n_sources : int
    Hillslope cells with an admissible daylight point.
n_unstable : int
    Of the sources, those with F <= 1.
rupture_slope : Raster
    Gradient of the rupture plane, (dz - h*) / dx (m/m).
safety_factor : Raster
    Factor of safety F at the critical daylight point (the minimum over the path); nodata on channels and cells with no admissible plane.
topographic_slope : Raster
    Gradient from the source to its daylight point, dz / dx (m/m).
unstable : Raster
    1 where F <= 1, 0 where F > 1.

Methods
-------
failure_probability(*, cohesion_distribution: str = 'uniform', cohesion_pa_min: float = 500.0, cohesion_pa_max: float = 1500.0, cohesion_pa_sd: float = 250.0, friction_coefficient_distribution: str = 'uniform', friction_coefficient_min: float = 0.6, friction_coefficient_max: float = 0.8, friction_coefficient_sd: float = 0.05, uncertainty_samples_per_parameter: int = 5) -> dict
    Probability of failure of the rupture planes by equal-density sampling of cohesion and friction.
landslide_catalogues(*, min_daylight_neighbours: int = 2) -> dict
    The potential-landslide catalogues: by daylight point and amalgamated.
slip_surfaces(*, family: str = 'circular', method: str = 'auto', search: str = 'seeded', back_scarp: str = 'vertical', janbu_correction: bool = True, structural_dip_gradient: float | Raster = 0.3, structural_dip_direction_degrees: float | Raster = 0.0, column_method: str = 'xie_revised', ellipsoid_search: str = 'seeded', ellipsoid_random_draws: int = 200, ellipsoid_random_seed: int = 1, ellipsoid_min_half_axis_pixels: float = 3.0, layer_depths: list[Raster] | None = None, screen_safety_factor: float = 1.5, toe_window_pixels: int = 0, max_length_metres: float = 0.0, slices_per_cell: int = 4, search_points: int = 12) -> SlipSurfaces
    Critical slip surfaces seeded from the rupture planes.

Authors
-------
B. Gailleton

### Members

#### `channel`

The channel mask used (1 channel, 0 hillslope).

#### `daylight`

Flat index of each cell's daylight point (numpy int64, -1 where none).

#### `depth`

Optimal rupture depth h* below the source cell (m).

#### `failure_probability`

failure_probability(self, *, cohesion_distribution: str = 'uniform', cohesion_pa_min: float = 500.0, cohesion_pa_max: float = 1500.0, cohesion_pa_sd: float = 250.0, friction_coefficient_distribution: str = 'uniform', friction_coefficient_min: float = 0.6, friction_coefficient_max: float = 0.8, friction_coefficient_sd: float = 0.05, uncertainty_samples_per_parameter: int = 5) -> dict

Probability of failure of the rupture planes by equal-density sampling of cohesion and friction.

lsdtt-slope-stability ``parameter_uncertainty: equal_density``, the same routine. For each cell and
each (c, mu) sample the search is the full deterministic one, so the daylight point may differ
between samples; Pf = sum of the weights of the samples with F_min <= 1. Exact for this model
(no Monte Carlo error). Density and the water table are held fixed.

Parameters
----------
cohesion_distribution : {"uniform", "normal"}, default="uniform", keyword-only
    uniform on [cohesion_pa_min, cohesion_pa_max], or normal about the central cohesion with
    cohesion_pa_sd, truncated to [max(0, c - 3 sd), c + 3 sd].
cohesion_pa_min, cohesion_pa_max, cohesion_pa_sd : float, default=500.0, 1500.0, 250.0, keyword-only
    Bounds (uniform) or standard deviation (normal) of cohesion (Pa).
friction_coefficient_distribution : {"uniform", "normal"}, default="uniform", keyword-only
    As cohesion_distribution, for the friction coefficient.
friction_coefficient_min, friction_coefficient_max, friction_coefficient_sd : float, default=0.6, 0.8, 0.05, keyword-only
    Bounds or standard deviation of the friction coefficient.
uncertainty_samples_per_parameter : int, default=5, keyword-only
    Equal-density samples of each parameter; their tensor product is used.

A region's rupture_cohesion_pa_min (and the other bounds and SDs) in the
region_parameters given to ``FlowInfo.rupture_planes`` override these.

Returns
-------
dict
    ``failure_probability`` (Raster, P(F <= 1)), ``mean_safety_factor``,
    ``sd_safety_factor`` (Rasters, weighted over the samples), ``n_sources``,
    ``n_possible_failure`` (Pf > 0), ``n_certain_failure`` (Pf = 1), ``n_samples``.

```python
def failure_probability(self, *, cohesion_distribution: str = 'uniform', cohesion_pa_min: float = 500.0, cohesion_pa_max: float = 1500.0, cohesion_pa_sd: float = 250.0, friction_coefficient_distribution: str = 'uniform', friction_coefficient_min: float = 0.6, friction_coefficient_max: float = 0.8, friction_coefficient_sd: float = 0.05, uncertainty_samples_per_parameter: int = 5) -> dict
```

#### `landslide_catalogues`

landslide_catalogues(self, *, min_daylight_neighbours: int = 2) -> dict

The potential-landslide catalogues: by daylight point and amalgamated.

lsdtt-slope-stability ``write_landslide_catalogue`` (scaLr's MATLAB): the
unstable cells sharing a daylight point form one landslide (after
re-attaching each daylight point to its source), the 1-D rule drops daylight
points with fewer than ``min_daylight_neighbours`` daylight-point neighbours,
and 8-connected landslides are merged into the amalgamated catalogue. Areas,
depths, volumes, F_L (eq. 7) and the planform ellipse are the CLI's columns.

Parameters
----------
min_daylight_neighbours : int, default=2, keyword-only
    The 1-D rule (CLI ``rupture_min_daylight_neighbours``); 0 keeps every
    daylight point (the paper's catalogue), 0 to 8.

Returns
-------
dict
    ``by_daylight`` and ``amalgamated`` (DataFrames, the CLI's
    ``_landslide_catalogue*.csv`` columns), ``label`` and
    ``label_amalgamated`` (Rasters of the ids), ``n_daylight_points``,
    ``n_daylight_points_removed``, ``n_landslide_cells``.

```python
def landslide_catalogues(self, *, min_daylight_neighbours: int = 2) -> dict
```

#### `length`

Horizontal distance dx from the source to its daylight point (m).

#### `n_channel`

Channel cells (not sources).

#### `n_sources`

Hillslope cells with an admissible daylight point.

#### `n_unstable`

Of the sources, those with F <= 1.

#### `rupture_slope`

Gradient of the rupture plane, (dz - h*) / dx (m/m).

#### `safety_factor`

Factor of safety F at the critical daylight point (the minimum over the path); nodata on channels and cells with no admissible plane.

#### `slip_surfaces`

slip_surfaces(self, *, family: str = 'circular', method: str = 'auto', search: str = 'seeded', back_scarp: str = 'vertical', janbu_correction: bool = True, structural_dip_gradient: float | Raster = 0.3, structural_dip_direction_degrees: float | Raster = 0.0, column_method: str = 'xie_revised', ellipsoid_search: str = 'seeded', ellipsoid_random_draws: int = 200, ellipsoid_random_seed: int = 1, ellipsoid_min_half_axis_pixels: float = 3.0, layer_depths: list[Raster] | None = None, screen_safety_factor: float = 1.5, toe_window_pixels: int = 0, max_length_metres: float = 0.0, slices_per_cell: int = 4, search_points: int = 12) -> SlipSurfaces

Critical slip surfaces seeded from the rupture planes.

lsdtt-slope-stability ``slip_surface_family``, the same routines; the keyword
names are the CLI's less the ``slip_surface_`` prefix. Crowns whose rupture
F is above ``screen_safety_factor`` are screened out. The 2-D families
(method of slices on the crown's D8 section) have no side resistance;
``ellipsoid`` (method of columns) does.

Parameters
----------
family : {"plane_integrated", "circular", "bilinear", "structural", "ellipsoid"}, default="circular", keyword-only
    plane_integrated: scaLr's plane with the forces summed over the body;
    circular: arcs by sagitta; bilinear: a Rankine back segment onto a
    flatter basal plane; structural: a plane at the apparent dip of a weak
    layer, ending where it daylights; ellipsoid: 3-D ellipsoids through the
    crown and the daylight point.
method : {"auto", "ordinary", "bishop", "janbu", "spencer"}, default="auto", keyword-only
    Limit-equilibrium method of the 2-D families; auto is ordinary for
    plane_integrated, bishop for circular, janbu for bilinear and structural.
search : {"seeded", "full"}, default="seeded", keyword-only
    seeded: the toe at the rupture daylight point (+/- toe_window_pixels);
    full: every toe down the path (up to max_length_metres).
back_scarp : {"vertical", "rankine"}, default="vertical", keyword-only
    The crown's back scarp: a tension crack, or a Rankine active plane
    carrying cohesion and friction (not for bilinear).
janbu_correction : bool, default=True, keyword-only
    Apply Janbu's f0.
structural_dip_gradient : float or Raster, default=0.3, keyword-only
    structural: tan(dip) of the weak layer, constant or per cell (at the crown).
structural_dip_direction_degrees : float or Raster, default=0.0, keyword-only
    structural: dip direction, degrees clockwise from north.
column_method : {"xie_revised", "hovland"}, default="xie_revised", keyword-only
    ellipsoid: the method of columns.
ellipsoid_search : {"seeded", "random"}, default="seeded", keyword-only
    ellipsoid: a grid of shapes, or random draws.
ellipsoid_random_draws : int, default=200, keyword-only
    ellipsoid random search: shapes per crown.
ellipsoid_random_seed : int, default=1, keyword-only
    ellipsoid random search: seed.
ellipsoid_min_half_axis_pixels : float, default=3.0, keyword-only
    ellipsoid: smallest length half-axis (cells); shorter crowns are "too small".
layer_depths : list of Raster or None, default=None, keyword-only
    ellipsoid: depths (m) of weak layers that truncate the ellipsoids; each
    is searched separately.
screen_safety_factor : float, default=1.5, keyword-only
    Crowns with rupture F above this are not evaluated.
toe_window_pixels : int, default=0, keyword-only
    Seeded search: toes this many cells either side of the daylight point.
max_length_metres : float, default=0.0, keyword-only
    Full search: longest crown-to-toe flow distance (0 = no cap).
slices_per_cell : int, default=4, keyword-only
    Slices per section cell.
search_points : int, default=12, keyword-only
    Points per shape parameter in the search.

Returns
-------
SlipSurfaces
    The critical surface of each evaluated crown, its bodies and counts.

```python
def slip_surfaces(self, *, family: str = 'circular', method: str = 'auto', search: str = 'seeded', back_scarp: str = 'vertical', janbu_correction: bool = True, structural_dip_gradient: float | Raster = 0.3, structural_dip_direction_degrees: float | Raster = 0.0, column_method: str = 'xie_revised', ellipsoid_search: str = 'seeded', ellipsoid_random_draws: int = 200, ellipsoid_random_seed: int = 1, ellipsoid_min_half_axis_pixels: float = 3.0, layer_depths: list[Raster] | None = None, screen_safety_factor: float = 1.5, toe_window_pixels: int = 0, max_length_metres: float = 0.0, slices_per_cell: int = 4, search_points: int = 12) -> SlipSurfaces
```

#### `topographic_slope`

Gradient from the source to its daylight point, dz / dx (m/m).

#### `unstable`

1 where F <= 1, 0 where F > 1.

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

## `SlipSurfaces`

Critical slip surfaces of each crown, from ``RupturePlanes.slip_surfaces``.

Members: ``family``, ``method``, ``safety_factor``, ``crown_depth``,
``max_depth``, ``sagitta``, ``basal_gradient``, ``length``, ``sinuosity``,
``unstable``, ``evaluated``, ``body_depth``, ``cover_safety_factor``,
``n_evaluated``, ``n_screened``, ``n_unstable``, ``n_inadmissible``,
``n_too_small``, ``failure_probability``, ``spencer_check``,
``body_catalogues``.

Returned by ``RupturePlanes.slip_surfaces()``; not constructible directly.

Attributes
----------
basal_gradient : Raster
    bilinear: tan of the basal plane; structural: the apparent dip gradient.
body_depth : Raster
    Envelope of the unstable bodies' thickness (m).
cover_safety_factor : Raster
    ellipsoid: lowest F of the critical surfaces covering each cell.
crown_depth : Raster
    Depth h below the crown (m); ellipsoid: the depth half-axis.
evaluated : Raster
    1 evaluated, 0 screened out (crowns with a rupture result).
family : str
    The slip-surface family.
length : Raster
    Flow distance crown to toe (m); ellipsoid: straight distance.
max_depth : Raster
    Largest body thickness (m).
method : str
    The limit-equilibrium method used (auto resolved), or the column method for ellipsoids.
n_evaluated : int
    Crowns evaluated.
n_inadmissible : int
    Crowns with no admissible surface.
n_screened : int
    Crowns screened out by screen_safety_factor.
n_too_small : int
    ellipsoid: crowns below ellipsoid_min_half_axis_pixels (0 for the 2-D families).
n_unstable : int
    Crowns with F <= 1.
safety_factor : Raster
    Factor of safety of each evaluated crown's critical surface.
sagitta : Raster
    circular: the sagitta s (m); ellipsoid: the lateral half-axis (m).
sinuosity : Raster
    2-D: flow distance over straight distance crown to toe.
unstable : Raster
    1 where F <= 1.

Methods
-------
body_catalogues() -> dict
    Catalogues of the unstable slip bodies: by toe and amalgamated.
failure_probability(*, cohesion_distribution: str = 'uniform', cohesion_pa_min: float = 500.0, cohesion_pa_max: float = 1500.0, cohesion_pa_sd: float = 250.0, friction_coefficient_distribution: str = 'uniform', friction_coefficient_min: float = 0.6, friction_coefficient_max: float = 0.8, friction_coefficient_sd: float = 0.05, uncertainty_samples_per_parameter: int = 5) -> dict
    Probability of failure on each crown's fixed critical surface.
spencer_check() -> dict
    Spencer's factor of safety on each critical 2-D surface.

Authors
-------
B. Gailleton

### Members

#### `basal_gradient`

bilinear: tan of the basal plane; structural: the apparent dip gradient.

#### `body_catalogues`

body_catalogues(self) -> dict

Catalogues of the unstable slip bodies: by toe and amalgamated.

lsdtt-slope-stability ``write_slip_body_catalogue``. By toe: the union of the
bodies of the unstable crowns sharing a toe (groups may overlap, so their
volumes need not add up). Amalgamated: 8-connected components of the body
envelope, which do not overlap.

Parameters
----------
self : SlipSurfaces
    The slip surfaces.

Returns
-------
dict
    ``by_toe`` and ``amalgamated`` (DataFrames, the CLI's
    ``_slip_body_catalogue*.csv`` columns) and ``label_amalgamated`` (Raster).

```python
def body_catalogues(self) -> dict
```

#### `body_depth`

Envelope of the unstable bodies' thickness (m).

#### `cover_safety_factor`

ellipsoid: lowest F of the critical surfaces covering each cell.

#### `crown_depth`

Depth h below the crown (m); ellipsoid: the depth half-axis.

#### `evaluated`

1 evaluated, 0 screened out (crowns with a rupture result).

#### `failure_probability`

failure_probability(self, *, cohesion_distribution: str = 'uniform', cohesion_pa_min: float = 500.0, cohesion_pa_max: float = 1500.0, cohesion_pa_sd: float = 250.0, friction_coefficient_distribution: str = 'uniform', friction_coefficient_min: float = 0.6, friction_coefficient_max: float = 0.8, friction_coefficient_sd: float = 0.05, uncertainty_samples_per_parameter: int = 5) -> dict

Probability of failure on each crown's fixed critical surface.

lsdtt-slope-stability ``parameter_uncertainty`` with a slip family. The surface found with the
central values is held fixed and evaluated for each (c, mu) sample, so Pf is a lower bound on the
probability that some surface fails. Screened-out crowns get 0. For ellipsoids the result also
has ``cover_failure_probability``, the highest Pf of the ellipsoids covering each cell.

Parameters
----------
cohesion_distribution : {"uniform", "normal"}, default="uniform", keyword-only
    uniform on [cohesion_pa_min, cohesion_pa_max], or normal about the central cohesion with
    cohesion_pa_sd, truncated to [max(0, c - 3 sd), c + 3 sd].
cohesion_pa_min, cohesion_pa_max, cohesion_pa_sd : float, default=500.0, 1500.0, 250.0, keyword-only
    Bounds (uniform) or standard deviation (normal) of cohesion (Pa).
friction_coefficient_distribution : {"uniform", "normal"}, default="uniform", keyword-only
    As cohesion_distribution, for the friction coefficient.
friction_coefficient_min, friction_coefficient_max, friction_coefficient_sd : float, default=0.6, 0.8, 0.05, keyword-only
    Bounds or standard deviation of the friction coefficient.
uncertainty_samples_per_parameter : int, default=5, keyword-only
    Equal-density samples of each parameter; their tensor product is used.

A region's rupture_cohesion_pa_min (and the other bounds and SDs) in the
region_parameters given to ``FlowInfo.rupture_planes`` override these.

Returns
-------
dict
    ``failure_probability`` (Raster, P(F <= 1)), ``mean_safety_factor``,
    ``sd_safety_factor`` (Rasters, weighted over the samples), ``n_sources``,
    ``n_possible_failure`` (Pf > 0), ``n_certain_failure`` (Pf = 1), ``n_samples``.

```python
def failure_probability(self, *, cohesion_distribution: str = 'uniform', cohesion_pa_min: float = 500.0, cohesion_pa_max: float = 1500.0, cohesion_pa_sd: float = 250.0, friction_coefficient_distribution: str = 'uniform', friction_coefficient_min: float = 0.6, friction_coefficient_max: float = 0.8, friction_coefficient_sd: float = 0.05, uncertainty_samples_per_parameter: int = 5) -> dict
```

#### `family`

The slip-surface family.

#### `length`

Flow distance crown to toe (m); ellipsoid: straight distance.

#### `max_depth`

Largest body thickness (m).

#### `method`

The limit-equilibrium method used (auto resolved), or the column method for ellipsoids.

#### `n_evaluated`

Crowns evaluated.

#### `n_inadmissible`

Crowns with no admissible surface.

#### `n_screened`

Crowns screened out by screen_safety_factor.

#### `n_too_small`

ellipsoid: crowns below ellipsoid_min_half_axis_pixels (0 for the 2-D families).

#### `n_unstable`

Crowns with F <= 1.

#### `safety_factor`

Factor of safety of each evaluated crown's critical surface.

#### `sagitta`

circular: the sagitta s (m); ellipsoid: the lateral half-axis (m).

#### `sinuosity`

2-D: flow distance over straight distance crown to toe.

#### `spencer_check`

spencer_check(self) -> dict

Spencer's factor of safety on each critical 2-D surface.

lsdtt-slope-stability ``slip_surface_spencer_check``: force and moment
equilibrium with a constant interslice inclination, on the surface the
family's own method found, as a check of the simplified methods. Not for
ellipsoids.

Parameters
----------
self : SlipSurfaces
    A 2-D family's result.

Returns
-------
dict
    ``safety_factor`` (Raster, nodata where Spencer has no solution) and
    ``n_no_solution``.

```python
def spencer_check(self) -> dict
```

#### `unstable`

1 where F <= 1.

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
bearings_relative_to_true_north : bool, default=False
    Convert ``regional_slope_bearing`` and ``flow_path_bearing`` to true
    north at each pixel, as the CLI's ``bearings_relative_to_true_north``;
    needs a DEM with a CRS.

Returns
-------
RoadFlowResult

Authors
-------
B. Gailleton
S. M. Mudd

```python
def find_road_flow_discrimination(dem: lsdtt3.Raster, regional_smoothing_metres: float = 100.0, threshold_contributing_pixels: int = 1000, flow_path_half_window: int = 5, intermediate_window_metres: float = 15.0, roughness_radius_metres: float = 2.5, raster_is_filled: bool = False, boundary_conditions: str = 'oooo', n_threads: int = 0, bearings_relative_to_true_north: bool = False) -> lsdtt3.RoadFlowResult
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
bearings_relative_to_true_north : bool, default=False
    Convert ``feature_orientation_degrees`` to true north at each pixel, as
    the CLI's ``bearings_relative_to_true_north``; needs a DEM with a CRS.

Returns
-------
RoadResult

Authors
-------
B. Gailleton
S. M. Mudd

```python
def find_roads(dem: lsdtt3.Raster, smoothing_window_metres: float = 8.0, score_threshold: float = 0.5, use_bench_test: bool = True, bench_half_width_metres: float = 6.0, n_threads: int = 0, bearings_relative_to_true_north: bool = False) -> lsdtt3.RoadResult
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

## `landslide_frequency_area`

Frequency - area of a landslide catalogue and its power-law tail.

lsdtt-slope-stability's ``_frequency_area.csv`` and ``_frequency_area_fit.csv``
for one catalogue (landslides or slip bodies). Areas are whole numbers of
cells, so the fit and the bins are discrete in cells. The tail p(A) ~ A^eta
is fitted by maximum likelihood with A_min chosen by the Kolmogorov -
Smirnov distance (Clauset, Shalizi and Newman 2009); ``binned_eta`` is the
least-squares slope of log p(A) against log A over the bins at or above
A_min holding at least 3 landslides.

Parameters
----------
areas_m2 : sequence of float
    Landslide areas (m^2), e.g. a catalogue's ``area_m2`` column.
cell_area_m2 : float
    Area of one cell (m^2).
bins_per_decade : int, default=10, keyword-only
    Logarithmic bins per decade of area (CLI ``frequency_area_bins_per_decade``).
min_tail : int, default=50, keyword-only
    Fewest landslides in the fitted tail (at least 2; CLI ``frequency_area_min_tail``).

Returns
-------
tuple of (pandas.DataFrame, dict)
    The density (``area_lower_m2, area_upper_m2, area_centre_m2,
    n_landslides, probability_density_per_m2``) and the fit (``n_landslides,
    method, area_min_m2, eta, eta_standard_error, n_tail, ks_distance,
    binned_eta, binned_eta_standard_error, n_bins_fitted``), the CLI's columns
    less ``catalogue``.

Authors
-------
B. Gailleton

```python
def landslide_frequency_area(areas_m2: Sequence[float], cell_area_m2: float, *, bins_per_decade: int = 10, min_tail: int = 50) -> tuple[pandas.DataFrame, dict]
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

## `sinmap`

SINMAP stability index, class, probability and factor-of-safety bounds.

lsdtt-slope-stability ``run_sinmap`` (Pack et al. 1998), the same routine;
the keyword names are the CLI's less the ``sinmap_`` prefix. The infinite
slope with steady wetness w = min(1, R a / (T sin theta)) is evaluated over
uniform ranges of dimensionless cohesion C, friction coefficient tan phi and
T/R. FS_min uses (C_min, t_min, (T/R)_min), FS_max the other bounds; the
stability index is FS_min where it is above 1, P(FS > 1) where FS_min < 1 <
FS_max (closed form), and 0 where FS_max < 1. Flat cells are class 1 and
nodata in the continuous outputs. The CLI uses the D-infinity slope and area
of the carved-then-filled DEM: ``FlowGraph(dem, method="dinf")``,
``.slope()`` and ``.drainage_area()``.

Parameters
----------
slope : Raster
    Gradient S (m/m), e.g. ``FlowGraph.slope()``.
drainage_area : Raster
    Contributing area (m^2) on the same grid, e.g. ``FlowGraph.drainage_area()``.
dimensionless_cohesion_min, dimensionless_cohesion_max : float, default=0.0, 0.25, keyword-only
    C = (C_r + C_s) / (h rho_s g), bounds.
friction_coefficient_min, friction_coefficient_max : float, default=tan 30, 1.0, keyword-only
    tan phi, bounds.
transmissivity_over_recharge_min_metres, transmissivity_over_recharge_max_metres : float, default=2000.0, 3000.0, keyword-only
    T / R (m), bounds.
water_to_soil_density_ratio : float, default=0.5, keyword-only
    rho_w / rho_s.
region : Raster or None, default=None, keyword-only
    Integer region ids selecting parameter sets; nodata cells are left out.
region_parameters : dict or None, default=None, keyword-only
    ``{region_id: {name: value}}`` with the CLI's region-table column names
    (``sinmap_dimensionless_cohesion_min``, ...); a missing name takes the
    keyword value. ``rupture_*`` names are accepted and ignored, so one table
    can drive both ``sinmap`` and ``FlowInfo.rupture_planes``. Needed with
    ``region``.
exclude : Raster or None, default=None, keyword-only
    Nonzero cells are left out, e.g. ``FlowGraph.dinf_edge_contaminated_mask()``.
contour_width : Raster or None, default=None, keyword-only
    Contour width (m) for a = A / w, e.g. ``FlowGraph.dinf_contour_width()``
    (the CLI's default); None uses the cell size, as SINMAP.
n_threads : int, default=0, keyword-only
    OpenMP threads (0 = default); the result does not depend on it.

Returns
-------
dict
    Rasters ``stability_index``, ``stability_class`` (1 stable to 6 defended),
    ``probability_stable``, ``fs_min``, ``fs_max``, ``saturation`` (SINMAP's
    coding: 3 always saturated, 2 partly, else the wettest w),
    ``specific_catchment_area`` (m); counts ``n_evaluated``, ``n_flat``,
    ``n_excluded``; ``class_areas``, a DataFrame of the CLI's class-area
    table (``region_id, stability_class, n_cells, area_m2,
    fraction_of_region``).

Authors
-------
B. Gailleton

```python
def sinmap(slope: Raster, drainage_area: Raster, *, dimensionless_cohesion_min: float = 0.0, dimensionless_cohesion_max: float = 0.25, friction_coefficient_min: float = 0.5773502691896257, friction_coefficient_max: float = 1.0, transmissivity_over_recharge_min_metres: float = 2000.0, transmissivity_over_recharge_max_metres: float = 3000.0, water_to_soil_density_ratio: float = 0.5, region: Raster | None = None, region_parameters: dict | None = None, exclude: Raster | None = None, contour_width: Raster | None = None, n_threads: int = 0) -> dict
```

---

## `sinmap_slope_area_curves`

The class boundaries of SINMAP's slope - specific catchment area plot.

lsdtt-slope-stability ``write_sinmap_slope_area_curves_table`` for one
parameter set: the saturation lines (``always_saturated_above``,
``never_saturated_below``) and the FS curves ``fs_min_1.5``, ``fs_min_1.0``
and ``fs_max_1.0`` as points. The CLI uses the largest evaluated slope,
the cell size and the largest specific catchment area as the extent.

Parameters
----------
slope_max : float
    Largest gradient (m/m) of the plot.
area_min, area_max : float
    Range of specific catchment area (m); points outside it are dropped.
dimensionless_cohesion_min, dimensionless_cohesion_max, friction_coefficient_min, friction_coefficient_max, transmissivity_over_recharge_min_metres, transmissivity_over_recharge_max_metres, water_to_soil_density_ratio
    As ``sinmap``, keyword-only.
n_slopes : int, default=400, keyword-only
    Slopes sampled along each curve.

Returns
-------
pandas.DataFrame
    ``curve, slope, specific_catchment_area_m``, the CLI table's columns
    less ``region_id``.

Authors
-------
B. Gailleton

```python
def sinmap_slope_area_curves(slope_max: float, area_min: float, area_max: float, *, dimensionless_cohesion_min: float = 0.0, dimensionless_cohesion_max: float = 0.25, friction_coefficient_min: float = 0.5773502691896257, friction_coefficient_max: float = 1.0, transmissivity_over_recharge_min_metres: float = 2000.0, transmissivity_over_recharge_max_metres: float = 3000.0, water_to_soil_density_ratio: float = 0.5, n_slopes: int = 400) -> pandas.DataFrame
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
