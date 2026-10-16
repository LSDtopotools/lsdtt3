# Landscape And Hillslopes API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

LandscapeModel thresholds use specific catchment area A/dx in metres. Flat relief windows use a full side length in pixels.

## `HillslopeCellularResult`

Conservative cellular hillslope-transport result and mass-transfer diagnostics.

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
    Conservative cellular hillslope-transport result and mass-transfer diagnostics.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


### Members

#### `completed`

(self) -> bool

#### `elevation`

(self) -> lsdtt3.Raster

#### `exported_volume_m3`

(self) -> float

#### `max_slope_ratio`

(self) -> float

#### `message`

(self) -> str

#### `n_face_events`

(self) -> int

#### `n_substeps`

(self) -> int

---

## `HillslopeSteadyStateResult`

Steady hillslope-surface solution and convergence diagnostics.

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
    Steady hillslope-surface solution and convergence diagnostics.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


### Members

#### `converged`

(self) -> bool

#### `elevation`

(self) -> lsdtt3.Raster

#### `final_residual`

(self) -> float

#### `n_iterations`

(self) -> int

---

## `HillslopeTransportResult`

Implicit hillslope-transport result and adaptive-solver diagnostics.

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
    Implicit hillslope-transport result and adaptive-solver diagnostics.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


### Members

#### `elevation`

(self) -> lsdtt3.Raster

#### `final_dt`

(self) -> float

#### `limiter_fired`

(self) -> bool

#### `max_slope_ratio`

(self) -> float

#### `n_rejected_substeps`

(self) -> int

#### `n_substeps`

(self) -> int

#### `total_picard_iterations`

(self) -> int

---

## `LandscapeModel`

Stateful coupled fluvial-incision and hillslope-transport landscape model.

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
    Stateful coupled fluvial-incision and hillslope-transport landscape model.

The model owns its evolving elevation raster and reuses native hillslope solver
workspaces between calls to ``step`` or ``run``.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


**Deprecated compatibility names:**
- `incision_threshold_area` -> `incision_threshold_specific_area`; earliest removal: following minor.

```python
def __init__(self, dem: lsdtt3.Raster, timestep_years: float = 100.0, erodibility: float = 0.0010000000474974513, m_exponent: float = 0.44999998807907104, n_exponent: float = 1.0, uplift_rate: float = 9.999999747378752e-05, hillslope_diffusivity: float = 0.0010000000474974513, hillslope_model: str = 'cidre', nonlinear_hillslopes: bool = False, critical_hillslope_slope: float = 0.800000011920929, incision_threshold_specific_area: float = 0.0, coupling_regime: str = 'automatic', boundary_conditions: str = 'oooo', fluvial_method: str = 'd8', channel_mask: object | None = None, hold_channel_mask_fixed: bool = False, fixed_pixels: object | None = None, n_threads: int = 0, baselevel_fall_rate: float = 0.0, cellular_r_match: float = 0.85, cellular_cfl_safety: float = 0.8, cellular_nonlocal_fraction: float = 0.0, cellular_nonlocal_activation_ratio: float = 0.8, cellular_nonlocal_deposition_probability: float = 0.5, cellular_nonlocal_max_hops: int = 4, erodibility_field: object | None = None, erodibility_m_exponent: object | None = None, erodibility_area_convention: object | None = None) -> None
def __init__(self, dem: lsdtt3.Raster, timestep_years: float = 100.0, erodibility: float = 0.0010000000474974513, m_exponent: float = 0.44999998807907104, n_exponent: float = 1.0, uplift_rate: float = 9.999999747378752e-05, hillslope_diffusivity: float = 0.0010000000474974513, hillslope_model: str = 'cidre', nonlinear_hillslopes: bool = False, critical_hillslope_slope: float = 0.800000011920929, incision_threshold_area: float = 0.0, coupling_regime: str = 'automatic', boundary_conditions: str = 'oooo', fluvial_method: str = 'd8', channel_mask: object | None = None, hold_channel_mask_fixed: bool = False, fixed_pixels: object | None = None, n_threads: int = 0, baselevel_fall_rate: float = 0.0, cellular_r_match: float = 0.85, cellular_cfl_safety: float = 0.8, cellular_nonlocal_fraction: float = 0.0, cellular_nonlocal_activation_ratio: float = 0.8, cellular_nonlocal_deposition_probability: float = 0.5, cellular_nonlocal_max_hops: int = 4, erodibility_field: object | None = None, erodibility_m_exponent: object | None = None, erodibility_area_convention: object | None = None) -> None
```

### Members

#### `__init__`

__init__(self, dem: lsdtt3.Raster, timestep_years: float = 100.0, erodibility: float = 0.0010000000474974513, m_exponent: float = 0.44999998807907104, n_exponent: float = 1.0, uplift_rate: float = 9.999999747378752e-05, hillslope_diffusivity: float = 0.0010000000474974513, hillslope_model: str = 'cidre', nonlinear_hillslopes: bool = False, critical_hillslope_slope: float = 0.800000011920929, incision_threshold_specific_area: float = 0.0, coupling_regime: str = 'automatic', boundary_conditions: str = 'oooo', fluvial_method: str = 'd8', channel_mask: object | None = None, hold_channel_mask_fixed: bool = False, fixed_pixels: object | None = None, n_threads: int = 0, baselevel_fall_rate: float = 0.0, cellular_r_match: float = 0.85, cellular_cfl_safety: float = 0.8, cellular_nonlocal_fraction: float = 0.0, cellular_nonlocal_activation_ratio: float = 0.8, cellular_nonlocal_deposition_probability: float = 0.5, cellular_nonlocal_max_hops: int = 4, erodibility_field: object | None = None, erodibility_m_exponent: object | None = None, erodibility_area_convention: object | None = None) -> None
__init__(self, dem: lsdtt3.Raster, timestep_years: float = 100.0, erodibility: float = 0.0010000000474974513, m_exponent: float = 0.44999998807907104, n_exponent: float = 1.0, uplift_rate: float = 9.999999747378752e-05, hillslope_diffusivity: float = 0.0010000000474974513, hillslope_model: str = 'cidre', nonlinear_hillslopes: bool = False, critical_hillslope_slope: float = 0.800000011920929, incision_threshold_area: float = 0.0, coupling_regime: str = 'automatic', boundary_conditions: str = 'oooo', fluvial_method: str = 'd8', channel_mask: object | None = None, hold_channel_mask_fixed: bool = False, fixed_pixels: object | None = None, n_threads: int = 0, baselevel_fall_rate: float = 0.0, cellular_r_match: float = 0.85, cellular_cfl_safety: float = 0.8, cellular_nonlocal_fraction: float = 0.0, cellular_nonlocal_activation_ratio: float = 0.8, cellular_nonlocal_deposition_probability: float = 0.5, cellular_nonlocal_max_hops: int = 4, erodibility_field: object | None = None, erodibility_m_exponent: object | None = None, erodibility_area_convention: object | None = None) -> None

Overloaded function.

1. ``__init__(self, dem: lsdtt3.Raster, timestep_years: float = 100.0, erodibility: float = 0.0010000000474974513, m_exponent: float = 0.44999998807907104, n_exponent: float = 1.0, uplift_rate: float = 9.999999747378752e-05, hillslope_diffusivity: float = 0.0010000000474974513, hillslope_model: str = 'cidre', nonlinear_hillslopes: bool = False, critical_hillslope_slope: float = 0.800000011920929, incision_threshold_specific_area: float = 0.0, coupling_regime: str = 'automatic', boundary_conditions: str = 'oooo', fluvial_method: str = 'd8', channel_mask: object | None = None, hold_channel_mask_fixed: bool = False, fixed_pixels: object | None = None, n_threads: int = 0, baselevel_fall_rate: float = 0.0, cellular_r_match: float = 0.85, cellular_cfl_safety: float = 0.8, cellular_nonlocal_fraction: float = 0.0, cellular_nonlocal_activation_ratio: float = 0.8, cellular_nonlocal_deposition_probability: float = 0.5, cellular_nonlocal_max_hops: int = 4, erodibility_field: object | None = None, erodibility_m_exponent: object | None = None, erodibility_area_convention: object | None = None) -> None``

Create a coupled landscape-evolution model from an in-memory DEM.

Parameters
----------
dem : Raster
    Initial, depression-conditioned elevation surface.
timestep_years : float, optional
    Duration of each call to ``step``.
erodibility, m_exponent, n_exponent : float, optional
    Stream-power incision parameters. ``erodibility`` must be > 0; it is also the
    reference K for ``erodibility_field``.
uplift_rate, hillslope_diffusivity : float, optional
    Uniform uplift in m yr⁻¹ and hillslope diffusivity in m² yr⁻¹.
hillslope_model : {"none", "linear", "roering", "roering_implicit", "cidre", "adi", "cellular"}, optional
    Hillslope process model. ``none`` disables hillslope transport.
nonlinear_hillslopes, critical_hillslope_slope : bool, float, optional
    Roering nonlinear-flux controls.
incision_threshold_specific_area : float, optional
    Channel threshold in specific catchment area A/dx (metres); zero supports fully coupled runs.
coupling_regime : {"automatic", "coupled", "channels_as_boundaries"}, optional
    Whether hillslope channels emerge or are imposed.
boundary_conditions : str, optional
    Four NESW open, no-flux, periodic, or base-level codes.
fluvial_method : {"d8", "dinf", "none"}, optional
    Fluvial incision routing method, or ``none`` for hillslope-only evolution.
channel_mask : Raster, optional
    Same-grid binary imposed channel raster.
hold_channel_mask_fixed : bool, optional
    Add imposed channel pixels to the fixed-elevation set.
fixed_pixels : array-like of int, optional
    Explicit N-by-2 ``(row, col)`` fixed-elevation cells. Defaults to open/base-level edges.
n_threads : int, optional
    OpenMP thread count; zero uses the runtime default.
baselevel_fall_rate : float, optional
    Downward motion of fixed base-level cells in m yr⁻¹.
cellular_r_match, cellular_cfl_safety, cellular_nonlocal_fraction : float, optional
    Cellular-hillslope calibration, stability, and optional nonlocal-transfer controls.
erodibility_field : Raster or numpy.ndarray, optional
    Per-pixel absolute erodibility K replacing ``erodibility`` in the fluvial step: a
    Raster strictly on the DEM grid or a ``(n_rows, n_cols)`` array (NaN only on DEM
    nodata). Every valid pixel needs a finite K > 0.
erodibility_m_exponent : float, optional
    The m the field was calibrated for; a different m raises ``ValueError``.
erodibility_area_convention : {"map_area", "specific_area"}, optional
    Area convention of the field; ``map_area`` K is converted by dx**m.

Returns
-------
LandscapeModel
    Ready-to-step evolving model.

Authors
-------
B. Gailleton


2. ``__init__(self, dem: lsdtt3.Raster, timestep_years: float = 100.0, erodibility: float = 0.0010000000474974513, m_exponent: float = 0.44999998807907104, n_exponent: float = 1.0, uplift_rate: float = 9.999999747378752e-05, hillslope_diffusivity: float = 0.0010000000474974513, hillslope_model: str = 'cidre', nonlinear_hillslopes: bool = False, critical_hillslope_slope: float = 0.800000011920929, incision_threshold_area: float = 0.0, coupling_regime: str = 'automatic', boundary_conditions: str = 'oooo', fluvial_method: str = 'd8', channel_mask: object | None = None, hold_channel_mask_fixed: bool = False, fixed_pixels: object | None = None, n_threads: int = 0, baselevel_fall_rate: float = 0.0, cellular_r_match: float = 0.85, cellular_cfl_safety: float = 0.8, cellular_nonlocal_fraction: float = 0.0, cellular_nonlocal_activation_ratio: float = 0.8, cellular_nonlocal_deposition_probability: float = 0.5, cellular_nonlocal_max_hops: int = 4, erodibility_field: object | None = None, erodibility_m_exponent: object | None = None, erodibility_area_convention: object | None = None) -> None``

```python
def __init__(self, dem: lsdtt3.Raster, timestep_years: float = 100.0, erodibility: float = 0.0010000000474974513, m_exponent: float = 0.44999998807907104, n_exponent: float = 1.0, uplift_rate: float = 9.999999747378752e-05, hillslope_diffusivity: float = 0.0010000000474974513, hillslope_model: str = 'cidre', nonlinear_hillslopes: bool = False, critical_hillslope_slope: float = 0.800000011920929, incision_threshold_specific_area: float = 0.0, coupling_regime: str = 'automatic', boundary_conditions: str = 'oooo', fluvial_method: str = 'd8', channel_mask: object | None = None, hold_channel_mask_fixed: bool = False, fixed_pixels: object | None = None, n_threads: int = 0, baselevel_fall_rate: float = 0.0, cellular_r_match: float = 0.85, cellular_cfl_safety: float = 0.8, cellular_nonlocal_fraction: float = 0.0, cellular_nonlocal_activation_ratio: float = 0.8, cellular_nonlocal_deposition_probability: float = 0.5, cellular_nonlocal_max_hops: int = 4, erodibility_field: object | None = None, erodibility_m_exponent: object | None = None, erodibility_area_convention: object | None = None) -> None
def __init__(self, dem: lsdtt3.Raster, timestep_years: float = 100.0, erodibility: float = 0.0010000000474974513, m_exponent: float = 0.44999998807907104, n_exponent: float = 1.0, uplift_rate: float = 9.999999747378752e-05, hillslope_diffusivity: float = 0.0010000000474974513, hillslope_model: str = 'cidre', nonlinear_hillslopes: bool = False, critical_hillslope_slope: float = 0.800000011920929, incision_threshold_area: float = 0.0, coupling_regime: str = 'automatic', boundary_conditions: str = 'oooo', fluvial_method: str = 'd8', channel_mask: object | None = None, hold_channel_mask_fixed: bool = False, fixed_pixels: object | None = None, n_threads: int = 0, baselevel_fall_rate: float = 0.0, cellular_r_match: float = 0.85, cellular_cfl_safety: float = 0.8, cellular_nonlocal_fraction: float = 0.0, cellular_nonlocal_activation_ratio: float = 0.8, cellular_nonlocal_deposition_probability: float = 0.5, cellular_nonlocal_max_hops: int = 4, erodibility_field: object | None = None, erodibility_m_exponent: object | None = None, erodibility_area_convention: object | None = None) -> None
```

#### `channel_mask`

channel_mask(self) -> lsdtt3.Raster

Return the imposed or currently threshold-derived channel mask.

Returns
-------
Raster
    Binary channel raster on the current model grid.

Authors
-------
B. Gailleton

```python
def channel_mask(self) -> lsdtt3.Raster
```

#### `coupling_regime`

Resolved coupling experiment: coupled or channels_as_boundaries.

#### `elevation`

(self) -> lsdtt3.Raster

#### `erodibility_field`

The per-pixel erodibility in use (run convention), or None.

#### `hillslope_model`

(self) -> str

#### `iteration`

(self) -> int

#### `run`

run(self, n_steps: int, compute_erosion_rate: bool = False) -> lsdtt3.LandscapeStepResult

Advance the coupled model repeatedly and return the final step result.

Parameters
----------
n_steps : int
    Number of timesteps to execute.
compute_erosion_rate : bool, optional
    Include a per-pixel, signed erosion-rate raster for the final step.

Returns
-------
LandscapeStepResult
    Final surface and diagnostics.

Authors
-------
B. Gailleton

```python
def run(self, n_steps: int, compute_erosion_rate: bool = False) -> lsdtt3.LandscapeStepResult
```

#### `set_erodibility_field`

set_erodibility_field(self, field: object | None) -> None

Replace the per-pixel erodibility used by subsequent steps.

Parameters
----------
field : Raster, numpy.ndarray or None
    Absolute K strictly on the model grid (Raster), or a ``(n_rows, n_cols)`` array
    with NaN only on DEM nodata; ``None`` returns to the scalar ``erodibility``.
    Every valid pixel needs a finite K > 0.

Returns
-------
None

Authors
-------
B. Gailleton

```python
def set_erodibility_field(self, field: object | None) -> None
```

#### `set_forcing`

set_forcing(self, uplift_rate: float, baselevel_fall_rate: float = 0.0) -> None

Change the uniform tectonic/base-level forcing used by subsequent steps.

Parameters
----------
uplift_rate, baselevel_fall_rate : float
    Non-negative uplift and fixed-base-level lowering rates in m yr⁻¹.

Returns
-------
None

Authors
-------
B. Gailleton

```python
def set_forcing(self, uplift_rate: float, baselevel_fall_rate: float = 0.0) -> None
```

#### `step`

step(self, compute_erosion_rate: bool = False) -> lsdtt3.LandscapeStepResult

Advance the coupled model by one configured timestep.

Parameters
----------
compute_erosion_rate : bool, optional
    Include a per-pixel, signed erosion-rate raster in the result.

Returns
-------
LandscapeStepResult
    Updated surface and diagnostics.

Authors
-------
B. Gailleton

```python
def step(self, compute_erosion_rate: bool = False) -> lsdtt3.LandscapeStepResult
```

#### `time_years`

(self) -> float

---

## `LandscapeStepResult`

Diagnostics and surface produced by one or more landscape-evolution steps.

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
    Diagnostics and surface produced by one or more landscape-evolution steps.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


### Members

#### `delta_max_elevation`

(self) -> float

#### `delta_mean_elevation`

(self) -> float

#### `elevation`

(self) -> lsdtt3.Raster

#### `erosion_rate`

(self) -> object

#### `iteration`

(self) -> int

#### `max_erosion_rate`

(self) -> float

#### `max_local_elevation_change`

(self) -> float

#### `mean_erosion_rate`

(self) -> float

#### `p999_local_elevation_change`

(self) -> float

#### `time_years`

(self) -> float

---

## `hillslope_steady_state`

hillslope_steady_state(dem: lsdtt3.Raster, channels: lsdtt3.Raster, model: str = 'roering', diffusivity: float = 0.0010000000474974513, uplift_rate: float = 9.999999747378752e-05, critical_slope: float = 0.800000011920929, nonlinear: bool = False, boundary_conditions: str = 'oooo', tolerance: float = 9.999999747378752e-05, max_iterations: int = 10000, sor_omega: float = 0.0, n_threads: int = 0) -> lsdtt3.HillslopeSteadyStateResult

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
    hillslope_steady_state(dem: lsdtt3.Raster, channels: lsdtt3.Raster, model: str = 'roering', diffusivity: float = 0.0010000000474974513, uplift_rate: float = 9.999999747378752e-05, critical_slope: float = 0.800000011920929, nonlinear: bool = False, boundary_conditions: str = 'oooo', tolerance: float = 9.999999747378752e-05, max_iterations: int = 10000, sor_omega: float = 0.0, n_threads: int = 0) -> lsdtt3.HillslopeSteadyStateResult

Solve a channel-bounded linear or nonlinear steady hillslope surface.

Parameters
----------
dem : Raster
    Initial surface; elevations at channel cells are fixed Dirichlet values.
channels : Raster
    Same-grid binary raster: one fixes a channel cell and zero solves it as hillslope.
model : {"auto", "roering", "cidre"}, optional
    Newton with safe SOR fallback, direct SOR/Picard, or pseudo-transient routed
    Cidre solution.
diffusivity : float, optional
    Hillslope diffusivity in m² yr⁻¹.
uplift_rate : float, optional
    Uniform non-channel uplift rate in m yr⁻¹.
critical_slope : float, optional
    Critical gradient used when ``nonlinear=True``.
nonlinear : bool, optional
    Use the Roering et al. (1999) nonlinear flux law.
boundary_conditions : str, optional
    Four NESW open or paired-periodic codes.
tolerance, max_iterations, sor_omega, n_threads
    Solver convergence and execution controls.

Returns
-------
HillslopeSteadyStateResult
    Surface and convergence diagnostics.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


```python
def hillslope_steady_state(dem: lsdtt3.Raster, channels: lsdtt3.Raster, model: str = 'roering', diffusivity: float = 0.0010000000474974513, uplift_rate: float = 9.999999747378752e-05, critical_slope: float = 0.800000011920929, nonlinear: bool = False, boundary_conditions: str = 'oooo', tolerance: float = 9.999999747378752e-05, max_iterations: int = 10000, sor_omega: float = 0.0, n_threads: int = 0) -> lsdtt3.HillslopeSteadyStateResult
```

---

## `hillslope_transport_adi`

hillslope_transport_adi(dem: lsdtt3.Raster, diffusivity: float, duration_years: float, boundary_conditions: str = 'oooo') -> lsdtt3.Raster

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
    hillslope_transport_adi(dem: lsdtt3.Raster, diffusivity: float, duration_years: float, boundary_conditions: str = 'oooo') -> lsdtt3.Raster

Advance unconstrained linear hillslope diffusion with Crank-Nicolson ADI sweeps.

Parameters
----------
dem : Raster
    Current elevation surface.
diffusivity : float
    Diffusivity in m² yr⁻¹.
duration_years : float
    Integration duration in years.
boundary_conditions : str, optional
    Four NESW open, no-flux, or paired-periodic codes.

Returns
-------
Raster
    New elevation surface.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


```python
def hillslope_transport_adi(dem: lsdtt3.Raster, diffusivity: float, duration_years: float, boundary_conditions: str = 'oooo') -> lsdtt3.Raster
```

---

## `hillslope_transport_cellular`

hillslope_transport_cellular(dem: lsdtt3.Raster, channels: lsdtt3.Raster, duration_years: float, diffusivity: float = 0.001, uplift_rate: float = 0.0, critical_slope: float = 0.8, r_match: float = 0.85, cfl_safety: float = 0.8, nonlocal_fraction: float = 0.0, nonlocal_activation_ratio: float = 0.8, nonlocal_deposition_probability: float = 0.5, nonlocal_max_hops: int = 4, boundary_conditions: str = 'oooo', n_threads: int = 0) -> lsdtt3.HillslopeCellularResult

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
    hillslope_transport_cellular(dem: lsdtt3.Raster, channels: lsdtt3.Raster, duration_years: float, diffusivity: float = 0.001, uplift_rate: float = 0.0, critical_slope: float = 0.8, r_match: float = 0.85, cfl_safety: float = 0.8, nonlocal_fraction: float = 0.0, nonlocal_activation_ratio: float = 0.8, nonlocal_deposition_probability: float = 0.5, nonlocal_max_hops: int = 4, boundary_conditions: str = 'oooo', n_threads: int = 0) -> lsdtt3.HillslopeCellularResult

Advance conservative face-transfer cellular hillslope transport.

Parameters
----------
dem, channels : Raster
    Current elevation and same-grid binary fixed-channel mask.
duration_years, diffusivity, uplift_rate, critical_slope : float
    Transport duration and physical flux-law parameters.
r_match, cfl_safety, nonlocal_fraction : float, optional
    Cellular stability and optional nonlocal-transfer controls.

Returns
-------
HillslopeCellularResult
    Updated elevation and conservative transfer diagnostics.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


```python
def hillslope_transport_cellular(dem: lsdtt3.Raster, channels: lsdtt3.Raster, duration_years: float, diffusivity: float = 0.001, uplift_rate: float = 0.0, critical_slope: float = 0.8, r_match: float = 0.85, cfl_safety: float = 0.8, nonlocal_fraction: float = 0.0, nonlocal_activation_ratio: float = 0.8, nonlocal_deposition_probability: float = 0.5, nonlocal_max_hops: int = 4, boundary_conditions: str = 'oooo', n_threads: int = 0) -> lsdtt3.HillslopeCellularResult
```

---

## `hillslope_transport_explicit`

hillslope_transport_explicit(dem: lsdtt3.Raster, channels: lsdtt3.Raster, diffusivity: float, duration_years: float, nonlinear: bool = False, critical_slope: float = 0.800000011920929, boundary_conditions: str = 'oooo', n_threads: int = 0) -> lsdtt3.Raster

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
    hillslope_transport_explicit(dem: lsdtt3.Raster, channels: lsdtt3.Raster, diffusivity: float, duration_years: float, nonlinear: bool = False, critical_slope: float = 0.800000011920929, boundary_conditions: str = 'oooo', n_threads: int = 0) -> lsdtt3.Raster

Advance linear or Roering nonlinear hillslope transport with explicit sub-cycling.

Parameters
----------
dem : Raster
    Current elevation surface.
channels : Raster
    Same-grid binary fixed-channel mask.
diffusivity : float
    Diffusivity in m² yr⁻¹.
duration_years : float
    Integration duration in years.
nonlinear : bool, optional
    Use nonlinear Roering transport.
critical_slope : float, optional
    Critical gradient for nonlinear transport.
boundary_conditions : str, optional
    Four NESW open, no-flux, or paired-periodic codes.
n_threads : int, optional
    OpenMP thread count; zero uses the runtime default.

Returns
-------
Raster
    New elevation surface; channels and nodata are preserved.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


```python
def hillslope_transport_explicit(dem: lsdtt3.Raster, channels: lsdtt3.Raster, diffusivity: float, duration_years: float, nonlinear: bool = False, critical_slope: float = 0.800000011920929, boundary_conditions: str = 'oooo', n_threads: int = 0) -> lsdtt3.Raster
```

---

## `hillslope_transport_implicit`

hillslope_transport_implicit(dem: lsdtt3.Raster, channels: lsdtt3.Raster, duration_years: float, diffusivity: float = 0.0010000000474974513, uplift_rate: float = 0.0, critical_slope: float = 0.800000011920929, nonlinear: bool = True, boundary_conditions: str = 'oooo', picard_tolerance: float = 9.999999747378752e-05, picard_max_iterations: int = 50, slope_ratio_cap: float = 0.9900000095367432, n_threads: int = 0) -> lsdtt3.HillslopeTransportResult

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
    hillslope_transport_implicit(dem: lsdtt3.Raster, channels: lsdtt3.Raster, duration_years: float, diffusivity: float = 0.0010000000474974513, uplift_rate: float = 0.0, critical_slope: float = 0.800000011920929, nonlinear: bool = True, boundary_conditions: str = 'oooo', picard_tolerance: float = 9.999999747378752e-05, picard_max_iterations: int = 50, slope_ratio_cap: float = 0.9900000095367432, n_threads: int = 0) -> lsdtt3.HillslopeTransportResult

Advance hillslope transport with implicit backward-Euler integration.

Parameters
----------
dem, channels : Raster
    Current surface and same-grid binary fixed-channel mask.
duration_years : float
    Integration duration in years.
diffusivity, uplift_rate, critical_slope, nonlinear
    Roering transport-law parameters in m² yr⁻¹, m yr⁻¹, and m/m.
boundary_conditions : str, optional
    Four NESW open, no-flux, or paired-periodic codes.
picard_tolerance, picard_max_iterations, slope_ratio_cap, n_threads
    Nonlinear-solver controls.

Returns
-------
HillslopeTransportResult
    New elevation surface and adaptive/Picard diagnostics.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


```python
def hillslope_transport_implicit(dem: lsdtt3.Raster, channels: lsdtt3.Raster, duration_years: float, diffusivity: float = 0.0010000000474974513, uplift_rate: float = 0.0, critical_slope: float = 0.800000011920929, nonlinear: bool = True, boundary_conditions: str = 'oooo', picard_tolerance: float = 9.999999747378752e-05, picard_max_iterations: int = 50, slope_ratio_cap: float = 0.9900000095367432, n_threads: int = 0) -> lsdtt3.HillslopeTransportResult
```
