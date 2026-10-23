# Landscape And Hillslopes API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

LandscapeModel thresholds use specific catchment area A/dx in metres. Flat relief windows use a full side length in pixels.

## `HillslopeCellularResult`

Conservative cellular hillslope-transport result and mass-transfer diagnostics.

Returned by ``hillslope_transport_cellular()``; not constructible directly.

Attributes
----------
completed : bool
elevation : Raster
    Elevation Raster (m) after the run.
exported_volume_m3 : float
    Net volume (m^3) exported through fixed-elevation (Dirichlet) cells.
max_slope_ratio : float
    Largest S / S_c reached.
message : str
n_face_events : int
n_substeps : int

Authors
-------
B. Gailleton

### Members

#### `completed`

(self) -> bool

#### `elevation`

Elevation Raster (m) after the run.

#### `exported_volume_m3`

Net volume (m^3) exported through fixed-elevation (Dirichlet) cells.

#### `max_slope_ratio`

Largest S / S_c reached.

#### `message`

(self) -> str

#### `n_face_events`

(self) -> int

#### `n_substeps`

(self) -> int

---

## `HillslopeSteadyStateResult`

Steady hillslope-surface solution and convergence diagnostics.

Returned by ``hillslope_steady_state()``; not constructible directly.

Attributes
----------
converged : bool
    True if the solver tolerance was met.
elevation : Raster
    Steady-state elevation Raster (m).
final_residual : float
n_iterations : int
    Solver iterations used.

Authors
-------
B. Gailleton

### Members

#### `converged`

True if the solver tolerance was met.

#### `elevation`

Steady-state elevation Raster (m).

#### `final_residual`

(self) -> float

#### `n_iterations`

Solver iterations used.

---

## `HillslopeTransportResult`

Implicit hillslope-transport result and adaptive-solver diagnostics.

Returned by ``hillslope_transport_implicit()``; not constructible directly.

Attributes
----------
elevation : Raster
    Elevation Raster (m) after the run.
final_dt : float
    Last substep length (years).
limiter_fired : bool
max_slope_ratio : float
    Largest S / S_c in the final surface.
n_rejected_substeps : int
    Substeps rejected by the limiter or divergence.
n_substeps : int
    Accepted substeps.
total_picard_iterations : int

Authors
-------
B. Gailleton

### Members

#### `elevation`

Elevation Raster (m) after the run.

#### `final_dt`

Last substep length (years).

#### `limiter_fired`

(self) -> bool

#### `max_slope_ratio`

Largest S / S_c in the final surface.

#### `n_rejected_substeps`

Substeps rejected by the limiter or divergence.

#### `n_substeps`

Accepted substeps.

#### `total_picard_iterations`

(self) -> int

---

## `LandscapeModel`

Stateful coupled fluvial-incision and hillslope-transport landscape model.

The model owns its evolving elevation raster and reuses native hillslope solver
workspaces between calls to ``step`` or ``run``. Each step applies uplift,
stream-power fluvial incision ``E = K a**m S**n`` (``a = A / dx`` the
specific catchment area in m, ``S`` the gradient in m/m, ``E`` in m/yr) and
hillslope transport. Typical use::

    model = lsdtt3.LandscapeModel(dem, timestep_years=100, uplift_rate=1e-4)
    result = model.run(1000)          # 1000 steps of 100 years
    result.elevation.write("evolved.tif")

Construct it as

    LandscapeModel(dem: Raster, timestep_years: float = 100.0, erodibility: float = 0.001, m_exponent: float = 0.45, n_exponent: float = 1.0, uplift_rate: float = 0.0001, hillslope_diffusivity: float = 0.001, hillslope_model: str = 'cidre', nonlinear_hillslopes: bool = False, critical_hillslope_slope: float = 0.8, incision_threshold_specific_area: float = 0.0, coupling_regime: str = 'automatic', boundary_conditions: str = 'oooo', fluvial_method: str = 'd8', channel_mask: object | None = None, hold_channel_mask_fixed: bool = False, fixed_pixels: object | None = None, n_threads: int = 0, baselevel_fall_rate: float = 0.0, cellular_r_match: float = 0.85, cellular_cfl_safety: float = 0.8, cellular_nonlocal_fraction: float = 0.0, cellular_nonlocal_activation_ratio: float = 0.8, cellular_nonlocal_deposition_probability: float = 0.5, cellular_nonlocal_max_hops: int = 4, erodibility_field: object | None = None, erodibility_m_exponent: object | None = None, erodibility_area_convention: object | None = None)
    LandscapeModel(dem: Raster, timestep_years: float = 100.0, erodibility: float = 0.001, m_exponent: float = 0.45, n_exponent: float = 1.0, uplift_rate: float = 0.0001, hillslope_diffusivity: float = 0.001, hillslope_model: str = 'cidre', nonlinear_hillslopes: bool = False, critical_hillslope_slope: float = 0.8, incision_threshold_area: float = 0.0, coupling_regime: str = 'automatic', boundary_conditions: str = 'oooo', fluvial_method: str = 'd8', channel_mask: object | None = None, hold_channel_mask_fixed: bool = False, fixed_pixels: object | None = None, n_threads: int = 0, baselevel_fall_rate: float = 0.0, cellular_r_match: float = 0.85, cellular_cfl_safety: float = 0.8, cellular_nonlocal_fraction: float = 0.0, cellular_nonlocal_activation_ratio: float = 0.8, cellular_nonlocal_deposition_probability: float = 0.5, cellular_nonlocal_max_hops: int = 4, erodibility_field: object | None = None, erodibility_m_exponent: object | None = None, erodibility_area_convention: object | None = None)

Parameters
----------
dem : Raster
    Initial, depression-conditioned elevation surface (m) in a projected
    metre CRS.
timestep_years : float, default=100.0
    Duration (years) of each call to ``step``.
erodibility : float, default=0.001
    Stream-power erodibility K (> 0), in the specific-area convention above;
    also the reference K for ``erodibility_field``.
m_exponent, n_exponent : float, default=0.45, 1.0
    Stream-power exponents on specific area and slope.
uplift_rate : float, default=0.0001
    Uniform rock uplift rate in m/yr.
hillslope_diffusivity : float, default=0.001
    Hillslope diffusivity D in m^2/yr.
hillslope_model : {"none", "linear", "roering", "roering_implicit", "cidre", "adi", "cellular"}, default="cidre"
    Hillslope process model: ``linear`` explicit diffusion, ``roering``
    explicit nonlinear flux (Roering et al. 1999), ``roering_implicit``
    backward-Euler nonlinear solver, ``cidre`` Carretier et al. (2016),
    ``adi`` implicit linear diffusion, ``cellular`` cellular transport;
    ``none`` disables hillslope transport.
nonlinear_hillslopes : bool, default=False
    Use the Roering nonlinear flux law where the model supports it.
critical_hillslope_slope : float, default=0.8
    Critical gradient S_c (m/m) of the nonlinear flux law.
incision_threshold_specific_area : float, default=0.0
    Channel threshold in specific catchment area A/dx (metres); zero supports fully coupled runs.
incision_threshold_area : float, default=0.0
    Deprecated alias of ``incision_threshold_specific_area`` (emits a
    ``DeprecationWarning``).
coupling_regime : {"automatic", "coupled", "channels_as_boundaries"}, default="automatic"
    ``coupled``: hillslope and fluvial processes act everywhere, so valleys
    and drainage density emerge; ``channels_as_boundaries``: channel pixels
    (above the threshold or in ``channel_mask``) are fixed boundaries of the
    hillslope solve; ``automatic`` chooses ``channels_as_boundaries`` when a
    positive threshold or a channel mask is given, else ``coupled``.
boundary_conditions : str, default='oooo'
    Four NESW codes: ``o`` open, ``n`` no-flux, ``p`` periodic, ``b`` base
    level.
fluvial_method : {"d8", "dinf", "none"}, default="d8"
    Fluvial incision routing method, or ``none`` for hillslope-only evolution.
channel_mask : Raster or None, default=None
    Same-grid 0/1 imposed channel raster.
hold_channel_mask_fixed : bool, default=False
    Add imposed channel pixels to the fixed-elevation set.
fixed_pixels : array-like of int or None, default=None
    Explicit N x 2 ``(row, col)`` fixed-elevation cells. None uses the open
    and base-level edges.
n_threads : int, default=0
    OpenMP thread count; zero uses the runtime default.
baselevel_fall_rate : float, default=0.0
    Downward motion of fixed base-level cells in m/yr.
cellular_r_match : float, default=0.85
    ``cellular`` only: slope ratio r = S / S_c (0 < r < 1) above which the
    nonlinear mobility ``1 / (1 - r**2)`` is continued linearly, avoiding
    the singularity at S_c.
cellular_cfl_safety : float, default=0.8
    ``cellular`` only: stability safety factor on the substep (0 < f <= 1).
cellular_nonlocal_fraction : float, default=0.0
    ``cellular`` only: fraction of near-critical excess flux moved
    non-locally (0 keeps the local solver exactly).
cellular_nonlocal_activation_ratio : float, default=0.8
    ``cellular`` only: slope / S_c ratio (0-1) above which non-local
    transfer starts.
cellular_nonlocal_deposition_probability : float, default=0.5
    ``cellular`` only: probability (0-1] that non-local material deposits at
    each hop.
cellular_nonlocal_max_hops : int, default=4
    ``cellular`` only: largest number of cells non-local material travels.
erodibility_field : Raster or numpy.ndarray or None, default=None
    Per-pixel absolute erodibility K replacing ``erodibility`` in the fluvial
    step: a Raster strictly on the DEM grid or a ``(n_rows, n_cols)`` array
    (NaN only on DEM nodata). Every valid pixel needs a finite K > 0.
erodibility_m_exponent : float or None, default=None
    The m the field was calibrated for; a different m raises ``ValueError``.
erodibility_area_convention : {"map_area", "specific_area"} or None, default=None
    Area convention of the field; ``map_area`` K is converted by dx**m.

Attributes
----------
coupling_regime : str
    Resolved coupling experiment: coupled or channels_as_boundaries.
elevation : Raster
    Current elevation Raster (m).
erodibility_field : object
    The per-pixel erodibility in use (run convention), or None.
hillslope_model : str
    Name of the hillslope model in use.
iteration : int
    Number of steps taken so far.
time_years : float
    Model time (years).

Methods
-------
channel_mask() -> Raster
    Return the imposed or currently threshold-derived channel mask.
run(n_steps: int, compute_erosion_rate: bool = False) -> LandscapeStepResult
    Advance the coupled model repeatedly and return the final step result.
set_erodibility_field(field: object | None) -> None
    Replace the per-pixel erodibility used by subsequent steps.
set_forcing(uplift_rate: float, baselevel_fall_rate: float = 0.0) -> None
    Change the uniform tectonic/base-level forcing used by subsequent steps.
step(compute_erosion_rate: bool = False) -> LandscapeStepResult
    Advance the coupled model by one configured timestep.

Authors
-------
B. Gailleton
S. M. Mudd

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

Current elevation Raster (m).

#### `erodibility_field`

The per-pixel erodibility in use (run convention), or None.

#### `hillslope_model`

Name of the hillslope model in use.

#### `iteration`

Number of steps taken so far.

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

Model time (years).

---

## `LandscapeStepResult`

Diagnostics and surface produced by one or more landscape-evolution steps.

Returned by ``LandscapeModel.run()``, ``LandscapeModel.step()``; not constructible directly.

Attributes
----------
delta_max_elevation : float
    Absolute change of the maximum elevation over the last step (m).
delta_mean_elevation : float
    Absolute change of the mean elevation over the last step (m).
elevation : Raster
    Elevation (m) after the step(s).
erosion_rate : Raster | None
    Per-pixel erosion rate Raster (m/yr, positive = erosion, uplift removed) if requested, else None.
iteration : int
    Total number of steps the model has taken.
max_erosion_rate : float
    Largest erosion rate over the last step (m/yr, before uplift).
max_local_elevation_change : float
    Largest absolute per-pixel elevation change over the last step (m).
mean_erosion_rate : float
    Mean erosion rate over the last step (m/yr, before uplift).
p999_local_elevation_change : float
    99.9th percentile of the absolute per-pixel elevation change (m); a robust steady-state measure.
time_years : float
    Model time (years).

Authors
-------
B. Gailleton

### Members

#### `delta_max_elevation`

Absolute change of the maximum elevation over the last step (m).

#### `delta_mean_elevation`

Absolute change of the mean elevation over the last step (m).

#### `elevation`

Elevation (m) after the step(s).

#### `erosion_rate`

Per-pixel erosion rate Raster (m/yr, positive = erosion, uplift removed) if requested, else None.

#### `iteration`

Total number of steps the model has taken.

#### `max_erosion_rate`

Largest erosion rate over the last step (m/yr, before uplift).

#### `max_local_elevation_change`

Largest absolute per-pixel elevation change over the last step (m).

#### `mean_erosion_rate`

Mean erosion rate over the last step (m/yr, before uplift).

#### `p999_local_elevation_change`

99.9th percentile of the absolute per-pixel elevation change (m); a robust steady-state measure.

#### `time_years`

Model time (years).

---

## `hillslope_steady_state`

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

```python
def hillslope_steady_state(dem: lsdtt3.Raster, channels: lsdtt3.Raster, model: str = 'roering', diffusivity: float = 0.0010000000474974513, uplift_rate: float = 9.999999747378752e-05, critical_slope: float = 0.800000011920929, nonlinear: bool = False, boundary_conditions: str = 'oooo', tolerance: float = 9.999999747378752e-05, max_iterations: int = 10000, sor_omega: float = 0.0, n_threads: int = 0) -> lsdtt3.HillslopeSteadyStateResult
```

---

## `hillslope_transport_adi`

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

```python
def hillslope_transport_adi(dem: lsdtt3.Raster, diffusivity: float, duration_years: float, boundary_conditions: str = 'oooo') -> lsdtt3.Raster
```

---

## `hillslope_transport_cellular`

Advance conservative face-transfer cellular hillslope transport.

Parameters
----------
dem, channels : Raster
    Current elevation and same-grid binary fixed-channel mask.
duration_years, diffusivity, uplift_rate, critical_slope : float
    Transport duration and physical flux-law parameters.
r_match : float, optional
    Slope ratio S / S_c (0-1) above which the nonlinear mobility
    ``1 / (1 - r**2)`` is continued linearly.
cfl_safety : float, optional
    Stability safety factor on the substep (0 < f <= 1).
nonlocal_fraction : float, optional
    Fraction of near-critical excess flux moved non-locally (0 = local only).
nonlocal_activation_ratio : float, optional
    Slope / S_c ratio above which non-local transfer starts.
nonlocal_deposition_probability : float, optional
    Probability (0-1] that non-local material deposits at each hop.
nonlocal_max_hops : int, optional
    Largest number of cells non-local material travels.
boundary_conditions : str, optional
    Four NESW codes ``o``, ``n``, ``p`` or ``b``.
n_threads : int, optional
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
HillslopeCellularResult
    Updated elevation and conservative transfer diagnostics.

Authors
-------
B. Gailleton

```python
def hillslope_transport_cellular(dem: lsdtt3.Raster, channels: lsdtt3.Raster, duration_years: float, diffusivity: float = 0.001, uplift_rate: float = 0.0, critical_slope: float = 0.8, r_match: float = 0.85, cfl_safety: float = 0.8, nonlocal_fraction: float = 0.0, nonlocal_activation_ratio: float = 0.8, nonlocal_deposition_probability: float = 0.5, nonlocal_max_hops: int = 4, boundary_conditions: str = 'oooo', n_threads: int = 0) -> lsdtt3.HillslopeCellularResult
```

---

## `hillslope_transport_explicit`

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

```python
def hillslope_transport_explicit(dem: lsdtt3.Raster, channels: lsdtt3.Raster, diffusivity: float, duration_years: float, nonlinear: bool = False, critical_slope: float = 0.800000011920929, boundary_conditions: str = 'oooo', n_threads: int = 0) -> lsdtt3.Raster
```

---

## `hillslope_transport_implicit`

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

```python
def hillslope_transport_implicit(dem: lsdtt3.Raster, channels: lsdtt3.Raster, duration_years: float, diffusivity: float = 0.0010000000474974513, uplift_rate: float = 0.0, critical_slope: float = 0.800000011920929, nonlinear: bool = True, boundary_conditions: str = 'oooo', picard_tolerance: float = 9.999999747378752e-05, picard_max_iterations: int = 50, slope_ratio_cap: float = 0.9900000095367432, n_threads: int = 0) -> lsdtt3.HillslopeTransportResult
```
