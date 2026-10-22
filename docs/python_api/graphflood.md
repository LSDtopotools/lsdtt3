# GraphFlood And Hydraulics API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

Canonical convergence names are used throughout. Legacy tol_* keywords warn and old/new conflicts fail, including accelerator backend options.

## `GraphFloodInducedResult`

Reduced-domain result and prepass diagnostics.

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
Reduced-domain result and prepass diagnostics.

Parameters
----------
result, prepass, active_mask, active_fraction, expansion_results
    Values returned by :func:`run_graphflood_induced`.

Attributes
----------
result : object
    Final reduced-domain backend result.
prepass : object
    One-step full-domain routing result.
active_mask : numpy.ndarray
    Final retained hydraulic domain.
active_fraction : float
    Retained fraction of the original active cells.
expansion_results : tuple
    Result from each domain-expansion solve.
source_discharge, outlet_qin, relative_mass_error : float
    Full forcing, final outlet accumulation, and their relative mismatch.

Authors
-------
B. Gailleton

```python
(result: Any, prepass: Any, active_mask: numpy.ndarray, active_fraction: float, expansion_results: tuple[typing.Any, ...], source_discharge: float, outlet_qin: float, relative_mass_error: float) -> None
```

### Members

#### `__init__`

Initialize self.  See help(type(self)) for accurate signature.

```python
(self, result: Any, prepass: Any, active_mask: numpy.ndarray, active_fraction: float, expansion_results: tuple[typing.Any, ...], source_discharge: float, outlet_qin: float, relative_mass_error: float) -> None
```

---

## `GraphFloodMultigridResult`

Final result plus every coarse-to-fine level.

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
Final result plus every coarse-to-fine level.

Parameters
----------
result, level_results, factors
    Values returned by :func:`run_graphflood_multigrid`.

Attributes
----------
result : object
    Result from the full-resolution level.
level_results : tuple
    Backend result at each requested factor.
factors : tuple of int
    Coarsening factors in execution order.

Authors
-------
B. Gailleton

```python
(result: Any, level_results: tuple[typing.Any, ...], factors: tuple[int, ...]) -> None
```

### Members

#### `__init__`

Initialize self.  See help(type(self)) for accurate signature.

```python
(self, result: Any, level_results: tuple[typing.Any, ...], factors: tuple[int, ...]) -> None
```

---

## `GraphFloodOptimisedPerformance`

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

#### `boundary_routed_cells`

(self) -> int

#### `convergence_seconds`

(self) -> float

#### `depth_update_seconds`

(self) -> float

#### `dirty_tiles_checked`

(self) -> int

#### `friction_seconds`

(self) -> float

#### `full_topology_checks`

(self) -> int

#### `interior_routed_cells`

(self) -> int

#### `persistent_working_bytes`

(self) -> int

#### `projected_cells`

(self) -> int

#### `projection_seconds`

(self) -> float

#### `radix_redistributions`

(self) -> int

#### `routed_cells`

(self) -> int

#### `routing_seconds`

(self) -> float

#### `scratch_reallocations`

(self) -> int

#### `topology_check_seconds`

(self) -> float

#### `topology_rebuilds`

(self) -> int

#### `topology_reuses`

(self) -> int

---

## `GraphFloodOptimisedResult`

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

#### `bed`

(self) -> lsdtt3.Raster

#### `converged`

(self) -> bool

#### `dhw`

(self) -> lsdtt3.Raster

#### `hydraulic_slope`

(self) -> lsdtt3.Raster

#### `performance`

(self) -> lsdtt3.GraphFloodOptimisedPerformance

#### `projection_volume`

(self) -> float

#### `qin`

(self) -> lsdtt3.Raster

#### `qout`

(self) -> lsdtt3.Raster

#### `specific_discharge`

(self) -> lsdtt3.Raster

#### `step_count`

(self) -> int

#### `total_discharge`

(self) -> lsdtt3.Raster

#### `total_water_volume`

(self) -> float

#### `water_depth`

(self) -> lsdtt3.Raster

#### `water_surface`

(self) -> lsdtt3.Raster

---

## `HydraulicResult`

Result of stationary GraphFlood-style hydraulic routing.

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
Result of stationary GraphFlood-style hydraulic routing.

All raster properties use the routed hydraulic domain. ``projection_volume``
is the explicitly non-conservative volume added while priority-flood projects
the evolving ``z + h`` surface.

Authors
-------
B. Gailleton

### Members

#### `bed`

(self) -> lsdtt3.Raster

#### `converged`

(self) -> bool

#### `current_time`

(self) -> float

#### `dhw`

(self) -> lsdtt3.Raster

#### `hydraulic_slope`

(self) -> lsdtt3.Raster

#### `projection_volume`

(self) -> float

#### `qin`

(self) -> lsdtt3.Raster

#### `qout`

(self) -> lsdtt3.Raster

#### `specific_discharge`

(self) -> lsdtt3.Raster

#### `step_count`

(self) -> int

#### `total_discharge`

(self) -> lsdtt3.Raster

#### `total_water_volume`

(self) -> float

#### `water_depth`

(self) -> lsdtt3.Raster

#### `water_surface`

(self) -> lsdtt3.Raster

---

## `HydraulicSurfaceFill`

Priority-flood projection of a hydraulic surface towards selected outlets.

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
Priority-flood projection of a hydraulic surface towards selected outlets.

``added_depth`` and ``added_volume`` quantify the explicitly non-conservative
water required to remove local minima in the supplied surface.

Authors
-------
B. Gailleton

### Members

#### `added_depth`

(self) -> lsdtt3.Raster

#### `added_volume`

(self) -> float

#### `filled_surface`

(self) -> lsdtt3.Raster

---

## `fill_hydraulic_surface`

fill_hydraulic_surface(surface: lsdtt3.Raster, node_state: lsdtt3.Raster, boundary_conditions: str = 'nnnn', slope_epsilon: float = 9.99999993922529e-09, check_valid_outlet: bool = True, priority_flood_backend: str = 'binary') -> lsdtt3.HydraulicSurfaceFill

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
fill_hydraulic_surface(surface: lsdtt3.Raster, node_state: lsdtt3.Raster, boundary_conditions: str = 'nnnn', slope_epsilon: float = 9.99999993922529e-09, check_valid_outlet: bool = True, priority_flood_backend: str = 'binary') -> lsdtt3.HydraulicSurfaceFill

Priority-flood a hydraulic ``z + h`` surface towards selected outlets.

The node-state raster uses 0 for excluded cells, 1 for active cells, and 3
for outlets. Unlike terrain pit filling, excluded/nodata cells are walls, not
drains. Code-3 cells are prescribed sinks wherever they occur.

Authors
-------
B. Gailleton

```python
def fill_hydraulic_surface(surface: lsdtt3.Raster, node_state: lsdtt3.Raster, boundary_conditions: str = 'nnnn', slope_epsilon: float = 9.99999993922529e-09, check_valid_outlet: bool = True, priority_flood_backend: str = 'binary') -> lsdtt3.HydraulicSurfaceFill
```

---

## `run_graphflood`

run_graphflood(bed: lsdtt3.Raster, node_state: lsdtt3.Raster, rainfall_rate: lsdtt3.Raster, friction_coefficient: lsdtt3.Raster, method: str = 'mfd', boundary_conditions: str = 'oooo', mfd_exponent: float = 1.0, condition_bed: str = 'none', fill_slope_epsilon: float = 9.99999993922529e-09, check_valid_outlet: bool = True, friction_depth_exponent: float = 1.6666666666666667, friction_slope_exponent: float = 0.5, dt: float = 0.001, min_slope: float = 0.0001, point_inflows: object | None = None, max_steps: int = 100, run_until_converged: bool = True, convergence_check_every: int = 10, median_dhw_convergence_tolerance_metres: float = 1e-06, median_dhw_relative_convergence_tolerance: float = 0.0001, n_threads: int = 0, initial_water_depth: object | None = None) -> lsdtt3.HydraulicResult

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
run_graphflood(bed: lsdtt3.Raster, node_state: lsdtt3.Raster, rainfall_rate: lsdtt3.Raster, friction_coefficient: lsdtt3.Raster, method: str = 'mfd', boundary_conditions: str = 'oooo', mfd_exponent: float = 1.0, condition_bed: str = 'none', fill_slope_epsilon: float = 9.99999993922529e-09, check_valid_outlet: bool = True, friction_depth_exponent: float = 1.6666666666666667, friction_slope_exponent: float = 0.5, dt: float = 0.001, min_slope: float = 0.0001, point_inflows: object | None = None, max_steps: int = 100, run_until_converged: bool = True, convergence_check_every: int = 10, median_dhw_convergence_tolerance_metres: float = 1e-06, median_dhw_relative_convergence_tolerance: float = 0.0001, n_threads: int = 0, initial_water_depth: object | None = None) -> lsdtt3.HydraulicResult

Run stationary graph-based hydraulic routing on an evolving ``z + h`` surface.

``node_state`` is a raster on the bed grid: 0 means excluded/impermeable, 1
means active, and 3 is a selected outlet/sink. Open edges do not drain unless
a cell is marked 3. With ``check_valid_outlet=True``, every connected active
component must contain a code-3 cell.
Closed edges use ``n`` and periodic edges use paired ``p`` codes. By default,
every connected active component must contain an outlet.

``rainfall_rate`` and ``friction_coefficient`` are required same-grid rasters.
The friction law is ``q / width = coefficient * h**a * S**b``. Point inflows
are an N by 3 array of ``(row, col, discharge_m3_s)``.

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `tol_median_abs_dhw` -> `median_dhw_convergence_tolerance_metres`; earliest removal: following minor.
- `tol_median_rel_dhw` -> `median_dhw_relative_convergence_tolerance`; earliest removal: following minor.

```python
def run_graphflood(bed: lsdtt3.Raster, node_state: lsdtt3.Raster, rainfall_rate: lsdtt3.Raster, friction_coefficient: lsdtt3.Raster, method: str = 'mfd', boundary_conditions: str = 'oooo', mfd_exponent: float = 1.0, condition_bed: str = 'none', fill_slope_epsilon: float = 9.99999993922529e-09, check_valid_outlet: bool = True, friction_depth_exponent: float = 1.6666666666666667, friction_slope_exponent: float = 0.5, dt: float = 0.001, min_slope: float = 0.0001, point_inflows: object | None = None, max_steps: int = 100, run_until_converged: bool = True, convergence_check_every: int = 10, median_dhw_convergence_tolerance_metres: float = 1e-06, median_dhw_relative_convergence_tolerance: float = 0.0001, n_threads: int = 0, initial_water_depth: object | None = None) -> lsdtt3.HydraulicResult
```

---

## `run_graphflood_induced`

Run GraphFlood on a discharge-selected domain and floodplain halo.

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
Run GraphFlood on a discharge-selected domain and floodplain halo.

Parameters
----------
bed, node_state, rainfall_rate, friction_coefficient : Raster
    Full-domain hydraulic inputs.
qin_threshold : float
    Minimum full-domain volumetric ``Qin`` retained as channel domain.
backend : str or callable, default="optimised_v2"
    Backend used by the prepass and reduced solves.
halo_cells : int, default=3
    Initial dilation and wet-edge expansion distance.
edge_depth_threshold : float, default=1e-4
    Depth that activates expansion at a retained-domain edge.
max_expansions : int, default=2
    Maximum number of domain expansions after the first solve.
point_inflows : array-like, optional
    Rows of ``(row, col, discharge_m3_s)``.
**backend_options
    Options forwarded to the backend.

Returns
-------
GraphFloodInducedResult
    Reduced result, prepass, active mask, and expansion history.

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `tol_median_abs_dhw` -> `median_dhw_convergence_tolerance_metres`; earliest removal: following minor.
- `tol_median_rel_dhw` -> `median_dhw_relative_convergence_tolerance`; earliest removal: following minor.

```python
(bed: Any, node_state: Any, rainfall_rate: Any, friction_coefficient: Any, qin_threshold: float, *, backend: str | Callable[..., Any] = 'optimised_v2', halo_cells: int = 3, edge_depth_threshold: float = 0.0001, max_expansions: int = 2, point_inflows: Any = None, median_dhw_convergence_tolerance_metres: float | None = None, median_dhw_relative_convergence_tolerance: float | None = None, **backend_options: Any) -> lsdtt3._graphflood_accelerators.GraphFloodInducedResult
```

---

## `run_graphflood_multigrid`

Run a compatible GraphFlood backend from coarse to fine grids.

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
Run a compatible GraphFlood backend from coarse to fine grids.

Parameters
----------
bed, node_state, rainfall_rate, friction_coefficient : Raster
    Full-resolution hydraulic inputs.
backend : str or callable, default="optimised_v2"
    Bundled backend name or compatible GraphFlood function.
factors : sequence of int, default=(8, 4, 2, 1)
    Strictly descending integer refinement factors ending at one.
steps_per_level : int or sequence of int, default=100
    Solver iterations allocated to each level.
point_inflows : array-like, optional
    Rows of ``(row, col, discharge_m3_s)``.
initial_water_depth : Raster, optional
    Full-resolution initial depth.
**backend_options
    Options forwarded to every backend call.

Returns
-------
GraphFloodMultigridResult
    Final result and per-level results.

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `tol_median_abs_dhw` -> `median_dhw_convergence_tolerance_metres`; earliest removal: following minor.
- `tol_median_rel_dhw` -> `median_dhw_relative_convergence_tolerance`; earliest removal: following minor.

```python
(bed: Any, node_state: Any, rainfall_rate: Any, friction_coefficient: Any, *, backend: str | Callable[..., Any] = 'optimised_v2', factors: Sequence[int] = (8, 4, 2, 1), steps_per_level: int | Sequence[int] = 100, point_inflows: Any = None, initial_water_depth: Any = None, median_dhw_convergence_tolerance_metres: float | None = None, median_dhw_relative_convergence_tolerance: float | None = None, **backend_options: Any) -> lsdtt3._graphflood_accelerators.GraphFloodMultigridResult
```

---

## `run_graphflood_optimised`

run_graphflood_optimised(bed: lsdtt3.Raster, node_state: lsdtt3.Raster, rainfall_rate: lsdtt3.Raster, friction_coefficient: lsdtt3.Raster, boundary_conditions: str = 'oooo', mfd_exponent: float = 1.0, condition_bed: str = 'none', fill_slope_epsilon: float = 9.99999993922529e-09, check_valid_outlet: bool = True, friction_depth_exponent: float = 1.6666666666666667, friction_slope_exponent: float = 0.5, dt: float = 0.001, min_slope: float = 0.0001, point_inflows: object | None = None, max_steps: int = 100, run_until_converged: bool = True, convergence_check_every: int = 10, median_dhw_convergence_tolerance_metres: float = 1e-06, median_dhw_relative_convergence_tolerance: float = 0.0001, reuse_topology: bool = True, initial_water_depth: object | None = None) -> lsdtt3.GraphFloodOptimisedResult

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
run_graphflood_optimised(bed: lsdtt3.Raster, node_state: lsdtt3.Raster, rainfall_rate: lsdtt3.Raster, friction_coefficient: lsdtt3.Raster, boundary_conditions: str = 'oooo', mfd_exponent: float = 1.0, condition_bed: str = 'none', fill_slope_epsilon: float = 9.99999993922529e-09, check_valid_outlet: bool = True, friction_depth_exponent: float = 1.6666666666666667, friction_slope_exponent: float = 0.5, dt: float = 0.001, min_slope: float = 0.0001, point_inflows: object | None = None, max_steps: int = 100, run_until_converged: bool = True, convergence_check_every: int = 10, median_dhw_convergence_tolerance_metres: float = 1e-06, median_dhw_relative_convergence_tolerance: float = 0.0001, reuse_topology: bool = True, initial_water_depth: object | None = None) -> lsdtt3.GraphFloodOptimisedResult

Run the separate CPU-optimised stationary MFD GraphFlood implementation.

The reference ``run_graphflood`` path is unchanged. This implementation uses
direct raster MFD routing and validates a cached topology before reuse.

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `tol_median_abs_dhw` -> `median_dhw_convergence_tolerance_metres`; earliest removal: following minor.
- `tol_median_rel_dhw` -> `median_dhw_relative_convergence_tolerance`; earliest removal: following minor.

```python
def run_graphflood_optimised(bed: lsdtt3.Raster, node_state: lsdtt3.Raster, rainfall_rate: lsdtt3.Raster, friction_coefficient: lsdtt3.Raster, boundary_conditions: str = 'oooo', mfd_exponent: float = 1.0, condition_bed: str = 'none', fill_slope_epsilon: float = 9.99999993922529e-09, check_valid_outlet: bool = True, friction_depth_exponent: float = 1.6666666666666667, friction_slope_exponent: float = 0.5, dt: float = 0.001, min_slope: float = 0.0001, point_inflows: object | None = None, max_steps: int = 100, run_until_converged: bool = True, convergence_check_every: int = 10, median_dhw_convergence_tolerance_metres: float = 1e-06, median_dhw_relative_convergence_tolerance: float = 0.0001, reuse_topology: bool = True, initial_water_depth: object | None = None) -> lsdtt3.GraphFloodOptimisedResult
```

---

## `run_graphflood_optimised_v2`

run_graphflood_optimised_v2(bed: lsdtt3.Raster, node_state: lsdtt3.Raster, rainfall_rate: lsdtt3.Raster, friction_coefficient: lsdtt3.Raster, boundary_conditions: str = 'oooo', mfd_exponent: float = 1.0, condition_bed: str = 'none', fill_slope_epsilon: float = 9.99999993922529e-09, check_valid_outlet: bool = True, friction_depth_exponent: float = 1.6666666666666667, friction_slope_exponent: float = 0.5, dt: float = 0.001, min_slope: float = 0.0001, point_inflows: object | None = None, max_steps: int = 100, run_until_converged: bool = True, convergence_check_every: int = 10, median_dhw_convergence_tolerance_metres: float = 1e-06, median_dhw_relative_convergence_tolerance: float = 0.0001, reuse_topology: bool = True, dirty_tile_size: int = 64, full_validation_every: int = 128, initial_water_depth: object | None = None, prescribed_qin: object | None = None) -> lsdtt3.GraphFloodOptimisedResult

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
run_graphflood_optimised_v2(bed: lsdtt3.Raster, node_state: lsdtt3.Raster, rainfall_rate: lsdtt3.Raster, friction_coefficient: lsdtt3.Raster, boundary_conditions: str = 'oooo', mfd_exponent: float = 1.0, condition_bed: str = 'none', fill_slope_epsilon: float = 9.99999993922529e-09, check_valid_outlet: bool = True, friction_depth_exponent: float = 1.6666666666666667, friction_slope_exponent: float = 0.5, dt: float = 0.001, min_slope: float = 0.0001, point_inflows: object | None = None, max_steps: int = 100, run_until_converged: bool = True, convergence_check_every: int = 10, median_dhw_convergence_tolerance_metres: float = 1e-06, median_dhw_relative_convergence_tolerance: float = 0.0001, reuse_topology: bool = True, dirty_tile_size: int = 64, full_validation_every: int = 128, initial_water_depth: object | None = None, prescribed_qin: object | None = None) -> lsdtt3.GraphFloodOptimisedResult

Run GraphFlood with the incremental v2 CPU kernels.

V2 uses a monotone radix heap in Priority-Flood, a bounds-free interior MFD
loop, specialised Manning exponents, and exact dirty-tile topology validation.
It falls back to a full Priority-Flood whenever a checked edge invalidates the
cached order. ``full_validation_every=0`` disables periodic full audits.
``node_state=2`` is a flux-only ghost cell: it routes its fixed
``prescribed_qin`` downslope while its depth and capacity-limited ``qout``
evolve normally. ``prescribed_qin=None`` means zero at ghost cells.

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `tol_median_abs_dhw` -> `median_dhw_convergence_tolerance_metres`; earliest removal: following minor.
- `tol_median_rel_dhw` -> `median_dhw_relative_convergence_tolerance`; earliest removal: following minor.

```python
def run_graphflood_optimised_v2(bed: lsdtt3.Raster, node_state: lsdtt3.Raster, rainfall_rate: lsdtt3.Raster, friction_coefficient: lsdtt3.Raster, boundary_conditions: str = 'oooo', mfd_exponent: float = 1.0, condition_bed: str = 'none', fill_slope_epsilon: float = 9.99999993922529e-09, check_valid_outlet: bool = True, friction_depth_exponent: float = 1.6666666666666667, friction_slope_exponent: float = 0.5, dt: float = 0.001, min_slope: float = 0.0001, point_inflows: object | None = None, max_steps: int = 100, run_until_converged: bool = True, convergence_check_every: int = 10, median_dhw_convergence_tolerance_metres: float = 1e-06, median_dhw_relative_convergence_tolerance: float = 0.0001, reuse_topology: bool = True, dirty_tile_size: int = 64, full_validation_every: int = 128, initial_water_depth: object | None = None, prescribed_qin: object | None = None) -> lsdtt3.GraphFloodOptimisedResult
```
