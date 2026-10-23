# GraphFlood And Hydraulics API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

Canonical convergence names are used throughout. Legacy tol_* keywords warn and old/new conflicts fail, including accelerator backend options.

## `GraphFloodInducedResult`

Reduced-domain result and prepass diagnostics.

Construct it as

    GraphFloodInducedResult(result: Any, prepass: Any, active_mask: numpy.ndarray, active_fraction: float, expansion_results: tuple[typing.Any, ...], source_discharge: float, outlet_qin: float, relative_mass_error: float)

It is also returned by ``run_graphflood_induced()``.

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

Construct it as

    GraphFloodMultigridResult(result: Any, level_results: tuple[typing.Any, ...], factors: tuple[int, ...])

It is also returned by ``run_graphflood_multigrid()``.

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

Timing (seconds) and work counters of an optimised GraphFlood run.

The ``*_seconds`` attributes are wall-clock time spent in each solver phase;
the other attributes count operations, for profiling only.

Returned by ``GraphFloodOptimisedResult.performance``; not constructible directly.

Attributes
----------
boundary_routed_cells : int
convergence_seconds : float
depth_update_seconds : float
dirty_tiles_checked : int
friction_seconds : float
full_topology_checks : int
interior_routed_cells : int
persistent_working_bytes : int
projected_cells : int
projection_seconds : float
radix_redistributions : int
routed_cells : int
routing_seconds : float
scratch_reallocations : int
topology_check_seconds : float
topology_rebuilds : int
topology_reuses : int

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

Result of ``run_graphflood_optimised`` and ``run_graphflood_optimised_v2``.

Same hydraulic fields as ``HydraulicResult`` plus solver performance counters.

Returned by ``run_graphflood_optimised()``, ``run_graphflood_optimised_v2()``; not constructible directly.

Attributes
----------
bed : Raster
    Bed elevation used for routing (m), after any conditioning.
converged : bool
    True if the convergence tolerances were met.
dhw : Raster
    Depth change over the last step (m).
hydraulic_slope : Raster
    Water-surface slope used in the flow law (m/m).
performance : GraphFloodOptimisedPerformance
    Timing and work counters (GraphFloodOptimisedPerformance).
projection_volume : float
    Water volume (m^3) added by priority-flood projection of the surface (not conserved).
qin : Raster
    Discharge into each cell (m^3/s).
qout : Raster
    Discharge out of each cell from the flow law (m^3/s).
specific_discharge : Raster
    Discharge per unit width (m^2/s).
step_count : int
    Number of time steps taken.
total_discharge : Raster
    Total discharge through each cell (m^3/s).
total_water_volume : float
    Water stored on the grid (m^3).
water_depth : Raster
    Water depth h (m).
water_surface : Raster
    Water surface elevation z + h (m).

Authors
-------
B. Gailleton

### Members

#### `bed`

Bed elevation used for routing (m), after any conditioning.

#### `converged`

True if the convergence tolerances were met.

#### `dhw`

Depth change over the last step (m).

#### `hydraulic_slope`

Water-surface slope used in the flow law (m/m).

#### `performance`

Timing and work counters (GraphFloodOptimisedPerformance).

#### `projection_volume`

Water volume (m^3) added by priority-flood projection of the surface (not conserved).

#### `qin`

Discharge into each cell (m^3/s).

#### `qout`

Discharge out of each cell from the flow law (m^3/s).

#### `specific_discharge`

Discharge per unit width (m^2/s).

#### `step_count`

Number of time steps taken.

#### `total_discharge`

Total discharge through each cell (m^3/s).

#### `total_water_volume`

Water stored on the grid (m^3).

#### `water_depth`

Water depth h (m).

#### `water_surface`

Water surface elevation z + h (m).

---

## `HydraulicResult`

Result of stationary GraphFlood-style hydraulic routing.

All raster properties use the routed hydraulic domain. ``projection_volume``
is the explicitly non-conservative volume added while priority-flood projects
the evolving ``z + h`` surface.

Returned by ``run_graphflood()``; not constructible directly.

Attributes
----------
bed : Raster
    Bed elevation used for routing (m), after any conditioning.
converged : bool
    True if the convergence tolerances were met.
current_time : float
    Model time at the end of the run (s).
dhw : Raster
    Depth change over the last step (m).
hydraulic_slope : Raster
    Water-surface slope used in the flow law (m/m).
projection_volume : float
    Water volume (m^3) added by priority-flood projection of the surface (not conserved).
qin : Raster
    Discharge into each cell (m^3/s).
qout : Raster
    Discharge out of each cell from the flow law (m^3/s).
specific_discharge : Raster
    Discharge per unit width (m^2/s).
step_count : int
    Number of time steps taken.
total_discharge : Raster
    Total discharge through each cell (m^3/s).
total_water_volume : float
    Water stored on the grid (m^3).
water_depth : Raster
    Water depth h (m).
water_surface : Raster
    Water surface elevation z + h (m).

Authors
-------
B. Gailleton

### Members

#### `bed`

Bed elevation used for routing (m), after any conditioning.

#### `converged`

True if the convergence tolerances were met.

#### `current_time`

Model time at the end of the run (s).

#### `dhw`

Depth change over the last step (m).

#### `hydraulic_slope`

Water-surface slope used in the flow law (m/m).

#### `projection_volume`

Water volume (m^3) added by priority-flood projection of the surface (not conserved).

#### `qin`

Discharge into each cell (m^3/s).

#### `qout`

Discharge out of each cell from the flow law (m^3/s).

#### `specific_discharge`

Discharge per unit width (m^2/s).

#### `step_count`

Number of time steps taken.

#### `total_discharge`

Total discharge through each cell (m^3/s).

#### `total_water_volume`

Water stored on the grid (m^3).

#### `water_depth`

Water depth h (m).

#### `water_surface`

Water surface elevation z + h (m).

---

## `HydraulicSurfaceFill`

Priority-flood projection of a hydraulic surface towards selected outlets.

``added_depth`` and ``added_volume`` quantify the explicitly non-conservative
water required to remove local minima in the supplied surface.

Returned by ``fill_hydraulic_surface()``; not constructible directly.

Attributes
----------
added_depth : Raster
    Depth added by filling (m): filled_surface - surface.
added_volume : float
    Volume added by filling (m^3).
filled_surface : Raster
    Filled surface elevation (m).

Authors
-------
B. Gailleton

### Members

#### `added_depth`

Depth added by filling (m): filled_surface - surface.

#### `added_volume`

Volume added by filling (m^3).

#### `filled_surface`

Filled surface elevation (m).

---

## `fill_hydraulic_surface`

Priority-flood a hydraulic ``z + h`` surface towards selected outlets.

The node-state raster uses 0 for excluded cells, 1 for active cells, and 3
for outlets. Unlike terrain pit filling, excluded/nodata cells are walls, not
drains. Code-3 cells are prescribed sinks wherever they occur.

Parameters
----------
surface : Raster
    Water-surface (or bed) elevation (m) to fill.
node_state : Raster
    Same grid: 0 excluded, 1 active, 3 outlet.
boundary_conditions : str, default='nnnn'
    Four NESW codes ``o``, ``n`` or ``p`` (periodic edges wrap).
slope_epsilon : float, default=1e-8
    Small gradient (m/m) imposed across filled flats so they drain.
check_valid_outlet : bool, default=True
    Require an outlet (code 3) in every connected active region.
priority_flood_backend : {"binary", "radix"}, default="binary"
    Priority-queue implementation; both give the same result.

Returns
-------
HydraulicSurfaceFill
    Filled surface, added depth (m) and added volume (m^3).

Authors
-------
B. Gailleton

```python
def fill_hydraulic_surface(surface: lsdtt3.Raster, node_state: lsdtt3.Raster, boundary_conditions: str = 'nnnn', slope_epsilon: float = 9.99999993922529e-09, check_valid_outlet: bool = True, priority_flood_backend: str = 'binary') -> lsdtt3.HydraulicSurfaceFill
```

---

## `run_graphflood`

Run stationary graph-based hydraulic routing on an evolving ``z + h`` surface.

GraphFlood-style solver (Gailleton et al. 2024): rain falls on every active
cell, discharge is routed downslope over the water surface ``z + h`` (bed
plus water depth) and the depth ``h`` is updated each step until inflow and
outflow balance. Water leaves the domain ONLY through cells marked 3 in
``node_state``, whatever the edge codes in ``boundary_conditions``; to
drain across all four raster edges, mark them as outlets::

    edges = lsdtt3.make_outlet_mask(bed, True, True, True, True).to_numpy()
    node_state = lsdtt3.Raster.from_numpy(1 + 2 * edges, bed.mask_numpy(),
                                          bed.metadata)

With ``check_valid_outlet=True`` (the default) every connected group of
active cells must contain a code-3 cell, otherwise ``ValueError`` is raised.

Parameters
----------
bed : Raster
    Bed (terrain) elevation (m) in a projected metre CRS.
node_state : Raster
    Same grid as ``bed``: 0 excluded (wall), 1 active, 3 outlet (water
    leaving here is removed). Nodata pixels are excluded.
rainfall_rate : Raster
    Same grid: rainfall (effective runoff) rate in m/s, for example
    ``10 mm/h = 10 / 1000 / 3600 = 2.78e-6`` m/s.
friction_coefficient : Raster
    Same grid: coefficient ``K`` of the flow law ``q = K * h**a * S**b``
    (``q`` discharge per unit width, m^2/s). For Manning, ``K = 1 / n``
    (n in s m^(-1/3), for example n = 0.033 gives K = 30).
method : {"mfd", "dinf"}, default="mfd"
    Routing over the water surface: multiple-flow-direction or D-infinity.
boundary_conditions : str, default='oooo'
    Four NESW codes ``o``, ``n`` or ``p`` setting the neighbour topology of
    the edges (``p`` wraps to the opposite edge and must be paired). Edges
    do not drain by themselves: see ``node_state``.
mfd_exponent : float, default=1.0
    Exponent on slope used to split MFD discharge among receivers
    (positive).
condition_bed : {"none", "fill"}, default="none"
    ``"fill"`` priority-floods the bed towards the code-3 cells before
    routing.
fill_slope_epsilon : float, default=1e-8
    Small gradient (m/m) imposed across flats when the water surface (or
    bed) is priority-flooded towards the outlets.
check_valid_outlet : bool, default=True
    Require an outlet (code 3) in every connected active region.
friction_depth_exponent : float, default=5/3
    Exponent ``a`` on depth in the flow law (5/3 for Manning).
friction_slope_exponent : float, default=0.5
    Exponent ``b`` on hydraulic slope in the flow law (1/2 for Manning).
dt : float, default=0.001
    Time step in seconds of the depth update.
min_slope : float, default=0.0001
    Floor (m/m) on the hydraulic slope used in the flow law, so flat water
    surfaces still discharge.
point_inflows : array-like of float or None, default=None
    N x 3 array of ``(row, col, discharge_m3_s)``: extra inflow (m^3/s)
    added at given pixels, for example a river entering the domain.
max_steps : int, default=100
    Maximum number of time steps.
run_until_converged : bool, default=True
    Stop early when converged; if False run exactly ``max_steps`` steps.
convergence_check_every : int, default=10
    Check convergence every this many steps (positive).
median_dhw_convergence_tolerance_metres : float, default=1e-6
    Converged when the median absolute depth change per step (m) is at most
    this ...
median_dhw_relative_convergence_tolerance : float, default=0.0001
    ... and the median relative depth change is at most this.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
initial_water_depth : Raster or None, default=None
    Starting water depth (m) on the bed grid; None starts dry.

Returns
-------
HydraulicResult
    Water depth (m), water surface, discharges and convergence state.

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
median_dhw_convergence_tolerance_metres : float or None, default=None
    Stop a solve when the median absolute water-depth change per step
    (m) is at most this value (and the relative test below also holds).
    None uses the backend's own default.
median_dhw_relative_convergence_tolerance : float or None, default=None
    Relative counterpart: the median depth change divided by the depth
    must also be at most this value (dimensionless). None uses the
    backend's own default.
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
median_dhw_convergence_tolerance_metres : float or None, default=None
    Stop a solve when the median absolute water-depth change per step
    (m) is at most this value (and the relative test below also holds).
    None uses the backend's own default.
median_dhw_relative_convergence_tolerance : float or None, default=None
    Relative counterpart: the median depth change divided by the depth
    must also be at most this value (dimensionless). None uses the
    backend's own default.
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

Run the separate CPU-optimised stationary MFD GraphFlood implementation.

The reference ``run_graphflood`` path is unchanged. This implementation uses
direct raster MFD routing and validates a cached topology before reuse. The
physics and inputs are those of ``run_graphflood`` (MFD routing only): water
leaves only through ``node_state`` 3 cells.

Parameters
----------
bed : Raster
    Bed (terrain) elevation (m) in a projected metre CRS.
node_state : Raster
    Same grid as ``bed``: 0 excluded (wall), 1 active, 3 outlet (water
    leaving here is removed). Nodata pixels are excluded.
rainfall_rate : Raster
    Same grid: rainfall (effective runoff) rate in m/s, for example
    ``10 mm/h = 10 / 1000 / 3600 = 2.78e-6`` m/s.
friction_coefficient : Raster
    Same grid: coefficient ``K`` of the flow law ``q = K * h**a * S**b``
    (``q`` discharge per unit width, m^2/s). For Manning, ``K = 1 / n``
    (n in s m^(-1/3), for example n = 0.033 gives K = 30).
boundary_conditions : str, default='oooo'
    Four NESW codes ``o``, ``n`` or ``p`` setting the neighbour topology of
    the edges (``p`` wraps to the opposite edge and must be paired). Edges
    do not drain by themselves: see ``node_state``.
mfd_exponent : float, default=1.0
    Exponent on slope used to split MFD discharge among receivers
    (positive).
condition_bed : {"none", "fill"}, default="none"
    ``"fill"`` priority-floods the bed towards the code-3 cells before
    routing.
fill_slope_epsilon : float, default=1e-8
    Small gradient (m/m) imposed across flats when the water surface (or
    bed) is priority-flooded towards the outlets.
check_valid_outlet : bool, default=True
    Require an outlet (code 3) in every connected active region.
friction_depth_exponent : float, default=5/3
    Exponent ``a`` on depth in the flow law (5/3 for Manning).
friction_slope_exponent : float, default=0.5
    Exponent ``b`` on hydraulic slope in the flow law (1/2 for Manning).
dt : float, default=0.001
    Time step in seconds of the depth update.
min_slope : float, default=0.0001
    Floor (m/m) on the hydraulic slope used in the flow law, so flat water
    surfaces still discharge.
point_inflows : array-like of float or None, default=None
    N x 3 array of ``(row, col, discharge_m3_s)``: extra inflow (m^3/s)
    added at given pixels, for example a river entering the domain.
max_steps : int, default=100
    Maximum number of time steps.
run_until_converged : bool, default=True
    Stop early when converged; if False run exactly ``max_steps`` steps.
convergence_check_every : int, default=10
    Check convergence every this many steps (positive).
median_dhw_convergence_tolerance_metres : float, default=1e-6
    Converged when the median absolute depth change per step (m) is at most
    this ...
median_dhw_relative_convergence_tolerance : float, default=0.0001
    ... and the median relative depth change is at most this.
reuse_topology : bool, default=True
    Reuse the cached routing order between steps when it is still valid.
initial_water_depth : Raster or None, default=None
    Starting water depth (m) on the bed grid; None starts dry.

Returns
-------
GraphFloodOptimisedResult
    Water depth (m), discharges, convergence state and timing counters.

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

Run GraphFlood with the incremental v2 CPU kernels.

V2 uses a monotone radix heap in Priority-Flood, a bounds-free interior MFD
loop, specialised Manning exponents, and exact dirty-tile topology validation.
It falls back to a full Priority-Flood whenever a checked edge invalidates the
cached order. ``full_validation_every=0`` disables periodic full audits.
``node_state=2`` is a flux-only ghost cell: it routes its fixed
``prescribed_qin`` downslope while its depth and capacity-limited ``qout``
evolve normally. ``prescribed_qin=None`` means zero at ghost cells.

Parameters
----------
bed : Raster
    Bed (terrain) elevation (m) in a projected metre CRS.
node_state : Raster
    Same grid as ``bed``: 0 excluded (wall), 1 active, 3 outlet (water
    leaving here is removed). Nodata pixels are excluded.
rainfall_rate : Raster
    Same grid: rainfall (effective runoff) rate in m/s, for example
    ``10 mm/h = 10 / 1000 / 3600 = 2.78e-6`` m/s.
friction_coefficient : Raster
    Same grid: coefficient ``K`` of the flow law ``q = K * h**a * S**b``
    (``q`` discharge per unit width, m^2/s). For Manning, ``K = 1 / n``
    (n in s m^(-1/3), for example n = 0.033 gives K = 30).
boundary_conditions : str, default='oooo'
    Four NESW codes ``o``, ``n`` or ``p`` setting the neighbour topology of
    the edges (``p`` wraps to the opposite edge and must be paired). Edges
    do not drain by themselves: see ``node_state``.
mfd_exponent : float, default=1.0
    Exponent on slope used to split MFD discharge among receivers
    (positive).
condition_bed : {"none", "fill"}, default="none"
    ``"fill"`` priority-floods the bed towards the code-3 cells before
    routing.
fill_slope_epsilon : float, default=1e-8
    Small gradient (m/m) imposed across flats when the water surface (or
    bed) is priority-flooded towards the outlets.
check_valid_outlet : bool, default=True
    Require an outlet (code 3) in every connected active region.
friction_depth_exponent : float, default=5/3
    Exponent ``a`` on depth in the flow law (5/3 for Manning).
friction_slope_exponent : float, default=0.5
    Exponent ``b`` on hydraulic slope in the flow law (1/2 for Manning).
dt : float, default=0.001
    Time step in seconds of the depth update.
min_slope : float, default=0.0001
    Floor (m/m) on the hydraulic slope used in the flow law, so flat water
    surfaces still discharge.
point_inflows : array-like of float or None, default=None
    N x 3 array of ``(row, col, discharge_m3_s)``: extra inflow (m^3/s)
    added at given pixels, for example a river entering the domain.
max_steps : int, default=100
    Maximum number of time steps.
run_until_converged : bool, default=True
    Stop early when converged; if False run exactly ``max_steps`` steps.
convergence_check_every : int, default=10
    Check convergence every this many steps (positive).
median_dhw_convergence_tolerance_metres : float, default=1e-6
    Converged when the median absolute depth change per step (m) is at most
    this ...
median_dhw_relative_convergence_tolerance : float, default=0.0001
    ... and the median relative depth change is at most this.
reuse_topology : bool, default=True
    Reuse the cached routing order between steps when it is still valid.
dirty_tile_size : int, default=64
    Side (pixels) of the tiles whose routing is re-checked after a change.
full_validation_every : int, default=128
    Run a full topology audit every this many steps; 0 disables it.
initial_water_depth : Raster or None, default=None
    Starting water depth (m) on the bed grid; None starts dry.
prescribed_qin : Raster or None, default=None
    Fixed inflow discharge (m^3/s) of ``node_state`` 2 ghost cells; None
    means zero.

Returns
-------
GraphFloodOptimisedResult
    Water depth (m), discharges, convergence state and timing counters.

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
