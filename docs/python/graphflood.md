# GraphFlood hydraulic routing

`run_graphflood` performs stationary hydraulic routing on a bed DEM. It uses a
separate node-state raster: `0` is nodata, `1` is an ordinary active node, and
`3` is a permitted outlet. Only state `3` can export water from the model.
Boundary conditions specify neighbour connectivity only; they never create an
outlet by themselves.

```python
import numpy as np
import lsdtt3

bed = lsdtt3.Raster.read("dem.tif")
state = np.ones(bed.shape, dtype=np.float32)
state[-1, :] = 3.0  # selected outlets on the open south edge
node_state = lsdtt3.Raster.from_numpy(state, bed.mask_numpy(), bed.metadata)
rainfall_rate = lsdtt3.Raster.from_numpy(
    np.full(bed.shape, 1.0e-6, dtype=np.float32), bed.mask_numpy(), bed.metadata
)
friction_coefficient = lsdtt3.Raster.from_numpy(
    np.full(bed.shape, 1.0 / 0.03, dtype=np.float32), bed.mask_numpy(), bed.metadata
)

result = lsdtt3.run_graphflood(
    bed,
    node_state,
    rainfall_rate,
    friction_coefficient,
    method="mfd",
    boundary_conditions="nnon",
    condition_bed="fill",
)

result.water_depth.write("water_depth.tif")
result.total_discharge.write("discharge.tif")
```

## Separate optimised MFD solver

`run_graphflood_optimised` leaves `run_graphflood` unchanged. It uses a direct
eight-neighbour MFD traversal, persistent compact work arrays, a validated
topology cache, and cumulative stage timings.

```python
fast = lsdtt3.run_graphflood_optimised(
    bed,
    node_state,
    rainfall_rate,
    friction_coefficient,
    boundary_conditions="nnon",
)

print(fast.performance.routing_seconds)
print(fast.performance.topology_rebuilds, fast.performance.topology_reuses)
```

## Incremental v2 solver

`run_graphflood_optimised_v2` retains the v1 interface and adds a radix-heap
Priority-Flood, a bounds-free interior MFD loop, specialised Manning kernels,
and dirty-tile topology checks. A local topology failure still triggers a full
flood, so tile reuse does not bypass the drainage invariant.

The radix heap can resolve equal-elevation cells in a different deterministic
order from v1's binary heap. Depression fills therefore satisfy the same
drainage and minimum-slope constraints but are not guaranteed bitwise-identical
where several equal-cost spill paths exist.

```python
fast_v2 = lsdtt3.run_graphflood_optimised_v2(
    bed, node_state, rainfall_rate, friction_coefficient,
    boundary_conditions="nnon", dirty_tile_size=64,
)
```

## Coarse-to-fine and induced-domain drivers

The multigrid driver can wrap any bundled GraphFlood version, or a compatible
callable. Each factor must divide the raster shape. The built-in v2 path keeps
bed, friction, rainfall, depth transfer, and solver hand-off in C++; other
backends use the transparent `RasterHierarchy` driver. Hydraulic fields must
share the bed grid and be valid on every active cell.
Coarse node states preserve any outlet in their source block, and point inflows
are regrouped while conserving volumetric discharge.

```python
multi = lsdtt3.run_graphflood_multigrid(
    bed, node_state, rainfall_rate, friction_coefficient,
    backend="optimised_v2", factors=(8, 4, 2, 1),
    steps_per_level=(50, 50, 100, 200),
)
```

The induced driver performs one full-domain routing pass, retains cells above
a volumetric `qin` threshold plus a halo, and transfers flux crossing from the
excluded hillslopes onto the retained domain. The retained mask is closed along
all downslope links before solving, preventing flow from leaving and re-entering
the induced domain. Wet domain edges may expand.

```python
induced = lsdtt3.run_graphflood_induced(
    bed, node_state, rainfall_rate, friction_coefficient,
    qin_threshold=1.0, halo_cells=4, max_expansions=2,
)
print(induced.active_fraction)
```

Use `fill_hydraulic_surface` separately when you want to inspect the
hydraulically filled surface before routing. `run_graphflood` returns the bed,
water depth/surface, discharge, slope, and convergence diagnostics in one
`HydraulicResult`.

See [GraphFlood API](../python_api/graphflood.md) for parameters and output
objects.
