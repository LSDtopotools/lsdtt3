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

Use `fill_hydraulic_surface` separately when you want to inspect the
hydraulically filled surface before routing. `run_graphflood` returns the bed,
water depth/surface, discharge, slope, and convergence diagnostics in one
`HydraulicResult`.

See [GraphFlood API](../python_api/graphflood.md) for parameters and output
objects.
