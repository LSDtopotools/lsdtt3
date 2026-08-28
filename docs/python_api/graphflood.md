# GraphFlood API

## Hydraulic surface conditioning

::: lsdtt3.fill_hydraulic_surface

::: lsdtt3.HydraulicSurfaceFill

## Stationary hydraulic routing

::: lsdtt3.run_graphflood

::: lsdtt3.HydraulicResult

`node_state` is explicit: valid active nodes use `1`, selected outlet nodes
use `3`, and excluded/nodata nodes use `0`. Open boundaries alone do not make a
cell drain; an outlet selection is required.
