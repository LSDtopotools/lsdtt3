# GraphFlood API

## Hydraulic surface conditioning

::: lsdtt3._api_docs.fill_hydraulic_surface

::: lsdtt3._api_docs.HydraulicSurfaceFill
    options:
      members: false

## Stationary hydraulic routing

::: lsdtt3._api_docs.run_graphflood

::: lsdtt3._api_docs.HydraulicResult
    options:
      members: false

`node_state` is explicit: valid active nodes use `1`, selected outlet nodes
use `3`, and excluded/nodata nodes use `0`. Open boundaries alone do not make a
cell drain; an outlet selection is required.
