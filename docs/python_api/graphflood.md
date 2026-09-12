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

## Optimised stationary MFD routing

::: lsdtt3._api_docs.run_graphflood_optimised

::: lsdtt3._api_docs.run_graphflood_optimised_v2

::: lsdtt3._api_docs.run_graphflood_multigrid

::: lsdtt3._api_docs.run_graphflood_induced

The canonical convergence keywords are
`median_dhw_convergence_tolerance_metres` and
`median_dhw_relative_convergence_tolerance`. The former `tol_median_abs_dhw`
and `tol_median_rel_dhw` spellings remain accepted with `DeprecationWarning`
during the compatibility window. Multigrid and induced routing normalize these
options before invoking built-in or custom backends; supplying old and new names
together is an error.

::: lsdtt3._api_docs.GraphFloodMultigridResult
    options:
      members: false

::: lsdtt3._api_docs.GraphFloodInducedResult
    options:
      members: false

::: lsdtt3._api_docs.GraphFloodOptimisedResult
    options:
      members: false

::: lsdtt3._api_docs.GraphFloodOptimisedPerformance
    options:
      members: false

`node_state` is explicit: valid active nodes use `1`, selected outlet nodes
use `3`, and excluded/nodata nodes use `0`. Open boundaries alone do not make a
cell drain; an outlet selection is required.
