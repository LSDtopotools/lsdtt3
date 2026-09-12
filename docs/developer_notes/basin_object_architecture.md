# Basin Object Architecture

This document records the design boundaries of the basin object hierarchy. User
instructions for the executable are in the
[basin extraction guide](../cli_reference/cli_reference_lsdtt_basin_extraction.md).

## `lsdtt_basin`

`lsdtt_basin` represents the pixels upstream of a validated outlet `flat_index`
from `lsdtt_flow_base`. It holds an `lsdtt_flat_index_slice` into the
Braun-Willett depth-first flow stack rather than copying basin pixels. This
contiguous slicing depends on D8 depth-first ordering and does not apply directly
to breadth-first multiple-flow-direction graphs.

The backend-neutral interface is `lsdtt_basin_base` in `src/core/`; the CPU
implementation is `lsdtt_basin_cpu` in `src/cpu/`. Basin behavior includes
completeness checks against open boundaries and nodata, outline extraction,
raster-value sampling, and pixel counts. Statistical algorithms from
`src/utils/` remain free functions rather than basin members unless they need
private object state.

## `lsdtt_junction_basin`

`lsdtt_junction_basin` associates a basin with an `lsdtt_junction_network`. It
records the first upstream junction and stream order for the supplied outlet.
The outlet itself need not be exactly at a junction.

## `lsdtt_cosmo_basin`

`lsdtt_cosmo_basin` adds cosmogenic information to a basin selected from a
sample location. Samples need not occur at junctions, so this type is separate
from `lsdtt_junction_basin`.

## Basin Management

Basin-manager analysis belongs in `src/cpu/`. A manager may construct ordinary
basins from outlet indices using flow information alone, or junction basins from
a junction network. Selection and filtering operations cover completeness,
pixel-area bounds, stream order, nesting, duplicate outlets, and retaining the
largest basin. Per-basin statistics should compose the manager's raster-value
access with functions from `src/utils/`, avoiding one forwarding member per
statistic.
