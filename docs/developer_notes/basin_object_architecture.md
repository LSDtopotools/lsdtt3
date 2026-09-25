# Basin Object Architecture

This document records the design boundaries of the basin object hierarchy. User
instructions for the executable are in the
[basin extraction guide](../cli_reference/cli_reference_lsdtt_basin_extraction.md).

## `lsdtt_basin_cpu`

`lsdtt_basin_cpu` represents the pixels upstream of a validated outlet `flat_index`
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

## `lsdtt_junction_basin_cpu`

`lsdtt_junction_basin_cpu` associates a basin with an `lsdtt_junction_network`. It
records the first upstream junction and stream order for the supplied outlet.
The outlet itself need not be exactly at a junction.

## Cosmogenic sites

There is no `lsdtt_cosmo_basin`. Attaching cosmogenic information to a basin
selected from a sample location is handled by `lsdtt_cosmo_site_cpu`, which
covers everything between "I have a DEM and some parameters" and "I have a
production raster and three shielding rasters", with `lsdtt_cosmo_raster_cpu`
underneath.

The design constraint to know about before touching it is **capabilities, not
a union**. The cosmogenic CLIs are not the same shape — `lsdtt-cosmo-tool` has
no shielding-raster parameters and no production-from-disk path — so a caller
states which capabilities it has, and anything outside that set is
unavailable rather than merely unused. Unifying the CLIs naively would have
added an untested code path to an inversion tool while passing every numeric
check, because the unused prefixes default to NULL. Changing a CLI's
capability set is a feature decision with its own review, never a side effect
of sharing code.

Two related rules: `validate_cosmo_site_config()` does no I/O and is meant to
be called immediately after parsing, so a bad `scaling_scheme` fails before a
DEM is read; and production is re-computable rather than merely readable,
because `lsdtt-cosmo-tool`'s self-consistent-age mode rebuilds production over
one basin's pixels per iterate. The full rationale is in the header comment of
`src/cpu/lsdtt_cosmo_site_cpu.hpp`. The science is in
[`cosmogenic_nuclides.md`](../design_and_algorithms/cosmogenic_nuclides.md).

## Basin Management

Basin-manager analysis belongs in `src/cpu/`. Two managers exist:
`lsdtt_junction_basin_manager_cpu`, which constructs junction basins from a
junction network, and `lsdtt_zero_order_basin_manager_cpu`, for zero-order
(unchannelised) basins. Selection and filtering operations cover completeness,
pixel-area bounds, stream order, nesting, duplicate outlets, and retaining the
largest basin. Per-basin statistics should compose the manager's raster-value
access with functions from `src/utils/`, avoiding one forwarding member per
statistic.

## Other modules in the basin layer

- `lsdtt_basin_outline_cpu` — extracts the basin outline as a polygon.
- `lsdtt_basin_outline_metrics_cpu` — planform shape metrics computed on that
  outline (elongation, circularity, Feret diameters and the rest). The
  derivations are in `docs/tex_technical/basin_shape_metrics.tex`.
- `lsdtt_basin_stats_cpu` — per-basin statistics composed over the manager's
  raster-value access.
- `lsdtt_side_tributary_basins_cpu` — tributary basins entering a trunk channel.
- `lsdtt_outlet_basin_verification_cpu` — checks that a resolved outlet yields
  the basin the caller expected. Part of the outlet-resolution stack
  (`lsdtt_outlet_bounds_discovery`, `lsdtt_outlet_clusters`,
  `lsdtt_outlet_resolution_cpu`) used by `lsdtt-outlet-matching`.
