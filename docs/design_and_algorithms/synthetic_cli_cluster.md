# Choosing a Synthetic Landscape Tool

lsdtt3 has three command-line tools that can create topography. They are not
three interfaces to the same operation. Choose between them according to what
should determine the resulting landscape.

## Quick Decision Guide

| You need | Use | Why |
|---|---|---|
| A controlled analytical or stochastic raster | `lsdtt-synthetic-surfaces` | The surface is constructed directly in one step |
| Ground-truth coefficients, labels, or masks | `lsdtt-synthetic-surfaces` | Truth is known and written by construction |
| An optimised drainage topology | `lsdtt-synthetic-networks` | The OCN receiver graph is the experiment |
| A synthetic island or polygon-domain network | `lsdtt-synthetic-networks` | Domain and outlet geometry are first-class inputs |
| Transient landscape evolution | `lsdtt-lem` | Erosion, transport, uplift, and time determine the result |
| A process steady state | `lsdtt-lem` | The selected process laws determine the topography |
| Evolution from a real DEM | `lsdtt-lem` | It is the only tool in this cluster that reads an initial DEM |
| A noisy initial condition that will immediately evolve | `lsdtt-lem` | Its internal generator avoids an unnecessary intermediate file |

The shortest distinction is:

- **Synthetic surfaces:** specified one-shot fixtures, often with an answer key.
- **Synthetic networks:** optimised drainage topology with elevations added.
- **LEM:** process-evolved topography through time or at steady state.

!!! note "Executable rename"
    `lsdtt-synthetic-surfaces` replaces the former
    `lsdtt-raster-synthetic-tests` executable. The old command is not retained;
    update scripts and reinstall lsdtt3 before using the new name.

## Capability Matrix

| Capability | Synthetic surfaces | Synthetic networks | LEM |
|---|---:|---:|---:|
| Read an initial DEM | No | No | Yes |
| Direct analytical surfaces | Yes | No | No |
| Seeded stochastic generation | Surface values | OCN graph and optional random domain | Initial condition only |
| Ground-truth coefficient/label/mask output | Yes | No | No |
| OCN topology optimisation | No | Yes | No |
| Chi-based fluvial steady state | No | Yes | Yes |
| Transient evolution | No | No | Yes |
| Optional landscape equilibration | No | Yes | Yes |
| D8 routing workflows | Fixture only | Yes | Yes |
| D-infinity process routing | No | No | Yes |
| Main product | Raster fixture and truth layers | Optimised graph represented by a DEM and area raster | Evolved or process-steady landscape |

## `lsdtt-synthetic-surfaces`

Use this CLI when the surface itself is the controlled input to another test or
analysis. It writes one or more selected surface families on a common grid:

- uniform, red, and blue noise;
- diamond-square terrain;
- periodic waves and global quadratics with coefficient truth;
- terrace staircases with tread IDs and channel masks;
- alluvial fans with lobe IDs;
- hillslope strips with channel masks.

The periodic and quadratic surfaces can also exercise wavelet quadratic fitting
and write true, fitted, and residual coefficients. Those modes are diagnostics:
they report errors but do not replace assertion-based CTest coverage.

The CLI does not expose every synthetic generator in the C++ library. In
particular, do not expect a general `surface_type` selector, Perlin noise, cones,
or mountain-range generation from this executable.

[Full parameter reference](../cli_reference/cli_reference_lsdtt_synthetic_surfaces.md)

## `lsdtt-synthetic-networks`

Use this CLI when drainage geometry is the quantity being generated. It builds a
rectangular, island, or polygon domain; selects one or more outlets; and changes
a spanning drainage graph through simulated-annealing edge flips. The objective
is an Optimal Channel Network energy, length-weighted by default on the D8
lattice.

After optimisation, the tool integrates chi upstream on the graph to create a
fluvial steady-state DEM. It can optionally reroute and equilibrate that DEM and
add hillslope steady state.

### Routing fidelity is a scientific result

The intended OCN receiver is always lower than its donor, but another neighbour
can be lower still. Routing the output DEM can therefore produce a different
graph. Always inspect `routing_fidelity` in `_ocn_summary.csv` before using the
DEM downstream.

The current CLI does not write an exact receiver raster, direction raster, edge
table, or network vector. Downstream tools must reroute the DEM. Equilibration
makes a DEM consistent with its new routing, but can move divides away from the
original optimised graph; `eq_ocn_agreement` measures that drift.

[Full parameter reference](../cli_reference/cli_reference_lsdtt_synthetic_networks.md)

## `lsdtt-lem`

Use the LEM when landscape form should emerge from geomorphic processes. It can
start from a real DEM or from an internal flat, white-noise, red-noise, or
diamond-square surface. In the latter case the synthetic surface is an initial
condition, not the experiment's final product.

Time-stepping mode combines stream-power incision, hillslope transport, uplift,
and boundary or base-level forcing. Steady-state modes compute fluvial and
optional hillslope equilibria, with an iterated landscape-equilibration path
available for D8 workflows. The LEM also supports D-infinity where documented.

[Full parameter reference](../cli_reference/cli_reference_lsdtt_lem.md)

## Where The Implementations Overlap

| Shared implementation | Used by |
|---|---|
| Uniform, red-noise, and diamond-square generators | synthetic-surfaces and LEM |
| Synthetic dimensions, cell size, and seed conventions | all three; surfaces currently declares them locally and LEM deliberately overrides grid defaults |
| Chi-based fluvial steady-state core | synthetic-networks and LEM |
| Landscape equilibration and hillslope steady-state solvers | synthetic-networks and LEM |
| Raster metadata and GDAL output | all three |

Shared implementation does not imply interchangeable experiments. An OCN and a
LEM can use the same steady-state equation on the same receiver graph, but only
the OCN tool optimises that graph and only the LEM represents process evolution.

## Common Workflows

### Validate a hillslope solver

1. Generate a hillslope strip and channel mask with
   `lsdtt-synthetic-surfaces`.
2. Load the strip and mask in `lsdtt-lem`.
3. Compare the evolved profile with the analytical or expected profile.

### Evolve an OCN-derived landscape

1. Generate a network and DEM with `lsdtt-synthetic-networks`.
2. Inspect routing fidelity and decide whether to use the raw or equilibrated
   DEM.
3. Load that DEM into `lsdtt-lem` for transient evolution.

### Test a classification or fitting algorithm

1. Generate a surface and companion truth layer with
   `lsdtt-synthetic-surfaces`.
2. Run the analysis CLI on the surface.
3. Score its result against the coefficient, tread, lobe, or channel truth.

### Evolve generic stochastic topography

Use `lsdtt-lem` directly with an internal red-noise or diamond-square initial
surface. Generate an external fixture only if another tool needs the unevolved
surface or controls not exposed by the LEM.

## Reproducibility

- Preserve the parameter file and executable version with every result.
- Set `random_seed`; the shared canonical seed is 33.
- Record dimensions and `cell_size_metres`, because scale changes the physical
  interpretation of every surface and process.
- Retain CRS metadata when chaining tools. Generated utility surfaces use the
  default synthetic projected georeference; polygon-domain networks can inherit
  a real polygon CRS and extent.
- For OCN comparisons, hold `n_iter`, cooling controls, domain, outlets, and
  energy weighting fixed. Record both energy and routing fidelity.
- For LEM comparisons, record timestep controls, process laws, forcing,
  boundaries, and whether a loaded or generated initial surface was used.

## Related Documentation

- [Landscape evolution modelling](landscape_evolution_modelling.md)
- [Surface metrics](surface_metrics.md)
- [Flow routing](flow_routing.md)
- [Testing and benchmarking](../developer_notes/testing_and_benchmarking.md)
