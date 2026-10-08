# CLI reference: lsdtt-advection-diffusion

`lsdtt-advection-diffusion` is a deliberately small transient model for
controlled numerical experiments. It combines D8 or D-infinity fluvial
incision, linear ADI hillslope diffusion, and uniform uplift.

It intentionally excludes nonlinear hillslopes, steady-state snapping,
adaptive timesteps, forcing histories, particle transport, and convergence-
based stopping. Those features remain in `lsdtt-lem`.

## Workflow

The fixed-timestep workflow is:

1. Parse and validate the parameter file.
2. Load the highest available intermediate surface, or generate the red-noise
   initial surface when no intermediate output exists.
3. Fill the initial surface and apply baselevel boundaries.
4. Enter the fixed-timestep loop.
5. Rebuild flow routing and apply fluvial incision.
6. Apply ADI hillslope diffusion.
7. Apply uniform uplift and restore baselevel boundaries.
8. Write scheduled outputs.
9. Write final elevation and statistics.

The fluvial incision implementation follows Braun and Willett (2013) for D8
and Anand et al. (2020) for D-infinity. The ADI hillslope diffusion follows
Perron et al. (2008).

## Usage

```text
lsdtt-advection-diffusion <path> <parameter_file>
lsdtt-advection-diffusion <parameter_file>
lsdtt-advection-diffusion -h
```

## Core parameters

```text
fluvial_incision_model: d8 | dinf
K: ...
m_exponent: ...
n_exponent: ...
use_specific_drainage_area: false
drainage_area_moving_average_kernel_size_pixels: 1
hillslope_model: adi
D: ...
uplift_rate: ...
```

`K` has units `m^(1 - 2m) yr^-1` when total drainage area is used. Note: if
`use_specific_drainage_area` is true, this parameter is interpreted as `K_a`,
an erodibility with units adjusted to be consistent with the specific drainage
area formulation (Bonetti et al., 2020, PNAS). `K_a` is equivalent to
`K*dx^m` (see Litwin et al., 2025, ESurf), so has units of
`m^(1 - m) yr^-1`. `D` has units `m^2 yr^-1`, and `uplift_rate` has units
`m yr^-1`.

## Input grids and continuation

The optional file-I/O parameters `read_path` and `read_fname` load an existing
grid as the initial surface for a **new run**. The raster metadata supplies
`n_rows`, `n_cols`, and `cell_size_metres` automatically. These parameters are
not a continuation or restart switch.

Continuation is separate and automatic: when no `read_path`/`read_fname` input
grid is supplied, the executable searches `write_path` for the highest
`<write_fname>_model_iteration_<N>.tif` and resumes from that iteration.

When `use_specific_drainage_area` is false, total drainage area `A` is used.
When true, the fluvial calculation uses `A / dx`. The former follows the
nondimensional scaling of Theodratos et al. (2018); the latter follows the
characteristic-scale convention of Litwin et al. (2025).

After flow routing, the drainage-area field used by stream-power incision can
be smoothed with an odd-sized square moving-average kernel. The default value
of `1` leaves the field unchanged. Larger positive odd values average valid
drainage-area cells in the local window before the new incised elevation is
calculated. Windows wrap across any side marked `p` in `boundary_conditions`
(`N,E,S,W` order); other sides are clipped at the DEM edge.

The initial surface is always red noise. `fill_initial_surface` defaults to
true. Set it to false only when initial internal drainage is intentional.

## Output checkpoints

Initial and final elevation fields are always written. If
`frequency_of_output_writing` is positive, elevation is also written every N
iterations. Optional hillshade and log10 drainage-area diagnostics use the
same initial, intermediate, and final checkpoints.

### Resuming a run

There is no parameter switch to enable restart; resuming is automatic. Before
generating a new surface, the executable scans `write_path` for files matching:

```text
<write_fname>_model_iteration_<N>.tif
```

It loads the file with the highest numeric `N` and starts the timestep loop at
`N + 1`, up to `number_of_iterations`. To extend a run, increase
`number_of_iterations` and rerun the same parameter file. If no matching
intermediate TIFF exists, normal red-noise initialisation is used.
