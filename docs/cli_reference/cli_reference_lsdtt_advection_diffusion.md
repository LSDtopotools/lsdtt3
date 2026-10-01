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
2. Generate the red-noise initial surface.
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
hillslope_model: adi
D: ...
uplift_rate: ...
```

`K` has units `m^(1 - 2m) yr^-1` when total drainage area is used. `D` has
units `m^2 yr^-1`, and `uplift_rate` has units `m yr^-1`.

When `use_specific_drainage_area` is false, total drainage area `A` is used.
When true, the fluvial calculation uses `A / dx`. The former follows the
nondimensional scaling of Theodratos et al. (2018); the latter follows the
characteristic-scale convention of Litwin et al. (2025).

The initial surface is always red noise. `fill_initial_surface` defaults to
true. Set it to false only when initial internal drainage is intentional.

## Output checkpoints

Initial and final elevation fields are always written. If
`frequency_of_output_writing` is positive, elevation is also written every N
iterations. Optional hillshade and log10 drainage-area diagnostics use the
same initial, intermediate, and final checkpoints.
