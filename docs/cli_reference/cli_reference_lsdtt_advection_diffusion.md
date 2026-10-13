# CLI reference: lsdtt-advection-diffusion

`lsdtt-advection-diffusion` is a deliberately small transient model for
controlled numerical experiments. It combines D8 or D-infinity fluvial
incision, linear ADI hillslope diffusion, and uniform uplift.

It intentionally excludes nonlinear hillslopes, steady-state snapping,
adaptive timesteps, forcing histories, particle transport, and convergence-
based stopping. Those features remain in `lsdtt-lem`.

## Usage

```text
lsdtt-advection-diffusion <path> <parameter_file>
lsdtt-advection-diffusion <parameter_file>
lsdtt-advection-diffusion -h
lsdtt-advection-diffusion -v
lsdtt-advection-diffusion -c
```

When two positional arguments are supplied, `<path>` is used as the parameter
file directory and `<parameter_file>` is the parameter filename. With one
argument, the parameter file is read from the current directory.

## Parameter file format

```text
# Lines beginning with # are comments.
key: value
```

Values are parsed according to the parameter type. Boolean values should be
written as `true` or `false`. Paths may be absolute or relative to the process
working directory unless otherwise stated.

## Pipeline overview

The fixed-timestep workflow is:

1. Parse and validate the parameter file.
2. Load an input grid, resume from the highest checkpoint, or generate the
   red-noise initial surface.
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

## Input/output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | empty | Directory containing an optional initial raster. |
| `read_fname` | string | empty | Optional initial raster filename. |
| `write_path` | string | — | Directory for outputs and continuation checkpoints. |
| `write_fname` | string | — | Prefix for output files. |
| `write_output_format` | string | `cog` | Elevation output format. |
| `frequency_of_output_writing` | int | — | Write an intermediate elevation raster every N iterations. |
| `number_of_iterations` | int | — | Final iteration number. |
| `intermediate_hillshade_format` | string | — | Optional intermediate hillshade format. |
| `intermediate_hillshade_write_path` | string | — | Hillshade output directory. |
| `intermediate_drainage_area_format` | string | — | Optional log10 drainage-area diagnostic format. |
| `intermediate_drainage_area_write_path` | string | — | Drainage-area diagnostic directory. |
| `write_intermediate_png_sidecars` | bool | `false` | Write PNG sidecars for raster checkpoints. |

`read_path` and `read_fname` load an existing grid as the initial surface for a
new run. Raster metadata supplies `n_rows`, `n_cols`, and `cell_size_metres`.
They are not a continuation switch.

Continuation is automatic when no input grid is supplied: the executable scans
`write_path` for the highest numeric file matching:

```text
<write_fname>_model_iteration_<N>.tif
```

It resumes from iteration `N + 1`. To extend a run, increase
`number_of_iterations` and rerun the same parameter file. Initial and final
elevation fields are always written; intermediate elevation, hillshade, and
drainage-area outputs are controlled by the output parameters above.

## Grid, surface, and boundary parameters

| Parameter | Type | Description |
|---|---|---|
| `n_rows` | int | Number of rows when no input raster is supplied. |
| `n_cols` | int | Number of columns when no input raster is supplied. |
| `cell_size_metres` | float | Grid spacing in metres. |
| `noise_target_stddev_metres` | float | Standard deviation target for the red-noise initial surface. |
| `random_seed` | int | Seed for reproducible initial noise. |
| `boundary_conditions` | string | Four-character `N,E,S,W` boundary code; `o` is open and `p` is periodic. |
| `baselevel_elevation` | float | Elevation imposed at open baselevel boundaries. |
| `fill_initial_surface` | bool | Fill the initial surface before starting the model. |

The initial surface is red noise. Set `fill_initial_surface: false` only when
initial internal drainage is intentional.

## Fluvial incision parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fluvial_incision_model` | string | — | `d8` or `dinf`. |
| `K` | float | — | Fluvial erodibility coefficient. |
| `m_exponent` | float | — | Drainage-area exponent. |
| `n_exponent` | float | — | Slope exponent. |
| `use_specific_drainage_area` | bool | `false` | Use `A / dx` instead of total area `A`. |
| `drainage_area_moving_average_kernel_size_pixels` | int | `1` | Odd square smoothing-window size. |

`K` has units `m^(1 - 2m) yr^-1` when total drainage area is used. If
`use_specific_drainage_area` is true, it is interpreted as `K_a`, with units
consistent with the specific-area formulation (Bonetti et al., 2020). The
coefficient conversion is:

```text
K_a = K_A dx^m
```

so `K_a` has units `m^(1 - m) yr^-1` (Litwin et al., 2025). When specific area
is enabled, the fluvial calculation uses `A / dx`; otherwise it uses total area
`A`.

After routing, drainage area can be smoothed with an odd-sized square moving
average. The default `1` leaves it unchanged. Windows wrap across sides marked
`p` in `boundary_conditions` and are clipped at open edges.

## Hillslope and uplift parameters

| Parameter | Type | Description |
|---|---|---|
| `hillslope_model` | string | Must be `adi` for this CLI. |
| `D` | float | Linear diffusivity in `m^2 yr^-1`. |
| `uplift_rate` | float | Uniform uplift in `m yr^-1`. |
| `timestep` | float | Fixed timestep in years. |
| `progress_reporting_interval` | int | Iteration interval for progress messages. |

## Threading

| Parameter | Type | Description |
|---|---|---|
| `n_threads` | int | Number of worker threads; `0` uses the implementation default. |
| `progress_reporting_interval` | int | Controls reporting frequency and does not change threading. |

## Example parameter file

```text
write_path: outputs/z
write_fname: example
write_output_format: cog
read_path:
read_fname:

n_rows: 400
n_cols: 500
cell_size_metres: 2
noise_target_stddev_metres: 4
random_seed: 33
boundary_conditions: opop
baselevel_elevation: 0
fill_initial_surface: true

fluvial_incision_model: dinf
hillslope_model: adi
use_specific_drainage_area: false
drainage_area_moving_average_kernel_size_pixels: 9
K: 0.00025
D: 0.025
uplift_rate: 0.001
m_exponent: 0.5
n_exponent: 1

timestep: 40
number_of_iterations: 1500000
frequency_of_output_writing: 20000
progress_reporting_interval: 500

intermediate_hillshade_format: png
intermediate_hillshade_write_path: outputs/hillshade
intermediate_drainage_area_format: png
intermediate_drainage_area_write_path: outputs/drainage_area
write_intermediate_png_sidecars: false
```
