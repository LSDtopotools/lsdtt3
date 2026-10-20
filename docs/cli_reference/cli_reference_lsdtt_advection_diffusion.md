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
| `write_hillshade_rasters` | bool | `false` | Write hillshade diagnostics at the initial, intermediate and final steps. |
| `intermediate_hillshade_write_path` | string | — | Hillshade output directory. |
| `write_drainage_area_rasters` | bool | `false` | Write log10 drainage-area diagnostics at the initial, intermediate and final steps. |
| `diagnostics_write_output_format` | string | `png` | Diagnostic format: `cog`, `png`, or `both`. |
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
| `multiply_stream_power_term_by_w_over_dx` | bool | `false` | Multiply the stream-power incision term by `w / dx` to account for channel width. |
| `channel_width_w` | float | `1` | Channel width `w` in the same length units as `dx`; used when the width correction is enabled. |
| `drainage_area_moving_average_kernel_size_pixels` | int | `1` | Odd square smoothing-window size. |

The optional `w / dx` multiplier is a flow-width correction proposed to reduce
grid-resolution dependence by Perron et al. (2008), following Howard (1994);
see also Hergarten (2020). It represents a fixed channel width `w` and is
independent of the choice to use total or specific drainage area.

This is distinct from the specific-area formulation discussed by Kargère et al.
(2025): using `A / dx` with `m_exponent: 1` and `n_exponent: 1` is an alternative
that implies channel concavity `theta = m / n = 1`. For `m_exponent: 0.5` and
`n_exponent: 1`, the channel concavity is `theta = 0.5`; the new width option
does not change the area or slope exponents.

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
multiply_stream_power_term_by_w_over_dx: false
channel_width_w: 1
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

write_hillshade_rasters: true
intermediate_hillshade_write_path: outputs/hillshade
write_drainage_area_rasters: true
diagnostics_write_output_format: png
intermediate_drainage_area_write_path: outputs/drainage_area
write_intermediate_png_sidecars: false
```

<!-- BEGIN GENERATED CLI PARAMETERS -->
## Parameters

This section is generated from the executable's `-h json` registry. Do not edit it manually.

**Status: experimental.** This CLI is not yet fully tested.

### File I/O

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `read_path` | `string` | - | `io` | - | Optional directory containing an existing input grid. | This loads a grid for a new run and automatically takes n_rows, n_cols and cell_size_metres from its raster metadata. It is distinct from automatic continuation: continuation instead discovers the highest &lt;write_fname&gt;_model_iteration_N.tif in write_path. |
| `read_fname` | `string` | - | `io` | - | Filename of the existing input grid in read_path. | Use with read_path to start a new run from an existing grid; this is not a continuation switch. |
| `write_path` | `string` | `./outputs` | `io` | - | Directory for outputs. | - |
| `write_fname` | `string` | `advection_diffusion` | `io` | - | Output filename prefix. | - |
| `write_output_format` | `string` | `cog` | `io` | - | Elevation output: cog, gtiff, or png. | - |

### Grid settings

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_rows` | `int` | `128` | `param` | - | Number of grid rows. | - |
| `n_cols` | `int` | `128` | `param` | - | Number of grid columns. | - |
| `cell_size_metres` | `float` | `1` | `param` | - | Grid cell size (m). | - |
| `boundary_conditions` | `string` | `opop` | `param` | - | Boundary conditions in N,E,S,W order. | o=open fixed-baselevel, n=no-flux, p=periodic, b=forced outlet. |
| `n_threads` | `int` | `0` | `io` | - | OpenMP threads; 0 uses all available. | - |

### Initial conditions

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `random_seed` | `int` | `33` | `param` | - | Red-noise random seed. | - |
| `noise_target_stddev_metres` | `float` | `1` | `param` | - | Red-noise standard deviation (m). | - |
| `fill_initial_surface` | `bool` | `true` | `switch` | - | Fill the red-noise surface before routing. | Set false only when initial internal drainage is intentional. |

### Fluvial incision

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fluvial_incision_model` | `string` | `dinf` | `param` | - | Fluvial model: d8 or dinf. | D8 follows the numerical scheme of Braun &amp; Willett (2013, doi:10.1016/j.geomorph.2012.10.008); whereas the dinf follows the numerical scheme of Anand et al. (2020, doi:10.1016/j.envsoft.2020.104804). |
| `K` | `float` | `1e-3` | `param` | - | Erodibility K (m^(1-2m) yr^-1). | Note: if use_specific_drainage_area is true, this parameter is interpreted as K_a, an erodibility with units adjusted to be consistent with the specific drainage area formulation (Bonetti et al., 2020, PNAS). K_a is equivalent to K*dx^m (see Litwin et al., 2025, ESurf), so has units of m^(1-m) yr^-1. |
| `m_exponent` | `float` | `0.45` | `param` | - | Drainage-area exponent m. | - |
| `n_exponent` | `float` | `1` | `param` | - | Slope exponent n. | - |
| `use_specific_drainage_area` | `bool` | `false` | `param` | - | Use A/dx rather than total drainage area A. | When false, normal drainage area is used and the model scales according to the nondimensionalised equations provided by Theodoratos et al., (2018, ESurf). When true, specific drainage area is used and the model scales according to the nondimensionalised equations in Litwin et al. (2025, ESurf). |
| `multiply_stream_power_term_by_w_over_dx` | `bool` | `false` | `switch` | - | Multiply the stream-power incision term by channel_width_w / dx. | This flow-width correction was proposed to reduce grid-resolution dependence by Perron et al. (2008, doi:10.1029/2007JF000977), following Howard (1994, doi:10.1029/94WR00757); see also Hergarten (2020, doi:10.5194/esurf-8-367-2020). The fixed channel-width approach is that of Perron et al. (2008). This correction is separate from use_specific_drainage_area. An alternative is to use specific area A/dx with m=1 and n=1 (Kargere et al., 2025, doi:10.5194/esurf-13-403-2025); that formulation implies channel concavity theta=m/n=1, rather than theta=0.5 for m=0.5 and n=1. |
| `channel_width_w` | `float` | `1` | `param` | - | Channel width w, in the same length units as dx. | Used only when multiply_stream_power_term_by_w_over_dx is true. The incision term is multiplied by w/dx. |
| `drainage_area_moving_average_kernel_size_pixels` | `int` | `1` | `param` | - | Odd-sized square moving-average kernel for drainage area used by stream-power incision; 1 disables smoothing. | Only positive odd values are allowed. |

### Hillslope diffusion

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `hillslope_model` | `string` | `adi` | `param` | - | Hillslope model; only adi is supported. | ADI linear diffusion follows Perron et al. (2008, doi:10.1029/2007JF000977). |
| `D` | `float` | `1e-3` | `param` | - | Hillslope diffusivity D (m^2 yr^-1). | - |

### Uplift

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `uplift_rate` | `float` | `1e-4` | `param` | - | Uniform uplift rate (m yr^-1). | - |
| `baselevel_elevation` | `float` | `0` | `param` | - | Fixed open-boundary elevation (m). | - |

### Timestepping

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `timestep` | `float` | `100` | `param` | - | Fixed timestep (yr). | Adaptive timesteps are not supported. |
| `number_of_iterations` | `int` | `100` | `param` | - | Number of fixed timesteps. | No convergence or steady-state stopping is performed. |
| `progress_reporting_interval` | `int` | `500` | `param` | - | Progress-report interval in iterations. | Set 0 to suppress progress reporting. |

### Outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `frequency_of_output_writing` | `int` | `0` | `param` | - | Write elevation every N iterations; 0 disables intermediate output. | Initial and final elevation are always written. |
| `write_hillshade_rasters` | `bool` | `false` | `switch` | - | Write hillshade diagnostic rasters. | Written for the initial, intermediate and final topography, in the format set by diagnostics_write_output_format. |
| `write_drainage_area_rasters` | `bool` | `false` | `switch` | - | Write log10 drainage-area diagnostic rasters. | Written for the initial, intermediate and final topography, in the format set by diagnostics_write_output_format. Uses the routing model selected by fluvial_incision_model. |
| `diagnostics_write_output_format` | `string` | `png` | `param` | - | Format of hillshade and drainage-area diagnostics: cog, png, or both. | Applies only to write_hillshade_rasters and write_drainage_area_rasters; elevation outputs use write_output_format. |
| `intermediate_hillshade_write_path` | `string` | - | `io` | - | Optional directory for hillshade outputs. | Empty uses write_path. |
| `intermediate_drainage_area_write_path` | `string` | - | `io` | - | Optional directory for drainage-area outputs. | Empty uses write_path. |
| `write_intermediate_png_sidecars` | `bool` | `false` | `switch` | - | Keep PNG georeferencing sidecars. | False removes .wld/.pgw and .aux.xml sidecars after writing. |

<!-- END GENERATED CLI PARAMETERS -->
