# CLI reference: lsdtt-graphflood

Runs stationary GraphFlood routing on a DEM. Routing is defined by the evolving
hydraulic surface, with node states selecting the hydraulic domain and its
outlets.

## Usage

```
lsdtt-graphflood <path> <parameter_file>
lsdtt-graphflood <parameter_file>
```

## Required inputs

| Parameter | Description |
|---|---|
| `read_fname` | Bed-elevation raster. |
| `hydraulic_node_state_fname` | Co-registered node-state raster: `0` wall/nodata, `1` active, `3` outlet. |
| `rainfall_rate_fname` | Co-registered rainfall/input-rate raster in m/s. |
| `friction_coefficient_fname` | Co-registered coefficient raster for `q / width = coefficient * h^a * S^b`. Use `1 / n` for Manning friction. |

`boundary_conditions` defaults to `oooo` (NESW). It controls whether an edge
can communicate beyond the raster (`o`), is closed (`n`), or is periodic (`p`).
It does not create an outlet: only a node-state value of `3` does that.

## Hydraulic parameters

| Parameter | Default | Description |
|---|---:|---|
| `hydraulic_dt` | `1e-3` | Fixed timestep in seconds. |
| `hydraulic_min_slope` | `1e-4` | Minimum hydraulic slope used in friction calculations. |
| `hydraulic_fill_slope_epsilon` | `1e-8` | Requested minimum slope while priority-flooding `z + h`; a representable float32 rise is always enforced. |
| `friction_depth_exponent` | `5/3` | Depth exponent `a`. |
| `friction_slope_exponent` | `1/2` | Slope exponent `b`. |
| `n_steps` | `100` | Iteration limit. |

There is no scalar rainfall fallback, scalar Manning parameter, or water-depth
rate limiter. Supply the forcing rasters explicitly.

## Outputs

`write_water_depth_raster` and `write_water_surface_raster` default to true.
Optional outputs include `write_qin_raster`, `write_qout_raster`,
`write_dhw_raster`, `write_hydraulic_slope_raster`,
`write_specific_discharge_raster`, `write_total_discharge_raster`, and
`write_effective_width_raster`.

<!-- BEGIN GENERATED CLI PARAMETERS -->
## Parameters

This section is generated from the executable's `-h json` registry. Do not edit it manually.

**Status: experimental.** This CLI is not yet fully tested.

### File I/O

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `read_path` | `string` | `./` | `io` | - | Directory from which input files are read. | Defaults to the directory holding the parameter file. |
| `read_fname` | `string` | - | `io` | - | Input filename, with or without extension (e.g. my_dem.tif). | Read from read_path. lsdtt3 uses GDAL, so any GDAL-readable raster format works. If empty, defaults to the parameter filename stem. |
| `write_path` | `string` | - | `io` | - | Directory to which output files are written. | Defaults to read_path when not set. |
| `write_fname` | `string` | - | `io` | - | Output filename prefix (no path, no extension). | All outputs are named &lt;write_fname&gt;_&lt;descriptor&gt;.&lt;ext&gt;. Defaults to the stem of read_fname when not set. |

### Hydraulic domain

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `hydraulic_node_state_fname` | `string` | - | `io` | - | Raster defining hydraulic nodes: 0 = outside, 1 = active, 3 = outlet. | Required. Every active connected component needs code 3 by default; code 3 must lie on an open boundary edge. |
| `check_valid_outlet` | `bool` | `true` | `switch` | - | Require every active connected component to contain a code-3 outlet. | Disable only for deliberate experimental no-drain configurations. |
| `bed_preprocessing` | `string` | `none` | `param` | - | Optional conditioning of the bed before hydraulic routing. | Allowed values: none or fill. Fill honours node states and boundary conditions. |
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope imposed by pit filling. | Used only for bed_preprocessing = fill. |
| `hydraulic_fill_slope_epsilon` | `float` | `1e-8` | `param` | - | Minimum slope imposed when filling the evolving z+h hydraulic surface. | A small positive value removes routing flats while avoiding artificial hydraulic gradients. |

### Flowgraph routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four-character string for north, east, south, west. o = potentially draining open edge, n = closed edge, p = periodic edge. Only selected code-3 cells drain through an o edge. |
| `flow_routing_method` | `string` | `mfd` | `param` | - | Flowgraph routing method. | Allowed values: mfd or dinf. |
| `mfd_exponent` | `float` | `1.0` | `param` | - | Slope-weighting exponent for MFD routing. | Used only when flow_routing_method = mfd. |
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all available hardware threads. |

### Hydraulic parameters

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `hydraulic_dt` | `float` | `1e-3` | `param` | - | Fixed hydraulic timestep in seconds. | Used for every hydraulic update. |
| `friction_coefficient_fname` | `string` | - | `io` | - | Raster of wide-flow friction coefficient. | Required. For Manning flow this is 1/n per cell. |
| `friction_depth_exponent` | `float` | `1.6666667` | `param` | - | Depth exponent in q/width = coefficient * h^a * S^b. | Manning uses 5/3; Chézy uses 3/2. |
| `friction_slope_exponent` | `float` | `0.5` | `param` | - | Slope exponent in q/width = coefficient * h^a * S^b. | Manning and Chézy both use 1/2. |
| `hydraulic_min_slope` | `float` | `1e-4` | `param` | - | Minimum hydraulic slope used internally. | Prevents zero-gradient numerical stagnation. |
| `rainfall_rate_fname` | `string` | - | `io` | - | Raster of rainfall/input rate in m/s. | Required and co-registered with the DEM. |
| `receiver_width` | `float` | `-1.0` | `param` | - | Effective receiver width in metres. | If &lt;= 0, the solver falls back to a geometry-based default (currently cell_size). |
| `n_steps` | `int` | `100` | `param` | - | Maximum stationary hydraulic iterations to run. | Used as the convergence limit or the fixed iteration count. |

### Point inflow

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `use_point_inflow` | `bool` | `false` | `switch` | - | Apply a single point inflow to the hydraulic solver. | If true, the solver adds a fixed discharge at point_inflow_flat_index every timestep. |
| `point_inflow_flat_index` | `int` | `0` | `param` | - | Flat raster index of the point inflow cell. | Flat index is row * n_cols + col in row-major order. |
| `point_inflow_discharge_m3_per_s` | `float` | `0.0` | `param` | - | Point inflow discharge in m^3/s. | Used only when use_point_inflow is true. |
| `run_until_converged` | `bool` | `true` | `switch` | - | Run until convergence instead of a fixed number of steps. | Used mainly for stationary mode. |
| `report_every` | `int` | `50` | `param` | - | Iteration interval for progress reporting. | Set &lt;= 0 to disable periodic progress output. |
| `convergence_check_every` | `int` | `10` | `param` | - | Iteration interval for convergence checking. | Used by run_until_converged. |
| `median_dhw_convergence_tolerance_metres` | `float` | `1e-6` | `param` | - | Absolute median dhw convergence threshold in metres. | Used by run_until_converged. |
| `median_dhw_relative_convergence_tolerance` | `float` | `1e-4` | `param` | - | Relative median dhw convergence threshold. | Used by run_until_converged. |
| `adaptive_dt_enabled` | `bool` | `false` | `param` | - | Enable simple adaptive timestep control in stationary convergence runs. | Adaptive dt only acts inside run_until_converged. |
| `adaptive_dt_min` | `float` | `1e-4` | `param` | - | Minimum adaptive timestep in seconds. | Used only when adaptive_dt_enabled is true. |
| `adaptive_dt_max` | `float` | `1.0` | `param` | - | Maximum adaptive timestep in seconds. | Used only when adaptive_dt_enabled is true. |
| `adaptive_target_median_abs_dhw` | `float` | `1e-3` | `param` | - | Target median absolute water-depth change per checked iteration in metres. | Adaptive dt grows or shrinks to keep median \|dhw\| near this target. |
| `adaptive_dt_growth_factor` | `float` | `1.5` | `param` | - | Multiplicative factor used to increase adaptive dt. | Must be &gt; 1. |
| `adaptive_dt_shrink_factor` | `float` | `0.5` | `param` | - | Multiplicative factor used to decrease adaptive dt. | Must be &gt; 0 and &lt; 1. |

### Output controls

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_water_depth_raster` | `bool` | `true` | `switch` | - | Write final water depth raster. | Output suffix: _water_depth.tif |
| `write_water_surface_raster` | `bool` | `true` | `switch` | - | Write final water surface raster. | Output suffix: _water_surface.tif |
| `write_qin_raster` | `bool` | `false` | `switch` | - | Write final qin raster. | Output suffix: _qin.tif |
| `write_qout_raster` | `bool` | `false` | `switch` | - | Write final qout raster. | Output suffix: _qout.tif |
| `write_dhw_raster` | `bool` | `false` | `switch` | - | Write final timestep water-depth increment raster. | Output suffix: _dhw.tif |
| `write_bed_raster` | `bool` | `false` | `switch` | - | Write bed elevation raster used by the hydraulic solver. | Output suffix: _bed.tif |
| `write_output_format` | `string` | `cog` | `io` | - | Output format for raster outputs. | Allowed values: cog, gtiff, same_as_input. |

<!-- END GENERATED CLI PARAMETERS -->
