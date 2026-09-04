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
