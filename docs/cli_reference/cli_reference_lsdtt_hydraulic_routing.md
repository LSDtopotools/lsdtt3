# CLI reference: lsdtt-hydraulic-routing

Performs graph-based hydraulic routing on a DEM using D-infinity or MFD flow routing. Computes water depth, water-surface elevation, and discharge flux through time. Implements Gailleton et al. (2024) hydraulic algorithm.

## Usage

```
lsdtt-hydraulic-routing  <path> <parameter_file>
lsdtt-hydraulic-routing  <parameter_file>        # path defaults to ./
lsdtt-hydraulic-routing  -h                      # print help and exit
lsdtt-hydraulic-routing  -v                      # print version and exit
lsdtt-hydraulic-routing  -c                      # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

Keys use `snake_case`. All keys and their defaults are listed below.

---

## Pipeline overview

Steps run in this order:

| Step | Description |
|---|---|
| 1. DEM preparation | Load DEM, verify projection, fill pits |
| 2. Flow graph construction | Build D-inf or MFD flow DAG (CSR layout) |
| 3. Initial condition | Set water depth at sources or inlets |
| 4. Time-stepping | Solve shallow-water equations for N timesteps |
| 5. Outputs | Write water depth, surface, and flux rasters |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input DEM filename. Any GDAL-supported format. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `upstream_discharge_csv` | string | _(empty)_ | Optional CSV with upstream inflow points. Columns: `latitude`, `longitude`, `discharge_m3s`. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` = all available hardware threads. |
| `boundary_conditions` | string | `oooo` | Four-character boundary condition string (N, E, S, W). `o` = open, `n` = no-flux. |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `raster_is_filled` | bool | `false` | Set `true` if DEM is already pit-filled. |
| `fill_min_slope` | float | `0.0001` | Minimum slope (m/m) imposed during pit filling. |

---

## Flow routing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `flow_routing_algorithm` | string | `dinf` | Flow routing method: `dinf` (D-infinity), `mfd` (multiple-flow-direction). |

---

## Initial conditions

| Parameter | Type | Default | Description |
|---|---|---|---|
| `initial_water_depth_metres` | float | `0.1` | Uniform initial water depth (metres) across domain. |
| `water_depth_at_channel_threshold_area_m2` | float | `10000.0` | Drainage area threshold (m²) defining channels. Channels get initial water depth. |

---

## Time-stepping

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_timesteps` | int | `10` | Number of time steps to solve. |
| `timestep_seconds` | float | `1.0` | Duration of each timestep (seconds). |
| `manning_coefficient` | float | `0.03` | Manning's n roughness coefficient. |
| `gravity_acceleration_m_s2` | float | `9.81` | Gravitational acceleration (m/s²). |

---

## Solver options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `solver_tolerance` | float | `1e-6` | Convergence tolerance for iterative solver. |
| `max_solver_iterations` | int | `100` | Maximum iterations per timestep. |

---

## Output: Water depth

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_water_depth` | bool | `true` | Write water depth raster for final timestep. |
| `write_water_depth_all_timesteps` | bool | `false` | Write water depth raster for every timestep. |

Output file: `<write_fname>_water_depth_final.tif` (or `_ts<N>.tif` for each timestep)

---

## Output: Water surface elevation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_water_surface` | bool | `true` | Write water-surface elevation raster. |

Output file: `<write_fname>_water_surface_final.tif`

---

## Output: Discharge flux

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_discharge_flux` | bool | `true` | Write discharge flux magnitude raster. |

Output file: `<write_fname>_discharge_flux_final.tif`

---

## Output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `cog` | Raster output format: `cog`, `gtiff`, `same_as_input`. |
