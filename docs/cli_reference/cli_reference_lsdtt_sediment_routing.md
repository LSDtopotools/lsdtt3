# CLI reference: lsdtt-sediment-routing

Routes sediment through a landscape using flow routing and sediment transport laws. Computes sediment flux, deposition, and erosion/incision across the DEM using either detachment-limited or transport-limited formulations.

## Usage

```
lsdtt-sediment-routing  <path> <parameter_file>
lsdtt-sediment-routing  <parameter_file>        # path defaults to ./
lsdtt-sediment-routing  -h                      # print help and exit
lsdtt-sediment-routing  -v                      # print version and exit
lsdtt-sediment-routing  -c                      # print citation and exit
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
| 2. Flow routing | Build D8 or D-inf flow network |
| 3. Accumulation | Compute drainage area and flow length |
| 4. Transport law | Solve transport equation (detachment or transport-limited) |
| 5. Output | Write sediment flux, deposition, and topographic change |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input DEM filename. Any GDAL-supported format. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` = all available hardware threads. |
| `boundary_conditions` | string | `oooo` | Four-character boundary condition string (N, E, S, W). |

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
| `flow_routing_algorithm` | string | `d8` | Flow routing method: `d8` (D8, preferred for sediment routing), `dinf` (D-infinity). |

---

## Sediment transport law

| Parameter | Type | Default | Description |
|---|---|---|---|
| `sediment_transport_model` | string | `stream_power` | Transport law: `stream_power` (U = K·A^m·S^n), `detachment_limited`, `transport_limited`. |
| `stream_power_K` | float | `1e-5` | Stream power erodibility coefficient (m/a). |
| `stream_power_m` | float | `0.5` | Area exponent. |
| `stream_power_n` | float | `1.0` | Slope exponent. |
| `slope_clamp_factor` | float | `0.0` | Maximum slope permissible (multiples of `dx`). `0` = no clamping. |

---

## Threshold area

| Parameter | Type | Default | Description |
|---|---|---|---|
| `threshold_area_m2` | float | `1000000.0` | Minimum drainage area (m²) for fluvial transport. Hillslope above this threshold. |
| `divide_area_by_dx` | bool | `true` | Divide area by pixel size (`cell_size`) in transport law (common in some formulations). |

---

## Sediment flux computation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `compute_sediment_flux` | bool | `true` | Compute and write sediment flux raster. |
| `compute_deposition_raster` | bool | `false` | Compute deposition locations. |

---

## Output: Sediment flux

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_sediment_flux_raster` | bool | `true` | Write sediment flux magnitude raster. |

Output file: `<write_fname>_sediment_flux.tif`

---

## Output: Erosion and deposition

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_erosion_rate_raster` | bool | `false` | Write erosion rate raster (m/a). |
| `write_deposition_rate_raster` | bool | `false` | Write deposition rate raster (m/a). |

Output files:
- `<write_fname>_erosion_rate.tif`
- `<write_fname>_deposition_rate.tif`

---

## Output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `cog` | Raster output format: `cog`, `gtiff`, `same_as_input`. |
