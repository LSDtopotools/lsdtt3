# CLI reference: lsdtt-basin-surface-metrics

Extracts predefined basins from outlet coordinates and computes multiscale surface metrics (slope, curvature, roughness, etc.) along basin outlines and within basin interiors. Outputs per-basin statistics and optional polygon layers with metrics as attributes.

## Usage

```
lsdtt-basin-surface-metrics  <path> <parameter_file>
lsdtt-basin-surface-metrics  <parameter_file>        # path defaults to ./
lsdtt-basin-surface-metrics  -h                      # print help and exit
lsdtt-basin-surface-metrics  -v                      # print version and exit
lsdtt-basin-surface-metrics  -c                      # print citation and exit
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
| 1. DEM preparation | Load DEM, verify projection, optionally fill pits |
| 2. D8 flow routing | Build the D8 flow-info object |
| 3. Basin outlet parsing | Read outlet coordinates from CSV |
| 4. Basin extraction | Extract basin boundaries and outlines |
| 5. Metric computation | Compute surface metrics at multiple scales |
| 6. Outputs | Write basin statistics and optional polygons |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input DEM filename (with or without extension). Any GDAL-supported format. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `basin_outlet_fname` | string | _(empty)_ | Path to CSV file with basin outlets. Required columns: `latitude`, `longitude`. Optional: `gauge_id`, `gauge_name`. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` = all available hardware threads. |
| `boundary_conditions` | string | `oooo` | Four-character boundary condition string (N, E, S, W). `o` = open, `n` = no-flux, `p` = periodic. |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `raster_is_filled` | bool | `false` | Set `true` if DEM is already pit-filled. If `false`, `fill_pits()` is called automatically. |
| `fill_min_slope` | float | `0.0001` | Minimum slope (m/m) imposed during pit filling. Used only when `raster_is_filled` is `false`. |
| `threshold_contributing_pixels` | int | `1000` | Contributing-pixel threshold defining channel sources used to snap basin outlets. |

---

## Basin extraction

| Parameter | Type | Default | Description |
|---|---|---|---|
| `snap_outlet_to_channel` | bool | `true` | Snap outlet coordinates to the nearest D8 channel pixel before basin extraction. |
| `outlet_snapping_search_radius_cells` | int | `10` | Search radius (pixels) for outlet snapping. |

---

## Surface metrics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `lambda_scales` | vector of floats | `[10, 20, 50, 100]` | Smoothing scales (metres) for polynomial surface fitting. |
| `polyfit_method` | string | `direct` | Polynomial fitting method: `direct` (rolling-moment, preferred) or `wavelet`. |
| `write_slope_raster` | bool | `false` | Write slope (gradient magnitude) raster. |
| `write_curvature_raster` | bool | `false` | Write curvature raster. |
| `write_roughness_raster` | bool | `false` | Write surface roughness raster. |

---

## Output: Basin statistics CSV

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_basin_statistics_csv` | bool | `true` | Write per-basin hydromorphic statistics CSV. |
| `calculate_shape_metrics` | bool | `false` | Compute basin planform shape metrics (elongation, circularity, etc.). |
| `calculate_centre_of_mass_bearing` | bool | `false` | Compute centre-of-mass bearing from outlet. |

Output file: `<write_fname>_basin_stats.csv`

---

## Output: Basin outlines polygon

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_basin_polygons` | bool | `false` | Write basin outlines as polygon features. |
| `write_vector_format` | string | `flatgeobuf` | Basin polygon format: `flatgeobuf` / `fgb` (preferred), `gpkg`. |

Output file: `<write_fname>_basin_polygons.<ext>`

---

## Output: Basin outline CSV

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_basin_outline_csv` | bool | `false` | Write basin outline as CSV with elevation and metrics sampled along perimeter. |

Output file: `<write_fname>_basin_outline.csv`

---

## Output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `cog` | Raster output format. `cog`: Cloud-Optimised GeoTIFF (preferred). `gtiff`: GeoTIFF. `same_as_input`: same as input DEM format. |
