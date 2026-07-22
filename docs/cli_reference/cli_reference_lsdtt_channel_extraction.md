# CLI reference: lsdtt-channel-extraction

Extracts a channel network from a DEM using D8 flow routing and either threshold-based or Wiener-based source extraction. Outputs the network as points, polylines, junctions, or a flow-direction raster. Optionally burns a user-supplied vector channel network into the DEM and writes elevation profiles.

## Usage

```
lsdtt-channel-extraction  <path> <parameter_file>
lsdtt-channel-extraction  <parameter_file>        # path defaults to ./
lsdtt-channel-extraction  -h                      # print help and exit
lsdtt-channel-extraction  -v                      # print version and exit
lsdtt-channel-extraction  -c                      # print citation and exit
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
| 1. Pit filling | Optional; always recommended unless DEM is pre-filled |
| 2. D8 flow routing | Build the D8 flow-info object |
| 3. Source extraction | Find channel heads using `source_extraction_algorithm_choice` (`threshold` or `wiener`) |
| 4. Junction network | Build compact channel link and junction objects |
| 5. Outputs | Write selected outputs |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input DEM filename (with or without extension). Any GDAL-supported format. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` = all available hardware threads. Applies to pit filling and flow-routing construction. |
| `boundary_conditions` | string | `oooo` | Four-character boundary condition string (N, E, S, W). `o` = open, `n` = no-flux, `p` = periodic, `b` = forced outlet. |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `raster_is_filled` | bool | `false` | Set `true` if the DEM is already pit-filled and ready for routing. If `false`, `fill_pits()` is called automatically. |
| `fill_min_slope` | float | `0.0001` | Minimum slope (m/m) imposed during pit filling. Used only when `raster_is_filled` is `false`. |

---

## Source extraction

| Parameter | Type | Default | Description |
|---|---|---|---|
| `source_extraction_algorithm_choice` | string | `threshold` | Source extraction method. Allowed values: `threshold`, `wiener`. |
| `source_pixels_threshold` | int | `1000` | Contributing-pixel count required to define a channel head. Reduce to extract a denser network. |

### Wiener parameters

Used only when `source_extraction_algorithm_choice: wiener`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `wiener_pruning_drainage_area_m2` | float | `1000.0` | Minimum drainage area for Wiener channel-mask pruning, in physical square metres (m^2). |
| `wiener_surface_fitting_radius_m` | float | `6.0` | Radius (metres) used when fitting surfaces for tangential-curvature estimation. |
| `wiener_connected_components_threshold_pixels` | int | `100` | Minimum connected-component size (pixels) retained in the Wiener channel mask before skeletonisation. |

### Wiener diagnostics (optional)

All diagnostics are opt-in and do not alter source extraction output.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `wiener_write_filtered_dem` | bool | `false` | Write Wiener-filtered DEM raster as `<write_fname>_wiener_filtered_dem.*`. |
| `wiener_write_tangential_curvature` | bool | `false` | Write tangential-curvature raster as `<write_fname>_wiener_tangential_curvature.*`. |
| `wiener_write_qq_table` | bool | `false` | Write QQ diagnostics table CSV. |
| `wiener_diagnostics_raster_format` | string | `same_as_input` | Raster format for Wiener diagnostics (`cog`, `gtiff`, `same_as_input`, `png`). |
| `wiener_qq_output_csv` | string | _(empty)_ | Explicit QQ CSV path. If empty and `wiener_write_qq_table` is true, output defaults to `<write_fname>_wiener_qq.csv`. |

---

## Output: general

| Parameter | Type | Default | Description |
|---|---|---|---|
| `calculate_spatial_index` | bool | `true` | Build a spatial index for vector outputs. Speeds up GIS spatial queries but requires all features to be buffered in memory before writing. Disable for large datasets. |

---

## Output: channel network points

One point per channel pixel, with drainage area as an attribute.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_channel_network_points` | bool | `false` | Write the channel network as a point layer (one feature per channel pixel). |
| `channel_network_points_output_format` | string | `flatgeobuf` | Output format. `csv`: latitude/longitude in EPSG:4326. `gpkg`: GeoPackage. `flatgeobuf` / `fgb`: FlatGeobuf (preferred on WSL2). |

Output file: `<write_fname>_channel_network.<ext>`

---

## Output: channel network lines

One polyline per channel link (source-to-junction or junction-to-junction). Smaller file than the point layer; does not record per-pixel drainage area.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_channel_network_lines` | bool | `false` | Write the channel network as a polyline layer. |
| `channel_network_lines_output_format` | string | `flatgeobuf` | Output format. `flatgeobuf` / `fgb`: FlatGeobuf. `gpkg`: GeoPackage. CSV is not supported for lines. |

Output file: `<write_fname>_channel_network_lines.<ext>`

---

## Output: junctions

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_channel_network_junctions` | bool | `false` | Write junction points with attributes `junction_number`, `receiver_junction`, and `stream_order`. |
| `channel_network_junctions_output_format` | string | `flatgeobuf` | Output format. `csv`, `gpkg`, `flatgeobuf` / `fgb`. |

Output file: `<write_fname>_junctions.<ext>`

---

## Output: D8 flow lines

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_flow_lines` | bool | `false` | Write one line from each valid pixel centre to its D8 receiver. |
| `flow_lines_output_format` | string | `flatgeobuf` | Output format. `flatgeobuf` / `fgb`, `gpkg`. CSV is not supported. |

Output file: `<write_fname>_flow_lines.<ext>`

---

## Output: ArcMap D8 flow direction raster

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_arcmap_d8_flowdirection_raster` | bool | `false` | Write a D8 flow-direction raster using ArcMap encoding: E=1, SE=2, S=4, SW=8, W=16, NW=32, N=64, NE=128. Outlet pixels are written as 0. |
| `flow_direction_output_format` | string | `cog` | Raster format. `cog`, `gtiff`, or `same_as_input`. |

Output file: `<write_fname>_d8_flow_direction.<ext>`

---

## Output: burned stream profiles

Burns a user-supplied vector polyline network into the DEM and writes an elevation profile showing original and burned elevations. Only `burn_in_place` mode is currently implemented.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_burned_stream_profiles` | bool | `false` | Burn vector channels into DEM and write elevation profiles. Requires `stream_vector_path`. |
| `stream_vector_path` | string | _(empty)_ | Path to the vector file containing polyline channel features. Must be in the same CRS as the DEM. Supported: GeoPackage, FlatGeobuf, Shapefile, and any GDAL/OGR format. |
| `stream_burn_mode` | string | `burn_in_place` | Burn algorithm. `burn_in_place`: samples DEM pixels beneath the input polyline and enforces monotone-decreasing elevation downstream. `snap_to_network`: not yet implemented. |
| `burned_stream_profiles_output_format` | string | `csv` | Output format for the profile table. `csv`: plain table with one row per channel pixel. `gpkg`: GeoPackage point layer. `flatgeobuf` / `fgb`: FlatGeobuf point layer. |

Output file: `<write_fname>_burned_stream_profiles.<ext>`

Profile columns: `feature_id`, `flow_distance`, `x`, `y`, `dem_elevation`, `burned_elevation`.

---

## Minimal example parameter file

```
# lsdtt-channel-extraction parameter file
read_path:  /data/dems/
read_fname: my_dem.tif

source_pixels_threshold:          500
write_channel_network_lines:      true
write_channel_network_junctions:  true
```

## Example: threshold source extraction

```
read_path: /data/dems/
read_fname: my_dem.tif

source_extraction_algorithm_choice: threshold
source_pixels_threshold: 750

write_channel_network_lines: true
```

## Example: Wiener source extraction with diagnostics

```
read_path: /data/dems/
read_fname: my_dem.tif

source_extraction_algorithm_choice: wiener
wiener_pruning_drainage_area_m2: 2500.0
wiener_surface_fitting_radius_m: 8.0
wiener_connected_components_threshold_pixels: 80

wiener_write_filtered_dem: true
wiener_write_tangential_curvature: true
wiener_write_qq_table: true
wiener_diagnostics_raster_format: gtiff

write_channel_network_points: true
```
