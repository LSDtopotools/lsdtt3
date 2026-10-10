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

!!! note "Multipart lines"
    `stream_vector_path` accepts multipart geometry (`MultiLineString` and friends). Each feature is split into its parts on read, with the attributes copied onto each, and the split is reported in the log.

    Shapefiles are the commonest source of multipart data, and **versions before August 2026 discarded it silently** — see [Vector data](../what_is_lsdtt3/vector_data.md).

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
| `threshold_contributing_pixels` | int | `1000` | Contributing-pixel count required to define a channel head. Reduce to extract a denser network. |

### Wiener parameters

Used only when `source_extraction_algorithm_choice: wiener`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `wiener_pruning_drainage_area_m2` | float | `1000.0` | Minimum drainage area for Wiener channel-mask pruning, in physical square metres (m^2). |
| `wiener_surface_fitting_radius_metres` | float | `6.0` | Radius (metres) used when fitting surfaces for tangential-curvature estimation. |
| `wiener_connected_components_threshold_pixels` | int | `100` | Minimum connected-component size (pixels) retained in the Wiener channel mask before skeletonisation. |

### Wiener diagnostics (optional)

All diagnostics are opt-in and do not alter source extraction output.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_wiener_filtered_dem_raster` | bool | `false` | Write Wiener-filtered DEM raster as `<write_fname>_wiener_filtered_dem.*`. |
| `write_wiener_tangential_curvature_raster` | bool | `false` | Write tangential-curvature raster as `<write_fname>_wiener_tangential_curvature.*`. |
| `wiener_write_qq_table` | bool | `false` | Write QQ diagnostics table CSV. |
| `wiener_diagnostics_write_output_format` | string | `same_as_input` | Raster-format override for Wiener diagnostics (`cog`, `gtiff`, `same_as_input`, `png`). |
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
| `channel_network_points_write_vector_format` | string | `same_as_write_vector_format` | Point-output override. `csv`: latitude/longitude in EPSG:4326. `gpkg`: GeoPackage. `flatgeobuf` / `fgb`: FlatGeobuf. |

!!! tip "Longitudinal profiles come from the POINTS, not the lines"
    Each channel point carries `flow_distance` (distance to the outlet along the
    routed flow path, in metres) and `elevation`. Plot one against the other for
    a longitudinal profile.

    **Do not derive distance from the line geometry.** `flow_distance` is
    accumulated from the D8 direction code — one cell size per cardinal step,
    `cell_size × √2` per diagonal — so it is correct on a periodic domain and is
    unaffected by how the lines are packaged or split. See
    [Vector data](../what_is_lsdtt3/vector_data.md).

Output file: `<write_fname>_channel_network.<ext>`

---

## Output: channel network lines

One polyline per channel link (source-to-junction or junction-to-junction). Smaller file than the point layer; does not record per-pixel drainage area.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_channel_network_lines` | bool | `false` | Optional channel-network polyline export; off because channel-extraction vector products are opt-in. |
| `write_vector_format` | string | `flatgeobuf` | Baseline channel-line format. `flatgeobuf` / `fgb` or `gpkg`; CSV is not supported for lines. |

Output file: `<write_fname>_channel_network_lines.<ext>`

!!! note "The lines are MultiLineString"
    Channel links are written as `MultiLineString`, uniformly. An ordinary link
    is a single-part one whose coordinates are exactly what a `LineString` would
    have carried. A link that crosses a **periodic** boundary (available in
    `lsdtt-lem`) is split at the wrap, so no renderer draws a line across the
    whole domain. `link_length` is measured along the flow path, so it is
    correct in either case. See [Vector data](../what_is_lsdtt3/vector_data.md).

---

## Output: junctions

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_channel_network_junctions` | bool | `false` | Write junction points with attributes `junction_number`, `receiver_junction`, and `stream_order`. |
| `channel_network_junctions_write_vector_format` | string | `same_as_write_vector_format` | Junction-point override: `csv`, `gpkg`, `flatgeobuf` / `fgb`. |

Output file: `<write_fname>_junctions.<ext>`

---

## Output: D8 flow lines

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_flow_lines` | bool | `false` | Write one line from each valid pixel centre to its D8 receiver. |
| `flow_lines_write_vector_format` | string | `same_as_write_vector_format` | Flow-line override. `flatgeobuf` / `fgb`, `gpkg`. CSV is not supported. |

Output file: `<write_fname>_flow_lines.<ext>`

---

## Output: ArcMap D8 flow direction raster

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_arcmap_d8_flowdirection_raster` | bool | `false` | Write a D8 flow-direction raster using ArcMap encoding: E=1, SE=2, S=4, SW=8, W=16, NW=32, N=64, NE=128. Outlet pixels are written as 0. |
| `write_output_format` | string | `cog` | Baseline raster format. `cog`, `gtiff`, or `same_as_input`. |

Output file: `<write_fname>_d8_flow_direction.<ext>`

---

## Output: burned stream profiles

Burns a user-supplied vector polyline network into the DEM and writes an elevation profile showing original and burned elevations. Only `burn_in_place` mode is currently implemented.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_burned_stream_profiles` | bool | `false` | Burn vector channels into DEM and write elevation profiles. Requires `stream_vector_path`. |
| `stream_vector_path` | string | _(empty)_ | Path to the vector file containing polyline channel features. Must be in the same CRS as the DEM. Supported: GeoPackage, FlatGeobuf, Shapefile, and any GDAL/OGR format. |
| `stream_burn_mode` | string | `burn_in_place` | Burn algorithm. `burn_in_place`: samples DEM pixels beneath the input polyline and enforces monotone-decreasing elevation downstream. `snap_to_network`: not yet implemented. |
| `burned_stream_profiles_write_vector_format` | string | `csv` | Output format for the profile table. `csv`: plain table with one row per channel pixel. `gpkg`: GeoPackage point layer. `flatgeobuf` / `fgb`: FlatGeobuf point layer. |

Output file: `<write_fname>_burned_stream_profiles.<ext>`

Profile columns: `feature_id`, `flow_distance`, `x`, `y`, `dem_elevation`, `burned_elevation`.

---

## Minimal example parameter file

```
# lsdtt-channel-extraction parameter file
read_path:  /data/dems/
read_fname: my_dem.tif

threshold_contributing_pixels: 500
write_channel_network_lines:      true
write_channel_network_junctions:  true
```

## Example: threshold source extraction

```
read_path: /data/dems/
read_fname: my_dem.tif

source_extraction_algorithm_choice: threshold
threshold_contributing_pixels: 750

write_channel_network_lines: true
```

## Example: Wiener source extraction with diagnostics

```
read_path: /data/dems/
read_fname: my_dem.tif

source_extraction_algorithm_choice: wiener
wiener_pruning_drainage_area_m2: 2500.0
wiener_surface_fitting_radius_metres: 8.0
wiener_connected_components_threshold_pixels: 80

write_wiener_filtered_dem_raster: true
write_wiener_tangential_curvature_raster: true
wiener_write_qq_table: true
wiener_diagnostics_write_output_format: gtiff

write_channel_network_points: true
```

<!-- BEGIN GENERATED CLI PARAMETERS -->
## Parameters

This section is generated from the executable's `-h json` registry. Do not edit it manually.

### File I/O

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `read_path` | `string` | `./` | `io` | - | Directory from which input files are read. | Defaults to the directory holding the parameter file. |
| `read_fname` | `string` | - | `io` | - | Input filename, with or without extension (e.g. my_dem.tif). | Read from read_path. lsdtt3 uses GDAL, so any GDAL-readable raster format works. If empty, defaults to the parameter filename stem. |
| `write_path` | `string` | - | `io` | - | Directory to which output files are written. | Defaults to read_path when not set. |
| `write_fname` | `string` | - | `io` | - | Output filename prefix (no path, no extension). | All outputs are named &lt;write_fname&gt;_&lt;descriptor&gt;.&lt;ext&gt;. Defaults to the stem of read_fname when not set. |

### Run control

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all available hardware threads. Applies to pit filling and flow-routing construction. |
| `calculate_spatial_index` | `bool` | `true` | `io` | - | Whether or not to calculate a spatial index for vector data writes. | A spatial index helps a GIS to speed up spatial queries, so it is often helpful to have a spatial index calculated. But calculating this requires the entire vector layer to sit in memory, preventing feature-by-feature writes and adding a lot of overhead. In some cases (e.g. large datasets) it may therefore be preferred to generate vectors without spatial index - these can always be calculated later on using e.g. gdal. |

### Flow routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four-character string for north, east, south, west. o = open, n = no-flux, p = periodic, b = forced outlet. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope (m/m) imposed by pit filling or flat resolution. | inline_gradient (default): imposed during the Priority-Flood fill, so used only when raster_is_filled is false. barnes2014b_gradient: the minimum gradient across flats, also on a raster_is_filled input. barnes2014b_float_step: ignored. |
| `fill_flat_resolution_choice` | `string` | `inline_gradient` | `param` | - | How flats are made drainable: inline_gradient, barnes2014b_gradient or barnes2014b_float_step. | inline_gradient (default, unchanged): the Priority-Flood fill imposes fill_min_slope during the flood. The barnes2014b options make flow converge across flats (Barnes et al. 2014, doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses fill_min_slope (m/m, must be &gt; 0) as the minimum gradient; flats are raised about 3-4x more, capped below their rims, and capped flats are reported. barnes2014b_float_step is equivalent to RichDEM ResolveFlats on float32 in its first pass: it ignores fill_min_slope and its slopes depend on elevation (of order ULP(z) / cell size). With either barnes2014b option the result is always routable: cells the first pass leaves without a lower neighbour are refilled and resolved again (up to 4 rounds), and any still left are raised by a final priority-flood pass by one float32 step each; the report gives the rounds and the cells changed after the first pass. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Channel sources

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `source_extraction_algorithm_choice` | `string` | `threshold` | `switch` | - | Source extraction algorithm. | Allowed values: threshold, wiener. Threshold uses threshold_contributing_pixels. Wiener uses filtered-DEM curvature, QQ thresholding, and area pruning. |
| `wiener_pruning_drainage_area_m2` | `float` | `1000.0` | `param` | `source_extraction_algorithm_choice` | Minimum drainage area for Wiener channel pixels (m^2). | Wiener-only parameter. Channel pixels with drainage area below this physical-area threshold are removed before skeleton endpoint extraction. |
| `wiener_surface_fitting_radius_metres` | `float` | `6.0` | `param` | `source_extraction_algorithm_choice` | Surface-fitting radius for Wiener tangential curvature (m). | Wiener-only parameter. Passed to local polynomial fitting used for tangential curvature computation. |
| `wiener_connected_components_threshold_pixels` | `int` | `100` | `param` | `source_extraction_algorithm_choice` | Minimum connected-component size in Wiener channel mask (pixels). | Wiener-only parameter. Components smaller than this are removed before skeletonisation. |
| `write_wiener_filtered_dem_raster` | `bool` | `false` | `param` | `source_extraction_algorithm_choice` | Write Wiener-filtered DEM diagnostic raster. | If true, writes &lt;write_fname&gt;_wiener_filtered_dem with the format set by wiener_diagnostics_write_output_format. |
| `write_wiener_tangential_curvature_raster` | `bool` | `false` | `param` | `source_extraction_algorithm_choice` | Write Wiener tangential-curvature diagnostic raster. | If true, writes &lt;write_fname&gt;_wiener_tangential_curvature with the format set by wiener_diagnostics_write_output_format. |
| `wiener_write_qq_table` | `bool` | `false` | `param` | `source_extraction_algorithm_choice` | Write Wiener QQ diagnostic table as CSV. | If true, writes QQ data to wiener_qq_output_csv, or to &lt;write_fname&gt;_wiener_qq.csv when that parameter is empty. |
| `wiener_diagnostics_write_output_format` | `string` | `same_as_input` | `param` | `source_extraction_algorithm_choice` | Raster format for Wiener diagnostic rasters. | Allowed values: cog, gtiff, same_as_input, png. |
| `wiener_qq_output_csv` | `string` | - | `param` | `source_extraction_algorithm_choice` | Output CSV path for Wiener QQ diagnostics. | Wiener-only parameter. If empty and wiener_write_qq_table is true, the CLI writes &lt;write_fname&gt;_wiener_qq.csv. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_output_format` | `string` | `cog` | `io` | - | Raster output format: cog, gtiff, png or same_as_input. | cog is a Cloud Optimised GeoTIFF and is the project default. It falls back to gtiff for small rasters, where the tiling overhead is not worth it. There is ONE raster format parameter per CLI: it governs every raster the tool writes. Vector outputs are governed separately by write_vector_format. |
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). CSV is available only through selectors for point/table outputs. Per-output overrides exist only where a container genuinely differs. Raster outputs are governed separately by write_output_format. |

### Channel network output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_channel_network_points` | `bool` | `false` | `switch` | - | Write the channel network as point data, where each channel pixel is a feature with attribute values, such as drainage_area. | Compared to write_channel_network_lines, this write will be larger because there are more features, but it does record the drainage area at each channel pixel. |
| `channel_network_points_write_vector_format` | `string` | `same_as_write_vector_format` | `param` | `write_channel_network_points` | Output format for the exported channel network. | csv: writes latitude/longitude point rows in EPSG:4326. gpkg: writes a GeoPackage in the raster CRS. flatgeobuf or fgb: writes FlatGeobuf in the raster CRS. |
| `write_channel_network_lines` | `bool` | `false` | `switch` | - | Write the channel network as vector lines, where each channel link (defined between source or upstream junction and outlet or downstream junction) is a feature. | Does not support csv writes, only gpkg and flatgeobuf/fgb. It is off by default because channel-extraction vector products are opt-in. Compared to write_channel_network as points, the lines filesize will be smaller (there are fewer features), at the cost of not recording drainage area at each channel pixel. |
| `write_channel_network_junctions` | `bool` | `false` | `switch` | - | Write the channel network junctions as point data | Contains feature attributes: junction_number, receiver_junctionand stream_order |
| `channel_network_junctions_write_vector_format` | `string` | `same_as_write_vector_format` | `param` | `write_channel_network_junctions` | Output format for the exported channel network. | csv: writes latitude/longitude point rows in EPSG:4326. gpkg: writes a GeoPackage in the raster CRS. flatgeobuf or fgb: writes FlatGeobuf in the raster CRS. |

### Flow routing output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_flow_lines` | `bool` | `false` | `switch` | - | Write D8 receiver links as a vector line layer. | Writes one line from each valid donor pixel centre to its receiver pixel centre, following the D8 directions. |
| `flow_lines_write_vector_format` | `string` | `same_as_write_vector_format` | `param` | `write_flow_lines` | Output format for D8 flow-line visualisation. | Allowed values: flatgeobuf, fgb, gpkg. FlatGeobuf is the default. GeoPackage may have file I/O issues on some filesystems. |
| `write_drainage_area_raster` | `bool` | `false` | `switch` | - | Write the D8 drainage-area raster. | Contributing area at every pixel, from the same flow routing used for channel extraction. Pair it with write_distance_to_farthest_ridge_raster for the length-area (Hack) scaling, or use it to weight channel lines by discharge. Output: &lt;write_fname&gt;_drainage_area. |
| `write_distance_to_farthest_ridge_raster` | `bool` | `false` | `switch` | - | Write the D8 flow-path length from the farthest upslope headwater to each pixel. | This is the L of Hack's law (L ~ A^h); pair it with the drainage-area raster to get the length-area scaling of the network. Headwater pixels are 0. Output: &lt;write_fname&gt;_distance_to_farthest_ridge. |
| `write_arcmap_d8_flowdirection_raster` | `bool` | `false` | `switch` | - | Write a D8 flow-direction raster using ArcMap conventions. | ArcMap D8 codes are used: E=1, SE=2, S=4, SW=8, W=16, NW=32, N=64, NE=128. Outlet pixels are written as 0. |

### Stream burning

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_burned_stream_profiles` | `bool` | `false` | `switch` | - | Burn a user-supplied vector polyline into the DEM and write an elevation profile showing original and burned elevations. | Requires stream_vector_path to be set. The routine finds pixels beneath each polyline feature, builds a flow-distance/elevation profile, and enforces a downhill-flow constraint on the burned elevation. Output columns: feature_id, flow_distance, x, y, dem_elevation, burned_elevation. |
| `stream_vector_path` | `string` | - | `param` | `write_burned_stream_profiles` | Path to the input vector file containing polyline channel features. | Required when write_burned_stream_profiles is true. The file must contain polyline or multilinestring geometries in the same CRS as the DEM. Supported formats: GeoPackage (.gpkg), FlatGeobuf (.fgb), Shapefile (.shp), and any other GDAL/OGR-supported format. |
| `stream_burn_mode` | `string` | `burn_in_place` | `param` | `write_burned_stream_profiles` | Method used to construct the burned channel profile. | burn_in_place: pixels directly beneath the input polyline are sampled from the DEM and a monotone-decreasing burned elevation is computed by linking low-lying points along the profile. snap_to_network: not yet implemented. |
| `burned_stream_profiles_write_vector_format` | `string` | `csv` | `param` | `write_burned_stream_profiles` | Output format for the burned stream profile table. | csv: writes a plain CSV table with one row per channel pixel. gpkg: writes a GeoPackage point layer in the raster CRS. flatgeobuf or fgb: writes a FlatGeobuf point layer in the raster CRS. |

<!-- END GENERATED CLI PARAMETERS -->
