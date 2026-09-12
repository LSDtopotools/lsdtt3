# CLI reference: lsdtt-cross-sections

Builds a channel network from user-supplied source points and extracts orthogonal cross sections along that network. Each cross section is drawn perpendicular to a locally fitted channel bearing direction. Optionally computes and samples vertical excess curvature.

`source_points_path` is required. All three output types (channel network lines, cross-section lines, and sampled points) are written by default.

## Usage

```
lsdtt-cross-sections <path> <parameter_file>
lsdtt-cross-sections <parameter_file>         # path defaults to ./
lsdtt-cross-sections -h                       # print help and exit
lsdtt-cross-sections -v                       # print version and exit
lsdtt-cross-sections -c                       # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Pipeline overview

| Step | Description |
|---|---|
| 1 | Pit filling (skipped when `raster_is_filled: true`) |
| 2 | D8 flow routing |
| 3 | Load source points and map to raster pixels |
| 4 | Build junction / channel network from source points |
| 5 | Compute vertical excess curvature (if requested) |
| 6 | Write channel network lines |
| 7 | Extract cross sections using bearing window |
| 8 | Write cross-section lines |
| 9 | Write sampled cross-section elevation points |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input raster filename (without extension). |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |

---

## Threading and routing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` uses all available hardware threads. Applies to pit filling, flow routing, and polyfit curvature calculation. |
| `boundary_conditions` | string | `oooo` | Boundary conditions for flow routing. Four-character string (N,E,S,W): `o` = open, `n` = no-flux, `p` = periodic, `b` = forced outlet. |
| `write_output_format` | string | `cog` | Raster output format for the vertical excess curvature raster. Options: `cog`, `gtiff`, `same_as_input`. |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `raster_is_filled` | bool | `false` | If true, skip pit filling and use the input raster directly for flow routing. |
| `fill_min_slope` | float | `0.0001` | Minimum gradient imposed during pit filling. Used only when `raster_is_filled: false`. |
| `sample_filled_dem` | bool | `false` | Sample cross-section elevations from the filled DEM instead of the original DEM. Routing always uses the filled DEM; by default cross sections are sampled from the original. |

---

## Source-point input

Source points define where channel sources are placed. CSV files are interpreted as latitude/longitude in EPSG:4326; GeoPackage and FlatGeobuf files are read through OGR and transformed to the raster CRS.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `source_points_path` | string | `` | **Required.** Path to the source-point input file. |
| `source_points_input_format` | string | `auto_detect` | Format hint. Options: `auto_detect`, `csv_latitude_longitude_epsg4326`, `ogr_vector_layer`. |
| `source_points_layer_name` | string | `` | Optional layer name for vector inputs. If empty, uses the only layer or the first layer. |
| `source_points_latitude_field_name` | string | `latitude` | Latitude field name for CSV inputs. |
| `source_points_longitude_field_name` | string | `longitude` | Longitude field name for CSV inputs. |

---

## Vertical excess curvature

Computed using direct polyfit. Required for the optional VEC raster output and for sampling VEC values along cross sections.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_vertical_excess_curvature_raster` | bool | `false` | Write vertical excess curvature raster. Output: `<write_fname>_vertical_excess_curvature.tif`. |
| `include_vertical_excess_curvature_in_cross_section_points` | bool | `false` | Sample vertical excess curvature at each cross-section point and add it to the point output. Adds `vertical_excess_curvature` and `vec_valid` attributes. |
| `vertical_excess_curvature_lambda_metres` | float | `50.0` | Polyfit half-wavelength (m). Controls the smoothing scale before curvature is derived. |
| `vertical_excess_curvature_z_factor` | float | `1.0` | Vertical scaling factor for polyfit derivatives. |
| `vertical_excess_curvature_flat_slope_tolerance` | float | `1e-12` | Slope tolerance below which direction-dependent metrics are nodata. |

---

## Cross-section geometry

| Parameter | Type | Default | Description |
|---|---|---|---|
| `cross_section_half_length_metres` | float | `50.0` | Half-length of each cross section (m). Each section extends this distance to either side of the channel pixel. Must be > 0. |
| `cross_section_sample_spacing_metres` | float | `5.0` | Spacing between sampled elevation points along each cross section (m). Must be > 0. |
| `bearing_window_length_metres` | float | `50.0` | Along-channel window length (m) used to estimate local channel bearing via orthogonal linear regression. Must be > 0. |
| `cross_section_thinning` | int | `1` | Keep every Nth channel pixel as a cross-section centre. `1` = every pixel. Must be > 0. |
| `skip_junction_pixels` | bool | `false` | Skip cross sections at source, confluence, or outlet junction pixels. |
| `min_bearing_window_pixels` | int | `3` | Minimum number of channel pixels required to fit a bearing. Cross sections with fewer pixels in the window are skipped. Must be ≥ 2. |

---

## Output control

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `calculate_spatial_index` | bool | `true` | Build a spatial index when writing vector files. Speeds up GIS use but requires buffering the layer in memory. |
| `write_channel_network_lines` | bool | `true` | One polyline per compact channel link. Enabled because the network locates and contextualises this workflow's cross sections. |
| `write_vector_format` | string | `flatgeobuf` | Format for channel-network and cross-section lines. Options: `flatgeobuf` / `fgb`, `gpkg`. |
| `cross_section_sample_points_write_vector_format` | string | `same_as_write_vector_format` | Optional point-output override; also accepts `csv`. |
| `write_cross_section_lines` | bool | `true` | `<write_fname>_cross_sections.<ext>` — one line per cross section. |
| `write_cross_section_sample_points` | bool | `true` | `<write_fname>_cross_section_sample_points.<ext>` — sampled elevation points with section id and optional VEC. |

---

## Example parameter file

```
read_path:  /data/dems/
read_fname: my_dem
write_path: /data/output/
write_fname: my_dem_xs

raster_is_filled:    false
fill_min_slope:      0.0001
boundary_conditions: oooo

source_points_path:         /data/sources/channel_sources.csv
source_points_input_format: csv_latitude_longitude_epsg4326

cross_section_half_length_metres:    100.0
cross_section_sample_spacing_metres: 5.0
bearing_window_length_metres:        50.0
cross_section_thinning:              5
skip_junction_pixels:                true

write_channel_network_lines:        true
write_vector_format: flatgeobuf
write_cross_section_lines:          true
write_cross_section_sample_points:  true
cross_section_sample_points_write_vector_format: gpkg
```
