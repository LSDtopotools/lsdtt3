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

### Threading / routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all available hardware threads. Applies to pit filling, flow routing, and polyfit curvature calculation. |
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four-character string for north, east, south, west. o = open, n = no-flux, p = periodic, b = forced outlet. |

### Raster output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_output_format` | `string` | `cog` | `io` | - | Raster output format. | Allowed values: cog, gtiff, same_as_input. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). CSV is available only through selectors for point/table outputs. Per-output overrides exist only where a container genuinely differs. Raster outputs are governed separately by write_output_format. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |
| `fill_flat_resolution_choice` | `string` | `inline_gradient` | `param` | - | How flats are made drainable: inline_gradient, barnes2014b_gradient or barnes2014b_float_step. | inline_gradient (default, unchanged): the Priority-Flood fill imposes fill_min_slope during the flood. The barnes2014b options make flow converge across flats (Barnes et al. 2014, doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses fill_min_slope (m/m, must be &gt; 0) as the minimum gradient; flats are raised about 3-4x more, capped below their rims, and capped flats are reported. barnes2014b_float_step is equivalent to RichDEM ResolveFlats on float32 in its first pass: it ignores fill_min_slope and its slopes depend on elevation (of order ULP(z) / cell size). With either barnes2014b option the result is always routable: cells the first pass leaves without a lower neighbour are refilled and resolved again (up to 4 rounds), and any still left are raised by a final priority-flood pass by one float32 step each; the report gives the rounds and the cells changed after the first pass. |

### DEM preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_min_slope` | `float` | `0.0001` | `param` | `sample_filled_dem` | Minimum slope (m/m) imposed by pit filling or flat resolution. | inline_gradient (default): imposed during the Priority-Flood fill, so used only when raster_is_filled is false. barnes2014b_gradient: the minimum gradient across flats, also on a raster_is_filled input. barnes2014b_float_step: ignored. |
| `sample_filled_dem` | `bool` | `false` | `switch` | - | Sample cross-section elevations from the filled DEM instead of the original DEM. | Routing always uses the routed DEM. By default the cross sections are sampled from the original DEM so that elevations reflect the observed topography. |

### Source-point input

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `source_points_path` | `string` | - | `param` | - | Path to the source-point input file. | CSV input is interpreted as latitude/longitude in EPSG:4326 using the configured latitude/longitude field names. Georeferenced vector input such as GeoPackage (.gpkg) and FlatGeobuf (.fgb) is read through OGR and transformed into the raster CRS before mapping to pixels. |
| `source_points_layer_name` | `string` | - | `param` | - | Optional layer name for source-point vector input. | Used only for georeferenced vector inputs. If empty, the point-input loader uses the only layer if present, otherwise the first layer. |
| `source_points_input_format` | `string` | `auto_detect` | `param` | - | Format hint for the source-point input. | Allowed values: auto_detect, csv_latitude_longitude_epsg4326, ogr_vector_layer. |
| `source_points_latitude_field_name` | `string` | `latitude` | `param` | - | Latitude field name for CSV source-point input. | Used only when the source-point input is CSV. |
| `source_points_longitude_field_name` | `string` | `longitude` | `param` | - | Longitude field name for CSV source-point input. | Used only when the source-point input is CSV. |

### Vertical excess curvature

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vertical_excess_curvature_raster` | `bool` | `false` | `switch` | - | Write the vertical excess curvature raster derived from the DEM. | Computed using direct polyfit followed by compute_polyfit_metrics. |
| `vertical_excess_curvature_lambda_metres` | `float` | `50.0` | `param` | `write_vertical_excess_curvature_raster` | Polyfit half-wavelength in metres used to compute vertical excess curvature. | Controls the smoothing scale of the fitted surface before curvature metrics are derived. |
| `vertical_excess_curvature_z_factor` | `float` | `1.0` | `param` | `write_vertical_excess_curvature_raster` | Vertical scaling factor used when interpreting polyfit derivatives. | - |

### General

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `include_vertical_excess_curvature_in_cross_section_points` | `bool` | `false` | - | - | Sample vertical excess curvature along each cross section and add it to the point output. | Adds vertical_excess_curvature and vec_valid attributes to each sampled cross-section point. |
| `vertical_excess_curvature_flat_slope_tolerance` | `float` | `1e-12` | - | - | Slope tolerance below which direction-dependent polyfit metrics are nodata. | Passed through to compute_polyfit_metrics. |

### Cross-section geometry

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `cross_section_half_length_metres` | `float` | `50.0` | `param` | - | Half-length of each cross section in metres. | Each cross-section line extends this distance to either side of the channel centre pixel. |
| `cross_section_sample_spacing_metres` | `float` | `5.0` | `param` | - | Spacing between sampled points along each cross section in metres. | Sample points are written with elevation attributes. |
| `bearing_window_length_metres` | `float` | `50.0` | `param` | - | Along-channel window length used to estimate local channel bearing. | A local channel segment of this total length is assembled around each channel pixel and fitted with orthogonal linear regression. The cross section is then drawn orthogonal to that fitted channel direction. |
| `cross_section_thinning` | `int` | `1` | `param` | - | Keep every Nth channel pixel as a cross-section centre. | 1 = every channel pixel. Larger values thin the output in row-major order. |
| `skip_junction_pixels` | `bool` | `false` | `param` | - | Skip cross sections centred on junction pixels. | If true, no cross sections are generated at source, confluence, or outlet junctions. |
| `min_bearing_window_pixels` | `int` | `3` | `param` | - | Minimum number of channel pixels required to fit a local bearing. | Cross sections are skipped where the bearing window contains too few channel pixels. |

### Output control

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `calculate_spatial_index` | `bool` | `true` | `io` | - | Whether or not to calculate a spatial index for vector data writes. | A spatial index can improve GIS performance but requires buffering the whole layer in memory. |
| `write_channel_network_lines` | `bool` | `true` | `switch` | - | Write the extracted channel network as vector lines. | Enabled by default because the network locates the cross sections and is part of this workflow's interpretation, not an optional extraction product. Each compact channel link between junctions is written as a line feature. |
| `write_cross_section_lines` | `bool` | `true` | `switch` | - | Write simple cross-section location lines. | Writes one line per cross section. This is the light-weight geometry output describing section position and section id. |
| `write_cross_section_sample_points` | `bool` | `true` | `switch` | - | Write sampled cross-section elevations as point data. | Each sampled point stores section id, location, sampled elevation, and optional sampled vertical excess curvature. |
| `cross_section_sample_points_write_vector_format` | `string` | `same_as_write_vector_format` | `param` | - | Output format for sampled cross-section points. | Allowed values: csv, flatgeobuf, fgb, gpkg. |

<!-- END GENERATED CLI PARAMETERS -->
