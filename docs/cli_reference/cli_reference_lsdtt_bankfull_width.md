# CLI reference: lsdtt-bankfull-width

Estimates bankfull channel width and depth from a high-resolution DEM with Manning's equation. It builds a single-source channel network from user-supplied source points, extracts cross sections along it, measures the along-channel slope and drainage area at each section, and then, per section, finds the water depth whose Manning discharge matches a target discharge.

`source_points_path` is required. Manning widths need `run_manning_widths: true` and a reach table (`reach_manning_parameters_path`) giving Manning's n and the target discharge for each stretch of channel. The tool does not derive either.

## Usage

```
lsdtt-bankfull-width <path> <parameter_file>
lsdtt-bankfull-width <parameter_file>         # path defaults to ./
lsdtt-bankfull-width -h                       # print help and exit
lsdtt-bankfull-width -v                       # print version and exit
lsdtt-bankfull-width -c                       # print citation and exit
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
| 3 | Load source points and map them to raster pixels |
| 4 | Build the channel network from the source points |
| 5 | Extract cross sections using a bearing window |
| 6 | Channel profile per cross section: distance from source, drainage area, thalweg elevation, along-channel slope |
| 7 | Sampled elevation profile per cross section |
| 8 | Manning width / depth solve per cross section, using the reach that contains it |

Every output comes from the same cross-section extraction, so `cross_id` / `cross_section_id` numbering matches across all files from one run.

---

## Inputs

**Reach table** (`reach_manning_parameters_path`, CSV). One row per reach:

```
start_distance_m,end_distance_m,manning_n,target_discharge_m3s
0,150,0.035,5
150,1000,0.04,8
```

A cross section belongs to a reach by its distance from the source, half-open: `start <= d < end`. `cross_section_id` is assigned by a row-major raster scan and has nothing to do with position along the channel, so it is not used for matching. A malformed row is reported with its line number.

**Channel profile** (`channel_profile_csv_path`, CSV). Read only when the profile is not computed in the same run (that is, when `write_channel_profile_csv` is false). It must come from the same DEM and cross-section geometry parameters.

---

## Outputs

| File | Switch | Contents |
|---|---|---|
| `<write_fname>_channel_profile.csv` | `write_channel_profile_csv` | `cross_id`, `distance_from_source_m`, `drainage_area_m2`, `thalweg_elevation_m`, `slope_m_per_m`, `slope_valid` |
| `<write_fname>_cross_section_sample_points.<ext>` | `write_cross_section_sample_points` | Sampled elevation profile per section: `cross_section_id`, `offset_metres`, `x`, `y`, `elevation`, `valid` |
| `<write_fname>_manning_bank_points.<ext>` | `write_manning_bank_points` | One point per bank per valid section |
| `<write_fname>_manning_widths.csv` | `write_manning_widths_csv` | Per section: solved depth, water-surface elevation, bank positions, width, flow area, wetted perimeter, hydraulic radius, discharge at the solution, `valid`, `quality_flag` |

Slopes are dimensionless gradients (m/m). Both point layers are written in `write_vector_format` (`flatgeobuf` / `fgb` or `gpkg`). CSV is not available for them, because a point CSV would have to carry EPSG:4326 latitude/longitude; use `write_manning_widths_csv` for tabular output.

### Renamed parameters

These older names are still read and print a deprecation warning:

| Old name | Use instead |
|---|---|
| `cross_section_sample_points_output_format` | `write_vector_format` (an old value is still applied to that layer for the run) |
| `manning_bank_points_output_format` | `write_vector_format` (an old value is still applied to that layer for the run) |
| `minimum_bearing_window_pixels` | `min_bearing_window_pixels` |

---

## Example parameter file

```
read_path:  /data/dems/
read_fname: my_dem
write_path: /data/output/
write_fname: bf

n_threads: 1
source_points_path:            /data/sources/channel_sources.csv
reach_manning_parameters_path: /data/reaches.csv

cross_section_half_length_metres:    160
cross_section_sample_spacing_metres: 2
bearing_window_length_metres:        20
cross_section_thinning:              10
thalweg_search_half_width_metres:    6
slope_window_half_length_metres:     30

run_manning_widths:       true
manning_max_depth_metres: 10
write_vector_format:      flatgeobuf
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

### Flow routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four-character string for north, east, south, west. o = open, n = no-flux, p = periodic, b = forced outlet. |
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all available hardware threads. A positive value is honoured verbatim and OVERRIDES OMP_NUM_THREADS, so setting 1 here runs serial on any machine. |

### DEM preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `vertical_units` | `string` | `metres` | `param` | - | Vertical unit of the input DEM elevations. | metres / feet / us_survey_feet. Slope, and therefore Manning discharge, is wrong by a factor of 3.28 if the DEM is in feet and this is left as metres. Nothing in the file records the unit. Horizontal units must already be metres. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope (m/m) imposed by pit filling or flat resolution. | inline_gradient (default): imposed during the Priority-Flood fill across filled depressions and flats, so used only when raster_is_filled is false. barnes2014b_gradient: the minimum gradient across flats, also on a raster_is_filled input. barnes2014b_float_step: ignored. |
| `fill_flat_resolution_choice` | `string` | `inline_gradient` | `param` | - | How flats are made drainable: inline_gradient, barnes2014b_gradient or barnes2014b_float_step. | inline_gradient (default, unchanged): the Priority-Flood fill imposes fill_min_slope during the flood. The barnes2014b options make flow converge across flats (Barnes et al. 2014, doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses fill_min_slope (m/m, must be &gt; 0) as the minimum gradient; flats are raised about 3-4x more, capped below their rims, and capped flats are reported. barnes2014b_float_step is equivalent to RichDEM ResolveFlats on float32 in its first pass: it ignores fill_min_slope and its slopes depend on elevation (of order ULP(z) / cell size). With either barnes2014b option the result is always routable: cells the first pass leaves without a lower neighbour are refilled and resolved again (up to 4 rounds), and any still left are raised by a final priority-flood pass by one float32 step each; the report gives the rounds and the cells changed after the first pass. |

### Source points

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `source_points_path` | `string` | - | `io` | - | Path to the source-point input file (exactly one point). | CSV: latitude/longitude in EPSG:4326. Vector (.gpkg/.fgb): read through OGR and transformed to the raster CRS. Required. When write_channel_profile_csv is true it must resolve to a single source (a single mainstem channel). |
| `source_points_layer_name` | `string` | - | `param` | - | Optional layer name for vector source-point input. | Only used for vector input. Empty means the first layer. |
| `source_points_input_format` | `string` | `auto_detect` | `param` | - | Source-point file format. | auto_detect / csv_latitude_longitude_epsg4326 / ogr_vector_layer. |
| `source_points_latitude_field_name` | `string` | `latitude` | `param` | - | Latitude column name for CSV source-point input. | Only used for CSV input. |
| `source_points_longitude_field_name` | `string` | `longitude` | `param` | - | Longitude column name for CSV source-point input. | Only used for CSV input. |

### Cross sections

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `cross_section_half_length_metres` | `float` | `50.0` | `param` | - | Half-length of each cross section [m]. | Must be greater than cross_section_sample_spacing_metres and comfortably larger than the expected bank-to-bank half-width, or banks will not be found (quality_flag bank_crossing_not_found). |
| `cross_section_sample_spacing_metres` | `float` | `5.0` | `param` | - | Spacing between sampled points along each cross section [m]. | Must be &gt; 0. Coarse spacing under-resolves the channel; width accuracy is limited to roughly this value. Use a few times the DEM cell size or finer. |
| `bearing_window_length_metres` | `float` | `50.0` | `param` | - | Along-channel window used to estimate local channel bearing [m]. | Sections are cut perpendicular to this bearing. Longer windows smooth sinuous channels but blur tight bends. |
| `cross_section_thinning` | `int` | `1` | `param` | - | Keep every Nth channel pixel as a cross-section centre. | Must be &gt;= 1. 1 keeps every channel pixel. |
| `skip_junction_pixels` | `bool` | `false` | `switch` | - | Skip cross sections centred on junction pixels. | Irrelevant for a single-source network, which has no junctions. |
| `min_bearing_window_pixels` | `int` | `3` | `param` | - | Minimum channel pixels required to fit a local bearing. | Sections with fewer pixels in the bearing window are dropped. |
| `thalweg_search_half_width_metres` | `float` | `25.0` | `param` | - | Half-width around each section centre searched for the thalweg [m]. | The thalweg is the lowest sample within +/- this distance of the section centre, so a ditch, gully or second channel elsewhere in the section cannot capture the solve. The centre is a D8 channel pixel, which can sit up to about half a channel width from the true thalweg: set this to roughly half the channel width or a little more, and check thalweg_s in the widths CSV. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `calculate_spatial_index` | `bool` | `true` | `switch` | - | Calculate a spatial index for vector writes. | Faster loading and queries in a GIS, but the whole layer is held in memory while writing. |
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). CSV is available only through selectors for point/table outputs. Per-output overrides exist only where a container genuinely differs. Raster outputs are governed separately by write_output_format. |

### Channel profile

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_channel_profile_csv` | `bool` | `true` | `switch` | - | Write distance from source, drainage area, thalweg elevation and slope per cross section. | Requires source_points_path to resolve to exactly one source point. Output: &lt;write_fname&gt;_channel_profile.csv. Its distance and slope columns drive the Manning solve. If run_manning_widths is also set in the same run, the in-memory profile is used and channel_profile_csv_path is ignored. |
| `slope_window_half_length_metres` | `float` | `25.0` | `param` | `write_channel_profile_csv` | Half-length of the along-channel slope-fit window [m]. | Slope is the least-squares slope of pit-filled elevation against along-channel distance over +/- this distance, so it has no cross-channel component. Must be at least two cell sizes; use roughly 5-10 channel widths. Near the source or outlet a one-sided window is accepted only if the valid cells span at least this distance, otherwise slope_valid is 0. |
| `write_cross_section_sample_points` | `bool` | `true` | `switch` | - | Write the full sampled elevation profile per cross section. | Fields: cross_section_id, offset_metres, x, y, elevation, valid. Output: &lt;write_fname&gt;_cross_section_sample_points.fgb (or .gpkg). |

### Manning widths

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_manning_widths` | `bool` | `false` | `switch` | - | Compute bankfull channel widths via Manning's equation. | Requires reach_manning_parameters_path (manning_n and target_discharge_m3s per reach, bounded by DISTANCE FROM SOURCE, not cross_id). For each cross section inside a reach, solves the depth that reproduces the reach's target discharge using THAT cross section's own slope: discharge and n are constant across a reach, but slope and geometry vary per section, so solved depth varies along the reach. Only water connected to the channel is counted. |
| `reach_manning_parameters_path` | `string` | - | `io` | `run_manning_widths` | CSV of reach parameters. | Columns: start_distance_m,end_distance_m,manning_n,target_discharge_m3s. Bounds are along-channel distance from source in metres. Reaches are half-open [start, end), must not overlap, and must have start &lt; end, manning_n &gt; 0 and target_discharge_m3s &gt; 0. Slope comes from the channel profile, not this file. |
| `channel_profile_csv_path` | `string` | - | `io` | `run_manning_widths` | Channel-profile CSV supplying per-section distance and slope. | Only read when write_channel_profile_csv is false. Must be the _channel_profile.csv from an earlier run against the SAME DEM and SAME cross-section geometry parameters; otherwise cross_ids do not correspond and results are wrong or empty (a coverage report is printed). |
| `manning_min_depth_metres` | `float` | `0.01` | `param` | `run_manning_widths` | Lower bracket of the depth solve [m]. | Must be &gt; 0 and &lt; manning_max_depth_metres. Sections whose target discharge lies outside the bracket are flagged target_discharge_not_bracketed. |
| `manning_max_depth_metres` | `float` | `15.0` | `param` | `run_manning_widths` | Upper bracket of the depth solve [m]. | Raise for very large rivers. Must exceed manning_min_depth_metres. |
| `manning_depth_tolerance_metres` | `float` | `0.005` | `param` | `run_manning_widths` | Bisection depth tolerance [m]. | Must be &gt; 0. Solved depth is within about this of the exact root. |
| `manning_max_discharge_relative_error` | `float` | `0.05` | `param` | `run_manning_widths` | Maximum \|Q_solved - Q_target\| / Q_target before a section is flagged. | Above this the section gets quality_flag spill_discontinuity: the water surface overtopped a bank crest into lower ground and no depth reproduces the target discharge. Must be &gt; 0. |
| `write_manning_bank_points` | `bool` | `true` | `switch` | `run_manning_widths` | Write Manning bank points as georeferenced points. | One point per bank per valid section. Output: &lt;write_fname&gt;_manning_bank_points.fgb (or .gpkg). |
| `write_manning_widths_csv` | `bool` | `true` | `switch` | `run_manning_widths` | Write per-section Manning results to CSV. | One row per section, valid or not, with quality_flag. Output: &lt;write_fname&gt;_manning_widths.csv. |

<!-- END GENERATED CLI PARAMETERS -->
