# CLI reference: lsdtt-water-centreline

Water/floodplain centreline extraction and tracing. Combines multiple methods to identify river centrelines and floodplain axes. Outputs centreline polylines and optional diagnostic rasters.

## Usage

```
lsdtt-water-centreline  <path> <parameter_file>
lsdtt-water-centreline  <parameter_file>        # path defaults to ./
lsdtt-water-centreline  -h                      # print help and exit
lsdtt-water-centreline  -v                      # print version and exit
lsdtt-water-centreline  -c                      # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

Keys use `snake_case`. All keys and their defaults are listed below.

---

!!! note "Multipart lines"
    `input_vector` accepts multipart geometry (`MultiLineString` and friends). Each feature is split into its parts on read, with the attributes copied onto each, and the split is reported in the log.

    Shapefiles are the commonest source of multipart data, and **versions before August 2026 discarded it silently** — see [Vector data](../what_is_lsdtt3/vector_data.md).

## Pipeline overview

Steps run in this order:

| Step | Description |
|---|---|
| 1. DEM preparation | Load DEM, verify projection, optionally fill pits |
| 2. Water mask detection | Identify water-related low-relief areas |
| 3. Skeletonisation | Extract 1-pixel centreline(s) |
| 4. Vectorisation | Convert centreline raster to polyline(s) |
| 5. Smoothing *(optional)* | Apply Douglas-Peucker or Visvalingam-Whyatt simplification |
| 6. Output | Write centreline polylines and diagnostics |

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

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `raster_is_filled` | bool | `false` | Set `true` if DEM is already pit-filled. |
| `fill_min_slope` | float | `0.0001` | Minimum slope (m/m) imposed during pit filling. |

---

## Water mask detection

| Parameter | Type | Default | Description |
|---|---|---|---|
| `detection_method` | string | `relief` | Method to identify water areas: `relief` (local relief), `slope` (slope-based), `curvature` (curvature-based). |
| `relief_threshold_metres` | float | `5.0` | Relief threshold (metres) for water mask. Pixels with relief < threshold are water. |
| `slope_threshold_degrees` | float | `3.0` | Slope threshold (degrees). Pixels below this are potential water. |
| `smoothing_scale_metres` | float | `100.0` | Smoothing scale (metres) for relief/curvature computation. |

---

## Skeleton extraction

| Parameter | Type | Default | Description |
|---|---|---|---|
| `skeleton_method` | string | `zhang_suen` | Skeletonisation algorithm: `zhang_suen`, `thinning`. |
| `prune_skeleton_length_pixels` | int | `5` | Remove skeleton spurs shorter than this (pixels). |
| `bridge_small_gaps_pixels` | int | `3` | Bridge isolated skeleton segments separated by gaps < this (pixels). |

---

## Centreline smoothing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `apply_smoothing` | bool | `false` | Simplify centreline using Douglas-Peucker or similar. |
| `smoothing_tolerance_metres` | float | `10.0` | Tolerance (metres) for centreline simplification. Larger = simpler line. |
| `smoothing_algorithm` | string | `douglas_peucker` | Simplification algorithm: `douglas_peucker`, `visvalingam_whyatt`. |

---

## Output: Centreline vector

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_centreline_vector` | bool | `true` | Write centreline(s) as polyline vector layer. |
| `centreline_vector_format` | string | `flatgeobuf` | Output format: `flatgeobuf` / `fgb` (preferred), `gpkg`. |

Output file: `<write_fname>_water_centreline.<ext>`

---

## Output: Diagnostic rasters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_water_mask` | bool | `false` | Write binary water mask raster. |
| `write_skeleton_raster` | bool | `false` | Write skeleton/centreline as raster. |
| `write_relief_raster` | bool | `false` | Write relief raster used for detection. |

Output files:
- `<write_fname>_water_mask.tif`
- `<write_fname>_skeleton.tif`
- `<write_fname>_relief.tif`

---

## Output vector format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_vector_format` | string | `flatgeobuf` | Output **vector** format for the centreline layer: `flatgeobuf`/`fgb`, `gpkg`. This is the vector key, not the raster key `write_output_format` used by the raster tools. |

`write_vector_format` was called `write_output_format` in earlier releases. The
old name still works as a deprecated alias: the value is carried across and a
warning naming the replacement is printed. Update your parameter files, because
`write_output_format` means a *raster* format in every other CLI.

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

### -

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `input_vector` | `string` | - | `param` | `clip_lines_to_areas` | Input vector dataset (e.g. a GeoPackage) holding both layers. | Required. Resolved against read_path if not absolute. |
| `water_area_layer` | `string` | - | `param` | `clip_lines_to_areas` | Name of the water-area POLYGON layer. | Required. e.g. vml-nr73sw_Water_Area. |
| `water_line_layer` | `string` | - | `param` | `clip_lines_to_areas` | Name of the water LINESTRING (tributary) layer. | Required. e.g. vml-nr73sw_Water_Line. |
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Output VECTOR format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. Distinct from the raster key write_output_format used by the raster tools: this CLI writes vectors only. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). GeoPackage is an SQLite database, and SQLite cannot take the file locks it needs on a drive mounted into WSL2, so a gpkg write there fails. |
| `epsg` | `int` | `27700` | `param` | `clip_lines_to_areas` | EPSG code for the synthesised grids and output CRS. | Must match the input layers' projected CRS (metres). |
| `cell_size_metres` | `float` | `1.0` | `param` | `clip_lines_to_areas` | Synthesised rasterisation grid resolution (m). | Finer resolves narrow necks but costs memory; auto-coarsens above a ~50M-cell cap per polygon. |
| `padding_metres` | `float` | `5.0` | `param` | `clip_lines_to_areas` | Padding (m) added around each polygon bbox in the synthesised grid. | - |
| `preserve_distance_pixels` | `float` | `1.0` | `param` | `clip_lines_to_areas` | Centreline preference for thin_distance_guided (pixels). | Higher protects the medial line more on blocky banks. |
| `max_spur_pixels` | `int` | `0` | `param` | `clip_lines_to_areas` | Maximum skeleton spur length pruned (pixels); 0 = skip (recommended). | The centreline is the skeleton's longest path, which already ignores spurs; pruning only trims the centreline ends back from the polygon tips, so leave at 0 unless you have a specific reason. |
| `max_loop_pixels` | `int` | `10` | `param` | `clip_lines_to_areas` | Maximum terminal loop ('lollipop') length pruned (pixels); 0 = skip. | - |
| `detect_channel_spurs` | `bool` | `true` | `param` | `clip_lines_to_areas` | Beyond the main long-axis line, follow narrow channel arms of the water body. | Uses the local half-width (distance to bank) from the EDT: a skeleton branch whose tip is narrow (a river channel) is retained; wide bays are ignored. Helps trace a channel running off a lake. |
| `channel_max_half_width_metres` | `float` | `12.0` | `param` | `clip_lines_to_areas` | Max half-width (m) at a branch tip for it to count as a channel spur. | ~half the widest river channel you want followed (full width ~2x this). Below this = channel; above = lake/bay. |
| `min_channel_spur_len_metres` | `float` | `20.0` | `param` | `clip_lines_to_areas` | Minimum length (m) of a channel spur to retain. | Filters short nubs; raise to keep only substantial channel arms. |
| `centreline_smooth_iterations` | `int` | `3` | `param` | `clip_lines_to_areas` | Laplacian smoothing passes applied to the centreline before DP. | Removes the pixel-skeleton staircase wobble (endpoints fixed). 0 = off. |
| `dp_tolerance_metres` | `float` | `1.0` | `param` | `clip_lines_to_areas` | Douglas-Peucker tolerance (m) for the output centreline. | Thins vertices after smoothing; ~1 cell size is a good default. |
| `tributary_snap_tolerance_metres` | `float` | `2.0` | `param` | `clip_lines_to_areas` | Max distance (m) a water-line endpoint may sit from the polygon boundary to count as a tributary contact. | - |
| `centreline_snap_max_distance_metres` | `float` | `100.0` | `param` | `clip_lines_to_areas` | Max connector (spur) length (m) from a tributary mouth to the centreline. | Mouth-to-centreline distance ~ the local half-width, so wide water bodies need a larger value than narrow rivers. |
| `min_centreline_component_metres` | `float` | `10.0` | `param` | `clip_lines_to_areas` | Drop centreline component pieces shorter than this (m). | Declutters the many tiny components of a fragmented lake skeleton. The per-polygon long-axis fallback is exempt, so every polygon keeps a line. |
| `max_centreline_aspect_for_loch` | `float` | `2.0` | `param` | `clip_lines_to_areas` | Feret aspect (d_max/d_min) below which a polygon is treated as a loch. | Lochs have no meaningful single long axis. |
| `fill_holes` | `bool` | `true` | `param` | `clip_lines_to_areas` | Fill interior holes/islands before skeletonising. | Avoids skeleton loops; disable to keep braided islands. |
| `skip_loch_centreline` | `bool` | `false` | `param` | `clip_lines_to_areas` | For lochs, skip the centreline and star-connect tributaries to the centroid. | When false, a loch still gets its long-axis path, flagged as a loch. |
| `clip_lines_to_areas` | `bool` | `true` | `switch` | - | Remove Water_Line portions that are the OUTLINE of a water area. | OS VectorMap Local draws water-body banks in the line layer. Segments whose midpoint is inside a Water_Area or within line_overlap_tolerance of its boundary are dropped; the outside runs are kept as tributaries. |
| `line_overlap_tolerance_metres` | `float` | `2.0` | `param` | `clip_lines_to_areas` | Distance (m) within which a line segment counts as a water-area outline. | Raise if bank lines are digitised a little off the polygon edge. |
| `min_tributary_length_metres` | `float` | `0.0` | `param` | `clip_lines_to_areas` | Drop clipped tributary runs shorter than this (m); 0 = keep all. | Use to remove tiny stub fragments left at water-body banks. |
| `emit_unchanged_tributaries` | `bool` | `true` | `param` | `clip_lines_to_areas` | Copy the CLIPPED tributary lines into the output as role='tributary'. | Gives a single connected network in one layer (outlines removed). |
| `calculate_spatial_index` | `bool` | `true` | `io` | - | Build a spatial index when writing the output layer. | - |

<!-- END GENERATED CLI PARAMETERS -->
