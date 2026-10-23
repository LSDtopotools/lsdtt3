# CLI reference: lsdtt-water-centreline

Turns mapped water bodies into a single connected line network. For each
water-**area** polygon (a wide river reach or a loch) it extracts a
centreline, then connects the surrounding water-**line** tributaries to it.
This is a pure vector tool: **no DEM is read.** It was written for Ordnance
Survey VectorMap Local, which supplies both layers in one GeoPackage, and it
reads that product's `Fid` and `FeatCode` attributes when they are present.

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

Keys use `snake_case`. All keys and their defaults are listed under
[Parameters](#parameters) below.

---

## Pipeline

1. Read the polygon layer (`water_area_layer`) and the line layer
   (`water_line_layer`) from `input_vector_fname`.
2. With `clip_lines_to_areas: true` (default), drop the parts of the water
   lines that are the outlines of water areas (VectorMap Local draws banks in
   the line layer). What is left outside the areas are the tributaries.
3. Rasterise each polygon on its own grid of `cell_size_metres` (default 1 m),
   fill its holes (`fill_holes`), and skeletonise it. The centreline is the
   skeleton's longest path, smoothed (`centreline_smooth_iterations`) and
   simplified (`dp_tolerance_metres`). With `detect_channel_spurs: true`,
   narrow channel arms (half-width below `channel_max_half_width_metres`) are
   followed as well.
4. A polygon whose Feret aspect ratio is below
   `max_centreline_aspect_for_loch` is a loch. It still gets its long-axis
   line, unless `skip_loch_centreline: true`, in which case its tributaries are
   joined to its centroid instead.
5. Join each tributary mouth (an endpoint within
   `tributary_snap_tolerance_metres` of a polygon boundary) to the nearest
   centreline with a connector, up to `centreline_snap_max_distance_metres`
   long.

## Example

```
read_path: /home/alice/os_vml/
input_vector_fname: vml-nr73sw.gpkg
water_area_layer: vml-nr73sw_Water_Area
water_line_layer: vml-nr73sw_Water_Line
epsg: 27700
```

`epsg` (default 27700, British National Grid) must match the CRS of the input
layers, which must be projected with metre units. The layers are not
reprojected.

## Output

`<write_fname>_water_centreline.fgb` (or `.gpkg` with
`write_vector_format: gpkg`): one line layer holding the whole network, with
attributes

| Field | Meaning |
|---|---|
| `role` | `centreline`, `connector` (tributary mouth to centreline) or `tributary` |
| `source_fid` | `Fid` of the source water-area polygon; empty for tributaries |
| `feat_code` | `FeatCode` of the source feature |
| `length_m` | Line length (m) |
| `segment_id` | Running id of the line |

Tributaries are copied into the output only with
`emit_unchanged_tributaries: true` (default).

`write_vector_format` was called `write_output_format` in earlier releases. The
old name still works as a deprecated alias: the value is carried across and a
warning naming the replacement is printed. Update your parameter files, because
`write_output_format` means a *raster* format in every other CLI.

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

### -

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `input_vector_fname` | `string` | - | `param` | `clip_lines_to_areas` | Input vector dataset (e.g. a GeoPackage) holding both layers. | Required. A relative path is read from read_path. |
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
