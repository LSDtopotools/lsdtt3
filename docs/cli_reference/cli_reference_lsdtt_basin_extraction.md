# lsdtt3 basin extraction guide

This document explains how to select and extract basins with
`lsdtt-basin-extraction`. Developer documentation for the underlying basin
objects is in [Basin object architecture](../developer_notes/basin_object_architecture.md).

---

!!! note "Multipart lines and points"
    `crossing_line_fname`, `single_channel_points_fname`, `side_tributary_points_fname` accept multipart geometry (`MultiLineString` and friends). Each feature is split into its parts on read, with the attributes copied onto each, and the split is reported in the log.

    Shapefiles are the commonest source of multipart data, and **versions before August 2026 discarded it silently** — see [Vector data](../what_is_lsdtt3/vector_data.md).
## Selecting all basins of one stream order

`select_all_basins_of_stream_order` is a basin-selection mode that extracts all
basins whose associated junction stream order is exactly a user-specified value.
Use this when you want a complete basin set at one hierarchy level (for example,
all 3rd-order basins) instead of a stream-order range filter.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `select_all_basins_of_stream_order` | bool | `false` | Enable exact stream-order basin selection. Mutually exclusive with other `select_*` modes. |
| `selected_stream_order` | int | `-1` | Target stream order. Required and must be `> 0` when `select_all_basins_of_stream_order: true`. |

Example parameter file:

```

## Selecting zero-order basins

Zero-order basins are selected separately from junction basins. They are the
complete upslope contributing areas of the channel-head source pixels returned
by `threshold_contributing_pixels`; they are not Strahler-order-1 junction
basins.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `select_all_zero_order_basins` | bool | `false` | Select one zero-order basin for every channel-head source pixel. |

Example:

```
select_all_zero_order_basins: true
write_basin_outlines: true
basin_outlines_write_vector_format: flatgeobuf
```

The output is `<write_fname>_zero_order_basin_outlines.fgb` with `zero_order_id`
and `source_flat` attributes.
read_path:   /data/nujang/
write_path:  /data/nujang/
read_fname:  nujang_COP30_fill.tif
write_fname: nujang_order3
raster_is_filled: true

select_all_basins_of_stream_order: true
selected_stream_order: 3

write_basin_outlines_as_polygon: true
```

## Selecting basins that cross a line (e.g. a fault)

`select_basins_from_line_crossings` is a basin-selection mode for studying basin
geometry relative to a baseline structure such as a fault or a valley axis. It
finds the channel pixels that lie on a supplied polyline and whose contributing
area exceeds a threshold, treats each as a basin outlet, and extracts the basin
upstream of it. Because the analysis is stream-order / geometry based, it
automatically keeps only the **largest complete basin per independent crossing**
(it forces the incomplete- and nested-basin filters on).

The supplied line need not be in the DEM CRS — it is reprojected into the DEM
CRS internally (its layer CRS is read and a coordinate transform applied), so
you can pass, for example, a WGS84 fault trace straight from the GEM database or
from `lsdtt-raster-preprocessing`'s `clip_vector_data_to_raster` output.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `select_basins_from_line_crossings` | bool | `false` | Enable the mode. Mutually exclusive with the other `select_*` modes. |
| `crossing_line_fname` | string | `NULL` | Polyline vector the channels must cross (any OGR format; any CRS). Required. |
| `crossing_minimum_contributing_pixels` | int | `500` | Minimum contributing pixels at a crossing to seed a basin outlet. |

Combine with `filter_basins_by_area` + `basins_max_area_pixels` to drop a trunk
valley that the line runs *along* (rather than across), keeping only the
fault-transverse catchments.

### Per-basin analysis output

In this mode the tool writes `<write_fname>_line_crossing_basins.csv` — one row
per retained basin, with every per-basin quantity computed in C++ (outlet, centre
of mass, area, the outlet→COM bearing, and the bearing relative to the local
fault orientation):

| Column | Description |
|---|---|
| `basin_id` | Index of the basin in the retained set |
| `area_m2` | Basin area (upslope of the crossing), m² |
| `outlet_x`, `outlet_y` | Crossing-pixel (outlet) cell-centre, raster CRS |
| `com_x`, `com_y` | Basin centre of mass, raster CRS |
| `outlet_to_com_bearing_deg` | Bearing from outlet to COM (deg from North, clockwise) |
| `local_fault_azimuth_deg` | Azimuth of the nearest line segment to the basin |
| `side` | Which side of the line the basin lies on (`east`/`west`, by COM easting) |
| `deflection_from_fault_orthogonal_deg` | Outlet→COM bearing minus the line-orthogonal pointing into the basin, in (−180, 180]. ~0 means the basin extends straight across the line; a systematic sign by side indicates shear-related bending. |

The polygons (`write_basin_outlines_as_polygon`) and this CSV are consumed by the
`lsdviztools3` `fault_basin_deformation.py` example for plotting (rose diagrams of
basin bearings, deflection-by-side distributions, etc.).

Example parameter file:

```
read_path:   /data/parkfield/
write_path:  /data/parkfield/
read_fname:  parkfield_COP30_fill.tif
write_fname: parkfield
raster_is_filled: true
threshold_contributing_pixels: 200

select_basins_from_line_crossings: true
crossing_line_fname: parkfield_COP30_clipped_vector.gpkg
crossing_minimum_contributing_pixels: 500

# focus on fault-transverse catchments (drop the along-fault trunk valley)
filter_basins_by_area: true
basins_max_area_pixels: 22222

write_basin_outlines_as_polygon: true
```

## Basins draining into a single channel

`extract_basins_from_single_channel` builds **one main channel** from a supplied
upstream point and extracts the side-tributary basins draining into it. It is the
single-trunk counterpart of `extract_side_tributary_basins_from_points`, but
identifies the channel endpoints explicitly (no reliance on input feature order)
and snaps the downstream limit onto the traced flow path.

The points file holds one or two points (any OGR format; reprojected to the DEM
CRS on read):

- **one point** — trace D8 from it to the DEM edge; that path is the single
  channel.
- **two points** — the **upstream** point is identified by D8 reachability: both
  points' D8 paths are traced and the upstream point is the one whose reach to
  the *common-downstream point* (where the two flows merge) is longer. The single
  channel runs from the upstream point to that confluence (so a downstream point
  that sits on a hillslope is snapped to where its flow joins the channel). If the
  two points are in different sub-basins (their flows never meet), a
  `downstream point not found` warning is printed and the longer of the two D8
  paths is used.

Side basins are then found exactly as for `extract_side_tributary_basins_from_points`
— a side outlet is the pixel just before a tributary joins the channel, kept if
its basin is complete (`reject_incomplete_side_tributary_basins`) and at least
`side_tributary_minimum_contributing_pixels` in size.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `extract_basins_from_single_channel` | bool | `false` | Enable the mode. |
| `single_channel_points_fname` | string | `NULL` | Point vector with 1 or 2 points (upstream, optional downstream). Required. |
| `side_tributary_minimum_contributing_pixels` | int | `1000` | Minimum side-basin size (pixels). |
| `reject_incomplete_side_tributary_basins` | bool | `true` | Drop basins touching nodata / the DEM edge. |

Outputs (shared with the side-tributary writers, selected by the
`write_side_tributary_*` flags and `side_tributary_write_vector_format`):
`<write_fname>_side_tributary_main_channel`, `..._side_tributary_outlets`,
`..._side_tributary_basin_outlines`.

The `lsdviztools3` `single_channel_basins.py` example drives this mode and plots
the single channel with the basins draining into it.

Example parameter file:

```
read_path:   /data/nujang/
write_path:  /data/nujang/
read_fname:  nujang_COP30_fill.tif
write_fname: nujang_sc
raster_is_filled: true

extract_basins_from_single_channel: true
single_channel_points_fname: nujang_mainstem_points.gpkg
side_tributary_minimum_contributing_pixels: 2000
reject_incomplete_side_tributary_basins: true

write_side_tributary_main_channel: true
write_side_tributary_outlets: true
write_side_tributary_basin_outlines: true
side_tributary_write_vector_format: gpkg
```

<!-- BEGIN GENERATED CLI PARAMETERS -->
## Parameters

This section is generated from the executable's `-h json` registry. Do not edit it manually.

**Status: stable.**

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
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all available hardware threads. Applies to pit filling and flow-routing construction. |
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four-character string for north, east, south, west. o = open, n = no-flux, p = periodic, b = forced outlet. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |
| `fill_flat_resolution_choice` | `string` | `inline_gradient` | `param` | - | How flats are made drainable: inline_gradient, barnes2014b_gradient or barnes2014b_float_step. | inline_gradient (default, unchanged): the Priority-Flood fill imposes fill_min_slope during the flood. The barnes2014b options make flow converge across flats (Barnes et al. 2014, doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses fill_min_slope (m/m, must be &gt; 0) as the minimum gradient; flats are raised about 3-4x more, capped below their rims, and capped flats are reported. barnes2014b_float_step is equivalent to RichDEM ResolveFlats on float32 in its first pass: it ignores fill_min_slope and its slopes depend on elevation (of order ULP(z) / cell size). With either barnes2014b option the result is always routable: cells the first pass leaves without a lower neighbour are refilled and resolved again (up to 4 rounds), and any still left are raised by a final priority-flood pass by one float32 step each; the report gives the rounds and the cells changed after the first pass. |

### DEM preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope (m/m) imposed by pit filling or flat resolution. | inline_gradient (default): imposed during the Priority-Flood fill, so used only when raster_is_filled is false. barnes2014b_gradient: the minimum gradient across flats, also on a raster_is_filled input. barnes2014b_float_step: ignored. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). CSV is available only through selectors for point/table outputs. Per-output overrides exist only where a container genuinely differs. Raster outputs are governed separately by write_output_format. |
| `calculate_spatial_index` | `bool` | `true` | `io` | - | Whether or not to calculate a spatial index for vector data writes. | A spatial index helps a GIS to speed up spatial queries, so it is often helpful to have a spatial index calculated. But calculating this requires the entire vector layer to sit in memory, preventing feature-by-feature writes and adding a lot of overhead. In some cases (e.g. large datasets) it may therefore be preferred to generate vectors without spatial index - these can always be calculated later on using e.g. gdal. |
| `write_channel_network_points` | `bool` | `false` | `switch` | - | Write the channel network as point data, where each channel pixel is a feature with attribute values, such as drainage_area. | Compared to write_channel_network_lines, this write will be larger because there are more features, but it does record the drainage area at each channel pixel. |
| `write_channel_network_lines` | `bool` | `false` | `switch` | - | Write the channel network as vector lines, where each channel link (defined between source or upstream junction and outlet or downstream junction) is a feature. | Does not support csv writes, only gpkg and flatgeobuf/fgb. It is off by default because basin extraction does not require a channel-line deliverable. Compared to write_channel_network as points, the lines filesize will be smaller (there are fewer features), at the cost of not recording drainage area at each channel pixel. |
| `write_channel_network_junctions` | `bool` | `false` | `switch` | - | Write the channel network junctions as point data. | Contains feature attributes: junction_number, receiver_junction, and stream_order. |

### Single-channel side-basin extraction

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `extract_basins_from_single_channel` | `bool` | `false` | `switch` | - | Build one main channel from an upstream point and extract the side basins draining into it. | Reads single_channel_points_fname (one or two points; any CRS). With one point the D8 path is traced to the DEM edge; with two, the upstream point is identified by D8 reachability and the channel is truncated at the common-downstream point where the other point's flow joins it. Side basins are then filtered by side_tributary_minimum_contributing_pixels and reject_incomplete_side_tributary_basins, and written with the side-tributary writers (main channel, outlets, basin outlines). |
| `single_channel_points_fname` | `string` | `NULL` | `param` | `extract_basins_from_single_channel` | Point vector file (1 or 2 points) for extract_basins_from_single_channel. | First/only point is the upstream seed; an optional second point is the downstream limit (identified and snapped automatically). Any OGR format; reprojected to the DEM CRS on read. |

### Side-tributary basin extraction from seed points

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `extract_side_tributary_basins_from_points` | `bool` | `false` | `switch` | - | Extract side-tributary basins draining into D8 flowlines traced from input points. | If true, reads a point vector file, traces D8 downstream from each point, finds side tributaries entering those traced channels, filters them by minimum contributing pixels, rejects incomplete basins if requested, and writes side-tributary outlet points and basin outlines. |
| `side_tributary_points_fname` | `string` | `NULL` | `param` | - | Input point vector file for side-tributary extraction. | Path to a GeoPackage or FlatGeobuf containing upstream seed points. If no path separator is present, read_path is prepended. Supported formats are those readable by GDAL/OGR. |
| `side_tributary_points_layer` | `string` | - | `param` | - | Layer name in the side-tributary input point vector. | If empty, the first layer in the vector dataset is used. |
| `side_tributary_minimum_contributing_pixels` | `int` | `1000` | `param` | - | Minimum contributing-pixel count for a side tributary basin. | Candidate side tributary outlets are only retained if their D8 upslope pixel count, including the outlet pixel itself, is at least this value. |
| `reject_incomplete_side_tributary_basins` | `bool` | `true` | `param` | - | Reject side basins touching nodata or the raster edge. | A side basin is incomplete if any pixel in its upslope set is adjacent to nodata or an open/out-of-domain edge according to the flow object's boundary conditions. |
| `include_side_tributary_stopping_point` | `bool` | `false` | `param` | - | Treat the last point in the input vector as a stopping point. | If true, the final feature in side_tributary_points_fname is not a source but a stopping point. Each D8 downstream trace along which side tributaries are found is truncated when it first enters a square window (half-width side_tributary_stopping_point_window_pixels) around this point. This limits side-tributary extraction to the reach upstream of the stopping point. The stopping point need not lie exactly on the main channel. |
| `side_tributary_stopping_point_window_pixels` | `int` | `5` | `param` | - | Half-width in pixels of the square stopping-point window. | Used only when include_side_tributary_stopping_point is true. A trace is truncated at the first pixel within this many pixels (Chebyshev distance) of the stopping point. Default 5 gives an 11x11 window. |
| `write_side_tributary_main_channel` | `bool` | `true` | `switch` | - | Write the traced main channel as vector lines. | The main channel is the set of D8 flow paths from the source points down to the outlet or stopping point, from which side tributaries are extracted. Output file: &lt;write_fname&gt;_side_tributary_main_channel.fgb or .gpkg. Uses side_tributary_write_vector_format. |
| `write_side_tributary_outlets` | `bool` | `true` | `switch` | - | Write accepted side-tributary outlet points. | The outlets are one pixel upstream from the D8 flowlines traced from the input points. Output file: &lt;write_fname&gt;_side_tributary_outlets.fgb or .gpkg. |
| `write_side_tributary_basin_outlines` | `bool` | `true` | `switch` | - | Write accepted side-tributary basin outlines. | Output file: &lt;write_fname&gt;_side_tributary_basin_outlines.fgb or .gpkg. |
| `side_tributary_write_vector_format` | `string` | `same_as_write_vector_format` | `param` | - | Vector output format for side-tributary outputs. | Valid options: flatgeobuf, fgb, gpkg. FlatGeobuf is recommended because it is a single-file binary GIS format and is more robust than GeoPackage on some filesystems. |

### Basin outlet input

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_outlet_fname` | `string` | `NULL` | `param` | - | File of requested basin outlet points (CSV, FlatGeobuf, or GeoPackage). | Required when select_basins_from_outlets or select_basins_from_outlets_nested is true. Format is auto-detected from the file extension: *.csv is read as a text file with latitude and longitude columns (WGS84 / EPSG:4326 decimal degrees); *.fgb or *.gpkg is read as a point vector layer in its own CRS, which is then transformed to the raster CRS. For CSV inputs the file must contain columns named exactly 'latitude' and 'longitude'. For vector inputs all features must be point geometry. |

### Basin outlet snapping

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_outlet_snapping_method` | `string` | `downstream` | `param` | `write_basin_outlet_points` | Optional snapping of basin outlet CSV points to pixels that exceed a contributing-pixel threshold. | Allowed values: none, expanding_square, downstream, largest_area. none: use the mapped raster pixel directly. expanding_square: search outward ring-by-ring until a qualifying pixel is found. downstream: follow the D8 receiver chain until a qualifying pixel is found. largest_area: take the pixel with the most contributing pixels in the square of half-width basin_outlet_snap_max_search_radius_pixels, if it is a qualifying pixel (a large radius can jump onto a bigger neighbouring river). |
| `basin_outlet_contributing_pixels_threshold` | `int` | `1000` | `param` | `write_basin_outlet_points` | Contributing-pixel threshold used for optional basin outlet snapping. | This is the total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. |
| `basin_outlet_snap_min_search_radius_pixels` | `int` | `0` | `param` | `write_basin_outlet_points` | Minimum expanding-square snapping search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square. 0 means the mapped input pixel is tested first. |
| `basin_outlet_snap_max_search_radius_pixels` | `int` | `25` | `param` | `write_basin_outlet_points` | Maximum expanding-square snapping search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square or largest_area. For expanding_square the search stops at the first radius containing one or more qualifying pixels; for largest_area it is the half-width of the square searched. |
| `write_basin_outlet_points` | `bool` | `false` | `switch` | - | Write requested basin outlet points to FlatGeobuf in the raster CRS. | Useful for debugging outlet locations after latitude/longitude input has been transformed into the DEM local/projected coordinate system. Points that fail coordinate transformation are omitted because they do not have valid raster-CRS coordinates. |

### Basin selection

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `select_basins_from_outlets` | `bool` | `false` | `switch` | - | Select basins from supplied basin outlet points. | Reads basin_outlet_fname, optionally snaps the mapped pixels, and then constructs one lsdtt_junction_basin per supplied basin outlet. |
| `select_basins_from_outlets_nested` | `bool` | `false` | `switch` | - | Select supplied basin outlets and all nested upstream basins. | Reads basin_outlet_fname, optionally snaps the mapped pixels, constructs one lsdtt_junction_basin per supplied basin outlet, and then expands selection to include all upstream nested basins from the junction-network s_vector. |
| `select_basins_from_line_crossings` | `bool` | `false` | `switch` | - | Select basins whose channels cross a supplied line (e.g. a fault). | Reads crossing_line_fname (a polyline vector that MUST be in the same CRS as the DEM). Channel pixels lying on the line whose contributing pixel count is at least crossing_minimum_contributing_pixels become basin outlets. Incomplete and nested basins are removed automatically, leaving the largest complete basin for each independent crossing. Outputs basin outlines and outlet points. |
| `crossing_line_fname` | `string` | `NULL` | `param` | - | Polyline vector file the channels must cross (same CRS as the DEM). | Required when select_basins_from_line_crossings is true. Any OGR-readable polyline layer (gpkg, fgb, shp, geojson). |
| `crossing_minimum_contributing_pixels` | `int` | `500` | `param` | - | Minimum contributing pixels at a line crossing to seed a basin outlet. | Channel-line crossings with fewer upslope pixels are ignored. |
| `select_all_basins` | `bool` | `false` | `switch` | - | Select one basin for every junction in the junction network. | Uses the junction network to derive a basin outlet flat index for each junction and constructs one lsdtt_junction_basin per junction. |
| `select_all_zero_order_basins` | `bool` | `false` | `switch` | - | Select one basin for every channel-head source pixel. | Zero-order basins are the complete upslope contributing areas of source pixels identified by threshold_contributing_pixels. This is separate from Strahler stream-order basin selection. |
| `select_all_basins_of_stream_order` | `bool` | `false` | `switch` | - | Select all basins whose junction stream order equals selected_stream_order. | This is equivalent to selecting all basins then filtering to one exact stream order. |
| `selected_stream_order` | `int` | `-1` | `param` | - | Exact stream order used by select_all_basins_of_stream_order. | Must be &gt; 0 when select_all_basins_of_stream_order is true. |

### Basin filtering

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `filter_basins_by_stream_order` | `bool` | `false` | `param` | - | Keep only basins within a stream-order range. | Only basins between basins_min_stream_order and basins_max_stream_order are kept. Use -1 for either bound to ignore that bound. |
| `basins_min_stream_order` | `int` | `-1` | `param` | - | Minimum stream order to keep when filtering basins. | Inclusive. Defaults to -1, which means ignore minimum. |
| `basins_max_stream_order` | `int` | `-1` | `param` | - | Maximum stream order to keep when filtering basins. | Inclusive. Defaults to -1, which means ignore maximum. |
| `filter_basins_by_area` | `bool` | `false` | `param` | - | Keep only basins within an area range in pixels. | Only basins between basins_min_area_pixels and basins_max_area_pixels are kept. Use -1 for either bound to ignore that bound. |
| `basins_min_area_pixels` | `int` | `-1` | `param` | - | Minimum basin area in pixels to keep when filtering basins by area. | Inclusive. Defaults to -1, which means ignore minimum. |
| `basins_max_area_pixels` | `int` | `-1` | `param` | - | Maximum basin area in pixels to keep when filtering basins by area. | Inclusive. Defaults to -1, which means ignore maximum. |
| `filter_incomplete_basins` | `bool` | `false` | `param` | - | Remove all incomplete basins. | An incomplete basin has at least one DEM pixel that is adjacent to an invalid pixel. Off-grid neighbours are invalid, unless an n or p boundary is supplied. |
| `filter_all_but_largest_basin` | `bool` | `false` | `param` | - | Keep only the largest basin. | Applied later in the filtering order so it can be used after completeness or other filters. |
| `filter_nested_basins` | `bool` | `false` | `param` | - | Remove nested basins. | Only the largest non-nested basins in the current selection are retained. |
| `filter_nested_by_stream_order` | `bool` | `false` | `param` | - | Remove nested basins of the same stream order. | In a consecutive chain of nested basins with the same order, only the downstream-most is preserved. Nested basins of different orders are allowed. |

### Basin printing to vector

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_basin_outlines` | `bool` | `false` | `switch` | - | Write selected basin outlines as vector features. | Features are labelled by their junction number. Use basin_outlines_write_vector_format to specify an override. To get basin outlines as polygons rather than polylines, use write_basin_outlines_as_polygon. To write basin outlines per stream order, use write_basin_outlines_per_order. |
| `write_basin_outlines_as_polygon` | `bool` | `false` | `switch` | - | Write basin outlines as polygons rather than polylines. | Set to true when you want basins printed as areas. |
| `write_basin_outlines_per_order` | `bool` | `false` | `switch` | - | Write basin outlines per stream order. | For GeoPackage output, this writes individual layers per stream order. For FlatGeobuf output, this writes multiple files, one per stream order. |
| `basin_outlines_write_vector_format` | `string` | `same_as_write_vector_format` | `param` | - | Output format for basin outline visualisation. | Allowed values: flatgeobuf, fgb, gpkg. FlatGeobuf is the default. GeoPackage may have file I/O issues on some filesystems. CSV is not allowed for basin outline lines. |
| `write_basin_outlines_to_csv` | `bool` | `false` | `switch` | - | Write basin outlines as a point CSV with distance and elevation. | Outputs a _basin_outlines.csv file containing points for all basins. |

### Basin outlet writing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_selected_basin_outlets_csv` | `bool` | `false` | `switch` | - | Write the outlets of the final selected basins to a CSV file. | Outputs a _selected_basin_outlets.csv file with longitude, latitude, basin_junction_number, basin_outlet_flat_index, and the basin size as basin_n_pixels, basin_area_km2 and basin_fraction_of_domain. This file can be directly read back into lsdtt as basin_outlet_fname. Note: write_basin_outlet_points writes the requested outlets, but this option writes the outlets of the selected basins. |

### Basin printing to raster

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_basins_raster` | `bool` | `false` | `switch` | - | Write the selected basins as a labelled raster. | Basins are labelled by junction number since that uniquely identifies extracted basins. Nested basins are printed on top of parent basins. |
| `write_basins_datamap` | `bool` | `false` | `switch` | - | Write selected basin datamap as CSV. | For each basin, there is a row containing the junction number and outlet flat index. This can be helpful for debugging. |

### Interior nodata

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `interior_nodata_method` | `string` | `sink` | `param` | - | How nodata pixels inside the DEM are treated: sink or fill. | Interior nodata are nodata pixels not connected to the raster edge through other nodata (isolated voids, moulins, masked lakes). sink (the default, unchanged behaviour): they are outlets, so flow that reaches them ends there, as at the DEM edge; use this for real internal drainage (e.g. moulins on an ice sheet). fill: they are interpolated with nodata_fill_method_choice before any filling or routing, so flow passes across them; use this for data voids. Nodata connected to the edge (e.g. sea removed by min_elevation) is never filled and always acts as an open boundary. |
| `nodata_fill_method_choice` | `string` | `laplacian` | `param` | `interior_nodata_method=fill` | Method for filling interior nodata pixels. | laplacian: smooth Dirichlet interpolation (recommended). iterative_weighted_average: fast, well-suited to thin voids. Used when interior_nodata_method is fill. |
| `nodata_fill_max_iter` | `int` | `10000` | `param` | `interior_nodata_method=fill` | Maximum iterations for iterative nodata fill methods. | Used when interior_nodata_method is fill. |

### Orientations

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `bearings_relative_to_true_north` | `bool` | `false` | `param` | - | Write bearings, aspects and orientations clockwise from true north. | Off (default): they are clockwise from grid north (grid up), so in an oblique frame with the feature along +x, 90 degrees means along the feature. On: each value is converted with the grid convergence at its own location (each cell of a raster, each row of a table), which also removes the up to about 3 degree difference between UTM grid north and true north. Angles between two directions (junction angles, deflections, contrasts) are unchanged, and hillshade illumination always stays relative to grid up. Needs a raster with a CRS. |

<!-- END GENERATED CLI PARAMETERS -->
