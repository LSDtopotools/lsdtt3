# CLI reference: lsdtt-fetch-river-lines

**Status: experimental.** This page lists the parameters only. A full guide
will be written once the tool has been tested on real catchments.

Reads river centrelines for an area, from OS Open Rivers (Great Britain,
downloaded once and cached) or from a line file you already have, and
prepares them for `lsdtt-stream-burning`: single-part lines, oriented
downstream, joined at junctions, clipped to the area plus a buffer, in the
DEM's CRS. It reports the network's topology (tips, outlets, loops,
bifurcations, direction conflicts) and does not repair it.

## Usage

```
lsdtt-fetch-river-lines <path> <parameter_file>
lsdtt-fetch-river-lines <parameter_file>         # path defaults to ./
lsdtt-fetch-river-lines -h                       # print help and exit
```

## Outputs

| File | Contents |
|---|---|
| `<write_fname>_river_lines.fgb` | Prepared lines: `line_id`, `source`, `source_feature_id`, `name`, `feature_class`, `direction_source`, `on_loop`, `component_id`, `from_node`, `to_node`, `n_splits`, `culverted`, `strahler_order`, `length_m` |
| `<write_fname>_river_nodes.fgb` | Every node: `node_id`, `node_type`, `degree_in`, `degree_out`, `component_id` |
| `<write_fname>_river_topology.csv` | One row per connected component |
| `<write_fname>_ATTRIBUTION.txt` | Licence and citation of the source |

## Parameters

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

### Bounding box

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `west` | `string` | - | `param` | - | Western edge of the bounding box, WGS84 decimal degrees. | Required. Must be less than east. Range -180 to 180. |
| `east` | `string` | - | `param` | - | Eastern edge of the bounding box, WGS84 decimal degrees. | Required. Must be greater than west. Range -180 to 180. |
| `south` | `string` | - | `param` | - | Southern edge of the bounding box, WGS84 decimal degrees. | Required. Must be less than north. Range -90 to 90. |
| `north` | `string` | - | `param` | - | Northern edge of the bounding box, WGS84 decimal degrees. | Required. Must be greater than south. Range -90 to 90. |
| `bounds_fname` | `string` | - | `param` | - | Vector file (or raster) defining the area, instead of west/south/east/north. | FlatGeobuf (.fgb), GeoPackage (.gpkg) or any other OGR-readable format. The layer must contain either one or more polygons (the envelope of all of them is used), or exactly two points marking opposite corners. For two points either diagonal works -- lower-left/upper-right and lower-right/upper-left give the same box, so the orientation is worked out automatically. If the layer is in a projected CRS such as UTM it is reprojected to WGS84 for you, so a box drawn in a GIS can be used directly. A raster with a georeference (e.g. an existing DEM GeoTIFF) works too: its footprint is reprojected to WGS84. Set this or the four edge parameters, not both. |
| `bounds_path` | `string` | - | `param` | - | Optional directory for bounds_fname. | If empty, read_path is used. Ignored when bounds_fname is an absolute path. |
| `bounds_layer` | `string` | - | `param` | - | Layer name within bounds_fname. | If empty, the first layer is used (with a warning when the file has more than one). |

### River line source

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `river_line_source` | `string` | - | `switch` | - | Where the river lines come from: local_file or gb_os_open_rivers. | local_file: any OGR-readable line file you already have (river_lines_fname), for example OS MasterMap Water Network or an Open Rivers extract. gb_os_open_rivers: OS Open Rivers (Great Britain, OGL v3), downloaded once (about 52 MB) into river_cache_dir and read from there afterwards. Datasets derived from a DEM (HydroRIVERS, MERIT Hydro, CCM2, REC2) are deliberately not offered: burning one into a DEM tests the DEM against itself. |
| `river_lines_fname` | `string` | - | `io` | `river_line_source` | The line file for river_line_source: local_file. | Relative to read_path. Any format OGR reads (FlatGeobuf, GeoPackage, shapefile, ...), including /vsizip/ paths. Lines and multilines are read; other geometries are counted and skipped. |
| `river_lines_layer` | `string` | - | `io` | `river_line_source` | Layer of river_lines_fname to read. Empty = the first layer. | - |
| `river_cache_dir` | `string` | - | `io` | `river_line_source` | Directory where national downloads are kept and reused. Empty = ~/.cache/lsdtt3/river_lines. | Open Rivers is stored as &lt;dir&gt;/gb_os_open_rivers/oprvrs_gpkg_gb.zip and read in place. Any oprvrs_gpkg*.zip already in that folder is used without a download, so a copy can be placed there by hand. Delete it to fetch a newer release. |

### Local file attributes

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `feature_class_field_name` | `string` | - | `param` | `river_line_source` | local_file: field giving each line's kind. Empty = every line is a river. | Values are mapped to river, stream, lake_connector, canal, drain, tidal or other (Open Rivers 'form' and OSM 'waterway' values are recognised: inlandRiver, lake, tidalRiver, canal, stream, ditch, ...). See river_feature_classes. |
| `direction_field_name` | `string` | - | `param` | `river_line_source` | local_file: field saying whether a line is digitised downstream. | Needs direction_field_mapping. A value not in the mapping leaves that line's direction unknown (direction_source = unknown). |
| `direction_field_mapping` | `string` | - | `param` | `direction_field_name` | local_file: value=1 (digitised downstream) or value=-1 (digitised upstream), comma separated. | Write a space inside a value as an underscore (a parameter value cannot contain a space); attribute values are compared the same way. Example for Open Rivers, whose values are 'in direction' and 'in opposite direction': in_direction=1,in_opposite_direction=-1. |
| `digitised_downstream` | `bool` | `false` | `param` | `river_line_source` | local_file without a direction field: trust the digitising order as downstream. | Sets direction_source = digitised. Otherwise direction_source is unknown and the lines are kept as digitised. |
| `start_node_field_name` | `string` | - | `param` | `river_line_source` | local_file: field holding the id of the node at each line's digitised start. | With end_node_field_name, lines whose ends share an id are joined exactly, even when the coordinates differ by more than node_snap_tolerance_metres (Open Rivers: start_node). |
| `end_node_field_name` | `string` | - | `param` | `river_line_source` | local_file: field holding the id of the node at each line's digitised end. | See start_node_field_name (Open Rivers: end_node). |
| `culvert_field_name` | `string` | - | `param` | `river_line_source` | local_file: field marking culverted or tunnelled lines (OSM: tunnel). | culvert, yes, 1 or true gives culverted = 1, anything else 0. Empty = culverted is -1 (unknown) on every line. |

### Preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `output_epsg` | `int` | `0` | `param` | - | EPSG code of the outputs. 0 = the CRS of bounds_fname when it is a raster, else the WGS 84 UTM zone of the box centre. | Must be a projected CRS in metres. Writing the lines in the DEM's own CRS means the tolerances here mean the same as in lsdtt-stream-burning, and the lines sit on the DEM in a GIS without reprojection. Web Mercator is refused: it distorts distances. |
| `clip_buffer_metres` | `float` | `-1` | `param` | - | Buffer beyond the area before clipping [m]. Negative = 5% of the shorter side, at least 200 m. | Lines are clipped to the area plus this buffer, not to the exact edge, so lines that leave and re-enter near the edge stay whole and lsdtt-stream-burning sees complete edges. Ends cut by the clip are node_type clip_edge. |
| `river_feature_classes` | `string` | `river,stream,lake_connector` | `param` | - | Comma list of feature classes to keep: river, stream, lake_connector, canal, drain, tidal, other. | Lake connectors are kept by default because dropping them breaks the network at every lake. Canals and drains are dropped because they make loops and cross divides; tidal reaches because they put outlets in the wrong place. Counts of every class read are written to the log. |
| `node_snap_tolerance_metres` | `float` | `-1` | `param` | - | Line ends this close are joined, and a line is split where another line's end lies this close to it [m]. Negative = 1 m. | No vertex moves further than this. The vector is used by stream burning for topology only, so its geometry is not otherwise changed. |
| `min_line_length_metres` | `float` | `0` | `param` | - | Drop isolated single lines shorter than this [m]. 0 = keep everything. | A line connected to any other line is never dropped. |

### Network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `http_timeout_seconds` | `int` | `60` | `param` | - | Maximum seconds for an individual HTTP request. | Bounds how long an unresponsive server can stall the fetch. |
| `http_connect_timeout_seconds` | `int` | `15` | `param` | - | Maximum seconds to establish a connection. | Keep this short so an offline machine fails promptly. |
| `http_max_retry` | `int` | `3` | `param` | - | Number of retries on a transient HTTP failure. | Cloud storage occasionally returns 5xx responses; retrying absorbs those without failing the run. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). CSV is available only through selectors for point/table outputs. Per-output overrides exist only where a container genuinely differs. Raster outputs are governed separately by write_output_format. |
| `calculate_spatial_index` | `bool` | `true` | `io` | - | Whether to calculate a spatial index for vector data writes. | A spatial index helps a GIS speed up spatial queries, so it is often helpful to have one. But calculating it requires the entire vector layer to sit in memory, preventing feature-by-feature writes and adding a lot of overhead. For large datasets it may therefore be preferable to write without an index; one can always be built later with e.g. gdal. |

<!-- END GENERATED CLI PARAMETERS -->
