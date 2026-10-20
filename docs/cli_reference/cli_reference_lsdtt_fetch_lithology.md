# CLI reference: lsdtt-fetch-lithology

Downloads geological map units from [Macrostrat](https://macrostrat.org) for a
bounding box and writes them as a polygon layer in the UTM zone of the box,
ready for `lsdtt-lithology-erodibility`. Each unit carries `lith_primary`, its
first lithology, which is the field `lsdtt-lithology-erodibility` reads by
default.

Macrostrat data are licensed CC BY 4.0. Every run writes an attribution file
naming Macrostrat and each source map used; cite them in work that uses the
output.

## Usage

```
lsdtt-fetch-lithology <path> <parameter_file>
lsdtt-fetch-lithology <parameter_file>         # path defaults to ./
lsdtt-fetch-lithology -h                       # print help and exit
```

## What it does

1. Fetches Macrostrat's vector tiles (`https://tiles.macrostrat.org/carto/{z}/{x}/{y}.mvt`,
   layer `units`) covering the box at `macrostrat_tile_zoom`. The `carto` layer
   blends source maps of different scales and serves more detailed ones at higher
   zoom. Tile geometry is quantised to 1/4096 of a tile: about 10 m at zoom 10.
2. Keeps the polygons whose extent overlaps the box. Whole polygons are kept, so
   the output reaches beyond the box; tiles also overlap slightly at their edges,
   so a unit can appear as more than one polygon.
3. With `lookup_map_units: true` (the default), asks the Macrostrat API for each
   map unit's linked lithologies, the lithology vocabulary and the source maps.
   If the API cannot be reached the run continues with a warning:
   `lith_primary` then comes from the tiles' text and the attribution lists
   source ids only.
4. Reprojects to the WGS 84 UTM zone of the box centre (or `output_epsg`) and
   writes the layer.

The tiles need GDAL's MVT driver (`ogrinfo --formats | grep MVT`).

## lith_primary

Macrostrat's lithology text has no single format (`Major:{granite,diorite}, Minor:{...}`,
`plutonic: granite`, prose, or nothing). `lith_primary` is, in order:

| `lith_primary_source` | `lith_primary` |
|---|---|
| `major` | the first term of a `Major:{...}` list |
| `liths` | the name of the unit's first linked Macrostrat lithology |
| `text` | the whole lithology text |
| `none` | empty |

Values are lower case with whitespace collapsed. Free-text values are usually
too varied to match a property table by name; map them onto your units with a
crosswalk (below).

## Outputs

| File | Description |
|---|---|
| `<write_fname>_macrostrat_units.fgb` (or `.gpkg`) | Layer `units`: `map_id`, `source_id`, `name`, `lith` (Macrostrat's text), `lith_primary`, `lith_primary_source`, `lith_ids`, `lith_names`, `lith_type`, `lith_class`, `best_age_top`, `best_age_bottom` (Ma), `ref_name`. |
| `<write_fname>_lithology_vocabulary.csv` | One row per `lith_primary`: its source, type, class, number of polygons and area (km², overlaps counted twice), largest first. |
| `<write_fname>_crosswalk_template.csv` | `value,code,lith_class,area_km2` with `code` blank: fill it in to make a `lithology_crosswalk_fname` for `lsdtt-lithology-erodibility`. |
| `<write_fname>_ATTRIBUTION.txt` | The citation for Macrostrat and every source map. |

## Example: fetch, then make a K raster

```
# fetch.param
west: -3.25
south: 55.90
east: -3.15
north: 55.97
write_fname: edinburgh
```

Fill in `edinburgh_crosswalk_template.csv` with the codes of your property table
and save it as `edinburgh_crosswalk.csv`, then:

```
# lithology.param
lithology_fname: edinburgh_macrostrat_units.fgb
lithology_field_name: lith_primary
lithology_crosswalk_fname: edinburgh_crosswalk.csv
lithology_property_table_fname: my_lithology_k.csv
cell_size_metres: 30
write_fname: edinburgh
```

## Example: fetch, then k_sn by lithology

The same units file can feed `lsdtt-chi-analysis` directly, without a property
table. The crosswalk template's columns are `value,code,lith_class,area_km2`; for the
chi crosswalk, keep `value` and `code` and add a `name` column (the chi reader uses
the name, and the other columns are ignored):

```
# fetch.param (the Allegheny Front box used in the West Virginia examples)
west: -79.45
south: 38.85
east: -79.15
north: 39.05
write_fname: wvfront
```

```
value,code,name
sandstone,10,sandstone
siltstone,20,siltstone
shale,30,shale
limestone,40,limestone
fine alluvium,90,alluvium
```

```
# chi.param (after lsdtt-fetch-raster has written wvfront_DEM.tif)
fit_chi_networks: true
lithology_fname: wvfront_macrostrat_units.fgb
lithology_field_name: lith_primary
lithology_crosswalk_fname: wvks_crosswalk_lith_primary.csv
lithology_contact_buffer_metres: 0
write_lithology_code_raster: true
```

Each unit's k_sn, its reach statistics and a bootstrap interval are written per m/n
value. See the [chi analysis reference](cli_reference_lsdtt_chi_analysis.md#k_sn-by-lithologic-unit)
for the outputs and the caveats on reading the contrasts.

## Working offline

`macrostrat_tile_dir` reads tiles from a directory laid out as
`<dir>/{z}/{x}/{y}.mvt` (or `.pbf`) instead of downloading them. With
`lookup_map_units: false` nothing is downloaded at all.

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

### Macrostrat

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `macrostrat_tile_zoom` | `int` | `10` | `param` | - | Zoom of the Macrostrat vector tiles fetched (0-16). | The blended 'carto' layer serves more detailed source maps at higher zoom, and the tile geometry is quantised to about 1/4096 of a tile (about 10 m at zoom 10, 40 m at zoom 8). Higher zoom means more tiles for the same box. |
| `max_tiles` | `int` | `256` | `param` | - | Refuse a request needing more tiles than this. | A guard against an accidental continental request. Raise it deliberately for large areas, or lower macrostrat_tile_zoom. |
| `lookup_map_units` | `bool` | `true` | `param` | - | Look up each map unit's linked lithologies and the source maps in the Macrostrat API. | Adds lith_ids, lith_names, lith_type and lith_class and the full attribution. Without it lith_primary comes from the tiles' lith text only. |
| `request_delay_milliseconds` | `int` | `200` | `param` | - | Pause after every request to Macrostrat. | Macrostrat publishes no rate limit; this keeps a large fetch polite. |
| `output_epsg` | `int` | `0` | `param` | - | EPSG code of the output layer. 0 (default) = the WGS 84 UTM zone of the box centre. | Must be a projected CRS in metres for lsdtt-lithology-erodibility; Web Mercator (EPSG:3857, the tiles' own) is refused by lsdtt because it distorts distances. |
| `macrostrat_tile_dir` | `string` | - | `param` | - | Read the tiles from a local directory instead of downloading them. | Tiles must be at &lt;dir&gt;/{z}/{x}/{y}.mvt (or .pbf), e.g. a previously downloaded set. Combine with lookup_map_units: false to work entirely offline. |

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

<!-- END GENERATED CLI PARAMETERS -->
