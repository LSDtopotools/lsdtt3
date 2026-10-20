# CLI reference: lsdtt-lithology-erodibility

Turns a lithology map and a lithology → erodibility table into a stream-power
erodibility (K) raster for `lsdtt-lem` (`erodibility_fname`). The map can be a
polygon vector file with a lithology attribute (for example a Macrostrat
download) or a raster of integer lithology codes. The K raster is made on the
grid the LEM run will use: a DEM's grid, a synthetic grid (with the same
placements `lsdtt-lem` offers), or the map's own grid.

## Usage

```
lsdtt-lithology-erodibility <path> <parameter_file>
lsdtt-lithology-erodibility <parameter_file>         # path defaults to ./
lsdtt-lithology-erodibility -h                       # print help and exit
```

## The property table

```
# lsdtt lithology property table v1
# erodibility_convention: map_area
# m_exponent: 0.45
# n_exponent: 1
# calibration_cell_size_metres: 30
# provenance: <where the values come from>
code,name,erodibility,reference
10,sandstone,2.0e-5,"..."
30,limestone,1.0e-5,"..."
```

* `erodibility_convention`, `m_exponent` and `n_exponent` are **required**: a K
  value means nothing without its area convention, m and n.
  `map_area` (K in m^(1−2m) yr⁻¹) is independent of the cell size and is the
  recommended form; `specific_area` K is tied to `calibration_cell_size_metres`.
* `code` is an integer in 0–65534, unique. `name` is unique (case-insensitive).
  `erodibility` must be finite and > 0. `reference` is optional.
* Hillslope columns (`critical_hillslope_slope`, `hillslope_diffusivity`) are not
  applied (a warning says so); other extra columns are ignored.

No default table ships with lsdtt3; supply your own values and their sources.

## From map values to codes

Each map value must become a table code. In order:

1. With `lithology_crosswalk_fname` (CSV `value,code`), every map value present
   must be in the crosswalk. Use it to map many source terms onto a few units.
2. Otherwise a value that is an integer is a code, and any other value must
   match a table `name` (case-insensitive).

Every value that cannot be mapped is listed with its pixel count, and the run
stops. Model-domain pixels with no lithology at all (a hole in the map, or the
map's edge inside the domain) are refused unless `lithology_unmapped_code`
names a table unit for them.

For polygons each pixel takes the polygon containing its centre. A pixel whose
centre falls in a sliver gap between polygons takes the value covering most of
it; a pixel no polygon touches stays unmapped. Where two polygons with
different values claim a pixel centre, the value that sorts first wins and the
count is reported.

## Choosing the grid

| `model_grid_source` | The K raster's grid | Run `lsdtt-lem` with |
|---|---|---|
| `dem` (with `dem_fname`) | the DEM's | `read_initial_dem: true` on the same DEM |
| `parameters` | `n_rows` × `n_cols` at `cell_size_metres` on the synthetic georeference (Watopia unless `synthetic_crs_wkt`) | the same `n_rows`, `n_cols`, `cell_size_metres` |
| `lithology` | a code raster's own grid, or a vector's extent at `cell_size_metres` (required); the map's footprint is the domain | `model_grid_source: erodibility_raster` |
| `auto` (default) | `dem` when `dem_fname` is given, else `lithology` | |

`input_raster_placement` places the map on the grid exactly as in `lsdtt-lem`
(see its reference page for the worked examples): `as_georeferenced` (the
default; a vector is rasterised in the grid's CRS), `warp_to_model_grid` (code
rasters, `near` or `mode`), `relabel_whole_raster` (a code raster with the grid's
dimensions), and `relabel_centre_clip` / `relabel_lower_left_clip`, which take an
`n_rows` × `n_cols` window at `cell_size_metres` from the map and put it on the
synthetic grid. For a vector map the clip window is centred on (or starts at) the
clip point exactly. Clips never rescale the map; a map in a geographic CRS, Web
Mercator, or any CRS distorting distances by more than 1 % is refused.

## Outputs

| File | Description |
|---|---|
| `<write_fname>_erodibility.tif` | K per pixel, in the table's convention. Outside the domain: nodata. |
| `<write_fname>_lithology_code.tif` | The lithology code per pixel (`write_lithology_code_raster`). |
| `<write_fname>_lithology_legend.csv` | `code,name,erodibility,n_pixels,map_values`: each unit present, its pixel count and the map values that became it. |
| `<write_fname>_lem_erodibility_parameters.txt` | The `lsdtt-lem` lines for this raster: `erodibility_fname` and the table's stated m, n and area convention (and calibration cell size), plus how to set the grid. |

The same lines are printed at the end of the run. `lsdtt-lem` then converts the
table's K into its own convention and prints the factor (for a map-area table
on a specific-area run, dx^m).

## Example: a Macrostrat map on its own grid, then the LEM

```
# lithology.param
lithology_fname: macrostrat_units.gpkg
lithology_field_name: lith_primary
lithology_property_table_fname: my_lithology_k.csv
lithology_crosswalk_fname: macrostrat_to_units.csv
cell_size_metres: 30
write_fname: area
```

```
# lem.param (the erodibility lines come from area_lem_erodibility_parameters.txt)
model_grid_source: erodibility_raster
erodibility_fname: area_erodibility.tif
erodibility_raster_m_exponent: 0.45
erodibility_raster_n_exponent: 1
erodibility_raster_area_convention: map_area
erodibility: 2e-5
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

### Lithology input

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `lithology_fname` | `string` | - | `param` | - | Lithology map (in read_path): a polygon vector file or a categorical code raster. | Vectors: any OGR polygon layer (GeoPackage, FlatGeobuf, GeoJSON, Shapefile), e.g. a Macrostrat download. Rasters: integer lithology codes. Must be in a projected CRS in metres for clips and for taking the grid from it. |
| `lithology_source` | `string` | `auto` | `param` | - | auto, vector_map or code_raster. | auto decides from the file extension (.tif, .tiff, .asc, .bil, .img, .vrt, .flt, .nc are rasters). |
| `lithology_layer_name` | `string` | - | `param` | `lithology_fname` | Vector layer to read. Empty (default) = the first layer. | - |
| `lithology_field_name` | `string` | `lith_primary` | `param` | `lithology_fname` | Attribute of the vector layer that holds the lithology. | Values are matched to the property table by integer code, by table name (case-insensitive), or through lithology_crosswalk_fname. lith_primary is the first lithology of each unit in lsdtt-fetch-lithology output. |
| `lithology_crosswalk_fname` | `string` | - | `param` | `lithology_fname` | Optional CSV (in read_path) with columns value,code mapping map values to table codes. | When given, every map value present must be in it. Use it to map many source names (e.g. Macrostrat lithology terms) onto a few table units. |
| `lithology_property_table_fname` | `string` | - | `param` | - | Lithology property table (in read_path): code, name, erodibility per unit. | CSV with '#' header lines stating erodibility_convention (map_area or specific_area), m_exponent and n_exponent (all required) and optionally calibration_cell_size_metres and provenance, then a header row code,name,erodibility[,reference]. K must be &gt; 0. See the CLI reference for the format. |
| `lithology_unmapped_code` | `int` | `-1` | `param` | - | Table code given to model-domain pixels with no lithology. -1 (default) = refuse. | Pixels a map does not cover (holes, the edge of a map) inside the domain are refused unless a table unit is named here. With model_grid_source: lithology the map's footprint is the domain, so uncovered pixels are nodata instead. |

### Grid

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `model_grid_source` | `string` | `auto` | `param` | - | What defines the output grid: auto, dem, parameters or lithology. | dem: the grid of dem_fname (use the same DEM in lsdtt-lem). parameters: n_rows, n_cols, cell_size_metres on the synthetic georeference (the same grid lsdtt-lem builds from the same keys). lithology: the code raster's own grid, or a vector's extent at cell_size_metres (required), with the map's footprint as the domain; run lsdtt-lem with model_grid_source: erodibility_raster on the result. auto = dem when dem_fname is given, else lithology. |
| `dem_fname` | `string` | - | `param` | - | DEM (in read_path) whose grid the K raster is made on (model_grid_source: dem). | - |
| `n_rows` | `int` | `200` | `param` | `model_grid_source` | Rows of the synthetic grid (model_grid_source: parameters). | - |
| `n_cols` | `int` | `400` | `param` | `model_grid_source` | Columns of the synthetic grid. | - |
| `cell_size_metres` | `float` | `30` | `param` | `model_grid_source` | Cell size of the synthetic grid, or of the grid taken from a vector map. | - |
| `synthetic_crs_wkt` | `string` | - | `param` | `model_grid_source` | CRS of the synthetic grid. Empty = Watopia (EPSG:32758), as in lsdtt-lem. | - |
| `input_raster_placement` | `string` | `as_georeferenced` | `param` | - | How the lithology map is placed on the grid (the lsdtt-lem placements). | as_georeferenced: in place (vectors are rasterised in the grid's CRS; rasters must already be on the grid). warp_to_model_grid: rasters warped (near or mode). relabel_whole_raster: a code raster with the grid's rows, cols and cell size. relabel_centre_clip / relabel_lower_left_clip: an n_rows x n_cols window at cell_size_metres around / from the clip point in the map's CRS, relabelled onto the synthetic grid. relabel_* need model_grid_source: parameters. |
| `input_raster_clip_easting_metres` | `double` | `nan` | `param` | `input_raster_placement` | Easting of the clip centre or lower-left corner, in the map's CRS. | - |
| `input_raster_clip_northing_metres` | `double` | `nan` | `param` | `input_raster_placement` | Northing of the clip centre or lower-left corner, in the map's CRS. | - |
| `input_raster_resampling_method` | `string` | `near` | `param` | `input_raster_placement` | Resampling for a warped code raster: near or mode (clips always use near). | - |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_lithology_code_raster` | `bool` | `true` | `switch` | - | Write the lithology codes the K raster was made from. | Output: &lt;write_fname&gt;_lithology_code.tif. |
| `write_output_format` | `string` | `cog` | `io` | - | Raster output format: cog or gtiff. | - |

<!-- END GENERATED CLI PARAMETERS -->
