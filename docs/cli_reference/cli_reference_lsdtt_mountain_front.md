# CLI reference: lsdtt-mountain-front

Identifies mountain masses in a DEM via adaptive local-relief thresholding. The pipeline classifies each pixel as either mountain or flat/piedmont, then optionally extracts the points where drainage crosses the mountain front boundary.

## Usage

```
lsdtt-mountain-front <path> <parameter_file>
lsdtt-mountain-front <parameter_file>         # path defaults to ./
lsdtt-mountain-front -h                       # print help and exit
lsdtt-mountain-front -v                       # print version and exit
lsdtt-mountain-front -c                       # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Pipeline overview

The pipeline runs all stages unconditionally. Output writing is controlled by the `write_*` parameters.

| Stage | Description |
|---|---|
| 1 | Compute local relief (moving square window) |
| 2 | Build flat mask: adaptive CDF divergence or fixed threshold |
| 3 | Remove small flat components |
| 4 | Morphological closing of flat mask (seal narrow gaps) |
| 5 | Morphological opening of flat mask (sever thin peninsulas) |
| 6 | Invert flat mask → initial mountain mask |
| 7 | Fill holes in mountain mask (enclosed non-mountain pixels become mountain) |
| 8 | Remove small mountain components |
| 9 | Valley infilling via morphological closing (optional) |
| 10 | Identify mountain front outlet points (optional) |

The `_relief_cdf.csv` file is always written when `relief_threshold_algorithm_choice` is `cdf_divergence` (the default). It contains the 500-bin empirical CDF and the best-fit parametric CDF, and is the primary tool for choosing a manual threshold if the adaptive result is not satisfactory.

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input raster filename. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `write_output_format` | string | `cog` | Output raster format. `cog`, `gtiff`, or `same_as_input`. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` uses all available hardware threads. |

---

## Stage 1: Local relief

| Parameter | Type | Default | Description |
|---|---|---|---|
| `relief_window_radius_metres` | float | `0` | **Radius** (m) of the square window local relief (max − min) is measured over, converted to the nearest whole number of cells. `0` = 2 cells, a 5×5 window. The same name and meaning as in `lsdtt-flat-surface-finder` and `lsdtt-fluvial-surfaces`. Replaces `relief_window_side_pixels` (a full side in cells), which is still accepted and converted. |

---

## Stage 2: Flat mask threshold

| Parameter | Type | Default | Description |
|---|---|---|---|
| `relief_threshold_algorithm_choice` | string | `cdf_divergence` | `cdf_divergence` computes a 500-bin empirical CDF, fits the best parametric distribution, and picks the x-value where the empirical CDF most exceeds the fitted CDF. `fixed` uses `max_relief_metres`. Replaces the boolean `adaptive_relief_threshold` (`true` = `cdf_divergence`). |
| `max_relief_metres` | float | `1.0` | Fixed local-relief bound in metres below which a pixel is classed flat. Only used when `relief_threshold_algorithm_choice` is `fixed`. Read a good value from the `_relief_cdf.csv` of a `cdf_divergence` run. Replaces `relief_threshold_metres`. |

---

## Stages 3–8: Flat and mountain mask filtering

| Parameter | Type | Default | Description |
|---|---|---|---|
| `min_flat_component_pixels` | int | `500` | Minimum 8-connected component size (pixels) to retain in the flat mask (Stage 3). Removes small isolated flat patches before the mask is inverted. |
| `flat_closing_radius_metres` | float | `500` | Radius (m) for morphological closing of the flat mask (Stage 4). Seals narrow background gaps between adjacent flat patches, such as riser pixels between floodplain segments. Set to `0` to skip. |
| `erosion_radius_metres` | float | `500` | Erosion radius (m) for morphological opening of the flat mask (Stage 5). Severs thin flat peninsulas narrower than twice this radius; broad flat bodies are restored to their exact original extent. Set to `0` to skip. |
| `min_component_pixels` | int | `1000` | Minimum 8-connected component size (pixels) to retain in the mountain mask after inverting and hole-filling (Stage 8). |

---

## Stage 9: Valley infilling

| Parameter | Type | Default | Description |
|---|---|---|---|
| `infill_valley_radius_metres` | float | `0` | Radius (m) of the structuring element for morphological closing of the mountain mask. Seals valley notches cut by incised rivers into the mountain front. Set to approximately half the valley width to be bridged. `0` = disabled. |

---

## Stage 10: Mountain front outlets

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_mountain_front_outlets` | bool | `false` | Compute and write mountain front outlet points. Mountain mask edge pixels whose upstream contributing area meets `basin_minimum_contributing_pixels` are identified as outlets. Upstream duplicates (outlets that have another outlet downstream on the same flow path) are pruned. |
| `basin_minimum_contributing_pixels` | int | `1000` | Minimum upstream contributing pixels for a mountain edge pixel to qualify as an outlet. The former `min_contributing_pixels` spelling is an unchanged-value alias. |
| `contributing_pixels_fname` | string | _(empty)_ | Filename of a pre-computed `n_contributing_pixels` raster. If empty, the CLI fills pits and builds D8 flow routing from the DEM automatically. Providing a pre-computed raster skips the fill and routing steps. |

---

## Output parameters

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_relief_raster` | bool | `false` | `<write_fname>_relief.tif` — local relief raster. Useful for inspecting the threshold. |
| `write_flat_mask_raster` | bool | `false` | `<write_fname>_flat_mask_closing.tif` and `<write_fname>_flat_mask_opening.tif` — flat mask after Stage 4 and Stage 5 respectively. |
| `write_mountain_mask_raster` | bool | `true` | `<write_fname>_mountain_mask.tif` — final binary mountain mask. 1 = mountain, 0 = flat / piedmont. |
| `write_mountain_labels_raster` | bool | `false` | `<write_fname>_mountain_labels.tif` — each 8-connected mountain component labelled with a unique integer (UInt32). |

When `write_mountain_front_outlets` is true, two files are always written:

- `<write_fname>_mountain_front_outlets.fgb` — outlet points as FlatGeobuf with attributes `row`, `col`, `contributing_pixels`
- `<write_fname>_mountain_front_outlets.csv` — same outlets in WGS84 latitude/longitude

---

## Example parameter file

```
read_path:  /data/dems/
read_fname: my_dem
write_path: /data/output/
write_fname: my_dem_mountain

# Relief window: radius in metres (0 = 2 cells, a 5 x 5 window)
relief_window_radius_metres: 0

# Bound chosen from the relief CDF (inspect _relief_cdf.csv to tune)
relief_threshold_algorithm_choice: cdf_divergence

# Flat mask filtering
min_flat_component_pixels:    500
flat_closing_radius_metres:   500
erosion_radius_metres:        500
min_component_pixels:         1000

# Seal valley notches 500 m wide
infill_valley_radius_metres:  250

# Outputs
write_mountain_mask_raster:       true
write_relief_raster:              true
write_mountain_front_outlets:     true
basin_minimum_contributing_pixels: 1000
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

### -

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | OpenMP thread count. | 0 = all hardware threads. |
| `write_output_format` | `string` | `cog` | `io` | - | Raster output format. | cog, gtiff, or same_as_input. |

### Admission

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `relief_window_radius_metres` | `float` | `0` | `param` | - | Radius (m) of the square window local relief (max - min) is measured over. 0 = 2 cells. | Converted to the nearest whole number of cells, at least one; the run log states the radius actually used. A radius r cells is a (2r+1) x (2r+1) window, so the default is a 5 x 5 window. Replaces relief_window_side_pixels (a full side) and fallback_relief_window_radius_pixels (a radius in cells), which spelled the same idea two ways. |
| `relief_threshold_algorithm_choice` | `string` | `cdf_divergence` | `param` | - | How the relief bound is chosen: fixed or cdf_divergence. | fixed uses max_relief_metres. cdf_divergence computes the relief CDF, fits the best of a family of parametric distributions (normal, exponential, Weibull, Pareto) and takes the relief at which the empirical CDF most exceeds the fit: flat pixels inflate the empirical distribution at low relief. The CDF and the chosen value are written out. Replaces the boolean adaptive_relief_threshold. |
| `max_relief_metres` | `float` | `1` | `param` | `relief_threshold_algorithm_choice` | Local relief (m) below which a pixel is admitted as flat. | Read a value off the relief CDF the tool writes. Where the tool has relief_threshold_algorithm_choice, this is the value used by `fixed`. |

### Stages 4-7: flat mask post-processing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `min_flat_component_pixels` | `int` | `500` | `param` | - | Minimum 8-connected component size (pixels) to keep in the flat mask. | Small isolated flat patches (e.g. individual floodplain pixels) are removed before the flat mask is used to define mountain pixels. |
| `flat_closing_radius_metres` | `float` | `500` | `param` | - | Radius (m) for morphological closing of the flat mask (Stage 4). | Seals narrow background gaps between adjacent flat patches (e.g. riser pixels between floodplain segments). 0 = skip. |
| `erosion_radius_metres` | `float` | `500` | `param` | - | Erosion radius (m) for morphological opening of the flat mask (Stage 5). | Severs thin flat peninsulas narrower than 2x this radius; broad flat bodies are restored to their exact original extent. 0 = skip. |
| `min_component_pixels` | `int` | `1000` | `param` | - | Minimum 8-connected component size (pixels) to keep in the mountain mask. | Applied to the eroded mountain mask after neck severing. |

### Stage 9: valley infilling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `infill_valley_radius_metres` | `float` | `0` | `param` | - | Radius (m) of the circular SE for morphological closing (valley infilling). | Seals valley notches cut by incised rivers into the mountain front. Set to ~half the valley width to be filled. 0 = disabled. |

### Stage 10: outlets

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `contributing_pixels_fname` | `string` | - | `param` | - | Filename of a pre-computed n_contributing_pixels raster. | Leave empty to auto-compute via fill + D8 flow routing. |
| `basin_minimum_contributing_pixels` | `int` | `1000` | `param` | - | Minimum upstream contributing pixels for a mountain edge pixel to be an outlet. | Only used when write_mountain_front_outlets is true. |

### Outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_relief_raster` | `bool` | `false` | `switch` | - | Write the local relief raster. | Useful for threshold inspection. |
| `write_flat_mask_raster` | `bool` | `false` | `switch` | - | Write the processed flat mask (after component filter, dilation, hole fill). | The mountain mask is the inverse of this. |
| `write_mountain_mask_raster` | `bool` | `true` | `switch` | - | Write the final mountain mask raster. | 1 = mountain, 0 = flat / piedmont. |
| `write_mountain_labels_raster` | `bool` | `false` | `switch` | - | Write a labelled mountain components raster (UInt32). | Each 8-connected mountain component receives a unique integer label. |
| `write_mountain_front_outlets` | `bool` | `false` | `switch` | - | Write mountain front outlet points. | Mountain mask edge pixels whose upstream area meets basin_minimum_contributing_pixels. The container is set by write_vector_format. |
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Output VECTOR format: flatgeobuf, gpkg or csv. | Governs the outlet point layer. Distinct from write_output_format, which is the RASTER format. csv writes latitude/longitude in WGS84 instead of geometry in the raster CRS. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). GeoPackage is an SQLite database, and SQLite cannot take the file locks it needs on a drive mounted into WSL2, so a gpkg write there fails. |

### Interior nodata

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `interior_nodata_method` | `string` | `sink` | `param` | - | How nodata pixels inside the DEM are treated: sink or fill. | Interior nodata are nodata pixels not connected to the raster edge through other nodata (isolated voids, moulins, masked lakes). sink (the default, unchanged behaviour): they are outlets, so flow that reaches them ends there, as at the DEM edge; use this for real internal drainage (e.g. moulins on an ice sheet). fill: they are interpolated with nodata_fill_method_choice before any filling or routing, so flow passes across them; use this for data voids. Nodata connected to the edge (e.g. sea removed by min_elevation) is never filled and always acts as an open boundary. |
| `nodata_fill_method_choice` | `string` | `laplacian` | `param` | `interior_nodata_method=fill` | Method for filling interior nodata pixels. | laplacian: smooth Dirichlet interpolation (recommended). iterative_weighted_average: fast, well-suited to thin voids. Used when interior_nodata_method is fill. |
| `nodata_fill_max_iter` | `int` | `10000` | `param` | `interior_nodata_method=fill` | Maximum iterations for iterative nodata fill methods. | Used when interior_nodata_method is fill. |

<!-- END GENERATED CLI PARAMETERS -->
