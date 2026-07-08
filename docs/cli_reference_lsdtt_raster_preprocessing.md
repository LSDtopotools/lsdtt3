# CLI reference: lsdtt-raster-preprocessing

Runs a configurable preprocessing pipeline on a DEM raster. Steps are applied in order; each step is optional and controlled by parameters in a plain-text parameter file.

## Usage

```
lsdtt-raster-preprocessing <path> <parameter_file>
lsdtt-raster-preprocessing <parameter_file>         # path defaults to ./
lsdtt-raster-preprocessing -h                       # print help and exit
lsdtt-raster-preprocessing -v                       # print version and exit
lsdtt-raster-preprocessing -c                       # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments and are ignored.
key: value
```

Keys use `snake_case`. Values have no quoting requirements. All keys and their defaults are listed below.

---

## Pipeline overview

Steps run in this order when enabled:

| Step | Trigger parameter |
|---|---|
| 0. Mosaic input tiles | `mosaic_file_list` |
| 1. Projection check / reprojection | always runs |
| 2. Threshold masking | `minimum_elevation`, `maximum_elevation` |
| 2b. Binary mask | `apply_binary_mask: true` |
| 2c. Clip to polygon | `clip_raster_to_polygon: true` |
| 3. Interior nodata fill | `write_filled_nodata_raster: true` |
| 4. Clip to max-valid rectangle | `write_clipped_rectangle_raster: true` |
| 5. Depression breaching | `write_breached_raster: true` or `breach_before_fill: true` |
| 6. Pit filling | `write_fill_raster: true` or `write_fill_diff_raster: true` |
| 7. Depression polygons | `write_depression_polygons: true` |
| 8. Hillshade | `write_hillshade_raster: true` |
| 8b. Multidirectional hillshade | `write_multidirectional_hillshade_raster: true` |
| 9. Vector clipping | `clip_vector_data_to_raster: true` |

---

## Input / output parameters

These are universal across all lsdtt3 CLIs.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input raster. |
| `read_fname` | string | parameter file stem | Input raster filename (with or without extension). |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `write_output_format` | string | `cog` | Output raster format. `cog` = Cloud Optimised GeoTIFF. `gtiff` = plain compressed GeoTIFF. `same_as_input` = match the input format. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` uses all available hardware threads. `1` disables parallelism. Applies to filling, breaching, nodata repair, and wavelet steps. |

---

## Step 0: Mosaic

Merge multiple input rasters into one before any other processing. When `mosaic_file_list` is set, `read_fname` is ignored — a warning is printed if both are present. `write_fname` must be set explicitly when using a mosaic.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `mosaic_file_list` | string | _(empty)_ | Path to a text file listing rasters to merge, one absolute path per line. Blank lines and lines beginning with `#` are ignored. Rasters with mismatched CRS, cell size, or misaligned origin are skipped with a warning. In overlap regions the last file in the list takes priority. Output: `<write_fname>_mosaic.tif`. |

---

## Step 1: Projection

The raster must be in a projected (metre-unit) CRS before analysis. This step checks the input and reprojects it if needed.

Output files produced by reprojection receive a suffix indicating the chosen CRS: `_utm`, `_polar` (NSIDC polar stereographic, used when |lat| > 75°), or `_lambert` (Lambert Azimuthal Equal Area, used for very wide rasters when `allow_lambert_projection` is true).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enforce_projection` | bool | `true` | If true and the raster is geographic or Web Mercator, reproject it automatically. If false and the raster is unprojected, the program exits with an error. |
| `force_reproject` | bool | `false` | Reproject even when the input is already in a projected CRS. Useful for converting an existing projection to the preferred UTM/polar/LAEA analysis CRS. |
| `resampling_method` | string | `cubic` | GDAL resampling algorithm used during reprojection. Options: `near`, `bilinear`, `cubic`, `cubicspline`, `lanczos`, `average`, `mode`, `max`, `min`, `med`, `q1`, `q3`, `sum`. |
| `grid_spacing` | float | `-1` | Target grid spacing in metres after reprojection. Negative (default) = auto-estimate from input pixel size, rounded to the nearest 5 m (≥20 m inputs), 2 m (5–20 m), or 1 m (<5 m). Set a positive value to enforce a specific spacing. |
| `force_resampling` | bool | `false` | If true, use `forced_grid_spacing` as the target spacing instead of the auto-estimated or user-supplied `grid_spacing`. |
| `forced_grid_spacing` | float | `-1` | Target grid spacing in metres when `force_resampling` is true. Must be positive. |
| `allow_lambert_projection` | bool | `false` | Use Lambert Azimuthal Equal Area for rasters spanning three or more UTM zones. |
| `warp_mem_mb` | int | `512` | GDALWarp working memory in MB. Larger values reduce source tile re-reads for large rasters. Minimum ~64 MB. |

---

## Step 2: Threshold masking

Pixels outside the elevation range are set to nodata. The masked raster is written as `<write_fname>_remove_seas.tif` when any pixels are changed.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `minimum_elevation` | float | `0.0` | Pixels at or below this elevation (metres) are set to nodata. The default of 0.0 removes ocean and below-sea-level artefacts common in global DEMs. Set to a large negative number (e.g. `-9999`) to disable. |
| `maximum_elevation` | float | `30000` | Pixels above this elevation (metres) are set to nodata. Reduce to remove artefacts in low-relief study areas. |

---

## Step 2b: Binary mask

Apply an external binary raster to keep only pixels that overlap the foreground of the mask. Runs after threshold masking.

Output: `<write_fname>_binary_masked.tif`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `apply_binary_mask` | bool | `false` | Enable binary mask application. |
| `binary_mask_fname` | string | _(empty)_ | Path to the binary mask raster. Any GDAL-readable single-band raster is accepted. Pixels > 0.5 are treated as foreground (keep); all others set the DEM pixel to nodata. Must cover the same extent and have the same cell size and origin as the (possibly reprojected) DEM. Required when `apply_binary_mask` is true. |

---

## Step 2c: Clip raster to polygon

Crop the DEM to the bounding box of a polygon, and optionally mask pixels that fall outside the polygon shape.

Output: `<write_fname>_polygon_clipped.tif`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `clip_raster_to_polygon` | bool | `false` | Enable polygon clipping. |
| `polygon_clip_fname` | string | _(empty)_ | Path to the polygon vector file. Any OGR-readable polygon format is accepted (GeoPackage, FlatGeobuf, Shapefile, GeoJSON, etc.). The first layer and first polygon are used. A warning is printed if the file contains more than one polygon. Required when `clip_raster_to_polygon` is true. |
| `force_rectangular_raster` | bool | `false` | If false (default), pixels inside the bounding box but outside the polygon shape are set to nodata. If true, all pixels within the bounding box are kept — a plain rectangular crop. |
| `polygon_clip_rule` | string | `majority_coverage` | Rule for deciding whether a pixel on the polygon boundary is kept. `any_coverage`: pixel is kept if the polygon touches any part of it. `majority_coverage`: pixel is kept if its centre falls inside the polygon. `full_coverage`: pixel is kept only if entirely inside. Ignored when `force_rectangular_raster` is true. |

---

## Step 3: Interior nodata fill

Fill nodata pixels that are completely surrounded by valid data (voids that do not touch the raster edge). Recommended before clipping to a rectangle.

Output: `<write_fname>_nodata_filled.tif`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_filled_nodata_raster` | bool | `false` | Enable interior nodata filling and write the result. |
| `nodata_fill_method_choice` | string | `laplacian` | Filling algorithm. `laplacian`: smooth Dirichlet interpolation (recommended). `iterative_weighted_average`: faster; well-suited to thin voids. |
| `nodata_fill_max_iter` | int | `10000` | Maximum iterations for iterative methods. |

---

## Step 4: Clip to maximum valid rectangle

Clip the raster to the largest axis-aligned rectangle that contains no nodata pixels. Running step 3 first maximises the size of this rectangle.

Output: `<write_fname>_clipped_rect.tif`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_clipped_rectangle_raster` | bool | `false` | Enable rectangle clipping and write the result. |

---

## Step 5: Depression breaching

Carve flow paths through topographic barriers to eliminate closed depressions. Breaching preserves more of the original surface than filling; it is particularly effective on lidar DEMs.

Output files:

| File | Trigger |
|---|---|
| `<write_fname>_breached.tif` | `write_breached_raster: true` |
| `<write_fname>_carve_diff.tif` | `write_carve_diff_raster: true` |

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_breached_raster` | bool | `false` | Write the breached DEM. |
| `write_carve_diff_raster` | bool | `false` | Write the difference raster (breached minus original). Negative values show where the surface was lowered. |
| `breach_before_fill` | bool | `false` | Run breaching before pit filling. The breached DEM is then used as input to the fill step. Recommended for lidar DEMs where channelised flow paths must be carved before residual pits are filled. |
| `fill_remainder_after_breach` | bool | `true` | Fill pits that could not be fully breached. Used with `breach_before_fill` or least-cost breaching. |
| `breach_algorithm_choice` | string | `lindsay2016` | Breaching algorithm. `lindsay2016`: Lindsay (2016) efficient breaching. `least_cost`: Lindsay (2020) pit-centric least-cost breaching. |
| `breach_lc_max_dist` | int | `1000` | Maximum search radius in pixels for least-cost breaching. Only used when `breach_algorithm_choice` is `least_cost`. |
| `breach_lc_max_cost` | float | unlimited | Maximum total elevation lowering allowed for a breach path. Paths requiring more lowering are skipped. Only used when `breach_algorithm_choice` is `least_cost`. |
| `breach_lc_flat_increment` | float | `0.0` | Small elevation increment added along flat breach paths to ensure drainage. Only used when `breach_algorithm_choice` is `least_cost`. |

---

## Step 6: Pit filling

Fill closed depressions so that every pixel drains to the raster edge. Run after breaching if `breach_before_fill` is true.

Output files:

| File | Trigger |
|---|---|
| `<write_fname>_fill.tif` | `write_fill_raster: true` |
| `<write_fname>_fill_diff.tif` | `write_fill_diff_raster: true` |

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_fill_raster` | bool | `false` | Write the pit-filled DEM. |
| `write_fill_diff_raster` | bool | `false` | Write the difference raster (filled minus original). Positive values show where the surface was raised. For `barnesetal2021` this represents water depth standing in depressions. |
| `fill_algorithm_choice` | string | `barnes2014` | Filling algorithm. `barnes2014`: Barnes et al. (2014) Improved Priority-Flood. `barnesetal2021`: Barnes et al. (2021) Fill-Spill-Merge — note this algorithm does not guarantee that all depressions are filled (see doi:10.5194/esurf-9-105-2021). |
| `fill_min_slope` | float | `0.0001` | Minimum gradient imposed between adjacent filled cells to prevent flat regions that produce ambiguous flow directions. Increase slightly for very coarse DEMs. Not used with `barnesetal2021`. |
| `fill_spill_merge_rainfall_rate` | float | `0.1` | Uniform rainfall depth in elevation units per cell for the Fill-Spill-Merge algorithm. Only used when `fill_algorithm_choice` is `barnesetal2021`. |
| `fsm_boundary_conditions` | string | `oooo` | Boundary conditions for the depression hierarchy used by Fill-Spill-Merge, given as four characters in north–east–south–west order. `o` = open (edge drains freely to the ocean). `n` = no-flux (edge is not an ocean boundary). `p` = periodic (edge wraps; not ocean). Example: `oooo` for a typical DEM tile; `nnnn` for a closed basin. Only used when `fill_algorithm_choice` is `barnesetal2021`. |

---

## Step 7: Depression polygons

Trace closed depressions identified by the fill difference raster into vector polygons. Requires that `write_fill_raster` or `write_fill_diff_raster` is also true (the fill diff must be computed).

Output: `<write_fname>_depressions.<ext>` (extension depends on `depression_polygons_output_format`)

Each polygon carries three attributes: `depression_id`, `area_m2`, `volume_m3`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_depression_polygons` | bool | `false` | Enable depression polygon output. Automatically enables `write_fill_diff_raster`. |
| `depression_min_pixels` | int | `1` | Minimum depression size in pixels. Depressions smaller than this are removed before polygon tracing. |
| `depression_min_depth` | float | `0.0` | Minimum fill depth in metres for a cell to be included in a depression. A small positive value (e.g. `0.01`) suppresses near-zero fill artefacts. |
| `depression_polygons_output_format` | string | `flatgeobuf` | Output vector format. `flatgeobuf` or `fgb`: FlatGeobuf (recommended). `gpkg`: GeoPackage. |

---

## Step 8: Hillshade

Compute a hillshade raster using the Zevenbergen and Thorne (1987) algorithm. When pit filling is also enabled, the hillshade is computed from the filled DEM.

Output: `<write_fname>_hillshade.tif`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_hillshade_raster` | bool | `false` | Enable hillshade computation and write the result. |
| `hillshade_azimuth` | float | `315` | Sun azimuth in degrees clockwise from north. 315° (northwest) is the standard cartographic convention. |
| `hillshade_altitude` | float | `45` | Sun altitude above the horizon in degrees. |
| `hillshade_z_factor` | float | `1.0` | Vertical exaggeration applied before computing the hillshade. Values greater than 1 emphasise relief in flat terrain. |

---

## Step 8b: Multidirectional hillshade

Compute a multidirectional hillshade by averaging hillshades from evenly spaced azimuth directions. This reduces directional bias and improves visual continuity of features that are only strongly illuminated from certain directions.

When pit filling is also enabled, this step uses the filled DEM (same source choice as Step 8).

Output: `<write_fname>_hillshade_multidirectional.tif`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_multidirectional_hillshade_raster` | bool | `false` | Enable multidirectional hillshade computation and write the result. |
| `multidirectional_hillshade_n_directions` | int | `8` | Number of evenly spaced azimuth directions. Must be `>= 2`. Typical values: `8`–`16`. |
| `multidirectional_hillshade_start_azimuth` | float | `0.0` | Starting azimuth in degrees clockwise from north. Subsequent azimuths are sampled uniformly around 360°. |
| `hillshade_altitude` | float | `45` | Sun altitude used for each direction. |
| `hillshade_z_factor` | float | `1.0` | Vertical exaggeration used for each direction. |

---

## Step 9: Vector clipping

Clip a vector file to the bounding box of the output raster. The vector CRS need not match the raster CRS — the raster bounding box is reprojected into the vector's native CRS before clipping.

Output: `<write_fname>_clipped_vector.<ext>` (extension depends on `vector_clip_output_format`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `clip_vector_data_to_raster` | bool | `false` | Enable vector clipping. |
| `vector_data_fname` | string | _(empty)_ | Path to the input vector file. Any OGR-readable format is accepted (Shapefile, GeoPackage, FlatGeobuf, GeoJSON, etc.). Required when `clip_vector_data_to_raster` is true. |
| `vector_clip_output_format` | string | `flatgeobuf` | Output vector format. `flatgeobuf` or `fgb`: FlatGeobuf (preferred, especially on WSL2). `gpkg`: GeoPackage. |

---

## Example parameter file

```
read_path:  /data/dems/
read_fname: my_dem
write_path: /data/output/
write_fname: my_dem_preprocessed

# Projection
enforce_projection:       true
resampling_method:        cubic
grid_spacing:             -1
allow_lambert_projection: false

# Masking
minimum_elevation:        0.0
maximum_elevation:        30000

# Fill interior nodata voids first, then clip to a clean rectangle
write_filled_nodata_raster:    true
write_clipped_rectangle_raster: true

# Breach lidar-style DEMs, then fill residual pits
breach_before_fill:        true
breach_algorithm_choice:   lindsay2016
write_fill_raster:         true
write_fill_diff_raster:    true

# Hillshade from the filled DEM
write_hillshade_raster:    true
hillshade_azimuth:         315
hillshade_altitude:        45
hillshade_z_factor:        1.0
```
