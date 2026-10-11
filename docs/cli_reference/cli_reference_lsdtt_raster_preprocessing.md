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

!!! note "Multipart points"
    `point_data_fname` accepts multipart geometry (`MultiPoint` and friends). Each feature is split into its parts on read, with the attributes copied onto each, and the split is reported in the log.

    Shapefiles are the commonest source of multipart data, and **versions before August 2026 discarded it silently** — see [Vector data](../what_is_lsdtt3/vector_data.md).

## Pipeline overview

Steps run in this order when enabled:

| Step | Trigger parameter |
|---|---|
| 0. Mosaic input tiles | `mosaic_file_list` |
| 1. Projection check / reprojection | always runs |
| 2. Threshold masking | `min_elevation`, `maximum_elevation` |
| 2b. Binary mask | `run_binary_mask: true` |
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
| `resampling_method` | string | `cubic` | GDAL resampling algorithm used during reprojection. Cubic is the default for continuous elevations; use `near` or `mode` for categorical rasters. Other options include `bilinear`, `cubicspline`, `lanczos`, `average`, `max`, `min`, `med`, `q1`, `q3`, `sum`. |
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
| `min_elevation` | float | `0.0` | Pixels at or below this elevation (metres) are set to nodata. The default of 0.0 removes ocean and below-sea-level artefacts common in global DEMs. Set to a large negative number (e.g. `-9999`) to disable. |
| `maximum_elevation` | float | `30000` | Pixels above this elevation (metres) are set to nodata. Reduce to remove artefacts in low-relief study areas. |

---

## Step 2b: Binary mask

Apply an external binary raster to keep only pixels that overlap the foreground of the mask. Runs after threshold masking.

Output: `<write_fname>_binary_masked.tif`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `run_binary_mask` | bool | `false` | Enable binary mask application. |
| `binary_mask_fname` | string | _(empty)_ | Path to the binary mask raster. Any GDAL-readable single-band raster is accepted. Pixels > 0.5 are treated as foreground (keep); all others set the DEM pixel to nodata. Must cover the same extent and have the same cell size and origin as the (possibly reprojected) DEM. Required when `run_binary_mask` is true. |

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

## Step 2d: Pixel replacement

Overwrite individual cells with new values, listed in a point file. Use it to
blank a handful of known-bad pixels, or to correct spot elevations, without
building a mask raster — `run_binary_mask` is the raster-driven equivalent.

Runs **before** breaching, filling and interior nodata filling: it is a
correction to the source data, so everything downstream sees the corrected DEM.
The corrected DEM is only written on its own if you ask for it.

Output: `<write_fname>_pixels_replaced.tif` (with `write_replaced_pixels_raster`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `replace_pixels_fname` | string | _(empty)_ | Point file listing pixels to overwrite. CSV, FlatGeobuf or GeoPackage. Empty disables the step. |
| `replace_pixels_value_field` | string | `value` | Name of the column holding each pixel's new value. |
| `write_replaced_pixels_raster` | bool | `false` | Write the DEM after replacement. Without it the corrected DEM still feeds every later stage but is not written separately. |

**Positions** come from the geometry for vector formats. For CSV, either
`latitude` / `longitude` (read as EPSG:4326) or `easting` / `northing` — `x` and
`y` also work — read as being in the raster's own coordinate system.

**Setting a cell to the raster's nodata value blanks it properly**: the cell is
marked invalid, not merely given the nodata number. That distinction matters,
because a cell holding `-9999` while still flagged as valid data would be treated
as a very deep pit by every routing routine downstream.

Pixels that fall outside the raster, or whose value column cannot be read, are
**skipped with a warning** rather than failing the run — so one bad row in a long
list does not cost you the whole job.

Example CSV:

```
easting,northing,value
-121503,-2502711,-99999
 486210,  4321870,  1204.5
```

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
| `fill_algorithm_choice` | string | `barnes2014` | Filling algorithm. `barnes2014`: Improved Priority-Flood (Barnes et al. 2014) — the default, and the only one that imposes `fill_min_slope` across flats. `vincent1993`: grayscale morphological reconstruction (Vincent 1993) — bit-identical to `barnes2014` and 2–3× faster, but requires `fill_min_slope: 0` (see below). `barnesetal2021`: Fill-Spill-Merge — note this algorithm does not guarantee that all depressions are filled (see doi:10.5194/esurf-9-105-2021). |
| `fill_min_slope` | float | `0.0001` | Minimum gradient imposed between adjacent filled cells to prevent flat regions that produce ambiguous flow directions. Increase slightly for very coarse DEMs. Not used with `barnesetal2021`. See the note below on how this differs from RichDEM. |

### Note: `vincent1993` is faster but cannot impose a gradient

`vincent1993` computes the same fixed point as `barnes2014` by a different
route — two raster scans plus a FIFO drain (O(N)) instead of a priority queue
(O(N log N)). It is verified **bit-identical** to `barnes2014` across all test
fixtures and degenerate cases, and is 2–3× faster (up to ~5× on real terrain
with few depressions).

The catch: it fills depressions to exactly their spill level, leaving **true
flats**. It has no mechanism for the gradient `fill_min_slope` imposes. Since
lsdtt3's D8 router rejects a surface containing true flats, `vincent1993` is
only usable with `fill_min_slope: 0`, and the CLI will refuse to run otherwise
rather than silently ignore the parameter.

Use it when you want a filled surface and are not routing flow across it, or
when a downstream step will resolve the flats itself. For the normal
fill-then-route workflow, stay on `barnes2014`.

Note that `fill_pits()` in the library already dispatches to reconstruction
automatically whenever `min_slope` is zero, so the speed-up is applied without
you asking for it — this parameter is for pinning the backend explicitly.

### Note: `fill_min_slope` is a slope, and differs from RichDEM's "epsilon fill"

lsdtt3 is benchmarked against [RichDEM](https://github.com/r-barnes/richdem)
(`docs/developer_notes/testing_and_benchmarking.md`). Both tools offer what is loosely called
an "epsilon fill", but **they are not the same operation**, and comparing their
outputs without knowing this is misleading.

- **`fill_min_slope: 0`** is equivalent to RichDEM's
  `FillDepressions(epsilon=false)`. Verified bit-exact against it on both a
  high-relief and a low-relief real DEM.
- **`fill_min_slope > 0`** imposes a true **slope**: each cell in a flat is
  raised by `fill_min_slope × cell_size` (× √2 diagonally). On a 30 m grid the
  default gives 3 mm of rise per cell. The gradient is uniform across the
  raster, scales with resolution, and does not depend on elevation or on the
  raster's data type.
- **RichDEM's `epsilon=true`** instead raises each cell to the next
  representable floating-point value (one ULP). That step depends on elevation
  and data type — about 1.2e-4 m at 2000 m in float32, 7.6e-6 m at 100 m, and
  effectively zero near sea level. RichDEM refuses integer rasters for this
  reason.

**What this means in practice.** The two approaches agree closely in high relief
and diverge in low relief, because the flat-gradient term is a negligible share
of the fill in steep terrain and a dominant one in flat terrain. Measured on
256² real DEMs:

| Terrain | Depression fill volume | Flat-gradient term | Its share |
|---|---|---|---|
| Alpine (Val Venosta) | 6091 m | 0.56 m | 0.01 % |
| Lakeland (Finland) | 188 m | 78.8 m | 29.5 % |

In the low-relief case lsdtt3 moves ~40× more material than RichDEM's epsilon
fill. Neither is wrong — they answer the question differently — but if you work
in **low-relief terrain** (coastal plains, glaciated shield, marsh, floodplain)
be aware that `fill_min_slope` is doing a substantial amount of work there, and
consider reducing it. The default is tuned for the steeper landscapes lsdtt3 is
most often applied to.

Setting `fill_min_slope: 0` removes the imposed gradient entirely, but leaves
true flats whose flow directions are ambiguous. That is a reasonable choice only
if the downstream analysis can tolerate them.
| `fill_spill_merge_rainfall_rate` | float | `0.1` | Uniform rainfall depth in elevation units per cell for the Fill-Spill-Merge algorithm. Only used when `fill_algorithm_choice` is `barnesetal2021`. |
| `fsm_boundary_conditions` | string | `oooo` | Boundary conditions for the depression hierarchy used by Fill-Spill-Merge, given as four characters in north–east–south–west order. `o` = open (edge drains freely to the ocean). `n` = no-flux (edge is not an ocean boundary). `p` = periodic (edge wraps; not ocean). Example: `oooo` for a typical DEM tile; `nnnn` for a closed basin. Only used when `fill_algorithm_choice` is `barnesetal2021`. |

---

## Step 7: Depression polygons

Trace closed depressions identified by the fill difference raster into vector polygons. Requires that `write_fill_raster` or `write_fill_diff_raster` is also true (the fill diff must be computed).

Output: `<write_fname>_depressions.<ext>` (extension depends on `write_vector_format`)

Each polygon carries three attributes: `depression_id`, `area_m2`, `volume_m3`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_depression_polygons` | bool | `false` | Enable depression polygon output. Automatically enables `write_fill_diff_raster`. |
| `depression_min_pixels` | int | `1` | Minimum depression size in pixels. Depressions smaller than this are removed before polygon tracing. |
| `depression_min_depth` | float | `0.0` | Minimum fill depth in metres for a cell to be included in a depression. A small positive value (e.g. `0.01`) suppresses near-zero fill artefacts. |
| `write_vector_format` | string | `flatgeobuf` | Output format for all vector writes. `flatgeobuf` or `fgb`: FlatGeobuf (recommended). `gpkg`: GeoPackage. |

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

Output: `<write_fname>_clipped_vector.<ext>` (extension depends on `write_vector_format`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `clip_vector_data_to_raster` | bool | `false` | Enable vector clipping. |
| `vector_data_fname` | string | _(empty)_ | Path to the input vector file. Any OGR-readable format is accepted (Shapefile, GeoPackage, FlatGeobuf, GeoJSON, etc.). Required when `clip_vector_data_to_raster` is true. |

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
min_elevation:            0.0
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

### Run control

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all hardware threads. 1 = single-threaded. Applies to filling, breaching, nodata repair, and wavelet steps. |
| `write_output_format` | `string` | `cog` | `io` | - | Output raster format. | cog: Cloud Optimised GeoTIFF (default). gtiff: plain compressed GeoTIFF. same_as_input: match the input raster format. |

### Projection

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `enforce_projection` | `bool` | `true` | `switch` | - | Check and reproject the raster to a suitable projected CRS. | If true and the raster is geographic or Web Mercator, it is reprojected automatically before any analysis. If false and the raster is not projected, the program exits with an error. |
| `resampling_method` | `string` | `cubic` | `param` | `enforce_projection` | GDAL resampling algorithm used during reprojection. | Cubic is the default here because this workflow reprojects continuous elevation surfaces. Use near or mode for categorical rasters. Accepted values: near, bilinear, cubic (default), cubicspline, lanczos, average, mode, max, min, med, q1, q3, sum. |
| `grid_spacing` | `float` | `-1` | `param` | `enforce_projection` | Target grid spacing in metres after reprojection. | If negative (default), spacing is inferred from the input pixel size and rounded: &gt;=20 m -&gt; nearest 5 m; 5-20 m -&gt; nearest 2 m; &lt;5 m -&gt; nearest 1 m. Set a positive value to enforce a specific spacing. |
| `allow_lambert_projection` | `bool` | `false` | `param` | `enforce_projection` | Use Lambert Azimuthal Equal Area for rasters spanning 3+ UTM zones. | If true and the raster is very wide (&gt;=3 UTM zones), a LAEA projection centred on the raster is used instead of UTM. Output files receive the suffix _lambert. |
| `warp_mem_mb` | `int` | `512` | `param` | `enforce_projection` | GDALWarp working memory in MB. | Larger values reduce source tile re-reads for large rasters. Reduce on memory-constrained systems (minimum ~64 MB). |
| `force_reproject` | `bool` | `false` | `switch` | - | Force reprojection to the recommended analysis CRS. | If true, the raster is reprojected to the recommended CRS even when the input raster is already projected. This allows conversion from an existing projected CRS to the preferred UTM, polar stereographic, or Lambert Azimuthal Equal Area analysis projection. |
| `force_resampling` | `bool` | `false` | `switch` | - | Force a user-specified grid spacing during reprojection. | If true, forced_grid_spacing is used as the target grid spacing instead of the automatically estimated or user-supplied grid_spacing. |
| `forced_grid_spacing` | `float` | `-1` | `param` | `force_resampling` | Target grid spacing in metres when force_resampling is true. | Must be positive when force_resampling is true. Ignored otherwise. |

### Vertical units

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `vertical_units` | `string` | `metres` | `param` | - | Vertical unit of the input elevations: metres, feet or us_survey_feet. | NO raster format records the vertical unit, so this CANNOT be detected and is NEVER applied automatically. Much US state lidar is delivered as NAVD88 FEET on a metre UTM grid; left uncorrected every slope, curvature and relief is too large by 3.28. Establish the vertical datum from the data provider before setting this - a wrong setting silently corrupts correct data. The tool warns when slopes look implausible. |
| `check_vertical_units_against_usgs` | `bool` | `true` | `param` | - | When the slope-plausibility warning fires, query USGS 3DEP to settle whether the elevations are metres or feet. | Plausibility can only FLAG a suspicious DEM; an external reference SETTLES it, because nothing but a unit error puts a DEM 3.28x off a national elevation dataset. Runs ONLY if the warning fires AND the DEM sits in US coverage - which is not the limitation it appears, since feet-on-a-metre-grid is a US delivery convention. Sends 5 pixel coordinates to epqs.nationalmap.gov. Set false to keep the run offline. Nothing is ever converted automatically. |

### Threshold masking

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `min_elevation` | `float` | `0.0` | `param` | - | Pixels at or below this elevation (metres) are set to nodata. | Default of 0.0 removes ocean and below-sea-level artefacts common in global DEMs. Set to a large negative number to disable. |
| `maximum_elevation` | `float` | `30000` | `param` | - | Pixels above this elevation (metres) are set to nodata. | Default covers the full range of terrestrial elevations. Reduce to remove artefacts in low-relief study areas. |

### Binary mask

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_binary_mask` | `bool` | `false` | `switch` | - | Apply an external binary mask raster to the DEM. | Loads the raster at binary_mask_fname and calls raster_mask_by_binary_keep: pixels where the mask &gt; 0.5 are kept; all others are set to nodata. The mask must cover the same spatial extent and have the same cell size and origin as the (possibly reprojected) DEM. Runs after threshold masking. Output suffix: _binary_masked.tif |
| `binary_mask_fname` | `string` | - | `param` | `run_binary_mask` | Full path to the binary mask raster. | Any GDAL-readable single-band raster is accepted. Pixels &gt; 0.5 are treated as keep (foreground); all others set the DEM pixel to nodata. Required when run_binary_mask = true. |

### Interior nodata filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_filled_nodata_raster` | `bool` | `false` | `switch` | - | Write a DEM with interior nodata pixels filled. | Only fills nodata pixels completely surrounded by valid data. Does not fill nodata that touches the raster edge. |
| `nodata_fill_method_choice` | `string` | `laplacian` | `param` | `write_filled_nodata_raster` | Method for filling interior nodata pixels. | laplacian: smooth Dirichlet interpolation (recommended). iterative_weighted_average: fast, well-suited to thin voids. |
| `nodata_fill_max_iter` | `int` | `10000` | `param` | `write_filled_nodata_raster` | Maximum iterations for iterative nodata fill methods. | Only used when nodata_fill_method_choice is laplacian or iterative_weighted_average. |

### Clip to rectangle

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_clipped_rectangle_raster` | `bool` | `false` | `switch` | - | Write a DEM clipped to the largest rectangle of all-valid pixels. | Run write_filled_nodata_raster first to maximise the rectangle size. |

### Depression breaching

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_breached_raster` | `bool` | `false` | `switch` | - | Write a depression-breached DEM. | Use breach_algorithm_choice to select the algorithm. |
| `write_carve_diff_raster` | `bool` | `false` | `switch` | - | Write the difference raster: breached minus original. | Negative values show where breaching lowered the surface. |
| `breach_before_fill` | `bool` | `false` | `switch` | - | Breach depressions before filling pits. | Recommended for lidar DEMs where channelised flow paths must be carved before residual pits are filled. |
| `fill_remainder_after_breach` | `bool` | `true` | `param` | `breach_before_fill` | Fill pits that could not be fully breached. | Only used when breach_algorithm_choice is least_cost. If cells without a downslope neighbour remain after the breach, they are filled with fill_min_slope and fill_flat_resolution_choice, as the barnes2014 pit fill would, so the breached DEM is routable. |
| `breach_algorithm_choice` | `string` | `lindsay2016` | `param` | `write_breached_raster` | Algorithm for depression breaching. | lindsay2016: Lindsay (2016) efficient breaching. least_cost: Lindsay (2020) pit-centric least-cost breaching. |
| `breach_lc_max_dist` | `int` | `1000 (outlet) / 100 (whitebox)` | `param` | `write_breached_raster` | Maximum breach path length in cells for least-cost breaching. | Counted in steps from the pit (a diagonal step counts as one), not as a radius. Pits that need a longer path are left unbreached (and filled if fill_remainder_after_breach is true). 1 to 32767. If not set, 1000 with breach_lc_edge_policy outlet (long enough to bridge large water bodies, so a far-off pit is not carved to the nearest raster edge instead) and 100 with whitebox (the WhiteboxTools default; longer paths there can carve across divides). Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_max_cost` | `float` | `unlimited` | `param` | `write_breached_raster` | Maximum cost of a least-cost breach path. | The cost of a path is the depth each cell on it is lowered, summed: metres when breach_lc_minimize_dist is false; with breach_lc_minimize_dist true each step is weighted by its length in cells (1 or 1.414), so the units are metres x cells. A pit whose cheapest path costs more is left unbreached. Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_minimize_dist` | `bool` | `false` | `param` | `write_breached_raster` | Weight least-cost breach paths by their length. | false: a path costs the total depth lowered (m), so the search prefers shallow cuts however long. true: each step is weighted by its length in cells, favouring shorter paths; breach_lc_max_cost is then in metres x cells. Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_edge_policy` | `string` | `outlet` | `param` | `write_breached_raster` | Whether a least-cost breach may drain off the raster edge. | outlet: the raster edge is an outlet, like a nodata cell, so a pit near the edge can breach straight off the grid. whitebox: the edge is not an outlet (as in WhiteboxTools), so a breach must reach an edge cell that is already low enough, which can carve deep cuts across divides near the edge. Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_flat_increment` | `float` | `0.0` | `param` | `write_breached_raster` | Elevation drop (m) per cell along a least-cost breach channel. | 0 or less selects the automatic value: twice the float32 spacing at the largest absolute elevation in the DEM (e.g. 2.4e-4 m for elevations of 1024-2047 m), the smallest drop that keeps the channel strictly downhill. Only used when breach_algorithm_choice is least_cost. |

### Smoothing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `smooth_algorithm_choice` | `string` | `none` | `param` | - | Edge-preserving smoothing applied BEFORE breaching and filling: none, perona_malik or coherence_diffusion. | For work that reads heights off surfaces - terrace treads, valley floors, anything the fluvial-surface ladder bands - the riser is the signal, not the noise. A mean or Gaussian filter removes DEM noise and the tread/riser break in equal measure, which moves every height break downstream. Perona-Malik (1990) diffuses within a facet and stops at a gradient above K, so the riser survives. It is a CONDITIONING step: the smoothed raster is written as its own file and everything downstream should read that, so which DEM an analysis ran on is visible rather than inferred. |
| `perona_malik_iterations` | `int` | `10` | `param` | - | Explicit diffusion timesteps. | More iterations remove more within-facet noise. The riser is preserved at every count, so this trades runtime against how quiet the treads become rather than against how much structure survives. |
| `perona_malik_lambda` | `float` | `0.20` | `param` | - | Diffusion timestep. Must be in (0, 0.25]. | 0.25 is the stability limit of the explicit 4-neighbour stencil. A value above it is REJECTED rather than warned about: an unstable run does not give a worse DEM, it gives a meaningless one. |
| `perona_malik_conductance` | `string` | `exponential` | `param` | - | Conductance function: exponential or quadratic. | exponential is exp(-(s/K)^2), Perona &amp; Malik eq. 3, and falls away faster above K, so it preserves a sharp riser more strictly. quadratic is 1/(1+(s/K)^2), eq. 4, which keeps diffusing weakly across strong edges and is the gentler choice on a noisy DEM. |
| `perona_malik_k_quantile` | `float` | `0.90` | `param` | - | Quantile of the DEM's own neighbour-difference distribution used to set the edge threshold K. 0 = state K directly instead. | K is in metres of rise between adjacent cells, the same units as the differences the stencil takes - not a central-difference gradient, which averages across two cells and so cannot see the one-cell noise the diffusion does. A K carried over from another landscape means nothing, so following Black et al. (1998) it is read off the distribution of the DEM in hand. The value actually used is printed. |
| `perona_malik_k_metres_per_cell` | `float` | `0` | `param` | - | Edge threshold K stated directly, in metres of rise per cell. Only read when perona_malik_k_quantile is 0. | Use this to pin K across a set of DEMs that must be conditioned identically; the quantile would otherwise give each one its own. |
| `write_stencil_difference_raster` | `bool` | `false` | `output` | - | Write the neighbour-difference raster (metres) K was read from. | Per cell, the largest absolute elevation difference to a valid 4-neighbour - the exact quantity the conductance is evaluated on. Write it the first time a DEM is conditioned: a quantile is only interpretable against the distribution it was taken from. |
| `coherence_iterations` | `int` | `10` | `param` | - | Coherence-enhancing diffusion steps. | More iterations bridge longer breaks along a feature and smooth further within a facet. |
| `coherence_time_step` | `float` | `0.15` | `param` | - | Diffusion timestep tau. Must be in (0, 0.25]. | Same explicit-stencil stability limit as perona_malik_lambda, and rejected the same way. |
| `coherence_derivative_radius_metres` | `float` | `0` | `param` | - | Gradient pre-smoothing radius before the structure tensor is formed. 0 = one cell. | Defaults are stated in CELLS, not metres, because a radius that is sensible on 1 m lidar is sub-pixel on 30 m COP30 and the structure tensor then has nothing to integrate over. |
| `coherence_integration_radius_metres` | `float` | `0` | `param` | - | Structure-tensor integration radius. 0 = three cells. | This sets the length scale over which a feature has to look linear before the filter will smooth along it. Too small and every noisy streak reads as coherent; too large and a curving riser stops registering as one feature. See the note on cells vs metres under coherence_derivative_radius_metres. |
| `coherence_contrast` | `float` | `0.2` | `param` | - | Coherence contrast C: how selective the along-feature smoothing is. | The coherence it is compared against is normalised to [0, 1], so C is scale-invariant and transfers between DEMs. Raise it to smooth only the most coherent threads; lower it to smooth more of the field. |
| `coherence_min_diffusivity` | `float` | `0.001` | `param` | - | Across-feature diffusivity floor, in (0, 1). | The eigenvalue perpendicular to the feature. Near zero means no smoothing at all across a riser, which is the point; raise it only to deliberately soften the across-feature behaviour. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_fill_raster` | `bool` | `false` | `switch` | - | Write a pit-filled DEM. | Use fill_algorithm_choice to select the algorithm. |
| `write_fill_diff_raster` | `bool` | `false` | `switch` | - | Write the difference raster: filled minus original. | Positive values show where filling raised the surface.If fill_agorithm_choice is barnesetal2021, this represents the water depth standing in topographic depressions. |
| `priority_flood_queue` | `string` | `binary` | `param` | `write_fill_raster` | Spill-queue backend for Barnes 2014 filling. | binary (default) uses the historic heap. radix uses the equivalent hybrid radix queue. Only used when fill_algorithm_choice is barnes2014. |
| `fill_min_slope` | `float` | `0.0001` | `param` | `write_fill_raster` | Minimum gradient imposed between adjacent filled cells. | Prevents flat regions that cause ambiguous flow direction. Increase slightly for very coarse DEMs. This does not apply to fill-spill-merge. Under fill_flat_resolution_choice barnes2014b_gradient it is the minimum gradient across flats; barnes2014b_float_step ignores it. |
| `fill_flat_resolution_choice` | `string` | `inline_gradient` | `param` | - | How flats are made drainable: inline_gradient, barnes2014b_gradient or barnes2014b_float_step. | inline_gradient (default, unchanged): the Priority-Flood fill imposes fill_min_slope during the flood. The barnes2014b options make flow converge across flats (Barnes et al. 2014, doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses fill_min_slope (m/m, must be &gt; 0) as the minimum gradient; flats are raised about 3-4x more, capped below their rims, and capped flats are reported. barnes2014b_float_step is equivalent to RichDEM ResolveFlats on float32 in its first pass: it ignores fill_min_slope and its slopes depend on elevation (of order ULP(z) / cell size). With either barnes2014b option the result is always routable: cells the first pass leaves without a lower neighbour are refilled and resolved again (up to 4 rounds), and any still left are raised by a final priority-flood pass by one float32 step each; the report gives the rounds and the cells changed after the first pass. |
| `fill_spill_merge_rainfall_rate` | `float` | `0.1` | `param` | `write_fill_raster` | Uniform rainfall depth in elevation units per cell for fill-spill-merge flow routing (Barnes et al., 2021). | Spatial patterns in rainfall to route into and over depressions are not currently supported. |
| `fsm_boundary_conditions` | `string` | `oooo` | `param` | `write_fill_raster` | Boundary conditions for the depression hierarchy (north east south west). | Four characters, one per edge: o = open (cells on this edge drain freely to the ocean, default). b = forced baselevel outlet (same as o for hierarchy purposes). n = no-flux (edge is NOT an ocean boundary). p = periodic (edge wraps; NOT an ocean boundary). Example: oooo = all four edges are open ocean (typical for a DEM tile). Example: nnnn = closed basin (no edge is ocean; requires interior pits). |

### Filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_algorithm_choice` | `string` | `barnes2014` | `param` | - | Depression filling algorithm. | barnes2014: Improved Priority-Flood (Barnes et al. 2014). The default. Imposes fill_min_slope across flats so the filled surface is routable. vincent1993: grayscale morphological reconstruction (Vincent 1993). Produces a BIT-IDENTICAL result to barnes2014 and is 2-3x faster, but fills to exactly the spill level and so requires fill_min_slope: 0 -- it cannot impose a gradient across flats, and lsdtt3's D8 router will reject a surface with true flats. barnesetal2021: Fill-Spill-Merge (https://doi.org/10.5194/esurf-9-105-2021); note this distributes a finite rainfall volume and does not guarantee all depressions are filled. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). CSV is available only through selectors for point/table outputs. Per-output overrides exist only where a container genuinely differs. Raster outputs are governed separately by write_output_format. |

### Pixel replacement

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `replace_pixels_fname` | `string` | - | `param` | - | Point file listing pixels to overwrite with new values. | CSV, FlatGeobuf or GeoPackage. Positions come from geometry, or from latitude/longitude, or from easting/northing in the raster CRS. The replacement value is read from the column named by replace_pixels_value_field. Setting a cell to the raster's nodata value marks it nodata rather than merely writing the number. Applied before breaching, filling and nodata filling. |
| `write_replaced_pixels_raster` | `bool` | `false` | `switch` | `replace_pixels_fname` | Write the DEM after pixel replacement. | Without this the corrected DEM is still used by every later stage, but is not written on its own. Set it when the corrected DEM is the thing you want. |
| `replace_pixels_value_field` | `string` | `value` | `param` | `replace_pixels_fname` | Column holding the replacement value. | Named field in replace_pixels_fname giving the new value for each listed pixel. |

### Sink candidates

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_sink_candidates` | `bool` | `false` | `switch` | - | Write a point layer of candidate internally drained sinks. | Finds closed depressions that never drain to the DEM edge and writes them for review in a GIS. Pass the reviewed file to lsdtt-internal-drainage to route flow into the accepted sinks. Output format follows write_vector_format. |
| `sink_candidate_min_pixels` | `int` | `1` | `param` | `write_sink_candidates` | Minimum sink size in pixels. | Closed depressions smaller than this are not reported. Raise it to suppress DEM noise on rough surfaces. |

### Depression polygons

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_depression_polygons` | `bool` | `false` | `switch` | - | Write a vector polygon layer of topographic depressions, based on the fill diff raster. | Each polygon carries depression_id, area_m2, and volume_m3 attributes. Output format is controlled by write_vector_format. |
| `depression_min_pixels` | `int` | `1` | `param` | `write_depression_polygons` | Minimum depression size in pixels. | Depressions with fewer pixels than this are removed before polygon tracing. Useful for suppressing single-cell fill artefacts. |
| `depression_min_depth` | `float` | `0.0` | `param` | `write_depression_polygons` | Minimum fill depth in metres for a cell to be considered a depression. | Cells where fill_diff &lt;= this value are treated as background. A small positive value (e.g. 0.01) suppresses near-zero fill artefacts. |

### Hillshade

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_hillshade_raster` | `bool` | `false` | `switch` | - | Write a hillshade raster. | Uses the Zevenbergen and Thorne (1987) algorithm. Computed from the filled DEM when write_fill_raster is also true. |
| `hillshade_azimuth` | `float` | `315` | `param` | `write_hillshade_raster` | Sun azimuth in degrees clockwise from north. | 315 (northwest) is the cartographic convention. |
| `hillshade_altitude` | `float` | `45` | `param` | `write_hillshade_raster` | Sun altitude above the horizon in degrees. | 45 degrees is a typical default. |
| `hillshade_z_factor` | `float` | `1.0` | `param` | `write_hillshade_raster` | Vertical exaggeration factor applied before computing hillshade. | Values greater than 1 emphasise relief in flat terrain. |
| `write_multidirectional_hillshade_raster` | `bool` | `false` | `switch` | - | Write a multidirectional hillshade raster. | Computes hillshade from multiple azimuths and writes the mean illumination raster to reduce directional bias. |
| `multidirectional_hillshade_n_directions` | `int` | `8` | `param` | `write_multidirectional_hillshade_raster` | Number of evenly spaced azimuth directions used for multidirectional hillshade. | Must be &gt;= 2. Typical values are 8 to 16. |
| `multidirectional_hillshade_start_azimuth` | `float` | `0.0` | `param` | `write_multidirectional_hillshade_raster` | Starting azimuth in degrees clockwise from north for multidirectional hillshade. | Subsequent azimuths are sampled uniformly around 360 degrees. |

### Mosaic

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `mosaic_file_list` | `string` | - | `switch` | - | Path to a text file listing rasters to mosaic (one path per line). | If set, the listed rasters are merged into a single output raster before any other preprocessing step. Rasters with a different CRS, cell size, or pixel-misaligned origin are skipped with a warning. Blend semantics: last valid pixel wins (later files take priority over earlier ones in overlap regions, matching gdal_merge.py default). The mosaic is written to write_path/write_fname_mosaic.tif and used as the input DEM for all subsequent pipeline steps. Lines starting with # and blank lines are ignored in the file list. |

### Vector clipping

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `clip_vector_data_to_raster` | `bool` | `false` | `switch` | - | Clip a vector file to the bounding box of the output raster. | The vector CRS need not match the raster CRS - the raster bbox is reprojected into the vector's native CRS before clipping. Set vector_data_fname to the path of the input vector file. |
| `vector_data_fname` | `string` | - | `param` | `clip_vector_data_to_raster` | Path to the input vector file to clip. | Any OGR-readable format is accepted (Shapefile, GeoPackage, FlatGeobuf, GeoJSON, etc.). Required when clip_vector_data_to_raster is true. |

### Clip raster to polygon

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `clip_raster_to_polygon` | `bool` | `false` | `switch` | - | Crop the DEM to the minimum bounding box of a polygon. | The DEM is cropped to the smallest axis-aligned rectangle that contains the polygon, in the raster's own CRS. Pixels outside the polygon are also set to nodata unless force_rectangular_raster is true. If the file contains more than one polygon, a warning is issued and only the first polygon is used. Set polygon_clip_fname to the path of the polygon file. Output suffix: _polygon_clipped.tif |
| `polygon_clip_fname` | `string` | - | `param` | `clip_raster_to_polygon` | Path to the polygon vector file used to clip the DEM. | Any OGR-readable polygon format is accepted (GeoPackage, FlatGeobuf, Shapefile, GeoJSON, etc.). The first layer and first polygon are used. Required when clip_raster_to_polygon is true. |
| `force_rectangular_raster` | `bool` | `false` | `param` | `clip_raster_to_polygon` | When clipping to a polygon, keep all data inside the bounding box. | If false (default), pixels that fall inside the bounding box but outside the polygon shape are set to nodata. If true, all pixels within the bounding box are kept regardless of whether they fall inside the polygon - the output is a plain rectangular crop. Only used when clip_raster_to_polygon is true. |
| `polygon_clip_rule` | `string` | `majority_coverage` | `param` | `clip_raster_to_polygon` | Pixel inclusion rule when masking by the polygon shape. | Only used when clip_raster_to_polygon is true and force_rectangular_raster is false. any_coverage: pixel is kept if the polygon touches any part of it. majority_coverage: pixel is kept if the polygon centre falls inside (approximately 50% coverage, default). full_coverage: pixel is kept only if the polygon covers it fully. |

### Point sampling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `sample_rasters_to_points` | `bool` | `false` | `switch` | - | Sample a list of rasters onto a point dataset. | Reads a point dataset and a CSV file with columns raster,attribute. Each raster is sampled at each point and written as a new attribute in an output point file. |
| `point_data_fname` | `string` | - | `param` | `sample_rasters_to_points` | Path to the input point dataset. | Recommended input is an OGR-readable point vector file such as FlatGeobuf, GeoPackage, Shapefile, or GeoJSON. CSV latitude/longitude support can also be added, but vector input is the cleanest first implementation because existing attributes are preserved. |
| `point_data_layer_name` | `string` | - | `param` | `sample_rasters_to_points` | Optional layer name for the input point dataset. | If empty, the first layer is used, or the only layer if the file has one. |
| `raster_attribute_fname` | `string` | - | `param` | `sample_rasters_to_points` | CSV file listing rasters to sample, with columns raster,attribute. | The raster column gives the raster filename or path. The attribute column gives the name of the output attribute to append to the point file. |

<!-- END GENERATED CLI PARAMETERS -->
