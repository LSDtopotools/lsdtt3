# CLI reference: lsdtt-raster-calculator

Two-raster map algebra calculator. Performs element-wise arithmetic operations (addition, subtraction, multiplication, division) on pairs of rasters. Handles automatic reprojection, coregistration, and nodata harmonisation.

## Usage

```
lsdtt-raster-calculator  <path> <parameter_file>
lsdtt-raster-calculator  <parameter_file>        # path defaults to ./
lsdtt-raster-calculator  -h                      # print help and exit
lsdtt-raster-calculator  -v                      # print version and exit
lsdtt-raster-calculator  -c                      # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

Keys use `snake_case`. All keys and their defaults are listed below.

---

## Pipeline overview

Steps run in this order:

| Step | Description |
|---|---|
| 1. Load first raster | Read primary input DEM or raster |
| 2. Load second raster | Read secondary raster |
| 3. Reproject second raster | If CRS mismatch, warp to first raster grid |
| 4. Coregister *(optional)* | Align second raster to first raster if requested |
| 5. Harmonise nodata | Ensure consistent nodata values and mask |
| 6. Compute overlap mask | Identify valid pixels in both rasters |
| 7. Calculate outputs | Perform requested algebra operations |
| 8. Write outputs | Save result rasters |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input rasters. |
| `read_fname` | string | parameter file stem | First raster filename. Any GDAL-supported format. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `second_raster_read_fname` | string | _(required)_ | Filename of second raster. |
| `second_raster_read_path` | string | same as `read_path` | Directory containing second raster. |

---

## Reprojection

| Parameter | Type | Default | Description |
|---|---|---|---|
| `reproject_second_raster_to_first` | bool | `true` | Automatically reproject second raster to first raster grid if CRS differs. |
| `resampling_method` | string | `bilinear` | GDAL method for warping the second raster. Bilinear avoids cubic overshoot before raster algebra; use `near` for categorical values. |

---

## Coregistration

| Parameter | Type | Default | Description |
|---|---|---|---|
| `coregister_second_raster` | bool | `false` | Coregister second raster to first using cross-correlation or other methods. |
| `coregistration_window_pixels` | int | `64` | Window size (pixels) for coregistration search. |
| `coregistration_search_distance_pixels` | int | `10` | Search distance (pixels) for offset detection. |

---

## Operations

| Parameter | Type | Default | Description |
|---|---|---|---|
| `calculate_sum` | bool | `false` | Output: first + second |
| `calculate_difference` | bool | `true` | Output: first - second |
| `calculate_absolute_difference` | bool | `false` | Output: abs(first - second) |
| `calculate_multiply` | bool | `false` | Output: first × second |
| `calculate_divide` | bool | `false` | Output: first ÷ second |
| `calculate_ratio` | bool | `false` | Output: first / second (alias for divide) |
| `calculate_first_on_overlap` | bool | `false` | Output: first raster where both are valid |
| `calculate_second_on_overlap` | bool | `false` | Output: second raster where both are valid |

---

## Output handling

| Parameter | Type | Default | Description |
|---|---|---|---|
| `set_nodata_to_zero` | bool | `false` | Replace nodata with `0` in output (careful: may introduce false signal). |
| `preserve_first_raster_nodata` | bool | `true` | Use first raster's nodata value for outputs. |

---

## Output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `cog` | Raster output format: `cog` (Cloud-Optimised GeoTIFF, preferred), `gtiff`, `same_as_input`. |

---

## Output filenames

Operations produce outputs named according to the operation:

- `<write_fname>_sum.tif`
- `<write_fname>_difference.tif`
- `<write_fname>_absolute_difference.tif`
- `<write_fname>_product.tif`
- `<write_fname>_ratio.tif`
- `<write_fname>_first_on_overlap.tif`
- `<write_fname>_second_on_overlap.tif`
