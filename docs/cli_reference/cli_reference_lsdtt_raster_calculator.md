# CLI reference: lsdtt-raster-calculator

Two-raster algebra: sum, product, ratio, difference, absolute difference, and
each raster masked to where both are valid. A typical use is a DEM of
difference. The second raster is warped onto the first raster's grid when the
two differ, and can optionally be coregistered to the first (Nuth and Kääb
2011) before any algebra.

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

Keys use `snake_case`. All keys and their defaults are listed under
[Parameters](#parameters) below.

---

## Pipeline

1. Load the first raster (`read_fname`) and the second (`second_raster_read_fname`).
2. If the grids differ, warp the second raster onto the first raster's grid
   with `resampling_method` (default `bilinear`; use `near` for categorical
   values). With `resample_second_to_first_grid: false` the tool stops instead.
3. Give the second raster the first raster's nodata value.
4. Optionally coregister the second raster to the first
   (`run_dem_coregistration: true`).
5. Build the overlap mask: pixels valid in both rasters.
6. Compute and write the selected outputs. Pixels outside the overlap are
   nodata, and so is a division by zero in the ratio.

## Example: a DEM of difference

```
read_fname: dem_2020.tif
second_raster_read_fname: dem_2010.tif
write_difference_raster: true
run_dem_coregistration: true
write_dem_coreg_diagnostics_csv: true
```

`write_difference_raster` is on by default; the other outputs are off.

## Outputs

`<pair>` below is `<write_fname>_and_<second_raster_write_fname>`, with
`_coreg` appended when `run_dem_coregistration: true`.
`second_raster_write_fname` defaults to the stem of
`second_raster_read_fname`.

| File | Written when |
|---|---|
| `<pair>_difference.tif` (first − second) | `write_difference_raster: true` (default) |
| `<pair>_absolute_difference.tif` | `write_absolute_difference_raster: true` |
| `<pair>_sum.tif` | `write_sum_raster: true` |
| `<pair>_multiply.tif` | `write_multiply_raster: true` |
| `<pair>_ratio.tif` (first / second) | `write_ratio_raster: true` |
| `<pair>_first_on_overlap.tif`, `<pair>_second_on_overlap.tif` | `write_first_on_overlap_raster` / `write_second_on_overlap_raster` |
| `<pair>_overlap_mask.tif` (1 = both valid) | `write_overlap_mask_raster: true` |
| `<second_raster_write_fname>_resampled_to_<write_fname>.tif` | whenever the second raster had to be warped |
| `<second_raster_write_fname>_coregistered_to_<write_fname>.tif` | `write_second_coregistered_raster: true` |
| `<pair>_diagnostics.csv` (`<pair>` ends in `_coreg` here) | `write_dem_coreg_diagnostics_csv: true` |

The output extension follows `write_output_format` (default `cog`).

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

### Second raster and alignment

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `second_raster_read_fname` | `string` | - | `param` | - | Second raster filename including extension. | Required. This is the right-hand raster operand in all two-raster operations. |
| `second_raster_read_path` | `string` | - | `param` | - | Optional directory for second_raster_read_fname. | If empty, read_path is used. If second_raster_read_fname is an absolute path, this parameter is ignored. |
| `second_raster_write_fname` | `string` | - | `param` | - | Output stem used when naming second-raster-derived outputs. | If empty, defaults to second_raster_read_fname stem (extension removed). Used in output naming such as &lt;second_raster_write_fname&gt;_resampled_to_&lt;write_fname&gt;. |
| `resample_second_to_first_grid` | `bool` | `true` | `param` | - | Warp second raster to first raster grid when needed. | If true, second raster is resampled and aligned to read_fname grid before algebra. If false, the CLI errors when grids differ. |
| `resampling_method` | `string` | `bilinear` | `param` | - | GDAL resampling method used when warping second raster. | Bilinear is the conservative default for aligning an arbitrary second raster because it does not introduce cubic overshoot before raster algebra. Use near for categorical values. Allowed values follow GDAL warp methods: near, bilinear, cubic, cubicspline, lanczos, average, mode, max, min, med, q1, q3, sum. |
| `verbose_alignment_report` | `bool` | `true` | `param` | - | Print grid-compatibility and overlap diagnostics. | If true, reports whether rasters already share a grid and whether resampling was required. |

### DEM coregistration

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_dem_coregistration` | `bool` | `false` | `switch` | - | Enable DEM coregistration before algebra outputs. | If true, the second raster is coregistered to the first raster after optional resampling and before all algebra operations. |
| `dem_coregistration_method` | `string` | `nuth_kaab` | `param` | `run_dem_coregistration` | DEM coregistration method. | Allowed values: nuth_kaab. Method reference: Nuth and Kaab (2011), https://doi.org/10.5194/tc-5-271-2011. This implementation was ported from xdem. |
| `write_second_coregistered_raster` | `bool` | `false` | `switch` | `run_dem_coregistration` | Write the coregistered second raster. | Output file: &lt;second_raster_write_fname&gt;_coregistered_to_&lt;write_fname&gt;. |
| `write_dem_coreg_diagnostics_csv` | `bool` | `false` | `switch` | `run_dem_coregistration` | Write per-iteration DEM coregistration diagnostics CSV. | Output file: &lt;write_fname&gt;_and_&lt;second_raster_write_fname&gt;_coreg_diagnostics.csv. |
| `nuth_kaab_coreg_max_iterations` | `int` | `8` | `param` | `run_dem_coregistration` | Maximum iterations for nuth_kaab DEM coregistration. | Used only when dem_coregistration_method is nuth_kaab. |
| `nuth_kaab_coreg_offset_threshold_pixels` | `float` | `0.001` | `param` | `run_dem_coregistration` | Convergence threshold in pixels for horizontal offset updates. | Matches Nuth and Kaab style stopping criterion: stop when estimated horizontal residual offset per iteration falls below this value. |
| `nuth_kaab_coreg_slope_min_degrees` | `float` | `0.0` | `param` | `run_dem_coregistration` | Minimum slope [degrees] used in nuth_kaab sample selection. | Used only when dem_coregistration_method is nuth_kaab. |
| `nuth_kaab_coreg_slope_max_degrees` | `float` | `90.0` | `param` | `run_dem_coregistration` | Maximum slope [degrees] used in nuth_kaab sample selection. | Used only when dem_coregistration_method is nuth_kaab. |
| `nuth_kaab_coreg_aspect_bin_count` | `int` | `72` | `param` | `run_dem_coregistration` | Number of aspect bins used before cosine fitting. | Nuth and Kaab fitting is performed on binned dh/tan(slope) by aspect. |
| `nuth_kaab_coreg_min_samples` | `int` | `500` | `param` | `run_dem_coregistration` | Minimum samples required after filtering. | Used only when dem_coregistration_method is nuth_kaab. |
| `nuth_kaab_coreg_subsample_count` | `int` | `500000` | `param` | `run_dem_coregistration` | Maximum number of eligible pixels used for iterative fitting. | Set -1 to use all eligible pixels; otherwise random subsample up to this count. |
| `nuth_kaab_coreg_apply_vertical_shift` | `bool` | `true` | `param` | `run_dem_coregistration` | Allow vertical shift term in nuth_kaab fitting. | If false, only horizontal shift is estimated. |
| `nuth_kaab_coreg_max_horizontal_step_metres` | `float` | `0.0` | `param` | `run_dem_coregistration` | Maximum horizontal step per iteration [m]. | &lt;= 0 disables per-iteration step capping. |
| `nuth_kaab_coreg_initial_dx_metres` | `float` | `0.0` | `param` | `run_dem_coregistration` | Initial easting shift guess [m] applied before iterations. | Useful to seed the optimizer when large known offsets are present. |
| `nuth_kaab_coreg_initial_dy_metres` | `float` | `0.0` | `param` | `run_dem_coregistration` | Initial northing shift guess [m] applied before iterations. | Useful to seed the optimizer when large known offsets are present. |
| `nuth_kaab_coreg_initial_dz_metres` | `float` | `0.0` | `param` | `run_dem_coregistration` | Initial vertical shift guess [m]. | If vertical-shift solving is enabled, iterations overwrite this with the final estimate. |

### Output format and diagnostics

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_output_format` | `string` | `cog` | `io` | - | Output raster format. | Allowed values: cog, gtiff, same_as_input, png. |
| `write_overlap_mask_raster` | `bool` | `false` | `switch` | - | Write overlap mask raster. | Writes &lt;write_fname&gt;_and_&lt;second_raster_write_fname&gt;[_coreg]_overlap_mask where _coreg is added when run_dem_coregistration=true; 1 means both rasters are valid and 0 means at least one raster is nodata. |

### Selected outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_sum_raster` | `bool` | `false` | `switch` | - | Write sum raster. | Output file: &lt;write_fname&gt;_and_&lt;second_raster_write_fname&gt;[_coreg]_sum. |
| `write_multiply_raster` | `bool` | `false` | `switch` | - | Write multiply raster. | Output file: &lt;write_fname&gt;_and_&lt;second_raster_write_fname&gt;[_coreg]_multiply. |
| `write_ratio_raster` | `bool` | `false` | `switch` | - | Write ratio raster (first / second). | Output file: &lt;write_fname&gt;_and_&lt;second_raster_write_fname&gt;[_coreg]_ratio. Division by zero is set to nodata. |
| `write_difference_raster` | `bool` | `true` | `switch` | - | Write difference raster (first - second). | Output file: &lt;write_fname&gt;_and_&lt;second_raster_write_fname&gt;[_coreg]_difference. |
| `write_absolute_difference_raster` | `bool` | `false` | `switch` | - | Write absolute difference raster. | Output file: &lt;write_fname&gt;_and_&lt;second_raster_write_fname&gt;[_coreg]_absolute_difference, computed as abs(first - second). |
| `write_first_on_overlap_raster` | `bool` | `false` | `switch` | - | Write first raster where overlap with second is valid. | Output file: &lt;write_fname&gt;_and_&lt;second_raster_write_fname&gt;[_coreg]_first_on_overlap. |
| `write_second_on_overlap_raster` | `bool` | `false` | `switch` | - | Write second raster where overlap with first is valid. | Output file: &lt;write_fname&gt;_and_&lt;second_raster_write_fname&gt;[_coreg]_second_on_overlap. |

<!-- END GENERATED CLI PARAMETERS -->
