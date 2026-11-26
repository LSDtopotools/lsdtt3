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
   With `reproject_to_metric_grid: true`, a first raster in degrees defines a UTM
   grid, and each DEM is warped onto it once, from its native grid.
2. If the grids or CRSs differ, warp the second raster onto the first raster's
   grid with `resampling_method`. The result is Float32, so an interpolated Int16
   DEM keeps its fractional heights. Cells whose resampling kernel touched
   nodata or left the source are nodata. With
   `resample_second_to_first_grid: false` the tool stops instead.
3. Convert both rasters to metres (`vertical_units`,
   `second_raster_vertical_units`) and, if declared, put the second on the first
   raster's geoid.
4. Optionally coregister the second raster to the first
   (`run_dem_coregistration: true`). The shift is fitted, then the original second
   raster is warped again with that shift, so it is interpolated only once.
5. Build the overlap mask: pixels valid in both rasters. No valid overlap stops
   the run.
6. Report the median and NMAD of first minus second, over stable terrain if a mask
   was given, and write the selected outputs. Pixels outside the overlap are
   nodata, and so is a division by zero in the ratio.

## Grids that do not line up

**Resampling.** Interpolating one DEM of a pair changes it, and the difference
then contains false change. `resampling_method: auto` (the default) uses `cubic`,
or `average` when the target cells are at least 1.5 times wider than the source
cells. Bilinear smooths the interpolated DEM, so ridges appear to lose height and
valleys to gain it. GDAL's `lanczos` can shift the result by a source cell when
downsampling, so use it only between grids of equal resolution.

**Which raster comes first.** The second raster is the one interpolated. If the
two differ in resolution, make the coarser DEM the first raster, so the finer one
is averaged rather than the coarser one being interpolated to detail it never had.

### Global 30 m DEMs

| DEM | Geoid | Cell centres at low latitude | Surface |
|---|---|---|---|
| SRTM, NASADEM | `egm96` | On whole arc-seconds | C-band radar, February 2000 |
| ALOS AW3D30 | `egm96` | Half a cell from whole arc-seconds | Optical stereo |
| Copernicus GLO-30 | `egm2008` | On whole arc-seconds | X-band radar |

- AW3D30 is therefore half a cell from the other two. This is handled by the warp.
- Above 50 degrees latitude GLO-30 cells are wider east-west than north-south
  (1.5 x 1 arc-seconds at first, and wider towards the poles).
- EGM96 and EGM2008 differ smoothly across a DEM, so a constant shift cannot remove
  the difference. Declare both geoids. This needs PROJ's geoid grids. If they are
  missing the run stops and names the command that installs them
  (`projsync --file us_nga_egm96_15 --file us_nga_egm08_25`), or set the environment
  variable `PROJ_NETWORK=ON`.
- These are surface models from different sensors. Radar penetration into snow and
  canopy, and vegetation growth, remain in the difference after coregistration, and
  a nonzero median over stable ground is common. Acquisition dates vary within a
  product, so check them from each product's metadata before reading a difference
  as change between two dates.

```
read_fname: nasadem.tif
second_raster_read_fname: cop30.tif
first_raster_vertical_datum: egm96
second_raster_vertical_datum: egm2008
reproject_to_metric_grid: true
run_dem_coregistration: true
nuth_kaab_coreg_slope_min_degrees: 5
dem_coregistration_stable_mask_fname: stable_terrain.tif
write_difference_raster: true
```

### Gridded lidar

- Two lidar DEMs on offset grids need only the default resampling.
- Lidar is often in feet on a metre grid, and no raster format records this. Set
  `vertical_units` and `second_raster_vertical_units`.
- Check that the surveys share a vertical datum and geoid model. `first_raster_vertical_datum`
  covers only EGM96 and EGM2008. A difference between other geoid models is not
  corrected.

### Stable terrain

Coregistration fits a shift by assuming most of the ground did not change. Where
landslides, lakes or ice cover much of a DEM that fails, so give
`dem_coregistration_stable_mask_fname`, with non-zero values on stable ground. The
reported NMAD over stable terrain is a measure of the noise in the difference.

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
| `resampling_method` | `string` | `auto` | `param` | - | GDAL resampling method used when warping a raster onto the analysis grid. | auto uses cubic, or average when the target cells are at least 1.5 times wider than the source cells. For a DEM of difference the method sets how much false change an offset grid creates. Bilinear smooths, so the difference follows curvature: ridges appear to lose height and valleys to gain it. GDAL's lanczos can shift the result by a source cell when downsampling, so use it only between grids of equal resolution. Use near for categorical values. Other allowed values follow GDAL warp methods: near, bilinear, cubic, cubicspline, lanczos, average, mode, max, min, med, q1, q3, sum. |
| `reproject_to_metric_grid` | `bool` | `false` | `param` | - | Put a geographic first raster on a UTM grid, and warp the second to it. | Global DEMs are distributed in degrees, and coregistration needs metres. Each DEM is warped once, directly from its native grid. The spacing comes from the north-south pixel size (30 m for 1 arc-second DEMs). Writes &lt;write_fname&gt;_&lt;projection suffix&gt;. Ignored if the first raster is already projected. |
| `verbose_alignment_report` | `bool` | `true` | `param` | - | Print grid-compatibility and overlap diagnostics. | If true, reports whether rasters already share a grid and whether resampling was required. |

### Vertical reference

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `vertical_units` | `string` | `metres` | `param` | - | Vertical unit of the first raster: metres, feet or us_survey_feet. | No raster format records this, so it cannot be detected. Lidar is often in feet on a metre grid. Elevations are converted to metres before any output. |
| `second_raster_vertical_units` | `string` | `metres` | `param` | - | Vertical unit of the second raster: metres, feet or us_survey_feet. | See vertical_units. |
| `first_raster_vertical_datum` | `string` | `unspecified` | `param` | - | Geoid of the first raster: unspecified, egm96 or egm2008. | SRTM, NASADEM and AW3D30 use EGM96; Copernicus GLO-30 uses EGM2008. Set both datums to convert the second raster to the first raster's geoid. The difference varies across a DEM, so a constant vertical shift does not remove it. Needs PROJ's geoid grids: EGM96 is bundled, EGM2008 (77 MB) is not; a missing grid stops the run before any raster is read, naming the command that installs it (tools/fetch_proj_grids.py; see the note on datums and PROJ grids). |
| `second_raster_vertical_datum` | `string` | `unspecified` | `param` | - | Geoid of the second raster: unspecified, egm96 or egm2008. | See first_raster_vertical_datum. Give both or neither. |

### DEM coregistration

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_dem_coregistration` | `bool` | `false` | `switch` | - | Enable DEM coregistration before algebra outputs. | If true, the second raster is coregistered to the first raster after optional resampling and before all algebra operations. |
| `dem_coregistration_method` | `string` | `nuth_kaab` | `param` | `run_dem_coregistration` | DEM coregistration method. | Allowed values: nuth_kaab. Method reference: Nuth and Kaab (2011), https://doi.org/10.5194/tc-5-271-2011. This implementation was ported from xdem. |
| `dem_coregistration_stable_mask_fname` | `string` | - | `param` | `run_dem_coregistration` | Raster of stable terrain used to fit the shift (non-zero = stable). | Exclude landslides, lakes, ice, forest change and other real change, so that they do not bias the shift. Read from read_path if the path is relative, and warped to the first raster's grid with nearest neighbour. Also used for the reported difference statistics. |
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
