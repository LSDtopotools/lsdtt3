# CLI reference: lsdtt-flat-surface-finder

Identifies broadly flat surfaces (floodplains, alluvial fans, terraces) in a DEM using local relief and optional slope-break detection. Outputs labelled flat masks, slope-break centrelines, and statistical summaries.

## Usage

```
lsdtt-flat-surface-finder  <path> <parameter_file>
lsdtt-flat-surface-finder  <parameter_file>        # path defaults to ./
lsdtt-flat-surface-finder  -h                      # print help and exit
lsdtt-flat-surface-finder  -v                      # print version and exit
lsdtt-flat-surface-finder  -c                      # print citation and exit
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
| 1. Relief computation | Compute local relief at multiple radii |
| 2. Relief thresholding | Classify low-relief pixels as candidates |
| 3. Morphological filtering | Opening/closing to clean mask |
| 4. Component filtering | Remove small isolated patches |
| 5. Slope-break detection *(optional)* | Identify edges using curvature and skeleton |
| 6. Channel breaking *(optional)* | Separate flats at channel boundaries |
| 7. Labelling | Assign unique IDs to each flat surface |
| 8. Plane analysis *(optional)* | Fit planes and split/merge surfaces |
| 9. Outputs | Write masks, centrelines, and statistics |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input DEM filename. Any GDAL-supported format. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` = all available hardware threads. |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `raster_is_filled` | bool | `false` | Set `true` if DEM is already pit-filled. |
| `fill_min_slope` | float | `0.0001` | Minimum slope (m/m) imposed during pit filling. |

---

## Relief thresholding

| Parameter | Type | Default | Description |
|---|---|---|---|
| `relief_radius_metres` | float | `100.0` | Radius (metres) for local relief calculation. |
| `relief_threshold_metres` | float | `10.0` | Relief threshold (metres). Pixels with relief < threshold are flat candidates. |
| `adaptive_relief_threshold` | bool | `false` | Adaptively select relief threshold using CDF divergence. If `true`, ignores `relief_threshold_metres`. |
| `relief_threshold_multiplier` | float | `1.0` | Scale factor applied to adaptively-selected relief threshold. |

---

## Morphological processing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `flat_opening_radius_metres` | float | `0.0` | Morphological opening radius (metres). `0` = skip. Removes thin connections. |
| `flat_closing_radius_metres` | float | `0.0` | Morphological closing radius (metres). `0` = skip. Seals small holes. |
| `open_then_close` | bool | `true` | If `true`, apply opening then closing. If `false`, apply closing then opening. |
| `minimum_component_pixels` | int | `100` | Minimum connected-component size (pixels). Smaller components removed. |

---

## Slope-break detection

Used only when OpenCV is available. Requires `LSDTT_HAS_OPENCV=1`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `detect_slope_breaks` | bool | `false` | Enable slope-break detection using curvature and skeletonisation *(requires OpenCV)*. |
| `curvature_threshold` | float | `0.01` | Vertical excess curvature threshold for slope-break detection. |
| `skeleton_prune_length_pixels` | int | `5` | Prune isolated skeleton spurs shorter than this (pixels). |
| `bridge_endpoints_by_elevation` | bool | `false` | Connect isolated slope-break endpoints if elevation difference is small. |
| `elevation_bridge_threshold_metres` | float | `2.0` | Maximum elevation difference (metres) for bridging isolated endpoints. |

---

## Channel breaking

| Parameter | Type | Default | Description |
|---|---|---|---|
| `break_flats_with_channels` | bool | `false` | Use D8 drainage-area threshold to separate flat surfaces at channels. |
| `channel_threshold_area_m2` | float | `1000000.0` | Minimum drainage area (m²) to define a channel. |
| `inject_drainage_area` | bool | `false` | Increase drainage area in flat areas to favour channel detection. |
| `drainage_area_injection_multiplier` | float | `2.0` | Factor to multiply drainage area in flat pixels. |

---

## Plane analysis

| Parameter | Type | Default | Description |
|---|---|---|---|
| `run_plane_analysis` | bool | `false` | Fit planes to flat surfaces and refine labelling. |
| `plane_residual_threshold_metres` | float | `1.0` | Maximum residual (metres) from plane fit to retain pixel in surface. |

---

## Output: Labelled flat mask

Always written.

| Parameter | Type | Default | Description |
|---|---|---|---|
| _(always written)_ | — | — | Write labelled flat-surface mask with unique ID per surface. |

Output file: `<write_fname>_labelled_flat_mask.tif`

---

## Output: Additional rasters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_combined_flat_mask` | bool | `false` | Write binary flat mask (union of all flat surfaces). |
| `write_relief_raster` | bool | `false` | Write local relief raster. |
| `write_slope_break_raster` | bool | `false` | Write binary slope-break mask. |
| `write_slope_break_centreline` | bool | `false` | Write slope-break centreline (skeleton). |

---

## Output: Statistics CSV

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_relief_statistics_csv` | bool | `true` | Write relief CDF and distribution-fit statistics. |
| `write_flat_statistics_csv` | bool | `false` | Write per-flat-surface area and elevation statistics. |

---

## Output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `cog` | Raster output format: `cog`, `gtiff`, `same_as_input`. |
