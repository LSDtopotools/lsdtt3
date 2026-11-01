# CLI reference: lsdtt-flat-surface-finder

**Pure morphometry: which ground is flat?** Identifies broadly flat surfaces —
floodplains, terrace treads, fans, piedmonts — from local relief and/or gradient
alone, and writes them as a binary mask and a labelled mask.

It does **no flow routing, no channel network, no height reference and no
floodplain/terrace classification.** Everything that needs a channel moved to
[`lsdtt-fluvial-surfaces`](cli_reference_lsdtt_fluvial_surfaces.md), which takes
this tool's `_combined_flat_mask.tif` as its `flat_mask_fname`. For the reasoning
behind the split and the algorithms, see
[Floodplains, terraces and valleys](../design_and_algorithms/fluvial_surfaces_and_valleys.md).

## Usage

```bash
lsdtt-flat-surface-finder <path> <parameter_file>
lsdtt-flat-surface-finder <parameter_file>      # path defaults to ./
lsdtt-flat-surface-finder -h                    # parameter help, then exit
lsdtt-flat-surface-finder -h <parameter_name>   # full guidance on one parameter
lsdtt-flat-surface-finder -h html               # write an HTML/CSV reference
lsdtt-flat-surface-finder --help-json           # machine-readable parameter list
lsdtt-flat-surface-finder -v                    # version
lsdtt-flat-surface-finder -c                    # citation
```

Parameter file lines are `key: value`; `#` starts a comment. Keys use
`snake_case`.

---

## Pipeline

| Step | What happens | Controlled by |
|---|---|---|
| 1 | Load DEM, check it is projected in metres | — |
| 2 | Local relief (max − min) in a square window | `relief_window_radius_metres` |
| 3 | Bound relief, fixed or chosen from the data | `max_relief_metres`, `relief_threshold_algorithm_choice` |
| 4 | Optional gradient bound from a fitted polynomial surface | `max_admission_gradient`, `admission_gradient_window_radius_metres` |
| 5 | Combine the two tests | `admission_algorithm_choice` |
| 6 | Morphological opening / closing | `flat_opening_radius_metres`, `flat_closing_radius_metres`, `open_then_close` |
| 7 | Optional directional closing of the non-flat mask | `directional_close_line_length_pixels` |
| 8 | (Riser / break-of-slope search moved to `lsdtt-feature-finder`, `feature_type: find_risers`) | — |
| 9 | Component filtering and labelling | `min_component_pixels` |
| 10 | Optional plane fitting per label | `run_plane_analysis` |
| 11 | Write masks, rasters and CSVs | the `write_*` switches |

---

## File I/O and run control

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | `./` | Directory from which input files are read. |
| `read_fname` | string | | Input filename, with or without extension (e.g. `my_dem.tif`). |
| `write_path` | string | same as `read_path` | Directory to which output files are written. |
| `write_fname` | string | stem of `read_fname` | Output filename prefix (no path, no extension). |
| `n_threads` | int | `0` | OpenMP thread count. `0` = all available. |
| `write_output_format` | string | `cog` | Raster output format: `cog`, `gtiff`, `same_as_input`. |

---

## Relief (always on)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `relief_window_radius_metres` | float | `0` | **Radius** (m) of the square window local relief (max − min) is measured over, converted to the nearest whole number of cells. `0` = 2 cells, a 5×5 window. The same name and meaning as in `lsdtt-fluvial-surfaces` and `lsdtt-mountain-front`. |
| `max_relief_metres` | float | `1.0` | Local relief (m) below which a pixel is admitted as flat. |
| `relief_threshold_algorithm_choice` | string | `fixed` | `fixed` uses `max_relief_metres`; `cdf_divergence` chooses the bound at the point of maximum divergence between the empirical CDF and the best-fitting parametric CDF; `qq_gaussian` chooses it from a quantile-quantile plot against a Gaussian (Clubb et al. 2017). |
| `relief_threshold_multiplier` | float | `1.0` | Multiplier applied to a bound chosen by `cdf_divergence`. |

`relief_threshold_algorithm_choice: cdf_divergence` writes `<write_fname>_relief_adaptive_cdf.csv` — the
whole curve, so you can see the break the run relied on rather than only the
number it picked.

---

## Admission

| Parameter | Type | Default | Description |
|---|---|---|---|
| `admission_algorithm_choice` | string | `relief` | Which test admits a pixel: `relief`, `gradient`, `relief_and_gradient`, `relief_or_gradient`. |
| `max_admission_gradient` | float | `0` | Gradient (m/m) below which a pixel is admitted. **Required** by every choice that names gradient. |
| `admission_gradient_window_radius_metres` | float | `0` | Window radius (m) of the polynomial fit the gradient is measured over. `0` = match the relief window. |
| `write_admission_gradient_raster` | bool | `false` | Write the gradient raster (m/m) the bound was applied to. |

**Relief and gradient are different tests and reject different ground.** Relief
charges a surface for its own tilt, so a smooth tread that slopes downvalley
fails it while a rough level surface passes; a gradient bound does the opposite.
On steep terrain the two admit nearly disjoint sets. Use `relief_or_gradient` when
the point is to find every candidate surface and let `lsdtt-fluvial-surfaces` sort
them.

The gradient comes from a polynomial surface fitted over
`admission_gradient_window_radius_metres`, not a 3×3 difference: a three-cell
operator on a 30 m DEM measures pixel noise. Gradients are m/m, not angles. There
is no default bound worth guessing — about 0.087 m/m is a valley floor,
0.27 m/m admits fan and pediment — so run the census at more than one
bound and report the sensitivity.

---

## Morphological post-processing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `flat_closing_radius_metres` | float | `0` | Closing radius (m) on the relief flat mask. `0` = skip. |
| `flat_opening_radius_metres` | float | `0` | Opening radius (m) on the relief flat mask. `0` = skip. |
| `open_then_close` | bool | `false` | Apply opening before closing instead of closing before opening. |

## Directional closing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_healed_flat_mask` | bool | `true` | Write the healed flat mask (after directional closing, if enabled). |
| `directional_close_line_length_pixels` | int | `0` | Line structuring-element length (pixels) for directional closing of the **non-flat** mask. `0` = disabled. |
| `directional_close_n_angles` | int | `18` | Number of line orientations sampled from [0°, 180°). |

## Component filtering

| Parameter | Type | Default | Description |
|---|---|---|---|
| `min_component_pixels` | int | `100` | Minimum 8-connected component size for inclusion in the labelled flat mask. |
| `cdf_n_bins` | int | `500` | Histogram bins for the empirical CDF output. |

Note that `min_component_pixels` filters components of the **whole mask** here.
`lsdtt-fluvial-surfaces` declares the same name but applies it **per band**.

---

## Slope-break skeleton — moved to `lsdtt-feature-finder`

The break-of-slope (VEC) skeleton and its endpoint bridging are now
`lsdtt-feature-finder` `feature_type: find_risers`, which writes open riser
segments with a step height. Its centreline can split surfaces in
`lsdtt-fluvial-surfaces` (`segmentation_algorithm_choice: whole_mask`,
`riser_mask_fname`).

Old keys in a parameter file here:

| Old key | Now |
|---|---|
| `write_slope_break_centreline_raster`, `write_slope_break_centreline_connected_raster`, `write_slope_break_bridge_mask_raster`, `write_slope_break_endpoints` (and their older spellings), `apply_slope_breaks_to_labels` | **Stop the run**, naming the replacement: they turned an output or a mask cut on. |
| `polyfit_lambda_metres`, `vec_adaptive_method`, `vec_adaptive_block_size_pixels`, `vec_adaptive_C`, `vec_adaptive_gaussian_sigma`, `slope_break_closing_radius_metres`, `slope_break_binary_min_component_pixels`, `slope_break_centreline_max_spur_pixels`, `slope_break_centreline_max_terminal_loop_pixels`, `slope_break_centreline_min_component_pixels`, `slope_break_connect_search_radius_metres`, `slope_break_connect_elevation_tolerance_metres`, `slope_break_connect_corridor_half_width_pixels`, `slope_break_connect_same_component` | Warn and are ignored; each message names its `riser_*` key in `lsdtt-feature-finder`. |

## Plane analysis

Fits a plane to each labelled surface and reports splits and merges. The fitting
and its diagnostics are current; the **re-break pass is deprecated**, because
banding the residual field in `lsdtt-fluvial-surfaces` achieves the same
separation without a second labelling pass over the whole raster.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `run_plane_analysis` | bool | `false` | Fit planes to each labelled flat surface and test for splits/merges. |
| `plane_min_pixels` | int | `10` | Minimum pixel count for a label to enter plane analysis. |
| `plane_merge_tolerance_rmse` | float | `2.0` | Merge threshold, in units of RMSE. |
| `plane_split_histogram_bins` | int | `25` | Histogram bins for bimodal split detection per label. |
| `plane_split_bimodal_threshold` | float | `0.4` | Valley/peak ratio threshold for bimodal split detection. |
| `write_plane_on_label_rmse` | bool | `false` | Write the dense plane-on-label RMSE matrix to CSV. |
| `write_plane_residual_raster` | bool | `true` | Write per-pixel DEM-minus-plane residuals for each label. |
| `write_residual_relief_raster` | bool | `true` | Write local relief computed on the plane-residual raster. |

### Residual-relief re-break — deprecated

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enable_residual_relief_rebreak` | bool | `false` | Run the residual-relief skeleton re-break pass. Gates every row below. |
| `residual_relief_threshold_multiplier` | float | `0.5` | Threshold for the residual-relief filter, as a multiplier of the flat-mask relief threshold. |
| `residual_relief_spur_max_pixels` | int | `10` | Maximum spur length (pixels) pruned from the residual-relief skeleton. |
| `merged_residual_relief_min_label_pixels` | int | `10` | Minimum pixel count for a label retained in the merged mask. |
| `write_residual_relief_skeleton_raster` | bool | `false` | Write the 1-pixel residual-relief skeleton used as break lines. |
| `write_residual_relief_labelled_mask_raster` | bool | `true` | Write the labelled flat mask after the residual-relief filter. |
| `write_residual_relief_merge_candidates` | bool | `false` | Write merge candidate pairs to CSV. |

---

## Validation input

| Parameter | Type | Default | Description |
|---|---|---|---|
| `test_surface_polygon_fname` | string | | Polygon vector file delineating a known flat surface, for comparison against the detected mask. |

## Output toggles

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_combined_flat_mask_raster` | bool | `true` | Write the combined binary flat mask. **This is the file `lsdtt-fluvial-surfaces` consumes.** |
| `write_metric_flat_mask_rasters` | bool | `false` | Write binary flat masks for each individual enabled metric. |
| `write_texture_rasters` | bool | `false` | Write per-metric texture rasters (IQR residuals, IQR texture). |

---

## Outputs

| File | Written when |
|---|---|
| `<write_fname>_combined_flat_mask.tif` | `write_combined_flat_mask_raster` (default on) — **the input to `lsdtt-fluvial-surfaces`** |
| `<write_fname>_labelled_flat_mask.tif` | always |
| `<write_fname>_relief.tif`, `_relief_flat_mask.tif` | `write_metric_flat_mask_rasters` |
| `<write_fname>_admission_gradient.tif` | `write_admission_gradient_raster` |
| `<write_fname>_relief_cdf.csv`, `_dist_fits.csv` | always |
| `<write_fname>_relief_adaptive_cdf.csv` | `relief_threshold_algorithm_choice: cdf_divergence` |
| `<write_fname>_flat_vertical_excess_curvature.tif`, `_slope_break_*` | the slope-break switches |
| `<write_fname>_plane_fits.csv`, `_plane_residual.tif`, `_residual_relief*.tif` | `run_plane_analysis` |
| `<write_fname>_README.csv` / `.html` | every successful run — the resolved parameters for that run |

---

## Renamed and retired parameters

The parameters below were renamed or moved when the tool was split. The run
prints a message naming the replacement; keys marked **removed** stop the run
rather than being ignored, so a stale parameter file cannot silently do nothing.

| Old key | Now |
|---|---|
| `find_floodplain` | **removed** — floodplain and terrace classification is `lsdtt-fluvial-surfaces` |
| `enable_leak_split`, `leak_split_min_side_fraction` | **removed** — band-then-label makes a HAND leak structurally impossible |
| `enable_hand_banding` | **removed** — the ladder is the method, not an option |
| `band_split_min_label_spread_metres` | **removed** — there is no parent label to split |
| `apply_plane_merges` | **removed** — it was never read |
| `write_combined_flat_mask` | `write_combined_flat_mask_raster` |
| `write_metric_flat_masks` | `write_metric_flat_mask_rasters` |
| `write_residual_relief_labelled_mask` | `write_residual_relief_labelled_mask_raster` |
| `write_residual_relief_skeleton` | `write_residual_relief_skeleton_raster` |
| `test_surface_vector_file` | `test_surface_polygon_fname` |
| `relief_window_side_pixels`, `relief_window_pixels` | `relief_window_radius_metres` — the full side in cells is converted to a radius in metres with the cell size |
| `relief_threshold_metres` | `max_relief_metres` |
| `adaptive_relief_threshold` | `relief_threshold_algorithm_choice` — `true` is `cdf_divergence`, `false` is `fixed` |
| `admission_mode` | `admission_algorithm_choice`; its values `slope`, `relief_and_slope`, `relief_or_slope` are now `gradient`, `relief_and_gradient`, `relief_or_gradient` |
| `admission_slope_cap_degrees` | `max_admission_gradient` — converted with gradient = tan(angle) |
| `admission_slope_window_radius_metres` | `admission_gradient_window_radius_metres` |
| `write_admission_slope_raster` | `write_admission_gradient_raster`; the file is now `_admission_gradient.tif`, in m/m |
| `floodplain_max_height_metres`, `terrace_min_height_metres`, `channel_adjacency_radius_metres`, `terrace_max_hand_spread_metres`, `terrace_corridor_per_pixel`, `band_reach_window_metres`, `terrace_max_distance_to_channel_metres` | moved to `lsdtt-fluvial-surfaces` (where the first two are now `max_floodplain_hand_metres` and `min_terrace_hand_metres`) |
| every `source_points_*`, `outlet_points_*`, `write_channel_*`, `break_flats_with_channels`, `fill_before_routing`, `connect_flat_surfaces_to_channels` | moved to `lsdtt-fluvial-surfaces` — this tool no longer does flow routing |

---

## See also

- [Floodplains, terraces and valleys](../design_and_algorithms/fluvial_surfaces_and_valleys.md) — the
  algorithms, and where they come from
- [CLI reference: lsdtt-fluvial-surfaces](cli_reference_lsdtt_fluvial_surfaces.md)
- [CLI reference: lsdtt-valley-centreline](cli_reference_lsdtt_valley_centreline.md)

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
| `max_relief_metres` | `float` | `1` | `param` | - | Local relief (m) below which a pixel is admitted as flat. | Read a value off the relief CDF the tool writes. Where the tool has relief_threshold_algorithm_choice, this is the value used by `fixed`. |
| `relief_threshold_algorithm_choice` | `string` | `fixed` | `param` | - | How the relief bound is chosen: fixed, cdf_divergence or qq_gaussian. | fixed uses max_relief_metres. cdf_divergence computes the relief CDF, fits the best of a family of parametric distributions (normal, exponential, Weibull, Pareto) and takes the relief at which the empirical CDF most exceeds the fit: flat pixels inflate the empirical distribution at low relief. qq_gaussian is the Clubb et al. (2017) rule: the lowest relief at which the quantiles rejoin a Gaussian drawn through two percentiles of the QQ plot (relief_qq_*) and stay on it. The table and the chosen value are written out. Replaces the boolean adaptive_relief_threshold. |
| `relief_threshold_multiplier` | `float` | `1.0` | `param` | `relief_threshold_algorithm_choice` | Multiplier applied to an automatically chosen relief bound. | Allows leeway around the auto-detected value. Values &gt; 1 loosen the criterion (more flat pixels); values &lt; 1 tighten it. Not used when relief_threshold_algorithm_choice is fixed. |
| `admission_algorithm_choice` | `string` | `relief` | `param` | - | Which test admits a pixel as flat: relief, gradient, relief_and_gradient or relief_or_gradient. | Local relief is not a gradient test and the two reject different ground. Relief charges a surface for its own tilt, so a smooth tread that slopes downvalley fails it while a rough level surface passes; a gradient bound does the opposite. On a Himalayan gorge the two select nearly disjoint sets. `relief` is the historical behaviour. `gradient` uses max_admission_gradient alone. `relief_and_gradient` requires both, which is the strictest reading of flat. `relief_or_gradient` admits either, which is the one to use when the point is to find every candidate surface and let the ladder sort them. Replaces admission_mode (whose `slope` values are now `gradient`). |
| `max_admission_gradient` | `float` | `0` | `param` | `admission_algorithm_choice` | Gradient (m/m) below which a pixel is admitted as flat. | A gradient, not an angle: 0.087 is 5 degrees, 0.176 is 10, 0.268 is 15 (gradient = tan(angle)). Measured from a fitted polynomial surface rather than a 3x3 difference. Required by every admission_algorithm_choice that names gradient; there is no default worth guessing - 0.087 (5 degrees) is a valley floor, 0.27 (15 degrees) admits fan and pediment - so run the census at more than one bound and report the sensitivity rather than defending one number. Replaces admission_slope_cap_degrees; an old value is converted with gradient = tan(angle). |
| `admission_gradient_window_radius_metres` | `float` | `0` | `param` | `admission_algorithm_choice` | Window radius (m) of the polynomial fit the admission gradient is measured over. 0 = match the relief window. | Defaulting to the relief window keeps the two admission tests commensurate: they then disagree about the ground, not about the scale they are looking at. |
| `write_admission_gradient_raster` | `bool` | `false` | `output` | - | Write the gradient raster (m/m) the admission bound was applied to. | Output file: &lt;stem&gt;_admission_gradient.tif. Write it the first time you set a bound on a new DEM. The bound is only meaningful against the gradient distribution it is cutting, and that distribution is not knowable in advance. |

### Threshold selection

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `relief_qq_lower_percentile` | `float` | `25` | `param` | `relief_threshold_algorithm_choice` | Lower percentile of the QQ cloud the Gaussian reference line is drawn through (qq_gaussian). | 25 and 75 (the IQR) are the values Clubb et al. (2017) used. |
| `relief_qq_upper_percentile` | `float` | `75` | `param` | `relief_threshold_algorithm_choice` | Upper percentile of the QQ cloud the Gaussian reference line is drawn through (qq_gaussian). | Must be above the lower percentile. |
| `relief_qq_threshold_fraction` | `float` | `0.01` | `param` | `relief_threshold_algorithm_choice` | Fraction of the data range within which a quantile counts as on the Gaussian line (qq_gaussian). | 0.01 (1%) is the paper's value. The threshold is the lowest quantile that is on the line and stays on it for the next 200 quantiles of 10000. Larger values admit more ground. |

### Morphological post-processing of relief mask

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `flat_closing_radius_metres` | `float` | `0` | `param` | - | Radius (m) for morphological closing of the relief flat mask. | Closing (dilation then erosion) seals small holes and narrow gaps inside flat blobs without expanding their outer boundary. 0 = disabled. |
| `flat_opening_radius_metres` | `float` | `0` | `param` | - | Radius (m) for morphological opening of the relief flat mask. | Opening (erosion then dilation) removes thin connections between blobs that should be separate flat surfaces, while restoring the size of the surviving blobs. 0 or negative = disabled. |
| `open_then_close` | `bool` | `false` | `param` | - | Apply opening before closing instead of closing before opening. | false (default): closing first, then opening. true: opening first, then closing. Has no effect when either radius is 0. |

### Component filtering

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `min_component_pixels` | `int` | `100` | `param` | - | Minimum 8-connected component size for inclusion in the labelled flat mask. | Smaller components are excluded. Set to 1 to include all. |
| `cdf_n_bins` | `int` | `500` | `param` | - | Number of histogram bins for the empirical CDF output. | - |

### Directional closing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `directional_close_line_length_pixels` | `int` | `0` | `param` | `write_healed_flat_mask` | Line-SE length (pixels) for directional closing of the non-flat mask. 0 = disabled. | Applied after mask combination, before component labelling. |
| `directional_close_n_angles` | `int` | `18` | `param` | `write_healed_flat_mask` | Number of line orientations sampled from [0 deg, 180 deg) for directional closing. | 18 -&gt; 10 deg spacing. |
| `write_healed_flat_mask` | `bool` | `true` | `switch` | - | Write the healed flat mask (after directional closing, if enabled). | Only written when directional_close_line_length_pixels &gt;= 3. |

### Plane analysis

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_plane_analysis` | `bool` | `false` | `switch` | - | Fit planes to each labelled flat surface and test for splits/merges. | Runs after labelling. Writes _plane_fits.csv and _merge_candidates.csv. |
| `plane_min_pixels` | `int` | `10` | `param` | - | Minimum pixel count for a label to be included in plane analysis. | - |
| `plane_merge_tolerance_rmse` | `float` | `2.0` | `param` | - | Merge threshold in units of RMSE. | - |
| `plane_split_histogram_bins` | `int` | `25` | `param` | - | Number of histogram bins for bimodal split detection per label. | - |
| `plane_split_bimodal_threshold` | `float` | `0.4` | `param` | - | Valley/peak ratio threshold for bimodal split detection. | - |
| `write_plane_on_label_rmse` | `bool` | `false` | `switch` | - | Write the dense plane-on-label RMSE matrix to CSV. | Output suffix: _plane_on_label_rmse.csv Requires run_plane_analysis = true. |
| `enable_residual_relief_rebreak` | `bool` | `false` | `switch` | `run_plane_analysis` | Run the residual-relief skeleton re-break pass (deprecated). | Default false. Gates the second labelling pass that skeletonises the relief of the plane residual and re-breaks the labels with it. Deprecated: lsdtt-fluvial-surfaces bands the residual field directly, which achieves the same separation without a second pass over the whole raster. Scheduled for deletion after one release. |
| `write_residual_relief_merge_candidates` | `bool` | `false` | `switch` | - | Write merge candidate pairs from the residual-relief plane analysis to CSV. | Output suffix: _residual_relief_merge_candidates.csv Requires run_plane_analysis = true. |
| `write_plane_residual_raster` | `bool` | `true` | `switch` | - | Write per-pixel residuals between DEM elevation and fitted plane for each label. | Output suffix: _plane_residual.tif Requires run_plane_analysis = true. |
| `write_residual_relief_raster` | `bool` | `true` | `switch` | - | Write local relief computed on the plane-residual raster. | Uses the same window as the original flat mask. Output suffix: _residual_relief.tif Requires run_plane_analysis = true. |
| `residual_relief_threshold_multiplier` | `float` | `0.5` | `param` | - | Threshold for the residual relief filter, as a multiplier of the flat-mask relief threshold. | Pixels where residual_relief &gt; multiplier * relief_threshold are removed from the flat mask before the residual-relief labelling step. |
| `write_residual_relief_labelled_mask_raster` | `bool` | `true` | `switch` | - | Write the labelled flat mask after applying the residual-relief filter. | Output suffix: _residual_relief_labelled_flat_mask.tif Requires run_plane_analysis = true. |
| `residual_relief_spur_max_pixels` | `int` | `10` | `param` | - | Maximum spur length (pixels) pruned from the residual-relief skeleton. | 0 = skip spur pruning. |
| `write_residual_relief_skeleton_raster` | `bool` | `false` | `switch` | - | Write the 1-pixel residual-relief skeleton used as break lines. | Output suffix: _residual_relief_skeleton.tif Requires run_plane_analysis = true. |
| `merged_residual_relief_min_label_pixels` | `int` | `10` | `param` | - | Minimum pixel count for a label to be retained in the merged residual-relief mask. | Labels smaller than this are set to 0 before writing _merged_residual_relief_labelled_flat_mask.tif. |

### Test surface vector

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `test_surface_polygon_fname` | `string` | - | `param` | - | Path to a polygon vector file delineating a known flat surface. | When set, the maximum value of each enabled metric within the polygon pixels replaces the percentile threshold for that metric. |

### Output toggles

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_texture_rasters` | `bool` | `false` | `switch` | - | Write per-metric texture rasters (IQR residuals, IQR texture). | Useful for threshold selection and visual inspection. |
| `write_metric_flat_mask_rasters` | `bool` | `false` | `switch` | - | Write binary flat masks for each individual enabled metric. | - |
| `write_combined_flat_mask_raster` | `bool` | `true` | `switch` | - | Write the combined binary flat mask. | 1 = flat, 0 = not flat. |

### Orientations

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `bearings_relative_to_true_north` | `bool` | `false` | `param` | - | Write bearings, aspects and orientations clockwise from true north. | Off (default): they are clockwise from grid north (grid up), so in an oblique frame with the feature along +x, 90 degrees means along the feature. On: each value is converted with the grid convergence at its own location (each cell of a raster, each row of a table), which also removes the up to about 3 degree difference between UTM grid north and true north. Angles between two directions (junction angles, deflections, contrasts) are unchanged, and hillshade illumination always stays relative to grid up. Needs a raster with a CRS. |

<!-- END GENERATED CLI PARAMETERS -->
