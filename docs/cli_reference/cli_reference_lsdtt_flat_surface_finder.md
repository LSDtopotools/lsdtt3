# CLI reference: lsdtt-flat-surface-finder

**Pure morphometry: which ground is flat?** Identifies broadly flat surfaces —
floodplains, terrace treads, fans, piedmonts — from local relief and/or gradient
alone, and writes them as a binary mask and a labelled mask.

It does **no flow routing, no channel network, no height reference and no
floodplain/terrace classification.** Everything that needs a channel moved to
[`lsdtt-fluvial-surfaces`](cli_reference_lsdtt_fluvial_surfaces.md), which takes
this tool's `_combined_flat_mask.tif` as its `flat_mask_fname`. For the reasoning
behind the split and the algorithms, see
[Floodplains, terraces and valleys](../fluvial_surfaces_and_valleys.md).

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
| 2 | Local relief (max − min) in a square window | `relief_window_side_pixels` |
| 3 | Threshold relief, fixed or adaptive | `relief_threshold_metres`, `adaptive_relief_threshold` |
| 4 | Optional gradient cap from a fitted polynomial surface | `admission_slope_cap_degrees`, `admission_slope_window_metres` |
| 5 | Combine the two tests | `admission_mode` |
| 6 | Morphological opening / closing | `flat_opening_radius_metres`, `flat_closing_radius_metres`, `open_then_close` |
| 7 | Optional directional closing of the non-flat mask | `directional_close_line_length_pixels` |
| 8 | Optional VEC break-of-slope skeleton *(requires OpenCV)* | `vec_adaptive_*`, `slope_break_*` |
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
| `relief_window_side_pixels` | int | `5` | Full window **side length** in pixels for local relief (max − min). Not a radius: must be odd and ≥ 3, so `5` is a 5×5 neighbourhood. `lsdtt-fluvial-surfaces` spells its equivalent `fallback_relief_window_radius_pixels`, which is a **radius**, so the same number means a different window there. The old name `relief_window_pixels` is accepted as a deprecated alias. |
| `relief_threshold_metres` | float | `1.0` | Local relief (m) below which a pixel is classed flat. |
| `adaptive_relief_threshold` | bool | `false` | Choose the relief threshold automatically, at the point of maximum divergence between the empirical CDF and the best-fitting parametric CDF. Ignores `relief_threshold_metres`. |
| `relief_threshold_multiplier` | float | `1.0` | Multiplier applied to the adaptively chosen threshold. |

`adaptive_relief_threshold` writes `<write_fname>_relief_adaptive_cdf.csv` — the
whole curve, so you can see the break the run relied on rather than only the
number it picked.

---

## Admission

| Parameter | Type | Default | Description |
|---|---|---|---|
| `admission_mode` | string | `relief` | Which test admits a pixel: `relief`, `slope`, `relief_and_slope`, `relief_or_slope`. |
| `admission_slope_cap_degrees` | float | `0` | Gradient (degrees) below which a pixel is admitted. **Required** by every mode that names slope. |
| `admission_slope_window_metres` | float | `0` | Window diameter (m) of the polynomial fit the gradient is measured over. `0` = match the relief window. |
| `write_admission_slope_raster` | bool | `false` | Write the gradient raster the cap was applied to. |

**Relief and slope are different tests and reject different ground.** Relief
charges a surface for its own tilt, so a smooth tread that slopes downvalley
fails it while a rough level surface passes; a gradient cap does the opposite. On
steep terrain the two admit nearly disjoint sets. Use `relief_or_slope` when the
point is to find every candidate surface and let `lsdtt-fluvial-surfaces` sort
them.

The gradient comes from a polynomial surface fitted over
`admission_slope_window_metres`, not a 3×3 difference: a three-cell operator on a
30 m DEM measures pixel noise. There is no default cap worth guessing — about
5 degrees is a valley floor, 15 degrees admits fan and pediment — so run the
census at more than one cap and report the sensitivity.

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

## Slope-break skeleton (VEC) — diagnostic, requires OpenCV

A break-of-slope skeleton derived by adaptive thresholding of vertical excess
curvature. It is a morphometric product worth having in its own right; it is
**not** part of the flat-surface admission and does not break the mask unless you
ask it to.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `apply_slope_breaks_to_labels` | bool | `false` | Break the flat mask with the skeleton before labelling. |
| `polyfit_lambda_metres` | float | `100` | Window half-wavelength (m) for the direct polyfit curvature computation. |
| `vec_adaptive_method` | string | `gaussian` | Adaptive thresholding method applied to vertical excess curvature. |
| `vec_adaptive_block_size_pixels` | int | `51` | Neighbourhood block size (pixels) for adaptive thresholding. |
| `vec_adaptive_C` | float | `0.0` | Constant subtracted from the local neighbourhood mean before thresholding. |
| `vec_adaptive_gaussian_sigma` | float | `0.0` | Sigma for an optional Gaussian pre-blur. |
| `slope_break_closing_radius_metres` | float | `0` | Closing radius (m) on the adaptive binary mask before skeletonisation. |
| `slope_break_binary_min_component_pixels` | int | `1` | Minimum component size retained in the adaptive binary mask. |
| `slope_break_centreline_max_spur_pixels` | int | `10` | Maximum endpoint-rooted spur length pruned from the skeleton. |
| `slope_break_centreline_max_terminal_loop_pixels` | int | `10` | Maximum endpoint-rooted terminal loop length pruned. |
| `slope_break_centreline_min_component_pixels` | int | `1` | Minimum component size retained in the skeleton. |
| `write_slope_break_centreline` | bool | `true` | Write the 1-pixel skeleton. |

### Endpoint bridging — deprecated

Bridging existed to re-close gaps in the riser network so that two treads at
different heights could not merge into one label. `lsdtt-fluvial-surfaces` labels
**within** height bands, which makes that merge impossible rather than
detectable, so the job is gone. Bridging is also the only route by which the
skeleton ever broke the flat mask, and it is off by default
(`slope_break_connect_search_radius_metres: 0`). **These five are scheduled for
deletion one release after the release that first ships them as deprecated.**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `slope_break_connect_search_radius_metres` | float | `0` | Search radius (m) for elevation-based endpoint connection. `0` = off. |
| `slope_break_connect_elevation_tolerance_metres` | float | `2` | Maximum elevation difference (m) between two bridged endpoints. |
| `slope_break_connect_corridor_half_width_pixels` | float | `12` | Half-width (pixels) of the A\* corridor for contour-routed bridges. |
| `slope_break_connect_same_component` | bool | `false` | Allow bridging between endpoints on the same skeleton component. |
| `write_slope_break_endpoints` | bool | `false` | Write skeleton endpoints and their elevations to FlatGeobuf. |
| `write_slope_break_bridge_mask` | bool | `true` | Write the raster of newly inserted bridge pixels. |
| `write_slope_break_centreline_connected` | bool | `true` | Write the connected skeleton after accepted bridge routing. |

---

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
| `<write_fname>_admission_slope_degrees.tif` | `write_admission_slope_raster` |
| `<write_fname>_relief_cdf.csv`, `_dist_fits.csv` | always |
| `<write_fname>_relief_adaptive_cdf.csv` | `adaptive_relief_threshold` |
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
| `floodplain_max_height_metres`, `terrace_min_height_metres`, `channel_adjacency_radius_metres`, `terrace_max_hand_spread_metres`, `terrace_corridor_per_pixel`, `band_reach_window_metres`, `terrace_max_distance_to_channel_metres` | moved to `lsdtt-fluvial-surfaces` |
| every `source_points_*`, `outlet_points_*`, `write_channel_*`, `break_flats_with_channels`, `fill_before_routing`, `connect_flat_surfaces_to_channels` | moved to `lsdtt-fluvial-surfaces` — this tool no longer does flow routing |

---

## See also

- [Floodplains, terraces and valleys](../fluvial_surfaces_and_valleys.md) — the
  algorithms, and where they come from
- [CLI reference: lsdtt-fluvial-surfaces](cli_reference_lsdtt_fluvial_surfaces.md)
- [CLI reference: lsdtt-valley-centreline](cli_reference_lsdtt_valley_centreline.md)
