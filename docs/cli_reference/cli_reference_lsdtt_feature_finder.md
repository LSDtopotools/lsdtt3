# CLI reference: lsdtt-feature-finder

Identifies geomorphic features on DEMs including terraces, piedmonts, hillslope morphology transitions, and other terrain classifications. Uses surface morphology analysis, relief thresholding, and optional vector data integration.

## Usage

```
lsdtt-feature-finder  <path> <parameter_file>
lsdtt-feature-finder  <parameter_file>        # path defaults to ./
lsdtt-feature-finder  -h                      # print help and exit
lsdtt-feature-finder  -v                      # print version and exit
lsdtt-feature-finder  -c                      # print citation and exit
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
| 1. DEM preparation | Load DEM, verify projection, optionally fill pits |
| 2. Surface analysis | Compute surface morphology metrics (curvature, relief, texture) |
| 3. Feature detection | Classify pixels based on surface properties and thresholds |
| 4. Post-processing | Morphological operations, smoothing, vectorisation |
| 5. Outputs | Write feature masks and optional vector layers |

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
| `boundary_conditions` | string | `oooo` | Four-character boundary condition string (N, E, S, W). |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `raster_is_filled` | bool | `false` | Set `true` if DEM is already pit-filled. |
| `fill_min_slope` | float | `0.0001` | Minimum slope (m/m) imposed during pit filling. |

---

## Feature detection

| Parameter | Type | Default | Description |
|---|---|---|---|
| `feature_type` | string | `terrace` | Feature to detect. Options include `terrace`, `piedmont`, `hillslope_break`, and `flow_road_discrimination`. |
| `relief_threshold_metres` | float | `5.0` | Local relief threshold (metres) for feature detection. |
| `slope_threshold_degrees` | float | `5.0` | Slope threshold (degrees) for classification. |
| `smoothing_scale_metres` | float | `50.0` | Smoothing scale (metres) for surface fitting prior to analysis. |

For `feature_type: flow_road_discrimination`, `threshold_contributing_pixels`
(int, default `1000`) defines the D8 channel network used to distinguish flow
pathways from roads. `big_channel_threshold_contributing_pixels` (int, default
`5000`) is the larger contributing-pixel threshold used by the valley-channel
point export. The former `flow_big_channel_threshold_pixels` spelling remains
an unchanged-value alias.

---

## Post-processing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `minimum_component_size_pixels` | int | `100` | Minimum connected component size (pixels). Smaller components removed. |
| `apply_morphological_opening` | bool | `true` | Apply morphological opening to remove thin features. |
| `morphological_radius_pixels` | int | `3` | Radius (pixels) for morphological operations. |

---

## Output: Feature mask raster

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_feature_mask` | bool | `true` | Write binary feature mask raster. |

Output file: `<write_fname>_<feature_type>_mask.tif`

---

## Output: Feature polygons

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_feature_polygons` | bool | `false` | Vectorise feature mask and write as polygon layer. |
| `feature_polygon_format` | string | `flatgeobuf` | Output format: `flatgeobuf` / `fgb`, `gpkg`. |

Output file: `<write_fname>_<feature_type>_polygons.<ext>`

---

## Output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `cog` | Raster output format: `cog`, `gtiff`, `same_as_input`, `png`. `png` writes an 8-bit scaled image with a worldfile and PAM sidecar for quick visualisation only, not for analysis reuse. Output files take the extension of the chosen format. |

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

### General

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `feature_type` | `string` | `channel_candidates` | `param` | - | Feature mode to run. | Supported values: channel_candidates, find_drains, relief_ratio, relative_path_impact, local_median_relative_elevation, find_valleys, find_roads, flow_road_discrimination, path_opening. |
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all hardware threads. Applies to wavelet fitting and polyfit metric computation. |
| `write_output_format` | `string` | `cog` | `io` | - | Output raster format. | Allowed values: cog, gtiff, same_as_input, png. cog: Cloud Optimised GeoTIFF (default). gtiff: plain compressed GeoTIFF. same_as_input: match the input raster format. png: 8-bit scaled image with a worldfile and PAM sidecar for georeferencing - for quick visualisation only, not for analysis reuse. Output files take the extension of the chosen format. |
| `channel_lambda_min_metres` | `float` | `3.0` | `param` | - | Minimum wavelength for multiscale channel-candidate curvature analysis, in metres. | Use a value near the narrowest expected natural channel width or bottom width. |
| `channel_lambda_max_metres` | `float` | `18.0` | `param` | - | Maximum wavelength for multiscale channel-candidate curvature analysis, in metres. | Use a broader value that captures larger natural channel forms. |
| `channel_n_scales` | `int` | `5` | `param` | - | Number of scales used in the multiscale channel-candidate detector. | Scales are log-spaced between channel_lambda_min_metres and channel_lambda_max_metres. |
| `polyfit_curvature_metric` | `string` | `plan_curvature` | `param` | - | Polyfit curvature metric used for channel-candidate detection. | Suggested starting values: plan_curvature or tangential_curvature. |
| `polyfit_curvature_sign` | `string` | `positive` | `param` | - | Whether channel-like pixels correspond to positive or negative curvature. | Allowed values: positive (default), negative. Positive curvature = concave = convergent (channel); negative curvature = convex = divergent (ridge). |
| `polyfit_method` | `string` | `direct` | `param` | - | Polyfit algorithm: direct (rolling-moment, O(n), default) or wavelet (FFT-based). | direct is faster and better at resolving fine-scale features and is the recommended choice. |
| `polyfit_z_factor` | `float` | `1.0` | `param` | - | Vertical scaling applied before interpreting polyfit derivatives. | Usually leave at 1.0 unless experimenting. |
| `polyfit_flat_slope_tolerance` | `float` | `1e-12` | `param` | - | Flat-slope tolerance passed to the polyfit metric engine. | Very small default. |
| `channel_depth_window_radius_metres` | `float` | `12.0` | `param` | - | Radius of moving-window mean used to derive local depth support, in metres. | Depth signal is positive_part(local_mean - dem). |
| `depth_window_shape` | `string` | `circular` | `param` | - | Shape of moving window used for local depth support. | Allowed values: square, circular. |
| `channel_weight_persistent` | `float` | `0.7` | `param` | - | Weight for multiscale curvature persistence in the final score. | Usually the dominant term. |
| `channel_weight_depth` | `float` | `0.3` | `param` | - | Weight for local depth support in the final score. | Helps favour locally low pixels. |
| `channel_weight_slope` | `float` | `0.0` | `param` | - | Optional slope weight in the final score. | Recommended to keep small or zero at first. |
| `channel_score_threshold` | `float` | `0.5` | `param` | - | Threshold applied to final channel-candidate score to make a binary mask. | Lower is more inclusive. |
| `write_channel_candidate_score_raster` | `bool` | `true` | `switch` | - | Write final channel-candidate score raster. | Continuous score in [0,1]. |
| `write_channel_candidate_mask_raster` | `bool` | `true` | `switch` | - | Write thresholded binary channel-candidate mask. | 1 indicates plausible channel candidate. |
| `write_channel_best_scale_raster` | `bool` | `true` | `switch` | - | Write best persistent scale raster for channel candidates. | Useful later for cross-road matching. |
| `write_channel_candidate_diagnostic_rasters` | `bool` | `true` | `switch` | - | Write intermediate diagnostic rasters for the channel-candidate detector. | Includes persistent response, depth signal, slope signal, and Hessian elongation. |
| `write_channel_scale_diagnostic_rasters` | `bool` | `false` | `switch` | - | Write per-scale diagnostic rasters (verbose). | Includes raw curvature, normalised channelness signal, and pairwise persistence for every scale. Disabled by default to reduce output volume. |

### Hessian elongation filter

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `use_hessian_filter` | `bool` | `true` | `switch` | - | Apply a Hessian eigenanalysis filter to retain only line-like features. | Pixels where the dominant curvature eigenvalue does not substantially exceed the cross-track eigenvalue are removed from the mask. |
| `hessian_elongation_threshold` | `float` | `0.6` | `param` | `use_hessian_filter` | Minimum Hessian elongation to retain in the channel mask. | elongation = \|lambda_major\| / (\|lambda_major\| + \|lambda_minor\|). 0.5 is vacuous (always true); 0.6 removes near-isotropic pits/peaks; higher values are more selective. |

### Relief-ratio detector

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_relief_ratio_raster` | `bool` | `true` | `switch` | - | Write the relief-ratio raster. | The output is local_relief(short_radius) / local_relief(long_radius), using circular moving windows. |
| `relief_ratio_short_radius_metres` | `float` | `1.8` | `param` | `write_relief_ratio_raster` | Short moving-window radius for local relief, in metres. | Use a radius matched to the narrow drain scale. This is a radius, not a diameter. |
| `relief_ratio_long_radius_metres` | `float` | `8.0` | `param` | `write_relief_ratio_raster` | Long moving-window radius for local relief, in metres. | Use a broader radius that captures channel-scale topographic relief. |

### Relative Path Impact Index (Tarolli-style, D8 adaptation)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `rpii_polyfit_lambda_metres` | `float` | `30.0` | `param` | - | Surface-fitting window (lambda, metres) for the smoothed reference DEM. | Computed using direct polyfit fitted_surface. |
| `rpii_fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope used when pit-filling original and smoothed DEMs. | Passed to fill_pits() before D8 routing. |
| `rpii_boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for D8 routing (N,E,S,W). | o=open, n=no-flux. |
| `rpii_flow_routing_method` | `string` | `d8` | `param` | - | Flow routing method for drainage-area rasters in RPII. | Allowed values: d8, dinf (d-infinity). |
| `rpii_max_steps` | `int` | `400` | `param` | - | Legacy parameter retained for backwards compatibility. | Not used by the current algebraic RPII definition. |
| `write_rpii_raster` | `bool` | `true` | `switch` | - | Write Relative Path Impact Index raster. | After Tarolli et al. (2013, https://doi.org/10.5721/EuJRS20134610): RPII = ln((A_r - A_sm) / A_sm), where A_r is drainage area on the filled original DEM and A_sm is drainage area on the filled smoothed DEM. Output: &lt;write_fname&gt;_relative_path_impact_index.&lt;ext&gt;. |
| `write_rpii_log10_area_raw_raster` | `bool` | `true` | `switch` | - | Write log10(A_r) where A_r is drainage area on the filled original DEM. | Output: &lt;write_fname&gt;_rpii_log10_area_raw.&lt;ext&gt;. |
| `write_rpii_log10_area_smoothed_raster` | `bool` | `true` | `switch` | - | Write log10(A_sm) where A_sm is drainage area on the filled smoothed DEM. | Output: &lt;write_fname&gt;_rpii_log10_area_smoothed.&lt;ext&gt;. |
| `write_rpii_smoothed_surface_raster` | `bool` | `false` | `switch` | - | Write the surface-fitted DEM used by RPI. | Output: &lt;write_fname&gt;_rpii_smoothed_surface.&lt;ext&gt;. |

### Local median + relative elevation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `local_median_window_radius_metres` | `float` | `6.0` | `param` | - | Circular moving-window radius for local median, in metres. | Increase to compare each pixel against broader neighbourhood context. |
| `write_local_median_raster` | `bool` | `true` | `switch` | - | Write local median elevation raster. | Output: &lt;write_fname&gt;_local_median.&lt;ext&gt;. |
| `write_local_median_relative_elevation_raster` | `bool` | `true` | `switch` | - | Write relative elevation with respect to local median. | relative = dem - local_median. Negative = below local median; positive = above local median. Output: &lt;write_fname&gt;_local_median_relative_elevation.&lt;ext&gt;. |
| `local_median_drain_negative_threshold_metres` | `float` | `0.15` | `param` | - | Magnitude threshold (m) for drain-candidate trough pixels. | Drain candidate mask keeps pixels where relative_elevation &lt;= -threshold. |
| `local_median_top_positive_threshold_metres` | `float` | `0.15` | `param` | - | Positive threshold (m) for top-surface candidate pixels. | Top mask keeps pixels where relative_elevation &gt;= threshold. |
| `local_median_bottom_negative_threshold_metres` | `float` | `0.15` | `param` | - | Magnitude threshold (m) for bottom-surface candidate pixels. | Bottom mask keeps pixels where relative_elevation &lt;= -threshold. |
| `local_median_hag_intermix_radius_metres` | `float` | `2.0` | `param` | - | Radius (m) for local top/bottom intermix detection. | A pixel is intermix where both top and bottom masks are present within this circular neighbourhood. |
| `local_median_hag_gradient_threshold` | `float` | `0.25` | `param` | - | Gradient threshold on relative elevation for hag-transition candidates. | Transition mask keeps intermix pixels where gradient magnitude &gt;= this value. |
| `write_local_median_gradient_raster` | `bool` | `true` | `switch` | - | Write gradient magnitude of local-median relative elevation. | Output: &lt;write_fname&gt;_local_median_relative_elevation_gradient.&lt;ext&gt;. |
| `write_local_median_drain_candidate_mask_raster` | `bool` | `true` | `switch` | - | Write binary drain-candidate mask from negative relative elevation. | Output: &lt;write_fname&gt;_local_median_drain_candidate_mask.&lt;ext&gt;. |
| `write_local_median_top_mask_raster` | `bool` | `true` | `switch` | - | Write binary top-surface candidate mask. | Output: &lt;write_fname&gt;_local_median_top_mask.&lt;ext&gt;. |
| `write_local_median_bottom_mask_raster` | `bool` | `true` | `switch` | - | Write binary bottom-surface candidate mask. | Output: &lt;write_fname&gt;_local_median_bottom_mask.&lt;ext&gt;. |
| `write_local_median_hag_intermix_mask_raster` | `bool` | `true` | `switch` | - | Write binary top/bottom intermix mask. | Output: &lt;write_fname&gt;_local_median_hag_intermix_mask.&lt;ext&gt;. |
| `write_local_median_hag_transition_mask_raster` | `bool` | `true` | `switch` | - | Write binary hag-transition candidate mask (intermix + high gradient). | Output: &lt;write_fname&gt;_local_median_hag_transition_mask.&lt;ext&gt;. |

### Drain finder

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `frangi_sigma_min_metres` | `float` | `1.0` | `param` | - | Minimum Frangi Gaussian scale in metres. | Use small values for narrow drain bottoms in high-resolution DEMs. |
| `frangi_sigma_max_metres` | `float` | `2.5` | `param` | - | Maximum Frangi Gaussian scale in metres. | Set this near the expected width of the drain-bottom signal rather than the spacing between drains. |
| `frangi_n_scales` | `int` | `4` | `param` | - | Number of Frangi scales. | Scales are sampled logarithmically between the minimum and maximum Frangi sigma values. |
| `frangi_beta` | `float` | `0.4` | `param` | - | Frangi beta parameter. | Lower values more strongly prefer line-like over blob-like responses. |
| `frangi_c` | `float` | `1.0` | `param` | - | Frangi c parameter. | Controls suppression of weak-structure responses. Smaller values make the filter respond more selectively to subtle narrow drains. |
| `drain_max_frangi_scale_metres` | `float` | `1.8` | `param` | - | Maximum allowed Frangi best scale, in metres, for retained drain candidates. | Use this to suppress broader features such as channels while still running Frangi over a wider diagnostic scale range. |
| `drain_log10_abs_accumulation_curvature_blur_radius_metres` | `float` | `20.0` | `param` | - | Radius of the moving-window mean used to smooth log10(abs(accumulation_curvature)), in metres. | The blurred log10 accumulation-curvature raster is then used for accumulation-curvature masking in the drain finder. Set to 0 to disable smoothing. |
| `drain_response_threshold` | `float` | `0.15` | `param` | - | Threshold applied to the Frangi response before shape filtering. | Lower values retain more candidate pixels; higher values are more conservative. |
| `drain_min_component_pixels` | `int` | `6` | `param` | - | Minimum connected-component size in pixels. | Small spotty detections below this size are removed. |
| `drain_min_major_axis_pixels` | `float` | `5.0` | `param` | - | Minimum estimated major-axis length in pixels. | Retains only components with a sufficiently line-like span. |
| `drain_min_elongation` | `float` | `3.0` | `param` | - | Minimum component elongation = major_axis / minor_axis. | Higher values more aggressively reject leopard-spot style detections. |
| `write_drain_mask_raster` | `bool` | `true` | `switch` | - | Write the final filtered binary drain mask. | Output pixels are 1.0 for retained drain detections and 0.0 otherwise. |
| `write_drain_diagnostic_rasters` | `bool` | `true` | `switch` | - | Write intermediate drain-finder diagnostic rasters. | Includes Frangi response, Frangi best scale, and the thresholded candidate mask before shape filtering. |
| `drain_accumulation_curvature_lambda` | `float` | `4.0` | `param` | - | Wavelength used for the polyfit accumulation-curvature veto, in metres. | Use a scale that captures channel and broken-ground curvature structure without being so large that subtle drains disappear. |
| `drain_accumulation_curvature_z_factor` | `float` | `1.0` | `param` | - | Vertical scaling applied before interpreting accumulation curvature. | Usually leave at 1.0 unless you deliberately want to amplify or suppress curvature magnitudes. |
| `drain_accumulation_curvature_flat_slope_tolerance` | `float` | `1e-12` | `param` | - | Flat-slope tolerance passed to the polyfit metric engine for the accumulation-curvature veto. | Very small default so only nearly flat terrain is treated as directionless. |
| `drain_log10_abs_accumulation_curvature_max` | `float` | `-5.0` | `param` | - | Maximum allowed value of log10(abs(accumulation_curvature)) for retained drain candidates. | Pixels above this value are removed from the thresholded Frangi candidates. |

### Valley finder (SMAP shape classifier on a smoothed DEM)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `valley_polyfit_lambda_metres` | `float` | `50.0` | `param` | - | Smoothing window (polyfit lengthscale) applied before shape classification, in metres. | The quadratic-fit window IS the smoothing scale. Larger values (e.g. 50 m) suppress fine channels and retain broad valley forms; smaller values resolve narrower valleys. |
| `valley_curvature_threshold` | `float` | `0.0` | `param` | - | Curvature threshold (kt) passed to the SMAP classifier, in 1/m. | Principal curvatures with magnitude &lt;= kt are set to zero before classification, collapsing weak/long-wavelength structure toward planar or cylindrical classes. 0 keeps every pixel in a non-planar class. |
| `valley_include_basin` | `bool` | `true` | `param` | - | Include the basin class (SMAP -3, concave bowl) in the valley mask. | Concave-up class with mean curvature &gt; 0. |
| `valley_include_synformal_saddle` | `bool` | `true` | `param` | - | Include the synformal-saddle class (SMAP -2) in the valley mask. | Concave-up class with mean curvature &gt; 0 and negative Gaussian curvature. |
| `valley_include_synform` | `bool` | `true` | `param` | - | Include the synform class (SMAP -1, cylindrical valley) in the valley mask. | The canonical broad-valley form: concave across the contour, straight along it. |
| `write_valley_mask_raster` | `bool` | `true` | `switch` | - | Write the binary valley mask. | 1 indicates a pixel in one of the selected concave-up SMAP classes. |
| `write_valley_shape_class_raster` | `bool` | `true` | `switch` | - | Write the full SMAP shape-class raster. | Integer codes -4..3 stored as float32 (see surface_shape_class). |
| `write_valley_diagnostic_rasters` | `bool` | `false` | `switch` | - | Write valley-finder diagnostic rasters. | Includes the smoothed fitted surface, the principal curvatures, and the post-threshold Gaussian and mean curvatures used by the classifier. |

### Road finder (along-contour bench detector)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `road_polyfit_lambda_metres` | `float` | `8.0` | `param` | - | Smoothing window (polyfit lengthscale) for road detection, in metres. | Roads are narrow, so use a small window matched to the bench width (cut-to-fill cross-section). Larger values blur the bench edges. |
| `road_flat_slope_tolerance` | `float` | `1e-4` | `param` | - | Gradient magnitude below which the local aspect is treated as undefined. | On near-flat ground the contour direction is meaningless, so such pixels are gated out of the road score. |
| `road_min_elongation` | `float` | `0.3` | `param` | - | Minimum Hessian line-likeness for retained road pixels. | elongation = (\|lam_major\| - \|lam_minor\|) / (\|lam_major\| + \|lam_minor\|). Higher values demand a more strongly linear curvature signature. |
| `road_min_contour_alignment` | `float` | `0.7` | `param` | - | Minimum alignment between the feature's long axis and the contour direction. | \|along_feature . contour_unit\| in [0,1]. ~1 for along-contour features (roads/benches); ~0 for down-slope features (channels). Raise to be stricter. |
| `road_curvature_half_saturation` | `float` | `0.1` | `param` | - | Half-saturation curvature (1/m) for the edge-strength weight in the road score. | edge_strength = \|cross_curvature\| / (\|cross_curvature\| + this). Smaller values let gentler curvature breaks score highly; larger values demand sharper edges. Set to 0 to disable the edge-strength weighting entirely. |
| `road_min_cross_curvature` | `float` | `0.05` | `param` | - | Minimum cross-feature curvature magnitude (1/m) for retained road pixels. | 0 disables the floor. Use to suppress featureless slopes with no real edge. |
| `road_max_cross_curvature` | `float` | `0.0` | `param` | - | Maximum cross-feature curvature magnitude (1/m) for retained road pixels. | 0 disables the cap. Use to exclude deeply incised channels and cliffs that have very high cross-curvature. |
| `road_score_threshold` | `float` | `0.5` | `param` | - | Threshold applied to the road score to make a binary mask. | score = elongation * contour_alignment * edge_strength. Lower is more inclusive. |
| `road_use_bench_test` | `bool` | `true` | `param` | - | Require road pixels to be a narrow flat ribbon (bench) on a slope. | Samples slope a short distance to each side along the cross-feature direction; both flanks must be steeper than the tread. This is what rejects floodplains and valley floors, which are flat on one flank. |
| `road_bench_half_width_metres` | `float` | `6.0` | `param` | - | Cross-feature sampling distance for the bench test, in metres. | Set a little larger than the road half-width plus its cut/fill batter so the samples land on the flanking slopes, not on the tread. |
| `road_bench_flank_slope_ratio` | `float` | `1.5` | `param` | - | Bench test: required ratio of flank slope to tread slope. | Both flanks must have slope &gt;= ratio * tread_slope (and &gt;= road_bench_min_flank_slope). Higher demands a flatter tread relative to its flanks. |
| `road_bench_min_flank_slope` | `float` | `0.15` | `param` | - | Bench test: minimum absolute flank slope (tan of slope angle). | Both flanks must be at least this steep (0.15 ~ 8.5 degrees). Guards against near-flat terrain passing the ratio test when the tread slope is ~0. |
| `write_road_score_raster` | `bool` | `true` | `switch` | - | Write the continuous road score raster. | score = elongation * contour_alignment, gated by slope/curvature. |
| `write_road_mask_raster` | `bool` | `true` | `switch` | - | Write the thresholded binary road mask. | 1 indicates a plausible along-contour bench/road pixel. |
| `write_road_diagnostic_rasters` | `bool` | `false` | `switch` | - | Write road-finder diagnostic rasters. | Includes contour alignment, Hessian elongation, cross-feature curvature, along-feature orientation (azimuth from north), and slope. |

### Flow-pathway road/channel discrimination

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `flow_regional_polyfit_lambda_metres` | `float` | `100.0` | `param` | - | Smoothing window (polyfit lengthscale) for the regional slope vector, in metres. | A large window (e.g. 100 m) gives the regional fall-line bearing that flow paths are compared against. |
| `flow_path_half_window` | `int` | `5` | `param` | - | Number of pixels traced up- and downstream when measuring the flow-path bearing. | The bearing is taken over a window of (2*this + 1) pixels along the D8 path (default 5 gives an 11-pixel segment). |
| `flow_short_curvature_radius_metres` | `float` | `0.0` | `param` | - | Radius for the short-wavelength tangential curvature, in metres. | The short polyfit half-wavelength is 2 * this radius; 0 uses ~3 pixels. Output: _flow_short_tangential_curvature.&lt;ext&gt;. |
| `flow_intermediate_polyfit_lambda_metres` | `float` | `15.0` | `param` | - | Polyfit window (lengthscale) for the intermediate curvature rasters, in metres. | Produces tangential, minimum, and maximum curvature at this scale. Outputs: _flow_intermediate_tangential_curvature.&lt;ext&gt;, _flow_intermediate_minimum_curvature.&lt;ext&gt;, _flow_intermediate_maximum_curvature.&lt;ext&gt;. |
| `flow_roughness_radius_metres` | `float` | `2.5` | `param` | - | Window radius for the local roughness (elevation standard deviation), in metres. | Roads are smoother than the surrounding landscape, so they show low roughness. |
| `write_flow_discrimination_rasters` | `bool` | `true` | `switch` | - | Write the flow-pathway discrimination diagnostic rasters. | Outputs: regional slope bearing and gradient, contributing pixels, channel mask, flow-path bearing, flow-vs-slope angle, short gradient, short total curvature, roughness, bearing contrast, and the valley-axis buffer mask. |
| `flow_valley_curvature_metric` | `string` | `tangential_curvature` | `param` | - | Regional curvature used to define the valley zone. | tangential_curvature (default): cross-slope/contour curvature (flow convergence) - picks valleys but NOT footslope slope breaks (those are concave down-slope only). plan_curvature: similar but normalised differently. laplacian or mean_curvature: total concavity - also flags slope breaks. Positive = concave; ridges (convex) are always excluded by sign. |
| `flow_valley_curvature_threshold` | `float` | `0.0` | `param` | - | Regional (smoothed) curvature above which a pixel is in the valley zone, in 1/m. | Uses flow_valley_curvature_metric of the heavily smoothed surface: concave-up (&gt;= this) is the valley zone and convex-up (ridges) is excluded. 0 keeps every concave pixel; raise to restrict to stronger valley convergence. |
| `flow_valley_require_bearing_contrast` | `bool` | `false` | `param` | - | Also require the bearing-contrast gate when building the valley-axis mask. | Off by default. On its own, bearing contrast only fires on sharp bearing flips (missing broad valleys) and cannot separate valleys from ridges; the regional curvature does that. Enable only to further restrict the valley zone. |
| `flow_bearing_contrast_window_metres` | `float` | `10.0` | `param` | - | Window radius for the regional-bearing-contrast metric, in metres. | contrast = 1 - \|mean unit-vector of the regional bearing\| over this window. High where the bearing flips, e.g. two valley sides facing each other. |
| `flow_bearing_contrast_threshold` | `float` | `0.4` | `param` | - | Bearing-contrast value above which a pixel is flagged as a valley/ridge axis. | 0 (no contrast) to 1 (opposing bearings). 0.4 captures clear facing-slope axes. |
| `flow_valley_buffer_radius_metres` | `float` | `10.0` | `param` | - | Dilation buffer applied around bearing-contrast axes, in metres. | Marks the valley-axis zone where flow-path bearings are expected to be ambiguous (not road misdirection). |
| `big_channel_threshold_contributing_pixels` | `int` | `5000` | `param` | - | Contributing-pixel threshold for the valley-channel point export. | Channel pixels with more than this many upslope pixels that fall inside the valley-axis buffer are exported to a FlatGeobuf point layer. |
| `flow_divergent_export_min_angle` | `float` | `0.0` | `param` | - | Minimum flow-vs-slope angle (deg) for a channel pixel to be exported as a point. | 0 exports every channel pixel (with its angle attribute) so it can be styled in a GIS; raise to export only the divergent (road-like) pixels. |
| `write_flow_divergent_points_fgb` | `bool` | `true` | `switch` | - | Export channel pixels (with flow-vs-slope angle) to a FlatGeobuf point layer. | Output: &lt;write_fname&gt;_flow_divergent_points.fgb. Attributes: row, col, contributing_px, flow_vs_slope_angle, flow_bearing, regional_bearing. |
| `write_flow_valley_channel_points_fgb` | `bool` | `true` | `switch` | - | Export large channel pixels inside the valley-axis buffer to a FlatGeobuf point layer. | Output: &lt;write_fname&gt;_flow_valley_channel_points.fgb. Selects channel pixels with contributing pixels &gt; big_channel_threshold_contributing_pixels that lie within the bearing-contrast buffer. Attributes: row, col, contributing_px, flow_vs_slope_angle, bearing_contrast. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |

### Path opening (generic curvilinear-feature filter)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `path_min_length_pixels` | `int` | `15` | `param` | - | Minimum path length L, in pixels, for the grayscale path opening. | A pixel is kept only if it lies on a connected, near-straight path of at least L pixels in the input response. Set near the shortest real road/channel segment to keep; larger L removes more isolated speckle but also drops short fragments. |
| `path_quantization_levels` | `int` | `256` | `param` | - | Number of grayscale threshold levels for the path opening. | Higher values preserve the input value range more finely at higher cost. 256 is usually plenty for a continuous response raster. |
| `path_trim_border_pixels` | `int` | `0` | `param` | - | Trim this many pixels in from the valid-data boundary before the path opening (0 = no trim). | Removes the boundary-margin artefacts that polyfit/flow/FFT response rasters carry. Because those artefacts are laterally continuous they otherwise survive the path opening and masquerade as features. Trims from both the grid edge and any irregular nodata footprint. |
| `path_do_closing` | `bool` | `false` | `param` | - | Run a path closing instead of a path opening. | A closing fills short dark breaks across an otherwise continuous bright thread (e.g. a one-pixel gap in a road). Opening (default) suppresses isolated bright speckle that lies on no long path. |
| `write_path_opened_raster` | `bool` | `true` | `switch` | - | Write the path-opened (or path-closed) response raster. | Output: &lt;write_fname&gt;_path_opened.&lt;ext&gt;. |
| `write_path_tophat_raster` | `bool` | `true` | `switch` | - | Write the white top-hat residual (input minus path opening). | The path OPENING (_path_opened.&lt;ext&gt;) is the detector: it retains the long, laterally continuous bright threads (roads/channels) and suppresses isolated speckle. The top-hat is the complement: high where bright detail was REMOVED, i.e. the short-scale speckle and blobs that do NOT lie on a length-L path. Useful as a rejected-noise map or to flatten background before thresholding. Output: &lt;write_fname&gt;_path_tophat.&lt;ext&gt;. Only meaningful for a path opening (not a closing). |

### Structure-tensor coherence (orientation / line-likeness)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `coherence_derivative_radius_metres` | `float` | `0.0` | `param` | - | Pre-smoothing radius for the gradient in the structure tensor, in metres. | 0 = no pre-smoothing (central difference on the raw field). Use a small value (~1 pixel) to suppress gradient noise on speckly inputs. |
| `coherence_integration_radius_metres` | `float` | `2.5` | `param` | - | Tensor-averaging (integration) radius for the structure tensor, in metres. | THE key scale: the neighbourhood over which the gradient outer product is pooled. Set near the feature half-width plus a few pixels. Larger values give smoother, more confident coherence but blur close features. |
| `write_coherence_raster` | `bool` | `true` | `switch` | - | Write the structure-tensor coherence raster (line-likeness, [0,1]). | Output: &lt;write_fname&gt;_coherence.&lt;ext&gt;. |
| `write_orientation_raster` | `bool` | `true` | `switch` | - | Write the structure-tensor orientation raster (feature azimuth, [0,180) deg). | Azimuth clockwise from north along which the feature runs. Meaningful only where coherence is high. Output: &lt;write_fname&gt;_orientation_deg.&lt;ext&gt;. |
| `write_energy_raster` | `bool` | `false` | `switch` | - | Write the structure-tensor energy raster (gradient energy, lam1+lam2). | Distinguishes strong from faint oriented structure. Output: &lt;write_fname&gt;_energy.&lt;ext&gt;. |

### Track spacing (local directional periodicity)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `spacing_window_pixels` | `int` | `64` | `param` | - | Analysis window side length for the spacing detector, in pixels. | Must span several periods of the largest spacing of interest (a few times spacing_max_metres / cell_size). |
| `spacing_step_pixels` | `int` | `32` | `param` | - | Window stride for the spacing detector, in pixels. | Smaller = smoother and slower; typically half the window size. |
| `spacing_min_metres` | `float` | `5.0` | `param` | - | Shortest repeat spacing to search for, in metres. | Must be &gt;= 2 * cell_size (Nyquist). |
| `spacing_max_metres` | `float` | `30.0` | `param` | - | Longest repeat spacing to search for, in metres. | Must be &lt;= spacing_window_pixels * cell_size. |
| `spacing_min_valid_fraction` | `float` | `0.6` | `param` | - | Minimum fraction of valid pixels for a window to be analysed. | Windows with more nodata than this are skipped. |
| `write_spacing_strength_raster` | `bool` | `true` | `switch` | - | Write the periodicity strength (peak prominence) raster. | High where a sharp regular spacing dominates. Output: &lt;write_fname&gt;_spacing_strength.&lt;ext&gt;. |
| `write_spacing_raster` | `bool` | `true` | `switch` | - | Write the dominant spacing raster, in metres. | Use to classify feature families (ditches vs tracks vs wave dams). Output: &lt;write_fname&gt;_spacing_metres.&lt;ext&gt;. |
| `write_spacing_orientation_raster` | `bool` | `false` | `switch` | - | Write the spacing orientation raster (line azimuth, [0,180) deg). | Output: &lt;write_fname&gt;_spacing_orientation_deg.&lt;ext&gt;. |

### Ridge lines (oriented NMS + hysteresis)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `ridge_high_threshold` | `float` | `40.0` | `param` | - | Strong-seed threshold (in input response units) for ridge hysteresis. | Thinned-response pixels at or above this start a retained centreline. Set relative to the input raster's value range. |
| `ridge_low_threshold` | `float` | `20.0` | `param` | - | Weak-link threshold (in input response units) for ridge hysteresis. | Thinned pixels at or above this are kept only if connected to a strong seed. Must be &lt;= ridge_high_threshold. |
| `ridge_min_coherence` | `float` | `0.0` | `param` | - | Minimum structure-tensor coherence to admit a pixel to the ridge response. | Pixels with coherence below this are zeroed before NMS, removing ridge maxima that sit in incoherent speckle. 0 disables the gate. Uses the coherence_*_radius_metres parameters for the tensor scale. |
| `write_ridge_thinned_raster` | `bool` | `true` | `switch` | - | Write the NMS-thinned ridge response (1-pixel centrelines, response-valued). | Output: &lt;write_fname&gt;_ridge_thinned.&lt;ext&gt;. |
| `write_ridge_mask_raster` | `bool` | `true` | `switch` | - | Write the binary hysteresis-linked ridge centreline mask. | Enabled by default because this is a primary product of the ridge-finder stage, not the optional hilltop diagnostic exposed by the coupling CLI. Output: &lt;write_fname&gt;_ridge_mask.&lt;ext&gt;. |

### Thin-line extractor (separate lines from blobs/speckle)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `thin_threshold` | `float` | `20.0` | `param` | - | Threshold to binarise the input response (keep response &gt;= this). | - |
| `thin_directional_closing_length_pixels` | `int` | `0` | `param` | - | Line length (pixels) for directional closing to bridge gaps along lines (0 = skip). | Reconnects fragmented features before component analysis. |
| `thin_path_length_pixels` | `int` | `40` | `param` | - | Binary path-opening length L (pixels); 0 = skip. | THE primary thin-vs-blob/speckle tool: keeps pixels on a connected thin path of at least L pixels (handles curves); removes speckle and any blob lacking such a path. Set near the shortest road to keep. |
| `thin_blob_removal_radius_metres` | `float` | `0.0` | `param` | - | Disk-opening top-hat radius (metres) to remove blobs wider than this (0 = skip). | THE thin-vs-blob filter: features wider than ~2x this radius (blobs, wide patches) are removed; thinner ribbons (roads) survive. Set a little above the widest real road half-width. |
| `thin_min_component_pixels` | `int` | `40` | `param` | - | Minimum connected-component size (pixels). | - |
| `thin_min_major_axis_pixels` | `float` | `0.0` | `param` | - | Minimum component major-axis length (pixels); 0 = skip. | - |
| `thin_min_elongation` | `float` | `1.0` | `param` | - | Minimum component elongation = major/minor axis; 1 = skip. | Raise (e.g. 3) to keep only clearly linear components and drop compact residue. |
| `thin_trim_border_pixels` | `int` | `0` | `param` | - | Trim this many pixels from the final mask border (0 = none). | Removes linear edge artefacts that survive the elongation filter. |
| `write_thin_mask_raster` | `bool` | `true` | `switch` | - | Write the final thin-line mask. | Output: &lt;write_fname&gt;_thin_mask.&lt;ext&gt;. |
| `write_thin_diagnostic_rasters` | `bool` | `false` | `switch` | - | Write the thresholded binary and the blob-removed intermediate. | - |

### Generic roughness / consistency of an input raster

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `roughness_input_radius_metres` | `float` | `5.0` | `param` | - | Moving-window radius for the roughness/consistency metric, in metres. | Larger windows pool more context (smoother result). |
| `roughness_circular` | `bool` | `false` | `param` | - | Treat the input as a circular angle in degrees (e.g. a flow bearing). | When true, computes circular consistency R = \|mean unit-vector\| over the window (1 = perfectly consistent direction, e.g. flow channelled along a road; 0 = directions cancel). Writes consistency (R) and roughness (1-R). When false, computes the linear moving-window standard deviation. |
| `write_roughness_raster` | `bool` | `true` | `switch` | - | Write the roughness raster (stddev, or 1-R for circular input). | Output: &lt;write_fname&gt;_input_roughness.&lt;ext&gt;. |
| `write_consistency_raster` | `bool` | `true` | `switch` | - | Write the circular consistency raster R (only for roughness_circular=true). | High where the angle is locally consistent. Output: &lt;write_fname&gt;_input_consistency.&lt;ext&gt;. |

### Smoothness road extractor (roads = smooth ribbons)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `smooth_roughness_radius_metres` | `float` | `1.5` | `param` | - | Window radius for the local-relief roughness, in metres. | Match to the hillslope speckle scale. Roads are smooth (low roughness) ribbons; the surrounding terrain is rough. |
| `smooth_roughness_half_saturation` | `float` | `0.1` | `param` | - | Roughness half-saturation (metres) for smoothness = 1/(1+roughness/h). | Roughness equal to this gives smoothness 0.5. Lower = more selective for very smooth surfaces; set near the road-vs-hillslope roughness boundary. |
| `smooth_trim_border_pixels` | `int` | `10` | `param` | - | Border trim applied before extraction (0 = none). | - |
| `smooth_path_min_length_pixels` | `int` | `40` | `param` | - | Minimum path length L (pixels) for the smoothness path opening. | Keeps smooth ribbons at least this long; drops short smooth patches. |
| `smooth_threshold` | `float` | `0.5` | `param` | - | Threshold on the path-opened smoothness to make the binary ribbon mask. | Higher = only the smoothest ribbons; lower = more inclusive. |
| `smooth_min_component_pixels` | `int` | `40` | `param` | - | Minimum connected-component size (pixels) for a retained road. | - |
| `smooth_min_major_axis_pixels` | `float` | `0.0` | `param` | - | Minimum component major-axis length (pixels); 0 = skip. | Use to demand a minimum linear span. |
| `smooth_min_elongation` | `float` | `1.0` | `param` | - | Minimum component elongation = major/minor axis; 1 = skip. | Raise (e.g. 3) to reject smooth blobs (fields) and keep thin roads. |
| `write_smooth_road_mask_raster` | `bool` | `true` | `switch` | - | Write the final smoothness-based road mask. | Output: &lt;write_fname&gt;_smooth_road_mask.&lt;ext&gt;. |
| `write_smooth_diagnostic_rasters` | `bool` | `false` | `switch` | - | Write roughness, smoothness, and path-opened smoothness diagnostics. | - |

### Road extractor (composed pipeline)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `road_trim_border_pixels` | `int` | `15` | `param` | - | Border trim applied to the response before road extraction (0 = none). | Removes boundary-margin artefacts that otherwise survive as linear features. |
| `road_min_coherence` | `float` | `0.5` | `param` | - | Minimum structure-tensor coherence to admit a pixel to the ridge response. | Pixels below this are zeroed before NMS, removing ridge maxima in incoherent speckle. 0 disables the gate. |
| `road_min_component_pixels` | `int` | `30` | `param` | - | Minimum connected-component size (pixels) for a retained road. | Drops short speckle fragments left by the ridge detector while keeping long road chains. Set near the shortest real road length in pixels. |
| `road_valley_curvature_fname` | `string` | - | `param` | - | Optional regional-curvature raster for the valley veto (relative to read path). | Typically &lt;prefix&gt;_flow_regional_valley_curvature.&lt;ext&gt; from flow_road_discrimination. Empty = no valley veto. Centreline pixels where this raster &gt;= road_valley_curvature_veto_threshold are removed. |
| `road_valley_curvature_veto_threshold` | `float` | `0.0` | `param` | - | Curvature threshold for the valley veto (sign per the supplied raster). | Centreline pixels with regional curvature &gt;= this (concave-up valley zone) are dropped, removing the in-valley channel false positive. |
| `road_classify_spacing` | `bool` | `false` | `param` | - | Classify each extracted component as road (1) or regularly-spaced track (2). | Runs the directional-periodicity detector on the response and labels a component a track when at least road_track_min_fraction of its pixels lie in a field with periodicity strength &gt;= road_track_min_strength. Uses the spacing_* parameters. Output: &lt;write_fname&gt;_road_extract_classified.&lt;ext&gt;. |
| `road_track_min_strength` | `float` | `6.0` | `param` | - | Periodicity strength above which a pixel counts as 'in a regular field'. | See track_spacing strength (whitened directional anisotropy). Raise to demand stronger regularity before calling a component a track. |
| `road_track_min_fraction` | `float` | `0.5` | `param` | - | Fraction of a component's pixels that must be 'in a regular field' to label it a track. | Higher = stricter (only components mostly inside a regular field become tracks). Ignored when road_classify_per_pixel is true. |
| `road_classify_per_pixel` | `bool` | `false` | `param` | - | Classify per pixel instead of per connected component. | false (default): whole component voted road/track (consistent per feature but diluted where roads and tracks connect). true: each pixel is a track where its local periodicity strength &gt;= road_track_min_strength (no dilution; road_track_min_fraction is ignored). |
| `write_road_extract_mask_raster` | `bool` | `true` | `switch` | - | Write the final extracted road mask. | Output: &lt;write_fname&gt;_road_extract_mask.&lt;ext&gt;. |
| `write_road_extract_diagnostic_rasters` | `bool` | `false` | `switch` | - | Write road-extractor intermediate rasters. | Diffused response, coherence, NMS-thinned ridge, and the pre-cleanup ridge mask. |

### Coherence-enhancing diffusion

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `ced_iterations` | `int` | `10` | `param` | `write_diffused_raster` | Number of coherence-enhancing diffusion iterations. | More iterations bridge longer gaps and smooth more strongly along feature orientation. |
| `ced_time_step` | `float` | `0.15` | `param` | `write_diffused_raster` | Explicit diffusion time step (tau). | Must be in (0, 0.25] for stability. 0.15 is a safe default. |
| `ced_derivative_radius_metres` | `float` | `0.5` | `param` | `write_diffused_raster` | Gradient pre-smoothing radius for the diffusion structure tensor, in metres. | 0 = none. A small value steadies the orientation estimate on noisy input. |
| `ced_integration_radius_metres` | `float` | `3.0` | `param` | `write_diffused_raster` | Structure-tensor integration radius for the diffusion, in metres. | Sets the scale over which a coherent orientation is detected before smoothing along it. |
| `ced_contrast` | `float` | `0.2` | `param` | `write_diffused_raster` | Coherence contrast (C) controlling along-feature smoothing selectivity. | Uses normalised coherence in [0,1], so C is scale-invariant. Higher = smooth only the most coherent threads; lower = smooth more of the field. |
| `ced_min_diffusivity` | `float` | `0.001` | `param` | `write_diffused_raster` | Floor (across-feature) diffusivity for the diffusion tensor. | Small value in (0,1); controls the tiny amount of smoothing allowed across the feature. |
| `write_diffused_raster` | `bool` | `true` | `switch` | - | Write the coherence-enhanced (diffused) raster. | Output: &lt;write_fname&gt;_diffused.&lt;ext&gt;. |

<!-- END GENERATED CLI PARAMETERS -->
