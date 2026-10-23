# CLI reference: lsdtt-surface-metrics

Extracts surface statistics and derived rasters from a DEM without requiring flow routing. All analyses operate directly on the input raster. The DEM must already be in a projected (metre-unit) CRS.

## Usage

```
lsdtt-surface-metrics <path> <parameter_file>
lsdtt-surface-metrics <parameter_file>         # path defaults to ./
lsdtt-surface-metrics -h                       # print help and exit
lsdtt-surface-metrics -v                       # print version and exit
lsdtt-surface-metrics -c                       # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Pipeline overview

| Step | Trigger |
|---|---|
| 1. Hillshade | `write_hillshade_raster: true` |
| 2. Local relief | `write_local_relief_raster: true` |
| 3. Polynomial surface fit + derived metrics | `write_polyfit_rasters: true` or `polyfit_metrics` is non-empty |
| 3b. Surface shape classification (SMAP) | `write_shape_classification_raster: true` |
| 4. MADk2 roughness | `write_madk2_rasters: true` |
| 5. Frangi filter | `write_frangi_raster: true` or `write_frangi_scale_raster: true` |
| 6. Otsu threshold segmentation | `write_otsu_raster: true` |
| 7. Hough line extraction | `write_hough_lines_raster: true` |
| 8. Curvature skeletons | `write_curvature_skeleton_rasters: true` |
| 9. VEC adaptive thresholding | `write_vec_adaptive_threshold_rasters: true` |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input raster filename. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `write_output_format` | string | `cog` | Output raster format. `cog` = Cloud Optimised GeoTIFF. `gtiff` = plain compressed GeoTIFF. `same_as_input` = match the input format. `png` = 8-bit scaled image with a worldfile and PAM sidecar, for quick visualisation only, not for analysis reuse. Output files take the extension of the chosen format. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` uses all hardware threads (default); `1` is serial. Applies to wavelet/direct polynomial fits, the polyfit metric engine, Frangi, and the 8-bit image-operation passes. Set an explicit limit no greater than your physical-core allocation in shared environments. |

---

## Step 1: Hillshade

Output: `<write_fname>_hillshade.tif`

Computed directly from the input DEM — no pit fill is applied first.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_hillshade_raster` | bool | `false` | Enable hillshade computation. Uses Zevenbergen and Thorne (1987). |
| `hillshade_azimuth` | float | `315` | Sun azimuth in degrees clockwise from north. 315° (northwest) is the standard cartographic convention. |
| `hillshade_altitude` | float | `45` | Sun altitude above the horizon in degrees. |
| `hillshade_z_factor` | float | `1.0` | Vertical exaggeration applied before computing the hillshade. Values greater than 1 emphasise relief in flat terrain. |

---

## Step 2: Local relief

Output: `<write_fname>_local_relief.tif`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_local_relief_raster` | bool | `false` | Enable local relief computation (local max − local min within a moving window over valid pixels). |
| `local_relief_radius_metres` | float | `10.0` | Radius of the moving window in metres. |
| `local_relief_window_shape` | string | `circular` | Window shape. `circular` or `square`. |

---

## Step 3: Polynomial surface fit

Fits a quadratic polynomial `z = ax² + by² + cxy + dx + ey + f` over a moving window and optionally derives surface metrics from the coefficients.

### Fit algorithm

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_polyfit_rasters` | bool | `false` | Write the six raw coefficient rasters (`a`–`f`) for the fitted surface. |
| `polyfit_lambda_metres` | float | `70.0` | Characteristic lengthscale in metres controlling the moving window size. Default (70 m) is tuned for Copernicus 30 m data. Reduce for lidar (e.g. 10–30 m). |
| `polyfit_method` | string | `direct` | Algorithm for the fit. `direct` (default, preferred): separable rolling-moment O(n) fit (no FFTW needed); faster and better at resolving fine-scale features. `wavelet`: FFT-based Gaussian wavelet fit (requires FFTW). `both`: run both algorithms and print a timing comparison. Output files are prefixed `_wavelet_polyfit_*` or `_direct_polyfit_*` so both sets coexist. |
| `polyfit_z_factor` | float | `1.0` | Vertical scaling applied before interpreting derivatives. Scales first and second derivatives consistently in the metric engine. |
| `polyfit_flat_slope_tolerance` | float | `1e-12` | Slope magnitude below which direction-dependent metrics (aspect, plan curvature, profile curvature, tangential curvature, rotor) are set to nodata. |

### Derived metrics

Set `polyfit_metrics` to a comma-separated list of metric names, or `all` to write every metric.

> **No spaces in the list.** Write `slope,aspect_degrees,mean_curvature`, not
> `slope, aspect_degrees, mean_curvature`. The parameter parser terminates a
> value at the first space, so a spaced list is silently truncated to its first
> entry — the run succeeds and quietly computes only that one metric. The same
> applies to `slope_thresholds` and to comma-separated lists in other CLIs
> (`m_over_n`, `gsd_fractions`, `psi_mm`).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `polyfit_metrics` | string | _(empty)_ | Comma-separated list of derived metrics to compute and write. Empty = none. |

Available metric names and their output suffixes:

> **Aspect convention.** `aspect_radians` and `aspect_degrees` report the
> **down-slope** azimuth — the direction the slope *faces* — clockwise from
> north, so 0 = faces north, 90 = faces east. This matches GDAL, ArcGIS and the
> geomorphology literature. **Changed 2026-07-22**: before that release the
> value was the *up-slope* azimuth, i.e. 180° opposite, so rasters produced by
> earlier versions must be rotated by 180° before being compared with new ones.


| Metric name | Output suffix |
|---|---|
| `fitted_surface` | `_<method>_polyfit_fitted_surface.tif` |
| `dzdx` | `_<method>_polyfit_dzdx.tif` |
| `dzdy` | `_<method>_polyfit_dzdy.tif` |
| `d2zdx2` | `_<method>_polyfit_d2zdx2.tif` |
| `d2zdxdy` | `_<method>_polyfit_d2zdxdy.tif` |
| `d2zdy2` | `_<method>_polyfit_d2zdy2.tif` |
| `slope` | `_<method>_polyfit_slope.tif` |
| `slope_degrees` | `_<method>_polyfit_slope_degrees.tif` |
| `aspect_radians` | `_<method>_polyfit_aspect_radians.tif` |
| `aspect_degrees` | `_<method>_polyfit_aspect_degrees.tif` |
| `laplacian` | `_<method>_polyfit_laplacian.tif` |
| `hessian_determinant` | `_<method>_polyfit_hessian_determinant.tif` |
| `mean_curvature` | `_<method>_polyfit_mean_curvature.tif` |
| `gaussian_curvature` | `_<method>_polyfit_gaussian_curvature.tif` |
| `maximum_curvature` | `_<method>_polyfit_maximum_curvature.tif` |
| `minimum_curvature` | `_<method>_polyfit_minimum_curvature.tif` |
| `curvedness` | `_<method>_polyfit_curvedness.tif` |
| `shape_index` | `_<method>_polyfit_shape_index.tif` |
| `profile_curvature` | `_<method>_polyfit_profile_curvature.tif` |
| `plan_curvature` | `_<method>_polyfit_plan_curvature.tif` |
| `tangential_curvature` | `_<method>_polyfit_tangential_curvature.tif` |
| `difference_curvature` | `_<method>_polyfit_difference_curvature.tif` |
| `horizontal_excess_curvature` | `_<method>_polyfit_horizontal_excess_curvature.tif` |
| `vertical_excess_curvature` | `_<method>_polyfit_vertical_excess_curvature.tif` |
| `accumulation_curvature` | `_<method>_polyfit_accumulation_curvature.tif` |
| `ring_curvature` | `_<method>_polyfit_ring_curvature.tif` |
| `rotor` | `_<method>_polyfit_rotor.tif` |

Where `<method>` is `wavelet` or `direct` depending on `polyfit_method`.

---

## Step 3b: Surface shape classification (SMAP)

Output: `<write_fname>_<method>_shape_class.tif`

A morphological surface-shape classification ported from the differential-geometry approach of TopoCurve (`CurveObj.m` / `CurveCalc`). Each pixel is labelled by the signs of Gaussian curvature (`KG = k₁·k₂`) and mean curvature (`KM = ½(k₁+k₂)`) formed from the two principal curvatures, after an optional threshold `kt` is applied to each principal curvature. This step builds on the polynomial fit, so it runs whenever `write_shape_classification_raster: true` (the polyfit is triggered automatically and the required curvature metrics are added to `polyfit_metrics` for you).

The classification uses the lsdtt3 curvature sign convention: **convex-up terrain (ridge, peak) → negative mean curvature; concave-up terrain (channel, basin) → positive mean curvature.**

### Class codes

The `shape_class` raster stores these integer codes (as float32; nodata elsewhere):

| Code | Class | KG | KM | Character |
|---|---|---|---|---|
| `3` | Peak / dome | > 0 | < 0 | Convex in all directions |
| `2` | Antiformal saddle | < 0 | < 0 | Saddle, net convex |
| `1` | Antiform | = 0 | < 0 | Cylindrical ridge |
| `0` | Plane | = 0 | = 0 | Flat |
| `-1` | Synform | = 0 | > 0 | Cylindrical valley |
| `-2` | Synformal saddle | < 0 | > 0 | Saddle, net concave |
| `-3` | Basin | > 0 | > 0 | Concave in all directions |
| `-4` | Perfect saddle | < 0 | = 0 | Symmetric saddle |

With `kt = 0` only the four generic classes (`3`, `2`, `-2`, `-3`) appear, because the cylindrical/planar classes require a principal curvature of exactly zero. Increasing `kt` zeroes weak principal curvatures and so introduces the `1`, `-1`, and `0` classes — use it to suppress low-amplitude, long-wavelength structure (e.g. a background basin).

### Curvature method

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_shape_classification_raster` | bool | `false` | Enable the shape classification. Requires the polynomial fit; `maximum_curvature`/`minimum_curvature` (polyfit method) or `fitted_surface` (fundamental_form method) are added to `polyfit_metrics` automatically. |
| `shape_curvature_method` | string | `polyfit` | Source of the principal curvatures. `polyfit`: analytic principal curvatures from the polynomial fit. `fundamental_form`: differential-geometry curvatures ported from TopoCurve `CurveCalc`, computed from the polyfit `fitted_surface` via the first and second fundamental forms (gradients of the surface and of its unit normal field). The smoothing scale is set by `polyfit_lambda_metres`. |
| `shape_curvature_threshold` | float | `0.0` | Curvature threshold `kt` (units: 1/m). Principal curvatures with magnitude ≤ `kt` are set to zero before `KG`/`KM` are formed. |
| `write_shape_curvature_rasters` | bool | `false` | Also write the post-threshold Gaussian and mean curvature rasters used for the classification. These differ from the raw `gaussian_curvature`/`mean_curvature` polyfit metrics because they reflect `kt`. |

### Outputs

| File | When |
|---|---|
| `_<method>_shape_class.tif` | Always (SMAP integer codes) |
| `_<method>_ff_maximum_curvature.tif`, `_<method>_ff_minimum_curvature.tif` | `shape_curvature_method: fundamental_form` (the fundamental-form principal curvatures) |
| `_<method>_shape_gaussian_curvature.tif`, `_<method>_shape_mean_curvature.tif` | `write_shape_curvature_rasters: true` (post-threshold KG/KM) |

Here `<method>` is `wavelet` or `direct` (the `polyfit_method`); for `fundamental_form` it denotes the polyfit method used to produce the smoothed `fitted_surface`.

---

## Step 3c: Slope threshold layers

Writes one raster per threshold showing where `slope_degrees` exceeds that
threshold. Requires the polynomial fit; `slope_degrees` is added to
`polyfit_metrics` automatically.

Outputs: `<write_fname>_<method>_slope_threshold_<threshold>.tif`, one per value.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_slope_threshold_layers` | bool | `false` | Enable the slope-threshold layers. |
| `slope_thresholds` | string | `3,7,12` | Comma-separated threshold values in **degrees**. Defaults to the categories used by the Environment Agency's ALERT tool. **No spaces** — see the warning under Derived metrics. |
| `clamp_slope_to_threshold` | bool | `true` | `true`: pixels above the threshold take the threshold value. `false`: they keep their original `slope_degrees` value. |

---

## Step 3d: Cloth-drape envelope

Runs a Verlet cloth simulation: a flat sheet starting at the maximum DEM
elevation falls under gravity and comes to rest on ridges and summits. The
resting surface is an upper envelope of the topography, and the gap between it
and the DEM is a measure of valley incision.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_cloth_envelope_raster` | bool | `false` | Write the resting cloth surface. Output `_cloth_envelope.tif`. |
| `write_cloth_incision_raster` | bool | `false` | Write envelope minus DEM. Positive = how far terrain lies below the envelope (incision depth). Output `_cloth_incision.tif`. |
| `write_cloth_crest_mask_raster` | bool | `false` | Write the contact mask: `1` where a cloth particle was pinned by terrain (ridge or summit), `0` where it settled freely. Useful for diagnosing whether the cloth is too rigid or too slack. Output `_cloth_crest_mask.tif`. |
| `cloth_rigidity` | float | `1.0` | Stiffness of the sheet: higher spans wider valleys without sagging. |
| `cloth_gravity` | float | `0.3` | Step size for the fall. |
| `cloth_max_iterations` | int | `500` | Iteration cap. |
| `cloth_convergence_tolerance_metres` | float | `0.005` | Movement below which the cloth is considered settled. |
| `cloth_class_threshold` | float | `0.5` | Contact distance for the crest mask. |

Run `lsdtt-surface-metrics -h <parameter>` for the current defaults and the
full guidance on each.

---

## Step 4: MADk2 roughness

Computes omnidirectional roughness (OR), direction of maximum continuity (DMC), and roughness anisotropy (RA) following Trevisani et al. (2023) using second-order directional differences.

Outputs: `<write_fname>_madk2_or.tif`, `_madk2_dmc.tif`, `_madk2_ra.tif`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_madk2_rasters` | bool | `false` | Enable MADk2 roughness computation. |
| `madk2_window_radius_metres` | float | `10.0` | Circular moving-window radius in metres. |
| `madk2_lag_pixels` | int | `1` | Lag in pixels for the second-order directional differences. Must be ≥ 1. Typical values are 1 or 2. |

---

## Step 5: Frangi filter

Enhances line-like terrain features (ridges or valleys) over multiple Gaussian smoothing scales using the eigenvalues of the Hessian matrix. The maximum response across all scales is retained per pixel.

Outputs: `<write_fname>_frangi.tif`, `<write_fname>_frangi_scale.tif`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_frangi_raster` | bool | `false` | Write the Frangi filter response raster. |
| `write_frangi_scale_raster` | bool | `false` | Write the scale-at-maximum-response raster (physical scale in metres). |
| `frangi_sigma_min_metres` | float | `1.0` | Minimum Gaussian smoothing scale in metres. |
| `frangi_sigma_max_metres` | float | `5.0` | Maximum Gaussian smoothing scale. This general terrain workflow searches broader lineaments up to 5 m; the drain-specific feature finder stops at 2.5 m. |
| `frangi_n_scales` | int | `4` | Number of scales sampled logarithmically between min and max. |
| `frangi_beta` | float | `0.5` | Controls suppression of blob-like structure. The general terrain default is 0.5; the drain finder uses 0.4 to prefer narrow lines more strongly. |
| `frangi_c` | float | `1.0` | Controls suppression of low-structureness responses. Phase 5 intentionally changes the scientific default from 15 to 1 under approved decision D7, matching `lsdtt-feature-finder`; set `15` explicitly to reproduce the former response. |
| `frangi_target` | string | `ridges` | Feature polarity. `ridges` or `valleys`. |

---

## Step 6: Otsu threshold segmentation

Output: `<write_fname>_otsu.tif`

Normalises the DEM to 8-bit and applies Otsu's method to select the threshold that minimises intra-class greyscale variance. Output is a binary raster: 1.0 = foreground (above threshold), 0.0 = background.

The 8-bit scaling uses the minimum and maximum of the **valid** pixels only, and the Otsu histogram counts valid pixels only, so nodata does not affect the threshold. (Before these operations were reimplemented natively, nodata entered both as zeros; on a DEM with a nodata border that made every valid pixel foreground.)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_otsu_raster` | bool | `false` | Enable Otsu segmentation. |

---

## Step 7: Hough line extraction

Output: `<write_fname>_hough_lines.tif`

Applies the progressive probabilistic Hough transform (Matas et al. 2000) to a Canny edge image derived from the DEM, and draws the detected segments 1 pixel wide. Output is a binary raster marking detected line pixels.

The algorithm and parameters are those of OpenCV's `HoughLinesP`. The order in which edge points are visited is random; lsdtt3 uses its own fixed-seed generator, so the result is identical on every run and platform, but the exact segments differ from OpenCV's.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_hough_lines_raster` | bool | `false` | Enable Hough line extraction. |
| `hough_rho` | float | `1.0` | Distance resolution of the Hough accumulator in pixels. |
| `hough_theta_degrees` | float | `1.0` | Angle resolution of the accumulator in degrees. |
| `hough_threshold` | int | `50` | Minimum accumulator votes required to accept a line. Increase to suppress short or faint lines. |
| `hough_min_line_length` | float | `30.0` | Minimum line segment length in pixels. |
| `hough_max_line_gap` | float | `10.0` | Maximum gap in pixels between collinear points to bridge into one segment. |
| `canny_threshold` | float | `50.0` | Lower Canny hysteresis threshold in 8-bit greyscale space [0, 255]. |
| `canny_ratio` | float | `2.0` | Ratio of upper to lower Canny threshold. Canny recommends 2.0–3.0. |
| `canny_gaussian_sigma` | float | `1.0` | Sigma for Gaussian blur applied before Canny edge detection. `0.0` = no blur. |

---

## Step 8: Curvature skeletons

Extracts slope-break skeletons from the vertical excess curvature (VEC) and the local relief of profile curvature, using a `direct` polynomial fit at `curvature_skeleton_lambda_metres`.

Outputs:

| File | Description |
|---|---|
| `<write_fname>_vec_skeleton_curvature.tif` | Raw VEC raster |
| `<write_fname>_vec_skeleton_binary.tif` | Binary VEC threshold mask |
| `<write_fname>_vec_skeleton.tif` | Thinned VEC skeleton |
| `<write_fname>_profile_curvature.tif` | Raw profile curvature raster |
| `<write_fname>_profile_curvature_relief.tif` | Local relief of profile curvature |
| `<write_fname>_profile_skeleton_binary.tif` | Binary profile curvature relief mask |
| `<write_fname>_profile_skeleton.tif` | Thinned profile curvature skeleton |
| `<write_fname>_vec_skeleton_cdf.csv` | VEC CDF (when `vec_skeleton_write_cdf: true`) |

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_curvature_skeleton_rasters` | bool | `false` | Enable curvature skeleton extraction. |
| `curvature_skeleton_lambda_metres` | float | `100.0` | Polyfit half-wavelength in metres for the direct fit. Larger values suppress noise; smaller values resolve finer features. |
| `vec_skeleton_threshold` | float | `0.0` | VEC value above which pixels are classified as slope-break candidates. Ignored when `vec_skeleton_use_cdf_threshold` is true. Run once and read the printed CDF divergence threshold to tune this. |
| `vec_skeleton_use_cdf_threshold` | bool | `false` | Auto-select the VEC threshold using CDF divergence. Overrides `vec_skeleton_threshold`. |
| `vec_skeleton_write_cdf` | bool | `false` | Write the VEC CDF to `_vec_skeleton_cdf.csv`. |
| `profile_relief_window_metres` | float | `50.0` | Circular window radius (m) for computing local relief of profile curvature. |
| `profile_relief_threshold` | float | `0.1` | Threshold on profile curvature local relief for binarisation. |
| `curvature_skeleton_spur_max_pixels` | int | `10` | Maximum endpoint-rooted spur length (pixels) pruned from skeletons after thinning. `0` = skip pruning. |
| `curvature_skeleton_min_component_pixels` | int | `1` | Minimum component size retained in skeletons after thinning and pruning. |

---

## Step 9: VEC adaptive thresholding

Applies adaptive mean and Gaussian thresholding to the vertical excess curvature. Each pixel is compared to the local neighbourhood mean (or Gaussian-weighted mean) minus the constant `vec_adaptive_C`.

Outputs: `<write_fname>_vec_adaptive_curvature.tif`, `_vec_adaptive_mean.tif`, `_vec_adaptive_gaussian.tif`, and optionally their skeletons.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_vec_adaptive_threshold_rasters` | bool | `false` | Enable VEC adaptive thresholding. Uses `curvature_skeleton_lambda_metres` for the polyfit scale. |
| `vec_adaptive_block_size_pixels` | int | `51` | Neighbourhood block size in pixels. Must be odd and ≥ 3. Should exceed the typical feature width but be smaller than the between-feature spacing. |
| `vec_adaptive_C` | float | `0.0` | Constant subtracted from the local mean before thresholding. `0` = threshold at the local mean. Negative (e.g. `-5`) = stricter; positive = looser. |
| `vec_adaptive_gaussian_sigma` | float | `0.0` | Sigma for optional Gaussian pre-blur before thresholding. `0` = no blur. |
| `write_vec_adaptive_skeleton_rasters` | bool | `false` | Thin the adaptive binary masks to 1-pixel-wide Zhang-Suen skeletons. Uses `curvature_skeleton_spur_max_pixels` and `curvature_skeleton_min_component_pixels`. Only active when `write_vec_adaptive_threshold_rasters` is true. |

---

## Example parameter file

```
read_path:  /data/dems/
read_fname: my_dem
write_path: /data/output/
write_fname: my_dem_metrics

n_threads: 0

# Hillshade
write_hillshade_raster: true
hillshade_azimuth:      315
hillshade_altitude:     45
hillshade_z_factor:     1.0

# Polynomial surface fit — slope and curvatures at 30 m scale
polyfit_lambda_metres: 70.0
polyfit_method:        direct
polyfit_metrics:       slope,plan_curvature,profile_curvature,vertical_excess_curvature

# Surface shape classification (SMAP)
write_shape_classification_raster: true
shape_curvature_method:        fundamental_form
shape_curvature_threshold:     0.0
write_shape_curvature_rasters: true

# MADk2 roughness
write_madk2_rasters:          true
madk2_window_radius_metres:   20.0
madk2_lag_pixels:             1

# Frangi ridge filter
write_frangi_raster:      true
frangi_sigma_min_metres:  50.0
frangi_sigma_max_metres:  200.0
frangi_n_scales:          4
frangi_target:            ridges
```

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

### Threading

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all hardware threads (default). 1 = serial. In shared environments, set an explicit limit no greater than the physical-core allocation. Applies to the direct and wavelet quadratic fits, polyfit metric engine, and the 8-bit image-operation passes (Otsu, adaptive thresholds, Canny, Hough). |

### Output format

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_output_format` | `string` | `cog` | `io` | - | Output raster format. | Allowed values: cog, gtiff, same_as_input, png. cog: Cloud Optimised GeoTIFF (default). gtiff: plain compressed GeoTIFF. same_as_input: match the input raster format. png: 8-bit scaled image with a worldfile and PAM sidecar for georeferencing - for quick visualisation only, not for analysis reuse. Output files take the extension of the chosen format. |

### Hillshade

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_hillshade_raster` | `bool` | `false` | `switch` | - | Write a hillshade raster. | Uses Zevenbergen and Thorne (1987). Computed directly from the input DEM - no fill is applied. |
| `hillshade_azimuth` | `float` | `315` | `param` | `write_hillshade_raster` | Sun azimuth in degrees clockwise from north. | 315 (northwest) is the cartographic convention. |
| `hillshade_altitude` | `float` | `45` | `param` | `write_hillshade_raster` | Sun altitude above the horizon in degrees. | 45 degrees is a typical default. |
| `hillshade_z_factor` | `float` | `1.0` | `param` | `write_hillshade_raster` | Vertical exaggeration factor applied before computing hillshade. | Values greater than 1 emphasise relief in flat terrain. |

### Local relief

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_local_relief_raster` | `bool` | `false` | `switch` | - | Write a local relief raster. | Local relief is computed as local_max - local_min within a moving window over valid pixels only. |
| `local_relief_radius_metres` | `float` | `10.0` | `param` | `write_local_relief_raster` | Radius of the local-relief window in metres. | This is a physical radius in the raster CRS. The default window shape is circular. |
| `local_relief_window_shape` | `string` | `circular` | `param` | `write_local_relief_raster` | Shape of the local-relief window. | Valid options are circular or square. Default is circular. |

### Polynomial surface fit

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_polyfit_rasters` | `bool` | `false` | `switch` | - | Write the six quadratic polynomial fit coefficient rasters. | Computes a, b, c, d, e, f where z = ax^2 + by^2 + cxy + dx + ey + f over a moving window controlled by polyfit_lambda. |
| `polyfit_lambda_metres` | `float` | `70.0` | `param` | - | Characteristic lengthscale (metres) for the quadratic polynomial fit. | Controls the size of the moving window. Default (70 m) is tuned for Copernicus 30 m data. Reduce for lidar (e.g. 10-30 m). |
| `polyfit_metrics` | `string` | - | `param` | - | Comma-separated list of derived polyfit metrics to write. | Examples: slope,aspect_degrees,plan_curvature,profile_curvature. Use all to write every supported metric. Available metrics: fitted_surface, dzdx, dzdy, d2zdx2, d2zdxdy, d2zdy2, slope, slope_degrees,aspect_radians, aspect_degrees, laplacian, hessian_determinant, mean_curvature, gaussian_curvature, maximum_curvature, minimum_curvature, curvedness, shape_index, profile_curvature, plan_curvature, tangential_curvature, difference_curvature, horizontal_excess_curvature, vertical_excess_curvature, accumulation_curvature, ring_curvature, rotor. |
| `polyfit_z_factor` | `float` | `1.0` | `param` | - | Vertical scaling applied before interpreting polyfit derivatives. | Scales first and second derivatives consistently in the polyfit metric engine. |
| `polyfit_flat_slope_tolerance` | `float` | `1e-12` | `param` | - | Slope tolerance below which direction-dependent polyfit metrics are nodata. | Applies to aspect, plan curvature, tangential curvature, profile curvature, and rotor. |
| `polyfit_method` | `string` | `direct` | `param` | - | Algorithm used for the quadratic polynomial fit. | direct: separable rolling-moment O(n) fit (no FFTW needed). Default and preferred: faster and better at resolving fine-scale features. wavelet: FFT-based Gaussian wavelet fit (requires FFTW). both: run both methods and write outputs from each; also prints a timing comparison. Outputs are prefixed _wavelet_polyfit_* or _direct_polyfit_* so results from both methods can coexist in the same directory. |
| `compute_backend` | `string` | `auto` | `param` | - | Preferred compute backend for kernels that have a GPU port. | auto: use the GPU when a device and GPU kernels are both available, otherwise CPU (default). cpu: force CPU. gpu: prefer the GPU, falling back to CPU per-function when no GPU kernel exists or no device is present. Currently only direct_polyfit has a GPU port; all other metrics run on the CPU regardless of this setting. Has no effect in a build without GPU support (-DLSDTT_GPU=OFF). |
| `write_slope_threshold_layers` | `bool` | `false` | `switch` | - | Writes layers for each slope threshold where the polyfit metric slope_degrees exceeds the threshold | Thresholds are to be supplied in a comma-separated string under the parameter slope_tresholds. |
| `clamp_slope_to_threshold` | `bool` | `true` | `param` | - | If true, the slope threshold layers will have the value of the threshold where the original slope_degrees layer exceeds the threshold. If false, slopes from slope_degrees layer are preserved where slopes exceed threshold. | In both cases, nodatavalue below threshold. Defaults to true. |
| `slope_thresholds` | `string` | `3,7,12` | `param` | - | The threshold values for which slope threshold layers are to be written. | Defaults to the categories used by the Environment Agency's ALERT tool. |

### Surface shape classification (SMAP, TopoCurve-style)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_shape_classification_raster` | `bool` | `false` | `switch` | - | Write the morphological surface-shape classification raster (SMAP). | Reproduces the differential-geometry shape classification of the TopoCurve approach. Each pixel is labelled by the signs of Gaussian and mean curvature derived from the polyfit principal curvatures: peak (3), antiformal saddle (2), antiform (1), plane (0), synform (-1), synformal saddle (-2), basin (-3), perfect saddle (-4). Requires the polynomial fit; maximum_curvature and minimum_curvature are added to the requested metrics automatically. Output: _&lt;method&gt;_shape_class.&lt;ext&gt;; write_output_format selects the extension. |
| `shape_curvature_method` | `string` | `polyfit` | `param` | - | Source of the principal curvatures used for the shape classification. | polyfit: principal curvatures from the polynomial fit (compute_polyfit_metrics maximum/minimum_curvature). Default. fundamental_form: the differential-geometry method ported from TopoCurve CurveCalc, computed from the fitted surface via the first and second fundamental forms (gradients of the surface and of its unit normal field). Both feed the same classifier. With fundamental_form the polyfit fitted_surface is used as the smoothed input (smoothing scale set by polyfit_lambda_metres), and the fundamental-form principal curvatures are written as _&lt;method&gt;_ff_maximum_curvature.&lt;ext&gt; / _&lt;method&gt;_ff_minimum_curvature.&lt;ext&gt;. |
| `shape_curvature_threshold` | `float` | `0.0` | `param` | - | Curvature threshold (kt) applied to each principal curvature before shape classification. | Principal curvatures with magnitude &lt;= kt are set to zero, so KG and KM collapse toward planar/cylindrical classes. Increase to remove low-amplitude, long-wavelength structure (e.g. a background basin) from the classification. Units: 1/m. Default 0 keeps every pixel in a non-planar class. |
| `write_shape_curvature_rasters` | `bool` | `false` | `switch` | - | Write the post-threshold Gaussian and mean curvature rasters used for the shape classification. | These differ from the raw gaussian_curvature/mean_curvature polyfit metrics because they reflect shape_curvature_threshold (kt). Outputs: _&lt;method&gt;_shape_gaussian_curvature.&lt;ext&gt;, _&lt;method&gt;_shape_mean_curvature.&lt;ext&gt;. Only active when write_shape_classification_raster = true. |

### MADk2 roughness

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_madk2_rasters` | `bool` | `false` | `switch` | - | Write MADk2 roughness rasters. | Writes omnidirectional roughness (OR), direction of maximum continuity (DMC), and roughness anisotropy (RA) following Trevisani et al. (2023). |
| `madk2_window_radius_metres` | `float` | `10.0` | `param` | `write_madk2_rasters` | Circular MADk2 window radius in metres. | This is the radius of the moving circular support window in projected map units. |
| `madk2_lag_pixels` | `int` | `1` | `param` | `write_madk2_rasters` | Lag in pixels for the MADk2 second-order directional differences. | Must be &gt;= 1. Typical values are 1 or 2. |

### Frangi filter

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_frangi_raster` | `bool` | `false` | `switch` | - | Write a Frangi filter response raster. | Enhances line-like terrain features over multiple scales using the Hessian of a Gaussian-smoothed DEM. |
| `write_frangi_scale_raster` | `bool` | `false` | `switch` | - | Write the scale-at-maximum-Frangi-response raster. | Stores the physical scale in metres at which the maximum Frangi response occurred for each pixel. |
| `frangi_sigma_min_metres` | `float` | `1.0` | `param` | - | Minimum Frangi scale in metres. | This is the smallest Gaussian smoothing scale tested by the Frangi filter. |
| `frangi_sigma_max_metres` | `float` | `5.0` | `param` | - | Maximum Frangi scale in metres. | This general terrain-metrics workflow searches broader lineaments up to 5 m; the drain-specific feature finder deliberately stops at 2.5 m. |
| `frangi_n_scales` | `int` | `4` | `param` | - | Number of Frangi scales. | Scales are sampled logarithmically between frangi_sigma_min_metres and frangi_sigma_max_metres. |
| `frangi_beta` | `float` | `0.5` | `param` | - | Frangi beta parameter. | Controls suppression of blob-like structure relative to line-like structure. The neutral 0.5 default supports general ridge/valley exploration; the drain finder uses 0.4 to prefer narrow lines more strongly. |
| `frangi_c` | `float` | `1.0` | `param` | - | Frangi c parameter. | Controls suppression of low-structureness responses. Phase 5 decision D7 intentionally changed the scientific default from 15 to 1, matching lsdtt-feature-finder; set 15 explicitly to reproduce the former response. |
| `frangi_target` | `string` | `ridges` | `param` | - | Feature polarity for the Frangi filter. | Valid options are ridges or valleys. |

### Otsu threshold segmentation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_otsu_raster` | `bool` | `false` | `switch` | - | Write a binary raster produced by Otsu threshold segmentation. | The DEM is normalised to 8-bit and Otsu's method selects the threshold that minimises intra-class greyscale variance. Output: 1.0 = foreground (above threshold), 0.0 = background. |

### Hough line extraction

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_hough_lines_raster` | `bool` | `false` | `switch` | - | Write a binary raster of detected straight line segments. | Uses the progressive probabilistic Hough transform (Matas et al. 2000) on a Canny edge image derived from the DEM. The algorithm and parameters are those of OpenCV's HoughLinesP; the random order in which edge points are visited uses lsdtt3's own fixed-seed generator, so results are reproducible on every platform but not identical to OpenCV's. Output: 1.0 = line pixel, 0.0 = other. |
| `hough_rho` | `float` | `1.0` | `param` | `write_hough_lines_raster` | Distance resolution of the Hough accumulator in pixels. | Smaller values give finer angular resolution but increase memory and compute cost. |
| `hough_theta_degrees` | `float` | `1.0` | `param` | `write_hough_lines_raster` | Angle resolution of the Hough accumulator in degrees. | Converted to radians internally. Typical range: 0.5-2.0 degrees. |
| `hough_threshold` | `int` | `50` | `param` | `write_hough_lines_raster` | Minimum number of accumulator votes required to accept a line. | Increase to suppress short or faint lines. Decrease to detect more lines in noisy data. |
| `hough_min_line_length` | `float` | `30.0` | `param` | `write_hough_lines_raster` | Minimum line segment length in pixels. | Segments shorter than this are discarded. Scale this relative to your pixel resolution. |
| `hough_max_line_gap` | `float` | `10.0` | `param` | `write_hough_lines_raster` | Maximum gap in pixels between collinear points to bridge into one segment. | Increase to connect dashed or broken linear features. |

### Shared Canny parameters (used by Hough; also available standalone)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `canny_threshold` | `float` | `50.0` | `param` | - | Lower Canny hysteresis threshold in 8-bit greyscale space [0,255]. | Used internally by the Hough line step. Pixels with gradient magnitude below this are definitely not edges. |
| `canny_ratio` | `float` | `2.0` | `param` | - | Ratio of upper to lower Canny threshold (high = canny_threshold * ratio). | Canny recommends 2.0-3.0. Increase to suppress weak edges. |
| `canny_gaussian_sigma` | `float` | `1.0` | `param` | - | Sigma for the Gaussian blur applied before Canny edge detection. | 0.0 = no blur. Increase (e.g. 1.5-3.0) for noisy or coarse DEMs. |

### Curvature skeleton (VEC + plan curvature relief)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_curvature_skeleton_rasters` | `bool` | `false` | `switch` | - | Extract slope-break skeletons directly from curvature rasters. | Computes vertical_excess_curvature and plan_curvature via direct polyfit (curvature_skeleton_lambda_metres), then: (1) thresholds VEC at vec_skeleton_threshold (or CDF divergence auto-threshold) and skeletonises; (2) computes local relief of plan_curvature over plan_relief_window_metres, thresholds at plan_relief_threshold, and skeletonises. Outputs use the selected raster extension: _vec_skeleton_curvature.&lt;ext&gt;, _vec_skeleton_binary.&lt;ext&gt;, _vec_skeleton.&lt;ext&gt;, _profile_curvature.&lt;ext&gt;, _profile_curvature_relief.&lt;ext&gt;, _profile_skeleton_binary.&lt;ext&gt;, _profile_skeleton.&lt;ext&gt;. |
| `curvature_skeleton_lambda_metres` | `float` | `100.0` | `param` | `write_curvature_skeleton_rasters` | Polyfit half-wavelength (m) for the direct fit used by curvature skeleton. | Controls the smoothing scale of the fitted surface. Larger values suppress noise; smaller values resolve finer features. |
| `vec_skeleton_threshold` | `float` | `0.0` | `param` | `write_curvature_skeleton_rasters` | VEC threshold for binarising vertical_excess_curvature. | Pixels with VEC &gt;= this value are classified as slope-break candidates. Run once with write_curvature_skeleton_rasters=true and read the printed CDF divergence threshold to choose an appropriate value. Ignored when vec_skeleton_use_cdf_threshold=true. |
| `vec_skeleton_use_cdf_threshold` | `bool` | `false` | `param` | `write_curvature_skeleton_rasters` | Auto-select the VEC threshold using CDF divergence. | Overrides vec_skeleton_threshold when true. |
| `vec_skeleton_write_cdf` | `bool` | `false` | `param` | `write_curvature_skeleton_rasters` | Write the VEC CDF to a CSV file alongside the skeleton outputs. | Output suffix: _vec_skeleton_cdf.csv. |
| `profile_relief_window_metres` | `float` | `50.0` | `param` | `write_curvature_skeleton_rasters` | Circular moving-window radius (m) for the profile curvature local relief. | Larger windows smooth sign-change boundaries but blur fine features. |
| `profile_relief_threshold` | `float` | `0.1` | `param` | `write_curvature_skeleton_rasters` | Threshold on profile curvature local relief for binarisation. | Pixels whose within-window profile curvature range &gt;= this value are classified as slope-break candidates. |
| `curvature_skeleton_spur_max_pixels` | `int` | `10` | `param` | `write_curvature_skeleton_rasters` | Maximum endpoint-rooted spur length (pixels) pruned from curvature skeletons. | 0 = skip spur pruning. |
| `curvature_skeleton_min_component_pixels` | `int` | `1` | `param` | `write_curvature_skeleton_rasters` | Minimum component size (pixels) retained in curvature skeletons. | Applied after thinning and spur pruning. |

### VEC adaptive thresholding

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vec_adaptive_threshold_rasters` | `bool` | `false` | `switch` | - | Apply adaptive mean and Gaussian thresholding to vertical_excess_curvature. | Uses curvature_skeleton_lambda_metres for the polyfit scale. Each pixel is compared to the mean (or Gaussian-weighted mean) of its local neighbourhood minus the constant vec_adaptive_C: a pixel is foreground if its normalised 8-bit value &gt; local_mean - vec_adaptive_C. Negative C is stricter (retains only locally brightest pixels). Outputs: _vec_adaptive_mean.&lt;ext&gt;, _vec_adaptive_gaussian.&lt;ext&gt;, and optionally their Zhang-Suen skeletons. |
| `vec_adaptive_block_size_pixels` | `int` | `51` | `param` | - | Neighbourhood block size (pixels) for adaptive VEC thresholding. | Must be odd and &gt;= 3. Should be larger than the typical feature width but smaller than the between-feature spacing. Larger values make the threshold surface smoother and less sensitive to local intensity variation. |
| `vec_adaptive_C` | `float` | `0.0` | `param` | - | Constant subtracted from the local neighbourhood mean before VEC thresholding. | Threshold = local_mean - C. 0 = threshold at the local mean (foreground = above-average pixels). Negative (e.g. -5) = stricter, retains only pixels well above the local mean. Positive = looser. |
| `vec_adaptive_gaussian_sigma` | `float` | `0.0` | `param` | - | Sigma for optional Gaussian pre-blur before adaptive VEC thresholding. | 0 = no blur. A small value (e.g. 1.0) smooths pixel-scale noise before the local threshold is applied. |
| `write_vec_adaptive_skeleton_rasters` | `bool` | `false` | `switch` | - | Thin the adaptive-threshold VEC binary masks to 1-pixel skeletons. | Applies Zhang-Suen thinning (and curvature_skeleton_spur_max_pixels pruning) to both the mean and Gaussian adaptive outputs. Outputs: _vec_adaptive_mean_skeleton.&lt;ext&gt;, _vec_adaptive_gaussian_skeleton.&lt;ext&gt;.Only active when write_vec_adaptive_threshold_rasters = true. |

### Cloth drape

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_cloth_envelope_raster` | `bool` | `false` | `switch` | - | Write the cloth-drape upper envelope surface. | Runs a Verlet cloth simulation: a flat sheet at max DEM elevation falls under gravity and comes to rest on ridges and summits. The resting position is the upper landscape envelope. Output: _cloth_envelope.&lt;ext&gt;. For best results supply a pit-filled or breached DEM. |
| `write_cloth_incision_raster` | `bool` | `false` | `switch` | - | Write the cloth-drape incision depth raster (envelope minus DEM). | Positive values indicate how far the terrain surface lies below the envelope at each pixel (valley incision depth). Negative values indicate residual pits or nodata artefacts above the envelope. Output: _cloth_incision.&lt;ext&gt;. |
| `write_cloth_crest_mask_raster` | `bool` | `false` | `switch` | - | Write the cloth contact (crest / summit) pixel mask. | 1 = cloth particle was pinned by terrain contact (ridge or summit); 0 = particle settled freely between ridges. Useful for diagnosing which features govern the envelope shape. Output: _cloth_crest_mask.&lt;ext&gt;. |
| `cloth_rigidity` | `float` | `1.0` | `param` | - | Spring stiffness of the cloth (controls envelope smoothness). | Higher values produce a smoother, longer-wavelength envelope that bridges across narrow valleys and ridges. Equivalent to the Tension parameter in the Carretier and Valenti EnvelopSurface script. Typical range: 0.5 (closely follows terrain) to 20 (broad smooth drape). |
| `cloth_gravity` | `float` | `0.3` | `param` | - | Downward gravitational force per Verlet step (m/step^2). | Controls how quickly the cloth descends. Reducing gravity gives the cloth more time to equilibrate spring forces before contact, which can improve convergence on complex terrain. Default (0.3) follows Zhang et al. (2016) CSF. |
| `cloth_class_threshold` | `float` | `0.5` | `param` | - | Contact tolerance in metres: cloth is pinned when within this distance of the DEM. | Increase (1-5 m) for coarse-resolution DEMs, noisy data, or DEMs with residual artefact pits. Decrease (0.1 m) for smooth high-resolution lidar where precise ridge contact is desired. |
| `cloth_max_iterations` | `int` | `500` | `param` | - | Maximum number of Verlet integration steps. | The simulation stops early when max_displacement &lt; cloth_convergence_tolerance_metres. Increase for large or complex DEMs where convergence is slow. |
| `cloth_convergence_tolerance_metres` | `float` | `0.005` | `param` | - | Convergence tolerance: maximum cloth displacement (m) per step to declare settled. | Tighten to 0.001 m for precise incision mapping on smooth terrain. Loosen to 0.01 m for faster runs on coarse DEMs where sub-centimetre precision is not needed. |

<!-- END GENERATED CLI PARAMETERS -->
