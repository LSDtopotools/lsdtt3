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
| 4. MADk2 roughness | `write_madk2_rasters: true` |
| 5. Frangi filter | `write_frangi_raster: true` or `write_frangi_scale_raster: true` |
| 6. Otsu threshold segmentation *(requires OpenCV)* | `write_otsu_raster: true` |
| 7. Hough line extraction *(requires OpenCV)* | `write_hough_lines_raster: true` |
| 8. Curvature skeletons | `write_curvature_skeletons: true` |
| 9. VEC adaptive thresholding *(requires OpenCV)* | `write_vec_adaptive_threshold: true` |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input raster filename. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `write_output_format` | string | `cog` | Output raster format. `cog` = Cloud Optimised GeoTIFF. `gtiff` = plain compressed GeoTIFF. `same_as_input` = match the input format. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `1` | Number of OpenMP threads. `0` uses all hardware threads. `1` = serial (default). Applies to the wavelet and direct polynomial fits, the polyfit metric engine, and OpenCV normalisation passes. |

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
| `polyfit_method` | string | `wavelet` | Algorithm for the fit. `wavelet`: FFT-based Gaussian wavelet fit (requires FFTW). `direct`: separable rolling-moment O(n) fit (no FFTW needed). `both`: run both algorithms and print a timing comparison. Output files are prefixed `_wavelet_polyfit_*` or `_direct_polyfit_*` so both sets coexist. |
| `polyfit_z_factor` | float | `1.0` | Vertical scaling applied before interpreting derivatives. Scales first and second derivatives consistently in the metric engine. |
| `polyfit_flat_slope_tolerance` | float | `1e-12` | Slope magnitude below which direction-dependent metrics (aspect, plan curvature, profile curvature, tangential curvature, rotor) are set to nodata. |

### Derived metrics

Set `polyfit_metrics` to a comma-separated list of metric names, or `all` to write every metric.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `polyfit_metrics` | string | _(empty)_ | Comma-separated list of derived metrics to compute and write. Empty = none. |

Available metric names and their output suffixes:

| Metric name | Output suffix |
|---|---|
| `fitted_surface` | `_<method>_polyfit_fitted_surface.tif` |
| `dzdx` | `_<method>_polyfit_dzdx.tif` |
| `dzdy` | `_<method>_polyfit_dzdy.tif` |
| `d2zdx2` | `_<method>_polyfit_d2zdx2.tif` |
| `d2zdxdy` | `_<method>_polyfit_d2zdxdy.tif` |
| `d2zdy2` | `_<method>_polyfit_d2zdy2.tif` |
| `slope` | `_<method>_polyfit_slope.tif` |
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
| `frangi_sigma_max_metres` | float | `5.0` | Maximum Gaussian smoothing scale in metres. |
| `frangi_n_scales` | int | `4` | Number of scales sampled logarithmically between min and max. |
| `frangi_beta` | float | `0.5` | Controls suppression of blob-like structure relative to line-like structure. |
| `frangi_c` | float | `15.0` | Controls suppression of low-structureness responses. |
| `frangi_target` | string | `ridges` | Feature polarity. `ridges` or `valleys`. |

---

## Step 6: Otsu threshold segmentation *(requires OpenCV)*

Output: `<write_fname>_otsu.tif`

Normalises the DEM to 8-bit and applies Otsu's method to select the threshold that minimises intra-class greyscale variance. Output is a binary raster: 1.0 = foreground (above threshold), 0.0 = background.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_otsu_raster` | bool | `false` | Enable Otsu segmentation. Silently skipped if OpenCV was not compiled in. |

---

## Step 7: Hough line extraction *(requires OpenCV)*

Output: `<write_fname>_hough_lines.tif`

Applies probabilistic Hough transform (HoughLinesP) to a Canny edge image derived from the DEM. Output is a binary raster marking detected line pixels.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_hough_lines_raster` | bool | `false` | Enable Hough line extraction. Silently skipped if OpenCV was not compiled in. |
| `hough_rho` | float | `1.0` | Distance resolution of the Hough accumulator in pixels. |
| `hough_theta_deg` | float | `1.0` | Angle resolution of the accumulator in degrees. |
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
| `write_curvature_skeletons` | bool | `false` | Enable curvature skeleton extraction. |
| `curvature_skeleton_lambda_metres` | float | `100.0` | Polyfit half-wavelength in metres for the direct fit. Larger values suppress noise; smaller values resolve finer features. |
| `vec_skeleton_threshold` | float | `0.0` | VEC value above which pixels are classified as slope-break candidates. Ignored when `vec_skeleton_use_cdf_threshold` is true. Run once and read the printed CDF divergence threshold to tune this. |
| `vec_skeleton_use_cdf_threshold` | bool | `false` | Auto-select the VEC threshold using CDF divergence. Overrides `vec_skeleton_threshold`. |
| `vec_skeleton_write_cdf` | bool | `false` | Write the VEC CDF to `_vec_skeleton_cdf.csv`. |
| `profile_relief_window_metres` | float | `50.0` | Circular window radius (m) for computing local relief of profile curvature. |
| `profile_relief_threshold` | float | `0.1` | Threshold on profile curvature local relief for binarisation. |
| `curvature_skeleton_spur_max_pixels` | int | `10` | Maximum endpoint-rooted spur length (pixels) pruned from skeletons after thinning. `0` = skip pruning. |
| `curvature_skeleton_min_component_pixels` | int | `1` | Minimum component size retained in skeletons after thinning and pruning. |

---

## Step 9: VEC adaptive thresholding *(requires OpenCV)*

Applies adaptive mean and Gaussian thresholding to the vertical excess curvature. Each pixel is compared to the local neighbourhood mean (or Gaussian-weighted mean) minus the constant `vec_adaptive_C`.

Outputs: `<write_fname>_vec_adaptive_curvature.tif`, `_vec_adaptive_mean.tif`, `_vec_adaptive_gaussian.tif`, and optionally their skeletons.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_vec_adaptive_threshold` | bool | `false` | Enable VEC adaptive thresholding. Silently skipped if OpenCV was not compiled in. Uses `curvature_skeleton_lambda_metres` for the polyfit scale. |
| `vec_adaptive_block_size_pixels` | int | `51` | Neighbourhood block size in pixels. Must be odd and ≥ 3. Should exceed the typical feature width but be smaller than the between-feature spacing. |
| `vec_adaptive_C` | float | `0.0` | Constant subtracted from the local mean before thresholding. `0` = threshold at the local mean. Negative (e.g. `-5`) = stricter; positive = looser. |
| `vec_adaptive_gaussian_sigma` | float | `0.0` | Sigma for optional Gaussian pre-blur before thresholding. `0` = no blur. |
| `write_vec_adaptive_skeletons` | bool | `false` | Thin the adaptive binary masks to 1-pixel-wide Zhang-Suen skeletons. Uses `curvature_skeleton_spur_max_pixels` and `curvature_skeleton_min_component_pixels`. Only active when `write_vec_adaptive_threshold` is true. |

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
polyfit_method:        wavelet
polyfit_metrics:       slope, plan_curvature, profile_curvature, vertical_excess_curvature

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
