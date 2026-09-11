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
pathways from roads.

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
