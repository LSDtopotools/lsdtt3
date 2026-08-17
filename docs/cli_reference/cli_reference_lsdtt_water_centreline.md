# CLI reference: lsdtt-water-centreline

Water/floodplain centreline extraction and tracing. Combines multiple methods to identify river centrelines and floodplain axes. Outputs centreline polylines and optional diagnostic rasters.

## Usage

```
lsdtt-water-centreline  <path> <parameter_file>
lsdtt-water-centreline  <parameter_file>        # path defaults to ./
lsdtt-water-centreline  -h                      # print help and exit
lsdtt-water-centreline  -v                      # print version and exit
lsdtt-water-centreline  -c                      # print citation and exit
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
| 2. Water mask detection | Identify water-related low-relief areas |
| 3. Skeletonisation | Extract 1-pixel centreline(s) |
| 4. Vectorisation | Convert centreline raster to polyline(s) |
| 5. Smoothing *(optional)* | Apply Douglas-Peucker or Visvalingam-Whyatt simplification |
| 6. Output | Write centreline polylines and diagnostics |

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

## Water mask detection

| Parameter | Type | Default | Description |
|---|---|---|---|
| `detection_method` | string | `relief` | Method to identify water areas: `relief` (local relief), `slope` (slope-based), `curvature` (curvature-based). |
| `relief_threshold_metres` | float | `5.0` | Relief threshold (metres) for water mask. Pixels with relief < threshold are water. |
| `slope_threshold_degrees` | float | `3.0` | Slope threshold (degrees). Pixels below this are potential water. |
| `smoothing_scale_metres` | float | `100.0` | Smoothing scale (metres) for relief/curvature computation. |

---

## Skeleton extraction

| Parameter | Type | Default | Description |
|---|---|---|---|
| `skeleton_method` | string | `zhang_suen` | Skeletonisation algorithm: `zhang_suen`, `thinning`. |
| `prune_skeleton_length_pixels` | int | `5` | Remove skeleton spurs shorter than this (pixels). |
| `bridge_small_gaps_pixels` | int | `3` | Bridge isolated skeleton segments separated by gaps < this (pixels). |

---

## Centreline smoothing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `apply_smoothing` | bool | `false` | Simplify centreline using Douglas-Peucker or similar. |
| `smoothing_tolerance_metres` | float | `10.0` | Tolerance (metres) for centreline simplification. Larger = simpler line. |
| `smoothing_algorithm` | string | `douglas_peucker` | Simplification algorithm: `douglas_peucker`, `visvalingam_whyatt`. |

---

## Output: Centreline vector

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_centreline_vector` | bool | `true` | Write centreline(s) as polyline vector layer. |
| `centreline_vector_format` | string | `flatgeobuf` | Output format: `flatgeobuf` / `fgb` (preferred), `gpkg`. |

Output file: `<write_fname>_water_centreline.<ext>`

---

## Output: Diagnostic rasters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_water_mask` | bool | `false` | Write binary water mask raster. |
| `write_skeleton_raster` | bool | `false` | Write skeleton/centreline as raster. |
| `write_relief_raster` | bool | `false` | Write relief raster used for detection. |

Output files:
- `<write_fname>_water_mask.tif`
- `<write_fname>_skeleton.tif`
- `<write_fname>_relief.tif`

---

## Output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `cog` | Raster output format (for diagnostics): `cog`, `gtiff`, `same_as_input`. |
