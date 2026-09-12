# CLI reference: lsdtt-swath-profile

Extracts swath profile statistics perpendicular to a baseline polyline. Computes min/max/mean elevation and other metrics within a corridor along the polyline, outputting profiles as CSV and optional rasters.

## Usage

```
lsdtt-swath-profile  <path> <parameter_file>
lsdtt-swath-profile  <parameter_file>        # path defaults to ./
lsdtt-swath-profile  -h                      # print help and exit
lsdtt-swath-profile  -v                      # print version and exit
lsdtt-swath-profile  -c                      # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

Keys use `snake_case`. All keys and their defaults are listed below.

---

!!! note "Multipart lines"
    `baseline_vector_fname` accepts multipart geometry (`MultiLineString` and friends). Each feature is split into its parts on read, with the attributes copied onto each, and the split is reported in the log.

    Shapefiles are the commonest source of multipart data, and **versions before August 2026 discarded it silently** — see [Vector data](../what_is_lsdtt3/vector_data.md).

## Pipeline overview

Steps run in this order:

| Step | Description |
|---|---|
| 1. DEM preparation | Load DEM, verify projection |
| 2. Baseline input | Read polyline from vector file or compute from points |
| 3. Corridor creation | Create swath perpendicular to baseline |
| 4. Sampling | Extract elevation and metrics within swath |
| 5. Statistics | Compute min, max, mean, percentiles, count |
| 6. Output | Write swath CSV and optional raster |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input DEM filename. Any GDAL-supported format. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |

---

## Baseline specification

| Parameter | Type | Default | Description |
|---|---|---|---|
| `baseline_type` | string | `vector` | Source of baseline: `vector` (polyline file), `points` (endpoints, auto-connect). |
| `baseline_vector_fname` | string | _(empty)_ | Polyline vector file (GeoPackage, FlatGeobuf, Shapefile). Required if `baseline_type: vector`. |
| `baseline_vector_format` | string | `gpkg` | Vector format: `gpkg`, `flatgeobuf`, `shapefile`. |
| `baseline_start_latitude` | float | _(empty)_ | Latitude of baseline start point. Used if `baseline_type: points`. |
| `baseline_start_longitude` | float | _(empty)_ | Longitude of baseline start point. |
| `baseline_end_latitude` | float | _(empty)_ | Latitude of baseline end point. |
| `baseline_end_longitude` | float | _(empty)_ | Longitude of baseline end point. |

---

## Swath parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `swath_width_metres` | float | `1000.0` | Total swath width (metres). Half extends on each side of baseline. |
| `swath_spacing_metres` | float | `100.0` | Spacing (metres) between sample points along baseline. |

---

## Statistics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `compute_min` | bool | `true` | Compute minimum elevation in swath at each position. |
| `compute_max` | bool | `true` | Compute maximum elevation. |
| `compute_mean` | bool | `true` | Compute mean elevation. |
| `compute_median` | bool | `false` | Compute median elevation. |
| `compute_std_deviation` | bool | `false` | Compute standard deviation. |
| `compute_percentiles` | bool | `false` | Compute percentile values (25th, 75th). |
| `compute_sample_count` | bool | `true` | Count valid pixels in swath. |

---

## Output: Swath CSV

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_swath_csv` | bool | `true` | Write swath profile CSV. |

Output file: `<write_fname>_swath.csv`

Columns: `distance_m`, `elevation_m`, `min_elevation_m`, `max_elevation_m`, `mean_elevation_m`, `std_dev_m`, `sample_count`, `latitude`, `longitude`

---

## Output: Swath raster

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_swath_raster` | bool | `false` | Write swath corridor as binary raster mask. |

Output file: `<write_fname>_swath_mask.tif`

---

## Output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `cog` | Raster output format (for swath mask): `cog`, `gtiff`, `same_as_input`. |

<!-- BEGIN GENERATED CLI PARAMETERS -->
## Parameters

This section is generated from the executable's `-h json` registry. Do not edit it manually.

### File I/O

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `read_path` | `string` | `./` | `io` | - | Directory from which input files are read. | Defaults to the directory holding the parameter file. |
| `read_fname` | `string` | - | `io` | - | Input filename, with or without extension (e.g. my_dem.tif). | Read from read_path. lsdtt3 uses GDAL, so any GDAL-readable raster format works. If empty, defaults to the parameter filename stem. |
| `write_path` | `string` | - | `io` | - | Directory to which output files are written. | Defaults to read_path when not set. |
| `write_fname` | `string` | - | `io` | - | Output filename prefix (no path, no extension). | All outputs are named &lt;write_fname&gt;_&lt;descriptor&gt;.&lt;ext&gt;. Defaults to the stem of read_fname when not set. |

### Baseline / value input

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `baseline_vector_fname` | `string` | - | `param` | - | Vector file (FlatGeobuf/GeoPackage/Shapefile) with the swath baseline. | The first LineString feature is used as the ordered baseline polyline. Coordinates must be in the raster's projected CRS (metres). |
| `value_raster_fname` | `string` | `NULL` | `param` | - | Optional raster whose values are binned along the swath. | If NULL (default), the input DEM's own values (elevation) are binned. Must share the input DEM's geometry. |

### Swath geometry

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `half_width_metres` | `float` | `0.0` | `param` | - | Only pixels within this distance of the baseline are included. | Set to 0 to include all valid pixels (full-extent swath). |
| `bin_width_metres` | `float` | `1000.0` | `param` | - | Along-axis bin width in metres. Must be &gt; 0. | - |
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads (0 = all available). | - |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_distance_rasters` | `bool` | `true` | `switch` | - | Write the along-axis, perpendicular, and signed-perpendicular rasters. | Output: &lt;prefix&gt;_swath_along_axis.tif, _swath_perp_distance.tif, _swath_signed_perp.tif. |
| `write_output_format` | `string` | `cog` | `io` | - | Raster output format for the distance rasters. | Allowed values: cog, gtiff, same_as_input. |

<!-- END GENERATED CLI PARAMETERS -->
