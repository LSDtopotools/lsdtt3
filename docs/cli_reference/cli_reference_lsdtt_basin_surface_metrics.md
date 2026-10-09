# CLI reference: lsdtt-basin-surface-metrics

Extracts predefined basins from outlet coordinates and computes multiscale surface metrics (slope, curvature, roughness, etc.) along basin outlines and within basin interiors. Outputs per-basin statistics and optional polygon layers with metrics as attributes.

## Usage

```
lsdtt-basin-surface-metrics  <path> <parameter_file>
lsdtt-basin-surface-metrics  <parameter_file>        # path defaults to ./
lsdtt-basin-surface-metrics  -h                      # print help and exit
lsdtt-basin-surface-metrics  -v                      # print version and exit
lsdtt-basin-surface-metrics  -c                      # print citation and exit
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
| 2. D8 flow routing | Build the D8 flow-info object |
| 3. Basin outlet parsing | Read outlet coordinates from CSV |
| 4. Basin extraction | Extract basin boundaries and outlines |
| 5. Metric computation | Compute surface metrics at multiple scales |
| 6. Outputs | Write basin statistics and optional polygons |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input DEM filename (with or without extension). Any GDAL-supported format. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `basin_outlet_fname` | string | _(empty)_ | Path to CSV file with basin outlets. Required columns: `latitude`, `longitude`. Optional: `gauge_id`, `gauge_name`. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` = all available hardware threads. |
| `boundary_conditions` | string | `oooo` | Four-character boundary condition string (N, E, S, W). `o` = open, `n` = no-flux, `p` = periodic. |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `raster_is_filled` | bool | `false` | Set `true` if DEM is already pit-filled. If `false`, `fill_pits()` is called automatically. |
| `fill_min_slope` | float | `0.0001` | Minimum slope (m/m) imposed during pit filling. Used only when `raster_is_filled` is `false`. |
| `threshold_contributing_pixels` | int | `1000` | Contributing-pixel threshold defining channel sources used to snap basin outlets. |

---

## Basin extraction

| Parameter | Type | Default | Description |
|---|---|---|---|
| `snap_outlet_to_channel` | bool | `true` | Snap outlet coordinates to the nearest D8 channel pixel before basin extraction. |
| `outlet_snapping_search_radius_cells` | int | `10` | Search radius (pixels) for outlet snapping. |

---

## Surface metrics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `lambda_scales` | vector of floats | `[10, 20, 50, 100]` | Smoothing scales (metres) for polynomial surface fitting. |
| `polyfit_method` | string | `direct` | Polynomial fitting method: `direct` (rolling-moment, preferred) or `wavelet`. |
| `write_slope_raster` | bool | `false` | Write slope (gradient magnitude) raster. |
| `write_curvature_raster` | bool | `false` | Write curvature raster. |
| `write_roughness_raster` | bool | `false` | Write surface roughness raster. |

---

## Output: Basin statistics CSV

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_basin_statistics_csv` | bool | `true` | Write per-basin hydromorphic statistics CSV. |
| `calculate_shape_metrics` | bool | `false` | Compute basin planform shape metrics (elongation, circularity, etc.). |
| `calculate_centre_of_mass_bearing` | bool | `false` | Compute centre-of-mass bearing from outlet. |

Output file: `<write_fname>_basin_stats.csv`

---

## Output: Basin outlines polygon

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_basin_polygons` | bool | `false` | Write basin outlines as polygon features. |
| `write_vector_format` | string | `flatgeobuf` | Basin polygon format: `flatgeobuf` / `fgb` (preferred), `gpkg`. |

Output file: `<write_fname>_basin_polygons.<ext>`

---

## Output: Basin outline CSV

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_basin_outline_csv` | bool | `false` | Write basin outline as CSV with elevation and metrics sampled along perimeter. |

Output file: `<write_fname>_basin_outline.csv`

---

## Output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `cog` | Raster output format. `cog`: Cloud-Optimised GeoTIFF (preferred). `gtiff`: GeoTIFF. `same_as_input`: same as input DEM format. |

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

### Threading &amp; I/O

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all hardware threads. 1 = serial (default). |
| `write_output_format` | `string` | `cog` | `io` | - | Output raster format. | cog: Cloud Optimised GeoTIFF (default). gtiff: plain compressed. same_as_input: match input format. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |
| `fill_flat_resolution_choice` | `string` | `inline_gradient` | `param` | - | How flats are made drainable: inline_gradient, barnes2014b_gradient or barnes2014b_float_step. | inline_gradient (default, unchanged): the Priority-Flood fill imposes fill_min_slope during the flood. The barnes2014b options make flow converge across flats (Barnes et al. 2014, doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses fill_min_slope (m/m, must be &gt; 0) as the minimum gradient; flats are raised about 3-4x more, capped below their rims, and capped flats are reported. barnes2014b_float_step is equivalent to RichDEM ResolveFlats on float32: it ignores fill_min_slope, its slopes depend on elevation (of order ULP(z) / cell size), and rim violations are reported. |

### Input &amp; Routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_outlet_fname` | `string` | `NULL` | `param` | - | Required file of requested basin outlet points (CSV format). | Must contain columns exactly named 'latitude' and 'longitude'. |
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four characters (north, east, south, west). o = open, n = no-flux, p = periodic, b = forced baselevel outlet. |
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope (m/m) imposed by pit filling or flat resolution. | inline_gradient (default): imposed during the Priority-Flood fill, so used only when raster_is_filled is false. barnes2014b_gradient: the minimum gradient across flats, also on a raster_is_filled input. barnes2014b_float_step: ignored. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Snapping options

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_outlet_snapping_method` | `string` | `downstream` | `param` | - | Method to snap CSV points to the routing object. | none: uses mapped pixel. downstream: follows flow path. expanding_square: searches outward until channel hit. |
| `basin_outlet_contributing_pixels_threshold` | `int` | `1000` | `param` | - | Pixel threshold for confirming an outlet snap. | Defines the required D8 accumulation to arrest snapping. |
| `basin_outlet_snap_min_search_radius_pixels` | `int` | `0` | `param` | `basin_outlet_snapping_method` | Minimum expanding-square search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square. |
| `basin_outlet_snap_max_search_radius_pixels` | `int` | `25` | `param` | `basin_outlet_snapping_method` | Maximum expanding-square search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square. |

### Smoothing Scales &amp; Polyfit

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `lambda_scales` | `string` | `10,50,100` | `param` | `write_polyfit_rasters` | List of lambda scales (metres) separating multi-scale sweeps. | These characteristic lengthscales dictate the moving window for the polynomial fit and represent roughness at varying scopes. |
| `polyfit_metrics` | `string` | `slope,maximum_curvature` | `param` | `write_polyfit_rasters` | List of derived polyfit metrics to extract. | Supported metrics: slope, slope_degrees, aspect_radians, maximum_curvature, minimum_curvature, mean_curvature, laplacian, plan_curvature, profile_curvature, etc. |
| `polyfit_method` | `string` | `direct` | `param` | `write_polyfit_rasters` | Algorithm for surface fitting. | direct: rolling-moment O(n) fit. wavelet: FFT Gaussian fit. |
| `write_polyfit_rasters` | `bool` | `false` | `switch` | - | Output all generated polyfit representations directly to disk. | - |

### Outline Outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_basin_outlines_csv` | `bool` | `false` | `switch` | - | Write 1D outline geometries and sampled surface metrics to CSVs. | Outputs one file per basin tracking exact values atop the divide. |
| `write_outline_curvature_flips` | `bool` | `false` | `switch` | - | Count spatial sign flips of maximum_curvature. | Requires 'maximum_curvature' polyfit metric. Outputs _outline_curvature_flips.csv charting frequency and spacing. |

### Slope-Area Outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_dinf_area_raster` | `bool` | `false` | `switch` | - | Dump spatial D-Infinity flow accumulation (in m2) as a GeoTIFF. | - |
| `write_dinf_area_log10_raster` | `bool` | `false` | `switch` | - | Dump spatial log10(D-Inf Area) to disk as a GeoTIFF. | - |
| `write_basin_slope_area_csv` | `bool` | `false` | `switch` | - | Export log-binned slope-area distributions per basin. | Writes percentiles of multi-scalar slope across log area bins. |
| `write_basin_slope_area_stats` | `bool` | `false` | `switch` | - | Export global summary locating the Area of maximum median slope. | Requires 'slope' polyfit metric. Bins using user specified width. |
| `slope_area_log_bin_width` | `float` | `0.1` | `param` | - | Resolution of log10(Area) binning. | Default 0.1 creates logarithmic intervals equal to 1/10th decade. |

### Vector Polygon Outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_outline_polygons_with_flipstats` | `bool` | `false` | `switch` | - | Write FGB/GPKG polygons populated with summary statistics. | If enabled, forces computation of slope_area_stats and curvature_flips, appending results per lambda as attributes. |
| `write_vector_format` | `string` | `flatgeobuf` | `param` | `write_outline_polygons_with_flipstats` | Vector output format for basin polygons. | Valid formats: flatgeobuf (or fgb) and gpkg. |

<!-- END GENERATED CLI PARAMETERS -->
