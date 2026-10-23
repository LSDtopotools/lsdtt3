# CLI reference: lsdtt-basin-surface-metrics

Extracts basins from a CSV of outlet points and measures multiscale surface
metrics (polynomial-fit slope, maximum curvature and others) along each basin
outline and inside each basin. For each smoothing length in `lambda_scales`
it can report how often curvature changes sign along the divide, log-binned
slope–area distributions, and the drainage area at which the median slope
peaks.

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

Keys use `snake_case`. All keys and their defaults are listed under
[Parameters](#parameters) below.

---

## Pipeline

1. Load the DEM, fill pits (skipped with `raster_is_filled: true`) and build
   D8 routing.
2. Read the outlets from `basin_outlet_fname` (a CSV with `latitude` and
   `longitude` columns, WGS84), snap them with `basin_outlet_snapping_method`
   (default `downstream`; also `expanding_square`, `largest_area` and `none`)
   and extract each basin.
3. Trace each basin's outline.
4. If a slope–area output is requested, build D-infinity routing on the filled
   DEM for the drainage area.
5. For each smoothing length λ in `lambda_scales` (default `10,50,100` m), fit
   the polynomial surface to the unfilled DEM and compute the metrics in
   `polyfit_metrics` (default `slope,maximum_curvature`).
6. Write the requested outputs.

## Example

```
read_fname: my_dem.tif
basin_outlet_fname: outlets.csv
lambda_scales: 10,50,100
polyfit_metrics: slope,maximum_curvature
write_basin_slope_area_stats: true
write_outline_curvature_flips: true
write_outline_polygons_with_flipstats: true
```

## Outputs

Every output is off by default. `<outlet_flat>` is the flat index of a
basin's snapped outlet pixel, and `<metric>_lambda<λ>` is, for example,
`slope_lambda50`.

| File | Written when | Contents |
|---|---|---|
| `<write_fname>_basin_<outlet_flat>_outline.csv` | `write_basin_outlines_csv` | One file per basin: distance along the outline, elevation and each metric at each λ |
| `<write_fname>_outline_curvature_flips.csv` | `write_outline_curvature_flips` (needs `maximum_curvature`) | Per basin and λ, how often and how far apart maximum curvature changes sign along the outline |
| `<write_fname>_basin_<outlet_flat>_slope_area_histogram.csv` | `write_basin_slope_area_csv` | One file per basin: slope percentiles (p05, p25, p50, p75, p95) in bins of log10 area, `slope_area_log_bin_width` wide |
| `<write_fname>_slope_area_stats.csv` | `write_basin_slope_area_stats` (needs `slope`) | Per basin and λ, the drainage area (m²) at which the median slope is greatest |
| `<write_fname>_basin_polygons.fgb` (or `.gpkg`) | `write_outline_polygons_with_flipstats` | Basin polygons with the curvature-flip and slope–area statistics for each λ as attributes. Turning this on also computes those statistics. |
| `<write_fname>_<metric>_lambda<λ>.tif` | `write_polyfit_rasters` | Each metric raster at each λ |
| `<write_fname>_dinf_area.tif`, `<write_fname>_dinf_log10area.tif` | `write_dinf_area_raster`, `write_dinf_area_log10_raster` | D-infinity drainage area (m²) and its log10 |

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

### Threading &amp; I/O

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all hardware threads. 1 = serial (default). |
| `write_output_format` | `string` | `cog` | `io` | - | Output raster format. | cog: Cloud Optimised GeoTIFF (default). gtiff: plain compressed. same_as_input: match input format. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |
| `fill_flat_resolution_choice` | `string` | `inline_gradient` | `param` | - | How flats are made drainable: inline_gradient, barnes2014b_gradient or barnes2014b_float_step. | inline_gradient (default, unchanged): the Priority-Flood fill imposes fill_min_slope during the flood. The barnes2014b options make flow converge across flats (Barnes et al. 2014, doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses fill_min_slope (m/m, must be &gt; 0) as the minimum gradient; flats are raised about 3-4x more, capped below their rims, and capped flats are reported. barnes2014b_float_step is equivalent to RichDEM ResolveFlats on float32 in its first pass: it ignores fill_min_slope and its slopes depend on elevation (of order ULP(z) / cell size). With either barnes2014b option the result is always routable: cells the first pass leaves without a lower neighbour are refilled and resolved again (up to 4 rounds), and any still left are raised by a final priority-flood pass by one float32 step each; the report gives the rounds and the cells changed after the first pass. |

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
| `basin_outlet_snapping_method` | `string` | `downstream` | `param` | - | Method to snap CSV points to the routing object. | none: uses mapped pixel. downstream: follows flow path. expanding_square: searches outward until channel hit. largest_area: the pixel with the most contributing pixels within basin_outlet_snap_max_search_radius_pixels. |
| `basin_outlet_contributing_pixels_threshold` | `int` | `1000` | `param` | - | Pixel threshold for confirming an outlet snap. | Defines the required D8 accumulation to arrest snapping. |
| `basin_outlet_snap_min_search_radius_pixels` | `int` | `0` | `param` | `basin_outlet_snapping_method` | Minimum expanding-square search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square. |
| `basin_outlet_snap_max_search_radius_pixels` | `int` | `25` | `param` | `basin_outlet_snapping_method` | Maximum expanding-square search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square or largest_area (the half-width of the square searched). |

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
