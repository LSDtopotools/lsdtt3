# CLI reference: lsdtt-slope-area

`lsdtt-slope-area` calculates slope–area statistics from a DEM for one or more
direct-polynomial-fit smoothing lengths. It can route flow with D-infinity or
D8, optionally clip the analysis with a same-grid basin mask, write per-lambda
slope rasters, and write quantile statistics in logarithmic drainage-area
bins.

## Usage

```text
lsdtt-slope-area <path> <parameter_file>
lsdtt-slope-area <parameter_file>
lsdtt-slope-area -h
lsdtt-slope-area -v
lsdtt-slope-area -c
```

## Parameter file format

```text
# Lines beginning with # are comments.
key: value
```

## Pipeline overview

1. Parse the parameter file and load the DEM.
2. Optionally fill pits using the configured minimum slope and boundaries.
3. Build D-infinity or D8 flow information and drainage area.
4. Optionally write flow-information and drainage-area products.
5. For each requested polynomial-fit length, calculate slope.
6. Apply the optional basin mask and logarithmic drainage-area bins.
7. Write per-lambda statistics, slope rasters, and summary results.

## Input/output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | parameter-file directory | Directory containing the DEM. |
| `read_fname` | string | — | Input DEM filename. |
| `write_path` | string | `read_path` | Output directory. |
| `write_fname` | string | DEM stem | Output filename prefix. |
| `write_flowinfo_rasters` | bool | `true` | Write the drainage-area (and, for d8, flow-direction) rasters. |
| `write_slope_rasters` | bool | `true` | Write one slope raster per lambda. |
| `summary_csv_fname` | string | automatic | Summary CSV filename. |

## Routing and pit filling parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fill_before_routing` | bool | `true` | Fill pits before routing. |
| `fill_min_slope` | float | `0.0001` | Minimum slope used during pit filling. |
| `flow_algorithm` | string | `dinf` | `dinf` or `d8`. |
| `boundary_conditions` | string | `oooo` | Four-character boundary code, in `N,E,S,W` order. |

## Polynomial-fit and binning parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `polyfit_lambda_metres_list` | list | — | Comma-separated smoothing lengths in metres. |
| `n_log_bins` | int | — | Number of logarithmic drainage-area bins. |

## Basin clipping parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `basin_clip_mode` | string | `none` | `none` or `mask`. |
| `basin_mask_fname` | string | — | Same-grid raster mask when mode is `mask`. |

The basin mask must have exactly the same rows, columns, cell size, origin, and
coordinate reference system as the DEM. Nonzero values are inside the basin.
This CLI does not rasterize vector polygons; use `lsdtt-basin-extraction` first.

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of worker threads; `0` uses the implementation default. |

## Outputs

For each lambda, the tool writes:

```text
<write_fname>_lambda_XXm.csv
<write_fname>_lambda_XXm_slope.tif
```

The per-lambda CSV columns are:

```text
drainage_area,d05_slope,d25_slope,d50_slope,d75_slope,d95_slope,number_of_pixels
```

The summary CSV columns are:

```text
polyfit_lambda,drainage_area_with_max_d50_slope
```

Decimal lambda values use underscores in filenames, for example
`lambda_12_5m`.

When `write_slope_rasters: true`, an existing expected slope raster is opened
and reused instead of recalculating polynomial-fit slope. The raster is
validated against the DEM grid, after which CSV statistics are regenerated with
the current bin count and basin mask.

## Example parameter file

```text
read_path: /data/Moreno
read_fname: Moreno.tif
write_path: /results/Moreno/slope_area
write_fname: Moreno_slope_area
n_threads: 4
boundary_conditions: opop
fill_before_routing: true
fill_min_slope: 0.0001
flow_algorithm: dinf
write_flowinfo_rasters: true
polyfit_lambda_metres_list: 3,15,21
n_log_bins: 200
write_slope_rasters: true
basin_clip_mode: none
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

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_before_routing` | `bool` | `true` | `switch` | - | Choose whether to pit-fill the DEM before building the flow network. | ASKS: should the DEM be filled at all? Filling is OPTIONAL here and this is a genuine analytical choice. True (the default) routes over a pit-free surface. False builds the flow network on the UNFILLED surface, so flow paths terminate at sinks instead of being carried through them. Do not confuse this with raster_is_filled, which is an assertion that the input is already filled. |

### DEM preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_min_slope` | `float` | `0.0001` | `param` | `fill_before_routing` | Minimum slope imposed during pit filling. | Used when fill_before_routing is true. |
| `breach_before_fill` | `bool` | `true` | `switch` | - | Breach (carve) depressions before filling pits. | The carve comes first, then the fill (fill_before_routing) takes only the pits the carve could not solve, so large filled areas are avoided. false = fill alone, the routing surface of lsdtt-slope-area before 2026-11-18 (set it to reproduce those runs). |
| `breach_algorithm_choice` | `string` | `least_cost` | `param` | `breach_before_fill` | Algorithm for depression breaching. | lindsay2016: Lindsay (2016) efficient breaching. least_cost: Lindsay (2020) pit-centric least-cost breaching. |
| `breach_lc_max_dist` | `int` | `1000 (outlet) / 100 (whitebox)` | `param` | `breach_before_fill` | Maximum breach path length in cells for least-cost breaching. | Counted in steps from the pit (a diagonal step counts as one), not as a radius. Pits that need a longer path are left unbreached (and filled if fill_before_routing is true). 1 to 32767. If not set, 1000 with breach_lc_edge_policy outlet (long enough to bridge large water bodies, so a far-off pit is not carved to the nearest raster edge instead) and 100 with whitebox (the WhiteboxTools default; longer paths there can carve across divides). Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_max_cost` | `float` | `unlimited` | `param` | `breach_before_fill` | Maximum cost of a least-cost breach path. | The cost of a path is the depth each cell on it is lowered, summed: metres when breach_lc_minimize_dist is false; with breach_lc_minimize_dist true each step is weighted by its length in cells (1 or 1.414), so the units are metres x cells. A pit whose cheapest path costs more is left unbreached. Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_minimize_dist` | `bool` | `false` | `param` | `breach_before_fill` | Weight least-cost breach paths by their length. | false: a path costs the total depth lowered (m), so the search prefers shallow cuts however long. true: each step is weighted by its length in cells, favouring shorter paths; breach_lc_max_cost is then in metres x cells. Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_edge_policy` | `string` | `outlet` | `param` | `breach_before_fill` | Whether a least-cost breach may drain off the raster edge. | outlet: the raster edge is an outlet, like a nodata cell, so a pit near the edge can breach straight off the grid. whitebox: the edge is not an outlet (as in WhiteboxTools), so a breach must reach an edge cell that is already low enough, which can carve deep cuts across divides near the edge. Only used when breach_algorithm_choice is least_cost. |
| `breach_lc_flat_increment` | `float` | `0.0` | `param` | `breach_before_fill` | Elevation drop (m) per cell along a least-cost breach channel. | 0 or less selects the automatic value: twice the float32 spacing at the largest absolute elevation in the DEM (e.g. 2.4e-4 m for elevations of 1024-2047 m), the smallest drop that keeps the channel strictly downhill. Only used when breach_algorithm_choice is least_cost. |
| `exclude_carved_and_filled_pixels` | `bool` | `true` | `switch` | - | Leave carved and filled pixels out of the slope-area measurements. | A pixel is filled where the routing surface (after breach_before_fill and fill_before_routing) is above the input DEM, carved where it is below. Their gradient is the conditioning's, not the landscape's (a carved reach has the breach's own tiny slope), so they are left out of the bins and of the D8 flow-path fits. Drainage area is still routed through them. Polyfit slope rasters are written without this exclusion (so a saved raster does not depend on it); &lt;write_fname&gt;_conditioned_pixels.tif shows the excluded pixels. false keeps them, as runs before 2026-11-18 did (set it to reproduce those runs). |
| `write_conditioned_pixels_raster` | `bool` | `true` | `switch` | - | Write the conditioned-pixels raster. | 0 untouched, 1 filled (the routing surface is above the input DEM), 2 carved (below it). Output: &lt;write_fname&gt;_conditioned_pixels.tif. |

### Slope-area analysis

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `polyfit_on_routing_surface` | `bool` | `false` | `switch` | - | Fit the polyfit slope to the routing surface instead of the input DEM. | By default polyfit runs on the input DEM (after interior-nodata handling, before any carve or fill), since the conditioning is for routing only; its slope raster is &lt;write_fname&gt;_lambda_XXm_input_dem_slope.tif. true fits the carved / filled routing surface, as runs before 2026-11-18 did (set it to reproduce those runs); its raster keeps the old name &lt;write_fname&gt;_lambda_XXm_slope.tif. Each mode reuses only its own saved raster. |
| `slope_algorithm_choice` | `string` | `polyfit` | `param` | - | Gradient method: polyfit or d8_flow_path. | polyfit: direct polynomial fit over each lambda in polyfit_lambda_metres_list (outputs &lt;write_fname&gt;_lambda_XXm*). d8_flow_path: the least-squares gradient along the D8 flow path, over a window centred on each pixel of each length in d8_gradient_window_metres_list: downstream along the D8 receivers, upstream along the D8 donor of largest drainage area (drainage_area_algorithm_choice), on the routing surface; outputs &lt;write_fname&gt;_d8window_XXm* and the summary &lt;write_fname&gt;_maximum_d50_slope_area_d8window.csv. Carved and filled pixels are left out of the D8 fits when exclude_carved_and_filled_pixels is true. The valid cells must span half the window less one diagonal step. |
| `d8_gradient_window_metres_list` | `string` | `150` | `param` | `slope_algorithm_choice=d8_flow_path` | Comma-separated total lengths (m) of the D8 flow-path gradient window. | The window is centred on the pixel (half downstream, half upstream). Each must be at least four cell sizes (120 m on 30 m data). A pixel whose valid cells span less than half the window (e.g. near a divide and an outlet at once) is nodata. D8 gradient rasters are always recomputed, never read back. Example: 90,150,300. |
| `polyfit_lambda_metres_list` | `string` | `15` | `param` | - | Comma-separated direct-polyfit lambdas in metres. | Example: 3,5,7,15,30. Decimal filename labels use _ instead of . |
| `n_log_bins` | `int` | `20` | `param` | - | Number of logarithmic drainage-area bins. | Only non-empty bins are written to each lambda CSV. |
| `write_slope_rasters` | `bool` | `true` | `switch` | - | Write each slope raster, after optional clipping. | Files are named &lt;write_fname&gt;_lambda_XXm_input_dem_slope.tif (polyfit on the input DEM, the default) or &lt;write_fname&gt;_lambda_XXm_slope.tif (polyfit_on_routing_surface). An existing file of that name is reused. |

### Flow information

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `drainage_area_algorithm_choice` | `string` | `dinf` | `param` | - | Drainage-area method: dinf (D-infinity, Tarboton 1997) or d8. | Area is routed on the routing surface (after breach_before_fill and fill_before_routing). Independent of slope_algorithm_choice: the D8 flow-path gradient walks D8 receivers whichever area is chosen, and uses this area to pick the upstream path. Renamed from flow_algorithm on 2026-11-18 (same values; an old file must be edited). |
| `write_flowinfo_rasters` | `bool` | `true` | `switch` | - | Write the drainage-area and flow-direction rasters. | The drainage-area raster is used by the slope-area analysis. |
| `write_log10_drainage_area_raster` | `bool` | `false` | `switch` | - | Write the base-10 logarithm of the drainage-area raster. | Values are log10 of each positive drainage-area pixel; non-positive and nodata pixels remain nodata. |
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four characters in north, east, south, west order. |

### Threading

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 uses the runtime default. |

### Basin clipping

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_clip_mode` | `string` | `none` | `param` | - | Optional basin clipping mode. | Allowed values are none, mask, or polygon. |
| `basin_mask_fname` | `string` | - | `param` | `basin_clip_mode=mask` | Raster basin mask filename. | The mask must have exactly the DEM grid settings. Cells with value 0 or nodata are excluded; any other value is kept, so a labelled basins raster works. A relative path is read from read_path. |
| `basin_polygon_fname` | `string` | - | `param` | `basin_clip_mode=polygon` | Basin polygon filename. | The polygon must be supplied as a same-grid rasterised polygon mask; vector rasterisation is not performed by this CLI. It is read exactly like basin_mask_fname. A relative path is read from read_path. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `summary_csv_fname` | `string` | - | `param` | - | Summary CSV filename. | Empty uses &lt;write_fname&gt;_maximum_d50_slope_area.csv. |

### Low concavity

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_low_concavity_detection` | `bool` | `false` | `switch` | - | Find channel reaches whose gradient changes little with drainage area. | At each channel pixel (threshold_contributing_pixels), take it and the next low_concavity_window_pixels - 1 pixels DOWNSTREAM along the D8 receivers (downstream because the path upstream branches; the window runs through tributary junctions), and fit log10 S = a + b log10 A over those with a gradient that are not carved or filled (exclude_carved_and_filled_pixels). The concavity is theta = -b (S ~ A^-theta). A pixel is low concavity if \|theta\| &lt; low_concavity_max_abs_concavity, so a weakly convex reach counts too. S is always the centred D8 flow-path gradient over low_concavity_gradient_window_metres, whatever slope_algorithm_choice is; A is the drainage_area_algorithm_choice area. Aimed at threshold slopes, debris-flow channels and alluvial fans. Outputs: &lt;write_fname&gt;_low_concavity_theta.tif and _low_concavity_flag.tif (write_low_concavity_rasters) and _low_concavity_points.csv (write_low_concavity_point_table). |
| `low_concavity_gradient_window_metres` | `float` | `150` | `param` | `run_low_concavity_detection` | Total length (m) of the centred D8 flow-path gradient window used by the low-concavity fit. | As d8_gradient_window_metres_list, but one value, and used whatever slope_algorithm_choice is. At least four cell sizes. |
| `low_concavity_window_pixels` | `int` | `10` | `param` | `run_low_concavity_detection` | Pixels in the downstream fit window (the pixel and the next ones down). | Fewer at an outlet. At least 2. |
| `low_concavity_min_pixels` | `int` | `5` | `param` | `run_low_concavity_detection` | Fewest (area, gradient) pairs for a valid fit. | Carved, filled and gradient-less pixels in the window do not count. At least 3 and at most low_concavity_window_pixels. A pixel without a valid fit has flag -1. |
| `low_concavity_max_abs_concavity` | `float` | `0.15` | `param` | `run_low_concavity_detection` | Threshold on \|theta\| below which a channel pixel is low concavity. | Fluvial channels are typically 0.35-0.6. Must be &gt;= 0. |
| `write_low_concavity_rasters` | `bool` | `true` | `switch` | `run_low_concavity_detection` | Write the low-concavity theta and flag rasters. | &lt;write_fname&gt;_low_concavity_theta.tif (theta where the fit is valid) and _low_concavity_flag.tif (1 low concavity, 0 not, -1 no valid fit; nodata off the channel network or outside the basin clip). |
| `write_low_concavity_point_table` | `bool` | `false` | `switch` | `run_low_concavity_detection` | Write one CSV row per channel pixel with its low-concavity fit. | &lt;write_fname&gt;_low_concavity_points.csv: flat, x, y, drainage_area_m2, slope, theta, r2, n_pixels, log10_area_range, n_junctions_in_window, low_concavity. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Interior nodata

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `interior_nodata_method` | `string` | `sink` | `param` | - | How nodata pixels inside the DEM are treated: sink or fill. | Interior nodata are nodata pixels not connected to the raster edge through other nodata (isolated voids, moulins, masked lakes). sink (the default, unchanged behaviour): they are outlets, so flow that reaches them ends there, as at the DEM edge; use this for real internal drainage (e.g. moulins on an ice sheet). fill: they are interpolated with nodata_fill_method_choice before any filling or routing, so flow passes across them; use this for data voids. Nodata connected to the edge (e.g. sea removed by min_elevation) is never filled and always acts as an open boundary. |
| `nodata_fill_method_choice` | `string` | `laplacian` | `param` | `interior_nodata_method=fill` | Method for filling interior nodata pixels. | laplacian: smooth Dirichlet interpolation (recommended). iterative_weighted_average: fast, well-suited to thin voids. Used when interior_nodata_method is fill. |
| `nodata_fill_max_iter` | `int` | `10000` | `param` | `interior_nodata_method=fill` | Maximum iterations for iterative nodata fill methods. | Used when interior_nodata_method is fill. |

<!-- END GENERATED CLI PARAMETERS -->
