# CLI reference: lsdtt-slope-area

`lsdtt-slope-area` measures how channel and hillslope gradient change with
drainage area. It:

- conditions the DEM for routing (least-cost carve, then fill), and leaves the
  carved and filled pixels out of the measurements;
- routes drainage area with D-infinity or D8;
- computes gradient by a direct polynomial fit (one or more smoothing lengths)
  or along the D8 flow path (one or more window lengths);
- writes quantiles of gradient in logarithmic drainage-area bins;
- optionally finds **low-concavity reaches**: channel pixels where gradient
  changes little with area. It writes them as rasters, a point table and a
  segments vector layer.

The methods, and why the defaults are what they are, are in
[Slope-area analysis and low-concavity reaches](../design_and_algorithms/slope_area_and_low_concavity.md).
The figures are drawn by `lsdviz3-slope-area` in lsdviztools3. The same
routines are in Python (see [the Python guide](../python/slope_area.md)).

## Usage

```text
lsdtt-slope-area <path> <parameter_file>
lsdtt-slope-area <parameter_file>
lsdtt-slope-area -h
lsdtt-slope-area -v
lsdtt-slope-area -c
```

`-h` prints every parameter; `-h json` prints them as JSON.

## Parameter file format

```text
# Lines beginning with # are comments.
key: value
```

An unknown key stops the run. The removed key `flow_algorithm` stops it with a
message naming the new key, `drainage_area_algorithm_choice` (same values).

## Pipeline

1. Read the DEM and treat interior nodata (`interior_nodata_method`).
2. Build the routing surface:
   - carve depressions (`breach_before_fill`, default true;
     `breach_algorithm_choice`, default `least_cost`);
   - fill what is left (`fill_before_routing`, default true).
3. Compare the routing surface with the input DEM: a pixel above it is
   **filled**, below it **carved**. These pixels are left out of the bins and
   the fits (`exclude_carved_and_filled_pixels`, default true) and written to
   `_conditioned_pixels.tif`.
4. Route drainage area on the routing surface (`drainage_area_algorithm_choice`:
   `dinf` or `d8`).
5. Compute gradient (`slope_algorithm_choice`):
   - `polyfit`: a direct polynomial fit for each length in
     `polyfit_lambda_metres_list`, on the **input DEM** by default;
   - `d8_flow_path`: a least-squares line along the D8 flow path, centred on
     each pixel, for each length in `d8_gradient_window_metres_list`.
6. Apply the optional basin clip, bin by log drainage area, and write the
   quantiles and the summary.
7. With `run_low_concavity_detection: true`, fit log S against log A along a
   short window downstream of each channel pixel and flag the pixels with a
   small concavity.

## Outputs

`<write_fname>` is the prefix. `XX` is the lambda or window length in metres;
a decimal uses `_` (`lambda_12_5m`).

| File | When | Contents |
|---|---|---|
| `_conditioned_pixels.tif` | `write_conditioned_pixels_raster` (default) | 0 untouched, 1 filled, 2 carved |
| `_dinf_drainage_area.tif` or `_d8_drainage_area.tif`, `_d8_flow_direction.tif` | `write_flowinfo_rasters` (default) | drainage area (m²); the D8 receivers |
| `_lambda_XXm.csv` | polyfit | binned quantiles (below) |
| `_lambda_XXm_input_dem_slope.tif` | polyfit, `write_slope_rasters` | polyfit gradient of the input DEM (the default) |
| `_lambda_XXm_slope.tif` | polyfit, `polyfit_on_routing_surface: true` | polyfit gradient of the routing surface |
| `_maximum_d50_slope_area.csv` | polyfit | `polyfit_lambda,drainage_area_with_max_d50_slope` |
| `_d8window_XXm.csv` | `d8_flow_path` | binned quantiles |
| `_d8window_XXm_slope.tif` | `d8_flow_path` | D8 flow-path gradient (always recomputed) |
| `_maximum_d50_slope_area_d8window.csv` | `d8_flow_path` | as the polyfit summary, one row per window |
| `_low_concavity_theta.tif`, `_low_concavity_theta_se.tif`, `_low_concavity_flag.tif` | low concavity, `write_low_concavity_rasters` (default) | concavity of the window fit; its standard error; flag 1 low concavity, 0 not, -1 no valid fit |
| `_low_concavity_points.csv` | low concavity, `write_low_concavity_point_table` | one row per channel pixel (below) |
| `_low_concavity_segments.fgb` (or `.gpkg`) | low concavity, `write_low_concavity_segments` (default) | the channel network split into runs of one flag (below) |
| `_run_summary.json`, `_ingested_params.param` | always | run record; the parameters as read, re-runnable |

**Binned quantile CSVs** (`_lambda_XXm.csv`, `_d8window_XXm.csv`):

```text
drainage_area,d05_slope,d25_slope,d50_slope,d75_slope,d95_slope,number_of_pixels
```

`drainage_area` is the geometric centre of the bin (m²); only bins with
pixels are written. Gradient is dimensionless (m/m).

**Point table** (`_low_concavity_points.csv`):

```text
flat,x,y,drainage_area_m2,slope,theta,theta_se,r2,n_pixels,log10_area_range,n_junctions_in_window,low_concavity
```

- `slope` is the pixel's own D8 flow-path gradient.
- `theta`, `theta_se` and `r2` come from the fit over the window below the pixel.
- `n_pixels` is the number of pairs fitted.
- `log10_area_range` is the spread of log10 A in the window.
- `n_junctions_in_window` is the number of confluences the window crosses.

**Segments.** Each channel link (from a source or junction to the next junction
or outlet) is split where the flag changes. The attributes are `segment_id`, `link_id`,
`low_concavity`, `n_pixels`, `length_m`, `median_theta`, `median_theta_se`,
`median_slope`, `min_area_m2`, `max_area_m2`, `stream_order` (Strahler),
`upstream_flat` and `downstream_flat`. The segments partition the channel
network: every channel pixel is in exactly one.

**Reuse of polyfit rasters.** With `write_slope_rasters: true`, an existing
polyfit slope raster of the expected name is read back instead of being
recomputed, after checking that it is on the DEM's grid. The two polyfit modes
use different names, so a raster made in one mode is never reused in the
other. Delete the raster to force a recomputation.

## Choosing the methods

- **Gradient.**
  - `polyfit` measures the local surface gradient everywhere, hillslopes
    included. Use it for the hillslope-to-channel transition on a slope-area
    plot.
  - `d8_flow_path` measures gradient along the channel, averaged over the
    window, and is valid on fewer pixels (a pixel needs a long enough D8 path
    around it). Use it for channel slope-area plots.
- **Drainage area.**
  - `dinf` spreads flow over hillslopes and gives area that grows smoothly down
    a channel.
  - `d8` gives area that jumps at each confluence.

  For the low-concavity finder this matters; see the algorithms page.
- **Window lengths** (`d8_gradient_window_metres_list`,
  `low_concavity_gradient_window_metres`). These must be at least four cell
  sizes. The default of 150 m is five cells on 30 m data.

## Reproducing runs made before 2026-11-18

Runs before 2026-11-18 filled the DEM without carving, kept the filled pixels,
and fitted polyfit to the filled DEM. To get their lambda CSVs, slope rasters,
drainage area and summary CSV byte for byte, edit the old parameter file:

```text
# was: flow_algorithm: dinf
drainage_area_algorithm_choice: dinf
breach_before_fill: false
exclude_carved_and_filled_pixels: false
polyfit_on_routing_surface: true
```

A regression test pins these outputs by sha256.

## Example parameter files

Binned slope-area on a 30 m DEM with three polyfit lengths and the default
conditioning:

```text
read_path: /data/my_site
read_fname: my_site_DEM.tif
write_path: /data/my_site/slope_area
write_fname: my_site
polyfit_lambda_metres_list: 90,150,300
n_log_bins: 30
```

The D8 flow-path gradient with D8 area, two windows:

```text
read_fname: my_site_DEM.tif
write_fname: my_site
slope_algorithm_choice: d8_flow_path
d8_gradient_window_metres_list: 150,300
drainage_area_algorithm_choice: d8
```

Low-concavity reaches on channels from 100 contributing pixels, with the point
table:

```text
read_fname: my_site_DEM.tif
write_fname: my_site
threshold_contributing_pixels: 100
run_low_concavity_detection: true
write_low_concavity_point_table: true
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
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). CSV is available only through selectors for point/table outputs. Per-output overrides exist only where a container genuinely differs. Raster outputs are governed separately by write_output_format. |

### Low concavity

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_low_concavity_detection` | `bool` | `false` | `switch` | - | Find channel reaches whose gradient changes little with drainage area. | At each channel pixel (threshold_contributing_pixels), take it and the next low_concavity_window_pixels - 1 pixels DOWNSTREAM along the D8 receivers (downstream because the path upstream branches; the window runs through tributary junctions), and fit log10 S = a + b log10 A over those with a gradient that are not carved or filled (exclude_carved_and_filled_pixels). The concavity is theta = -b (S ~ A^-theta). A pixel is low concavity if \|theta\| &lt; low_concavity_max_abs_concavity, so a weakly convex reach counts too. S is always the centred D8 flow-path gradient over low_concavity_gradient_window_metres, whatever slope_algorithm_choice is; A is the drainage_area_algorithm_choice area. Aimed at threshold slopes, debris-flow channels and alluvial fans. Outputs: &lt;write_fname&gt;_low_concavity_theta.tif and _low_concavity_flag.tif (write_low_concavity_rasters) and _low_concavity_points.csv (write_low_concavity_point_table). |
| `low_concavity_gradient_window_metres` | `float` | `150` | `param` | `run_low_concavity_detection` | Total length (m) of the centred D8 flow-path gradient window used by the low-concavity fit. | As d8_gradient_window_metres_list, but one value, and used whatever slope_algorithm_choice is. At least four cell sizes. |
| `low_concavity_window_pixels` | `int` | `10` | `param` | `run_low_concavity_detection` | Pixels in the downstream fit window (the pixel and the next ones down). | Fewer at an outlet. At least 2. |
| `low_concavity_min_pixels` | `int` | `5` | `param` | `run_low_concavity_detection` | Fewest (area, gradient) pairs for a valid fit. | Carved, filled and gradient-less pixels in the window do not count. At least 3 and at most low_concavity_window_pixels. A pixel without a valid fit has flag -1. |
| `low_concavity_max_abs_concavity` | `float` | `0.15` | `param` | `run_low_concavity_detection` | Threshold on \|theta\| below which a channel pixel is low concavity. | Fluvial channels are typically 0.35-0.6. Must be &gt;= 0. |
| `write_low_concavity_rasters` | `bool` | `true` | `switch` | `run_low_concavity_detection` | Write the low-concavity theta and flag rasters. | &lt;write_fname&gt;_low_concavity_theta.tif (theta where the fit is valid), _low_concavity_theta_se.tif (its standard error) and _low_concavity_flag.tif (1 low concavity, 0 not, -1 no valid fit; nodata off the channel network or outside the basin clip). |
| `write_low_concavity_point_table` | `bool` | `false` | `switch` | `run_low_concavity_detection` | Write one CSV row per channel pixel with its low-concavity fit. | &lt;write_fname&gt;_low_concavity_points.csv: flat, x, y, drainage_area_m2, slope, theta, theta_se (standard error of theta), r2, n_pixels, log10_area_range, n_junctions_in_window, low_concavity. |
| `low_concavity_min_run_pixels` | `int` | `1` | `param` | `run_low_concavity_detection` | Shortest run of low-concavity pixels kept as a low-concavity segment. | A shorter flagged run (along one channel link) is relabelled 0 and merged with its neighbours in the segments file; the rasters and point table keep every pixel's own flag. 1 (default) keeps all runs. |
| `write_low_concavity_segments` | `bool` | `true` | `switch` | `run_low_concavity_detection` | Write the channel network as segments marked low concavity or not. | &lt;write_fname&gt;_low_concavity_segments (format from write_vector_format): each channel link (source or junction to the next junction or outlet) is split where the flag changes; attributes segment_id, link_id, low_concavity (1, 0, -1 no valid fit), n_pixels, length_m, median_theta, median_theta_se, median_slope, min_area_m2, max_area_m2, stream_order (Strahler), upstream_flat, downstream_flat. Each line ends at the first pixel below its last, so the lines join. |

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
