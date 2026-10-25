# CLI reference: lsdtt-basin-stats

Computes hydromorphic statistics for basins defined by user-supplied outlet points. Requires a CSV file with outlet coordinates in WGS84; the program builds D8 flow routing, delineates each basin, and writes a summary CSV with elevation, slope, flow-distance, and optionally channel-length, planform-shape, and network-wide junction area-ratio statistics.

## Usage

```
lsdtt-basin-stats  <path> <parameter_file>
lsdtt-basin-stats  <parameter_file>        # path defaults to ./
lsdtt-basin-stats  -h                      # print help and exit
lsdtt-basin-stats  -v                      # print version and exit
lsdtt-basin-stats  -c                      # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Pipeline overview

| Step | Description |
|---|---|
| 1. Sea masking | Pixels at or below `min_elevation` are set to nodata |
| 2. Pit filling | Optional (`fill_before_routing`) |
| 3. D8 flow routing | Build the D8 flow-info object |
| 4. Outlet snapping | Optional (`basin_outlet_snapping_method`) |
| 5. Basin delineation | Extract upslope area for each outlet |
| 6. Statistics | Compute elevation, slope, flow-distance stats per basin |
| 7. Channel stats | Optional Strahler-order channel-length calculation |
| 8. Shape metrics | Optional per-basin planform shape metrics |
| 9. Junction area ratios | Optional network-wide junction drainage-area-ratio analysis |
| 10. Write outputs | `_basin_stats.csv`, optional per-basin PDF/CDF CSVs, and optional area-ratio CSVs |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input DEM filename. Any GDAL-supported format. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `write_output_format` | string | `cog` | Raster output format. `cog`, `gtiff`, or `same_as_input`. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` = all available hardware threads. |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `min_elevation` | float | `0.0` | Pixels at or below this elevation (m) are masked to nodata. The default removes sea and ocean pixels. Negative values apply as given; `-9999` disables. |
| `fill_before_routing` | bool | `true` | Fill DEM pits before building flow routing. Recommended. |
| `fill_min_slope` | float | `0.0001` | Minimum slope (m/m) imposed during pit filling. |

---

## Basin outlet input

The outlet CSV **must** contain columns named `latitude` and `longitude` (decimal degrees, EPSG:4326). Optional columns `gauge_id` and `gauge_name` are passed through to the output.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `basin_outlet_fname` | string | `NULL` | **Required.** Path to the outlet CSV file. |

---

## Outlet snapping

| Parameter | Type | Default | Description |
|---|---|---|---|
| `basin_outlet_snapping_method` | string | `none` | Snap outlets to pixels exceeding a contributing-pixel threshold. Unlike extraction workflows, basin statistics defaults to `none` so an authoritative gauge/outlet pixel is not moved implicitly. `expanding_square` searches ring-by-ring; `downstream` follows the D8 receiver chain; `largest_area` takes the pixel with the most contributing pixels in the square of half-width `basin_outlet_snap_max_search_radius_pixels`, if it reaches the threshold (a wide square can jump onto a bigger neighbouring river). |
| `basin_outlet_contributing_pixels_threshold` | int | `1000` | Minimum contributing-pixel count for snapping. |
| `basin_outlet_snap_min_search_radius_pixels` | int | `0` | Minimum search radius for `expanding_square` (pixels). `0` tests the input pixel first. |
| `basin_outlet_snap_max_search_radius_pixels` | int | `25` | Maximum search radius for `expanding_square`, and the square's half-width for `largest_area` (pixels). |

---

## Basin quality control

| Parameter | Type | Default | Description |
|---|---|---|---|
| `reject_nodata_edge_basins` | bool | `true` | Skip basins whose upslope area touches any nodata or raster-edge pixel. Prevents incorrect statistics for incomplete basins. |

---

## Hillshade (optional QC output)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_hillshade_raster` | bool | `false` | Write a hillshade raster after sea masking. |
| `hillshade_azimuth` | float | `315` | Sun azimuth (degrees clockwise from north). |
| `hillshade_altitude` | float | `45` | Sun altitude above horizon (degrees). |
| `hillshade_z_factor` | float | `1.0` | Vertical exaggeration factor. |

---

## Channel statistics

When enabled, builds a Strahler-order channel network and computes total channel length for orders 1–6.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `calculate_channel_stats` | bool | `true` | Build a channel/junction network and compute channel length by Strahler order. |
| `threshold_contributing_pixels` | int | `1000` | Contributing-pixel threshold defining channel sources for channel statistics. Also defines the network used for junction area ratios. |

---

## Basin shape metrics

Opt-in planform shape metrics for each basin. Off by default so the original curated hydromorphic stat set is reproduced unless explicitly requested. When enabled, the program traces each basin outline polygon (on the pixel-corner divide grid), converts it to projected coordinates, and derives a set of dimensionless shape indices. Nine columns are appended to `_basin_stats.csv`. See `docs/tex_technical/basin_shape_metrics.tex` for the full derivations, references, and bibliography.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `calculate_shape_metrics` | bool | `false` | Trace each basin outline polygon and compute planform shape metrics. When false the outline is not traced and the shape columns are omitted. |
| `basin_length_definition` | string | `flow_path` | Definition of basin length `L_b` used in the elongation ratio and form factor. `flow_path`: longest flow-path distance from the outlet (hydrologic basin length; the usual geomorphic convention). `geometric`: maximum Feret diameter of the basin outline polygon (longest straight-line caliper). The two are **not** interchangeable; both `Re` and `Rf` scale with `L_b`. |
| `shape_perimeter_simplify_cells` | float | `1.0` | Douglas–Peucker simplification tolerance for the perimeter, in cells. Applied as `(value × cell_size)` metres. Set to `0` to disable simplification and use the raw staircased perimeter. Area-based metrics are unaffected. |

### Metric definitions

`A` = basin polygon area (m²), `P` = perimeter (m, on the simplified ring), `L_b` = basin length (m), `A_hull` = convex-hull area (m²), `λ1 ≥ λ2` = eigenvalues of the area inertia tensor.

| Column | Definition / formula | Reference |
|---|---|---|
| `basin_length_m` | `L_b` under the chosen `basin_length_definition` | — |
| `elongation_ratio` | `Re = (2/√π)·√A / L_b = 1.128·√A / L_b` | Schumm (1956) |
| `form_factor` | `Rf = A / L_b²` | Horton (1932) |
| `circularity_ratio` | `Rc = 4π·A / P²` | Miller (1953) |
| `gravelius_compactness` | `Kc = P / (2√(πA)) = 0.2821·P / √A` | Gravelius (1914) |
| `solidity` | `A / A_hull` (≤ 1; 1 = convex) | — |
| `eccentricity` | `√(1 − λ2/λ1)` from the area inertia tensor | — |
| `orientation_deg` | COMPASS axial orientation (degrees clockwise from north, [0, 180)) of the principal axis from the area inertia tensor. Changed 2026-07-22 from the mathematical convention; older rasters in mathematical frame need converting via `compass_axial = (90 - maths_deg) mod 180`. | — |
| `hypsometric_integral` | `HI = (mean_elev − min_elev) / (max_elev − min_elev)` | Pike & Wilson (1971) |

Two algebraic identities act as internal cross-checks (they hold to rounding): `Kc²·Rc = 1` (Gravelius and Miller circularity carry the same `{P, A}` information) and `Re = (2/√π)·√Rf` (elongation and form factor carry the same `{A, L_b}` information).

**Perimeter staircase caveat.** The outline is traced on the raster pixel-corner grid, so its raw boundary is staircased: a 45° divide of true length `L` traces as `√2·L`. This systematically **overestimates** the perimeter and therefore biases the two perimeter-based metrics (`circularity_ratio`, `gravelius_compactness`) — the coastline paradox. To make these comparable across basins the perimeter is measured on a Douglas–Peucker-simplified ring with tolerance `shape_perimeter_simplify_cells × cell_size` metres. Area-based metrics (`elongation_ratio`, `form_factor`, `solidity`) use the area directly and are unaffected by the staircase.

---

## Junction drainage-area ratios

Opt-in, network-wide analysis that is a **continuous, high-sample-count complement to Tokunaga analysis**. At every interior confluence two (or more) channels merge; the program records the ratio `r = A_trib / A_main` of the two largest incoming drainage areas (in (0, 1]) at *every* junction — far more samples than the few well-populated Strahler-order bins. The distribution of `r` encodes the network's area-accumulation self-similarity. This requires a channel network: if `calculate_channel_stats` did not already build one, it is built here using `threshold_contributing_pixels`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `calculate_junction_area_ratios` | bool | `false` | Compute per-junction drainage-area ratios and write the three output CSVs below (with branch-bootstrap confidence intervals). |
| `junction_area_ratio_bootstrap` | int | `1000` | Branch-bootstrap replicates for the confidence intervals. Resampling is over non-mainstem tributary branches. `0` disables it. |

### Metric definitions

- **`r = A_trib / A_main`** — ratio of the smaller to the larger merging tributary drainage area at a confluence, in (0, 1]. Areas are D8 contributing-pixel counts.
- **`lambda` tail exponent** — `x = −ln(r)` is exponentially distributed if `r` has a power-law tail. `lambda_all` is the exponential rate fitted over all junctions; `lambda_tail` is fitted over the small-`r` tail (`x` above its median). A power-law tail in `r` is an exponential tail in `x`.
- **`R_A`** — Horton area ratio from a per-order regression of mean downstream area.
- **`c_equivalent`** — a **model-dependent** Tokunaga-equivalent `c ≈ R_A^lambda`. This is an estimate only, *not* a substitute for the order-binned Tokunaga `(a, c)`.
- **Branch bootstrap** — junctions are nested and hence autocorrelated, so the correct unit of independence is a tributary branch, not an individual junction. The main stem is always retained and resampling is over the non-mainstem branches, giving honest confidence intervals.
- **Area-decade collapse** — summary statistics restricted to confluences whose downstream area falls in each log10 area decade. If the network is self-similar these rows collapse (similar `median_ratio` / `lambda_tail`); this is the self-similarity diagnostic.

**Caveats.** This measures area accumulation, **not** side-branching topology, so it does not uniquely recover Tokunaga `(a, c)`. The metric is sensitive to the channel-extraction threshold (which truncates the small-`r` tail) and to the drainage-area definition (here D8 contributing pixels).

### Output files

| File | Contents | Columns |
|---|---|---|
| `<write_fname>_junction_area_ratios.csv` | One row per interior confluence | `junction_id`, `x`, `y`, `area_main_cells`, `area_trib_cells`, `ratio`, `area_down_cells`, `order_main`, `order_trib`, `order_receiver`, `n_tributaries`, `branch_id` |
| `<write_fname>_junction_area_ratio_summary.csv` | Single-row distribution summary with branch-bootstrap CIs | `n_junctions`, `mean_ratio`, `median_ratio`, `lambda_all`, `lambda_tail`, `r_a`, `c_equivalent`, `n_bootstrap`, `n_branches`, `median_ratio_ci_low`, `median_ratio_ci_high`, `lambda_tail_ci_low`, `lambda_tail_ci_high`, `r_a_ci_low`, `r_a_ci_high`, `c_equivalent_ci_low`, `c_equivalent_ci_high` |
| `<write_fname>_junction_area_ratio_decades.csv` | Area-decade collapse table (one row per populated log10 area decade) | `area_low_cells`, `area_high_cells`, `n`, `median_ratio`, `lambda_tail` |

---

## Per-basin elevation distribution

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_basin_pdf_cdf_csv` | bool | `false` | Write one elevation PDF/CDF CSV per basin. |
| `n_histogram_bins` | int | `100` | Number of histogram bins for the elevation PDF/CDF. |

Per-basin output: `<write_fname>_basin_<id>_pdf_cdf.csv`

Columns: `basin_id`, `bin_index`, `bin_lower_edge`, `bin_upper_edge`, `bin_centre`, `bin_width`, `count`, `probability`, `pdf`, `cdf`.

---

## Output

Main output: `<write_fname>_basin_stats.csv`

| Column | Description |
|---|---|
| `basin_id` | Sequential integer assigned by the program |
| `gauge_id`, `gauge_name` | Pass-through from the input CSV if present |
| `latitude`, `longitude` | Requested outlet coordinates (EPSG:4326) |
| `outlet_x`, `outlet_y` | Outlet in projected CRS |
| `outlet_row`, `outlet_col`, `outlet_flat` | Raster cell indices |
| `n_pixels` | Number of pixels in the basin |
| `drainage_area_km2` | Basin drainage area (km²) |
| `nodata_on_edge` | 1 if any upslope pixel touches nodata or raster edge |
| `centroid_x`, `centroid_y` | Area-weighted centroid in projected CRS |
| `centroid_latitude`, `centroid_longitude` | Centroid in WGS84 |
| `elevation_min`, `elevation_max`, `elevation_range` | Elevation statistics (m) |
| `elevation_mean`, `elevation_25th`, `elevation_50th`, `elevation_75th` | |
| `elevation_iqr`, `elevation_stddev`, `elevation_mad` | |
| `slope_25th`, `slope_50th`, `slope_75th` | Gradient-magnitude percentiles (m/m) |
| `flow_distance_25th`, `flow_distance_50th`, `flow_distance_75th` | Flow distance from outlet (m) |
| `basin_length_m`, `elongation_ratio`, `form_factor`, `circularity_ratio`, `gravelius_compactness`, `solidity`, `eccentricity`, `orientation_deg`, `hypsometric_integral` | Planform shape metrics (when `calculate_shape_metrics: true`) |
| `channel_length_O1` … `channel_length_O6` | Channel length (m) by Strahler order (when `calculate_channel_stats: true`) |

The shape-metric columns appear before the channel-length columns and are present only when `calculate_shape_metrics` is true.

---

## Minimal example parameter file

```
# lsdtt-basin-stats parameter file
read_path:         /data/dems/
read_fname:        catchment_dem.tif

basin_outlet_fname:                     /data/outlets.csv
basin_outlet_snapping_method:           downstream
basin_outlet_contributing_pixels_threshold: 2000

reject_nodata_edge_basins:  true
calculate_channel_stats:    true
threshold_contributing_pixels: 1000
```

## Example: enable basin shape metrics

```
calculate_shape_metrics:         true
basin_length_definition:         flow_path
shape_perimeter_simplify_cells:  1.0
```

## Example: enable junction drainage-area ratios

```
calculate_junction_area_ratios:  true
threshold_contributing_pixels:   1000
junction_area_ratio_bootstrap:   1000
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
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all hardware threads. 1 = single-threaded. |

### Output format

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_output_format` | `string` | `cog` | `io` | - | Output raster format. | cog: Cloud Optimised GeoTIFF (default). gtiff: plain compressed GeoTIFF. same_as_input: match the input raster format. |

### Basin outlets

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_outlet_fname` | `string` | `NULL` | `param` | - | CSV file of requested basin outlet points. | Required. CSV must contain columns named latitude and longitude. Optional columns: gauge_id, gauge_name. Coordinates are interpreted as WGS84 (EPSG:4326) decimal degrees. |

### Basin outlet snapping

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_outlet_snapping_method` | `string` | `none` | `param` | - | Optional snapping of basin outlet CSV points to pixels that exceed a contributing-pixel threshold. | Unlike extraction workflows, basin statistics defaults to none so an authoritative gauge/outlet pixel is not moved implicitly. Allowed values: none, expanding_square, downstream, largest_area. none: use the mapped raster pixel directly. expanding_square: search outward ring-by-ring until a qualifying pixel is found. downstream: follow the D8 receiver chain until a qualifying pixel is found. largest_area: take the pixel with the most contributing pixels in the square of half-width basin_outlet_snap_max_search_radius_pixels, if it is a qualifying pixel (a large radius can jump onto a bigger neighbouring river). |
| `basin_outlet_contributing_pixels_threshold` | `int` | `1000` | `param` | - | Contributing-pixel threshold used for optional basin outlet snapping. | This is the total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. |
| `basin_outlet_snap_min_search_radius_pixels` | `int` | `0` | `param` | - | Minimum expanding-square snapping search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square. 0 means the mapped input pixel is tested first. |
| `basin_outlet_snap_max_search_radius_pixels` | `int` | `25` | `param` | - | Maximum expanding-square snapping search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square or largest_area. For expanding_square the search stops at the first radius containing one or more qualifying pixels; for largest_area it is the half-width of the square searched. |

### Sea removal

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `min_elevation` | `float` | `0.0` | `param` | - | Pixels at or below this elevation are set to nodata. -9999 disables. | Default of 0.0 removes sea and ocean pixels from many DEM products. Negative values apply as given: for a basin below sea level (Death Valley -86 m, Dead Sea -430 m) set it just below the lowest land. -9999 (or lower) disables masking. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_before_routing` | `bool` | `true` | `switch` | - | Choose whether to pit-fill the DEM before building the flow network. | ASKS: should the DEM be filled at all? Filling is OPTIONAL here and this is a genuine analytical choice. True (the default) routes over a pit-free surface. False builds the flow network on the UNFILLED surface, so flow paths terminate at sinks instead of being carried through them. Do not confuse this with raster_is_filled, which is an assertion that the input is already filled. |
| `fill_flat_resolution_choice` | `string` | `inline_gradient` | `param` | - | How flats are made drainable: inline_gradient, barnes2014b_gradient or barnes2014b_float_step. | inline_gradient (default, unchanged): the Priority-Flood fill imposes fill_min_slope during the flood. The barnes2014b options make flow converge across flats (Barnes et al. 2014, doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses fill_min_slope (m/m, must be &gt; 0) as the minimum gradient; flats are raised about 3-4x more, capped below their rims, and capped flats are reported. barnes2014b_float_step is equivalent to RichDEM ResolveFlats on float32 in its first pass: it ignores fill_min_slope and its slopes depend on elevation (of order ULP(z) / cell size). With either barnes2014b option the result is always routable: cells the first pass leaves without a lower neighbour are refilled and resolved again (up to 4 rounds), and any still left are raised by a final priority-flood pass by one float32 step each; the report gives the rounds and the cells changed after the first pass. |

### Routing preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope (m/m) imposed during pit filling or flat resolution. | Used only when fill_before_routing is true. inline_gradient (default): imposed during the Priority-Flood fill. barnes2014b_gradient: the minimum gradient across flats. barnes2014b_float_step: ignored. |

### Basin quality control

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `reject_nodata_edge_basins` | `bool` | `true` | `param` | - | Reject basins whose upslope area touches any nodata or raster-edge pixel. | When true, basins that drain from the edge of the DEM are skipped. This prevents incorrect drainage area and statistics for incomplete basins. |

### Hillshade

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_hillshade_raster` | `bool` | `false` | `switch` | - | Write a hillshade raster after sea masking. | Useful for quick visual QC of the masked DEM. |
| `hillshade_azimuth` | `float` | `315` | `param` | `write_hillshade_raster` | Sun azimuth in degrees clockwise from grid north (grid up). | 315 is the standard cartographic northwest illumination. It is relative to the raster grid, so in an oblique frame it lights from the frame's north-west; bearings_relative_to_true_north never changes it. |
| `hillshade_altitude` | `float` | `45` | `param` | `write_hillshade_raster` | Sun altitude above horizon in degrees. | 45 degrees is a common default. |
| `hillshade_z_factor` | `float` | `1.0` | `param` | `write_hillshade_raster` | Vertical exaggeration factor for hillshade. | Values greater than 1 emphasise relief. |

### Channel statistics

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `calculate_channel_stats` | `bool` | `true` | `switch` | - | Build a channel/junction network and compute channel length by stream order. | Requires threshold_contributing_pixels to define channel sources. |
| `calculate_junction_area_ratios` | `bool` | `false` | `switch` | - | Compute per-junction drainage-area ratios (complement to Tokunaga). | Opt-in. Builds a channel network (using threshold_contributing_pixels) if one is not already built, then writes the ratio r = A_trib/A_main at every interior confluence plus a distribution summary. Outputs: &lt;write_fname&gt;_junction_area_ratios.csv, &lt;write_fname&gt;_junction_area_ratio_summary.csv (with branch-bootstrap confidence intervals), and a _junction_area_ratio_decades.csv collapse table. |
| `junction_area_ratio_bootstrap` | `int` | `1000` | `param` | - | Branch-bootstrap replicates for junction area-ratio confidence intervals. | Resampling is over non-mainstem tributary branches. 0 disables it. |
| `write_basin_polygons` | `bool` | `false` | `switch` | - | Write basin outline polygons with the per-basin metrics as attributes. | Attaches the hydromorphic stats (and shape metrics / centre-of-mass bearing when those are enabled) to each basin's outline polygon as a vector layer. Output: &lt;write_fname&gt;_basin_polygons.&lt;ext&gt;. Format set by write_vector_format. |
| `basin_stats_write_vector_format` | `string` | `flatgeobuf` | `param` | - | Format of the basin statistics table: flatgeobuf, gpkg or csv. | flatgeobuf (default) writes &lt;write_fname&gt;_basin_stats.fgb, one point per basin at its outlet in the DEM's CRS with every statistic as an attribute; gpkg the same as GeoPackage; csv writes &lt;write_fname&gt;_basin_stats.csv, where latitude/longitude and outlet_x/outlet_y locate each basin. Numbers are written at full double precision. |
| `write_vector_format` | `string` | `flatgeobuf` | `param` | - | Vector format for the basin-polygon output. | flatgeobuf (default) or gpkg. |
| `calculate_centre_of_mass_bearing` | `bool` | `false` | `switch` | - | Compute the bearing from each basin outlet to its centre of mass. | Opt-in. The basin centre of mass is the planimetric centroid of all upslope pixels (already reported as centroid_x/centroid_y). This adds com_bearing_deg (bearing clockwise from grid north, or from true north with bearings_relative_to_true_north, outlet -&gt; centre of mass) and com_distance_m columns. The outlet is the user-supplied (optionally snapped) outlet pixel; the channel is not followed. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Basin shape metrics

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `calculate_shape_metrics` | `bool` | `false` | `switch` | - | Trace each basin outline polygon and compute planform shape metrics. | Opt-in; off by default so the original curated hydromorphic stat set is reproduced unless explicitly requested. When true, adds basin_length_m, elongation_ratio, form_factor, circularity_ratio, gravelius_compactness, solidity, eccentricity, orientation_deg (principal axis, [0, 180) clockwise from grid north, or from true north with bearings_relative_to_true_north), and hypsometric_integral columns to the basin statistics table. When false, the outline is not traced and those columns are omitted, reproducing the original basin-stats output. See docs/tex_technical/basin_shape_metrics.tex. |
| `basin_length_definition` | `string` | `flow_path` | `param` | `calculate_shape_metrics` | Definition of basin length L_b used in elongation ratio and form factor. | flow_path: the longest flow-path distance from the outlet (hydrologic basin length; the default and the usual geomorphic convention). geometric: the maximum Feret diameter of the basin outline polygon (longest straight-line caliper). The two are not interchangeable and Re, Rf both scale with L_b. See docs/tex_technical/basin_shape_metrics.tex. |
| `shape_perimeter_simplify_cells` | `float` | `1.0` | `param` | `calculate_shape_metrics` | Douglas-Peucker simplification tolerance for the perimeter, in cells. | The basin outline is traced on the pixel-corner grid and is staircased, which inflates the perimeter and biases circularity_ratio and gravelius_compactness (the coastline paradox). The perimeter for those two metrics is measured on a ring simplified with tolerance (this value * cell_size) in metres. Set to 0 to disable simplification and use the raw staircased perimeter. Area-based metrics (elongation, form factor, solidity) are unaffected. |

### Basin distribution output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_basin_pdf_cdf_csv` | `bool` | `false` | `switch` | - | Write PDF/CDF CSV output for each basin. | If true, writes one CSV per basin containing histogram counts, probabilities, PDF, and CDF of basin elevations. |
| `n_histogram_bins` | `int` | `100` | `param` | `write_basin_pdf_cdf_csv` | Number of bins used for basin elevation histogram/PDF/CDF output. | Used when write_basin_pdf_cdf_csv is true. Default is 100 bins. |

### Interior nodata

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `interior_nodata_method` | `string` | `sink` | `param` | - | How nodata pixels inside the DEM are treated: sink or fill. | Interior nodata are nodata pixels not connected to the raster edge through other nodata (isolated voids, moulins, masked lakes). sink (the default, unchanged behaviour): they are outlets, so flow that reaches them ends there, as at the DEM edge; use this for real internal drainage (e.g. moulins on an ice sheet). fill: they are interpolated with nodata_fill_method_choice before any filling or routing, so flow passes across them; use this for data voids. Nodata connected to the edge (e.g. sea removed by min_elevation) is never filled and always acts as an open boundary. |
| `nodata_fill_method_choice` | `string` | `laplacian` | `param` | `interior_nodata_method=fill` | Method for filling interior nodata pixels. | laplacian: smooth Dirichlet interpolation (recommended). iterative_weighted_average: fast, well-suited to thin voids. Used when interior_nodata_method is fill. |
| `nodata_fill_max_iter` | `int` | `10000` | `param` | `interior_nodata_method=fill` | Maximum iterations for iterative nodata fill methods. | Used when interior_nodata_method is fill. |

### Orientations

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `bearings_relative_to_true_north` | `bool` | `false` | `param` | - | Write bearings, aspects and orientations clockwise from true north. | Off (default): they are clockwise from grid north (grid up), so in an oblique frame with the feature along +x, 90 degrees means along the feature. On: each value is converted with the grid convergence at its own location (each cell of a raster, each row of a table), which also removes the up to about 3 degree difference between UTM grid north and true north. Angles between two directions (junction angles, deflections, contrasts) are unchanged, and hillshade illumination always stays relative to grid up. Needs a raster with a CRS. |

<!-- END GENERATED CLI PARAMETERS -->
