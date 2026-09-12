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
| `min_elevation` | float | `0.0` | Pixels at or below this elevation (m) are masked to nodata. The default removes sea and ocean pixels. |
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
| `basin_outlet_snapping_method` | string | `none` | Snap outlets to pixels exceeding a contributing-pixel threshold. Unlike extraction workflows, basin statistics defaults to `none` so an authoritative gauge/outlet pixel is not moved implicitly. `expanding_square` searches ring-by-ring; `downstream` follows the D8 receiver chain. |
| `basin_outlet_contributing_pixels_threshold` | int | `1000` | Minimum contributing-pixel count for snapping. |
| `basin_outlet_snap_min_search_radius_pixels` | int | `0` | Minimum search radius for `expanding_square` (pixels). `0` tests the input pixel first. |
| `basin_outlet_snap_max_search_radius_pixels` | int | `25` | Maximum search radius for `expanding_square` (pixels). |

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
