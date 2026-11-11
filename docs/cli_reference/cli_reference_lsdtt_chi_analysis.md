# CLI reference: lsdtt-chi-analysis

End-to-end chi (χ) analysis of river networks. Extracts chi coordinates, fits piecewise-linear chi-elevation networks, computes disorder metrics with optional jackknife/bootstrap uncertainty, and performs spatial k_sn / u* inversion with m/n model selection.

**Requires Eigen3** for the k_sn inversion (`run_inversion`). Chi extraction and disorder sweep run without Eigen3.

## Usage

```
lsdtt-chi-analysis  <path> <parameter_file>
lsdtt-chi-analysis  <parameter_file>        # path defaults to ./
lsdtt-chi-analysis  -h                      # print help and exit
lsdtt-chi-analysis  -v                      # print version and exit
lsdtt-chi-analysis  -c                      # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Pipeline overview

Steps are gated by booleans and run in this order:

| Step | Trigger parameter |
|---|---|
| 1. Pit filling | `raster_is_filled: false` |
| 2. D8 flow routing | always |
| 3. Source extraction | `threshold_contributing_pixels` |
| 4. Chi extraction + data map | `write_chi_data_maps: true` |
| 5. Chi coordinate raster | `write_chi_coordinate_raster: true` |
| 6. Chi network fitting + disorder sweep | `fit_chi_networks: true` |
| 7. Disorder uncertainty | `run_subset_sensitivity`, `run_jackknife`, `run_bootstrap` |
| 8. Spatial k_sn / u* inversion | `run_inversion: true` (requires Eigen3) |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input DEM filename. Any GDAL-supported format. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `write_output_format` | string | `cog` | Raster output format for the chi coordinate and inversion rasters. `cog`, `gtiff`, or `same_as_input`. |

---

## Threading and routing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` = all available hardware threads. |
| `boundary_conditions` | string | `oooo` | Four-character boundary condition string (N, E, S, W). `o` = open, `n` = no-flux, `p` = periodic, `b` = forced outlet. |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `raster_is_filled` | bool | `false` | Set `true` if the DEM is already pit-filled. If `false`, `fill_pits()` is called before routing. |
| `fill_min_slope` | float | `0.0001` | Minimum slope (m/m) imposed during pit filling. |

---

## Basin outlet input

If `basin_outlet_fname` is `NULL` the largest basin in the DEM is used automatically.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `basin_outlet_fname` | string | `NULL` | CSV, FlatGeobuf, or GeoPackage of outlet points. CSV must have `latitude` and `longitude` columns (EPSG:4326). |
| `basin_outlet_snapping_method` | string | `downstream` | `none`, `expanding_square`, `downstream`, or `largest_area`. `largest_area` takes the pixel with the most contributing pixels in the square of half-width `basin_outlet_snap_max_search_radius_pixels`, if it reaches the threshold (a wide square can jump onto a bigger neighbouring river). |
| `basin_outlet_contributing_pixels_threshold` | int | `1000` | Contributing-pixel count required at the snapped outlet. |
| `basin_outlet_snap_min_search_radius_pixels` | int | `0` | Minimum search radius for `expanding_square` (pixels). |
| `basin_outlet_snap_max_search_radius_pixels` | int | `25` | Maximum search radius for `expanding_square`, and the square's half-width for `largest_area` (pixels). |

---

## Source extraction

| Parameter | Type | Default | Description |
|---|---|---|---|
| `threshold_contributing_pixels` | int | `1000` | Contributing-pixel count used to define channel heads. Chi is computed only on pixels upstream of channel heads. |
| `min_channel_pixels` | int | `20` | Reject basins whose extracted channel has fewer than this many nodes. `0` or `1` keeps all basins. |
| `n_pixels_move_outlet_upstream` | int | `10` | Move interior (non-edge) auto-detected outlets this many pixels upstream so the basin clears a coastal nodata fringe. `0` disables. |
| `n_largest_basins` | int | `5` | Keep only the N largest complete basins, ranked by basin pixels. `0` = no size cap. |

---

## Chi data maps (diagnostic output)

Runs before `fit_chi_networks`. Use to verify basin and network geometry before committing to a full m/n sweep.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_chi_data_maps` | bool | `false` | Extract chi networks and write per-node data maps. Outputs: nodes (chi, drainage_area, flow_distance, elevation, source_key, basin_key), valid basin outlets, source heads, and basin polygons. |
| `chi_data_maps_write_vector_format` | string | `csv` | Override for chi data map **point** files. `csv`: latitude/longitude in EPSG:4326. `gpkg`: GeoPackage. `flatgeobuf` / `fgb`: FlatGeobuf. Basin polygons are always vector (FlatGeobuf fallback when `csv` is requested). |
| `write_chi_data_map_lines` | bool | `false` | Write chi channels as polylines (one feature per `source_key`). Attributes: `source_key`, `basin_key`. |
| `write_vector_format` | string | `flatgeobuf` | Output format for chi channel and fitted-segment polylines. `flatgeobuf` / `fgb` or `gpkg`. CSV not supported for lines. |

---

## Chi coordinate raster

A raster of chi for **every pixel** of every retained complete basin, hillslopes
included, rather than for channel pixels only. This is the chi map used for
divide-migration and cross-divide chi-contrast analysis.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_chi_coordinate_raster` | bool | `false` | Write `<write_fname>_chi_coordinate_mn_*.tif` in `write_output_format`. |

Behaviour:

- Chi is integrated upstream from each retained basin outlet with chi = 0 at the outlet.
- Basins are selected exactly as for the other chi outputs: `threshold_contributing_pixels`,
  then `min_channel_pixels`, then `n_largest_basins`. The completeness filter
  **always** runs here, so a basin whose upslope area touches nodata or the DEM edge
  is never written even when `n_largest_basins: 0`.
- Every pixel outside the retained basins is nodata.
- One m/n value per run: the first value of `m_over_n`, or the first value of the
  `mn_min`/`mn_max`/`mn_step` sweep. Set `m_over_n` explicitly to control it.
- Nested outlets that disagree on baselevel produce inconsistent chi in a single
  raster; the step then warns and is skipped rather than aborting the run.

---

## Chi network fitting and m/n sweep

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fit_chi_networks` | bool | `false` | Run the piecewise-linear chi-network fitting and disorder sweep over the m/n range. |
| `chi_A0` | float | `1.0` | Reference drainage area A₀ (m²) for chi integration. χ = ∫(A₀/A)^(m/n) ds from outlet. **Leave at 1.0** unless you know why: the gradient of a chi-elevation profile equals `k_sn` only when `A0 = 1 m²`; any other value silently rescales it by `A0^(m/n)` (≈631× at `A0 = 1e6`, `m/n = 0.45`). |
| `m_over_n` | string | _(empty)_ | Optional explicit m/n value(s) for chi extraction and fitting. Comma-separated list (e.g. `0.45` or `0.35,0.45,0.55`). When set, it overrides `mn_min`, `mn_max`, and `mn_step`. |
| `mn_min` | float | `0.3` | Minimum m/n ratio for the disorder sweep. |
| `mn_max` | float | `0.6` | Maximum m/n ratio for the disorder sweep. |
| `mn_step` | float | `0.02` | Step size for the m/n disorder sweep. Must be > 0. |
| `fit_critical_divergence_metres` | float | `5.0` | Critical elevation divergence (m) used to split chi-elevation branches. |

---

## Disorder uncertainty

These options are only meaningful when `fit_chi_networks: true`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `run_subset_sensitivity` | bool | `false` | Run subset sensitivity analysis on chi disorder. |
| `subset_n_side_branches_per_sample` | int | `3` | Number of side branches selected per subset sensitivity sample. |
| `run_jackknife` | bool | `false` | Run leave-one-out jackknife analysis on chi disorder. |
| `run_bootstrap` | bool | `false` | Run bootstrap analysis on chi disorder. |
| `bootstrap_replicates` | int | `1000` | Number of bootstrap replicates. |
| `bootstrap_sample_full_tributary_count` | bool | `true` | Each replicate samples the full number of side branches. |
| `bootstrap_n_side_branches_per_replicate` | int | `3` | Number of side branches per bootstrap replicate (when not sampling full count). |
| `bootstrap_with_replacement` | bool | `true` | Bootstrap resampling with replacement. |
| `bootstrap_random_seed` | int | `5489` | Random seed for bootstrap resampling. |

---

## Spatial k_sn / u* inversion

Requires Eigen3. Performs Tikhonov-regularised spatial inversion of channel steepness over a coarse grid. Runs for each m/n value in the sweep; the m/n with minimum MAE or RMSE is selected as the best-fit concavity.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `run_inversion` | bool | `false` | Run spatial chi k_sn / u* inversion. |
| `inversion_coarse_cell_size` | float | `2000.0` | Cell size (m) of the coarse grid used in inversion. |
| `inversion_roughness_weight` | float | `10000.0` | Tikhonov regularisation weight for spatial smoothness. |
| `inversion_diagonal_jitter` | float | `1e-8` | Small diagonal perturbation added to the normal equations for numerical stability. |
| `inversion_use_mae_for_model_selection` | bool | `true` | Use MAE (`true`) or RMSE (`false`) for m/n model selection. |
| `write_inversion_rasters` | bool | `false` | Write u* and support-count rasters for each basin and each m/n value. |

---

## k_sn by lithologic unit

With `lithology_fname` set (and `fit_chi_networks: true`), the fitted channel
network is split by lithologic unit: each channel node is given the unit of the
DEM pixel it sits on, and the k_sn of the fitted segments is summarised per unit,
for every m/n value the run fits. Relative erodibility is not implemented; the
per-unit k_sn values are the measurement, and their ratios are left to the reader.
The model behind a k_sn contrast is described in
[Lithology and channel steepness statistics](../design_and_algorithms/lithology_channel_statistics.md).

### Inputs

- **`lithology_fname`** is a polygon vector (any OGR polygon layer: GeoPackage,
  FlatGeobuf, GeoJSON, Shapefile; reprojected to the DEM's CRS) or an integer
  code raster (`lithology_source`: `auto` decides from the extension). A raster
  not on the DEM grid is warped onto it with nearest neighbour.
- **`lithology_field_name`** (default `lith_primary`) names the attribute that
  holds the unit. Each distinct value is a unit unless a crosswalk groups them.
  Output from `lsdtt-fetch-lithology` carries `lith_primary` (the first
  lithology of each map unit).
- **`lithology_crosswalk_fname`** is an optional CSV with columns
  `value,code[,name]`. When given, every value present in the map must be listed
  (a missing value is a FATAL error that names the value and its pixel count).
  The `code` column sets the unit codes; the optional `name` column names each
  code, and without it a code is named by the map values joined with `; `. Two
  values with the same code are one unit.
- **Without a crosswalk**, units are numbered 1, 2, ... within the run, so the
  codes are only meaningful within one output set. Compare runs by unit name, not
  by code (`_lithology_units.csv` gives the legend).
- **Unmapped** pixels (a valid DEM pixel that no polygon touches, or nodata in a
  code raster) get code -1 and the name `unmapped`. They form their own group and
  appear as an `unmapped` row in every per-unit table.

### How nodes get their unit

Every channel node takes the unit of its own pixel. Each fitted segment (one piece of
the piecewise chi-elevation fit) also gets a length-weighted majority unit and a purity
(`lith_majority_fraction`, with `n_lith_units`), but segments are **not split at
contacts**. A node's k_sn is the k_sn of its segment, so a segment that crosses a
contact gives every node in it the segment's k_sn. This is why the per-unit
statistics use the node's own unit and the reaches described below, and why the
contact buffer exists.

### Distance to contact and the contact buffer

For each node the distance to the nearest mapped unit boundary is computed
(`distance_to_contact_m`). With `lithology_contact_buffer_metres` (default `0`),
nodes closer than that distance to a contact are set aside from the per-unit
statistics and flagged `near_contact = 1` in the node table
(`n_nodes_excluded_near_contact` counts them per unit). `0` keeps every node.
A buffer is useful when the map's positional error is comparable to a reach's
length: a node next to a contact may be in the wrong unit. The buffer also
shortens reaches at a contact: the end of a reach at a buffer is recorded as
`buffer`, not `contact`.

### Lithologic reaches

A **lithologic reach** is a run of consecutive channel nodes in one unit on one
branch. Reaches are the sampling unit of the per-unit k_sn statistics.

- **Construction.** Reaches are built along each branch of the fitted network
  (`branch_id`) and break whenever the unit changes or a buffer is entered.
  Each reach is a contiguous piece, so a unit that crosses the same channel
  several times gives several reaches.
- **OLS k_sn.** For each reach, the elevation of its nodes is regressed on chi:
  `k_sn_ols` is the slope, `k_sn_ols_se` its standard error and `r_squared` the
  coefficient of determination. The fit is ordinary least squares on the reach's
  own nodes. `k_sn_fit_chi_weighted` is a second estimate: the mean of the nodes'
  segment k_sn values, each weighted by the chi interval it spans.
- **End types.** `upstream_end` and `downstream_end` are one of `source` (the
  channel head), `junction` (a confluence), `outlet` (the network outlet),
  `contact` (the unit changes) or `buffer` (the contact buffer).
- **Qualifies.** A reach enters the per-unit statistics when it has at least
  `lithology_min_reach_pixels` nodes (default 10, minimum 2), spans a chi
  interval greater than 0, and has a finite slope. `qualifies` is written as 1 or 0
  in the reaches file.

### Statistics per unit and the bootstrap interval

Per unit (and per basin, in the basin table) the summary gives:

- **Node statistics** over every node not set aside by the buffer: node count,
  segment and basin counts, channel length, the median and quartiles of node
  k_sn, mean, minimum, maximum, and length-weighted mean and quartiles (weights
  are each node's flow length to its receiver).
- **Reach statistics** over the qualifying reaches: the number of reaches and of
  qualifying reaches, the median and quartiles of reach k_sn, the reach-length
  weighted median, and the median drainage area of the unit's nodes.
- **`low_sample`** is 1 when fewer than 5 reaches qualify. No bootstrap interval
  is computed for a low-sample unit (its interval columns are NaN).

The interval on the per-unit median reach k_sn is a percentile bootstrap: the
resampling draws reaches (or whole basins) with replacement, pools them, takes the
median, repeats `lithology_bootstrap_replicates` times (default 1000), and reports
the 2.5 and 97.5 percentiles of the medians.

- **`lithology_bootstrap_cluster`** `reach` (default) resamples lithologic reaches
  within a unit. `basin` resamples whole basins, so the reaches on one basin stay
  together. Reaches on one trunk channel are correlated, so the reach interval is
  a lower bound on the uncertainty; the basin interval is wider when there are few
  basins, and the per-basin table shows the spread between basins.
- **Seeding.** Each m/n value gets its own random stream, a `std::mt19937_64`
  seeded with `std::seed_seq{lithology_bootstrap_random_seed, round(m/n * 10000)}`.
  A fixed seed reproduces the intervals, and a value's interval does not depend
  on which other m/n values are in the same run.

### Outputs

Written once per run and per m/n value. The per-m/n names carry the value with
its `.` replaced by `p` (m/n 0.45 is `mn_0p450`). Each file is written only when
`write_ksn_lithology_summary` is true (the default with `lithology_fname`).

| File | Columns and contents |
|---|---|
| `<write_fname>_chi_ksn_by_lithology_mn_*.csv` | One row per unit, plus `unmapped` and `all`. `lith_code, lith_name, n_nodes, n_segments, n_basins, channel_length_m, k_sn_median, k_sn_q25, k_sn_q75, k_sn_mean, k_sn_min, k_sn_max, k_sn_length_weighted_mean, low_sample, m_over_n, A_0, contact_buffer_metres, min_reach_pixels, n_nodes_excluded_near_contact, n_reaches, n_reaches_qualifying, median_drainage_area_m2, reach_k_sn_median, reach_k_sn_p25, reach_k_sn_p75, reach_k_sn_length_weighted_median, reach_k_sn_median_ci_low, reach_k_sn_median_ci_high, n_bootstrap, bootstrap_cluster, node_k_sn_length_weighted_median, node_k_sn_length_weighted_p25, node_k_sn_length_weighted_p75`. The reach interval is in the `reach_k_sn_median_ci_*` columns. |
| `<write_fname>_chi_ksn_by_lithology_and_basin_mn_*.csv` | Reach statistics per unit and basin, with no interval. `lith_code, lith_name, basin_index, n_nodes, n_nodes_excluded_near_contact, channel_length_m, k_sn_median, n_reaches, n_reaches_qualifying, median_drainage_area_m2, reach_k_sn_median, reach_k_sn_p25, reach_k_sn_p75, reach_k_sn_length_weighted_median, low_sample, m_over_n, A_0, contact_buffer_metres, min_reach_pixels`. |
| `<write_fname>_chi_ksn_histogram_by_lithology_mn_*.csv` | Counts per unit in common k_sn bins. `lith_code, lith_name, bin_index, bin_lower, bin_upper, n_nodes, channel_length_m, m_over_n, A_0`. The bins span the 1st to 99th percentile of all nodes (`lithology_ksn_histogram_n_bins`, default 30, optionally log-spaced with `lithology_ksn_histogram_log_bins`). Nodes outside the range are in a row with `bin_index` -1 (below) or `n_bins` (above). |
| `<write_fname>_chi_ksn_lithology_nodes_mn_*.csv` | Every channel node, whether or not it was set aside. `basin_index, node_flat, row, col, x, y, chi, elevation, segment_id, k_sn, flow_length_to_receiver_m, lith_code, lith_name, distance_to_contact_m, latitude, longitude, drainage_area_m2, branch_id, near_contact, reach_id`. |
| `<write_fname>_chi_lithology_reaches_mn_*.csv` | Every lithologic reach. `reach_id, basin_index, branch_id, lith_code, lith_name, n_nodes, length_metres, chi_span, dz_metres, median_drainage_area_m2, k_sn_ols, k_sn_ols_se, r_squared, k_sn_fit_chi_weighted, upstream_flat, downstream_flat, upstream_end, downstream_end, qualifies`. |
| `<write_fname>_chi_lithology_segments_mn_*.csv` | Every fitted segment. `basin_index, segment_id, branch_id, k_sn, length_metres, lith_code, lith_name, lith_majority_fraction, n_lith_units, n_nodes`. |
| `<write_fname>_lithology_units.csv` | Once per run: the unit legend. `code, name, map_values, n_pixels, area_km2, n_pixels_in_analysed_basins`. |
| `<write_fname>_lithology_code.tif` | With `write_lithology_code_raster: true`: the unit code of every DEM pixel (-1 unmapped). |

### Interpreting per-unit k_sn

- **Uniform uplift and steady state.** A k_sn contrast between units reflects
  lithology only if uplift and base level are uniform and the landscape is near
  steady state. The CLI prints this as a reminder at each run.
- **Knickzones.** A transient knickzone puts a high or low k_sn on whichever unit
  it is crossing.
- **Strike-parallel sampling.** Drainage parallel to strike puts a unit mostly at
  channel heads or near water gaps, where k_sn is not comparable with mid-network
  reaches. `median_drainage_area_m2` per unit shows this; compare units at similar
  drainage areas before reading a contrast as lithological.
- **Alluvium** is a surface deposit, not bedrock. The k_sn of channels crossing
  it is not a measure of the bedrock beneath, and it often has few reaches (see
  the `low_sample` flag).
- **Map-scale positional error.** A contact drawn at 1:50 000 can be 50 m or more
  from the true contact. Use `lithology_contact_buffer_metres` to set aside the
  nodes that may be in the wrong unit, and compare runs with and without it.
- **Non-independence.** Nodes are not independent samples, and reaches on one
  channel are correlated. The intervals are a lower bound (see above).

### Coverage checks

- A map that labels **no** valid DEM pixel is refused (FATAL). The message gives
  the map's and the DEM's extents and CRSs, so a wrong CRS or an offset map is
  visible at once.
- A **partial** cover gives a note in the run output with the number and fraction of
  valid DEM pixels that have no unit. Above 10 % of the valid pixels a WARNING is
  also printed: the map covers part of the DEM, and the unlabelled pixels are in
  the `unmapped` group.

### Worked example: fetch, basins, chi

This uses the West Virginia Allegheny Front data (COP30 DEM, Macrostrat SGMC
units, `lith_primary`). The steps assume the Macrostrat fetch writes
`<write_fname>_macrostrat_units.fgb` and `<write_fname>_crosswalk_template.csv`,
as described in the [fetch lithology page](cli_reference_lsdtt_fetch_lithology.md).

```
# 1. DEM (keyless COP30)
write_fname: wvaf
dem_source: cop30_aws
west: -79.45
east: -79.15
south: 38.85
north: 39.05
target_epsg: 32617
grid_spacing: 30
```

Run `lsdtt-fetch-raster` on this, then `lsdtt-fetch-lithology` for the same box
(see its page). Fill in the crosswalk template with codes 10 (sandstone), 20
(siltstone), 30 (shale), 40 (limestone) and 90 (alluvium), and save it as
`wvks_crosswalk_lith_primary.csv` with a `name` column:

```
value,code,name
sandstone,10,sandstone
siltstone,20,siltstone
shale,30,shale
limestone,40,limestone
fine alluvium,90,alluvium
```

The basins come from `lsdtt-basin-extraction` (threshold 200 contributing pixels,
complete and non-nested basins of at least 1000 pixels, `write_selected_basin_outlets_csv: true`),
and the chi run takes its outlets from that file (`basin_outlet_fname`).

```
# lsdtt-chi-analysis parameter file (West Virginia, k_sn by lithology)
write_fname: wvks_lithprimary
read_fname: wvfront_DEM.tif
raster_is_filled: false
threshold_contributing_pixels: 200
min_channel_pixels: 30
m_over_n: 0.45
chi_A0: 1
n_largest_basins: 0
basin_outlet_fname: wvaf_selected_basin_outlets.csv
basin_outlet_snapping_method: none
fit_chi_networks: true
write_ksn_basin_summary: true
lithology_fname: wvfront_macrostrat_units.fgb
lithology_field_name: lith_primary
lithology_crosswalk_fname: wvks_crosswalk_lith_primary.csv
lithology_min_reach_pixels: 10
lithology_bootstrap_replicates: 1000
lithology_bootstrap_random_seed: 5489
lithology_bootstrap_cluster: reach
lithology_contact_buffer_metres: 0
write_lithology_code_raster: true
write_fit_tables: true
write_chi_fit_segment_lines: true
```

Measured on that run (live data, not a CTest; the same values are recorded in the
developer testing notes). Reach median k_sn with the 95 % bootstrap interval, buffer
0 m, m/n 0.45, 1000 replicates by reach, seed 5489:

| Unit | Reaches (qualifying) | Reach median k_sn [95 % interval] |
|---|---|---|
| sandstone | 358 (190) | 65.34 [60.01, 69.59] |
| siltstone | 68 (44) | 54.09 [45.86, 60.09] |
| shale | 334 (177) | 41.14 [36.66, 48.97] |
| limestone | 196 (113) | 57.84 [51.62, 64.87] |
| alluvium | 11 (7) | 8.63 [4.65, 14.69] |
| all | 967 (531) | 54.58 [51.05, 58.22] |

The run had 15847 channel nodes and none unmapped. With a 120 m contact buffer,
4677 of those nodes are set aside, and the reach medians become sandstone 64.56,
siltstone 55.99, shale 36.23, limestone 44.31 and alluvium 7.02 (the buffer
changes limestone most). The basin table, the histogram and the node table are
in the same output set.

## Relative erodibility from k_sn by lithologic unit

With `write_relative_erodibility: true` (it needs `lithology_fname` and
`write_ksn_lithology_summary: true`) the CLI turns each unit's reach-median k_sn into
an erodibility relative to a reference unit, and, given an erosion rate, into an
absolute erodibility and a property table for `lsdtt-lithology-erodibility`.

### The relation

With the stream power law E = K A^m S^n, uniform erosion (E the same everywhere) and
steady state, the channel steepness of unit i is k_sn,i = (E / K_i)^(1/n), with
theta = m/n the m/n of the chi run. So

    K_i / K_ref = (k_sn,ref / k_sn,i)^n

for any reference unit, and with a known E (m/yr)

    K_i = E / (A0^m k_sn,i^n),   m = (m/n) n,

which is E / k_sn,i^n with `chi_A0: 1`. k_sn,i is the unit's `reach_k_sn_median`
(the median of `k_sn_ols` over its qualifying lithologic reaches). K is in the
map-area convention: E [m/yr] = K A^m S^n with A in m^2 and S dimensionless gives K
in m^(1 - 2m) yr^-1, for any n. The CSV writes this as text in `K_units`, for example
`m^(0.1) yr^-1` for m = 0.45 and `m^(-0.8) yr^-1` for m = 0.9.

`stream_power_n_values` takes a list (default `1`). Each n gets its own rows,
with m = (m/n) n. The same key sets n for the erodibility calibration below;
`relative_erodibility_n_values` is its old name and still works, with a warning.

### The reference unit

`relative_erodibility_reference_unit` names the reference by its `lith_name`. Left
empty, the reference is the unit with the most qualifying reaches; ties go to the
larger channel length, then to the lower code. The run prints which unit was chosen
and why. `unmapped` and `all` are never a reference and get no row. A name that is
not a unit of the map stops the run, and the message lists the units.

### The interval

The ratio is resampled; it is never formed from the ends of two separate intervals.
Each bootstrap replicate resamples the qualifying reaches of every unit and of the
reference, takes each one's replicate median, and forms the ratio (and K) inside that
replicate. The 95 % interval is the 2.5 and 97.5 % quantiles of the replicate values.
Every n uses the same replicates.

- `lithology_bootstrap_cluster: reach` draws each unit's reaches independently (a
  reach belongs to one unit, so there is nothing to draw jointly). The reference row
  is 1 with the interval [1, 1].
- `lithology_bootstrap_cluster: basin` draws basins once per replicate, from the basins
  that hold a qualifying reach of any unit, and each unit pools its reaches from the
  drawn basins: a basin draw brings the reference's and unit i's reaches together. A
  replicate in which a unit drew no reach is left out for that unit, and `n_bootstrap`
  gives the number of replicates used.
- `lithology_bootstrap_replicates` and `lithology_bootstrap_random_seed` apply. The
  random stream is per m/n (so an interval does not depend on the other m/n values of
  the run) and separate from the k_sn summary's (so switching relative erodibility on
  leaves every k_sn output unchanged).
- A unit with fewer than 5 qualifying reaches has `low_sample = 1`, its point values
  and NaN intervals. If the reference has fewer than 5, every ratio interval is NaN.

The intervals carry only the sampling spread of the reaches. They do not include the
uncertainty of the assumptions below, or of n itself.

### Outputs

| File | Contents |
|---|---|
| `<write_fname>_chi_relative_erodibility_by_lithology_mn_*.csv` | One row per unit and n: `lith_code, lith_name, m_over_n, A_0, n, m, reference_unit, n_reaches_qualifying, reach_k_sn_median, k_ratio, k_ratio_ci_low, k_ratio_ci_high, n_bootstrap, bootstrap_cluster, low_sample`, and with an erosion rate `erosion_rate_m_per_yr, K_abs, K_abs_ci_low, K_abs_ci_high, K_units`. `k_ratio` is K_i / K_ref. |
| `<write_fname>_chi_relative_erodibility_property_table_mn_*_n_*.txt` | With `relative_erodibility_erosion_rate_m_per_yr` > 0: one draft property table per n (`n_1`, `n_2`, `n_1p5`). |

The property table is in the format `lsdtt-lithology-erodibility` reads. For the
West Virginia run below (n = 1):

```
# lsdtt lithology property table v1
# erodibility_convention: map_area
# m_exponent: 0.45
# n_exponent: 1
# provenance: draft from lsdtt-chi-analysis relative erodibility; DEM wvfront_DEM.tif; lithology map wvfront_macrostrat_units.fgb (field lith_primary) with crosswalk wvks_crosswalk_lith_primary.csv; m/n 0.45, A0 1 m^2, n 1, m 0.45; erosion rate 3e-05 m/yr (given, applied everywhere); reference unit sandstone; K = E / (A0^m k_sn^n) from each unit's reach-median k_sn; assumes uniform erosion and steady state
# K units: m^(0.1) yr^-1 (map area in m^2)
# bootstrap: reach clusters; the interval in each reference column is the 95 % percentile interval of K
code,name,erodibility,reference
10,sandstone,4.59121057e-07,"reach median k_sn 65.342244 from 190 qualifying reaches; K/K_ref 1; K 95 % interval [4.3085158e-07, 4.9655754e-07]"
...
```

Give it to `lsdtt-lithology-erodibility` as `lithology_property_table_fname`, with the
same `lithology_fname`, `lithology_field_name` and `lithology_crosswalk_fname` as the
chi run, so the codes and names match. A unit with no qualifying reach gets no row; it
is listed in a `# not_estimated` line, and `lsdtt-lithology-erodibility` refuses a map
that contains it until a row is added by hand.

### Conditions

The CLI prints these once per run. The values hold only when:

- the erosion rate is the same everywhere (uniform uplift and base-level fall);
- the channels are at steady state: a transient knickzone puts its k_sn on whichever
  unit it is crossing;
- each unit's k_sn is its own: reaches next to a contact can carry the k_sn of the unit
  across it (`lithology_contact_buffer_metres` sets nodes near contacts aside);
- the unit has enough qualifying reaches (`low_sample`).

A ratio raised to the power n also raises the relative error of the medians by about
n: on the two-unit known-answer test the ratio is within 0.44 % for n = 1 and 0.89 %
for n = 2 (developer testing notes).

### Example: West Virginia

The run above with these lines added (an E of 3e-5 m/yr, 30 m/Myr, is an
illustrative uniform rate, not a measurement; basin-average erosion rates of order
10-30 m/Myr are reported for the Appalachians):

```
write_relative_erodibility: true
stream_power_n_values: 1,2
relative_erodibility_erosion_rate_m_per_yr: 3e-5
```

Reference sandstone (190 qualifying reaches, the most), 1000 replicates by reach,
seed 5489:

| Unit | Qualifying reaches | K/K_ref, n = 1 | K/K_ref, n = 2 | K (n = 1), m^(0.1) yr^-1 |
|---|---|---|---|---|
| sandstone | 190 | 1 [1, 1] | 1 [1, 1] | 4.591e-07 [4.309e-07, 4.966e-07] |
| siltstone | 44 | 1.208 [1.032, 1.488] | 1.459 [1.065, 2.214] | 5.546e-07 [4.965e-07, 6.536e-07] |
| shale | 177 | 1.588 [1.303, 1.864] | 2.523 [1.698, 3.476] | 7.292e-07 [6.126e-07, 8.370e-07] |
| limestone | 113 | 1.130 [0.9755, 1.281] | 1.276 [0.9516, 1.642] | 5.186e-07 [4.646e-07, 5.850e-07] |
| alluvium | 7 | 7.573 [4.196, 14.65] | 57.35 [17.61, 214.5] | 3.477e-06 [2.043e-06, 6.445e-06] |

All the k_sn-by-lithology outputs of this run are byte-identical to the run without
relative erodibility. The K values scale with E: halving E halves every K and leaves
every ratio unchanged.

## Erodibility calibration from basin erosion rates

With `erosion_rate_fname`, the CLI calibrates the stream-power erodibility K from
basin-averaged erosion rates E (for example cosmogenic Be-10 rates), basin by basin,
across basins, and per lithologic unit. `lsdtt-fetch-erosion-rates` writes the
OCTOPUS rates in this format. It needs `fit_chi_networks: true`. Every other output of
the run is unchanged: the chi fit tables, the k_sn summaries and the relative
erodibility are byte-identical with and without it.

### Inputs and units

The sample file is a CSV with `latitude` and `longitude` (EPSG:4326), or any OGR point
layer in any CRS. One row per sample: the outlet, the rate (`erosion_rate_column`),
its 1-sigma uncertainty (`erosion_rate_uncertainty_column`, same units) and,
optionally, the reported drainage area in km² (`erosion_rate_reported_area_column`).

`erosion_rate_units` is required and has no default: `mm_per_kyr` (= `m_per_myr`,
1e-6 m/yr), `cm_per_kyr` (1e-5), `mm_per_yr` (1e-3) or `m_per_yr`. Mass fluxes
(g cm⁻² yr⁻¹) are refused; divide them by the rock density first. When the file has
an `erosion_rate_units` column, every row must agree with the parameter, or the run
stops naming the first row that disagrees.

### Snapping, basins and nesting

Each sample is snapped to a channel cell (at least `threshold_contributing_pixels`)
within `erosion_rate_snap_max_distance_metres` (rounded up to whole pixels). With a
reported area, `reported_area` snapping takes the nearest channel cell whose drainage
area is within `erosion_rate_area_mismatch_tolerance` of it (a ratio, default 1.25),
otherwise the cell closest in area. Without one it falls back to the nearest channel
cell, and the run says that mis-snapping onto a neighbouring channel cannot then be
detected. `expanding_square`, `downstream`, `largest_area` and `none` are also accepted.

The basin of each snapped outlet must be complete (it may not reach the DEM edge or
nodata). The DEM area is compared with the reported area; a ratio outside 1/tol..tol
is `area_mismatch`, and the sample is excluded (`erosion_rate_area_mismatch_action:
exclude`) or kept and flagged (`flag`).

A sample whose outlet lies upslope of another sample's outlet is nested in it;
samples on one outlet are duplicates. Nesting groups are the connected components of
these relations. They are the bootstrap clusters, and the number of groups is the
effective number of independent basins. `erosion_rate_nested_basins` keeps all
samples (`keep`), only those nested in no other (`outermost`), or only those holding
no other (`innermost`).

### Basin k_sn

Each basin is extracted and fitted on its own chi network (the m/n and A0 of the run,
`fit_critical_divergence_metres`), with no `n_largest_basins` or `min_channel_pixels`
filter. Its k_sn is the median `k_sn_ols` of its qualifying reaches, where reaches end
only at confluences, sources and the outlet: the value does not depend on a lithology
map. A reach qualifies with at least `erosion_rate_min_reach_pixels` nodes. A basin
needs `erosion_rate_min_reaches_per_basin` qualifying reaches (default 5); fewer gives
`too_few_reaches`. The 95 % interval is a percentile bootstrap of the reach median,
and sigma_ln_k = (ln high - ln low) / 3.92. The chi-z slope through the origin (with
its R²), the segment median, the chi-weighted mean and the reach IQR / median are
written alongside as diagnostics; `erosion_rate_ksn_statistic` can choose the chi-z
slope or the segment median instead.

### K

E = K A0^m k_sn^n and K = E / (A0^m k_sn^n), with m = (m/n) n. K is in the map-area
convention, m^(1 - 2m) yr^-1 (written in `K_units`), and does not depend on A0. For
each n of `stream_power_n_values`:

- **Per basin**: K of each sample. The interval is a Monte Carlo over the bootstrap
  replicates. In each replicate, k_sn is the basin's replicate median and E is drawn
  from N(E, sigma_E), redrawn while it is <= 0. A sample without sigma_E
  (`missing_uncertainty`) gets the point value only and enters no fit.
- **Fixed n, across basins**: the weighted mean of ln K_i, with weights
  1 / (sigma_lnE² + n² sigma_lnk²). The fit file also gives `ln_K_dispersion` (the
  weighted SD of ln K_i) and `corr_lnK_lnk`. At the wrong n, K_i trends with k_sn, and
  both grow.
- **Fitted n**: a York regression of ln E on ln k_sn with errors in both, only with at
  least `erosion_rate_min_groups_for_n_fit` nesting groups (default 8) and a k_sn range
  (max / min) of at least `erosion_rate_min_ksn_range_ratio` (default 2). Otherwise
  the row says why (`not_fitted_too_few_groups`, `not_fitted_narrow_ksn_range`) and
  stdout gives a first-order SE(n) for the basins at hand. Both York errors are
  written: `n_se_york` and `ln_K_se` assume the stated sigmas are complete;
  `*_mswd_scaled` multiply them by sqrt(MSWD). n and ln K are strongly correlated
  (`corr_n_lnK`); `E_at_median_ksn_m_per_yr` is the better-constrained pivot.

The 95 % intervals of the fixed-n and fitted rows are percentile bootstraps over
nesting groups. A York replicate with fewer than 3 distinct groups, a zero k_sn range
or no convergence is skipped, and `n_bootstrap` counts the replicates used.
`erosion_rate_systematic_relative_uncertainty` (default 0) adds a common relative
error of E, such as the production-rate term, to every K interval in quadrature.

### K per lithologic unit

With `lithology_fname`, for `erosion_rate_unit_k_methods`:

- `mono`: samples whose dominant unit covers at least
  `erosion_rate_mono_lithology_min_area_fraction` (default 0.8) of the basin area. The
  fixed-n estimator is applied per unit.
- `unit_in_basin`: for each unit and used basin with at least
  `erosion_rate_min_unit_reaches_per_basin` (default 3) qualifying lithologic reaches
  (the map's reaches, as for relative erodibility), the unit's reach-median k_sn in
  that basin with the basin's E. The fixed-n estimator is applied per unit. This
  assumes E is uniform across the units of a basin.

`low_sample` is 1 below 3 nesting groups. `N_basins` counts sample entries, so
duplicates on one outlet count once each but form one group.
`erosion_rate_property_table_method` (default `unit_in_basin`) writes one draft
property table per n for `lsdtt-lithology-erodibility`. Units without K are listed in
its `# not_estimated` line. A mixing model over area fractions is not available.

### Outputs

Per m/n (`<w>` = `write_fname`):

| File | Columns |
|---|---|
| `<w>_erodibility_calibration_samples_mn_*.csv` | one row per input sample, in input order: `sample_id, input_row, input_latitude, input_longitude, snapped_latitude, snapped_longitude, snapping_method, snap_distance_m, outlet_flat, status, status_detail, reported_area_km2, dem_area_km2, area_ratio, area_check, nesting_group, nested_in, n_samples_inside, duplicate_outlet_of, erosion_rate_input, erosion_rate_uncertainty_input, erosion_rate_input_units, erosion_rate_m_per_yr, erosion_rate_uncertainty_m_per_yr, m_over_n, A_0, k_sn_statistic, k_sn, k_sn_ci_low, k_sn_ci_high, n_reaches, n_reaches_qualifying, k_sn_reach_iqr_over_median, k_sn_chi_z_slope, chi_z_r_squared, k_sn_segment_median, k_sn_chi_weighted_mean, dominant_lith_code, dominant_lith_name, dominant_area_fraction, dominant_channel_length_fraction, unmapped_area_fraction` |
| `<w>_erodibility_calibration_basin_K_mn_*.csv` | per sample with a K and per n: `sample_id, status, nesting_group, m_over_n, A_0, n, m, erosion_rate_m_per_yr, erosion_rate_uncertainty_m_per_yr, k_sn, k_sn_ci_low, k_sn_ci_high, K, K_ci_low, K_ci_high, n_bootstrap, K_units` |
| `<w>_erodibility_calibration_fit_mn_*.csv` | one `fixed_n` row per n and one `fitted_n` row: `fit_type, method, status, N_basins, N_groups, ksn_range_ratio, m_over_n, A_0, n, n_ci_low, n_ci_high, n_se_york, n_se_york_mswd_scaled, m, K, K_ci_low, K_ci_high, ln_K_se, ln_K_se_mswd_scaled, corr_n_lnK, E_at_median_ksn_m_per_yr, mswd, ln_K_dispersion, corr_lnK_lnk, n_bootstrap, bootstrap_cluster, K_units` |
| `<w>_erodibility_calibration_by_lithology_mn_*.csv` | with a map, one row per method, n and unit: `lith_code, lith_name, method, m_over_n, A_0, n, m, N_basins, N_groups, K, K_ci_low, K_ci_high, n_bootstrap, low_sample, status, purity_threshold, K_units` |
| `<w>_erodibility_calibration_property_table_mn_*_n_*.txt` | with a map: the draft property table |

Once per run: `<w>_erodibility_calibration_basins.fgb` (or `.gpkg`: the polygons of
every complete sample basin in the DEM CRS, with `sample_id, status, nesting_group,
dem_area_km2, reported_area_km2, area_ratio`; the status is the one known before k_sn,
`kept` for a basin that goes on to the k_sn step, whose final status is in the samples
file). With a map, also
`<w>_erodibility_calibration_basin_lithology.csv` (`sample_id, lith_code, lith_name,
area_km2, area_fraction, channel_length_m, channel_length_fraction`).

`status` is one of `outside_dem`, `snap_failed`, `incomplete_basin`, `bad_rate`,
`area_mismatch`, `nested_excluded`, `too_few_reaches`, `missing_uncertainty` or `used`
(the first that applies). `area_check` (`ok`, `mismatch`, `no_reported_area`) is
reported separately.

### Conditions

The calibration assumes steady state and an erosion rate that is uniform within each
basin. CRN rates are weighted by quartz flux, and this is not corrected (units poor in
quartz are under-represented). A basin above a migrating knickpoint carries a relict
k_sn on top of a recent E; MSWD well above 1 shows scatter beyond the stated sigmas.
Prefer K at a fixed n where possible. These conditions are printed with every run.

### Example: West Virginia

The OCTOPUS rates of the box -80.5/-78.5/38.3/39.7 (`lsdtt-fetch-erosion-rates`,
scheme `st`; 23 samples with outlets in the box), a COP30 DEM fetched with that run's
`_octopus_dem_bounds.fgb` (3602 x 5469 pixels at 30 m, EPSG:32617) and Macrostrat
units (`lith_primary`, every value its own unit):

```
read_fname: wvcal_DEM.tif
write_fname: wvcal_chi
raster_is_filled: false
threshold_contributing_pixels: 200
basin_outlet_fname: wvoct_large_erosion_rates.csv
basin_outlet_snapping_method: downstream
basin_outlet_contributing_pixels_threshold: 200
n_largest_basins: 0
m_over_n: 0.45
chi_A0: 1
fit_chi_networks: true
lithology_fname: wvcal_lith_macrostrat_units.fgb
lithology_field_name: lith_primary
lithology_crosswalk_fname: wvcal_crosswalk_lith_primary.csv
write_relative_erodibility: true
erosion_rate_fname: wvoct_large_erosion_rates.csv
erosion_rate_units: mm_per_kyr
stream_power_n_values: 1,2
```

Measured on 2026-10-19 (m/n 0.45, A0 1 m²):

- **Samples:** 23 read. 17 used in 14 nesting groups. Excluded: 2 `incomplete_basin`
  (the two basins larger than 1000 km², left out of the DEM box), 1 `area_mismatch`
  (S046WTS032, ratio 0.455) and 3 `too_few_reaches` (S046WTS044, 045 and 054, 0.9 to
  2.0 km²). S046WTS009 to 012 share one outlet (one group).
- **Fixed n:** n = 1 gives K = 3.081e-07 [2.372e-07, 4.153e-07] m^(0.1) yr^-1, with
  `ln_K_dispersion` 0.520 and `corr_lnK_lnk` -0.404. n = 2 gives K = 8.581e-09
  [5.253e-09, 1.436e-08] m^(-0.8) yr^-1, with `ln_K_dispersion` 0.893 and
  `corr_lnK_lnk` -0.805.
- **Fitted n:** the fit ran with 14 groups and a k_sn range of 5.12x. n = 0.658
  [0.130, 2.025] (York SE 0.124, x sqrt(MSWD) 0.260, MSWD 4.39) and K = 1.005e-06
  [9.14e-09, 7.75e-06] m^(0.408) yr^-1, with `corr_n_lnK` -0.991.

The run took 20.9 s wall time with a peak RSS of 1.67 GB, on a laptop (WSL2); laptop
timings are thermally limited.

## Outputs

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_fit_tables` | bool | `true` | Write chi network fit node and segment tables to CSV. |
| `write_chi_fit_segment_lines` | bool | `false` | Write fitted chi segments as polylines. One feature per piecewise-linear segment. Attributes: `basin_index`, `outlet_flat`, `segment_id`, `branch_id`, `k_sn`, `downstream_chi`, `upstream_chi`, `downstream_elevation`, `upstream_elevation`, `n_nodes`. Written once per m/n value. |
| `write_disorder_sample_tables` | bool | `false` | Write detailed disorder sample tables for each m/n value. |
| `calculate_spatial_index` | bool | `false` | Calculate a spatial index for vector outputs. Unlike smaller vector-export CLIs, chi analysis leaves this off because its per-m/n layers can be large; enabling it speeds GIS queries but buffers each layer in memory. |

---

## Minimal example parameter file

```
# lsdtt-chi-analysis parameter file — chi data maps only
read_path:  /data/dems/
read_fname: catchment.tif

raster_is_filled:           false
threshold_contributing_pixels: 1000

write_chi_data_maps:        true
chi_data_maps_write_vector_format: flatgeobuf
write_chi_data_map_lines:   true
m_over_n:                   0.45
```

```
# lsdtt-chi-analysis parameter file — m/n sweep with bootstrap
read_path:  /data/dems/
read_fname: catchment.tif

threshold_contributing_pixels: 1000
basin_outlet_fname:         outlets.csv

fit_chi_networks:   true
chi_A0:             1.0
mn_min:             0.2
mn_max:             0.7
mn_step:            0.02

run_bootstrap:      true
bootstrap_replicates: 500
```

---

## References

Perron, J. T., and Royden, L. (2013). An integral approach to bedrock river profile analysis. *Earth Surface Processes and Landforms*, 38(6), 570–576. <https://doi.org/10.1002/esp.3302> (the chi coordinate)

Mudd, S. M., Clubb, F. J., Gailleton, B., and Hurst, M. D. (2018). How concave are river channels? *Earth Surface Dynamics*, 6, 505–523. <https://doi.org/10.5194/esurf-6-505-2018> (the disorder method for m/n, which follows Hergarten et al., 2016)

How chi is computed, how the profiles are segmented, and how knickpoints and the disorder metric are found are described in [Chi analysis](../design_and_algorithms/chi_analysis.md).

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

### Threading / routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all available hardware threads. |
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four-character string: north, east, south, west. o = open, n = no-flux, p = periodic, b = forced outlet. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |
| `fill_flat_resolution_choice` | `string` | `inline_gradient` | `param` | - | How flats are made drainable: inline_gradient, barnes2014b_gradient or barnes2014b_float_step. | inline_gradient (default, unchanged): the Priority-Flood fill imposes fill_min_slope during the flood. The barnes2014b options make flow converge across flats (Barnes et al. 2014, doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses fill_min_slope (m/m, must be &gt; 0) as the minimum gradient; flats are raised about 3-4x more, capped below their rims, and capped flats are reported. barnes2014b_float_step is equivalent to RichDEM ResolveFlats on float32 in its first pass: it ignores fill_min_slope and its slopes depend on elevation (of order ULP(z) / cell size). With either barnes2014b option the result is always routable: cells the first pass leaves without a lower neighbour are refilled and resolved again (up to 4 rounds), and any still left are raised by a final priority-flood pass by one float32 step each; the report gives the rounds and the cells changed after the first pass. |

### DEM preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope (m/m) imposed by pit filling or flat resolution. | inline_gradient (default): imposed during the Priority-Flood fill, so used only when raster_is_filled is false. barnes2014b_gradient: the minimum gradient across flats, also on a raster_is_filled input. barnes2014b_float_step: ignored. |

### Basin outlet input

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_outlet_fname` | `string` | `NULL` | `param` | - | File of basin outlet points (CSV, FlatGeobuf, or GeoPackage). | If NULL (default), the largest basin outlet in the DEM is used. CSV files must have 'latitude' and 'longitude' columns (EPSG:4326). |
| `basin_outlet_snapping_method` | `string` | `downstream` | `param` | - | Optional snapping of supplied outlet points to qualifying channel pixels. | Allowed values: none, expanding_square, downstream, largest_area. largest_area takes the pixel with the most contributing pixels in the square of half-width basin_outlet_snap_max_search_radius_pixels, if it reaches basin_outlet_contributing_pixels_threshold. |
| `basin_outlet_contributing_pixels_threshold` | `int` | `1000` | `param` | - | Contributing-pixel threshold used for optional outlet snapping. | Total upslope contributing pixels required at the snapped cell. |
| `basin_outlet_snap_min_search_radius_pixels` | `int` | `0` | `param` | - | Minimum expanding-square snapping search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square. |
| `basin_outlet_snap_max_search_radius_pixels` | `int` | `25` | `param` | - | Maximum expanding-square snapping search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square or largest_area (the half-width of the square searched). |

### Diagnostic data maps

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_chi_data_maps` | `bool` | `false` | `switch` | - | Extract chi networks and print per-node data maps for diagnostic purposes. | Writes nodes (chi, drainage_area, flow_distance, elevation, source_key, basin_key), valid basin outlets, source-head points, and basin polygons. Basin keys are 0-based; source keys increment globally across all basins. Runs before fit_chi_networks so extraction can be verified independently. |
| `chi_data_maps_write_vector_format` | `string` | `csv` | `param` | - | Output format for chi data map point files. | csv: latitude/longitude in EPSG:4326 (requires raster CRS). gpkg: GeoPackage in raster CRS. flatgeobuf/fgb: FlatGeobuf in raster CRS. Basin polygons are always vector (flatgeobuf fallback when csv requested). |
| `write_chi_data_map_lines` | `bool` | `false` | `switch` | - | Write chi channels as polylines, one feature per source_key. | Each source_key represents the continuous channel from a headwater source to the basin outlet following the primary-donor (longest-path) convention. Attributes: source_key, basin_key. Only complete basins are written. Keys match those in the chi data map point output. |

### Chi coordinate raster

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_chi_coordinate_raster` | `bool` | `false` | `switch` | - | Write a raster of the chi coordinate for every pixel of every complete basin. | Chi is integrated upstream from each retained basin outlet (chi = 0 at the outlet) across the whole upslope domain, so hillslope pixels get a chi value as well as channel pixels: this is the chi map used for divide-migration and cross-divide chi-contrast analysis. Basins are selected exactly as for the other chi outputs (threshold_contributing_pixels, then min_channel_pixels, then n_largest_basins), except that the completeness filter always runs: basins whose upslope area touches nodata or the DEM edge are never written, and every pixel outside the retained basins is nodata. Uses one m/n value per run (the first value of m_over_n, or of the mn sweep). Output: &lt;write_fname&gt;_chi_coordinate_mn_*.tif in write_output_format. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Source extraction

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `min_channel_pixels` | `int` | `20` | `param` | - | Reject basins whose extracted channel has fewer than this many pixels. | After chi extraction, any basin with fewer than this number of chi network nodes (channel pixels) is discarded before fitting, disorder, and inversion. Very short channels cannot be segmented meaningfully and can destabilise fitting. Set to 0 or 1 to keep all basins. |
| `n_pixels_move_outlet_upstream` | `int` | `10` | `param` | - | Move interior (non-edge) auto-detected outlets this many pixels upstream. | When no basin_outlet_fname is supplied, outlets are auto-detected. An outlet not on the DEM bounding-box edge drains to interior nodata (e.g. an island coast where the sea is masked); it is walked this many pixels up the main stem so the resulting basin clears the coastal nodata fringe and can pass the completeness test. Set to 0 to disable. |
| `n_largest_basins` | `int` | `5` | `param` | - | Keep only the N largest complete basins for all chi outputs. | After the min_channel_pixels filter, incomplete basins (upslope area touches the DEM edge / nodata) are dropped and only the N largest remaining complete basins, ranked by total basin pixels (drainage area), are analysed. This applies consistently to data maps, fitting, disorder, and inversion so basin keys stay stable. Set to 0 to keep all complete basins (no size cap). |

### Chi extraction / fitting

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fit_chi_networks` | `bool` | `false` | `switch` | - | Run the piecewise-linear chi-network fitting and disorder sweep. | When false (default) only write_chi_data_maps (if enabled) and the inversion (if enabled) run. Set to true for the full m/n sweep. |
| `chi_A0` | `float` | `1.0` | `param` | `fit_chi_networks` | Reference drainage area A0 for chi integration (m2). | Chi = integral of (A0/A)^(m/n) ds from the outlet. LEAVE THIS AT 1.0 UNLESS YOU KNOW WHY YOU ARE CHANGING IT. The slope of a chi-elevation profile equals the channel steepness index k_sn ONLY when A0 = 1 m2; for any other A0 the gradient is k_sn scaled by A0^(m/n), so it is a steepness index referenced to that area rather than k_sn itself. The rescaling is silent - the profile looks the same and only the numbers on the axis change - and it is large: A0 = 1e6 m2 at m/n = 0.45 scales every gradient by about 631. Other packages default differently (TopoToolbox's chitransform uses a0 = 1e6), so state A0 explicitly whenever a value is compared across tools or quoted in a paper. |
| `mn_min` | `float` | `0.3` | `param` | `fit_chi_networks` | Minimum m/n ratio for the disorder sweep. | - |
| `mn_max` | `float` | `0.6` | `param` | `fit_chi_networks` | Maximum m/n ratio for the disorder sweep. | - |
| `mn_step` | `float` | `0.02` | `param` | `fit_chi_networks` | Step size for the m/n disorder sweep. | Must be &gt; 0. |
| `m_over_n` | `string` | - | `param` | `fit_chi_networks` | Optional explicit m/n value(s) for chi extraction and fitting. | Comma-separated list (e.g. 0.45 or 0.35,0.45,0.55). When set, this overrides mn_min/mn_max/mn_step for chi data maps, fitting, and inversion. |
| `fit_critical_divergence_metres` | `float` | `5.0` | `param` | `fit_chi_networks` | Critical elevation divergence (m) used to split chi-elevation branches. | - |

### Disorder

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_subset_sensitivity` | `bool` | `false` | `switch` | - | Run subset sensitivity analysis on chi disorder. | - |
| `subset_n_side_branches_per_sample` | `int` | `3` | `param` | - | Number of side branches selected per subset sensitivity sample. | - |
| `run_jackknife` | `bool` | `false` | `switch` | - | Run jackknife (leave-one-out) analysis on chi disorder. | - |
| `run_bootstrap` | `bool` | `false` | `switch` | - | Run bootstrap analysis on chi disorder. | - |
| `bootstrap_replicates` | `int` | `1000` | `param` | - | Number of bootstrap replicates for disorder uncertainty estimation. | - |
| `bootstrap_sample_full_tributary_count` | `bool` | `true` | `param` | - | Whether each bootstrap replicate samples the full number of side branches. | - |
| `bootstrap_n_side_branches_per_replicate` | `int` | `3` | `param` | - | Number of side branches per bootstrap replicate. | - |
| `bootstrap_with_replacement` | `bool` | `true` | `param` | - | Whether bootstrap resampling is done with replacement. | - |
| `bootstrap_random_seed` | `int` | `5489` | `param` | - | Random seed for bootstrap resampling. | - |

### Inversion

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_inversion` | `bool` | `false` | `switch` | - | Run spatial chi k_sn / u* inversion over the m/n sweep. | Requires Eigen3 at build time. If Eigen3 is not available, this flag is ignored and inversion is skipped with a warning. |
| `inversion_coarse_cell_size` | `float` | `2000.0` | `param` | - | Cell size (m) of the coarse grid used in chi k_sn inversion. | - |
| `inversion_roughness_weight` | `float` | `10000.0` | `param` | - | Tikhonov regularisation weight for spatial smoothness. | - |
| `inversion_diagonal_jitter` | `float` | `1e-8` | `param` | - | Small diagonal perturbation added to the normal equations. | - |
| `inversion_use_mae_for_model_selection` | `bool` | `true` | `param` | - | Use MAE (true) or RMSE (false) for m/n model selection. | - |
| `write_inversion_rasters` | `bool` | `false` | `switch` | - | Write u* and support-count rasters for each basin. | - |

### Outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_fit_tables` | `bool` | `true` | `switch` | - | Write chi network fit node and segment tables to CSV. | - |
| `write_chi_tributary_table` | `bool` | `true` | `switch` | `fit_chi_networks` | Write one row per tributary (branch) of the chi fit, with its drainage area where it joins. | A tributary is a branch of the fit: from its head down to the node where it joins its parent branch (branch 0 is the main stem). Columns: basin_index, basin_key and source_key (the chi data maps' keys; source_key is the key of its head and is also written on every fit segment and node), branch_id, parent_branch_id, branch_order (0 main stem, 1 its tributaries, ...), head and junction flats and x / y (DEM CRS), junction_drainage_area_m2 (at the tributary's last node before the junction), basin_drainage_area_m2, area_fraction, area_rank_in_basin (1 = largest tributary), length_m and chi_span (head to junction), n_segments, k_sn min / median / max with the steepest segment's id and mid-chi, and n_descendants (tributaries draining into it at any depth). Written once per m/n. Output: &lt;write_fname&gt;_chi_tributaries_mn_*.csv. |
| `write_chi_knickpoint_table` | `bool` | `true` | `switch` | `fit_chi_networks` | Write the nodes where fitted k_sn segments meet (knickpoint candidates) and threshold suggestions. | One row per boundary: within a branch, where a segment meets the next segment downstream; and, flagged at_junction = 1, where a tributary's first segment meets the parent-branch segment at the junction. Columns include x / y, chi, elevation, flow distance, drainage area, k_sn upstream and downstream, delta_ksn (upstream minus downstream: positive where the channel steepens upstream) and log_ratio_ksn (log10 of their ratio; empty if either k_sn &lt;= 0). A second file suggests thresholds for a major knickpoint over the within-branch rows: for \|delta_ksn\| and \|log_ratio_ksn\|, the 90th and 95th percentiles, the median plus 3 MAD (MAD scaled by 1.4826) and Otsu's two-class split; and an outlier k_sn for segments (Q3 + 3 IQR). Written once per m/n. Outputs: &lt;write_fname&gt;_chi_knickpoints_mn_*.csv, &lt;write_fname&gt;_chi_knickpoint_thresholds_mn_*.csv. |
| `write_ksn_basin_summary` | `bool` | `true` | `switch` | - | Write a per-basin channel-steepness (k_sn) summary CSV. | One row per fitted basin with the number of segments, number of main-stem segments, total chi length, and mean / median / min / max / chi-weighted-mean k_sn plus main-stem mean and median k_sn. Written once per m/n value tested. Output: &lt;write_fname&gt;_chi_ksn_basin_summary_mn_*.csv. |
| `write_ksn_lithology_summary` | `bool` | `true` | `switch` | `lithology_fname` | With lithology_fname: write k_sn by lithologic unit. | Per m/n: &lt;write_fname&gt;_chi_ksn_by_lithology_mn_*.csv (per unit, plus unmapped and all: node, segment, basin and reach counts, node k_sn statistics, the median and quartiles of lithologic-reach k_sn with a bootstrap interval of the median, and length-weighted node k_sn quantiles), _chi_ksn_by_lithology_and_basin_mn_*.csv (the reach statistics per unit and basin, no interval), _chi_ksn_histogram_by_lithology_mn_*.csv (common bins for every unit), _chi_ksn_lithology_nodes_mn_*.csv (every channel node with its segment k_sn, unit, distance to the nearest contact and reach), _chi_lithology_reaches_mn_*.csv (every lithologic reach: OLS k_sn of elevation on chi, its standard error and R^2, end types) and _chi_lithology_segments_mn_*.csv (every fitted segment with its length-weighted majority unit and purity). Once per run: &lt;write_fname&gt;_lithology_units.csv (the unit legend with pixel counts). A lithologic reach is a run of channel nodes in one unit on one branch; reaches are the sampling unit of the per-unit k_sn statistics. Node k_sn is its fitted segment's k_sn; segments are not split at contacts. Nodes are not independent samples and the interval does not account for correlation between reaches on one channel. |
| `lithology_ksn_histogram_n_bins` | `int` | `30` | `param` | `write_ksn_lithology_summary` | Number of k_sn histogram bins, shared by every unit. | The bins span the 1st-99th percentile of all nodes' k_sn; nodes outside are counted in a row below (bin_index -1) and above (bin_index = n_bins). |
| `lithology_ksn_histogram_log_bins` | `bool` | `false` | `param` | `write_ksn_lithology_summary` | Space the k_sn histogram bins evenly in log10(k_sn). | k_sn &lt;= 0 then falls in the row below the range. |
| `write_chi_fit_segment_lines` | `bool` | `false` | `switch` | - | Write fitted chi segments as a vector polyline file. | One polyline feature per piecewise-linear segment. Attributes: basin_index, outlet_flat, segment_id, branch_id, k_sn (dz/dchi), downstream_chi, upstream_chi, downstream_elevation, upstream_elevation, n_nodes. Written once per m/n value tested. |
| `calculate_spatial_index` | `bool` | `false` | `io` | - | Calculate a spatial index for vector outputs. | Unlike the smaller vector-export CLIs, chi analysis leaves this off by default because its per-m/n outputs can be large. Enabling it speeds up later GIS queries but buffers each layer in memory before writing. |
| `write_disorder_sample_tables` | `bool` | `false` | `switch` | - | Write detailed disorder sample tables for each m/n value. | - |
| `write_output_format` | `string` | `cog` | `io` | - | Raster output format for the chi coordinate and inversion rasters. | Allowed values: cog, gtiff, same_as_input. |

### Lithology input

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `lithology_fname` | `string` | - | `param` | - | Lithology map (in read_path) for k_sn by lithologic unit: a polygon vector file or a categorical code raster. Empty (default) = off. Needs fit_chi_networks: true. | Vectors: any OGR polygon layer (GeoPackage, FlatGeobuf, GeoJSON, Shapefile), e.g. lsdtt-fetch-lithology output; it is reprojected to the DEM's CRS. Rasters: integer lithology codes, warped onto the DEM grid with nearest neighbour if not on it already. |
| `lithology_source` | `string` | `auto` | `param` | - | auto, vector_map or code_raster. | auto decides from the file extension (.tif, .tiff, .asc, .bil, .img, .vrt, .flt, .nc are rasters). |
| `lithology_layer_name` | `string` | - | `param` | `lithology_fname` | Vector layer to read. Empty (default) = the first layer. | - |
| `lithology_field_name` | `string` | `lith_primary` | `param` | `lithology_fname` | Attribute of the vector layer that holds the lithology. | Each distinct value is a lithologic unit unless lithology_crosswalk_fname groups them. lith_primary is the first lithology of each unit in lsdtt-fetch-lithology output; name (the map unit) or lith_class group more finely or more coarsely. Units are numbered within a run; compare runs by unit name. |
| `lithology_crosswalk_fname` | `string` | - | `param` | `lithology_fname` | Optional CSV (in read_path) with columns value,code[,name] mapping map values to units. | When given, every map value present must be in it. Use it to group many map values into a few units; the optional name column names each code (else the values joined with '; '). |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_lithology_code_raster` | `bool` | `false` | `switch` | - | Write the lithologic unit code of every pixel. | Output: &lt;write_fname&gt;_lithology_code.tif. |
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). CSV is available only through selectors for point/table outputs. Per-output overrides exist only where a container genuinely differs. Raster outputs are governed separately by write_output_format. |

### Lithology summaries

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `lithology_contact_buffer_metres` | `float` | `0.0` | `param` | `lithology_fname` | Exclude channel nodes closer than this to a mapped lithologic contact from the per-unit k_sn statistics. | Excluded nodes are still written, flagged near_contact = 1. 0 keeps every node. Mapped contacts carry a positional error of about 0.5 mm at map scale (about 125 m at 1:250,000). |
| `lithology_min_reach_pixels` | `int` | `10` | `param` | `lithology_fname` | Minimum number of channel nodes for a lithologic reach to enter the per-unit statistics. | A lithologic reach is a run of channel nodes in one unit on one branch. It must also span a chi interval greater than 0. |
| `lithology_bootstrap_cluster` | `string` | `reach` | `param` | `lithology_fname` | Resampling unit of the bootstrap interval of the per-unit median k_sn: reach or basin. | reach resamples lithologic reaches within a unit; basin resamples whole basins. Reaches on one trunk channel are correlated, so the reach interval is a lower bound on the uncertainty. |
| `lithology_bootstrap_replicates` | `int` | `1000` | `param` | `lithology_fname` | Number of bootstrap replicates for the per-unit confidence intervals. | 0 = no bootstrap; the interval columns are written as NaN. |
| `lithology_bootstrap_random_seed` | `int` | `5489` | `param` | `lithology_fname` | Random seed for the per-unit bootstrap. | lsdtt-chi-analysis: each m/n gets its own std::mt19937_64 from seed_seq{seed, round(m/n * 10000)}, so an interval does not depend on the other m/n values in the list. lsdtt-channel-network-analysis: one std::mt19937_64 per run. A fixed seed reproduces the intervals. Same default as bootstrap_random_seed. |

### Relative erodibility

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_relative_erodibility` | `bool` | `false` | `switch` | `write_ksn_lithology_summary` | Write the erodibility of each lithologic unit relative to a reference unit, from the units' reach-median k_sn. | Needs lithology_fname and write_ksn_lithology_summary: true. Under uniform erosion at steady state with the stream power law E = K A^m S^n (m = (m/n) n), K_i / K_ref = (k_sn,ref / k_sn,i)^n, from each unit's reach-median k_sn (reach_k_sn_median of _chi_ksn_by_lithology_mn_*.csv). The 95 % interval resamples the ratio: each bootstrap replicate resamples the reaches of the unit and of the reference (lithology_bootstrap_cluster: reach draws each unit's reaches independently; basin draws whole basins, bringing both units' reaches together) and forms the ratio in that replicate; the interval is the 2.5 and 97.5 % quantiles of the replicate ratios. lithology_bootstrap_replicates and lithology_bootstrap_random_seed apply; the stream is per m/n and separate from the k_sn summary's. Units with fewer than 5 qualifying reaches get low_sample = 1 and NaN intervals. Unmapped and all get no row. Output per m/n: &lt;write_fname&gt;_chi_relative_erodibility_by_lithology_mn_*.csv, one row per unit and n: lith_code, lith_name, m_over_n, A_0, n, m, reference_unit, n_reaches_qualifying, reach_k_sn_median, k_ratio (K_i / K_ref), k_ratio_ci_low, k_ratio_ci_high, n_bootstrap, bootstrap_cluster, low_sample, and with an erosion rate erosion_rate_m_per_yr, K_abs, K_abs_ci_low, K_abs_ci_high, K_units. The ratios hold only where the erosion rate is uniform and the channels are at steady state; reaches near a contact mix the k_sn of two units, and a transient knickzone puts its k_sn on whichever unit it crosses. |
| `relative_erodibility_reference_unit` | `string` | - | `param` | `write_relative_erodibility` | Name of the reference lithologic unit (K_ref) for write_relative_erodibility. | Empty (default): the unit with the most qualifying reaches; ties go to the larger channel length, then the lower code. Give the unit name as in lith_name (the crosswalk name, or the map value). A name that is not a unit of the map stops the run, listing the units. unmapped and all are never a reference. |
| `relative_erodibility_erosion_rate_m_per_yr` | `double` | `0.0` | `param` | `write_relative_erodibility` | Uniform erosion rate E (m/yr) for absolute erodibility; 0 = not given. | When &gt; 0, also K_i = E / (A0^m k_sn,i^n) for every unit (map-area convention, units m^(1 - 2m) yr^-1, written in K_units), with its interval from the same replicates, and one draft property table per n, &lt;write_fname&gt;_chi_relative_erodibility_property_table_mn_*_n_*.txt, in the format lsdtt-lithology-erodibility reads (erodibility_convention map_area, m_exponent, n_exponent, and a provenance line). Use it there with the same lithology map, field and crosswalk to make a K raster. Units with no qualifying reach have no row and are listed in a # not_estimated line. |

### Stream power

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `stream_power_n_values` | `string` | `1` | `param` | - | Slope exponent(s) n of the stream power law E = K A^m S^n, a comma list (e.g. 1,2). | Used by write_relative_erodibility and the erodibility calibration (erosion_rate_fname): each n gives its own rows, with m = (m/n) n for the m/n of the chi run. Every value must be &gt; 0. Replaces relative_erodibility_n_values (still accepted, with a warning). |

### Erosion-rate samples

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `erosion_rate_fname` | `string` | - | `param` | - | Basin-averaged erosion-rate samples: a CSV with latitude/longitude (EPSG:4326) or an OGR point layer in any CRS. Empty = off. | One row per sample: the basin outlet, the erosion rate and its 1-sigma uncertainty, optionally the reported drainage area (lsdtt-fetch-erosion-rates writes this format from OCTOPUS). Needs fit_chi_networks: true. Each sample is snapped to the channel network (threshold_contributing_pixels), its basin delineated and its k_sn measured on its own fitted chi network; K is calibrated per basin, at each stream_power_n_values n across basins, as a fitted (n, K), and per lithologic unit with lithology_fname. Outputs per m/n: &lt;write_fname&gt;_erodibility_calibration_samples_mn_*.csv, _basin_K_mn_*.csv, _fit_mn_*.csv and, with a map, _by_lithology_mn_*.csv and _property_table_mn_*_n_*.txt; once: _erodibility_calibration_basins.&lt;fgb\|gpkg&gt; and, with a map, _erodibility_calibration_basin_lithology.csv. The main chi outputs do not change. |
| `erosion_rate_layer_name` | `string` | - | `param` | `erosion_rate_fname` | Layer of a vector erosion-rate file. | Empty takes the only (or first) layer. Not used for CSV. |
| `erosion_rate_sample_id_column` | `string` | `sample_id` | `param` | `erosion_rate_fname` | Column holding the sample identifier. | When the column is missing or a value is blank the sample is called row_&lt;n&gt; (n counts data rows from 1). |
| `erosion_rate_latitude_column` | `string` | `latitude` | `param` | `erosion_rate_fname` | CSV column holding the outlet latitude. | WGS84 decimal degrees. Not used for vector input. |
| `erosion_rate_longitude_column` | `string` | `longitude` | `param` | `erosion_rate_fname` | CSV column holding the outlet longitude. | WGS84 decimal degrees. Not used for vector input. |
| `erosion_rate_column` | `string` | `erosion_rate` | `param` | `erosion_rate_fname` | Column holding the basin-averaged erosion rate. | In the units of erosion_rate_units. A missing column is an error. |
| `erosion_rate_uncertainty_column` | `string` | `erosion_rate_uncertainty` | `param` | `erosion_rate_fname` | Column holding the 1-sigma uncertainty of the erosion rate. | Same units as the rate. A missing column is an error; a blank or non-positive value marks that sample missing_uncertainty. |
| `erosion_rate_units` | `string` | - | `param` | `erosion_rate_fname` | Units of the erosion rate and its uncertainty (required, no default): mm_per_kyr, m_per_myr, cm_per_kyr, mm_per_yr or m_per_yr. | mm_per_kyr equals m_per_myr (1e-6 m/yr). Mass fluxes (g cm^-2 yr^-1) are refused: divide by the rock density first. When the file has an erosion_rate_units column every row must agree with this parameter. |
| `erosion_rate_reported_area_column` | `string` | `reported_area_km2` | `param` | `erosion_rate_fname` | Column holding each sample's reported drainage area in km^2. | Used for area-matched snapping and the DEM-vs-reported area check. Empty, or a column not in the file, skips both. |
| `erosion_rate_snapping_method` | `string` | `reported_area` | `param` | `erosion_rate_fname` | How a sample is snapped to the channel network: reported_area, expanding_square, downstream, largest_area or none. | reported_area takes the nearest channel cell whose drainage area is within erosion_rate_area_mismatch_tolerance of the reported area, else the cell closest in area. Without a reported area it falls back to expanding_square (nearest channel cell). largest_area takes the cell with the most contributing pixels within erosion_rate_snap_max_distance_metres. |
| `erosion_rate_snap_max_distance_metres` | `float` | `300` | `param` | `erosion_rate_fname` | Search radius for snapping a sample, in metres. | Converted to pixels with the cell size and rounded up, so a parameter file stays valid across resolutions. |
| `erosion_rate_area_mismatch_tolerance` | `float` | `1.25` | `param` | `erosion_rate_fname` | Ratio tolerance (&gt; 1) between the DEM drainage area and the reported area. | A basin with A_dem / A_reported outside [1/tol, tol] is area_mismatch. Also the matching tolerance of reported_area snapping. |
| `erosion_rate_area_mismatch_action` | `string` | `exclude` | `param` | `erosion_rate_fname` | What an area_mismatch basin does: exclude (written, kept out of every fit) or flag (kept in the fits and flagged). | - |
| `erosion_rate_nested_basins` | `string` | `keep` | `param` | `erosion_rate_fname` | Nested sample basins: keep, outermost or innermost. | keep uses all basins and resamples nesting groups in the bootstrap; outermost keeps only basins nested in no other sample; innermost only basins containing no other sample. |

### Erodibility calibration

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `erosion_rate_ksn_statistic` | `string` | `reach_median` | `param` | `erosion_rate_fname` | Basin k_sn that drives the calibration: reach_median, chi_z_slope or segment_median. | reach_median is the median k_sn of the basin's qualifying reaches (split at confluences only), with a bootstrap interval. The others are written as diagnostics either way. |
| `erosion_rate_min_reach_pixels` | `int` | `10` | `param` | `erosion_rate_fname` | Minimum channel nodes for a basin reach to qualify. | Reaches run between confluences, sources and the outlet. |
| `erosion_rate_min_reaches_per_basin` | `int` | `5` | `param` | `erosion_rate_fname` | Minimum qualifying reaches for a basin to enter the calibration. | Basins with fewer are written with status too_few_reaches and enter no fit. |
| `erosion_rate_fit_n` | `bool` | `true` | `param` | `erosion_rate_fname` | Fit n and K together (York regression of ln E on ln k_sn). | Attempted only with enough independent basin groups and a wide enough k_sn range; K at each fixed n is always written. |
| `erosion_rate_min_groups_for_n_fit` | `int` | `8` | `param` | `erosion_rate_fname` | Minimum independent basin groups for fitting n. | Nested basins form one group. Below this only K at fixed n is reported. |
| `erosion_rate_min_ksn_range_ratio` | `float` | `2.0` | `param` | `erosion_rate_fname` | Minimum max/min k_sn ratio over the used basins for fitting n. | - |
| `erosion_rate_mono_lithology_min_area_fraction` | `float` | `0.8` | `param` | `erosion_rate_fname` | Minimum area fraction of a basin's dominant unit for the basin to count as mono-lithologic. | Needs lithology_fname. |
| `erosion_rate_min_unit_reaches_per_basin` | `int` | `3` | `param` | `erosion_rate_fname` | Minimum qualifying lithologic reaches of a unit in a basin for the unit-in-basin estimate. | Needs lithology_fname. |
| `erosion_rate_unit_k_methods` | `string` | `mono,unit_in_basin` | `param` | `erosion_rate_fname` | Per-unit K estimators, comma separated: mono, unit_in_basin. | Needs lithology_fname. |
| `erosion_rate_property_table_method` | `string` | `unit_in_basin` | `param` | `erosion_rate_fname` | Per-unit estimator written to the draft lithology property table. | Needs lithology_fname. |
| `erosion_rate_systematic_relative_uncertainty` | `float` | `0` | `param` | `erosion_rate_fname` | Relative systematic uncertainty of every erosion rate (e.g. the production-rate term), added to the K interval in quadrature. | 0 adds nothing. |
| `erosion_rate_bootstrap_replicates` | `int` | `1000` | `param` | `erosion_rate_fname` | Bootstrap replicates of the calibration intervals. | 0 = no intervals. |
| `erosion_rate_bootstrap_random_seed` | `int` | `5489` | `param` | `erosion_rate_fname` | Seed of the calibration bootstrap. | One stream per m/n, separate from the k_sn summary and relative-erodibility streams. |
| `write_erodibility_calibration_basin_polygons` | `bool` | `true` | `param` | `erosion_rate_fname` | Write the sample basins as polygons (write_vector_format). | For checking the snapping against the published basins in a GIS. |
| `write_erodibility_calibration_property_table` | `bool` | `true` | `param` | `erosion_rate_fname` | Write a draft lithology property table from the per-unit K. | Needs lithology_fname. |

### Interior nodata

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `interior_nodata_method` | `string` | `sink` | `param` | - | How nodata pixels inside the DEM are treated: sink or fill. | Interior nodata are nodata pixels not connected to the raster edge through other nodata (isolated voids, moulins, masked lakes). sink (the default, unchanged behaviour): they are outlets, so flow that reaches them ends there, as at the DEM edge; use this for real internal drainage (e.g. moulins on an ice sheet). fill: they are interpolated with nodata_fill_method_choice before any filling or routing, so flow passes across them; use this for data voids. Nodata connected to the edge (e.g. sea removed by min_elevation) is never filled and always acts as an open boundary. |
| `nodata_fill_method_choice` | `string` | `laplacian` | `param` | `interior_nodata_method=fill` | Method for filling interior nodata pixels. | laplacian: smooth Dirichlet interpolation (recommended). iterative_weighted_average: fast, well-suited to thin voids. Used when interior_nodata_method is fill. |
| `nodata_fill_max_iter` | `int` | `10000` | `param` | `interior_nodata_method=fill` | Maximum iterations for iterative nodata fill methods. | Used when interior_nodata_method is fill. |

<!-- END GENERATED CLI PARAMETERS -->
