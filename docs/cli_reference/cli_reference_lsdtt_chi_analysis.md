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
| `basin_outlet_snapping_method` | string | `downstream` | `none`, `expanding_square`, or `downstream`. |
| `basin_outlet_contributing_pixels_threshold` | int | `1000` | Contributing-pixel count required at the snapped outlet. |
| `basin_outlet_snap_min_search_radius_pixels` | int | `0` | Minimum search radius for `expanding_square` (pixels). |
| `basin_outlet_snap_max_search_radius_pixels` | int | `25` | Maximum search radius for `expanding_square` (pixels). |

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

Mudd, S.M., et al. (2018). A statistical framework to quantify spatial variation in channel gradients using the integral method of channel profile analysis. *Journal of Geophysical Research: Earth Surface*, 123(2), 344–364.

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

### Threading / routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all available hardware threads. |
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four-character string: north, east, south, west. o = open, n = no-flux, p = periodic, b = forced outlet. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |
| `fill_flat_resolution_choice` | `string` | `inline_gradient` | `param` | - | How flats are made drainable: inline_gradient, barnes2014b_gradient or barnes2014b_float_step. | inline_gradient (default, unchanged): the Priority-Flood fill imposes fill_min_slope during the flood. The barnes2014b options make flow converge across flats (Barnes et al. 2014, doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses fill_min_slope (m/m, must be &gt; 0) as the minimum gradient; flats are raised about 3-4x more, capped below their rims, and capped flats are reported. barnes2014b_float_step is equivalent to RichDEM ResolveFlats on float32: it ignores fill_min_slope, its slopes depend on elevation (of order ULP(z) / cell size), and rim violations are reported. |

### DEM preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope (m/m) imposed by pit filling or flat resolution. | inline_gradient (default): imposed during the Priority-Flood fill, so used only when raster_is_filled is false. barnes2014b_gradient: the minimum gradient across flats, also on a raster_is_filled input. barnes2014b_float_step: ignored. |

### Basin outlet input

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_outlet_fname` | `string` | `NULL` | `param` | - | File of basin outlet points (CSV, FlatGeobuf, or GeoPackage). | If NULL (default), the largest basin outlet in the DEM is used. CSV files must have 'latitude' and 'longitude' columns (EPSG:4326). |
| `basin_outlet_snapping_method` | `string` | `downstream` | `param` | - | Optional snapping of supplied outlet points to qualifying channel pixels. | Allowed values: none, expanding_square, downstream. |
| `basin_outlet_contributing_pixels_threshold` | `int` | `1000` | `param` | - | Contributing-pixel threshold used for optional outlet snapping. | Total upslope contributing pixels required at the snapped cell. |
| `basin_outlet_snap_min_search_radius_pixels` | `int` | `0` | `param` | - | Minimum expanding-square snapping search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square. |
| `basin_outlet_snap_max_search_radius_pixels` | `int` | `25` | `param` | - | Maximum expanding-square snapping search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square. |

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
| `write_ksn_basin_summary` | `bool` | `true` | `switch` | - | Write a per-basin channel-steepness (k_sn) summary CSV. | One row per fitted basin with the number of segments, number of main-stem segments, total chi length, and mean / median / min / max / chi-weighted-mean k_sn plus main-stem mean and median k_sn. Written once per m/n value tested. Output: &lt;write_fname&gt;_chi_ksn_basin_summary_mn_*.csv. |
| `write_chi_fit_segment_lines` | `bool` | `false` | `switch` | - | Write fitted chi segments as a vector polyline file. | One polyline feature per piecewise-linear segment. Attributes: basin_index, outlet_flat, segment_id, branch_id, k_sn (dz/dchi), downstream_chi, upstream_chi, downstream_elevation, upstream_elevation, n_nodes. Written once per m/n value tested. |
| `calculate_spatial_index` | `bool` | `false` | `io` | - | Calculate a spatial index for vector outputs. | Unlike the smaller vector-export CLIs, chi analysis leaves this off by default because its per-m/n outputs can be large. Enabling it speeds up later GIS queries but buffers each layer in memory before writing. |
| `write_disorder_sample_tables` | `bool` | `false` | `switch` | - | Write detailed disorder sample tables for each m/n value. | - |
| `write_output_format` | `string` | `cog` | `io` | - | Raster output format for the chi coordinate and inversion rasters. | Allowed values: cog, gtiff, same_as_input. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). CSV is available only through selectors for point/table outputs. Per-output overrides exist only where a container genuinely differs. Raster outputs are governed separately by write_output_format. |

<!-- END GENERATED CLI PARAMETERS -->
