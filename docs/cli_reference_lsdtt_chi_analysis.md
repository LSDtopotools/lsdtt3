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
| 3. Source extraction | `source_pixels_threshold` |
| 4. Chi extraction + data map | `print_chi_data_maps: true` |
| 5. Chi network fitting + disorder sweep | `fit_chi_networks: true` |
| 6. Disorder uncertainty | `run_subset_sensitivity`, `run_jackknife`, `run_bootstrap` |
| 7. Spatial k_sn / u* inversion | `run_inversion: true` (requires Eigen3) |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input DEM filename. Any GDAL-supported format. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `write_output_format` | string | `cog` | Raster output format for inversion rasters. `cog`, `gtiff`, or `same_as_input`. |

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

If `basin_outlet_file` is `NULL` the largest basin in the DEM is used automatically.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `basin_outlet_file` | string | `NULL` | CSV, FlatGeobuf, or GeoPackage of outlet points. CSV must have `latitude` and `longitude` columns (EPSG:4326). |
| `basin_outlet_snapping_method` | string | `downstream` | `none`, `expanding_square`, or `downstream`. |
| `basin_outlet_contributing_pixels_threshold` | int | `1000` | Contributing-pixel count required at the snapped outlet. |
| `basin_outlet_min_search_radius_pixels` | int | `0` | Minimum search radius for `expanding_square` (pixels). |
| `basin_outlet_max_search_radius_pixels` | int | `25` | Maximum search radius for `expanding_square` (pixels). |

---

## Source extraction

| Parameter | Type | Default | Description |
|---|---|---|---|
| `source_pixels_threshold` | int | `1000` | Contributing-pixel count used to define channel heads. Chi is computed only on pixels upstream of channel heads. |

---

## Chi data maps (diagnostic output)

Runs before `fit_chi_networks`. Use to verify basin and network geometry before committing to a full m/n sweep.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `print_chi_data_maps` | bool | `false` | Extract chi networks and write per-node data maps. Outputs: nodes (chi, drainage_area, flow_distance, elevation, source_key, basin_key), valid basin outlets, source heads, and basin polygons. |
| `chi_data_maps_output_format` | string | `csv` | Output format for chi data map **point** files. `csv`: latitude/longitude in EPSG:4326. `gpkg`: GeoPackage. `flatgeobuf` / `fgb`: FlatGeobuf. Basin polygons are always vector (FlatGeobuf fallback when `csv` is requested). |
| `write_chi_data_map_lines` | bool | `false` | Write chi channels as polylines (one feature per `source_key`). Attributes: `source_key`, `basin_key`. |
| `chi_data_map_lines_output_format` | string | `flatgeobuf` | Output format for chi channel polylines. `flatgeobuf` / `fgb` or `gpkg`. CSV not supported for lines. |

---

## Chi network fitting and m/n sweep

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fit_chi_networks` | bool | `false` | Run the piecewise-linear chi-network fitting and disorder sweep over the m/n range. |
| `chi_A0` | float | `1.0` | Reference drainage area A₀ (m²) for chi integration. χ = ∫(A₀/A)^(m/n) ds from outlet. |
| `mn_min` | float | `0.3` | Minimum m/n ratio for the disorder sweep. |
| `mn_max` | float | `0.6` | Maximum m/n ratio for the disorder sweep. |
| `mn_step` | float | `0.02` | Step size for the m/n disorder sweep. Must be > 0. |
| `fit_critical_divergence_m` | float | `5.0` | Critical elevation divergence (m) used to split chi-elevation branches. |

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
| `chi_fit_segment_lines_output_format` | string | `flatgeobuf` | Output format for fitted chi segment polylines. `flatgeobuf` / `fgb` or `gpkg`. CSV not supported. |
| `write_disorder_sample_tables` | bool | `false` | Write detailed disorder sample tables for each m/n value. |
| `calculate_spatial_index` | bool | `false` | Calculate a spatial index for vector outputs. Requires buffering all features in memory; disable for large datasets. |

---

## Minimal example parameter file

```
# lsdtt-chi-analysis parameter file — chi data maps only
read_path:  /data/dems/
read_fname: catchment.tif

raster_is_filled:           false
source_pixels_threshold:    1000

print_chi_data_maps:        true
chi_data_maps_output_format: flatgeobuf
write_chi_data_map_lines:   true
```

```
# lsdtt-chi-analysis parameter file — m/n sweep with bootstrap
read_path:  /data/dems/
read_fname: catchment.tif

source_pixels_threshold:    1000
basin_outlet_file:          outlets.csv

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
