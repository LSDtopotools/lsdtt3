# CLI reference: lsdtt-ridgeline-analysis

Analyses drainage divides and their elevation profiles (ridgelines). The tool builds D8 flow routing and a junction/channel network, selects basins, extracts ridgelines from the basin topology, then processes each ridgeline (elevation, detrending, snout segmentation) and writes point, statistics, and extrema CSVs. It can also fit a polynomial surface and sample its metrics along each ridgeline, including 1‑D power spectra and characteristic spacing.

## Usage

```
lsdtt-ridgeline-analysis <path> <parameter_file>
lsdtt-ridgeline-analysis <parameter_file>     # path defaults to ./
lsdtt-ridgeline-analysis -h                   # print help and exit
lsdtt-ridgeline-analysis -v                   # print version and exit
lsdtt-ridgeline-analysis -c                   # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Analysis modes

Set with `ridgeline_analysis_mode`. The mode controls how ridgelines are derived from the basin topology.

| Mode | Description | Basin selection |
|---|---|---|
| `basin_outline_ridgelines` *(default)* | Splits each selected basin outline into a left and a right ridgeline. The split point is chosen by `outline_splitting_method`; ridge snouts are segmented following Robinson et al. (2025). | Required (`select_all_basins`, `select_basins_from_outlets`, or `..._nested`). |
| `interfluvial_divide_ridgelines` | Extracts the divide shared between every pair of non‑nested selected basins. | Required. |
| `range_defined_ridgelines` | Detects range‑spanning basins draining to opposing open (`o`) boundaries (mountain fronts) and extracts transverse divides and the Main Drainage Divide. | Automatic (range‑spanning detection); needs 4 boundary segments with at least one opposing pair of open boundaries. |

To analyse **all ridges in a DEM**, use the default mode with `select_all_basins: true`.

---

## Pipeline overview

| Phase | Description |
|---|---|
| 1 | Open DEM; pit‑fill unless `raster_is_filled` |
| 2 | Build D8 flow routing and the junction/channel network |
| 3 | Select basins (from outlets, all, or range‑spanning detection) and apply filters |
| 4 | Optionally write basin outlines / rasters |
| 5 | Extract ridgelines according to the mode |
| 6 | Fill elevations, fit/detrend each profile, segment snouts, smooth |
| 7 | Write ridgeline point / stats / extrema CSVs |
| 9 | Optional polynomial surface fit, sample metrics along ridgelines, optional power spectra and characteristic spacing |

Nothing is written unless the relevant `print_*` / `write_*` switch is set.

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input raster filename (including extension, e.g. `my_dem.bil`). |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | input stem | Output filename stem. |

---

## Threading and flow routing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` uses all available hardware threads. |
| `boundary_conditions` | string | `oooo` | Four‑character N, E, S, W string. `o` = open, `p` = periodic. `range_defined_ridgelines` interprets open boundaries as mountain fronts. |
| `raster_is_filled` | bool | `false` | If false, `fill_pits()` is run before flow routing. |
| `fill_min_slope` | float | `0.0001` | Minimum slope imposed by pit filling (used when `raster_is_filled` is false). |
| `source_pixels_threshold` | int | `1000` | Contributing‑pixel threshold defining channel heads. Reduce for a denser network and more (smaller) basins. |

---

## Basin selection

Exactly one selection method must be enabled for `basin_outline_ridgelines` and `interfluvial_divide_ridgelines` (selection is automatic for `range_defined_ridgelines`).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `select_all_basins` | bool | `false` | Select one basin for every junction in the network. Use this to analyse all ridges. |
| `select_basins_from_outlets` | bool | `false` | Select basins from the outlet points in `basin_outlet_file`. |
| `select_basins_from_outlets_nested` | bool | `false` | As above, plus all nested upstream basins. |
| `basin_outlet_file` | string | `NULL` | Outlet points (CSV with `latitude`/`longitude` in EPSG:4326, or `.fgb`/`.gpkg` point layer). Required when selecting from outlets. |
| `basin_outlet_snapping_method` | string | `downstream` | Snap outlet points to a qualifying channel pixel: `none`, `expanding_square`, or `downstream`. |
| `basin_outlet_contributing_pixels_threshold` | int | `1000` | Contributing‑pixel threshold used for outlet snapping. |
| `basin_outlet_min_search_radius_pixels` | int | `0` | Minimum search radius (pixels) for `expanding_square` snapping. |
| `basin_outlet_max_search_radius_pixels` | int | `25` | Maximum search radius (pixels) for `expanding_square` snapping. |

### Basin filtering

| Parameter | Type | Default | Description |
|---|---|---|---|
| `filter_basins_by_stream_order` | bool | `false` | Keep only basins within a stream‑order range. |
| `basins_min_stream_order` | int | `-1` | Minimum stream order to keep (`-1` = ignore). |
| `basins_max_stream_order` | int | `-1` | Maximum stream order to keep (`-1` = ignore). |
| `filter_basins_by_area` | bool | `false` | Keep only basins within an area range in pixels. |
| `basins_min_area_pixels` | int | `-1` | Minimum basin area in pixels (`-1` = ignore). |
| `basins_max_area_pixels` | int | `-1` | Maximum basin area in pixels (`-1` = ignore). |
| `filter_incomplete_basins` | bool | `false` | Remove basins touching nodata or the DEM edge. |
| `filter_all_but_largest_basin` | bool | `false` | Keep only the largest basin. |
| `filter_nested_basins` | bool | `false` | Keep only the largest non‑nested basins. |

---

## Basin outputs

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_basin_outlines` | bool | `false` | `<write_fname>_junction_basin_outlines.<ext>` — basin outlines as polylines. |
| `write_basin_outlines_as_polygon` | bool | `false` | `<write_fname>_junction_basin_polygons.<ext>` — basin outlines as polygons. |
| `basin_outlines_output_format` | string | `flatgeobuf` | Format for basin outlines: `flatgeobuf`, `fgb`, or `gpkg`. |
| `write_basin_outlines_to_csv` | bool | `false` | `<write_fname>_basin_outlines.csv` — outline points with distance and elevation. |
| `write_basins_raster` | bool | `false` | `<write_fname>_junction_basins.tif` — basins labelled by junction number. |

---

## Ridgeline extraction and processing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `outline_splitting_method` | string | `longest_channel` | How a basin outline is split into two ridgelines (`basin_outline_ridgelines` only): `longest_channel` (point nearest the longest channel's head), `max_elevation`, or `equal_length`. |
| `ridgeline_smoothing_window_fraction` | float | `0.05` | LOWESS smoothing window, as a fraction (0–1) of total ridgeline length, used for max‑curvature snout segmentation (Robinson et al. 2025; Kneed algorithm, Satopää et al. 2011). |
| `ridgeline_pixels_smoothing_window` | int | `1` | Pixel window for the moving‑average smoothing of the detrended ridgeline before extracting extrema. `1` = no smoothing. Distinct from `ridgeline_smoothing_window_fraction`. |
| `use_polynomial_fit_for_ridgeline` | bool | `false` | Detrend ridgeline elevation with a polynomial fit. Slope/steepness are always from a linear fit; this only changes detrended‑elevation outputs. |
| `use_constant_fit_for_mdd_ridgeline` | bool | `false` | Fit `z = constant` through the Main Drainage Divide for detrending (overrides linear/polynomial). `range_defined_ridgelines` only. |
| `m_over_n` | float | `0.5` | Concavity index (θ) for normalised ridgeline steepness. Constrain with `lsdtt-chi-analysis` and supply the result here. |

---

## Ridgeline CSV outputs

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `print_ridgelines` | bool | `false` | `<write_fname><mode_prefix>.csv` — per‑point ridgeline data: lat/long, distance, elevation, fitted/detrended/smoothed elevation, smoothed detrended slope. |
| `print_ridgeline_stats` | bool | `false` | `<write_fname><mode_prefix>_stats.csv` — per‑ridgeline statistics: fitted slope, linear‑fit R², normalised steepness, polynomial R² and coefficients. |
| `print_ridgeline_extrema` | bool | `false` | `<write_fname><mode_prefix>_extrema.csv` — ridgeline extrema with distance/relief to the next extremum and distance to the triple divide. |

`<mode_prefix>` is `_basin_outline_ridgelines`, `_interfluvial_ridgelines`, or `_range_defined_ridgelines`.

---

## Polynomial surface fit and along‑ridgeline mapping (Phase 9)

Runs when `write_polyfit_rasters` is true **or** `polyfit_metrics` is non‑empty. Sampled metric values are appended to a `_with_polyfit` ridgeline CSV per method.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `polyfit_metrics` | string | _(empty)_ | Comma‑separated polyfit metrics to compute and sample (or `all`). E.g. `slope,tangential_curvature,plan_curvature`. |
| `polyfit_lambda_metres` | float | `70.0` | Characteristic lengthscale (m) for the quadratic fit window. Reduce for lidar (e.g. 10–30 m). |
| `polyfit_method` | string | `direct` | `direct` (rolling‑moment, no FFTW), `wavelet` (FFT, requires FFTW), or `both`. Outputs are prefixed `_direct_*` / `_wavelet_*`. |
| `polyfit_z_factor` | float | `1.0` | Vertical scaling applied before interpreting derivatives. |
| `polyfit_flat_slope_tolerance` | float | `1e-12` | Slope tolerance below which direction‑dependent metrics are nodata. |
| `write_polyfit_rasters` | bool | `false` | Write the six coefficient rasters and each requested metric raster (`_<method>_polyfit_*.tif`). |
| `write_ridgeline_power_spectra_csv` | bool | `false` | Write the 1‑D power spectrum of each sampled metric along each ridgeline (`_<method><mode_prefix>_power_spectra.csv`). Requires FFTW. |
| `write_ridgeline_autocorrelation_csv` | bool | `false` | Write the characteristic spacing (dominant wavelength) of each sampled metric per ridgeline (`_<method><mode_prefix>_characteristic_spacing.csv`). |

---

## Example: all ridges in a DEM (defaults)

```
read_path:  /data/IndianCreek_1m/
read_fname: indian_creek.bil
write_fname: indian_creek

# Default mode: split every basin outline into two ridgelines
ridgeline_analysis_mode: basin_outline_ridgelines
select_all_basins: true

# Channel network density
source_pixels_threshold: 1000

# Write the ridgeline CSVs
print_ridgelines:        true
print_ridgeline_stats:   true
print_ridgeline_extrema: true
```
