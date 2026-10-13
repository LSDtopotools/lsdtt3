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
| `range_defined_ridgelines` | Identifies range-spanning basins on both opposing open (`o`) boundaries, interpreted as mountain fronts, and constructs the transverse divides and Main Drainage Divide from those basin sequences. | Automatic (range-spanning detection); needs 4 boundary segments with at least one opposing pair of open boundaries. |

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
| 5b | Optional paired flow path statistics (lengthwise asymmetry) |
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
| `threshold_contributing_pixels` | int | `1000` | Contributing‑pixel threshold defining channel heads. Reduce for a denser network and more (smaller) basins. |

---

## Basin selection

Exactly one selection method must be enabled for `basin_outline_ridgelines` and `interfluvial_divide_ridgelines` (selection is automatic for `range_defined_ridgelines`).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `select_all_basins` | bool | `false` | Select one basin for every junction in the network. Use this to analyse all ridges. |
| `select_basins_from_outlets` | bool | `false` | Select basins from the outlet points in `basin_outlet_fname`. |
| `select_basins_from_outlets_nested` | bool | `false` | As above, plus all nested upstream basins. |
| `basin_outlet_fname` | string | `NULL` | Outlet points (CSV with `latitude`/`longitude` in EPSG:4326, or `.fgb`/`.gpkg` point layer). Required when selecting from outlets. |
| `basin_outlet_snapping_method` | string | `downstream` | Snap outlet points to a qualifying channel pixel: `none`, `expanding_square`, or `downstream`. |
| `basin_outlet_contributing_pixels_threshold` | int | `1000` | Contributing‑pixel threshold used for outlet snapping. |
| `basin_outlet_snap_min_search_radius_pixels` | int | `0` | Minimum search radius (pixels) for `expanding_square` snapping. |
| `basin_outlet_snap_max_search_radius_pixels` | int | `25` | Maximum search radius (pixels) for `expanding_square` snapping. |

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
| `write_vector_format` | string | `flatgeobuf` | Format for basin outlines: `flatgeobuf`, `fgb`, or `gpkg`. |
| `write_basin_outlines_to_csv` | bool | `false` | `<write_fname>_basin_outlines.csv` — outline points with distance and elevation. |
| `write_basins_raster` | bool | `false` | `<write_fname>_junction_basins.tif` — basins labelled by junction number. |

---

## Ridgeline extraction and processing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `outline_splitting_method` | string | `longest_channel` | How a basin outline is split (`basin_outline_ridgelines` only): `longest_channel`, `max_elevation`, or `equal_length`; `none` returns one ridgeline following the complete basin outline. |
| `ridgeline_smoothing_window_fraction` | float | `0.05` | LOWESS smoothing window, as a fraction (0–1) of total ridgeline length, used for max‑curvature snout segmentation (Robinson et al. 2025; Kneed algorithm, Satopää et al. 2011). |
| `ridgeline_pixels_smoothing_window` | int | `1` | Pixel window for the moving‑average smoothing of the detrended ridgeline before extracting extrema. `1` = no smoothing. Distinct from `ridgeline_smoothing_window_fraction`. |
| `use_polynomial_fit_for_ridgeline` | bool | `false` | Detrend ridgeline elevation with a polynomial fit. Slope/steepness are always from a linear fit; this only changes detrended‑elevation outputs. |
| `use_constant_fit_for_mdd_ridgeline` | bool | `false` | Fit `z = constant` through the Main Drainage Divide for detrending (overrides linear/polynomial). `range_defined_ridgelines` only. |
| `m_over_n` | float | `0.5` | Concavity index (θ) for normalised ridgeline steepness. Constrain with `lsdtt-chi-analysis` and supply the result here. |

---

## Ridgeline CSV outputs

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_ridgelines` | bool | `false` | Write `<write_fname><mode_prefix>.csv` point data and a matching FlatGeobuf line layer. |
| `write_ridgeline_stats` | bool | `false` | `<write_fname><mode_prefix>_stats.csv` — per‑ridgeline statistics: fitted slope, linear‑fit R², normalised steepness, polynomial R² and coefficients. |
| `write_ridgeline_extrema` | bool | `false` | `<write_fname><mode_prefix>_extrema.csv` — ridgeline extrema with distance/relief to the next extremum and distance to the triple divide. |

`<mode_prefix>` is `_basin_outline_ridgelines`, `_interfluvial_ridgelines`, or `_range_defined_ridgelines`.

The FlatGeobuf line layer is written to `<write_fname><mode_prefix>.fgb` and
contains ridgeline IDs, basin IDs, node counts, and total lengths.

When `write_slope_area_analysis` is enabled, ridgeline point CSVs and
polyfit-enriched node outputs also contain `elevation_nearest_channel_left` and
`elevation_nearest_channel_right`. These are the elevations of the first marked
D8 channel pixel reached by tracing from the pixels immediately on either side
of each outgoing divide edge. The final node has `NaN` values.

## Selected-ridgeline triple-divide outputs

These opt-in switches run after ridgeline selection and write one combined
point layer for all enabled triple-divide methods:

| Parameter | Description |
|---|---|
| `write_triple_divides_on_selected_ridgelines` | Write triple and quadruple divide nodes identified where three or more surrounding flow basins meet a selected ridgeline. |
| `write_zero_order_triple_divides_on_selected_ridglines` | Write the join and leave nodes where first-order source-basin outlines bound a selected ridgeline. The parameter spelling is retained exactly as the CLI name. |
| `write_topological_distance_triple_divides_on_selected_ridgelines` | Write selected-ridgeline nodes where either left or right flow-topological distance exceeds a supplied threshold. |
| `topological_distance_triple_divide_threshold` | Backward-compatible single threshold in metres. |
| `topological_distance_triple_divide_thresholds` | Comma-separated thresholds in metres. A node qualifies when `left > threshold OR right > threshold`; one output pair is written per threshold. |

The native surrounding-basin and zero-order methods are written to separate
outputs so they can be loaded and compared independently. Enabled topological
distance points are also written separately:

```text
<write_fname>_native_triple_divides.csv
<write_fname>_native_triple_divides.fgb
<write_fname>_zero_order_triple_divides.csv
<write_fname>_zero_order_triple_divides.fgb
<write_fname>_topological_distance_triple_divides_XXXm.csv
<write_fname>_topological_distance_triple_divides_XXXm.fgb
<write_fname>_topological_distance_paired_flowpaths_XXXm.fgb
```

Each non-empty layer uses the common schema containing ridgeline ID/type,
method, node index, distance along ridgeline, elevation, selected basin
junction IDs, and latitude/longitude in the CSV. The FlatGeobuf stores point
geometry in the DEM CRS. CSV files are written even when a selected method
finds no points. Threshold-specific topological FlatGeobufs are also created
when empty so threshold cases remain explicit.

For each topological-distance threshold, the paired-flowpath FlatGeobuf stores
the incoming and outgoing branches on every qualifying side of every TDJ. Its
attributes include ridgeline/node IDs, side, branch, threshold, pair distance,
and the common junction flat index.

Zero-order output uses the dedicated source-pixel basin view. Each source pixel
identified by `threshold_contributing_pixels` defines one zero-order basin from
the complete upslope contributing area. Its ordered divide-grid outline is
matched to each selected ridgeline; each contiguous contact interval contributes
its join and leave ridgeline nodes.

---

## Paired flow path statistics — lengthwise asymmetry (Phase 5b)

For the two flow paths that diverge across a divide and rejoin at a common
junction, the lengthwise asymmetry of Goren & Shelef (2024) is

```
delta_L   = 2 |L_i - L_j| / (L_i + L_j)
delta_chi = 2 |chi_i - chi_j| / (chi_i + chi_j)
```

where `L_i` is the along‑channel distance from where path `i` enters the channel
network down to the shared junction. Both are dimensionless and scale‑free.
Published medians run from **0.19** for sub‑parallel trunk streams to **1.00**
for tortuous networks.

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `write_paired_flowpath_statistics` | bool | `false` | Write `<write_fname>_paired_flowpaths.csv` (one row per pair) and `<write_fname>_paired_flowpaths_summary.csv` (medians plus the rejection counts). |
| `paired_flowpath_threshold_divide_distance` | float | `1000.0` | Reject divide edges within this distance, **along the divide**, of either end. See the warning below. |
| `paired_flowpath_unique_pairs` | bool | `false` | Keep one record per distinct channel pair. The default reproduces Goren & Shelef, who iterate every divide edge and so can count a pair twice. |
| `paired_flowpath_only_channel_heads` | bool | `false` | Keep only pairs whose nearest downstream channels are both channel heads. Goren & Shelef do **not** restrict to heads, so the default reproduces them. |

### Requires `interfluvial_divide_ridgelines`, and sibling sub‑basins

The two paths in a pair must **rejoin**, so the divide has to be internal — with
both flanks draining to a common junction. The CLI refuses to run in any other
mode: a basin *outline* separates a basin from its neighbours, so its flanking
pixels drain to different outlets and every pair would be discarded.

For the same reason, select **sibling sub‑basins**, not top‑level basins.
`filter_nested_basins` leaves only top‑level basins, whose neighbours drain to
different outlets — every pair is then rejected as "paths never rejoin", and the
output is empty. Use `select_all_basins` with `filter_basins_by_stream_order`
instead:

```
ridgeline_analysis_mode: interfluvial_divide_ridgelines
select_all_basins: true
filter_basins_by_stream_order: true
basins_min_stream_order: 2
basins_max_stream_order: 3
write_paired_flowpath_statistics: true
```

A large "paths never rejoin" count is normal even so — sibling sub‑basins in
different parents are correctly excluded. Read the rejection counts in the
summary CSV before reading the median: they are how you tell "no complexity"
from "no surviving pairs".

### Choosing the divide‑end threshold

The 1000 m default is Goren & Shelef's, and it is calibrated for real
landscapes. On small synthetic grids it can remove essentially everything: on a
400×200 grid of 30 m cells with stream order 2–3 basins, a 1000 m two‑ended
filter left 40 pairs out of 9755. Check `n_rejected_distance` in the summary
CSV, and lower or disable the threshold (`0`) when the divides are short.

Note this measurement is taken on the **raw** ridgeline geometry, before the
Phase 6 elevation fill and smoothing, which move the nodes.

---

## Polynomial surface fit and along-ridgeline mapping (Phase 9)

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

## Slope-area channel analysis

Slope-area analysis runs after basin selection and ridgeline extraction. It
uses a direct quadratic polynomial fit to the original DEM for topographic
slope, while both D8 and optional D-infinity routing use the filled DEM. The
maximum-median-slope logarithmic area bin supplies the channel-area threshold.
The binary channel mask always uses filled-D8 accumulation and sets pixels
outside the selected basin scope to zero. Nested basin selections are not
supported and produce a fatal error.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_slope_area_analysis` | bool | `false` | Run slope-area analysis and append nearest-channel elevation fields to node outputs. |
| `slope_area_flow_method` | string | `d8` | Area used for slope-area histograms: `d8` or `dinfinity`. Both route the filled DEM. |
| `slope_area_polyfit_lambda_metres` | float | `70.0` | Direct polynomial slope-fit lengthscale, independent of the general polyfit settings. |
| `slope_area_number_of_bins` | int | `100` | Number of evenly spaced bins in log10 drainage area. |
| `slope_area_clipping` | string | `all_basins` | Use `all_basins` for one union histogram/threshold or `per_basin` for one per selected basin. |
| `write_slope_area_raster` | bool | `true` | Write `<write_fname>_slope_area_channel_mask.tif`. Values are binary 0/1. |
| `write_slope_area_histogram` | bool | `true` | Write histogram and threshold CSV outputs. |

Outputs are:

```text
<write_fname>_slope_area_histogram.csv
<write_fname>_slope_area_thresholds.csv
<write_fname>_slope_area_channel_mask.tif
```

The histogram CSV records `d05`, `d25`, `d50`, `d75`, and `d95` slope
percentiles for each area bin and marks the selected maximum-median bin. The
threshold CSV records the selected area threshold for the union or each basin.

---

## Example: all ridges in a DEM (defaults)

```
read_path:  /data/IndianCreek_1m/
read_fname: indian_creek.bil
write_fname: indian_creek

# Default mode: split every basin outline into two ridgelines
ridgeline_analysis_mode: basin_outline_ridgelines
select_all_basins: true
filter_incomplete_basins:true
filter_nested_basins:true

# Channel network density
threshold_contributing_pixels: 1000

# Write the ridgeline CSVs
write_ridgelines:        true
write_ridgeline_stats:   true
write_ridgeline_extrema: true
```

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

### Threading / Routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all available hardware threads. |
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four-character string (N, E, S, W). o = open, p = periodic. Note that range_defined_ridgelines relies on open boundaries (o) being interpreted as mountain fronts. |
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope (m/m) imposed by pit filling or flat resolution. | inline_gradient (default): imposed during the Priority-Flood fill, so used only when raster_is_filled is false. barnes2014b_gradient: the minimum gradient across flats, also on a raster_is_filled input. barnes2014b_float_step: ignored. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |
| `fill_flat_resolution_choice` | `string` | `inline_gradient` | `param` | - | How flats are made drainable: inline_gradient, barnes2014b_gradient or barnes2014b_float_step. | inline_gradient (default, unchanged): the Priority-Flood fill imposes fill_min_slope during the flood. The barnes2014b options make flow converge across flats (Barnes et al. 2014, doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses fill_min_slope (m/m, must be &gt; 0) as the minimum gradient; flats are raised about 3-4x more, capped below their rims, and capped flats are reported. barnes2014b_float_step is equivalent to RichDEM ResolveFlats on float32 in its first pass: it ignores fill_min_slope and its slopes depend on elevation (of order ULP(z) / cell size). With either barnes2014b option the result is always routable: cells the first pass leaves without a lower neighbour are refilled and resolved again (up to 4 rounds), and any still left are raised by a final priority-flood pass by one float32 step each; the report gives the rounds and the cells changed after the first pass. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Workflow mode

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `ridgeline_analysis_mode` | `string` | `basin_outline_ridgelines` | `param` | `extract_only_mdd` | Main ridgeline analysis mode. | Allowed values: basin_outline_ridgelines, interfluvial_divide_ridgelines, range_defined_ridgelines. |
| `extract_only_mdd` | `bool` | `false` | `switch` | - | If true, skips the transverse divides and only extracts the main drainage divide. | Only applies when ridgeline_analysis_mode is range_defined_ridgelines. |

### Range-defined boundaries

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `corner_coordinates_fname` | `string` | `NULL` | `param` | - | CSV of DEM corner coordinates for boundary segmentation. | Required for range_defined_ridgelines if the DEM is not rectangular. |

### Basin Selection

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_outlet_fname` | `string` | `NULL` | `param` | - | File of requested basin outlet points (CSV, FlatGeobuf, or GeoPackage). | Required when select_basins_from_outlets or select_basins_from_outlets_nested is true. Format is auto-detected from the file extension: *.csv is read as a text file with latitude and longitude columns (WGS84 / EPSG:4326 decimal degrees); *.fgb or *.gpkg is read as a point vector layer in its own CRS, which is then transformed to the raster CRS. For CSV inputs the file must contain columns named exactly 'latitude' and 'longitude'. For vector inputs all features must be point geometry. |
| `select_basins_from_outlets` | `bool` | `false` | `switch` | - | Select basins from supplied basin outlet points. | Reads basin_outlet_fname, optionally snaps the mapped pixels, and then constructs one lsdtt_junction_basin per supplied basin outlet. |
| `select_basins_from_outlets_nested` | `bool` | `false` | `switch` | - | Select supplied basin outlets and all nested upstream basins. | Reads basin_outlet_fname, optionally snaps the mapped pixels, constructs one lsdtt_junction_basin per supplied basin outlet, and then expands selection to include all upstream nested basins from the junction-network s_vector. |
| `select_all_basins` | `bool` | `false` | `switch` | - | Select one basin for every junction in the junction network. | Uses the junction network to derive a basin outlet flat index for each junction and constructs one lsdtt_junction_basin per junction. |
| `select_sibling_basins` | `bool` | `false` | `switch` | - | Select basins as SIBLING sets: every junction with two or more direct donors becomes a parent, and its donor sub-basins are the selection. | This is the pair source for write_paired_flowpath_statistics. delta_L only exists for a divide whose two flanking basins share a parent, because otherwise the two flow paths never rejoin and the pair is discarded. Selecting by an area or stream-order window over all junctions leaves that to chance: on one snapping-sweep surface a 400-4000 pixel window selected 106 basins and every one of 10007 divide edges was rejected as non-rejoining, since no two of those basins shared a parent. Taking the donors of a parent directly makes the guarantee structural, and makes the interfluvial divide search linear in the number of parents rather than quadratic in the number of basins. It replaces the standard basin filters, which are SKIPPED when it is on, because dropping one member of a sibling group would remove the very pairing it exists to provide. |
| `sibling_parent_min_stream_order` | `int` | `0` | `param` | - | Minimum stream order of the PARENT junction in a sibling group. | 0 means no lower bound. This is the order of the parent, not of the donor sub-basins that make up the pair. |
| `sibling_parent_max_stream_order` | `int` | `0` | `param` | - | Maximum stream order of the PARENT junction in a sibling group. | 0 means no upper bound. |
| `sibling_donor_min_area_pixels` | `int` | `0` | `param` | - | Minimum size, in upslope pixels, of a donor sub-basin in a sibling group. | 0 means no lower bound. A donor outside the window is dropped from its group, and a group with fewer than two surviving donors is dropped whole, so the sibling guarantee is never weakened by this filter. |
| `sibling_donor_max_area_pixels` | `int` | `0` | `param` | - | Maximum size, in upslope pixels, of a donor sub-basin in a sibling group. | 0 means no upper bound. |
| `sibling_max_groups` | `int` | `0` | `param` | - | Keep at most this many sibling groups, the largest parent basins first. | 0 keeps every group. Groups are ranked by descending parent drainage area with ties broken by ascending junction number, truncated, and then restored to ascending junction order, so the cap is reproducible. It is how two treatments of the same landscape are compared at matched sample size, which the snapping sweep could not do: it produced 12671 pairs on one surface and 186 on another, and a median over 186 is not comparable with a median over 12671. Note what it matches directly: the number of PARENT BASINS sampled, not the number of pairs. Pairs per parent depend on how long the shared divide is, and snapped surfaces have much longer ones -- 200 groups gave 3761 pairs unsnapped and 39571 after Roering snapping on the same landscape. To match PAIR counts, bisect this value against n_pairs in the summary CSV: adding groups only ever adds pairs, so n_pairs is monotone in it (10 / 20 / 30 groups gave 1421 / 3384 / 5888 pairs), and the ranking is deterministic so the bisection is reproducible. |

### Basin Outline Snapping / Filtering

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_outlet_snapping_method` | `string` | `downstream` | `param` | - | Optional snapping of basin outlet CSV points to pixels that exceed a contributing-pixel threshold. | Allowed values: none, expanding_square, downstream. none: use the mapped raster pixel directly. expanding_square: search outward ring-by-ring until a qualifying pixel is found. downstream: follow the D8 receiver chain until a qualifying pixel is found. |
| `basin_outlet_contributing_pixels_threshold` | `int` | `1000` | `param` | - | Contributing-pixel threshold used for optional basin outlet snapping. | This is the total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. |
| `basin_outlet_snap_min_search_radius_pixels` | `int` | `0` | `param` | - | Minimum expanding-square snapping search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square. 0 means the mapped input pixel is tested first. |
| `basin_outlet_snap_max_search_radius_pixels` | `int` | `25` | `param` | - | Maximum expanding-square snapping search radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square. The search stops at the first radius containing one or more qualifying pixels. |
| `basins_min_stream_order` | `int` | `-1` | `param` | - | Minimum stream order to keep when filtering basins. | Inclusive. Defaults to -1, which means ignore minimum. |
| `basins_max_stream_order` | `int` | `-1` | `param` | - | Maximum stream order to keep when filtering basins. | Inclusive. Defaults to -1, which means ignore maximum. |
| `basins_min_area_pixels` | `int` | `-1` | `param` | - | Minimum basin area in pixels to keep when filtering basins by area. | Inclusive. Defaults to -1, which means ignore minimum. |
| `basins_max_area_pixels` | `int` | `-1` | `param` | - | Maximum basin area in pixels to keep when filtering basins by area. | Inclusive. Defaults to -1, which means ignore maximum. |
| `filter_basins_by_stream_order` | `bool` | `false` | `switch` | - | Apply the stream order filter to the selected basins. | Uses basins_min_stream_order and basins_max_stream_order, each inclusive, with -1 meaning ignore that bound. |
| `filter_basins_by_area` | `bool` | `false` | `switch` | - | Apply the area filter to the selected basins. | Uses basins_min_area_pixels and basins_max_area_pixels, each inclusive, with -1 meaning ignore that bound. |
| `filter_incomplete_basins` | `bool` | `false` | `param` | - | Remove all incomplete basins. | An incomplete basin has at least one DEM pixel that is adjacent to an invalid pixel. Off-grid neighbours are invalid, unless an n or p boundary is supplied. |
| `filter_all_but_largest_basin` | `bool` | `false` | `param` | - | Keep only the largest basin. | Applied later in the filtering order so it can be used after completeness or other filters. |
| `filter_nested_basins` | `bool` | `false` | `param` | - | Remove nested basins. | Only the largest non-nested basins in the current selection are retained. |

### Basin Output Prints

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_basin_outlines` | `bool` | `false` | `switch` | - | Write selected basin outlines as vector features. | Features are labelled by their junction number. |
| `write_basin_outlines_as_polygon` | `bool` | `false` | `switch` | - | Write basin outlines as polygons rather than polylines. | Set to true when you want basins printed as areas. |
| `write_vector_format` | `string` | `flatgeobuf` | `param` | - | Output format for basin outline visualisation. | Allowed values: flatgeobuf, fgb, gpkg. FlatGeobuf is the default. |
| `write_basin_outlines_to_csv` | `bool` | `false` | `switch` | - | Write basin outlines as a point CSV with distance and elevation. | Outputs a _basin_outlines.csv file containing points for all basins. |
| `write_basins_raster` | `bool` | `false` | `switch` | - | Write the selected basins as a labelled raster. | Basins are labelled by junction number. |

### Ridgeline Algorithms

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `outline_splitting_method` | `string` | `longest_channel` | `param` | - | Determines how basin outlines are split into two ridgelines: max_elevation, longest_channel or equal_length. | Only used if ridgeline_analysis_mode is basin_outline_ridgelines. |
| `ridgeline_smoothing_window_fraction` | `float` | `0.05` | `param` | - | When analysing the ridgelines along basin outlines, it is necessary to segment the ridge snouts. | When segmenting the snout based on max curvature (following Robinson et al. 2025), lowess smoothing is required first with a smoothing window as parameter. This defines the smoothing window as a fraction (0-1) of total ridgeline distance. |
| `ridgeline_pixels_smoothing_window` | `int` | `1` | `param` | - | The pixel window to smooth the detrended ridgeline with when extracting ridgeline extrema. Using a very simple moving average. | A pixel window of 1 (default) effectively means no smoothing happens. Note this is NOT the same as ridgeline_smoothing_window_fraction which is used for ridge-snout segmentation. |
| `use_polynomial_fit_for_ridgeline` | `bool` | `false` | `switch` | - | When set to true, detrends the ridgeline with a polynomial fit. | Ridgeline slope and steepness always calculated with linear fit, but this alters the detrended elevation outputs. |
| `use_constant_fit_for_mdd_ridgeline` | `bool` | `false` | `switch` | - | When set to true, it will fit a z = constant line through the main drainage divide for detrendeing. | Overrules the default linear fit and the optional polynomial fit. Only applicable if ridgeline_analysis_mode is range_defined_ridgelines. |
| `m_over_n` | `float` | `0.5` | `param` | - | The concavity index for calculating ridgeline steepness. Usually denoted as the Greek symbol theta. | Default is 0.5, can be constrained for the DEM using the lsdtt-chi-mapping program and supply the result here. |

### Ridgeline CSV Outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_ridgelines` | `bool` | `false` | `switch` | - | This triggers ridgeline printing to csv. | NOTE for any ridgeline to print, you must select an analysis mode. |
| `write_ridgeline_stats` | `bool` | `false` | `switch` | - | This triggers printing of ridgeline statistics (including ridgeline steepness) to csv. | Includes slope, R2, normalized steepness, and polynomial coeffs. |
| `write_ridgeline_extrema` | `bool` | `false` | `switch` | - | This triggers printing of ridgeline extrema to csv. | Includes distance to next extremum, relief, and distance to TDJ. |
| `write_ridgeline_flow_topological_distance` | `bool` | `false` | `switch` | - | Write the left and right flow topological distances for internal ridgeline nodes. | Outputs a csv capturing the topological flow distance through a common junction between adjacent DEM pixels on each side of the divide node. |

### Triple divide outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_triple_divides_on_selected_ridgelines` | `bool` | `false` | `switch` | - | Write triple and quadruple divide nodes on selected ridgelines. | The surrounding-basin test is used. |
| `write_zero_order_triple_divides_on_selected_ridglines` | `bool` | `false` | `switch` | - | Write triple divide join/leave nodes bounding zero-order basins. | Uses the existing basin-outline traversal; zero-order means a first-order source basin in the Strahler network. |
| `write_zero_order_basins_touching_selected_ridgelines` | `bool` | `false` | `switch` | - | Write source-based zero-order basin outlines that touch selected ridgelines. | Uses the same exact ridgeline-grid contact test as zero-order TDJ detection. |
| `write_topological_distance_triple_divides_on_selected_ridgelines` | `bool` | `false` | `switch` | - | Write nodes where either flanking flow topological distance exceeds the threshold. | Results are written with the other selected triple-divide modes in one metadata-compatible layer. |
| `topological_distance_triple_divide_threshold` | `float` | `0.0` | `param` | - | Flow topological distance threshold for topological triple divides. | A node qualifies when the left OR right distance is greater than this value. |
| `topological_distance_triple_divide_thresholds` | `string` | - | `param` | - | Comma-separated flow topological distance thresholds in metres. | When supplied, one topological triple-divide CSV and FlatGeobuf is written for each value. An empty value falls back to the singular topological_distance_triple_divide_threshold parameter. |
| `native_tdj_channel_thresholds` | `string` | - | `param` | - | Comma-separated channel-area thresholds in square metres for native TDJs. | Native TDJs are recalculated on the selected ridgelines for each value. |

### Paired flow paths

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_paired_flowpath_statistics` | `bool` | `false` | `switch` | - | Write the lengthwise asymmetry (delta_L) of paired flow paths across each divide. | For the two flow paths that diverge across a divide and rejoin at a common junction, delta_L = 2\|L_i - L_j\|/(L_i + L_j) measures planform complexity (Goren and Shelef 2024). Published medians run from 0.19 for sub-parallel trunk streams to 1.00 for tortuous networks. Requires ridgeline_analysis_mode = interfluvial_divide_ridgelines, because the two paths must rejoin: a basin outline separates a basin from its neighbours, so its flanking pixels are in different basins and every pair would be discarded. |
| `paired_flowpath_threshold_divide_distance` | `float` | `1000.0` | `param` | - | Reject divide edges within this distance (metres) ALONG THE DIVIDE of either end. | The ends of an interfluvial divide are where it meets a junction, and pairs there are trivially symmetric because the two paths have barely separated before rejoining. Including them drags the median down. Goren and Shelef use 1000 m. Set to 0 to disable the filter. |
| `paired_flowpath_unique_pairs` | `bool` | `false` | `switch` | - | Keep only one record per distinct channel pair. | Goren and Shelef iterate every divide edge and take data from the nearest downstream channel, so two edges can yield the same channel pair and be counted twice. The default (false) reproduces them; true de-duplicates. |
| `paired_flowpath_only_channel_heads` | `bool` | `false` | `switch` | - | Keep only pairs where both nearest downstream channels are channel heads. | Goren and Shelef do NOT restrict pairs to channel heads - a paired path may start further down a channel - so the default reproduces them. |
| `paired_flowpath_straight_link` | `bool` | `false` | `switch` | - | Also compute the straight-link variant of delta_L, which ignores channel sinuosity. | The straight-link distance accumulates the straight-line (Pythagorean) distance between successive junctions instead of the true flow distance, so delta_L_straight_link separates genuine planform asymmetry from channel sinuosity. Adds straight_link_a, straight_link_b and delta_L_straight_link to the paired flowpath csv and the medians to the summary csv. THIS IS NOT the Goren and Shelef (2024) statistic - it is an extension - so it must NOT be compared with their published medians of 0.19 and 1.00, which are for the flow-distance delta_L only. |

### Polynomial surface fit

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_polyfit_rasters` | `bool` | `false` | `switch` | - | Write the six quadratic polynomial fit coefficient rasters. | Computes a, b, c, d, e, f where z = ax^2 + by^2 + cxy + dx + ey + f over a moving window controlled by polyfit_lambda. |
| `polyfit_lambda_metres` | `float` | `70.0` | `param` | - | Characteristic lengthscale (metres) for the wavelet quadratic fit. | Controls the size of the moving window. Default (70 m) is tuned for Copernicus 30 m data. Reduce for lidar (e.g. 10-30 m). |
| `polyfit_lambda_metres_list` | `string` | - | `param` | - | Comma-separated polyfit lengthscales in metres. | When supplied, Phase 9 runs each lambda sequentially after the ridgelines are extracted, releasing each polyfit result before the next lambda is calculated. |
| `polyfit_metrics` | `string` | - | `param` | - | Comma-separated list of derived polyfit metrics to write. | Examples: slope,aspect_degrees,plan_curvature,profile_curvature. Use all to write every supported metric. Available metrics: fitted_surface, dzdx, dzdy, d2zdx2, d2zdxdy, d2zdy2, slope, slope_degrees,aspect_radians, aspect_degrees, laplacian, hessian_determinant, mean_curvature, gaussian_curvature, maximum_curvature, minimum_curvature, curvedness, shape_index, profile_curvature, plan_curvature, tangential_curvature, difference_curvature, horizontal_excess_curvature, vertical_excess_curvature, accumulation_curvature, ring_curvature, rotor. |
| `polyfit_z_factor` | `float` | `1.0` | `param` | - | Vertical scaling applied before interpreting polyfit derivatives. | Scales first and second derivatives consistently in the polyfit metric engine. |
| `polyfit_flat_slope_tolerance` | `float` | `1e-12` | `param` | - | Slope tolerance below which direction-dependent polyfit metrics are nodata. | Applies to aspect, plan curvature, tangential curvature, profile curvature, and rotor. |
| `polyfit_method` | `string` | `wavelet` | `param` | - | Algorithm used for the quadratic polynomial fit. | wavelet: FFT-based Gaussian wavelet fit (requires FFTW). direct: separable rolling-moment O(n) fit (no FFTW needed). both: run both methods and map outputs from each. Outputs are prefixed _wavelet_* or _direct_*. |
| `calculate_polyfit_metrics_directly_on_dem` | `bool` | `false` | `switch` | - | If true, skips the large-scale polynomial surface fit and calculates polyfit metrics directly on the raw DEM. | Internally forces a 3x3 direct Evans (1980) quadratic fit, overriding polyfit_method and polyfit_lambda_metres. Also supports periodic boundaries. |

<!-- END GENERATED CLI PARAMETERS -->
