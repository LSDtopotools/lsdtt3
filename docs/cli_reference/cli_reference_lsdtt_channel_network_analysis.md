# CLI reference: lsdtt-channel-network-analysis

Builds a D8 channel and junction network from a DEM and then runs one or more channel-network analyses: channel-link bearing, junction-angle geometry, Strahler/Tokunaga ordering, and junction drainage-area ratios.

At least one output must be requested or the CLI exits with an error.

## Usage

```
lsdtt-channel-network-analysis <path> <parameter_file>
lsdtt-channel-network-analysis <parameter_file>         # path defaults to ./
lsdtt-channel-network-analysis -h                       # print help and exit
lsdtt-channel-network-analysis -v                       # print version and exit
lsdtt-channel-network-analysis -c                       # print citation and exit
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
| 1 | Pit filling (skipped when `raster_is_filled: true`) |
| 2 | D8 flow routing |
| 3 | Source pixel extraction from contributing-area threshold |
| 4 | Junction / channel network construction |
| 5 | Junction selection for bearing / angle analyses |
| 6 | Channel-link bearing analysis |
| 7 | Junction-angle analysis |
| 8 | Junction drainage-area ratios, then Strahler / Tokunaga ordering and output |

Steps 6 and 7 operate on the **selected junction subset** (see junction selection below). Step 8 (junction area ratios, Strahler / Tokunaga) always operates on the **full extracted network**.

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input raster filename. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |

There is no `write_output_format` for this CLI — vector and tabular outputs each have their own format parameter.

---

## Threading and routing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` uses all available hardware threads. Applies to pit filling and flow-routing construction. |
| `boundary_conditions` | string | `oooo` | Four-character string specifying boundary conditions for flow routing, in north–east–south–west order. `o` = open (drains freely). `n` = no-flux. `p` = periodic. `b` = forced outlet. |

---

## Step 1: DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `raster_is_filled` | bool | `false` | If true, skip pit filling and use the input raster directly for flow routing. |
| `fill_min_slope` | float | `0.0001` | Minimum gradient imposed between adjacent filled cells. Only used when `raster_is_filled` is false. |

---

## Step 3: Source extraction

| Parameter | Type | Default | Description |
|---|---|---|---|
| `threshold_contributing_pixels` | int | `1000` | Minimum upstream contributing pixels required for a pixel to be a channel source. Reduce to extract denser networks. Must be > 0. |

---

## Junction selection for bearing / angle analyses

Exactly one selection mode must be active. If none is specified, `analyse_all_junctions` is used automatically. Setting more than one is a fatal error.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `analyse_all_junctions` | bool | `false` | Analyse every junction in the extracted network. This is the default when no selection mode is explicitly set. |
| `analyse_upstream_of_all_outlets` | bool | `false` | Analyse all junctions upstream of every outlet junction in the network. |
| `analyse_upstream_of_outlet_csv` | bool | `false` | Analyse all junctions upstream of outlet points supplied in a file. Requires `basin_outlet_fname`. |

### Outlet file and snapping (used with `analyse_upstream_of_outlet_csv`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `basin_outlet_fname` | string | `NULL` | Path to the outlet point file. Format is auto-detected: `.csv` files must have columns named `latitude` and `longitude` (EPSG:4326). `.fgb` and `.gpkg` files are read as point vector layers. Required when `analyse_upstream_of_outlet_csv` is true. |
| `basin_outlet_snapping_method` | string | `downstream` | How supplied outlet points are snapped to the channel network before junction lookup. `none` = use the pixel as-is. `expanding_square` = search outward in a square until a pixel meets `basin_outlet_contributing_pixels_threshold`. `downstream` = walk downstream until a pixel meets the threshold. |
| `basin_outlet_contributing_pixels_threshold` | int | `1000` | Contributing-pixel threshold used when snapping outlets to the channel network. Must be > 0. |
| `basin_outlet_snap_min_search_radius_pixels` | int | `0` | Minimum search radius in pixels for expanding-square snapping. Used only with `expanding_square`. |
| `basin_outlet_snap_max_search_radius_pixels` | int | `25` | Maximum search radius in pixels for expanding-square snapping. Must be ≥ `basin_outlet_snap_min_search_radius_pixels`. |

---

## Complete-basin filtering

Restricts the junction-based analyses to *complete* basins — basins whose entire
upslope area touches no nodata or DEM edge. Because the bearing, junction-angle,
and centre-of-mass analyses are stream-order based, edge-truncated basins (whose
true upstream network is clipped by the DEM boundary) would bias the statistics,
so they are dropped before analysis. The filter is applied after junction
selection and affects the `_channel_link_bearings`, `_junction_angles`, and
`_junction_com_bearings` outputs.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `only_complete_basins` | bool | `true` | Keep only junctions whose upslope basin is complete (no upslope pixel is influenced by nodata or the DEM edge). **Note:** the default is `true`, so output now excludes edge-truncated basins (a change from older versions). Set `false` to analyse every selected junction regardless of completeness. |
| `min_basin_pixels` | int | `0` | Drop junctions whose upslope basin has fewer than this many pixels. `0` (default) disables the size filter. Applied alongside `only_complete_basins`. |

---

## Link sampling (Steps 6 and 7)

Controls how channel link geometry is sampled when fitting a bearing direction.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `link_sampling_mode` | string | `full_link` | How link geometry is sampled. `full_link` = use all nodes in the compact downstream link (recommended). `fixed_length` = sample from the anchor end up to `link_sampling_fixed_length_metres`. `fixed_n_pixels` = sample from the anchor end up to `link_sampling_fixed_n_pixels` pixels. |
| `link_sampling_fixed_length_metres` | float | `1000.0` | Length threshold for `fixed_length` mode. |
| `link_sampling_fixed_n_pixels` | int | `10` | Pixel-count threshold for `fixed_n_pixels` mode. |
| `link_sampling_min_nodes` | int | `2` | Minimum number of sampled nodes required to fit a bearing. Links with fewer nodes are skipped. |
| `channel_link_bearing_sample_anchor` | string | `upstream_end` | Which end of the compact downstream link to anchor bearing sampling. `upstream_end` or `downstream_end`. |

---

## Step 6: Channel-link bearing analysis

One record is written per analysed compact downstream link. Attributes include: junction ID, bearing (degrees from north), link length, and contributing area at the junction.

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_channel_link_bearings` | bool | `false` | `<write_fname>_channel_link_bearings.<ext>` |
| `write_vector_format` | string | `flatgeobuf` | Record-output format: `flatgeobuf` / `fgb`, `gpkg`, or `csv`. |
| `write_channel_link_bearing_vectors` | bool | `false` | `<write_fname>_channel_link_bearing_vectors.<ext>` — short line features for visual QA in a GIS. CSV is not supported for vector outputs. |
| `qa_write_vector_format` | string | `same_as_write_vector_format` | Override for line-based QA layers: `flatgeobuf` / `fgb` or `gpkg`. |

---

## Step 7: Junction-angle analysis

One record is written per valid two-donor junction. Attributes include: junction ID, donor-1 bearing, donor-2 bearing, receiver bearing, junction angle (between the two donors), and asymmetry angle.

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_junction_angles` | bool | `false` | `<write_fname>_junction_angles.<ext>` |
| `write_junction_angle_bearing_vectors` | bool | `false` | `<write_fname>_junction_angle_bearing_vectors.<ext>` — three short line features per junction (donor 1, donor 2, receiver) for visual QA. CSV is not supported. |

---

## Step 8a: Junction drainage-area ratios

A **continuous, high-sample-count complement to order-binned Tokunaga analysis**. At every interior confluence two (or more) channels merge; the program records the ratio `r = A_trib / A_main` of the two largest incoming drainage areas (in (0, 1]) at *every* junction — far more samples than the few well-populated Strahler-order bins. The distribution of `r` encodes the network's area-accumulation self-similarity. This analysis uses the full extracted network.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_junction_area_ratio_data` | bool | `false` | Write per-junction drainage-area ratios. Output: `<write_fname>_junction_area_ratios.csv`. |
| `write_junction_area_ratio_summary` | bool | `false` | Write the distribution summary (mean/median `r`, tail exponent `lambda`, Horton `R_A`, model-dependent `c_equivalent`) with branch-bootstrap confidence intervals. **Also** writes the area-decade collapse table. Outputs: `<write_fname>_junction_area_ratio_summary.csv` and `<write_fname>_junction_area_ratio_decades.csv`. |
| `junction_area_ratio_bootstrap` | int | `1000` | Number of branch-bootstrap replicates for the area-ratio confidence intervals. Resampling is over non-mainstem tributary branches (the main stem is always retained), the correct unit of independence for nested junctions. `0` disables bootstrapping (CI columns left as nan). Used only when `write_junction_area_ratio_summary` is true. |
| `junction_area_ratio_seed` | int | `1` | RNG seed for the area-ratio branch bootstrap. Fix for reproducible confidence intervals. |

### Metric definitions

- **`r = A_trib / A_main`** — ratio of the smaller to the larger merging tributary drainage area at a confluence, in (0, 1]. Areas are D8 contributing-pixel counts.
- **`lambda` tail exponent** — `x = −ln(r)` is exponentially distributed when `r` has a power-law tail. `lambda_all` is fitted over all junctions; `lambda_tail` over the small-`r` tail (`x` above its median).
- **`R_A`** — Horton area ratio from a per-order regression of mean downstream area.
- **`c_equivalent`** — a **model-dependent** Tokunaga-equivalent `c ≈ R_A^lambda`; an estimate only, not a substitute for the order-binned Tokunaga `(a, c)`.
- **Branch bootstrap** — junctions are nested and autocorrelated, so the unit of independence is a tributary branch. The main stem is retained, non-mainstem branches are resampled, yielding honest confidence intervals.
- **Area-decade collapse** — summary statistics restricted to confluences whose downstream area falls in each log10 area decade. If the network is self-similar these rows collapse (similar `median_ratio` / `lambda_tail`); this is the self-similarity diagnostic.

**Caveats.** This measures area accumulation, **not** side-branching topology, so it does not uniquely recover Tokunaga `(a, c)`. It is sensitive to the channel-extraction threshold (which truncates the small-`r` tail) and to the drainage-area definition (here D8 contributing pixels).

### Output files

| File | Contents | Columns |
|---|---|---|
| `<write_fname>_junction_area_ratios.csv` | One row per interior confluence | `junction_id`, `x`, `y`, `area_main_cells`, `area_trib_cells`, `ratio`, `area_down_cells`, `order_main`, `order_trib`, `order_receiver`, `n_tributaries`, `branch_id` |
| `<write_fname>_junction_area_ratio_summary.csv` | Single-row distribution summary with branch-bootstrap CIs | `n_junctions`, `mean_ratio`, `median_ratio`, `lambda_all`, `lambda_tail`, `r_a`, `c_equivalent`, `n_bootstrap`, `n_branches`, `median_ratio_ci_low`, `median_ratio_ci_high`, `lambda_tail_ci_low`, `lambda_tail_ci_high`, `r_a_ci_low`, `r_a_ci_high`, `c_equivalent_ci_low`, `c_equivalent_ci_high` |
| `<write_fname>_junction_area_ratio_decades.csv` | Area-decade collapse table (one row per populated log10 area decade) | `area_low_cells`, `area_high_cells`, `n`, `median_ratio`, `lambda_tail` |

---

## Step 8b: Strahler / Tokunaga analysis

These analyses operate on the **full extracted network**, not the outlet-selected subset. A *Strahler channel* is a maximal downstream chain of compact links of the same Strahler order. For a Strahler channel of order `i` flowing into a receiving channel of order `j`, the **Tokunaga offset** is `k = j − i`.

### Strahler channels

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_strahler_channels` | bool | `false` | `<write_fname>_strahler_channels.<ext>` — one polyline per Strahler channel with order, length, and elevation attributes. |

### Tokunaga channels

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_tokunaga_channels` | bool | `false` | `<write_fname>_tokunaga_channels.<ext>` — same geometry as Strahler channels with additional Tokunaga attributes: donor order, receiver order, Tokunaga offset. |
| `write_tokunaga_channel_data` | bool | `false` | `<write_fname>_tokunaga_channel_data.csv` — table-only export of Tokunaga channel attributes with no geometry. |
| `write_tokunaga_summary_data` | bool | `false` | Tokunaga pair counts, offset statistics, and self-similarity fits (see below). |

### Tokunaga self-similarity fits

The classical Tokunaga law models the mean count of side-tributary channels per offset:

```
T_k = a · c^(k − 1)
```

so that `a` is the count of channels one order lower joining a channel (offset `k = 1`) and `c` is the geometric ratio by which side-branching abundance grows per additional offset. When `write_tokunaga_summary_data` is true the CLI writes three tables (using the prefix `<write_fname>_tokunaga_summary`):

| File | Contents | Columns |
|---|---|---|
| `<write_fname>_tokunaga_summary_pair_counts.csv` | Counts of observed order pairs | `donor_order`, `receiver_order`, `offset`, `count` |
| `<write_fname>_tokunaga_summary_offset_stats.csv` | Counts / mean `T_k` grouped by offset `k` | `offset`, `count`, `mean_t` |
| `<write_fname>_tokunaga_summary_fit.csv` | Fitted Tokunaga-law parameters, **two rows** | `method`, `valid`, `a`, `c`, `intercept`, `exponent`, `r_squared`, `n_offsets_used`, `rmse_log_space` |

The `_fit.csv` table now contains **two fits**, distinguished by the `method` column:

- **`ols`** — the original ordinary log-linear least squares: regress `log(T_k)` against `(k − 1)`. The slope is `log(c)` and the intercept is `log(a)`. Every offset bin is weighted equally in log space.
- **`poisson`** — a count-weighted Poisson GLM fit (log link, solved by iteratively reweighted least squares). Each offset `k` contributes its observed mean `T_k` with a prior weight equal to its pair count.

**Why count-weighting matters.** `T_k` is a count rate, and high-offset bins are sparsely populated (few channel pairs separated by several Strahler orders). In ordinary log-linear least squares these noisy, low-count high offsets carry the same weight as the well-populated low offsets and can dominate the fitted slope `c`. The Poisson GLM weights each offset by its count, so it is the statistically appropriate model and gives a more stable `(a, c)` for networks where the high-offset tail is thin. Both rows expose the same `(a, c)` interpretation; compare them to judge how sensitive the fit is to the sparse tail. (`rmse_log_space` is reported on the `ols` row only; it is nan for `poisson`.)

---

## General output options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `calculate_spatial_index` | bool | `true` | Build a spatial index when writing vector files. Speeds up GIS use but increases memory overhead. |

---

## Example parameter file: bearing and junction angles

```
read_path:  /data/dems/
read_fname: my_dem_fill
write_path: /data/output/
write_fname: my_dem_network

raster_is_filled:         true
threshold_contributing_pixels: 1000
boundary_conditions:      oooo

# Select junctions upstream of supplied basin outlets
analyse_upstream_of_outlet_csv:  true
basin_outlet_fname:              /data/outlets/my_outlets.csv
basin_outlet_snapping_method:    downstream
basin_outlet_contributing_pixels_threshold: 1000

# Link sampling
link_sampling_mode:     full_link
link_sampling_min_nodes: 2

# Outputs
write_channel_link_bearings:          true
write_vector_format:                  flatgeobuf
write_junction_angles:                true
```

## Example parameter file: Strahler and Tokunaga

```
read_path:  /data/dems/
read_fname: my_dem
write_path: /data/output/
write_fname: my_dem_strahler

raster_is_filled:         false
fill_min_slope:           0.0001
threshold_contributing_pixels: 500

# These analyses use the full network; no junction selection needed.
write_strahler_channels:        true
write_tokunaga_channels:         true
write_tokunaga_summary_data:     true     # writes the ols + poisson _fit.csv
```

## Example parameter file: junction drainage-area ratios

```
read_path:  /data/dems/
read_fname: my_dem
write_fname: my_dem_areas

raster_is_filled:         false
threshold_contributing_pixels: 500

write_junction_area_ratio_data:    true
write_junction_area_ratio_summary: true   # also writes the decades collapse table
junction_area_ratio_bootstrap:     1000
junction_area_ratio_seed:          1
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

### Threading / routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all available hardware threads. Applies to pit filling and flow-routing construction. |
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four-character string for north, east, south, west. o = open, n = no-flux, p = periodic, b = forced outlet. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |
| `fill_flat_resolution_choice` | `string` | `inline_gradient` | `param` | - | How flats are made drainable: inline_gradient, barnes2014b_gradient or barnes2014b_float_step. | inline_gradient (default, unchanged): the Priority-Flood fill imposes fill_min_slope during the flood. The barnes2014b options make flow converge across flats (Barnes et al. 2014, doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses fill_min_slope (m/m, must be &gt; 0) as the minimum gradient; flats are raised about 3-4x more, capped below their rims, and capped flats are reported. barnes2014b_float_step is equivalent to RichDEM ResolveFlats on float32 in its first pass: it ignores fill_min_slope and its slopes depend on elevation (of order ULP(z) / cell size). With either barnes2014b option the result is always routable: cells the first pass leaves without a lower neighbour are refilled and resolved again (up to 4 rounds), and any still left are raised by a final priority-flood pass by one float32 step each; the report gives the rounds and the cells changed after the first pass. |

### DEM preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope (m/m) imposed by pit filling or flat resolution. | inline_gradient (default): imposed during the Priority-Flood fill, so used only when raster_is_filled is false. barnes2014b_gradient: the minimum gradient across flats, also on a raster_is_filled input. barnes2014b_float_step: ignored. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Selection mode for bearing / junction-angle analyses

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `analyse_all_junctions` | `bool` | `false` | `param` | - | Analyse all junctions in the extracted junction network. | If no explicit analysis-selection mode is chosen, the CLI will default to this mode. |
| `analyse_upstream_of_all_outlets` | `bool` | `false` | `param` | - | Analyse all junctions upstream of all outlet junctions. | This expands the selected outlet junctions to their full upstream junction subnetworks using the junction-network s_vector. |
| `analyse_upstream_of_outlet_csv` | `bool` | `false` | `param` | - | Analyse all junctions upstream of outlet points supplied in a file. | The file (CSV, FlatGeobuf, or GeoPackage) is loaded, optionally snapped to the channel network, converted to junctions, and expanded to their full upstream junction subnetworks. Set basin_outlet_fname to the path of the file. |

### Basin outlet input / snapping

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_outlet_fname` | `string` | `NULL` | `param` | - | File of requested outlet points (CSV, FlatGeobuf, or GeoPackage). | Required when analyse_upstream_of_outlet_csv is true. Format is auto-detected from the file extension: *.csv is read as latitude and longitude columns (EPSG:4326); *.fgb or *.gpkg is read as a point vector layer. For CSV inputs columns must be named 'latitude' and 'longitude'. For vector inputs all features must be point geometry. |
| `basin_outlet_snapping_method` | `string` | `downstream` | `param` | - | Optional snapping of supplied outlet CSV points to channel pixels. | Allowed values: none, expanding_square, downstream. |
| `basin_outlet_contributing_pixels_threshold` | `int` | `1000` | `param` | - | Contributing-pixel threshold used for optional outlet snapping. | This is total upslope contributing pixels from the flow-routing object. |
| `basin_outlet_snap_min_search_radius_pixels` | `int` | `0` | `param` | - | Minimum expanding-square snapping radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square. |
| `basin_outlet_snap_max_search_radius_pixels` | `int` | `25` | `param` | - | Maximum expanding-square snapping radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square. |

### Output general

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `calculate_spatial_index` | `bool` | `true` | `io` | - | Whether or not to calculate a spatial index for vector writes. | Spatial indexes speed up GIS use, but increase memory overhead. |

### Link sampling for bearings / junction angles

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `link_sampling_mode` | `string` | `full_link` | `param` | - | How link geometry is sampled for bearing fitting. | Allowed values: full_link, fixed_length, fixed_n_pixels. full_link is the default and recommended mode. |
| `link_sampling_fixed_length_metres` | `float` | `1000.0` | `param` | - | Length threshold for fixed_length sampling mode. | Used only when link_sampling_mode is fixed_length. |
| `link_sampling_fixed_n_pixels` | `int` | `10` | `param` | - | Pixel-count threshold for fixed_n_pixels sampling mode. | Used only when link_sampling_mode is fixed_n_pixels. |
| `link_sampling_min_nodes` | `int` | `2` | `param` | - | Minimum sampled nodes required for fitting a link bearing. | If sampling returns fewer nodes than this, the bearing fit is skipped. |
| `channel_link_bearing_sample_anchor` | `string` | `upstream_end` | `param` | - | Which end of a compact downstream link to anchor link-bearing sampling. | Allowed values: upstream_end, downstream_end. For channel-link bearing output the upstream_end is a sensible default. |

### Bearing outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_channel_link_bearings` | `bool` | `false` | `switch` | - | Write channel-link bearing output. | Writes one record per analysed compact downstream link. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Format for record-oriented network outputs. | Allowed values: csv, gpkg, flatgeobuf or fgb. Line-based QA layers use qa_write_vector_format. |

### Junction-angle outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_junction_angles` | `bool` | `false` | `switch` | - | Write junction-angle output. | Writes one record per valid two-donor junction. |
| `write_channel_link_bearing_vectors` | `bool` | `false` | `switch` | - | Write short downstream-pointing visual bearing vectors for analysed links. | This writes one short line per analysed compact downstream link. It is intended as a visual QA layer rather than a full attribute export. |
| `write_junction_angle_bearing_vectors` | `bool` | `false` | `switch` | - | Write visual bearing vectors for junction-angle analysis. | This writes three short centred line features per analysed junction: donor_1, donor_2, and receiver. It is intended as a visual QA layer. |
| `qa_write_vector_format` | `string` | `same_as_write_vector_format` | `param` | - | Output format override for line-based QA layers. | Allowed values: gpkg, flatgeobuf, fgb. CSV is not supported. |

### Strahler outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_strahler_channels` | `bool` | `false` | `switch` | - | Write Strahler-channel output. | Writes one record per derived Strahler channel. |

### Tokunaga outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_tokunaga_channels` | `bool` | `false` | `switch` | - | Write Tokunaga-channel output. | Writes the same Strahler-channel geometry but with Tokunaga-focused attributes such as donor order, receiver order, and Tokunaga offset. |
| `write_tokunaga_channel_data` | `bool` | `false` | `switch` | - | Write Tokunaga channel attributes to CSV. | This is a table-only export with no geometry. |
| `write_tokunaga_summary_data` | `bool` | `false` | `switch` | - | Write Tokunaga pair counts, offset statistics, and fit summary to CSV. | This is a summary-of-network export, not a geometry export. |

### Junction area-ratio outputs (complement to Tokunaga)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_junction_area_ratio_data` | `bool` | `false` | `switch` | - | Write per-junction drainage-area ratios to CSV. | One row per interior confluence: the ratio r = A_trib/A_main of the two merging tributary drainage areas, plus their Strahler orders. A continuous, high-sample-count complement to order-binned Tokunaga analysis. Output: &lt;write_fname&gt;_junction_area_ratios.csv. |
| `write_area_scaling_summary` | `bool` | `false` | `switch` | - | Write the drainage-area exceedance exponent beta, P(A &gt; a) ~ a^-beta. | The standard description of how drainage area is partitioned across a network, complementing the junction area ratios: those describe how area splits at one confluence, this describes the resulting distribution. Literature values are ~0.43 for real networks and ~0.45 for Optimal Channel Networks (Rodriguez-Iturbe et al. 1992; Rigon et al. 1996). NOTE lsdtt3 currently measures ~0.54-0.58 on generated OCNs and the gap is unexplained -- read the reported a_min, n_decades and n_used before using the number, and see docs/developer_notes/testing_and_benchmarking.md 8b.1c. Outputs &lt;write_fname&gt;_area_scaling_summary.csv and _area_scaling_curve.csv. |
| `area_scaling_upper_exclude_decades` | `float` | `1.0` | `param` | `write_area_scaling_summary` | Decades of the largest drainage areas excluded before fitting beta. | The top of the distribution is the finite-size rollover -- no basin can exceed the domain -- and fitting it measures the truncation, not the network. Raise this on small grids, lower it on large ones. |
| `area_scaling_min_fit_decades` | `float` | `1.5` | `param` | `write_area_scaling_summary` | Minimum decades of drainage area the fitted range must span. | Cutoff selection on its own favours short ranges, because a short segment of any smooth curve looks like a power law. If this cannot be met the result is NaN, which is the correct answer when there is too little range to fit. |
| `write_junction_area_ratio_summary` | `bool` | `false` | `switch` | - | Write a junction area-ratio distribution summary to CSV. | Mean/median r, the exponential tail rate lambda of x=-ln(r), a Horton area ratio R_A, and a model-dependent Tokunaga-equivalent c=R_A^lambda, with branch-bootstrap confidence intervals. Also writes an area-decade collapse table (&lt;write_fname&gt;_junction_area_ratio_decades.csv). Output: &lt;write_fname&gt;_junction_area_ratio_summary.csv. |
| `junction_area_ratio_bootstrap` | `int` | `1000` | `param` | - | Number of branch-bootstrap replicates for area-ratio confidence intervals. | Resampling is over non-mainstem tributary branches (the main stem is always retained), the correct unit of independence for nested junctions. 0 disables bootstrapping (CI columns left as nan). |
| `junction_area_ratio_seed` | `int` | `1` | `param` | - | RNG seed for the area-ratio branch bootstrap. | Fix for reproducible confidence intervals. |

### Junction centre-of-mass bearing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_junction_centre_of_mass_bearings` | `bool` | `false` | `switch` | - | Write the basin centre-of-mass bearing for every junction. | For each junction the basin is defined by the pixel just before the next downstream junction (base level for the outlet junction); this writes the geographic bearing (degrees from North, clockwise) from that outlet to the basin centre of mass, plus the outlet and centre-of-mass coordinates and the basin pixel count. Output: &lt;write_fname&gt;_junction_com_bearings.csv. |

### Complete-basin filtering

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `only_complete_basins` | `bool` | `true` | `param` | - | Restrict junction/bearing/angle analyses to complete basins. | A basin is complete if no pixel in its upslope area is influenced by nodata or the DEM edge. Incomplete (edge-truncated) basins are dropped before computing channel-link bearings, junction angles, and centre-of-mass bearings. Needed for stream-order-based statistics. |
| `min_basin_pixels` | `int` | `0` | `param` | - | Drop junctions whose upslope basin has fewer than this many pixels. | 0 disables the size filter. Applied alongside only_complete_basins. |

<!-- END GENERATED CLI PARAMETERS -->
