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
| `analyse_upstream_of_outlet_csv` | bool | `false` | Analyse all junctions upstream of outlet points supplied in a file. Requires `basin_outlet_file`. |

### Outlet file and snapping (used with `analyse_upstream_of_outlet_csv`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `basin_outlet_file` | string | `NULL` | Path to the outlet point file. Format is auto-detected: `.csv` files must have columns named `latitude` and `longitude` (EPSG:4326). `.fgb` and `.gpkg` files are read as point vector layers. Required when `analyse_upstream_of_outlet_csv` is true. |
| `basin_outlet_snapping_method` | string | `downstream` | How supplied outlet points are snapped to the channel network before junction lookup. `none` = use the pixel as-is. `expanding_square` = search outward in a square until a pixel meets `basin_outlet_contributing_pixels_threshold`. `downstream` = walk downstream until a pixel meets the threshold. |
| `basin_outlet_contributing_pixels_threshold` | int | `1000` | Contributing-pixel threshold used when snapping outlets to the channel network. Must be > 0. |
| `basin_outlet_min_search_radius_pixels` | int | `0` | Minimum search radius in pixels for expanding-square snapping. Used only with `expanding_square`. |
| `basin_outlet_max_search_radius_pixels` | int | `25` | Maximum search radius in pixels for expanding-square snapping. Must be ≥ `basin_outlet_min_search_radius_pixels`. |

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
| `minimum_basin_pixels` | int | `0` | Drop junctions whose upslope basin has fewer than this many pixels. `0` (default) disables the size filter. Applied alongside `only_complete_basins`. |

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
| `channel_link_bearings_output_format` | string | `flatgeobuf` | Output format: `flatgeobuf` / `fgb`, `gpkg`, or `csv`. |
| `write_channel_link_bearing_vectors` | bool | `false` | `<write_fname>_channel_link_bearing_vectors.<ext>` — short line features for visual QA in a GIS. CSV is not supported for vector outputs. |
| `channel_link_bearing_vectors_output_format` | string | `flatgeobuf` | Output format for bearing vector QA layer: `flatgeobuf` / `fgb` or `gpkg`. |

---

## Step 7: Junction-angle analysis

One record is written per valid two-donor junction. Attributes include: junction ID, donor-1 bearing, donor-2 bearing, receiver bearing, junction angle (between the two donors), and asymmetry angle.

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_junction_angles` | bool | `false` | `<write_fname>_junction_angles.<ext>` |
| `junction_angles_output_format` | string | `flatgeobuf` | Output format: `flatgeobuf` / `fgb`, `gpkg`, or `csv`. |
| `write_junction_angle_bearing_vectors` | bool | `false` | `<write_fname>_junction_angle_bearing_vectors.<ext>` — three short line features per junction (donor 1, donor 2, receiver) for visual QA. CSV is not supported. |
| `junction_angle_bearing_vectors_output_format` | string | `flatgeobuf` | Output format for junction-angle QA vectors: `flatgeobuf` / `fgb` or `gpkg`. |

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
| `strahler_channels_output_format` | string | `flatgeobuf` | Output format: `flatgeobuf` / `fgb`, `gpkg`, or `csv`. |

### Tokunaga channels

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_tokunaga_channels` | bool | `false` | `<write_fname>_tokunaga_channels.<ext>` — same geometry as Strahler channels with additional Tokunaga attributes: donor order, receiver order, Tokunaga offset. |
| `tokunaga_channels_output_format` | string | `flatgeobuf` | Output format: `flatgeobuf` / `fgb`, `gpkg`, or `csv`. |
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
basin_outlet_file:               /data/outlets/my_outlets.csv
basin_outlet_snapping_method:    downstream
basin_outlet_contributing_pixels_threshold: 1000

# Link sampling
link_sampling_mode:     full_link
link_sampling_min_nodes: 2

# Outputs
write_channel_link_bearings:          true
channel_link_bearings_output_format:  flatgeobuf
write_junction_angles:                true
junction_angles_output_format:        flatgeobuf
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
strahler_channels_output_format: flatgeobuf
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
