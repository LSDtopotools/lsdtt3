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
| 7b | Bearing and junction-angle statistics by group (optional, `write_channel_bearing_statistics` / `write_junction_angle_statistics`) |
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
| `basin_outlet_snapping_method` | string | `downstream` | How supplied outlet points are snapped to the channel network before junction lookup. `none` = use the pixel as-is. `expanding_square` = search outward in a square until a pixel meets `basin_outlet_contributing_pixels_threshold`. `downstream` = walk downstream until a pixel meets the threshold. `largest_area` takes the pixel with the most contributing pixels in the square of half-width `basin_outlet_snap_max_search_radius_pixels`, if it reaches the threshold (a wide square can jump onto a bigger neighbouring river). |
| `basin_outlet_contributing_pixels_threshold` | int | `1000` | Contributing-pixel threshold used when snapping outlets to the channel network. Must be > 0. |
| `basin_outlet_snap_min_search_radius_pixels` | int | `0` | Minimum search radius in pixels for expanding-square snapping. Used only with `expanding_square`. |
| `basin_outlet_snap_max_search_radius_pixels` | int | `25` | Maximum search radius in pixels for expanding-square snapping, and the square's half-width for `largest_area`. Must be ≥ `basin_outlet_snap_min_search_radius_pixels`. |

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

The whole network can be written too, so maps can show analysed links and
junctions against the rest:

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_channel_network_lines` | bool | `false` | `<write_fname>_channel_network_lines.<ext>` — every link of the network as a line (`junction_upstream, junction_downstream, stream_order, length, drainage_area`), in `qa_write_vector_format`. A link was analysed when its `junction_upstream` is an `up_junction` of the link-bearing outputs. |
| `write_junction_network_points` | bool | `false` | `<write_fname>_junction_network_points.<ext>` — every junction as a point (`junction_number, receiver_junction, junction_stream_order`), in `write_vector_format`. A junction was analysed when its `junction_number` is in the `junction` column of `_junction_angles`. |

---

## Step 7: Junction-angle analysis

One record is written per valid two-donor junction. Attributes include: junction ID, donor-1 bearing, donor-2 bearing, receiver bearing, junction angle (between the two donors), and asymmetry angle.

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_junction_angles` | bool | `false` | `<write_fname>_junction_angles.<ext>` |
| `write_junction_angle_bearing_vectors` | bool | `false` | `<write_fname>_junction_angle_bearing_vectors.<ext>` — three short line features per junction (donor 1, donor 2, receiver) for visual QA. CSV is not supported. |

---

## Step 7b: Bearing and junction-angle statistics by group

Steps 6 and 7 write one record per channel link and per junction. Step 7b asks
whether their angles differ between groups: between Strahler orders, between
basins or regions, or between lithologic units. It summarises each group's
distribution, compares the groups pairwise (or one against the rest), and
writes the comparisons with two kinds of p-value. Two switches turn it on:
`write_channel_bearing_statistics` and `write_junction_angle_statistics`. They
add output files; the record files are unchanged except for the lithology
columns listed below.

The statistical model and its tolerances are in
[Channel-link bearing statistics](../design_and_algorithms/channel_bearing_statistics.md)
and [Lithology and channel steepness statistics](../design_and_algorithms/lithology_channel_statistics.md).

### Angular convention

- **Bearings** are in degrees, measured clockwise from grid north (the projected
  CRS's y axis). A link's bearing is the direction from its upstream junction to
  its downstream junction; the **directional** form uses this downstream bearing
  (period 360).
- The **axial** form folds a bearing and its reverse onto one value
  (theta and theta + 180, period 180), so its mean lies in [0, 180). Use it when
  the direction of flow along a link is not the quantity of interest.
  `channel_bearing_angular_type` chooses `directional`, `axial` or `both` (default).
- **Junction angles** are not directions. `alpha` is the donor opening angle
  (`donor_opening_angle_deg`), and the sector magnitudes `beta_min` and `beta_max`
  are the smaller and the larger donor-receiver angle. The beta values do not
  depend on which donor is labelled 1 or 2. They are compared with linear
  statistics (KS and linear Wasserstein), not circular ones.
- **D8 quantisation.** The flow direction is on a grid, so a link's bearing is
  quantised: a two-node link is exactly 0, 45, 90, ... degrees. The axial
  mean resultant length cancels the four-fold D8 pattern, so the axial statistics
  are the more robust of the two to D8 quantisation. The summaries report
  `fraction_on_d8_multiples` (within 0.5 degrees of a multiple of 45) so the
  quantisation can be seen directly.

### Grouping and stratifying keys

| Key | Applies to | Groups |
|---|---|---|
| `strahler_order` | bearings and junctions | the link's own order (bearings) or the order at the confluence pixel (junctions). `channel_bearing_max_order_group` = k pools every order >= k into one group labelled `k+` (default 0 = no pooling). High orders are small, so pooling is usually needed. |
| `donor_order_pair` | junctions only | the two donor orders, as `i_j`. Refused for bearings. |
| `basin` | bearings and junctions | the basin holding the record. `channel_bearing_basin_order` = 0 uses the analysis-selection outlets (nested outlets resolve so the innermost wins). k > 0 uses the Strahler order-k basin; records of order above k have no such basin and are excluded and counted. |
| `region` | bearings and junctions | the region in `channel_region_fname` (a polygon field or a code raster). A junction takes the region of its pixel; a link takes the length-weighted majority region of the nodes its bearing was fitted on, as for lithology, with its purity. |
| `lithology` | bearings and junctions | the unit of the record, with `lithology_fname` (below). |
| `lithology_interior` | junctions only | the unit of the junction's pixel, but only when all three sampled arms are in it. |

`*_stratify_by` takes a second key whose levels are analysed separately (or
`none`); it must differ from the group key. Groups are compared only within a
level of the stratum. A pooled level `all` is always written.

Keys for bearings are set with `channel_bearing_group_by` and
`channel_bearing_stratify_by`; keys for junctions with `junction_angle_group_by` and
`junction_angle_stratify_by`.

### Lithology grouping and the "mixed" and "arms disagree" groups

With `lithology_fname`, each link takes the **length-weighted majority unit** of the
nodes its bearing was fitted on (ties go to the lower code). Its **purity** is that
majority's share of the length. A link below `lithology_min_majority_fraction`
(default 0.8) goes into a group named `mixed (<80% one unit)`, with the percentage
following the parameter. Every record keeps its majority unit and its purity in
the record columns, so a record can be re-grouped later.

For junctions there are two selections:

- `lithology` (the group key): every junction, by the unit of its own pixel.
- `lithology_interior`: the unit only when the three sampled arms (donor 1, donor 2,
  receiver) all have that unit with purity at or above the threshold; otherwise the
  group is `arms disagree`.

**Descriptive only.** The `mixed` and `arms disagree` groups have summary and rose
rows (with the note `descriptive only: not compared`). They are never compared, they
are never part of `rest` in one-vs-rest mode, and they are not compared with a
reference. The reason is that they have no single unit: pooling them into `rest`
would put a contact effect into a lithology contrast. The run prints how many
records this leaves out of the comparisons.

Unmapped records form a group like any other unit.

### Comparison modes and the Holm adjustment

`channel_bearing_comparison_mode` (also used for junction angles) is one of:

- `pairwise`: every pair of groups in a stratum. Default.
- `one_vs_rest`: each group against the other (eligible) records of its stratum.
- `none`: summaries and roses only.
- `reference` comparisons are added whenever `channel_bearing_reference_fname` is
  given: each group is also compared with the bearings in that file (a
  `_channel_link_bearings.csv` from another run, read from the column
  `channel_bearing_reference_column`, default `bearing_deg`). With
  `channel_bearing_reference_order_filter` > 0 only reference rows of that Strahler
  order are used. Reference comparisons form a family of their own.

Groups smaller than `channel_bearing_min_group_size` (default 30) are summarised
with `meets_min_group_size` = 0 and take no part in any comparison. The Holm
adjustment is applied per family (per stratum and angular type). The comparisons
for bearings are circular Wasserstein distance (W1, in degrees) and Kuiper's V;
for junction angles they are the KS statistic D and linear W1. With
`channel_bearing_rotation_invariant`, the bearing W1 also allows one global
rotation; comparisons with more than 250000 pairs of records are then skipped
(nan, with a note).

### Count and length weighting

Each summary is given twice. The **count** version gives every link weight 1. The
**length** version weights each link by its `total_length_metres` (the whole link,
not the sampled part). Both are written in the group summary.

### Block order, block permutation and the p-values

Links in one valley are not independent, so the p-values that assume independence
are optimistic. Two things are reported beside them.

**Blocks.** `channel_bearing_block_order` = 0 (default) chooses the block order
automatically: the largest Strahler order k whose order-k basins, holding at least
one analysed record, number at least `channel_bearing_min_blocks` (default 30). A
block is one order-k basin, and `n_blocks` in every summary and comparison row is
the number of distinct blocks that contribute records to the group. `n` is the count
of records, `n_blocks` the count of quasi-independent units.

The reason for the default of 30 is a resolution argument. The smallest attainable
permutation p is
1 / (number of distinct label assignments + 1). With 30 blocks even a 5-vs-25
split has C(30, 5) = 142506 assignments, so the number of permutations, not the
block count, sets the resolution. With about 10 blocks an unbalanced split can have
fewer than 100 assignments and the p-value is coarse.

- A positive `channel_bearing_block_order` is used as given. If it has fewer blocks
  than `channel_bearing_min_blocks`, the run notes it.
- **Too small.** If no order reaches the minimum, no blocks are defined. Every
  permutation p is nan, with the note `network too small for block permutation`,
  the descriptive statistics are still written, and a WARNING is printed.
- **Records above the block order** (links of order above k, for the chosen k)
  have no block. Their summary rows say so in `note`, and `n_unblocked` counts them;
  a comparison row that leaves them out of its permutation says so in its `note`.
  They are summarised and enter the iid comparisons, but take no part in the
  block permutation.
- The counts per order are in `<write_fname>_channel_network_block_order_counts.csv`
  (`strahler_order, n_basins_with_records, chosen, note`), so the automatic choice
  can be checked; `note` marks the block order and the orders above it.

**Permutation.** `channel_bearing_permutations` (default 999; 0 = off) and
`channel_bearing_permutation_seed` (default 1) give each comparison a block
permutation p-value:

    permutation_p = (1 + #{T* >= T}) / (1 + permutations)

T is Kuiper's V for bearings and the KS D for junction angles, computed over the
records that have a block. What is permuted depends on the key:

- **Spatial keys** (`basin`, `region`, `lithology`): each block takes the majority
  group of its records in the comparison (ties go to group a), and the group labels
  are permuted among the blocks. `permutation_statistic` is T computed on these
  block-majority samples, so it is not the same number as the record-level T in
  `kuiper_V` or `ks_D`.
- **Order keys** (`strahler_order`, `donor_order_pair`): the group labels are
  permuted among the records within each block, so block membership is kept.

The p is nan, with a note, when a group has no blocked record, when every block
has the same majority group, or, for order keys, when no block holds both groups.
Each comparison row has its own random stream (`std::mt19937_64`, seeded from the
seed and the row's labels), so a row's p does not depend on which other rows were
run.

**Labels.** The iid p-values are named by the statistic: `kuiper_p_iid_approx` (bearing
comparisons), `ks_p_iid_approx` (junction comparisons) and `rayleigh_p_iid_approx`
(the per-group Rayleigh test). Their `_iid_approx` suffix means they assume
independent records and use an approximation to the null distribution; they are
given beside `n_blocks` so the two can be read together. The `kuiper_p_holm` and
`ks_p_holm` columns are the Holm-adjusted iid p-values. `permutation_p` is the
block permutation p-value. These are reported values, not verdicts: read them with
the block counts.

### Outputs

All files are written under `write_fname`. Bearing statistics write four files and
junction statistics two; the block-order counts file is written with either.

| File | Columns |
|---|---|
| `_channel_bearing_group_summary.csv` | `measurand, group_key, group_value, stratum_key, stratum_value, angular_type, period_deg, n, n_blocks, total_length_m`, the count-weighted `circular_mean_deg_count, mean_resultant_length_count, circular_variance_count, circular_sd_deg_count`, the length-weighted `circular_mean_deg_length, mean_resultant_length_length, n_eff_length`, the Rayleigh `rayleigh_z, rayleigh_p_iid_approx`, `fraction_on_d8_multiples, meets_min_group_size, block_order, block_order_mode, min_blocks, n_unblocked, note`. |
| `_channel_bearing_comparisons.csv` | `measurand, stratum_key, stratum_value, angular_type, comparison_mode, group_a, group_b, n_a, n_b, n_blocks_a, n_blocks_b, circular_wasserstein_deg, rotation_invariant, optimal_rotation_deg, kuiper_V, kuiper_p_iid_approx, kuiper_p_holm, n_comparisons_in_family, permutation_p, note, permutation_statistic`. |
| `_channel_bearing_rose.csv` | `measurand, group_key, group_value, stratum_key, stratum_value, angular_type, bin_centre_deg, bin_lower_deg, bin_upper_deg, count, fraction_count, length_m, fraction_length`. Bins are centred on 0, so D8 directions sit on bin centres (`channel_bearing_rose_bin_width_degrees`, default 10). |
| `_channel_bearing_group_membership.csv` | One row per record: its upstream and downstream junction, `stream_order, total_length_m, bearing_deg, axial_bearing_deg`, its `group_value` and `stratum_value`, `basin_id, block_id, region_value, region_purity`. |
| `_junction_angle_group_summary.csv` | `measurand, group_key, group_value, stratum_key, stratum_value, n, n_blocks, mean_deg, median_deg, sd_deg, q25_deg, q75_deg, iqr_deg, meets_min_group_size, median_ci_low_deg, median_ci_high_deg, n_bootstrap, block_order, block_order_mode, min_blocks, n_unblocked, note`. The median interval is given only for lithology-keyed rows (one bootstrap stream per run). |
| `_junction_angle_comparisons.csv` | `measurand, stratum_key, stratum_value, comparison_mode, group_a, group_b, n_a, n_b, n_blocks_a, n_blocks_b, ks_D, ks_p_iid_approx, ks_p_holm, n_comparisons_in_family, linear_wasserstein_deg, permutation_p, note, permutation_statistic`. For lithology rows the `note` carries `group_key=lithology` or `group_key=lithology_interior`. |
| `_channel_network_block_order_counts.csv` | `strahler_order, n_basins_with_records, chosen, note`, one row per order. `note` marks the block order and the orders above it (no block: not in permutations). |

With `lithology_fname` the record files gain columns at the end:

- `_channel_link_bearings.csv`: `lith_code, lith_name, lith_majority_fraction`.
- `_junction_angles.csv`: `junction_lith_code, junction_lith_name, donor_1_lith_code,
  donor_1_lith_fraction, donor_2_lith_code, donor_2_lith_fraction, receiver_lith_code,
  receiver_lith_fraction, lith_all_arms_agree`.

Without `lithology_fname` every record output is unchanged.

Two more outputs come from the lithology map: `<write_fname>_lithology_units.csv`
(`code, name, map_values, n_pixels, area_km2, n_pixels_in_analysed_basins`) and,
with `write_lithology_code_raster`, `<write_fname>_lithology_code.tif`.

### Worked example: West Virginia, the Allegheny Front

This run uses the 2185 x 1148 cell COP30 box of the Allegheny Front at 30 m (EPSG:32617),
with Macrostrat SGMC `lith_primary` grouped into five units (codes 10 sandstone,
20 siltstone, 30 shale, 40 limestone, 90 alluvium) by the crosswalk used in the
chi example. The parameter file is the key lines of the lithology run:

```
read_fname: wvbig_DEM.tif
threshold_contributing_pixels: 300
analyse_all_junctions: true
only_complete_basins: true
link_sampling_mode: full_link
link_sampling_min_nodes: 5
write_channel_link_bearings: true
write_channel_link_bearing_vectors: true
write_junction_angles: true
write_strahler_channels: true
write_channel_bearing_statistics: true
channel_bearing_angular_type: both
channel_bearing_max_order_group: 4
channel_bearing_rose_bin_width_degrees: 10
write_junction_angle_statistics: true
write_fname: wvbig_lith_by_unit
lithology_fname: wvbig_macrostrat_units.fgb
lithology_field_name: lith_primary
lithology_crosswalk_fname: wvks_crosswalk_lith_primary.csv
lithology_min_majority_fraction: 0.8
lithology_bootstrap_replicates: 1000
lithology_bootstrap_random_seed: 5489
write_lithology_code_raster: true
channel_bearing_comparison_mode: one_vs_rest
channel_bearing_group_by: lithology
channel_bearing_stratify_by: strahler_order
junction_angle_group_by: lithology
junction_angle_stratify_by: strahler_order
```

The block order is automatic: order 3 (108 basins with records); order 4 has 23,
below the minimum of 30. Of 2974 links, 732 are `mixed (<80% one unit)` (majority sandstone 266,
siltstone 40, shale 250, limestone 159, alluvium 17); of 953 junctions, 489 are
`arms disagree` in `lithology_interior`. No channel node is unmapped.

Axial bearings, all orders pooled, one-vs-rest against the other comparable records
(counts, n_blocks as in the file, 999 permutations, seed 1):

| Unit | n (rest) | n_blocks (rest) | Axial mean / R (count) | W1 (deg) | Kuiper V | kuiper_p_iid_approx | kuiper_p_holm | permutation_p |
|---|---|---|---|---|---|---|---|---|
| sandstone | 585 (1657) | 66 (102) | 110.6 / 0.047 | 2.68 | 0.069 | 0.215 | 0.464 | 0.562 |
| siltstone | 318 (1924) | 27 (101) | 114.8 / 0.138 | 2.97 | 0.083 | 0.291 | 0.464 | 0.959 |
| shale | 945 (1297) | 82 (91) | 113.4 / 0.124 | 2.47 | 0.065 | 0.154 | 0.464 | 0.678 |
| limestone | 343 (1899) | 46 (106) | 63.7 / 0.023 | 4.50 | 0.118 | 0.008 | 0.041 | 0.863 |
| alluvium | 51 (2191) | 5 (108) | 110.9 / 0.322 | 9.66 | 0.236 | 0.062 | 0.249 | 0.200 |

The `rest` count excludes the mixed records: 2974 - 732 mixed - 585 sandstone =
1657. Mixed has its own row (585, not compared): 732 links, 85 blocks, axial mean
110.7, R 0.262, and no comparison row.

Junction angles, by lithology (median alpha, with the bootstrap interval of the
median, n junctions in brackets):

| Unit | lithology | lithology_interior |
|---|---|---|
| sandstone | 62.3 [60.3, 68.2] (284) | 64.8 [60.1, 72.5] (126) |
| siltstone | 66.2 [58.9, 72.4] (113) | 70.7 [60.3, 76.9] (76) |
| shale | 60.3 [57.0, 65.5] (385) | 58.4 [54.5, 60.9] (206) |
| limestone | 63.9 [59.3, 69.6] (155) | 66.9 [61.2, 80.0] (53) |
| alluvium | 43.6 [28.0, 83.1] (16) | 19.8, no interval (3) |
| arms disagree | - | 62.0 [60.0, 66.2] (489) |

The order run (`wvbig_order`, the same DEM with no lithology map) gives the
Strahler comparisons for the same box. Axial, pairwise: order 1 against order 2 has
W1 1.9 degrees, Kuiper V 0.056, `kuiper_p_iid_approx` 0.38 and `permutation_p` 0.33;
order 1 against order 3 has W1 4.7 degrees, V 0.107, iid p 0.046 and permutation p
0.012. Order 4+ has no order-3 block (n_blocks 0), so its `permutation_p` is nan. The
isotropic-null run (`wvnull_order`, a diamond-square surface matched to the DEM)
chooses block order 4 and gives axial R of 0.063, 0.056, 0.074 and 0.137 for orders
1, 2, 3 and 4+.

The WV numbers are from a live-data run and are not a CTest. The test suite's
own numbers are in the developer testing notes.

### Removed command and reference comparisons

The separate `lsdtt-bearing-comparison` CLI has been removed; this step replaces
it. To compare the bearings of this run with another site, set
`channel_bearing_reference_fname` to the other run's `_channel_link_bearings.csv`.
For arbitrary CSV files, the Python function `lsdtt3.compare_bearings` takes any
bearing CSV.

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
| `basin_outlet_snapping_method` | `string` | `downstream` | `param` | - | Optional snapping of supplied outlet CSV points to channel pixels. | Allowed values: none, expanding_square, downstream, largest_area. largest_area takes the pixel with the most contributing pixels in the square of half-width basin_outlet_snap_max_search_radius_pixels, if it reaches basin_outlet_contributing_pixels_threshold. |
| `basin_outlet_contributing_pixels_threshold` | `int` | `1000` | `param` | - | Contributing-pixel threshold used for optional outlet snapping. | This is total upslope contributing pixels from the flow-routing object. |
| `basin_outlet_snap_min_search_radius_pixels` | `int` | `0` | `param` | - | Minimum expanding-square snapping radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square. |
| `basin_outlet_snap_max_search_radius_pixels` | `int` | `25` | `param` | - | Maximum expanding-square snapping radius in pixels. | Used only when basin_outlet_snapping_method is expanding_square or largest_area (the half-width of the square searched). |

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
| `write_channel_link_bearings` | `bool` | `false` | `switch` | - | Write channel-link bearing output. | Writes one record per analysed compact downstream link. bearing_rad and bearing_deg are clockwise from grid north, or from true north (at the link's representative point) with bearings_relative_to_true_north; direction_x and direction_y stay grid components. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Format for record-oriented network outputs. | Allowed values: csv, gpkg, flatgeobuf or fgb. Line-based QA layers use qa_write_vector_format. |
| `write_lithology_code_raster` | `bool` | `false` | `switch` | - | Write the lithologic unit code of every pixel. | Output: &lt;write_fname&gt;_lithology_code.tif. |

### Junction-angle outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_junction_angles` | `bool` | `false` | `switch` | - | Write junction-angle output. | Writes one record per valid two-donor junction. |
| `write_channel_link_bearing_vectors` | `bool` | `false` | `switch` | - | Write short downstream-pointing visual bearing vectors for analysed links. | This writes one short line per analysed compact downstream link. It is intended as a visual QA layer rather than a full attribute export. |
| `write_junction_angle_bearing_vectors` | `bool` | `false` | `switch` | - | Write visual bearing vectors for junction-angle analysis. | This writes three short centred line features per analysed junction: donor_1, donor_2, and receiver. It is intended as a visual QA layer. |
| `qa_write_vector_format` | `string` | `same_as_write_vector_format` | `param` | - | Output format override for line-based QA layers. | Allowed values: gpkg, flatgeobuf, fgb. CSV is not supported. |

### Network outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_channel_network_lines` | `bool` | `false` | `switch` | - | Write the complete channel network as one line per channel link. | Every link of the network, analysed or not, with junction_upstream, junction_downstream, stream_order, length and drainage_area (at the upstream junction). Join junction_upstream to the up_junction of the link-bearing outputs to tell analysed links from the rest. Written in qa_write_vector_format (no CSV). |
| `write_junction_network_points` | `bool` | `false` | `switch` | - | Write every junction of the channel network as a point. | One point per junction with junction_number, receiver_junction and junction_stream_order. Join junction_number to the junction column of the junction-angle outputs to tell analysed junctions from the rest. Written in write_vector_format (CSV gives latitude/longitude). |

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
| `write_junction_centre_of_mass_bearings` | `bool` | `false` | `switch` | - | Write the basin centre-of-mass bearing for every junction. | For each junction the basin is defined by the pixel just before the next downstream junction (base level for the outlet junction); this writes the bearing (degrees clockwise from grid north, or from true north with bearings_relative_to_true_north) from that outlet to the basin centre of mass, plus the outlet and centre-of-mass coordinates and the basin pixel count. Output: &lt;write_fname&gt;_junction_com_bearings.csv. |

### Complete-basin filtering

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `only_complete_basins` | `bool` | `true` | `param` | - | Restrict junction/bearing/angle analyses to complete basins. | A basin is complete if no pixel in its upslope area is influenced by nodata or the DEM edge. Incomplete (edge-truncated) basins are dropped before computing channel-link bearings, junction angles, and centre-of-mass bearings. Needed for stream-order-based statistics. |
| `min_basin_pixels` | `int` | `0` | `param` | - | Drop junctions whose upslope basin has fewer than this many pixels. | 0 disables the size filter. Applied alongside only_complete_basins. |

### Channel-link bearing statistics

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_channel_bearing_statistics` | `bool` | `false` | `switch` | - | (experimental) Summarise and compare channel-link bearing distributions by group. | Groups the analysed links by channel_bearing_group_by (optionally within each level of channel_bearing_stratify_by) and writes &lt;write_fname&gt;_channel_bearing_group_summary.csv (circular mean, mean resultant length, count- and length-weighted), _channel_bearing_comparisons.csv (circular Wasserstein and Kuiper between groups), _channel_bearing_rose.csv and _channel_bearing_group_membership.csv. Bearings are directional (downstream, period 360) and axial (period 180). p-values assume independent links, which links in one valley are not: read n_blocks. |
| `channel_bearing_group_by` | `string` | `strahler_order` | `param` | `write_channel_bearing_statistics` | (experimental) The key whose groups are compared. | strahler_order (the link's own order), basin (see channel_bearing_basin_order), region (channel_region_fname) or lithology (lithology_fname: the length-weighted majority unit of the nodes the bearing was fitted on; links below lithology_min_majority_fraction form the group 'mixed (&lt;x% one unit)', which is summarised but never compared and never part of 'rest'). |
| `channel_bearing_stratify_by` | `string` | `none` | `param` | `write_channel_bearing_statistics` | (experimental) A key whose levels are analysed separately. | none, strahler_order, basin, region or lithology; must differ from channel_bearing_group_by. Groups are compared only within a level; a pooled level 'all' is always written. |
| `channel_bearing_angular_type` | `string` | `both` | `param` | `write_channel_bearing_statistics` | (experimental) directional, axial or both. | Directional uses the downstream bearing (period 360). Axial folds theta and theta + 180 (period 180); a four-fold D8 artefact cancels in the axial mean resultant length, so axial is the more robust to D8 quantisation. |
| `channel_bearing_comparison_mode` | `string` | `pairwise` | `param` | `write_channel_bearing_statistics` | (experimental) pairwise, one_vs_rest or none. | pairwise compares every pair of groups in a stratum; one_vs_rest each group against the other records of its stratum. Also used for junction angles. Holm-adjusted p-values are given per family (stratum, angular type); comparisons against the reference file are a family of their own. |
| `channel_bearing_min_group_size` | `int` | `30` | `param` | `write_channel_bearing_statistics` | (experimental) Groups smaller than this are summarised but not compared. | Flagged meets_min_group_size = 0. Also used for junction angles. |
| `channel_bearing_max_order_group` | `int` | `0` | `param` | `write_channel_bearing_statistics` | (experimental) Pool Strahler orders &gt;= this into one group labelled '&lt;k&gt;+'. | 0 = no pooling. High-order groups are small, so pooling is usually needed. |
| `channel_bearing_basin_order` | `int` | `0` | `param` | `write_channel_bearing_statistics` | (experimental) Which basins the basin key uses. | 0 = the analysis-selection outlets (the outlet file's outlets, or the network outlets), nested outlets resolved so the innermost wins. k &gt; 0 = the Strahler order-k basin holding the record; records of order &gt; k have none and are excluded (counted). |
| `channel_bearing_block_order` | `int` | `0` | `param` | `write_channel_bearing_statistics` | (experimental) Strahler order of the basins counted as independence blocks; 0 = automatic. | Automatic: the largest order k whose order-k basins holding at least one analysed record number at least channel_bearing_min_blocks; if no order reaches it, no blocks are defined (permutation p-values nan, a WARNING). A positive value is used as given. n_blocks is the number of distinct blocks contributing records to a group: a count of quasi-independent units beside n. Records above the block order have no block: they are summarised but take no part in permutations. The counts per order are written to &lt;write_fname&gt;_channel_network_block_order_counts.csv. |
| `channel_bearing_min_blocks` | `int` | `30` | `param` | `write_channel_bearing_statistics` | (experimental) Minimum number of blocks for the automatic block order. | With B blocks the smallest attainable permutation p is 1 / (number of distinct label assignments + 1). With 30 blocks even a 5-vs-25 split has C(30, 5) = 142506 assignments, so the number of permutations, not the block count, sets the resolution; with about 10 blocks an unbalanced split can have fewer than 100 assignments and the p-value is coarse. |
| `channel_bearing_permutations` | `int` | `999` | `param` | `write_channel_bearing_statistics` | (experimental) Block permutations per comparison; 0 = off. | permutation_p = (1 + #{T* &gt;= T}) / (1 + permutations), T = Kuiper V for bearings, KS D for junction angles, over records with a block. Spatial keys (basin, region, lithology): each block takes the majority group of its records in the comparison and block labels are permuted among blocks (permutation_statistic is T on those block-majority samples). Orders (strahler_order, donor_order_pair): groups are permuted among the records of each block. Descriptive-only groups (lithology mixed, arms disagree) take no part. |
| `channel_bearing_permutation_seed` | `int` | `1` | `param` | `write_channel_bearing_statistics` | (experimental) Seed of the block permutations. | Each comparison row draws from its own std::mt19937_64, seeded from this and the row's labels, so a row's p does not depend on the other rows. |
| `channel_bearing_rose_bin_width_degrees` | `float` | `10.0` | `param` | `write_channel_bearing_statistics` | (experimental) Rose-histogram bin width in degrees. | Bins are centred on 0 so D8 directions sit on bin centres. 180 divided by the width must be an integer. |
| `channel_bearing_rotation_invariant` | `bool` | `false` | `param` | `write_channel_bearing_statistics` | (experimental) Also allow one global rotation in the Wasserstein distance. | Cost grows as n_a * n_b; comparisons with n_a * n_b &gt; 250000 are skipped (nan, with a note). |
| `channel_bearing_reference_fname` | `string` | `NULL` | `param` | `write_channel_bearing_statistics` | (experimental) _channel_link_bearings.csv of another channel-network-analysis run. | Each group is also compared against these bearings (comparison_mode 'reference'), e.g. to compare two sites. Relative to read_path. |
| `channel_bearing_reference_column` | `string` | `bearing_deg` | `param` | `channel_bearing_reference_fname` | (experimental) Column of the reference file holding the bearings. | - |
| `channel_bearing_reference_order_filter` | `int` | `-1` | `param` | `channel_bearing_reference_fname` | (experimental) If &gt; 0, use only reference rows of this order. | Read from channel_bearing_reference_order_column. -1 uses all rows. |
| `channel_bearing_reference_order_column` | `string` | `stream_order` | `param` | `channel_bearing_reference_fname` | (experimental) Column used by channel_bearing_reference_order_filter. | - |

### Junction-angle statistics

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_junction_angle_statistics` | `bool` | `false` | `switch` | - | (experimental) Summarise and compare junction angles by group. | Linear statistics (the angles are sector magnitudes, not directions) of alpha (donor_opening_angle_deg) and beta_min / beta_max (the smaller and larger donor-receiver angle, which do not depend on the arbitrary donor labels). Writes &lt;write_fname&gt;_junction_angle_group_summary.csv and _junction_angle_comparisons.csv (Kolmogorov-Smirnov, linear Wasserstein). Uses the channel_bearing_* comparison settings. |
| `junction_angle_group_by` | `string` | `strahler_order` | `param` | `write_junction_angle_statistics` | (experimental) The key whose junction groups are compared. | strahler_order (the order at the confluence pixel), donor_order_pair ('i_j' from the two donor orders), basin, region (the region of the junction pixel) or lithology (lithology_fname). Lithology writes two selections: group_key lithology (every junction, by the unit of its pixel) and lithology_interior (the same unit when the three sampled arms are all in it with purity &gt;= lithology_min_majority_fraction, else the group 'arms disagree', summarised but never compared). Lithology-keyed rows also get a bootstrap interval of the median (the junction is the resampling unit; lithology_bootstrap_replicates, _random_seed, one random stream per run; none below 5 junctions). |
| `junction_angle_stratify_by` | `string` | `none` | `param` | `write_junction_angle_statistics` | (experimental) A key whose levels are analysed separately. | none or one of the junction_angle_group_by keys; must differ from it. |

### Channel region input

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `channel_region_fname` | `string` | - | `param` | - | Region map (in read_path): a polygon vector file or an integer code raster. | User-drawn regions such as physiographic provinces, used to group channel statistics by region. Vectors: any OGR polygon layer, reprojected to the DEM's CRS. Rasters: integer region codes, warped onto the DEM grid with nearest neighbour if not on it already. Empty (default) = no region grouping. |
| `channel_region_source` | `string` | `auto` | `param` | `channel_region_fname` | auto, vector_map or code_raster. | auto decides from the file extension (.tif, .tiff, .asc, .bil, .img, .vrt, .flt, .nc are rasters). |
| `channel_region_layer_name` | `string` | - | `param` | `channel_region_fname` | Vector layer to read. Empty (default) = the first layer. | - |
| `channel_region_field_name` | `string` | - | `param` | `channel_region_fname` | Attribute of the vector layer whose values name the regions. | Required for a vector map; each distinct value is one region. Ignored for a code raster. |
| `channel_region_min_purity` | `float` | `0.0` | `param` | `channel_region_fname` | (experimental) Minimum length fraction of a link in its majority region. | A link takes the length-weighted majority region of the nodes its bearing was fitted on (as for lithology); links below this purity are excluded from region groups (counted). 0 keeps all; the purity is written to the membership CSV either way. |

### Lithology input

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `lithology_fname` | `string` | - | `param` | - | Lithology map (in read_path) for the grouping key lithology: a polygon vector file or a categorical code raster. Empty (default) = off. | With it, _channel_link_bearings gains lith_code, lith_name, lith_majority_fraction and _junction_angles gains the junction and arm units, and &lt;write_fname&gt;_lithology_units.csv is written. Unmapped records form a group named unmapped; they are never dropped. Vectors: any OGR polygon layer (GeoPackage, FlatGeobuf, GeoJSON, Shapefile), e.g. lsdtt-fetch-lithology output; it is reprojected to the DEM's CRS. Rasters: integer lithology codes, warped onto the DEM grid with nearest neighbour if not on it already. |
| `lithology_source` | `string` | `auto` | `param` | - | auto, vector_map or code_raster. | auto decides from the file extension (.tif, .tiff, .asc, .bil, .img, .vrt, .flt, .nc are rasters). |
| `lithology_layer_name` | `string` | - | `param` | `lithology_fname` | Vector layer to read. Empty (default) = the first layer. | - |
| `lithology_field_name` | `string` | `lith_primary` | `param` | `lithology_fname` | Attribute of the vector layer that holds the lithology. | Each distinct value is a lithologic unit unless lithology_crosswalk_fname groups them. lith_primary is the first lithology of each unit in lsdtt-fetch-lithology output; name (the map unit) or lith_class group more finely or more coarsely. Units are numbered within a run; compare runs by unit name. |
| `lithology_crosswalk_fname` | `string` | - | `param` | `lithology_fname` | Optional CSV (in read_path) with columns value,code[,name] mapping map values to units. | When given, every map value present must be in it. Use it to group many map values into a few units; the optional name column names each code (else the values joined with '; '). |

### Lithology summaries

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `lithology_min_majority_fraction` | `float` | `0.8` | `param` | `lithology_fname` | Minimum fraction of a channel link or junction arm in its majority lithologic unit for it to count as that unit. | The fraction is length-weighted over the sampled nodes. No record is dropped: links below it form the group 'mixed (&lt;x% one unit)' and junctions whose arms fall short the group 'arms disagree'; both are summarised but left out of every comparison (never compared, never part of 'rest'). Every record keeps its majority unit and fraction in the per-record columns. |
| `lithology_bootstrap_replicates` | `int` | `1000` | `param` | `lithology_fname` | Number of bootstrap replicates for the per-unit confidence intervals. | 0 = no bootstrap; the interval columns are written as NaN. |
| `lithology_bootstrap_random_seed` | `int` | `5489` | `param` | `lithology_fname` | Random seed for the per-unit bootstrap. | lsdtt-chi-analysis: each m/n gets its own std::mt19937_64 from seed_seq{seed, round(m/n * 10000)}, so an interval does not depend on the other m/n values in the list. lsdtt-channel-network-analysis: one std::mt19937_64 per run. A fixed seed reproduces the intervals. Same default as bootstrap_random_seed. |

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
