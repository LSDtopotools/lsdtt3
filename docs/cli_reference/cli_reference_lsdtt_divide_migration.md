# CLI reference: lsdtt-divide-migration

Quantifies drainage‑divide migration proxies from user‑supplied outlet points. The tool builds flow routing and a junction/channel network, delineates a basin from each outlet, groups those basins into comparison **units**, extracts the divides separating each pair of adjacent units, then measures every channel head against the nearest divide: elevation, chi, relief, and gradient.

The workflow follows [DivideTools](https://github.com/amforte/DivideTools) (Forte & Whipple), reimplemented natively in C++. Interactive selection (`ginput`) is deliberately replaced by outlet files and attribute fields, so a run is fully reproducible and testable.

The interpretation is the standard one: across a divide, the side whose channel heads sit at **lower** elevation, with **higher** gradient and relief and **lower** chi, is the aggressor — the divide is expected to migrate away from it.

## Usage

```
lsdtt-divide-migration <path> <parameter_file>
lsdtt-divide-migration <parameter_file>     # path defaults to ./
lsdtt-divide-migration -h                   # print help and exit
lsdtt-divide-migration -h <parameter>       # full guidance on one parameter
lsdtt-divide-migration -v                   # print version and exit
lsdtt-divide-migration -c                   # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Status of this tool

`lsdtt-divide-migration` is experimental. **Every parameter it accepts is live** — there are no declared‑but‑inert parameters. The four pair‑summary controls (`write_pair_summary_csv`, `pair_uncertainty_method`, `bootstrap_replicates`, `bootstrap_seed`) became active on 2026‑07‑31 and are documented under *Pair‑summary statistics* below.

Divide fragments **are** merged into continuous polylines (see *Outputs* below). Two units that meet along two genuinely separate stretches of divide still produce two features sharing one `pair_id` — that is correct, not a failure to merge.

Unit pairs are found from **basin adjacency** (2026‑07‑31): a single pass over the basin cells records which basins touch on the grid, and only those pairs are handed to the divide extraction, so the cost no longer grows with the square of the outlet count. Contact is tested eight‑connected, which is how the divide tracing itself sees the grid — basins that meet only at a cell corner can still share a divide — and it is a superset of the pairs that can yield a divide, so the results are exactly what the previous all‑vs‑all loop produced.

One design limit is worth knowing: `nearest_divide_method` has only one implementation (`euclidean`).

---

## Pipeline overview

| Phase | Description |
|---|---|
| 1 | Open DEM; pit‑fill unless `raster_is_filled` |
| 2 | Build flow routing and the junction/channel network |
| 3 | Load outlet points, reproject into the raster CRS, snap them to channels |
| 4 | Resolve each outlet to a **unit** (cluster label, or one unit per outlet) |
| 5 | Delineate a basin per outlet; apply `outlet_control_mode` |
| 6 | Extract the interfluvial divide between every adjacent pair of units |
| 7 | Extract channel heads and integrate chi with chi = 0 at each outlet |
| 8 | Match each head to the nearest point on a divide bounding **its own unit**; measure relief and gradient |
| 9 | Optionally summarise each pair side against side: means, A − B deltas, uncertainty |
| 10 | Write divide vectors, channel‑head vectors, the metrics CSV and the pair summary |

---

## Units, pairs, and what actually gets compared

A **unit** is the thing a divide is measured *between*:

- with `outlet_cluster_field` set, all outlets sharing a label form one unit;
- without it, every outlet is its own unit.

Pairs are formed only where two units genuinely share a divide in the DEM; units that do not touch produce no record. Each divide feature carries `pair_id` as `"<lower_unit>|<higher_unit>"`, e.g. `1|2`.

A channel head is only ever matched to a divide that bounds **its own unit**. A head in unit 3 will never be assigned to the `1|2` divide, even if that divide is far closer in space.

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | `./` | Directory containing the input DEM. |
| `read_fname` | string | _(required)_ | Input raster filename, with extension (e.g. `my_dem.tif`). |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | input stem | Output filename prefix (no path, no extension). |

---

## Threading and flow routing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` uses all available hardware threads. |
| `boundary_conditions` | string | `oooo` | Four‑character N, E, S, W string. `o` = open, `n` = no‑flux, `p` = periodic, `b` = forced outlet. |
| `raster_is_filled` | bool | `false` | If false, `fill_pits()` runs before flow routing. |
| `fill_min_slope` | float | `0.0001` | Minimum slope imposed by pit filling. Used only when `raster_is_filled` is false. |
| `threshold_contributing_pixels` | int | `1000` | Contributing‑pixel threshold defining channel heads. Must be > 0. Lower it for a denser network and more heads — this is the main control on how many rows the metrics CSV has. |

---

## Outlets and clusters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `basin_outlet_fname` | string | `NULL` | **Required.** Point dataset of outlet locations: CSV with latitude/longitude columns, or a point vector (FGB, GPKG, shapefile). Points are reprojected into the raster CRS on load. At least two must land inside the DEM. |
| `outlet_cluster_field` | string | `NULL` | Attribute field holding cluster labels. `NULL` or empty means one unit per outlet. |
| `outlet_cluster_mode` | string | `auto` | How labels are interpreted: `auto` groups numerically when *every* label parses as an integer and as text otherwise; `integer` requires integer labels and rejects the run if any is not; `string` always groups on the raw text. `auto` and `string` differ only when labels look numeric — e.g. `01` and `1` are one group under `auto`, two under `string`. |
| `basin_outlet_snapping_method` | string | `downstream` | `none` (use the point as given), `downstream` (walk downstream to the first cell above the threshold), `expanding_square` (search outward for the nearest such cell), or `largest_area` (the cell with the most contributing pixels in the square of half-width `basin_outlet_snap_max_search_radius_pixels`, if it is above the threshold). |
| `basin_outlet_contributing_pixels_threshold` | int | `1000` | Contributing‑pixel threshold used when snapping. Higher values snap outlets onto larger channels. Must be > 0. |
| `basin_outlet_snap_min_search_radius_pixels` | int | `0` | Minimum expanding‑square search radius. `expanding_square` only. |
| `basin_outlet_snap_max_search_radius_pixels` | int | `25` | Maximum expanding‑square search radius, and the square's half-width for `largest_area`. Must be ≥ the minimum. `expanding_square` and `largest_area` only. |

Outlets that fail to snap are dropped with a message naming the point and the reason; the run continues on the remainder.

---

## Chi and base‑level control

| Parameter | Type | Default | Description |
|---|---|---|---|
| `chi_A0` | float | `1.0` | Reference drainage area A₀ (m²) in `chi = ∫(A₀/A)^(m/n) ds`. **Leave at 1.0** unless you know why: the gradient of a chi-elevation profile equals `k_sn` only when `A0 = 1 m²`; any other value silently rescales it by `A0^(m/n)` (≈631× at `A0 = 1e6`, `m/n = 0.45`). |
| `m_over_n` | float | `0.45` | Reference concavity θ = m/n. Must be > 0. Constrain it with `lsdtt-chi-analysis` and supply the result here. |
| `outlet_control_mode` | string | `none` | Base‑level treatment — see below. |
| `outlet_control_min_elevation` | float | `-9999` | Elevation floor (m). Used only when `outlet_control_mode: elevation`. |

Chi is always integrated with **chi = 0 at each outlet**. Chi values are only comparable across a divide when the basins share a base level, which is what the control modes are for.

| Mode | Behaviour |
|---|---|
| `none` *(default)* | Every delineated basin is kept, including basins whose upslope area touches the DEM edge. |
| `elevation` | Outlets below `outlet_control_min_elevation` are dropped before basin selection. |
| `max_out_elevation` | The highest selected outlet elevation becomes a common datum: each network is truncated below it and chi is re‑zeroed at each retained fragment. A basin can legitimately break into several fragments. |
| `complete_only` | Basins whose upslope area touches nodata or the DEM edge are dropped. |

**Prefer `complete_only` or `max_out_elevation` for any cross‑divide chi comparison on a real DEM.** Under `none`, a basin clipped by the DEM edge has a truncated network and its chi is not comparable with a complete neighbour's.

---

## Divide metrics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `nearest_divide_method` | string | `euclidean` | How heads are linked to divide geometry. Only `euclidean` is implemented: straight‑line distance in map units to the nearest point on a divide segment (projected onto the segment, not snapped to a vertex). |

`divide_elevation_m` is sampled from the routed (pit‑filled) DEM, as is `head_elevation_m`, so `relief_m` never mixes surfaces. Both are in the raster's vertical units.

Edge cases are reported rather than hidden:

- a head lying **on** its divide gives `distance_to_divide_m` ≈ 0 and `gradient_to_divide` = NaN, never an infinity;
- **negative relief** (a head above the nearest divide point) is passed through unchanged — it is a real signal, not an error;
- a nearest point falling on **nodata** still reports `distance_to_divide_m` and `pair_id_nearest`, with the elevation‑derived fields NaN;
- a head whose unit borders no divide keeps NaN divide fields and an empty `pair_id_nearest`.

---

## Outputs

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_divide_vectors` | bool | `true` | `<write_fname>_divides.fgb` \| `.gpkg` — divide segments as polylines. |
| `write_channel_heads` | bool | `true` | `<write_fname>_channel_heads.fgb` \| `.gpkg` — channel heads as points, with the full metric set. |
| `write_metrics_csv` | bool | `true` | `<write_fname>_divide_migration_metrics.csv` — the same metrics as a flat table. |
| `write_pair_summary_csv` | bool | `false` | `<write_fname>_divide_pair_summary.csv` — one side‑vs‑side row per divide pair. See *Pair‑summary statistics* below. |
| `write_vector_format` | string | `flatgeobuf` | `flatgeobuf` (or `fgb`) or `gpkg`. Applies to both vector outputs. |

Three provenance files are always written regardless of the switches above: `<write_fname>_README.html` and `<write_fname>_README.csv` (the parameter reference for this run) and `<write_fname>_ingested_params.param` (every parameter as actually parsed, including defaults you did not set — useful for reproducing a run exactly).

### Divide vector attributes

| Field | Type | Description |
|---|---|---|
| `pair_id` | string | `"<unit_a>\|<unit_b>"`, units in ascending order. |
| `unit_a`, `unit_b` | int | The two units this divide separates. |
| `junc_a`, `junc_b` | int | Junction numbers of the two source basins, or `-1` when the feature was merged from more than one fragment (see below). |
| `segment_id` | int | Unique per segment within the run. |
| `n_fragments` | int | How many extracted fragments were merged into this feature. `1` for a divide between a single pair of basins. |
| `n_vertices` | int | Vertex count of the polyline. |
| `length_m` | double | Segment length in map units, measured on the polyline as written. |

A unit made of several basins borders a neighbouring unit through one fragment
per adjacent basin pair. Fragments that are contiguous in space — endpoints
within one cell diagonal of each other — are merged into a single polyline, so
one feature is one divide and `length_m` is the whole along‑divide length. Two
units that meet along two genuinely separate stretches of divide (wrapping
around a third unit, say) still get one feature per stretch: sharing a `pair_id`
does not merge anything on its own.

Because `junc_a` / `junc_b` describe one basin pair, they are only meaningful on
a feature with `n_fragments = 1`; a merged feature reports `-1` for both rather
than one contributing basin's junction as though it described the whole line.

### Channel‑head attributes (vector and CSV)

| Field | Description |
|---|---|
| `head_id` | Unique per head within the run. |
| `outlet_id` | `point_id` of the outlet whose basin contains this head. |
| `unit_id` | Unit that basin belongs to. |
| `pair_id_nearest` | `pair_id` of the nearest divide; empty if unmatched. |
| `head_x`, `head_y` | Head position in the raster CRS. |
| `head_elevation_m` | Head elevation, from the routed DEM. |
| `chi` | Chi at the head, zeroed at its outlet. |
| `nearest_divide_x`, `nearest_divide_y` | Nearest point on the divide. |
| `divide_elevation_m` | Elevation there, from the routed DEM. |
| `distance_to_divide_m` | Straight‑line distance, map units. |
| `relief_m` | `divide_elevation_m − head_elevation_m`. May be negative. |
| `gradient_to_divide` | `relief_m / distance_to_divide_m`. NaN at zero distance. |

---

## Pair‑summary statistics

Set `write_pair_summary_csv: true` and the run writes `<write_fname>_divide_pair_summary.csv`: **one row per divide pair**, comparing the channel heads on one side against those on the other. This is the table you interpret; the per‑head CSV is the evidence behind it.

The row set comes from the **divides**, not from the heads: every pair present in the `_divides` layer gets exactly one row, so the two outputs always list the same pairs. A pair whose divide was extracted but which attracted no channel heads on either side still gets a row, with `n_heads_a` and `n_heads_b` at 0 and every mean, delta and uncertainty `nan` — the same convention a single empty side follows. Two separate stretches of divide between the same unit pair share one `pair_id` and are summarised together, in one row.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_pair_summary_csv` | bool | `false` | Write the pair‑summary CSV. |
| `pair_uncertainty_method` | string | `std_err` | `std_dev`, `std_err`, `bootstrap` or `ttest` — see below. |
| `bootstrap_replicates` | int | `1000` | Resamples per side. Read only when the method is `bootstrap`. Must be > 0. |
| `bootstrap_seed` | int | `5489` | Seed for the bootstrap. Read only when the method is `bootstrap`. |

### The sign convention: side A is the lower‑numbered unit

**Side A is always the lower‑numbered unit of the pair** — the same ordering the `pair_id` itself uses, so pair `3|7` has `unit_a = 3` and `unit_b = 7`. Every delta is

```
delta_<metric> = mean_<metric>_a − mean_<metric>_b
```

This is fixed by the unit ids, never by the order heads happen to arrive in, so a column cannot silently change sign between runs. The run prints a reminder on stdout.

### Rows, grouping and sample sizes

Heads are grouped by `pair_id_nearest`. A head that matched no divide (empty `pair_id_nearest`) takes no part in any row, and a head is only ever counted on the side its own `unit_id` names.

`n_heads_a` / `n_heads_b` count the heads assigned to each side. They are **not** the sample size of every mean: `relief_m` and `gradient_to_divide` are legitimately NaN for some heads (the nearest divide point fell on nodata, or the head sits on its divide), and those heads are **skipped** rather than allowed to turn a side mean into NaN. Each metric therefore carries its own valid‑value count, `n_<metric>_a` / `n_<metric>_b`, which is the n its mean *and* its uncertainty were computed from. If `n_relief_a` is well below `n_heads_a`, treat that side's relief with caution.

A pair in which one side has **no** heads still gets a row, with that side's counts at 0 and its mean, the delta and the uncertainty NaN. The divide exists in the vector output, so a row saying one side is empty is more useful than a divide that vanishes from the table.

### What each uncertainty method reports

The column is `uncertainty_<metric>` — one number per metric — and the row's `uncertainty_method` field always records which of the four produced it.

DivideTools' `AcrossDivide.m` reports a dispersion for *each* side and then asks whether the two intervals overlap. This tool has one column, so each option is the equivalent statistic **of the delta**. The per‑side ingredients are DivideTools' (the n−1 sample standard deviation, `s/√n`, a bootstrap of the mean, and Welch's unequal‑variance t‑test); combining them into a single number per metric is this tool's choice, not something inherited from DivideTools.

| `pair_uncertainty_method` | The `uncertainty_*` columns hold | Compare it with |
|---|---|---|
| `std_dev` | Pooled standard deviation of the two sides, `sqrt(((na−1)sa² + (nb−1)sb²)/(na+nb−2))`. A spread of the **head populations**, not of the delta — it does not shrink as you add heads. | Nothing directly; use it to see how variable the two sides are. |
| `std_err` *(default)* | Standard error of the **difference of means**, `sqrt(sa²/na + sb²/nb)`. | The delta. `\|delta\| > 2 ×` this is roughly a 95 % separation of the two side means. |
| `bootstrap` | Half‑width of the 95 % confidence interval of the **bootstrap distribution of the delta**. Both sides are resampled with replacement `bootstrap_replicates` times and the delta recorded each time; the reported value is the mean of the distances from the bootstrap mean to the 2.5th and 97.5th percentiles — the statistic `AcrossDivide.m`'s `bootCI` helper returns. | The delta directly: `\|delta\| >` this is a ~95 % separation. Roughly 1.96 × `std_err` for a large, well‑behaved sample, so **never** compare a bootstrap column against a std_err column from another run. |
| `ttest` | The **two‑sided p‑value** of Welch's unequal‑variance t‑test on the difference of means, matching `ttest2(...,'Vartype','unequal')`. This column is a probability, not a dispersion. | 0.05 (or your own threshold). Small means the two sides genuinely differ. |

**Small samples.** Every method needs a within‑side spread, so all four require at least **two valid values on both sides** of the metric concerned and return **NaN** otherwise. A side with one head, or none, gives a NaN uncertainty rather than a meaningless zero.

**Reproducibility.** The bootstrap uses an explicitly seeded `std::mt19937`; the same seed and the same input always give bit‑identical output. The seed is combined with the pair id and the metric rather than with a loop counter, so adding an unrelated pair to a run cannot change an existing pair's numbers, and the result does not depend on thread count.

### Pair‑summary CSV columns

| Field | Description |
|---|---|
| `pair_id` | `"<unit_a>\|<unit_b>"`, matching the divide vector output. |
| `unit_a`, `unit_b` | The two units. `unit_a < unit_b` always. |
| `n_heads_a`, `n_heads_b` | Channel heads assigned to each side. |
| `mean_head_elevation_a`, `mean_head_elevation_b`, `delta_head_elevation` | Head elevation (m), and A − B. |
| `n_head_elevation_a`, `n_head_elevation_b` | Valid values behind those means. |
| `mean_gradient_a`, `mean_gradient_b`, `delta_gradient` | Gradient to divide, and A − B. |
| `n_gradient_a`, `n_gradient_b` | Valid values behind those means. |
| `mean_relief_a`, `mean_relief_b`, `delta_relief` | Relief to divide (m), and A − B. |
| `n_relief_a`, `n_relief_b` | Valid values behind those means. |
| `mean_chi_a`, `mean_chi_b`, `delta_chi` | Chi, and A − B. |
| `n_chi_a`, `n_chi_b` | Valid values behind those means. |
| `uncertainty_method` | `std_dev`, `std_err`, `bootstrap` or `ttest`. |
| `uncertainty_head_elevation`, `uncertainty_gradient`, `uncertainty_relief`, `uncertainty_chi` | As defined in the table above. |

### Reading a row

Across a divide, the **aggressor** is the side whose channel heads sit at *lower* elevation, with *higher* gradient and relief, and *lower* chi. The divide is expected to migrate *away* from the aggressor, into the victim.

With A − B deltas, side A is the aggressor when:

| Column | Points to A as aggressor when |
|---|---|
| `delta_head_elevation` | **negative** (A's heads are lower) |
| `delta_gradient` | **positive** (A is steeper towards the divide) |
| `delta_relief` | **positive** (more relief between A's heads and the divide) |
| `delta_chi` | **negative** (A's heads have lower chi) |

Reverse every sign and B is the aggressor.

Judge each metric against its uncertainty before believing it. Under `std_err` or `bootstrap`, a metric whose `|delta|` does not clearly exceed its uncertainty column is **not** evidence of anything; under `ttest`, one whose p‑value is above your threshold is not either. A divide with all four metrics small relative to their uncertainties is best read as **stable**, which is a real and common result, not a failed analysis.

The four metrics can and do disagree. That is informative rather than a problem — chi responds to the drainage‑network geometry and is sensitive to base‑level and to `m_over_n`, while gradient and relief respond to the present hillslope form. Cross‑divide chi is only comparable at all when the basins share a base level, so a chi delta from a run with `outlet_control_mode: none` on a DEM whose basins are clipped at the edge should be discounted. Agreement between elevation, gradient and relief, with chi consistent, is the strong case.

---

## Example: minimal run, no clusters

Every outlet is its own unit, so a divide is measured between each adjacent pair of basins.

```
read_path:   /data/santa_cruz/
read_fname:  santa_cruz.tif
write_path:  /data/santa_cruz/divides/
write_fname: sc

# Required: the outlet points
basin_outlet_fname: outlets.csv

# Channel network density: lower for more heads
threshold_contributing_pixels: 500

# Chi settings -- constrain m_over_n with lsdtt-chi-analysis first
chi_A0:   1.0
m_over_n: 0.45

# Only compare basins that are complete within the DEM
outlet_control_mode: complete_only
```

## Example: clustered outlets

Here `outlets.gpkg` carries a `range_side` attribute with values such as `west` and `east`. All outlets sharing a value collapse into one unit, so the run measures the divide **between the two sides of the range** rather than between individual catchments.

```
read_path:   /data/santa_cruz/
read_fname:  santa_cruz.tif
write_path:  /data/santa_cruz/divides/
write_fname: sc_sides

basin_outlet_fname:   outlets.gpkg
outlet_cluster_field: range_side
outlet_cluster_mode:  string

# Snap outlets onto trunk channels
basin_outlet_snapping_method: downstream
basin_outlet_contributing_pixels_threshold: 5000

threshold_contributing_pixels: 500
m_over_n: 0.45

# Impose a common base level so chi is comparable across the divide
outlet_control_mode: max_out_elevation

write_metrics_csv:    true
write_channel_heads:  true
write_vector_format: gpkg

# The side-vs-side table, with a bootstrap uncertainty on each delta
write_pair_summary_csv:  true
pair_uncertainty_method: bootstrap
bootstrap_replicates:    2000
bootstrap_seed:          5489
```

This writes `sc_sides_divide_pair_summary.csv` with one row per pair of range sides. If units are `1` (west) and `2` (east), then side A is west and a **negative** `delta_head_elevation` with a **positive** `delta_gradient` and `delta_relief` and a **negative** `delta_chi` says the west side is the aggressor and the divide should be migrating east — provided each `|delta|` clearly exceeds its `uncertainty_*` column. Because the uncertainty method here is `bootstrap`, those columns are 95 % interval half‑widths, so compare `|delta|` against them directly.

If you would rather do the statistics yourself, the per‑head CSV has everything: group by `unit_id`, restrict to the rows whose `pair_id_nearest` names the divide of interest, and take the means of `head_elevation_m`, `gradient_to_divide`, `relief_m` and `chi`.

---

## Troubleshooting

### "basin_outlet_fname must provide at least two points that map into the DEM extent"

Almost always a **CRS or extent mismatch**. Outlet points are reprojected into the raster CRS on load, so this means the reprojected points fall outside the DEM, not that the file is unreadable. Check:

- the DEM has a valid CRS (`gdalinfo my_dem.tif` — a missing or wrong projection is the usual cause);
- for a CSV, that the columns really are latitude/longitude in decimal degrees and are not swapped (a lat/long transposition puts points in the wrong hemisphere);
- for a vector, that the file's own CRS is declared correctly;
- that the outlets genuinely lie inside the DEM footprint — compare `gdalinfo` corner coordinates against `ogrinfo -al -so outlets.gpkg`.

### "fewer than two valid snapped outlets in two distinct units remain for unit mapping (outlets=N, units=M)"

Divide extraction needs **two units to compare**. The message gives both counts, which tells you which problem you have:

- `units=1` with `outlets` ≥ 2: every outlet resolved to the same cluster. Check `outlet_cluster_field` names the column you meant and that its values actually differ.
- both counts low: outlets were dropped during snapping. Look for the per‑point snap messages above the error and lower `basin_outlet_contributing_pixels_threshold`, or switch `basin_outlet_snapping_method` to `expanding_square` (or `none` if your outlets are already on channel cells).

### "could not load outlet clusters from field '<field>'"

The field name was not found. **Matching differs by file type**, which matters if the name looks right:

- **CSV**: an exact match on the header is tried first, then a case‑insensitive one. So `Cluster_ID` will find `cluster_id`.
- **Vector (FGB/GPKG/shapefile)**: matching is **exact and case‑sensitive**. `Cluster_ID` will *not* find `cluster_id`.

Other causes are leading/trailing whitespace in a CSV header, or pointing at a file that does not carry the attribute at all. For a vector file, `ogrinfo -al -so <file>` lists the fields with their exact spelling; for a CSV, check the header row. Note that shapefiles truncate field names to 10 characters.

### "outlet cluster value '<x>' is not an integer (point_id=N)"

Raised only under `outlet_cluster_mode: integer`. Either clean the labels or use `auto` / `string`, which accept arbitrary text.

### "empty outlet cluster value for point_id=N"

A blank cell in the cluster column. Every outlet must carry a label — there is no implicit "unassigned" group.

### "two outlets from different units snap to the same flat index"

Two outlets in *different* clusters snapped onto the same cell, so the basin cannot be attributed to a unit. Move one outlet, or lower the snapping threshold so they resolve to distinct channel cells.

### "fewer than two unique basins remain after selection"

Filtering removed too much. Under `outlet_control_mode: complete_only` this usually means the basins touch the DEM edge — fetch a larger DEM, or fall back to `none` and accept that chi is not strictly comparable. Under `elevation`, the floor is probably above your outlets.

### The metrics CSV has rows but every divide field is NaN

The heads were not matched to any divide. Either only one unit had basins that survived filtering, or the units are not adjacent in the DEM (units separated by a third basin share no divide). Check `divide segments` in the run summary — if it is `0`, no divides were extracted at all.

### The pair‑summary CSV has a row but the uncertainties are NaN

Every uncertainty method needs at least **two valid values on both sides** of the metric. Check `n_<metric>_a` and `n_<metric>_b` on the row:

- one of them is `0` or `1`: that side has too few heads. Lower `threshold_contributing_pixels` for a denser network, or check that side's basins survived `outlet_control_mode` filtering.
- `n_heads_a` is healthy but `n_relief_a` or `n_gradient_a` is small: most of that side's heads have NaN relief, which means their nearest divide point landed on nodata — see the all‑NaN case above.

A NaN here is deliberate: a standard error from one value would be meaningless, and reporting `0` would be worse than reporting nothing.

### The pair‑summary CSV is missing but the metrics CSV is there

`write_pair_summary_csv` defaults to **false**. Set it to `true`.

### The pair‑summary CSV has a header but no rows

No head was matched to any divide, so there was nothing to summarise. The file is still written so the absence is visible rather than silent. Check `divide segments` and `heads with divide` in the run summary — see *The metrics CSV has rows but every divide field is NaN* above for the usual causes.

### Too few or too many channel heads

`threshold_contributing_pixels` is the control. It is a contributing‑*pixel* count, so its meaning in area depends on cell size: at 30 m, `1000` pixels ≈ 0.9 km². Lower it for a denser network.

<!-- BEGIN GENERATED CLI PARAMETERS -->
## Parameters

This section is generated from the executable's `-h json` registry. Do not edit it manually.

**Status: experimental.** This CLI is not yet fully tested.

### File I/O

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `read_path` | `string` | `./` | `io` | - | Directory from which input files are read. | Defaults to the directory holding the parameter file. |
| `read_fname` | `string` | - | `io` | - | Input filename, with or without extension (e.g. my_dem.tif). | Read from read_path. lsdtt3 uses GDAL, so any GDAL-readable raster format works. If empty, defaults to the parameter filename stem. |
| `write_path` | `string` | - | `io` | - | Directory to which output files are written. | Defaults to read_path when not set. |
| `write_fname` | `string` | - | `io` | - | Output filename prefix (no path, no extension). | All outputs are named &lt;write_fname&gt;_&lt;descriptor&gt;.&lt;ext&gt;. Defaults to the stem of read_fname when not set. |

### Routing / DEM

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four-character string for north, east, south, west. o = open, n = no-flux, p = periodic, b = forced outlet. |
| `fill_min_slope` | `float` | `0.0001` | `param` | `raster_is_filled` | Minimum slope used when pit filling is enabled. | Only used when raster_is_filled is false. |
| `n_threads` | `int` | `0` | `param` | - | Number of OpenMP threads. | 0 = all available hardware threads. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Outlets / clusters

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_outlet_fname` | `string` | `NULL` | `param` | - | Point dataset of outlet locations. | Required. Supports CSV latitude/longitude or point vectors (FGB/GPKG/etc.). |
| `outlet_cluster_field` | `string` | `NULL` | `param` | - | Optional attribute field naming outlet clusters. | If NULL/empty, each outlet is treated as its own unit. |
| `outlet_cluster_mode` | `string` | `auto` | `param` | `outlet_cluster_field` | How to interpret cluster field values. | Options: auto (numeric grouping if every label is an integer, otherwise string grouping), integer (labels must be integers, grouped numerically), string (labels always grouped as raw text, so 01 and 1 are different units). |
| `basin_outlet_snapping_method` | `string` | `downstream` | `param` | - | Outlet snapping method. | Options: none, downstream, expanding_square, largest_area. largest_area takes the pixel with the most contributing pixels within basin_outlet_snap_max_search_radius_pixels. |
| `basin_outlet_contributing_pixels_threshold` | `int` | `1000` | `param` | `basin_outlet_snapping_method` | Contributing-pixel threshold used in outlet snapping. | Higher values snap outlets onto larger channels. |
| `basin_outlet_snap_min_search_radius_pixels` | `int` | `0` | `param` | `basin_outlet_snapping_method` | Minimum expanding-square snap radius (pixels). | Used only for basin_outlet_snapping_method=expanding_square. |
| `basin_outlet_snap_max_search_radius_pixels` | `int` | `25` | `param` | `basin_outlet_snapping_method` | Maximum expanding-square snap radius (pixels). | Used only for basin_outlet_snapping_method=expanding_square or largest_area (the half-width of the square searched). |

### Chi

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `chi_A0` | `float` | `1.0` | `param` | - | Reference area A0 for chi integration (m^2). | Used in chi = integral((A0/A)^(m/n) ds). |
| `m_over_n` | `float` | `0.45` | `param` | - | Reference m/n (theta) for chi integration. | Must be &gt; 0. |
| `outlet_control_mode` | `string` | `none` | `param` | - | Outlet/base-level control mode. | Options: none, elevation, max_out_elevation, complete_only. |
| `outlet_control_min_elevation` | `float` | `-9999` | `param` | `outlet_control_mode` | Minimum base-level elevation for outlet_control_mode=elevation. | Ignored for other outlet_control_mode values. |

### Divide metrics

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `nearest_divide_method` | `string` | `euclidean` | `param` | - | Method to link channel heads to divide geometry. | Current phase supports: euclidean. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_divide_vectors` | `bool` | `true` | `io` | - | Write divide segment vectors. | Output path stem: &lt;write_fname&gt;_divides. |
| `write_channel_heads` | `bool` | `true` | `io` | - | Write channel-head point vectors. | Output path stem: &lt;write_fname&gt;_channel_heads. |
| `write_metrics_csv` | `bool` | `true` | `io` | - | Write per-channel-head metrics CSV. | Output path: &lt;write_fname&gt;_divide_migration_metrics.csv. |
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format. | Options: flatgeobuf/fgb or gpkg. |
| `write_pair_summary_csv` | `bool` | `false` | `io` | - | Write the side-vs-side pair summary CSV. | One row per divide pair with side means, A-minus-B deltas (A is the LOWER-numbered unit) and an uncertainty per metric. Output path: &lt;write_fname&gt;_divide_pair_summary.csv. |
| `pair_uncertainty_method` | `string` | `std_err` | `param` | `write_pair_summary_csv` | Uncertainty method for pair summaries. | Options: std_dev (pooled standard deviation of the two sides), std_err (standard error of the difference of means), bootstrap (95% CI half-width of the bootstrapped delta), ttest (two-sided Welch p-value). Needs two valid values on both sides, else NaN. |
| `bootstrap_replicates` | `int` | `1000` | `param` | `pair_uncertainty_method` | Bootstrap replicate count for pair summary uncertainty. | Used when pair_uncertainty_method is bootstrap. Must be &gt; 0. |
| `bootstrap_seed` | `int` | `5489` | `param` | `pair_uncertainty_method` | Random seed for bootstrap uncertainty. | Used when pair_uncertainty_method is bootstrap. The same seed and the same input always give bit-identical output. |

<!-- END GENERATED CLI PARAMETERS -->
