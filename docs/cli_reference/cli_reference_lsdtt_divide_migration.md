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

`lsdtt-divide-migration` is at version 0.1.0 and some of its parameter surface is declared ahead of its implementation. **Four parameters are accepted and validated but currently have no effect**: `write_pair_summary_csv`, `pair_uncertainty_method`, `bootstrap_replicates`, `bootstrap_seed`. No pair‑summary CSV is written yet; the per‑head metrics CSV carries everything needed to compute side‑by‑side comparisons downstream.

One further limitation to be aware of: divide fragments are **not merged into continuous polylines**. Where a unit consists of several basins, each adjacent basin pair contributes its own segment, so more than one feature may share a `pair_id`.

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
| 9 | Write divide vectors, channel‑head vectors, and the metrics CSV |

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
| `source_pixels_threshold` | int | `1000` | Contributing‑pixel threshold defining channel heads. Must be > 0. Lower it for a denser network and more heads — this is the main control on how many rows the metrics CSV has. |

---

## Outlets and clusters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `basin_outlet_file` | string | `NULL` | **Required.** Point dataset of outlet locations: CSV with latitude/longitude columns, or a point vector (FGB, GPKG, shapefile). Points are reprojected into the raster CRS on load. At least two must land inside the DEM. |
| `outlet_cluster_field` | string | `NULL` | Attribute field holding cluster labels. `NULL` or empty means one unit per outlet. |
| `outlet_cluster_mode` | string | `auto` | How labels are interpreted: `auto` groups numerically when *every* label parses as an integer and as text otherwise; `integer` requires integer labels and rejects the run if any is not; `string` always groups on the raw text. `auto` and `string` differ only when labels look numeric — e.g. `01` and `1` are one group under `auto`, two under `string`. |
| `basin_outlet_snapping_method` | string | `downstream` | `none` (use the point as given), `downstream` (walk downstream to the first cell above the threshold), or `expanding_square` (search outward for the nearest such cell). |
| `basin_outlet_contributing_pixels_threshold` | int | `1000` | Contributing‑pixel threshold used when snapping. Higher values snap outlets onto larger channels. Must be > 0. |
| `basin_outlet_min_search_radius_pixels` | int | `0` | Minimum expanding‑square search radius. `expanding_square` only. |
| `basin_outlet_max_search_radius_pixels` | int | `25` | Maximum expanding‑square search radius. Must be ≥ the minimum. `expanding_square` only. |

Outlets that fail to snap are dropped with a message naming the point and the reason; the run continues on the remainder.

---

## Chi and base‑level control

| Parameter | Type | Default | Description |
|---|---|---|---|
| `chi_A0` | float | `1.0` | Reference drainage area A₀ (m²) in `chi = ∫(A₀/A)^(m/n) ds`. |
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
| `vector_output_format` | string | `flatgeobuf` | `flatgeobuf` (or `fgb`) or `gpkg`. Applies to both vector outputs. |

Three provenance files are always written regardless of the switches above: `<write_fname>_README.html` and `<write_fname>_README.csv` (the parameter reference for this run) and `<write_fname>_ingested_params.param` (every parameter as actually parsed, including defaults you did not set — useful for reproducing a run exactly).

### Divide vector attributes

| Field | Type | Description |
|---|---|---|
| `pair_id` | string | `"<unit_a>\|<unit_b>"`, units in ascending order. |
| `unit_a`, `unit_b` | int | The two units this divide separates. |
| `junc_a`, `junc_b` | int | Junction numbers of the two source basins. |
| `segment_id` | int | Unique per segment within the run. |
| `n_vertices` | int | Vertex count of the polyline. |
| `length_m` | double | Segment length in map units. |

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

## Example: minimal run, no clusters

Every outlet is its own unit, so a divide is measured between each adjacent pair of basins.

```
read_path:   /data/santa_cruz/
read_fname:  santa_cruz.tif
write_path:  /data/santa_cruz/divides/
write_fname: sc

# Required: the outlet points
basin_outlet_file: outlets.csv

# Channel network density: lower for more heads
source_pixels_threshold: 500

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

basin_outlet_file:    outlets.gpkg
outlet_cluster_field: range_side
outlet_cluster_mode:  string

# Snap outlets onto trunk channels
basin_outlet_snapping_method: downstream
basin_outlet_contributing_pixels_threshold: 5000

source_pixels_threshold: 500
m_over_n: 0.45

# Impose a common base level so chi is comparable across the divide
outlet_control_mode: max_out_elevation

write_metrics_csv:    true
write_channel_heads:  true
vector_output_format: gpkg
```

To compare the two sides, group the metrics CSV by `unit_id` and take the mean of `head_elevation_m`, `gradient_to_divide`, `relief_m`, and `chi` for the rows whose `pair_id_nearest` names the divide of interest. (An automated pair‑summary table is planned; see *Status* above.)

---

## Troubleshooting

### "basin_outlet_file must provide at least two points that map into the DEM extent"

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

### Too few or too many channel heads

`source_pixels_threshold` is the control. It is a contributing‑*pixel* count, so its meaning in area depends on cell size: at 30 m, `1000` pixels ≈ 0.9 km². Lower it for a denser network.
