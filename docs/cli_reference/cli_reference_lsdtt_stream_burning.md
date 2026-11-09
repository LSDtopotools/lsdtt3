# CLI reference: lsdtt-stream-burning

Conditions a DEM using a user-supplied vector channel network. The tool builds two parallel network views (vector topology vs D8-routed DEM channels), then applies one of three `stream_burn_mode` workflows:

- `burn_in_place` — directly lower DEM cells under the vector polylines to enforce downhill channel profiles.
- `snap_to_network` — compare vector and DEM junction topology, then optionally repair mismatches/gaps.
- `snap_to_accumulation` — for accurate vectors (for example OS VectorMap Local), trace D8 flow along the vector and, with `run_lcp_burn: true`, least-cost carve only the unresolved connectivity gaps between accepted vector reaches to recover a continuous routed network.

## Usage

```
lsdtt-stream-burning <path> <parameter_file>
lsdtt-stream-burning <parameter_file>         # path defaults to ./
lsdtt-stream-burning -h                       # print help and exit
lsdtt-stream-burning -v                       # print version and exit
lsdtt-stream-burning -c                       # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

The `stream_vector_path` parameter is required.

---

!!! note "Multipart lines"
    `stream_vector_path` accepts multipart geometry (`MultiLineString` and friends). Each feature is split into its parts on read, with the attributes copied onto each, and the split is reported in the log.

    Shapefiles are the commonest source of multipart data, and **versions before August 2026 discarded it silently** — see [Vector data](../what_is_lsdtt3/vector_data.md).

## Pipeline overview

All runs execute the common preparation and network-building stages. The mode-specific processing is selected by `stream_burn_mode`.

| Stage | Description |
|---|---|
| 1 | Open DEM; optional breach and pit-fill |
| 2 | Write filled DEM (`_filled.tif`) if filling was performed |
| 3 | Build D8 flow routing from the filled DEM |
| 4 | Load the vector channel file; build clipped vector junction network |
| 5 | Write `_vector_junctions.fgb` and `_vector_edges.fgb` |
| 6 | Extract source junctions from the vector network (degree-1 non-outlet tips) |
| 7 | Map each vector source to the DEM pixel that contains it |
| 8 | Write `_vector_endpoints` (optional, controlled by `write_vector_sources`) |
| 9 | Build DEM channel network seeded from the mapped source junctions |
| 10 | Write DEM network points / lines / junctions (optional) |
| 11 | **Mode: `burn_in_place`** — lower DEM pixels under vector polylines to enforce monotone profiles |
| 11 | **Mode: `snap_to_network`** — find candidate mismatched junctions and optionally repair them |
| 11 | **Mode: `snap_to_accumulation`** — join adjacent water-line segments that connect end-to-end (simple pass-through joins, not confluences) into longer lines (`_merged_lines.fgb`), then run D8 trace acceptance along those lines (`_accepted_segments.fgb` -> `_accepted_clean.fgb`), and if `run_lcp_burn` is enabled carve only the unresolved obstruction gaps; optionally write the carved DEM (`_lcp_burned.tif`) and the post-burn connected network (`_connected_network_lines.fgb`) |

The `snap_to_accumulation` workflow has two stages, controlled independently:

- **Initial stage** (always runs): breach/fill → flow → clipped+repaired vector network → merged vector lines → source-seeded D8 network (`_dem_network_lines.fgb`). Setting `run_stream_burn: false` stops here, so no DEM conditioning is applied.
- **LCP burn stage** (`run_lcp_burn: true`, requires `segmented_trace_test: true`): bridge only the unresolved gaps between accepted vector reaches with least-cost-path carving, then (if requested) write the conditioned DEM (`_lcp_burned.tif`) and the fixed connected network extracted from it (`_connected_network_lines.fgb`).

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input raster filename. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` uses all available hardware threads. |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `raster_is_filled` | bool | `false` | Set `true` if the DEM is already pit-filled. Skips fill and breach steps. |
| `carve_before_fill` | bool | `false` | Run depression breaching before pit-filling. Useful when road or embankment obstructions block channels. |
| `breach_algorithm_choice` | string | `lindsay2016` | Breaching algorithm: `lindsay2016` or `least_cost`. Used only when `carve_before_fill` is true. |
| `fill_min_slope` | float | `0.0001` | Minimum slope imposed during pit-filling. |

---

## Vector channel input

| Parameter | Type | Default | Description |
|---|---|---|---|
| `stream_vector_path` | string | _(required)_ | Path to the polyline vector channel file. Must be a linestring or multilinestring layer in the same CRS as the DEM. Supported: GeoPackage (`.gpkg`), FlatGeobuf (`.fgb`), Shapefile (`.shp`), and any GDAL/OGR format. |

---

## Routing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `boundary_conditions` | string | `oooo` | D8 flow-routing boundary conditions: four characters for north, east, south, west. `o` = open, `n` = no-flux, `p` = periodic, `b` = forced outlet. |

---

## Burn mode

| Parameter | Type | Default | Description |
|---|---|---|---|
| `stream_burn_mode` | string | `snap_to_network` | How the vector network conditions the DEM. `burn_in_place`: directly lowers DEM pixels under the vector polylines to enforce monotone-decreasing elevation profiles. `snap_to_network`: compares vector junctions against DEM-extracted junctions and produces diagnostic outputs showing where the two networks agree or diverge (best for less accurate layers, e.g. OS Rivers). `snap_to_accumulation`: for accurate vector layers (e.g. OS VectorMap Local) — traces D8 flow along the vector and bridges obstruction gaps with least-cost-path carving (see the LCP-burn section below). |

---

## burn_in_place outputs

Used only when `stream_burn_mode: burn_in_place`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_burned_profiles_csv` | bool | `true` | Write burned stream elevation profiles as CSV. Schema: `feature_id`, `flow_distance`, `x`, `y`, `dem_elevation`, `burned_elevation`. |
| `write_burned_profiles_points` | bool | `false` | Write burned stream elevation profiles as a vector point file. |
| `burned_profiles_write_vector_format` | string | `same_as_write_vector_format` | Burned-profile point override: `csv`, `flatgeobuf`, `fgb`, or `gpkg`. |

---

## snap_to_network: junction comparison

Used only when `stream_burn_mode: snap_to_network`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_junction_comparison` | bool | `true` | Write a point file pairing first-layer vector junctions with the nearest DEM-extracted junctions. Each pair appears twice: once for the vector junction (`is_dem_junction=0`) and once for the DEM junction (`is_dem_junction=1`). The `distance_m` attribute gives the Euclidean distance. |
| `junction_comparison_write_vector_format` | string | `same_as_write_vector_format` | Junction-comparison point override: `flatgeobuf`, `fgb`, `gpkg`, or `csv`. |
| `junction_tolerance_metres` | float | `50` | Junction-pair distance (m) below which a pair is considered correctly matched and excluded from repair candidates. |

---

## snap_to_network: junction gap repair

| Parameter | Type | Default | Description |
|---|---|---|---|
| `repair_junction_gaps` | bool | `false` | Attempt to repair DEM obstructions at mismatched first-layer junctions. For each mismatched pair, a least-cost terrain path is carved from the DEM junction toward the vector junction. Falls back to burning along the input vector where no terrain path is found. |
| `max_gap_search_pixels` | int | `60` | Search-window half-width in pixels. The Dijkstra search is bounded to a rectangle spanning both junction locations plus this many pixels on each side. |
| `gap_uphill_penalty` | float | `20` | Cost multiplier for uphill elevation steps during gap repair. Higher values force the path to cross obstructions at their lowest point. |
| `max_gap_climb_metres` | float | `10` | Maximum elevation the repair path may climb above the DEM junction (m). Paths exceeding this limit fall back to the vector-burn strategy. |
| `write_repaired_dem_raster` | bool | `false` | Write the gap-repaired DEM as a raster. Only used when `repair_junction_gaps` is true. |
| `write_repair_log` | bool | `true` | Write a CSV log of every junction repair. Columns: `repair_index`, `network_index`, `vector_x`, `vector_y`, `dem_x`, `dem_y`, `distance_m`, `fix_type`, `n_burned`. |

---

## snap_to_accumulation: D8 trace and LCP burn

Used when `stream_burn_mode: snap_to_accumulation`. This is the recommended workflow for accurate vector layers (e.g. OS VectorMap Local). The linework is first simplified into longer end-to-end segments, then D8 flow is traced along those segments; vector reaches whose traces stay within tolerance are accepted, and least-cost-path carving is applied on the DEM between anchor points derived from accepted vector reaches to close unresolved connectivity gaps.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `run_stream_burn` | bool | `true` | Run the burn stage. Set `false` to stop after the **initial stage** (fill/flow, clipped+repaired vector network, `_merged_lines.fgb`, source-seeded `_dem_network_lines.fgb`); no burning is done. |
| `write_merged_lines` | bool | `true` | Write `_merged_lines.fgb`: nearby end-to-end channel segments merged into longer lines before trace testing. Written in the initial stage. |
| `merge_snap_tolerance_metres` | float | `2` | Endpoint snap tolerance (m) used for that end-to-end merge. Keep **small** (smaller than the shortest repair bridge); too large can collapse nearby junction geometry and prevent sensible merges. |
| `segmented_trace_test` | bool | `false` | Diagnostic trace (no burning): from points every `segmented_trace_spacing_metres` along each merged line, trace D8 until the flow strays > `trace_max_divergence_metres` (or reaches the line end / outlet / sink). Writes `_segmented_traces.fgb`. Required for `run_lcp_burn`. |
| `segmented_trace_spacing_metres` | float | `20` | Spacing (m) of trace start points along each line. |
| `trace_max_divergence_metres` | float | `10` | Max distance (m) the D8 trace may stray from the vector before the segment is treated as diverged. Set above the typical vector-to-D8 offset so a channel merely beside the routed path is not flagged. |
| `segmented_trace_min_pixels` | int | `20` | A trace is **accepted** if it reaches the line end or stays within tolerance for more than this many cells. Accepted traces (overlaps merged) → `_accepted_segments.fgb`. |
| `segmented_trace_endpoint_fraction` | float | `0.5` | An accepted trace is trimmed from its downstream end until its last cell is within this fraction of `trace_max_divergence_metres` of the line. |
| `segmented_trace_spur_pixels` | int | `20` | During cleaning, a remaining trace is dropped as a spur if its unique (non-shared) cells number fewer than this. Cleaned linework → `_accepted_clean.fgb`. |
| `run_lcp_burn` | bool | `false` | Connect the cleaned accepted vector reaches into a continuous network by least-cost-path burning the gaps (source / intermediate / outlet). Requires `segmented_trace_test: true`. Writes `_burn_segments.fgb`, `_lcp_burn_paths.fgb`, and (with `write_burned_outputs`) the burned DEM plus re-extracted connected network. |
| `connect_tolerance_metres` | float | `0` | Gap tolerance (m): vector reaches whose ends are farther apart than this get a burn segment. Default `0` = **every** gap between consecutive accepted vector reaches is burned (no short gap is dropped). |
| `lcp_channel_alpha` | float | `4` | Drainage-area weighting for the least-cost path: higher values make the path prefer existing channels (high drainage area) over a shorter straight route. |
| `write_burned_outputs` | bool | `true` | Write the burned DEM and the network re-extracted from it. |

### Using OS VectorMap Local with `snap_to_accumulation`

Recommended sequence for a new user:

1. Repair the raw VML water layer with `lsdtt-water-line-repair` (see `docs/cli_reference/cli_reference_lsdtt_water_line_repair.md`).
2. Use `<repair_write_fname>_water_line_repaired.<ext>` as `stream_vector_path` in this CLI.
3. Run with:
   - `stream_burn_mode: snap_to_accumulation`
   - `segmented_trace_test: true`
   - `run_lcp_burn: true`
   - `write_burned_outputs: true` (if you want final conditioned outputs)
4. Inspect outputs in order:
   - `_accepted_clean.fgb` (vector reaches with accepted D8 traces retained)
   - `_burn_segments.fgb` and `_lcp_burn_paths.fgb` (where carving was applied)
   - `_lcp_burned.tif` (final conditioned DEM)
   - `_connected_network_lines.fgb` (fixed post-burn network consistent with the conditioned DEM)

### Inward carving of DEM-edge sources

Shared by `snap_to_network` and `snap_to_accumulation`. When a vector source sits on the DEM boundary and the DEM's D8 flow there drains back off the raster, the channel cannot be seeded inward; if the vector network confirms the source drains inward (its downstream junction is lower), the DEM is carved a few pixels along the vector path so D8 routes into the DEM.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `carve_inward_edge_sources` | bool | `true` | Enable inward carving at DEM-edge vector sources that the vector confirms drain inward. Runs before the LCP burn so the trace and the seeded network build on the inward-routed DEM. |
| `edge_source_carve_pixels` | int | `3` | Minimum pixels to carve inward along the vector path at edge sources. |
| `edge_source_carve_max_pixels` | int | `40` | Maximum pixels to carve inward; the carve extends adaptively up to this until the terminus reaches `edge_source_min_interior_metres` from the boundary. |
| `edge_source_min_interior_metres` | float | `10` | Target distance (m) from the nearest DEM boundary that the inward carve aims to reach. |

---

## QC / diagnostic outputs

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_vector_sources` | bool | `true` | Write vector source tip points. Includes an `is_outlet` attribute to distinguish the outlet from upstream sources. |
| `write_vector_outlet` | bool | `true` | Include the vector network outlet tip in the source-tip output. Only used when `write_vector_sources` is true. |
| `vector_sources_write_vector_format` | string | `same_as_write_vector_format` | Source-tip point override: `flatgeobuf`, `fgb`, `gpkg`, or `csv`. |

---

## DEM network outputs

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_dem_network_points` | bool | `false` | Write the DEM-extracted channel network as a point layer with `drainage_area` and `stream_order` attributes. |
| `write_dem_network_lines` | bool | `true` | Write the DEM-extracted channel network as a line layer. Each link between junctions (or source to junction) becomes a polyline feature. |
| `write_dem_network_junctions` | bool | `false` | Write DEM channel network junctions as a point layer with `junction_number`, `receiver_junction`, and `stream_order`. |
| `write_vector_format` | string | `flatgeobuf` | Baseline DEM-network format. `flatgeobuf`, `fgb`, or `gpkg`. |
| `calculate_spatial_index` | bool | `true` | Compute a spatial index on vector output files. Improves GIS query speed; buffers the entire layer in memory before writing. |

---

## D8-guided vector-network cleaning

| Parameter | Type | Default | Description |
|---|---|---|---|
| `threshold_contributing_pixels` | int | `1000` | Contributing-pixel threshold defining the D8 channel mask used to score loop arms. |

---

## Output files (always written)

| File | Description |
|---|---|
| `<write_fname>_filled.tif` | Pit-filled DEM (written when `raster_is_filled: false`). |
| `<write_fname>_vector_junctions.fgb` | Junctions of the clipped vector network with attributes `junction_id`, `is_outlet`, `donor_count`. |
| `<write_fname>_vector_edges.fgb` | Edges of the clipped vector network. |

---

## snap_to_accumulation / LCP burn output files

Written by the `snap_to_accumulation` workflow (depending on the flags above).

| File | Written when | Description |
|---|---|---|
| `<write_fname>_merged_lines.fgb` | `write_merged_lines` | Repaired vector segments merged end-to-end into longer continuous lines (without merging true confluences). |
| `<write_fname>_dem_network_lines.fgb` | initial stage | Source-seeded D8 network on the (pre-burn) filled DEM. |
| `<write_fname>_segmented_traces.fgb` | `segmented_trace_test` | All D8 traces from points along each line, with `stop_reason`. |
| `<write_fname>_accepted_segments.fgb` | `segmented_trace_test` | Accepted traces (overlaps merged), with `edge_id`, `strahler`. |
| `<write_fname>_accepted_clean.fgb` | `segmented_trace_test` | Accepted linework after spur/containment cleaning. |
| `<write_fname>_burn_segments.fgb` | `run_lcp_burn` | One feature per bridged gap (start→end, case, connected). |
| `<write_fname>_lcp_burn_paths.fgb` | `run_lcp_burn` | The carved least-cost paths. |
| `<write_fname>_lcp_burned.tif` | `run_lcp_burn` + `write_burned_outputs` | The final burned DEM. |
| `<write_fname>_connected_network_lines.fgb` | `run_lcp_burn` | The continuous network re-extracted from the burned DEM. |

---

## Example parameter file: `burn_in_place`

```
read_path:  /data/dems/
read_fname: my_dem
write_path: /data/output/
write_fname: my_burn_in_place

stream_vector_path: /data/vectors/channel_network.gpkg
stream_burn_mode:   burn_in_place

raster_is_filled: false
fill_min_slope:   0.0001

write_burned_outputs:      true
write_burned_profiles_csv: true
write_burned_profiles_points: false
```

## Example parameter file: `snap_to_network`

```
read_path:  /data/dems/
read_fname: my_dem
write_path: /data/output/
write_fname: my_snap_to_network

stream_vector_path: /data/vectors/channel_network.gpkg
stream_burn_mode:   snap_to_network

raster_is_filled: false
fill_min_slope:   0.0001

write_junction_comparison: true
junction_tolerance_metres: 50
repair_junction_gaps:      false

write_dem_network_lines:      true
write_dem_network_junctions:  true
write_vector_format:          gpkg
write_vector_sources:         true
calculate_spatial_index:      true
```

## Example parameter file: `snap_to_accumulation`

```
read_path:  /data/dems/
read_fname: my_lidar_dem
write_path: /data/output/
write_fname: my_snap_to_accum

# Prefer repaired VML lines from lsdtt-water-line-repair.
stream_vector_path: /data/vectors/vml_water_line_repaired.gpkg
stream_burn_mode:   snap_to_accumulation

raster_is_filled:        false
carve_before_fill:       true
breach_algorithm_choice: lindsay2016
fill_min_slope:          0.001

segmented_trace_test:       true
segmented_trace_spacing_metres:  20
trace_max_divergence_metres:     10
run_lcp_burn:               true
connect_tolerance_metres:   0

write_burned_outputs:       true
write_merged_lines:         true
merge_snap_tolerance_metres: 2
```

To inspect the starting point before any burning, set `run_stream_burn: false` (writes the merged vector lines and the source-seeded DEM network, no burn).

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

### Threading

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all available hardware threads. |

### Routing boundary conditions

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for D8 flow routing. | Four-character string: north, east, south, west. o = open, n = no-flux, p = periodic, b = forced outlet. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `raster_is_filled` | `bool` | `false` | `switch` | - | Assert that the input raster has ALREADY been pit-filled. | ASKS: has the input already been filled? Filling is REQUIRED for this analysis, so this switch does not choose whether to fill -- it only lets the tool skip filling a surface that is filled already. Setting it true for a DEM that still contains pits gives BAD routing, not different routing. Do not confuse this with fill_before_routing, which asks whether to fill at all. |

### DEM preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `carve_before_fill` | `bool` | `false` | `param` | - | Run depression breaching before pit filling. | Useful when the DEM has road or embankment obstructions across channels. |
| `breach_algorithm_choice` | `string` | `lindsay2016` | `param` | - | Breaching algorithm used when carve_before_fill is true. | Allowed values: lindsay2016 (unbounded), least_cost (supports bounded breaching via breach_max_cost / breach_max_dist_pixels). |
| `min_elevation` | `float` | `-9999` | `param` | - | Cells at or below this elevation (m) become nodata before any analysis. -9999 disables. | For removing the sea, which many DEMs store as valid cells at about 0 m, and nodata stored as 0 (a block of zeros beside real terrain is a cliff that captures flow paths). Removed cells that touch the DEM edge act as an open boundary, so a river reaching the coast drains off the DEM there. Negative values apply as given: where land lies below sea level (fens, polders, Death Valley -86 m, the Dead Sea -430 m) set it just below the lowest land, or below the sea surface of the DEM, so that flow is still routed over the sub-zero land. -9999 (or lower) disables. |

### Bounded least-cost breaching (breach_algorithm_choice = least_cost)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `breach_max_cost` | `float` | `0` | `param` | - | Cap on the total depth (m) carved along one least-cost breach path. 0 = unbounded. | The cost of a path is the depth each cell on it is lowered, summed (metres; with breach_minimize_dist true each step is weighted by its length in cells, so the units become metres x cells). Discriminates obstruction type by carve depth: shallow carves (lake spillways, small pits) stay below the cap and breach; deep carves through road/rail embankments exceed it and are left filled, so those crossings remain as gaps for the LCP burn (run_lcp_burn) to bridge under the vector line. Start ~2-5 (m). Only used when breach_algorithm_choice is least_cost. |
| `breach_max_dist_pixels` | `int` | `1000` | `param` | - | Maximum breach path length in cells for least-cost breaching. | Counted in steps from the pit (a diagonal step counts as one), not as a radius; 1 to 32767. Must be large enough to span the widest water body you want breached. Only used when breach_algorithm_choice is least_cost. |
| `breach_minimize_dist` | `bool` | `false` | `param` | - | Weight least-cost breach paths by their length. | false: a path costs the total depth carved (m), which is what breach_max_cost caps. true: each step is weighted by its length in cells, favouring shorter paths, and breach_max_cost is then in metres x cells. Only used when breach_algorithm_choice is least_cost. |
| `breach_edge_policy` | `string` | `outlet` | `param` | - | Whether a least-cost breach may drain off the raster edge. | outlet: the raster edge is an outlet, like a nodata cell. whitebox: the edge is not an outlet (as in WhiteboxTools), so a breach must reach an edge cell that is already low enough, which can carve deep cuts across divides near the edge. Only used when breach_algorithm_choice is least_cost. |
| `breach_flat_increment` | `float` | `0` | `param` | - | Per-step elevation drop (m) imposed along a breach channel. | 0 = auto: twice the float32 spacing at the largest absolute elevation in the DEM, the smallest drop that keeps the channel strictly downhill. Only used when breach_algorithm_choice is least_cost. |
| `breach_fill_remainder` | `bool` | `true` | `param` | - | Fill any depression that could not be breached within the cost/distance bounds (so flow still routes off it). | Keep true for bounded breaching: road-dammed depressions that exceed breach_max_cost are filled with fill_min_slope rather than carved. Only used when breach_algorithm_choice is least_cost. |
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope imposed during pit filling. | Used only when raster_is_filled is false. |

### Vector channel input

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `stream_vector_path` | `string` | - | `param` | - | Path to the input vector file containing the channel network polylines. | Required. Must be a polyline (linestring / multilinestring) layer in the same CRS as the DEM. Supported: GeoPackage (.gpkg), FlatGeobuf (.fgb), Shapefile (.shp), and any GDAL/OGR format. |

### Burn mode

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `vector_direction_algorithm_choice` | `string` | `dem_elevation` | `param` | - | How the vector network's flow direction is decided: dem_elevation or digitised. | dem_elevation: each vector junction drains through its steepest downhill edge by DEM elevation. Where the DEM rises along a channel (trees, an embankment, a culvert) an inflowing reach can be oriented backwards and take its junction's only downstream slot, cutting the reach below out of the network (counted as loop-labelled edges). digitised: drain through the edge whose digitised order leaves the junction, as written by lsdtt-fetch-river-lines, whose lines are oriented downstream and checked against their neighbours. A local rise along a reach is then followed, not reversed. Guard: a whole network whose outlets stand higher than the median of its channel heads is taken to be digitised upstream and is oriented by dem_elevation instead (reported in the log). |
| `stream_burn_mode` | `string` | `snap_to_network` | `param` | - | How the supplied vector network is used to condition the DEM. | burn_in_place: directly lower DEM pixels under the vector polylines to enforce a monotone-decreasing elevation profile (fast, no network extraction required beyond DEM preparation). snap_to_network: compare the supplied vector junctions against junctions extracted from the DEM flow network and produce a diagnostic point file showing where the two networks agree or diverge (first step toward guided carving). Best for less accurate vector layers (e.g. OS Rivers). snap_to_accumulation: for ACCURATE vector layers (e.g. OS VectorMap Local). Trace D8 along the vector (segmented_trace_test), accept the reaches that follow the LiDAR thalweg, and bridge the remaining obstruction gaps (roads/embankments) with least-cost-path carving (run_lcp_burn) to recover one continuous network. Preserves the LiDAR thalweg and writes the burned DEM plus the network extracted from it. burn_to_corridor: for a vector network that is the truth (e.g. OS Open Rivers prepared by lsdtt-fetch-river-lines, oriented downstream). The channel is the lowest path inside a corridor of vector_corridor_half_width_metres around each line (0 = the line's own cells), on the input DEM with its small pits filled (burn_profile_max_pit_pixels). Its elevations connect the minima of that DEM along the path (burn_profile_fit_algorithm_choice, burn_profile_max_cut_metres), so the burned channel follows the real valley floor, not a filled lake. The rest of the DEM is then carved (always least-cost breach, using the breach_* parameters) and filled to drain to the channel, and finally the whole DEM is conditioned as in the other modes (carve_before_fill, breach_algorithm_choice, fill_min_slope), so flow routing on the output follows the burned channel. |

### burn_to_corridor

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `vector_corridor_half_width_metres` | `float` | `0` | `param` | - | burn_to_corridor: how far the channel may be pulled away from the vector line into the DEM's valley [m]. | The channel is the least-cost path inside this distance of the line, on the input DEM with its small pits filled, where low ground is cheap (cost per step: 1 + height above the lowest corridor cell within this distance plus one cell). 0 (or under half a cell) = the line's own cells, i.e. burn in place, as a network. A few cells = follow the DEM valley floor but stay near the mapped river. Where a lower cell beside the channel would still take its flow (the valley lies farther away than this), that cell is raised into a bank. Each junction also gets a disc of max(this, 1.5 cells) so reaches can meet. |
| `burn_profile_max_pit_pixels` | `int` | `16` | `param` | - | burn_to_corridor: depressions of the input DEM with fewer cells than this are filled before the burned profile is sampled. | DEM noise makes small pits that are not real; a pure minima profile would follow them. Each depression (including nested ones, from the depression hierarchy) with fewer than this many cells is filled to its spill level, so noise pits inside a large real depression (e.g. the valley behind an embankment) are filled while that depression is not. 0 = no filling (pure minima). The count is in cells, so the same value covers a smaller area on a finer DEM. Used only to choose and sample the channel; the rest of the DEM is conditioned from the input DEM. |
| `burn_profile_max_cut_metres` | `float` | `2000` | `param` | - | burn_to_corridor: a minimum along the channel whose use would cut more than this below the DEM downstream is not used [m]. | Walking downstream, a cell lower than the last minimum used is a candidate minimum. The cut it causes is the cut below the DEM downstream of it with it, minus the cut there without it (and without the lower cells it leads into): the profile beyond a deep hole (a quarry, a cutting, a void) is held down until the DEM falls to the hole's level, so the hole causes that cut. Above this value the candidate is rejected and the profile passes over it (those cells are raised). The minimum just upstream of an obstruction is not blamed for cutting through it (the obstruction is cut either way); a path that still cuts more than this somewhere is reported. The default does nothing in practice. |
| `burn_profile_fit_algorithm_choice` | `string` | `linear` | `param` | - | burn_to_corridor: how the burned profile joins successive minima: linear or chi. | linear: a straight line in distance between minima. chi: a straight line in chi, chi = integral((chi_A0/A)^m_over_n ds), i.e. slope S = k_s A^-m_over_n with k_s set for each segment so the profile passes through both minima (concave up where the drainage area A grows). A is read along the channel from D8 on the DEM after a first linear pass, so it follows the burned network. Either way the profile passes through the minima, falls downstream and is never above the DEM except where a rejected minimum or a junction raises it. |
| `m_over_n` | `float` | `0.45` | `param` | - | Reference m/n (theta) for chi integration. | burn_to_corridor with burn_profile_fit_algorithm_choice chi. Must be &gt;= 0. |
| `chi_A0` | `float` | `1.0` | `param` | - | Reference area A0 for chi integration (m^2). | Used in chi = integral((A0/A)^(m/n) ds). burn_to_corridor with burn_profile_fit_algorithm_choice chi. |
| `vector_name_field` | `string` | - | `param` | - | burn_to_corridor: attribute of the vector lines holding the river name, copied to the outputs. | Empty = no names. The lines prepared by lsdtt-fetch-river-lines carry the name in 'name'. Written to the name column of _burned_profiles.csv and _corridor_channels. |
| `write_burn_depth_raster` | `bool` | `false` | `switch` | - | Write the change made to the DEM: original minus burned [m]. | Output file: &lt;write_fname&gt;_burn_depth.tif. Positive where the DEM was lowered, negative where it was raised (filled pits and hollows, banks, a channel raised to meet the one it joins), nodata where it is unchanged. Measured against the input DEM, not the filled one, so a cut through a filled hollow is not understated. Currently written by stream_burn_mode burn_to_corridor. |

### snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `trace_max_divergence_metres` | `float` | `10` | `param` | - | Max distance (m) the D8 flow trace may stray from the vector before the segment is burned. | Set above the typical offset between the vector lines and the LiDAR valley bottom, so a channel that merely sits beside the thalweg is not burned; only genuine divergences (roads, captures) are. |
| `merge_snap_tolerance_metres` | `float` | `2` | `param` | - | Endpoint snap tolerance (m) for merging degree-2 vector fragments. | Used by stream_burn_mode=snap_to_accumulation to write &lt;write_fname&gt;_merged_lines.fgb: fragmented vector lines whose endpoints fall within this distance and form a degree-2 junction (shared by exactly two fragments) are merged into one continuous line. The repaired vector layer fills road gaps with short bridge fragments, so endpoints are nearly coincident -- keep this SMALL (&lt; the shortest bridge fragment); too large collapses a short bridge and its two neighbours into one node (degree &gt;=3) and blocks the merge. Confluences (degree &gt;=3) and free ends (degree 1) remain breaks. |

### switch

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_stream_burn` | `bool` | `true` | - | - | Run the stream-burning stage. Set false to stop after the INITIAL stage: fill/flow, the clipped+repaired vector network, the merged vector lines (&lt;write_fname&gt;_merged_lines.fgb) and the source-seeded D8 network (&lt;write_fname&gt;_dem_network_lines.fgb) are written, but no burning is done. Useful for inspecting the starting point before any burn. | snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML) |
| `run_lcp_burn` | `bool` | `false` | - | - | After segmented_trace_test, connect the cleaned accepted reaches into a continuous network by least-cost-path burning the gaps (source / intermediate / outlet). Requires segmented_trace_test=true. Writes &lt;write_fname&gt;_burn_segments.fgb (start-&gt;end per gap), _lcp_burn_paths.fgb (the carved paths), _lcp_burn_raised.fgb (cells raised to block diversions), the burned DEM, and re-extracts the connected network. | snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML) |
| `write_merged_lines` | `bool` | `true` | - | - | Write &lt;write_fname&gt;_merged_lines.fgb: the repaired vector fragments merged across degree-2 junctions into continuous lines (see merge_snap_tolerance_metres). Written in the initial stage, independent of the burn. | snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML) |

### param

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `connect_tolerance_metres` | `float` | `0` | - | - | Gap tolerance (m) for run_lcp_burn: reaches whose ends are farther apart than this get a burn segment. Default 0 = EVERY gap between consecutive accepted reaches is burned (no short gap is ever dropped). | snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML) |
| `lcp_channel_alpha` | `float` | `4` | - | - | Drainage-area weighting for the run_lcp_burn least-cost path: higher = the path prefers existing channels (high drainage area) more strongly over a shorter straight route. | snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML) |
| `segmented_trace_test` | `bool` | `false` | - | - | Diagnostic (no burning): trace D8 from points every segmented_trace_spacing_metres along each merged+repaired vector line, stopping each trace where the flow strays &gt; trace_max_divergence_metres from the line (or reaches the line end / outlet / sink). Writes &lt;write_fname&gt;_segmented_traces.fgb (fields line_id, seg_id, start_s_m, max_divergence_m, stop_reason). Runs in the initial stage; pair with run_stream_burn=false to stop after it. | snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML) |
| `segmented_trace_spacing_metres` | `float` | `20` | - | - | Spacing (m) of trace start points along each line for segmented_trace_test. | snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML) |
| `segmented_trace_min_pixels` | `int` | `20` | - | - | For segmented_trace_test: a trace is ACCEPTED if it reaches the line end OR stays within trace_max_divergence_metres of the line for more than this many cells. The accepted D8 traces (overlaps merged) are written to &lt;write_fname&gt;_accepted_segments.fgb. | snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML) |
| `segmented_trace_spur_pixels` | `int` | `20` | - | - | For segmented_trace_test cleaning: after dropping contained traces, a remaining trace is removed as a SPUR if the cells unique to it (not shared with another kept trace) number fewer than this. The cleaned linework is written to &lt;write_fname&gt;_accepted_clean.fgb. | snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML) |
| `segmented_trace_endpoint_fraction` | `float` | `0.5` | - | - | For segmented_trace_test: an accepted trace is trimmed back from its downstream end until its last cell is within this FRACTION of trace_max_divergence_metres of the line. So a trace may wander up to the full divergence tolerance mid-way, but its endpoint sits close to the channel (default 0.5 = half the stop threshold). | snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML) |

### burn_in_place outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_burned_profiles_csv` | `bool` | `true` | `switch` | - | Write the burned stream elevation profiles as a CSV file. | burn_in_place: feature_id, flow_distance, x, y, dem_elevation, burned_elevation. burn_to_corridor: one row per channel cell, downstream order within each path (a chain of reaches from its head to where it joins another path or leaves the network): path_id, path_rank (fitting order, 0 = first and longest), reach_id, name, point_index, distance_m (along the path), flow_distance_m (to the network outlet along the burned network), x, y, z_original (the input DEM), z_cleaned (small pits filled), z_profile (the fitted profile), z_final (the output DEM), pit_filled, pit_id (the filled depression, -1 if none), is_record_low (a minimum the profile passes through), is_rejected_low (a minimum rejected by burn_profile_max_cut_metres), drainage_area_m2 and chi (chi fit only), joins_path_id and join_distance_m (where the path ends on another path, distance along that path), starts_on_path_id. The log reports any channel that rises downstream. |
| `write_burned_profiles_points` | `bool` | `false` | `switch` | - | Write the burned stream elevation profiles as a vector point file. | Only used when stream_burn_mode is burn_in_place. |
| `burned_profiles_write_vector_format` | `string` | `same_as_write_vector_format` | `param` | - | Output format for burned profile point output. | Only used when write_burned_profiles_points is true. Allowed values: csv, flatgeobuf, fgb, gpkg. |

### snap_to_network outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_junction_comparison` | `bool` | `true` | `switch` | - | Write a point file pairing first-layer junctions in the supplied vector network with the nearest junctions in the DEM-extracted network. | Only used when stream_burn_mode is snap_to_network. Each pair_id appears twice: once for the vector junction (is_dem_junction=0) and once for the DEM junction (is_dem_junction=1). Attribute distance_m gives the Euclidean distance between them. |
| `junction_comparison_write_vector_format` | `string` | `same_as_write_vector_format` | `param` | - | Output format for the junction comparison point file. | Allowed values: flatgeobuf, fgb, gpkg, csv. |
| `write_reach_divergence` | `bool` | `false` | `switch` | - | Diagnostic: compare each vector edge with the DEM channel walked between the matched DEM junctions. | For every vector edge whose both endpoints match a DEM junction, walks the DEM channel downstream and records geometric divergence metrics (detour ratio, lateral separation, gradient ratio) plus the reach status: ok / diverged / topology_break / unmatched_endpoint. Writes &lt;write_fname&gt;_reach_divergence_dem_reach + _vector_reach line layers and _reach_divergence.csv. Detection is purely geometric; mean channel score is logged for confirmation only (a road ditch scores like a channel). Only used when stream_burn_mode is snap_to_network. |
| `reach_divergence_separation_metres` | `float` | `25.0` | `param` | - | Lateral separation [m] above which a DEM channel pixel counts as off the vector edge. | Set above typical digitising offset so a vector that merely fails to sit on the DEM thalweg is not flagged. Drives divergence_length_m and the severity score. Only used when write_reach_divergence is true. |
| `reach_divergence_min_length_metres` | `float` | `100.0` | `param` | - | Sustained off-vector length [m] that flags a reach as diverged. | A reach is 'diverged' if its divergence_length_m exceeds this OR its detour ratio exceeds reach_divergence_detour_ratio. Only used when write_reach_divergence is true. |
| `reach_divergence_detour_ratio` | `float` | `1.5` | `param` | - | Detour ratio (DEM path length / vector length) that flags a reach as diverged. | Secondary trigger alongside reach_divergence_min_length_metres. Only used when write_reach_divergence is true. |
| `reach_divergence_match_radius_metres` | `float` | `50.0` | `param` | - | Max distance [m] for matching a vector junction to a DEM junction in the reach diagnostic. | A vector edge endpoint with no DEM junction within this radius makes the reach status 'unmatched_endpoint'. DEM confluences are offset from vector vertices, so this is generous by default (cf. junction_tolerance_metres in the repair). Only used when write_reach_divergence is true. |

### snap_to_network: junction gap repair

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `repair_junction_gaps` | `bool` | `false` | `param` | - | Attempt to repair DEM obstructions at mismatched first-layer junctions. | For each junction pair whose distance exceeds junction_tolerance_metres, runs a least-cost terrain path from the DEM junction to the vector junction and carves the DEM to allow correct D8 routing. Falls back to burning under the input vector where no terrain path can be found (tunnels). Only used when stream_burn_mode is snap_to_network. |
| `repair_reach_divergence` | `bool` | `false` | `param` | - | After topological junction repairs, repair reaches the junction matcher misses: a DEM channel captured by a road/ditch that detours a long way before rejoining the correct channel (status diverged or topology_break from the reach-divergence diagnostic). | Reuses the same channel-guided narrows repair as the junction gaps, anchored on the reach's downstream DEM confluence; the vector supplies topology only. Requires channel guidance; only used when repair_junction_gaps is true and stream_burn_mode is snap_to_network. |
| `junction_tolerance_metres` | `float` | `50` | `param` | - | Minimum junction-pair distance (m) to trigger gap repair. | Pairs closer than this are considered correctly matched and skipped. |
| `max_repair_gap_metres` | `float` | `100` | `param` | - | Maximum gap (m) the repair will attempt to carve. | A real road/embankment/culvert gap in a LiDAR DEM spans tens of metres. A much larger vector-to-DEM junction gap means the DEM flow is already correct and the vector network is a coarse local representation (e.g. a confluence whose branches genuinely reunite far downstream). Candidates beyond this are reported as topology mismatches and left unrepaired rather than carving a long spurious channel. Set 0 to disable the cap. |
| `max_gap_search_pixels` | `int` | `60` | `param` | - | Search-window half-width in pixels for gap repair. | The Dijkstra search is bounded to a rectangle that spans the vector and DEM junction locations expanded by this many pixels on each side. |
| `gap_uphill_penalty` | `float` | `20` | `param` | - | Cost multiplier applied to uphill elevation steps during gap repair. | Higher values force the repair path to cross obstructions at their lowest point; lower values allow more relaxed routing. |
| `max_gap_climb_metres` | `float` | `10` | `param` | - | Maximum elevation the repair path may climb above the DEM junction (m). | Paths exceeding this limit are rejected as likely tunnels or deep cuts and fall back to the vector-burn strategy. |
| `write_repaired_dem_raster` | `bool` | `false` | `switch` | - | Write the modified (gap-repaired) DEM as a raster. | Only used when repair_junction_gaps is true. |
| `write_burned_outputs` | `bool` | `true` | `switch` | - | Write the final burned DEM and the stream network extracted from it. | In snap_to_network mode, after all conditioning and gap repair, writes &lt;write_fname&gt;_burned.tif (the carved elevation surface) and &lt;write_fname&gt;_burned_network_lines.fgb plus _burned_network_junctions.fgb (the D8 stream network extracted from the burned DEM, for checking that flow paths route correctly under roads and other obstructions). |
| `write_channel_score_raster` | `bool` | `false` | `switch` | - | Write the channel candidate score raster used to guide gap repair. | Values in [0,1]: 1 = strong channel signal, 0 = off-channel. Useful for diagnosing why a repair path is or is not found. Only written when repair_junction_gaps is true and channel guidance is enabled. |
| `write_divergence_points` | `bool` | `false` | `switch` | - | Write the divergence point chosen for each repair attempt as a vector point file. | Shows exactly where the channel-guided LCP starts. Only written when repair_junction_gaps is true. |
| `write_repair_log` | `bool` | `true` | `switch` | - | Write a CSV log of every junction repair performed. | Only used when repair_junction_gaps is true. |
| `max_repair_phases` | `int` | `10` | `param` | - | Maximum number of candidate-detection/repair phases. | Each phase repairs at most one upstream-most candidate per incomplete basin, then re-routes the DEM and re-runs candidate detection. |
| `repair_min_slope` | `float` | `-1` | `param` | - | Minimum downhill slope imposed on carved repair paths (m/m). | Values &lt;= 0 default to fill_min_slope. Set 0 explicitly only if non-increasing flat burns are desired. |

### snap_to_network: inward carving of edge sources

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `carve_inward_edge_sources` | `bool` | `true` | `param` | - | Carve a few DEM pixels inward at vector source junctions on the DEM edge. | When a vector source sits on the DEM boundary and the DEM's D8 flow at that pixel drains back off the raster, the DEM channel cannot be seeded inward. If the vector network confirms the source drains inward (its downstream junction is lower), the DEM is carved a few pixels along the vector path so D8 routes into the DEM. Sources the vector confirms drain OFF the edge (downstream junction off-raster or higher) are left alone. Used by stream_burn_mode snap_to_network and snap_to_accumulation; for snap_to_accumulation it runs before the LCP burn so the segmented trace and the source-seeded network build on the inward-routed DEM. |
| `edge_source_carve_pixels` | `int` | `3` | `param` | - | Minimum pixels to carve inward along the vector path at edge sources. | The carve extends adaptively beyond this minimum (up to edge_source_carve_max_pixels) until the terminus reaches edge_source_min_interior_metres from the nearest DEM boundary. Only used when carve_inward_edge_sources is true. |
| `edge_source_carve_max_pixels` | `int` | `40` | `param` | - | Maximum pixels to carve inward along the vector path at edge sources. | Safety cap on the adaptive carve. The carve normally stops earlier, once the terminus reaches edge_source_min_interior_metres from the boundary; this cap only binds for channels that approach the edge at a very shallow angle. Only used when carve_inward_edge_sources is true. |
| `edge_source_min_interior_metres` | `float` | `10.0` | `param` | - | Perpendicular distance inward [m] the edge-source carve must reach. | The carve follows the downstream vector edge inward until the terminus pixel is at least this far from the nearest DEM boundary, so the seeded channel reaches genuine interior terrain rather than hugging the edge. A distance test is used rather than a D8-drains-off test because near the boundary the DEM flow runs along the edge and may never reach an outlet. Only used when carve_inward_edge_sources is true. |

### Channel-score guidance for gap repair

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `channel_score_alpha` | `float` | `2.0` | `param` | - | Strength of channel-score weighting in the repair Dijkstra cost function. | Cost per step = step_distance x (1 + alpha x (1 - channel_score)) + max(0, Deltaelev) x gap_uphill_penalty. alpha=0 disables channel guidance (terrain cost only). Only used when repair_junction_gaps is true. |
| `channel_lambda_min_metres` | `float` | `-1` | `param` | - | Minimum wavelength (m) for multiscale channel candidate detection. | Values &lt;= 0 default to 3 x DEM cell size. Only used when repair_junction_gaps is true. |
| `channel_lambda_max_metres` | `float` | `-1` | `param` | - | Maximum wavelength (m) for multiscale channel candidate detection. | Values &lt;= 0 default to 10 x DEM cell size. Only used when repair_junction_gaps is true. |
| `channel_n_scales` | `int` | `4` | `param` | - | Number of wavelength scales for channel candidate detection. | More scales improve robustness but increase compute time. Only used when repair_junction_gaps is true. |

### Output - QC/diagnostic

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vector_sources` | `bool` | `true` | `switch` | - | Write the vector network source tips as a point layer. | Useful for checking that the source detection found the correct upstream channel heads. The outlet tip is also included with an attribute flag to distinguish it. |
| `write_vector_outlet` | `bool` | `true` | `switch` | - | Include the vector network outlet tip in the source-tip output. | Only used when write_vector_sources is true. The outlet tip is labelled with is_outlet = 1. |
| `vector_sources_write_vector_format` | `string` | `same_as_write_vector_format` | `param` | - | Output format for the vector source tip points. | Allowed values: flatgeobuf, fgb, gpkg, csv. |

### Output - DEM network (seeded from vector sources)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_dem_network_points` | `bool` | `false` | `switch` | - | Write the DEM-extracted channel network as a point layer. | Each channel pixel becomes a feature with drainage_area and stream_order attributes. |
| `write_dem_network_lines` | `bool` | `true` | `switch` | - | Write the DEM-extracted channel network as a line layer. | Each channel link (between junctions or source and junction) becomes a polyline feature. |
| `write_dem_network_junctions` | `bool` | `false` | `switch` | - | Write the DEM-extracted channel network junctions as a point layer. | Includes junction_number, receiver_junction, and stream_order. |
| `calculate_spatial_index` | `bool` | `true` | `io` | - | Compute a spatial index when writing vector outputs. | Improves GIS spatial query speed at the cost of buffering the entire layer in memory before writing. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). CSV is available only through selectors for point/table outputs. Per-output overrides exist only where a container genuinely differs. Raster outputs are governed separately by write_output_format. |

### D8 cleaning of the clipped vector network (loops + dead ends)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `clean_vector_network_d8` | `bool` | `true` | `param` | - | Run D8-guided cleaning of the clipped vector network before junction finding. | Resolves loops (removes the arm farthest from a D8 channel) and dead-ending segments (links the downstream tip to the next segment along the D8 flow path). Writes a cleaned-network line layer. |
| `d8_dead_end_search_max_metres` | `float` | `750.0` | `param` | - | Maximum D8 walk distance (m) when reconnecting a dead-ending tip. | If the flow path does not reach another vector segment within this distance the dead end is left unresolved. Only used when clean_vector_network_d8 is true. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Interior nodata

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `interior_nodata_method` | `string` | `sink` | `param` | - | How nodata pixels inside the DEM are treated: sink or fill. | Interior nodata are nodata pixels not connected to the raster edge through other nodata (isolated voids, moulins, masked lakes). sink (the default, unchanged behaviour): they are outlets, so flow that reaches them ends there, as at the DEM edge; use this for real internal drainage (e.g. moulins on an ice sheet). fill: they are interpolated with nodata_fill_method_choice before any filling or routing, so flow passes across them; use this for data voids. Nodata connected to the edge (e.g. sea removed by min_elevation) is never filled and always acts as an open boundary. |
| `nodata_fill_method_choice` | `string` | `laplacian` | `param` | `interior_nodata_method=fill` | Method for filling interior nodata pixels. | laplacian: smooth Dirichlet interpolation (recommended). iterative_weighted_average: fast, well-suited to thin voids. Used when interior_nodata_method is fill. |
| `nodata_fill_max_iter` | `int` | `10000` | `param` | `interior_nodata_method=fill` | Maximum iterations for iterative nodata fill methods. | Used when interior_nodata_method is fill. |

<!-- END GENERATED CLI PARAMETERS -->
