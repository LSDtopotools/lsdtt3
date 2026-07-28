# CLI reference: lsdtt-stream-burning

Conditions a DEM using a user-supplied vector channel network. The tool builds two parallel network views (vector topology vs D8-routed DEM channels), then applies one of three `stream_burn_mode` workflows:

- `burn_in_place` — directly lower DEM cells under the vector polylines to enforce downhill channel profiles.
- `snap_to_network` — compare vector and DEM junction topology, then optionally repair mismatches/gaps.
- `snap_to_accumulation` — for accurate vectors (for example OS VectorMap Local), trace D8 flow along the vector and, with `do_lcp_burn: true`, least-cost carve only the unresolved connectivity gaps between accepted vector reaches to recover a continuous routed network.

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
| 11 | **Mode: `snap_to_accumulation`** — join adjacent water-line segments that connect end-to-end (simple pass-through joins, not confluences) into longer lines (`_merged_lines.fgb`), then run D8 trace acceptance along those lines (`_accepted_segments.fgb` -> `_accepted_clean.fgb`), and if `do_lcp_burn` is enabled carve only the unresolved obstruction gaps; optionally write the carved DEM (`_lcp_burned.tif`) and the post-burn connected network (`_connected_network_lines.fgb`) |

The `snap_to_accumulation` workflow has two stages, controlled independently:

- **Initial stage** (always runs): breach/fill → flow → clipped+repaired vector network → merged vector lines → source-seeded D8 network (`_dem_network_lines.fgb`). Setting `do_stream_burn: false` stops here, so no DEM conditioning is applied.
- **LCP burn stage** (`do_lcp_burn: true`, requires `segmented_trace_test: true`): bridge only the unresolved gaps between accepted vector reaches with least-cost-path carving, then (if requested) write the conditioned DEM (`_lcp_burned.tif`) and the fixed connected network extracted from it (`_connected_network_lines.fgb`).

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
| `burned_profiles_output_format` | string | `flatgeobuf` | Output format for burned profile points. `flatgeobuf`, `fgb`, or `gpkg`. |

---

## snap_to_network: junction comparison

Used only when `stream_burn_mode: snap_to_network`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_junction_comparison` | bool | `true` | Write a point file pairing first-layer vector junctions with the nearest DEM-extracted junctions. Each pair appears twice: once for the vector junction (`is_dem_junction=0`) and once for the DEM junction (`is_dem_junction=1`). The `distance_m` attribute gives the Euclidean distance. |
| `junction_comparison_output_format` | string | `flatgeobuf` | Output format for the junction comparison point file. `flatgeobuf`, `fgb`, `gpkg`, or `csv`. |
| `junction_tolerance_m` | float | `50` | Junction-pair distance (m) below which a pair is considered correctly matched and excluded from repair candidates. |

---

## snap_to_network: junction gap repair

| Parameter | Type | Default | Description |
|---|---|---|---|
| `repair_junction_gaps` | bool | `false` | Attempt to repair DEM obstructions at mismatched first-layer junctions. For each mismatched pair, a least-cost terrain path is carved from the DEM junction toward the vector junction. Falls back to burning along the input vector where no terrain path is found. |
| `max_gap_search_px` | int | `60` | Search-window half-width in pixels. The Dijkstra search is bounded to a rectangle spanning both junction locations plus this many pixels on each side. |
| `gap_uphill_penalty` | float | `20` | Cost multiplier for uphill elevation steps during gap repair. Higher values force the path to cross obstructions at their lowest point. |
| `max_gap_climb_m` | float | `10` | Maximum elevation the repair path may climb above the DEM junction (m). Paths exceeding this limit fall back to the vector-burn strategy. |
| `write_repaired_dem` | bool | `false` | Write the gap-repaired DEM as a raster. Only used when `repair_junction_gaps` is true. |
| `write_repair_log` | bool | `true` | Write a CSV log of every junction repair. Columns: `repair_index`, `network_index`, `vector_x`, `vector_y`, `dem_x`, `dem_y`, `distance_m`, `fix_type`, `n_burned`. |

---

## snap_to_accumulation: D8 trace and LCP burn

Used when `stream_burn_mode: snap_to_accumulation`. This is the recommended workflow for accurate vector layers (e.g. OS VectorMap Local). The linework is first simplified into longer end-to-end segments, then D8 flow is traced along those segments; vector reaches whose traces stay within tolerance are accepted, and least-cost-path carving is applied on the DEM between anchor points derived from accepted vector reaches to close unresolved connectivity gaps.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `do_stream_burn` | bool | `true` | Run the burn stage. Set `false` to stop after the **initial stage** (fill/flow, clipped+repaired vector network, `_merged_lines.fgb`, source-seeded `_dem_network_lines.fgb`); no burning is done. |
| `write_merged_lines` | bool | `true` | Write `_merged_lines.fgb`: nearby end-to-end channel segments merged into longer lines before trace testing. Written in the initial stage. |
| `merge_snap_tol_m` | float | `2` | Endpoint snap tolerance (m) used for that end-to-end merge. Keep **small** (smaller than the shortest repair bridge); too large can collapse nearby junction geometry and prevent sensible merges. |
| `segmented_trace_test` | bool | `false` | Diagnostic trace (no burning): from points every `segmented_trace_spacing_m` along each merged line, trace D8 until the flow strays > `trace_max_divergence_m` (or reaches the line end / outlet / sink). Writes `_segmented_traces.fgb`. Required for `do_lcp_burn`. |
| `segmented_trace_spacing_m` | float | `20` | Spacing (m) of trace start points along each line. |
| `trace_max_divergence_m` | float | `10` | Max distance (m) the D8 trace may stray from the vector before the segment is treated as diverged. Set above the typical vector-to-D8 offset so a channel merely beside the routed path is not flagged. |
| `segmented_trace_min_pixels` | int | `20` | A trace is **accepted** if it reaches the line end or stays within tolerance for more than this many cells. Accepted traces (overlaps merged) → `_accepted_segments.fgb`. |
| `segmented_trace_endpoint_fraction` | float | `0.5` | An accepted trace is trimmed from its downstream end until its last cell is within this fraction of `trace_max_divergence_m` of the line. |
| `segmented_trace_spur_pixels` | int | `20` | During cleaning, a remaining trace is dropped as a spur if its unique (non-shared) cells number fewer than this. Cleaned linework → `_accepted_clean.fgb`. |
| `do_lcp_burn` | bool | `false` | Connect the cleaned accepted vector reaches into a continuous network by least-cost-path burning the gaps (source / intermediate / outlet). Requires `segmented_trace_test: true`. Writes `_burn_segments.fgb`, `_lcp_burn_paths.fgb`, and (with `write_burned_outputs`) the burned DEM plus re-extracted connected network. |
| `connect_tol_m` | float | `0` | Gap tolerance (m): vector reaches whose ends are farther apart than this get a burn segment. Default `0` = **every** gap between consecutive accepted vector reaches is burned (no short gap is dropped). |
| `lcp_channel_alpha` | float | `4` | Drainage-area weighting for the least-cost path: higher values make the path prefer existing channels (high drainage area) over a shorter straight route. |
| `write_burned_outputs` | bool | `true` | Write the burned DEM and the network re-extracted from it. |

### Using OS VectorMap Local with `snap_to_accumulation`

Recommended sequence for a new user:

1. Repair the raw VML water layer with `lsdtt-water-line-repair` (see `docs/cli_reference/cli_reference_lsdtt_water_line_repair.md`).
2. Use `<repair_write_fname>_water_line_repaired.<ext>` as `stream_vector_path` in this CLI.
3. Run with:
   - `stream_burn_mode: snap_to_accumulation`
   - `segmented_trace_test: true`
   - `do_lcp_burn: true`
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
| `vector_sources_output_format` | string | `flatgeobuf` | Output format for vector source tip points. `flatgeobuf`, `fgb`, `gpkg`, or `csv`. |

---

## DEM network outputs

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_dem_network_points` | bool | `false` | Write the DEM-extracted channel network as a point layer with `drainage_area` and `stream_order` attributes. |
| `write_dem_network_lines` | bool | `true` | Write the DEM-extracted channel network as a line layer. Each link between junctions (or source to junction) becomes a polyline feature. |
| `write_dem_network_junctions` | bool | `false` | Write DEM channel network junctions as a point layer with `junction_number`, `receiver_junction`, and `stream_order`. |
| `dem_network_output_format` | string | `flatgeobuf` | Output format for DEM network points, lines, and junctions. `flatgeobuf`, `fgb`, or `gpkg`. Points also accept `csv` (lat/lon in EPSG:4326). |
| `calculate_spatial_index` | bool | `true` | Compute a spatial index on vector output files. Improves GIS query speed; buffers the entire layer in memory before writing. |

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
| `<write_fname>_burn_segments.fgb` | `do_lcp_burn` | One feature per bridged gap (start→end, case, connected). |
| `<write_fname>_lcp_burn_paths.fgb` | `do_lcp_burn` | The carved least-cost paths. |
| `<write_fname>_lcp_burned.tif` | `do_lcp_burn` + `write_burned_outputs` | The final burned DEM. |
| `<write_fname>_connected_network_lines.fgb` | `do_lcp_burn` | The continuous network re-extracted from the burned DEM. |

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
junction_tolerance_m:      50
repair_junction_gaps:      false

write_dem_network_lines:      true
write_dem_network_junctions:  true
dem_network_output_format:    gpkg
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
segmented_trace_spacing_m:  20
trace_max_divergence_m:     10
do_lcp_burn:                true
connect_tol_m:              0

write_burned_outputs:       true
write_merged_lines:         true
merge_snap_tol_m:           2
```

To inspect the starting point before any burning, set `do_stream_burn: false` (writes the merged vector lines and the source-seeded DEM network, no burn).
