# CLI reference: lsdtt-stream-burning

Conditions a DEM using a user-supplied vector channel network. The tool builds two parallel channel network representations — one from the vector file, one extracted from D8 flow routing — and, depending on `stream_burn_mode`, either burns the vector network directly into the DEM (`burn_in_place`), diagnoses/repairs where the two networks diverge (`snap_to_network`), or — for accurate vector layers such as OS VectorMap Local — traces the LiDAR thalweg along the vector and bridges the obstruction gaps (roads/embankments) with least-cost-path carving so the burned DEM routes one continuous network (`snap_to_accumulation` + `do_lcp_burn`).

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
| 11 | **Mode: `snap_to_accumulation`** — merge degree-2 vector fragments (`_merged_lines.fgb`); optional inward-carve of DEM-edge sources; with `segmented_trace_test`, trace D8 along the vector and accept the reaches that follow the LiDAR thalweg (`_accepted_segments.fgb`, cleaned to `_accepted_clean.fgb`); with `do_lcp_burn`, least-cost-path carve the remaining gaps and re-extract the connected network |

The `snap_to_accumulation` workflow has two stages, controlled independently:

- **Initial stage** (always runs): breach/fill → flow → clipped+repaired vector network → merged vector lines → source-seeded D8 network (`_dem_network_lines.fgb`). Setting `do_stream_burn: false` stops here, for inspecting the starting point before any burn.
- **LCP burn stage** (`do_lcp_burn: true`, requires `segmented_trace_test: true`): bridge the gaps between accepted reaches with least-cost-path carving, write the burned DEM, and re-extract the continuous connected network.

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
| `stream_burn_mode` | string | `snap_to_network` | How the vector network conditions the DEM. `burn_in_place`: directly lowers DEM pixels under the vector polylines to enforce monotone-decreasing elevation profiles. `snap_to_network`: compares vector junctions against DEM-extracted junctions and produces diagnostic outputs showing where the two networks agree or diverge (best for less accurate layers, e.g. OS Rivers). `snap_to_accumulation`: for accurate vector layers (e.g. OS VectorMap Local) — traces the LiDAR thalweg along the vector and bridges obstruction gaps with least-cost-path carving (see the LCP-burn section below). |

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

## snap_to_accumulation: thalweg trace and LCP burn

Used when `stream_burn_mode: snap_to_accumulation`. This is the recommended workflow for accurate vector layers (e.g. OS VectorMap Local). The vector fragments are merged across degree-2 junctions, the D8 flow is traced along the merged lines, the reaches that follow the LiDAR thalweg are accepted, and the remaining obstruction gaps (roads, embankments) are bridged with least-cost-path carving.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `do_stream_burn` | bool | `true` | Run the burn stage. Set `false` to stop after the **initial stage** (fill/flow, clipped+repaired vector network, `_merged_lines.fgb`, source-seeded `_dem_network_lines.fgb`); no burning is done. |
| `write_merged_lines` | bool | `true` | Write `_merged_lines.fgb`: the repaired vector fragments merged across degree-2 junctions into continuous lines. Written in the initial stage. |
| `merge_snap_tol_m` | float | `2` | Endpoint snap tolerance (m) for merging degree-2 fragments. Keep **small** (smaller than the shortest bridge fragment); too large collapses a bridge plus its two neighbours into a degree-≥3 node and blocks the merge. |
| `segmented_trace_test` | bool | `false` | Diagnostic trace (no burning): from points every `segmented_trace_spacing_m` along each merged line, trace D8 until the flow strays > `trace_max_divergence_m` (or reaches the line end / outlet / sink). Writes `_segmented_traces.fgb`. Required for `do_lcp_burn`. |
| `segmented_trace_spacing_m` | float | `20` | Spacing (m) of trace start points along each line. |
| `trace_max_divergence_m` | float | `10` | Max distance (m) the trace may stray from the vector before the segment is treated as diverged. Set above the typical vector-to-thalweg offset so a channel merely beside the thalweg is not flagged. |
| `segmented_trace_min_pixels` | int | `20` | A trace is **accepted** if it reaches the line end or stays within tolerance for more than this many cells. Accepted traces (overlaps merged) → `_accepted_segments.fgb`. |
| `segmented_trace_endpoint_fraction` | float | `0.5` | An accepted trace is trimmed from its downstream end until its last cell is within this fraction of `trace_max_divergence_m` of the line. |
| `segmented_trace_spur_pixels` | int | `20` | During cleaning, a remaining trace is dropped as a spur if its unique (non-shared) cells number fewer than this. Cleaned linework → `_accepted_clean.fgb`. |
| `do_lcp_burn` | bool | `false` | Connect the cleaned accepted reaches into a continuous network by least-cost-path burning the gaps (source / intermediate / outlet). Requires `segmented_trace_test: true`. Writes `_burn_segments.fgb`, `_lcp_burn_paths.fgb`, `_lcp_burn_raised.fgb`, the burned DEM, and re-extracts the connected network. |
| `connect_tol_m` | float | `0` | Gap tolerance (m): reaches whose ends are farther apart than this get a burn segment. Default `0` = **every** gap between consecutive accepted reaches is burned (no short gap is dropped). |
| `lcp_channel_alpha` | float | `4` | Drainage-area weighting for the least-cost path: higher values make the path prefer existing channels (high drainage area) over a shorter straight route. |
| `write_burned_outputs` | bool | `true` | Write the burned DEM and the network re-extracted from it. |

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
| `<write_fname>_merged_lines.fgb` | `write_merged_lines` | Repaired vector fragments merged across degree-2 junctions. |
| `<write_fname>_dem_network_lines.fgb` | initial stage | Source-seeded D8 network on the (pre-burn) filled DEM. |
| `<write_fname>_segmented_traces.fgb` | `segmented_trace_test` | All D8 traces from points along each line, with `stop_reason`. |
| `<write_fname>_accepted_segments.fgb` | `segmented_trace_test` | Accepted traces (overlaps merged), with `edge_id`, `strahler`. |
| `<write_fname>_accepted_clean.fgb` | `segmented_trace_test` | Accepted linework after spur/containment cleaning. |
| `<write_fname>_burn_segments.fgb` | `do_lcp_burn` | One feature per bridged gap (start→end, case, connected). |
| `<write_fname>_lcp_burn_paths.fgb` | `do_lcp_burn` | The carved least-cost paths. |
| `<write_fname>_lcp_burn_raised.fgb` | `do_lcp_burn` | Cells raised to block competing diversions. |
| `<write_fname>_lcp_burned.tif` | `do_lcp_burn` + `write_burned_outputs` | The final burned DEM. |
| `<write_fname>_connected_network_lines.fgb` | `do_lcp_burn` | The continuous network re-extracted from the burned DEM. |

---

## Example parameter file

```
read_path:  /data/dems/
read_fname: my_dem
write_path: /data/output/
write_fname: my_dem_burned

stream_vector_path: /data/vectors/channel_network.gpkg

# DEM preparation
carve_before_fill: true
fill_min_slope:    0.0001

# Compare networks; identify mismatched junctions
stream_burn_mode:          snap_to_network
junction_tolerance_m:      50
write_junction_comparison: true

# DEM network outputs
write_dem_network_lines:   true
write_dem_network_junctions: true
dem_network_output_format: gpkg

# QC points
write_vector_sources:      true
calculate_spatial_index:   true
```

### Example: snap_to_accumulation with LCP burn (accurate vector layer)

```
read_path:  /data/dems/
read_fname: my_lidar_dem
write_path: /data/output/
write_fname: my_burned

stream_vector_path: /data/vectors/vectormap_local_water.fgb

# DEM preparation (breach lakes/pits, leave road embankments as gaps)
raster_is_filled:        false
carve_before_fill:       true
breach_algorithm_choice: lindsay2016
fill_min_slope:          0.001

# Trace the thalweg along the vector, then LCP-burn the gaps
stream_burn_mode:        snap_to_accumulation
segmented_trace_test:    true
segmented_trace_spacing_m: 20
trace_max_divergence_m:  10
do_lcp_burn:             true
connect_tol_m:           0
write_burned_outputs:    true
write_merged_lines:      true
merge_snap_tol_m:        2
```

To inspect the starting point before any burning, set `do_stream_burn: false` (writes the merged vector lines and the source-seeded DEM network, no burn).
