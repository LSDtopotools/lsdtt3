# CLI reference: lsdtt-stream-burning

Conditions a DEM using a user-supplied vector channel network. The tool builds two parallel channel network representations — one from the vector file, one extracted from D8 flow routing — and either burns the vector network directly into the DEM or identifies and optionally repairs locations where the two networks diverge.

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
| `stream_burn_mode` | string | `snap_to_network` | How the vector network conditions the DEM. `burn_in_place`: directly lowers DEM pixels under the vector polylines to enforce monotone-decreasing elevation profiles. `snap_to_network`: compares vector junctions against DEM-extracted junctions and produces diagnostic outputs showing where the two networks agree or diverge. |

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
