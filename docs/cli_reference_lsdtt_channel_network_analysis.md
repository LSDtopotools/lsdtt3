# CLI reference: lsdtt-channel-network-analysis

Builds a D8 channel and junction network from a DEM and then runs one or more channel-network analyses: channel-link bearing, junction-angle geometry, and Strahler/Tokunaga ordering.

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
| 8 | Strahler / Tokunaga ordering and output |

Steps 6 and 7 operate on the **selected junction subset** (see junction selection below). Step 8 (Strahler / Tokunaga) always operates on the **full extracted network**.

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input raster filename. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |

There is no `write_output_format` for this CLI — all outputs are vector files.

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
| `source_pixels_threshold` | int | `1000` | Minimum upstream contributing pixels required for a pixel to be a channel source. Reduce to extract denser networks. Must be > 0. |

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

## Step 8: Strahler / Tokunaga analysis

These analyses operate on the **full extracted network**, not the outlet-selected subset.

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
| `write_tokunaga_summary_data` | bool | `false` | `<write_fname>_tokunaga_summary.csv` — Tokunaga pair counts, offset statistics, and self-similarity fit summary for the full network. |

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
source_pixels_threshold:  1000
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
source_pixels_threshold:  500

# These analyses use the full network; no junction selection needed.
write_strahler_channels:        true
strahler_channels_output_format: flatgeobuf
write_tokunaga_channels:         true
write_tokunaga_summary_data:     true
```
