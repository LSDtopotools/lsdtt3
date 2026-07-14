# CLI reference: lsdtt-valley-centreline

Identifies a broad per-channel valley mask and extracts a medial-axis valley centreline for each source channel. The tool can also export reach-adaptive HAND-break diagnostics, regularly spaced valley cross sections, and optional cross-section-derived products.

Unlike `lsdtt-flat-surface-finder`, this tool is intentionally broad: it aims to map "the valley" around each channel and one centreline per channel, not detailed floodplain/terrace label classes.

## Usage

```bash
lsdtt-valley-centreline <path> <parameter_file>
lsdtt-valley-centreline <parameter_file>      # path defaults to ./
lsdtt-valley-centreline -h                    # print help and exit
lsdtt-valley-centreline -v                    # print version and exit
lsdtt-valley-centreline -c                    # print citation and exit
```

## Parameter file format

```text
# Lines beginning with # are comments.
key: value
```

---

## Pipeline overview

| Step | Description |
|---|---|
| 1 | Load DEM and pit-fill (`fill_min_slope_for_routing`) for robust routing |
| 2 | Build D8 flow object (`boundary_conditions`) |
| 3 | Load and map source/outlet points to raster pixels |
| 4 | Resolve channel scope (source snapping, optional outlet truncation, optional auto source/outlet role detection) |
| 5 | For each channel: build a wide corridor and compute Euclidean HAND above that single channel |
| 6 | Detect HAND break globally and by along-channel reach; classify valley pixels below break |
| 7 | Clean valley mask (opening/channel-adjacency/detour checks) |
| 8 | Route centreline as least-cost path that maximizes distance to valley walls (medial-axis behavior) |
| 9 | Fit pinned monotonic baseline along centreline |
| 10 | Optionally extract cross sections, feature picks, and hydraulic diagnostics |
| 11 | Write rasters/tables/vectors (`.tif`, `.csv`, `.fgb`) |

---

## Required and core I/O parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the DEM. |
| `read_fname` | string | parameter-file stem | Input raster filename (with or without extension). |
| `write_path` | string | same as `read_path` | Output directory. |
| `write_fname` | string | same as `read_fname` | Output stem used for all products. |
| `source_points_path` | string | `` | **Required.** Source point file. One centreline is produced per resolved source. |
| `write_output_format` | string | `gtiff` | Raster format: `cog`, `gtiff`, `same_as_input`. |

Important: this CLI currently requires `source_points_path`; threshold-only all-channel mode is not supported here.

---

## Adding sources and endpoints

### Source and outlet files

Use either:
- a **source-only** file, or
- a source file plus an explicit **outlet** file, or
- a **single combined source file** with `auto_detect_outlet_points: true` (default).

Point files can be CSV (`latitude`/`longitude` in EPSG:4326 by default) or OGR vectors (`.gpkg`, `.fgb`, etc.).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `source_points_path` | string | `` | Source points file (required). |
| `source_points_input_format` | string | `auto_detect` | `auto_detect`, `csv_latitude_longitude_epsg4326`, `ogr_vector_layer`. |
| `source_points_layer_name` | string | `` | Optional vector layer name. |
| `source_points_latitude_field_name` | string | `latitude` | CSV latitude field. |
| `source_points_longitude_field_name` | string | `longitude` | CSV longitude field. |
| `outlet_points_path` | string | `` | Optional outlet points file to truncate channels. |
| `outlet_points_input_format` | string | `auto_detect` | Same choices as source input format. |
| `outlet_points_layer_name` | string | `` | Optional vector layer name for outlet file. |
| `outlet_points_latitude_field_name` | string | `latitude` | CSV latitude field for outlets. |
| `outlet_points_longitude_field_name` | string | `longitude` | CSV longitude field for outlets. |
| `auto_detect_outlet_points` | bool | `true` | If no explicit outlet file is given, classify supplied points: downstream points become outlets. |

### Snapping and scope controls

| Parameter | Type | Default | Description |
|---|---|---|---|
| `channel_threshold_contributing_pixels` | int | `1000` | Drainage threshold used to define channel cells for source snapping. |
| `snap_sources_to_channel_threshold` | bool | `true` | Snap source points onto threshold channel cells. |
| `source_snap_max_radius_pixels` | int | `25` | Max expanding-square snap search radius. |
| `outlet_snap_max_distance_metres` | float | `0` | Max raw-outlet to traced-path snap distance (`0` = unlimited). |

---

## Valley extraction controls

| Parameter | Type | Default | Description |
|---|---|---|---|
| `valley_corridor_width_metres` | float | `5000` | Corridor half-width around channel used for HAND populations and candidate valley pixels. |
| `valley_fallback_hand_metres` | float | `10` | Fallback HAND break if no clear break is detected. |
| `valley_force_hand_metres` | float | `0` | Force a constant HAND threshold for valley extraction (`>0` enables; bypasses break detection). |
| `valley_min_hand_metres` | float | `2` | Floor on detected break to avoid channel-hugging degenerate valleys. |
| `valley_max_hand_metres` | float | `0` | Optional cap on break (`0` = no cap). |
| `valley_break_reach_window_metres` | float | `0` | Reach window length for adaptive breaks (`0` auto, `<0` single global break). |
| `refine_valley_with_baseline` | bool | `true` | Re-test valley mask against nearest-centreline baseline to remove stale-reference wedges. |
| `valley_dilation_pixels` | int | `4` | Cleanup dilation/opening strength (`0` disables dilation). |
| `valley_max_detour_ratio` | float | `3.0` | Prune pixels with long geodesic detours from channel (`0` disables). |

---

## Reach HAND statistics and break detection

The valley threshold is derived from HAND populations in two levels:
- **global per channel** (all corridor pixels for that channel), and
- **per reach window** (corridor pixels grouped by nearest channel node window).

The core statistic path is:
1. Compute Euclidean HAND above the single target channel.
2. Build log10(HAND) histogram (128 bins; low HAND floor clamp).
3. Smooth histogram with a 2-pass moving average.
4. Segment modes at local minima.
5. Merge weak/noisy mode splits based on dip depth and mode mass.
6. Take the first surviving dip after the first surviving mode as the valley break.
7. For reach windows without their own break: inherit nearest window break, then global break, then fallback.

Controls:

| Parameter | Type | Default | Description |
|---|---|---|---|
| `valley_wall_peak_fraction` | float | `0.02` | Minimum mode mass/peak floor (keeps small valley-floor mode detectable next to high wall mode). |
| `valley_break_method` | string | `fixed` | `fixed` or `relief_adaptive`. |
| `valley_break_fraction` | float | `0.6` | Fixed dip-depth acceptance fraction. Used only in `fixed`. |
| `valley_break_fraction_min` | float | `0.5` | `relief_adaptive`: strict end for low-relief reaches. |
| `valley_break_fraction_max` | float | `0.7` | `relief_adaptive`: permissive end for high-relief reaches. |
| `valley_break_relief_window_metres` | float | `300` | Radius for local-relief raster driving adaptive fraction. |

`relief_adaptive` uses a saturating relief-to-fraction mapping per reach, then applies the same relative dip-depth logic with that local fraction.

Diagnostics:
- `*_valley_breaks.csv` (always written): applied break per window, own/inherited flag, plus applied fraction and relief.
- `write_valley_hand_histogram: true`: `*_valley_hand_hist.csv` and `*_valley_hand_window_hist.csv`.
- `write_valley_break_metric: true`: `*_valley_break_metric.csv` with smoothed decision metric, peaks/boundaries, and chosen break bin.
- `write_valley_reach_id_raster: true`: `*_valley_reach_id_ch<sid>.tif` showing which pixels fed each reach histogram.

---

## Centreline controls and outputs

Centreline extraction now uses least-cost medial-axis routing inside the valley mask. Some legacy trough parameters remain exposed for compatibility.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `min_centreline_nodes` | int | `10` | Minimum acceptable routed centreline length. |
| `trough_scaling_factor` | float | `0.1` | Legacy trough parameter (retained in interface). |
| `trough_along_channel_slope` | float | `0.02` | Legacy trough parameter (retained in interface). |
| `centreline_loops` | int | `5` | Legacy trough iteration parameter (retained in interface). |
| `min_bank_elevation_window_radius_metres` | float | `200` | Legacy trough bank filter radius (retained in interface). |
| `trough_nodata_buffer_pixels` | int | `2` | Dilates valley hole before wall-distance field. |
| `trough_fill_min_slope` | float | `0.0001` | Legacy trough-fill slope parameter (retained in interface). |

Primary centreline outputs:
- `write_valley_centreline: true` -> `*_valley_centreline.fgb`
  - fields: `source_id`, `n_nodes`, `valley_hand_break_m`, `routing_ok`
- `write_valley_centreline_nodes: true` -> `*_valley_centreline_nodes.csv`
  - columns: `source_id`, `node_seq`, `x`, `y`, `elevation`, `trough_elevation`, `baseline_elevation`, `is_pin`, `distance_m`

---

## Valley mask outputs

| Parameter | Type | Default | Output |
|---|---|---|---|
| `write_valley_mask_raster` | bool | `true` | `*_valley_mask.tif` (UInt32, value `source_id + 1`) and `*_valley_mask_ch<sid>.tif` (Byte 0/1 per channel). |
| `write_valley_reach_id_raster` | bool | `true` | `*_valley_reach_id_ch<sid>.tif` (Int16 reach window IDs in corridor). |

Optional context rasters:
- `write_hillshade_raster` -> `*_hillshade.tif`
- `write_multidirectional_hillshade_raster` -> `*_multidirectional_hillshade.tif`

---

## Cross sections

Enable with `write_valley_cross_sections: true` (default).

### Geometry and spacing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `cross_section_spacing_metres` | float | `500` | Along-centreline station spacing. |
| `cross_section_length_multiplier` | float | `1.5` | Total section length = multiplier x widest bounded valley width on channel. |
| `cross_section_tangent_window_nodes` | int | `8` | Tangent smoothing window for orthogonal direction. |
| `cross_section_min_half_length_metres` | float | `0` | Optional floor for half-length. |
| `cross_section_max_half_length_metres` | float | `5000` | Hard cap for side search and half-length. |

### Overlap control and curved sections

| Parameter | Type | Default | Description |
|---|---|---|---|
| `clip_overlapping_sections` | bool | `true` | Straight-section non-overlap via local Voronoi/bisector clipping. |
| `neighbour_clip_window` | int | `2` | Number of neighboring stations checked on each side for clipping. |
| `min_section_reach_fraction` | float | `0.35` | Drop sections clipped below this fraction of nominal length. |
| `curved_cross_sections` | bool | `false` | Use boundary-fitted curved sections (streamwise isocontours). |
| `curved_wall_clearance_fraction` | float | `0.10` | Curved mode minimum wall clearance margin. |
| `curved_streamwise_smoothing_passes` | int | `40` | Curved mode smoothing passes for streamwise field. |
| `curved_coarse_spacing_metres` | float | `700` | Curved mode coarse-to-fine station placement spacing. |

### Cross-section outputs

- `*_valley_cross_sections.fgb`: section lines with `source_id`, `station_seq`, `along_distance_m`, `valley_width_m`, endpoint coords.
- `*_valley_cross_section_sides.fgb`: left/right wall points.
- `*_valley_cross_section_boundaries.fgb`: `wall_left`, `wall_right`, and `hydraulic_edge` boundary points.
- `*_valley_cross_section_samples.csv`: sampled profile rows (`offset_m`, `elevation`, `in_valley`, optional feature rasters).

---

## Cross-section feature and hydraulic products

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_cross_section_features` | bool | `true` | Enables feature extraction and writes `*_valley_cross_section_features.csv`. |
| `cross_section_curvature_window_metres` | float | `60` | Polyfit window for maximum-curvature raster used in feature picks. |
| `cross_section_feature_prominence_fraction` | float | `0.15` | Minimum prominence for curvature/hillshade peak reporting. |
| `write_cross_section_hydraulics` | bool | `false` | Writes BankfullMapper-style stage table `*_valley_cross_section_hydraulics.csv`. |
| `write_cross_section_informed_valley_mask` | bool | `false` | Builds `*_cross_section_informed_valley_mask.tif` (requires feature path). |

Additional feature/hydraulic outputs:
- `*_valley_cross_section_features.csv`
- `*_valley_cross_section_hydraulic_breaks.csv`
- `*_cross_section_informed_valley_mask.tif` and `*_cross_section_informed_valley_mask_ch<sid>.tif`

---

## Threading and routing options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | OpenMP threads (`0` uses all available). |
| `boundary_conditions` | string | `oooo` | Flow-routing boundary conditions (N,E,S,W). |
| `fill_min_slope_for_routing` | float | `0.0001` | Minimum slope when pit-filling DEM used for routing. |
| `multidirectional_hillshade_n_directions` | int | `8` | Number of azimuths for multidirectional hillshade. |

---

## Example parameter files

### 1) Single channel from one source point

```text
read_path:  /data/royalton/
read_fname: ry_utm.tif
write_path: /data/royalton/
write_fname: ry_valley_single

source_points_path: /data/royalton/channel_sources.gpkg
source_points_input_format: ogr_vector_layer
source_points_layer_name: channel_sources

auto_detect_outlet_points: false
snap_sources_to_channel_threshold: true
channel_threshold_contributing_pixels: 1000

valley_break_method: relief_adaptive
valley_corridor_width_metres: 5000
valley_break_reach_window_metres: 0
valley_min_hand_metres: 2
valley_fallback_hand_metres: 10

write_valley_centreline: true
write_valley_centreline_nodes: true
write_valley_mask_raster: true
write_valley_cross_sections: false
```

### 2) Source + explicit endpoint bounding

```text
read_path:  /data/mayoco/
read_fname: mayoco_COP30_utm.tif
write_path: /data/mayoco/
write_fname: mayoco_valley_bounded

source_points_path: /data/mayoco/channel_sources.gpkg
source_points_input_format: ogr_vector_layer
source_points_layer_name: sources

outlet_points_path: /data/mayoco/channel_outlets.gpkg
outlet_points_input_format: ogr_vector_layer
outlet_points_layer_name: outlets

auto_detect_outlet_points: false
outlet_snap_max_distance_metres: 200

valley_break_method: fixed
valley_break_fraction: 0.6
valley_max_hand_metres: 25

write_valley_centreline: true
write_valley_centreline_nodes: true
write_valley_hand_histogram: true
write_valley_break_metric: true
```

### 3) Combined point file with automatic outlet detection

```text
read_path:  /data/salinas/
read_fname: seco_utm.tif
write_path: /data/salinas/
write_fname: salinas_valley_autorole

source_points_path: /data/salinas/sources_for_fp.gpkg
source_points_input_format: ogr_vector_layer

auto_detect_outlet_points: true
snap_sources_to_channel_threshold: true
source_snap_max_radius_pixels: 30

valley_break_method: relief_adaptive
valley_break_fraction_min: 0.5
valley_break_fraction_max: 0.7
valley_break_relief_window_metres: 300

write_valley_mask_raster: true
write_valley_reach_id_raster: true
write_valley_hand_histogram: true
write_valley_break_metric: true
```

### 4) Dense cross sections, curved mode, and hydraulic outputs

```text
read_path:  /data/royalton/
read_fname: ry_utm.tif
write_path: /data/royalton/
write_fname: ry_valley_xs_curved

source_points_path: /data/royalton/sources_for_fp.gpkg
source_points_input_format: ogr_vector_layer
auto_detect_outlet_points: true

write_valley_cross_sections: true
cross_section_spacing_metres: 150
cross_section_length_multiplier: 1.5
cross_section_max_half_length_metres: 5000

curved_cross_sections: true
curved_coarse_spacing_metres: 700
curved_streamwise_smoothing_passes: 40
curved_wall_clearance_fraction: 0.10

clip_overlapping_sections: true
neighbour_clip_window: 2
min_section_reach_fraction: 0.35

write_cross_section_features: true
cross_section_curvature_window_metres: 60
cross_section_feature_prominence_fraction: 0.15
write_cross_section_hydraulics: true
write_cross_section_informed_valley_mask: true
```
