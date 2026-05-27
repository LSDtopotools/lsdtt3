# CLI reference: lsdtt-marsh-topography

High-level marsh and tidal-channel analysis workflow. Three linked stages can be run independently or together: marsh-platform extraction, tidal-channel extraction, and tidal-skeleton metrics. At least one stage must be enabled. The DEM must be in a projected CRS (e.g. UTM).

## Pipeline overview

Each stage is independently toggled. Outputs from earlier stages are passed forward automatically.

| Stage | Toggle | Description |
|---|---|---|
| 1 | `run_marsh_platform` | Extract marsh platform via scarp detection and iterative platform growth |
| 2 | `run_tidal_channels` | Extract tidal channels via multi-scale median residual thresholding |
| 3 | `run_tidal_skeleton_metrics` | Compute section, station, and strip metrics from the tidal skeleton |

Stages 2 and 3 can use outputs from stage 1 (marsh mask in strip metrics). Stage 3 requires a channel mask and skeleton from stage 2, or pre-computed rasters supplied via `read_channel_mask_fname` and `read_skeleton_fname`.

## Usage

```
lsdtt-marsh-topography <path> <parameter_file>
lsdtt-marsh-topography <parameter_file>         # path defaults to ./
lsdtt-marsh-topography -h                       # print help and exit
lsdtt-marsh-topography -v                       # print version and exit
lsdtt-marsh-topography -c                       # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input raster filename (without extension). |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |

---

## Workflow selection

| Parameter | Type | Default | Description |
|---|---|---|---|
| `run_marsh_platform` | bool | `true` | Run marsh-platform extraction. |
| `run_tidal_channels` | bool | `true` | Run tidal-channel extraction. |
| `run_tidal_skeleton_metrics` | bool | `true` | Run tidal skeleton metrics. Requires a channel mask and skeleton (from the tidal-channel stage, or supplied via `read_channel_mask_fname` and `read_skeleton_fname`). |

---

## Threading and output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` uses all available hardware threads. |
| `write_output_format` | string | `cog` | Raster output format for all raster outputs. Options: `cog`, `gtiff`, `same_as_input`. |
| `polygon_vector_output_format` | string | `flatgeobuf` | Output format for polygon vector outputs. Options: `flatgeobuf` / `fgb`, `gpkg`. |
| `line_vector_output_format` | string | `flatgeobuf` | Output format for line vector outputs. Options: `flatgeobuf` / `fgb`, `gpkg`. |
| `point_vector_output_format` | string | `gpkg` | Output format for point vector outputs. Options: `csv`, `flatgeobuf` / `fgb`, `gpkg`. |

---

## Optional input rasters

These allow stages to be skipped or to use pre-computed inputs.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `marsh_use_supplied_slope` | bool | `false` | Use an externally supplied slope raster for marsh-platform extraction. Requires `read_slope_fname`. |
| `read_slope_fname` | string | `NULL` | Slope raster filename. Used only when `marsh_use_supplied_slope: true`. |
| `read_channel_mask_fname` | string | `NULL` | Tidal-channel mask raster filename. Used when `run_tidal_channels: false` but skeleton metrics are requested. |
| `read_skeleton_fname` | string | `NULL` | Skeleton raster filename. Used when `run_tidal_channels: false`, or when an externally prepared skeleton should replace the extracted one. |
| `read_outlet_mask_fname` | string | `NULL` | Outlet mask raster filename. Required when `tidal_outlet_mode: outlet_mask`. |

---

## Marsh-platform parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `marsh_histogram_n_bins` | int | `100` | Number of bins in marsh search/cutoff histograms. |
| `marsh_search_space_histogram_derivative_threshold` | float | `-2.0` | Histogram-derivative threshold defining the marsh search space. |
| `marsh_max_scarp_growth_iterations` | int | `50` | Maximum iterative scarp-growth passes. |
| `marsh_scarp_clean_relief_fraction` | float | `0.85` | Local-relief fraction for scarp cleaning. Higher values keep only scarps in stronger relief. |
| `marsh_scarp_clean_min_neighbour_pixels` | int | `8` | Minimum nearby scarp pixels required after relief cleaning. Removes isolated artefacts. |
| `marsh_max_platform_growth_iterations` | int | `100` | Maximum iterative marsh-platform growth passes. |
| `marsh_cutoff_histogram_run_length_bins` | int | `8` | Run length below mean histogram frequency needed to define the cutoff elevation. |
| `marsh_preserve_scarps_in_final_platform` | bool | `true` | Reinsert cleaned scarps into the final marsh-platform mask. Matches legacy LSDMarshPlatform behaviour. |
| `marsh_fill_high_gaps` | bool | `true` | Apply the legacy high-gap fill stage. |
| `marsh_fill_pools` | bool | `true` | Apply the legacy pool-filling stage. |
| `marsh_fill_small_holes` | bool | `true` | Fill isolated one-pixel holes in the final marsh-platform mask. |
| `marsh_apply_scarp_closing` | bool | `false` | Apply binary closing to the cleaned marsh scarp mask. |
| `marsh_scarp_closing_radius_metres` | float | `0.0` | Binary closing radius for the marsh scarp mask (m). Used only when `marsh_apply_scarp_closing: true`. |
| `marsh_label_use_scarp_barriers` | bool | `false` | Treat scarps as barriers when labelling marsh-platform components. |
| `marsh_label_use_d4_connectivity` | bool | `false` | Use D4 connectivity when labelling marsh-platform components. |
| `marsh_min_platform_component_pixels` | int | `1` | Minimum component size for the final marsh-platform mask. Set > 1 to remove isolated fragments. |

### Marsh outputs

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_marsh_search_space` | bool | `true` | `<write_fname>_marsh_search_space.tif` |
| `write_marsh_search_crossover` | bool | `false` | `<write_fname>_marsh_search_crossover.tif` |
| `write_marsh_scarp_seed_mask` | bool | `false` | `<write_fname>_marsh_scarp_seed_mask.tif` |
| `write_marsh_scarp_order_raster` | bool | `false` | `<write_fname>_marsh_scarp_order.tif` |
| `write_marsh_scarp_mask` | bool | `true` | `<write_fname>_marsh_scarp_mask.tif` |
| `write_marsh_platform_seed_mask` | bool | `false` | `<write_fname>_marsh_platform_seed_mask.tif` |
| `write_marsh_platform_provisional` | bool | `false` | `<write_fname>_marsh_platform_provisional.tif` |
| `write_marsh_platform_mask` | bool | `true` | `<write_fname>_marsh_platform_mask.tif` |
| `write_marsh_platform_labels` | bool | `true` | `<write_fname>_marsh_platform_labels.tif` |
| `write_marsh_polygons` | bool | `false` | `<write_fname>_marsh_platforms.<ext>` — vectorised platform components. |

---

## Tidal-channel parameters

### Channel detection

| Parameter | Type | Default | Description |
|---|---|---|---|
| `tidal_use_flood_mask` | bool | `false` | Use simple elevation-threshold channel mask instead of multi-scale median. When true, pixels at or below `tidal_flood_elevation_metres` are classified as channel. |
| `tidal_flood_elevation_metres` | float | `0.0` | Elevation threshold for the flood-mask method (m). Used only when `tidal_use_flood_mask: true`. |
| `tidal_median_window_radii_metres` | string | `5,10,20` | Comma-delimited median-window radii (m) for multi-scale channel extraction. One residual raster per scale. |
| `tidal_residual_thresholds_metres` | string | `0.1` | Comma-delimited residual thresholds (m). One value is used for all scales; otherwise one per scale. |
| `tidal_scale_mode` | string | `any_scale_pass` | How per-scale masks are combined. Options: `any_scale_pass`, `all_scales_pass`. |

### Channel mask cleanup

| Parameter | Type | Default | Description |
|---|---|---|---|
| `tidal_apply_binary_closing` | bool | `false` | Apply binary closing to the channel mask. |
| `tidal_binary_closing_radius_metres` | float | `0.0` | Closing radius for the channel mask (m). |
| `tidal_fill_holes` | bool | `false` | Fill holes in the channel mask. |
| `tidal_min_channel_component_pixels` | int | `1` | Minimum connected-component size for the channel mask. |

### Skeleton extraction

| Parameter | Type | Default | Description |
|---|---|---|---|
| `tidal_extract_skeleton` | bool | `true` | Extract a 1-pixel skeleton from the tidal-channel mask. |
| `tidal_skeleton_spur_max_pixels` | int | `10` | Maximum spur length removed from the tidal skeleton. |
| `tidal_skeleton_terminal_loop_max_pixels` | int | `10` | Maximum small terminal loop length removed from the tidal skeleton. |
| `tidal_skeleton_min_component_pixels` | int | `1` | Minimum connected-component size for the tidal skeleton. |

### Skeleton reconnection (optional)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `tidal_connect_skeleton_segments` | bool | `false` | Attempt to reconnect nearby skeleton endpoints. |
| `tidal_connect_tangent_fit_length_pixels` | float | `6.0` | Tangent support length for endpoint reconnection. |
| `tidal_connect_tangent_min_support_pixels` | int | `2` | Minimum support pixels for endpoint tangent estimation. |
| `tidal_connect_max_gap_pixels` | float | `20.0` | Maximum gap distance for reconnection candidates. |
| `tidal_connect_max_facing_angle_degrees` | float | `35.0` | Maximum endpoint-facing angle for reconnection candidates. |
| `tidal_connect_max_heading_misalignment_degrees` | float | `45.0` | Maximum tangent-heading mismatch for reconnection candidates. |
| `tidal_connect_max_candidate_score` | float | `1.0` | Maximum accepted candidate score for reconnection. |
| `tidal_connect_route_corridor_half_width_pixels` | float | `12.0` | Half-width of the least-cost routing corridor for accepted skeleton bridges. |

### UPL

| Parameter | Type | Default | Description |
|---|---|---|---|
| `tidal_compute_upl` | bool | `true` | Compute unchannelled path length (Euclidean distance to nearest channel pixel). |

### Tidal outputs

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_tidal_median_rasters` | bool | `false` | `<write_fname>_tidal_median_<i>.tif` for each scale. |
| `write_tidal_residual_rasters` | bool | `false` | `<write_fname>_tidal_residual_<i>.tif` for each scale. |
| `write_tidal_threshold_masks` | bool | `false` | `<write_fname>_tidal_threshold_mask_<i>.tif` for each scale. |
| `write_tidal_channel_mask` | bool | `true` | `<write_fname>_tidal_channel_mask.tif` |
| `write_tidal_upl` | bool | `true` | `<write_fname>_tidal_upl.tif` |
| `write_tidal_channel_skeleton` | bool | `true` | `<write_fname>_tidal_channel_skeleton.tif` |
| `write_tidal_connected_skeleton` | bool | `true` | `<write_fname>_tidal_connected_skeleton.tif` |
| `write_tidal_channel_polygons` | bool | `false` | `<write_fname>_tidal_channel_polygons.<ext>` — vectorised channel components. |
| `polygonise_min_component_pixels` | int | `1` | Minimum component size when polygonising binary masks. |

---

## Tidal skeleton metrics parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `tidal_outlet_mode` | string | `edge_nodes` | How outlet nodes are chosen. Options: `edge_nodes`, `endpoint_nodes`, `outlet_mask`, `lowest_endpoint`. |
| `tidal_station_spacing_metres` | float | `5.0` | Station spacing along tidal skeleton sections (m). |
| `tidal_compute_cross_sections` | bool | `true` | Compute width, depth, and area at tidal skeleton stations. |
| `tidal_tangent_support_points` | int | `2` | Number of neighbouring stations used to estimate local section tangent. |
| `tidal_cross_section_half_length_metres` | float | `25.0` | Half-length of tidal cross sections (m). |
| `tidal_cross_section_sample_spacing_metres` | float | `1.0` | Sampling spacing along tidal cross sections (m). |
| `tidal_max_width_ratio` | float | `3.0` | Maximum section-length to local-width ratio. Sections exceeding this are rejected. |
| `tidal_compute_station_strip_metrics` | bool | `true` | Assign each DEM pixel to its nearest tidal station and compute strip metrics (local and upstream strip area, optional strip-mean values). |
| `tidal_strip_include_dem` | bool | `true` | Include DEM elevation in station strip summaries. |
| `tidal_strip_include_upl` | bool | `false` | Include UPL in station strip summaries. Available only when the tidal-channel stage computed UPL. |
| `tidal_strip_include_marsh_mask` | bool | `false` | Include marsh-platform mask values in station strip summaries. Available only when the marsh-platform stage was run. |

### Tidal skeleton metric outputs

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_tidal_node_id_raster` | bool | `true` | `<write_fname>_tidal_node_id.tif` |
| `write_tidal_section_id_raster` | bool | `true` | `<write_fname>_tidal_section_id.tif` |
| `write_tidal_section_strahler_raster` | bool | `true` | `<write_fname>_tidal_section_strahler.tif` |
| `write_tidal_section_hack_raster` | bool | `true` | `<write_fname>_tidal_section_hack.tif` |
| `write_tidal_nearest_station_id_raster` | bool | `true` | `<write_fname>_tidal_nearest_station_id.tif` |
| `write_tidal_station_downstream_distance_raster` | bool | `true` | `<write_fname>_tidal_station_downstream_distance.tif` |
| `write_tidal_sections_vector` | bool | `true` | `<write_fname>_tidal_sections.<ext>` — oriented tidal section lines. |
| `write_tidal_stations_vector` | bool | `true` | `<write_fname>_tidal_stations.<ext>` — station points with metrics. |
| `write_tidal_stations_csv` | bool | `true` | `<write_fname>_tidal_stations.csv` — station metrics table. |

---

## Example parameter file: full workflow

```
read_path:  /data/dems/
read_fname: marsh_dem
write_path: /data/output/
write_fname: marsh_output

n_threads:           0
write_output_format: cog

run_marsh_platform:          true
run_tidal_channels:          true
run_tidal_skeleton_metrics:  true

# --- Marsh platform ---
marsh_histogram_n_bins:                          100
marsh_search_space_histogram_derivative_threshold: -2.0
marsh_max_scarp_growth_iterations:               50
marsh_scarp_clean_relief_fraction:               0.85
marsh_scarp_clean_min_neighbour_pixels:          8
marsh_max_platform_growth_iterations:            100

write_marsh_search_space:  true
write_marsh_scarp_mask:    true
write_marsh_platform_mask: true
write_marsh_platform_labels: true

# --- Tidal channels ---
tidal_median_window_radii_metres:   5,10,20
tidal_residual_thresholds_metres:   0.1
tidal_scale_mode:                   any_scale_pass
tidal_extract_skeleton:             true
tidal_compute_upl:                  true

write_tidal_channel_mask:       true
write_tidal_upl:                true
write_tidal_channel_skeleton:   true
write_tidal_connected_skeleton: true

# --- Skeleton metrics ---
tidal_outlet_mode:                   edge_nodes
tidal_station_spacing_metres:        5.0
tidal_compute_cross_sections:        true
tidal_cross_section_half_length_metres: 25.0
tidal_strip_include_dem:             true
tidal_strip_include_marsh_mask:      true

write_tidal_sections_vector:   true
write_tidal_stations_vector:   true
write_tidal_stations_csv:      true
```

## Example parameter file: tidal channels and metrics only

```
read_path:  /data/dems/
read_fname: tidal_dem
write_path: /data/output/
write_fname: tidal_channels

run_marsh_platform:          false
run_tidal_channels:          true
run_tidal_skeleton_metrics:  true

tidal_median_window_radii_metres:    10,20
tidal_residual_thresholds_metres:    0.15
tidal_scale_mode:                    any_scale_pass
tidal_extract_skeleton:              true
tidal_skeleton_spur_max_pixels:      15
tidal_compute_upl:                   true

tidal_outlet_mode:                   lowest_endpoint
tidal_station_spacing_metres:        5.0
tidal_compute_cross_sections:        true
tidal_strip_include_dem:             true
tidal_strip_include_upl:             true
```

## Example parameter file: skeleton metrics from pre-computed inputs

```
read_path:  /data/dems/
read_fname: dem
write_path: /data/output/
write_fname: skeleton_metrics

run_marsh_platform:          false
run_tidal_channels:          false
run_tidal_skeleton_metrics:  true

read_channel_mask_fname: /data/output/dem_tidal_channel_mask.tif
read_skeleton_fname:     /data/output/dem_tidal_connected_skeleton.tif

tidal_outlet_mode:            edge_nodes
tidal_station_spacing_metres: 10.0
tidal_compute_cross_sections: true
```
