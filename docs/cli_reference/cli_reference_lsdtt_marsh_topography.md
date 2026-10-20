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
| `write_vector_format` | string | `flatgeobuf` | Output format for polygon and line vector outputs. Options: `flatgeobuf` / `fgb`, `gpkg`. |
| `marsh_points_write_vector_format` | string | `gpkg` | Point-output override. Options: `csv`, `flatgeobuf` / `fgb`, `gpkg`. |

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
| `write_marsh_search_space_raster` | bool | `true` | `<write_fname>_marsh_search_space.tif` |
| `write_marsh_search_crossover_raster` | bool | `false` | `<write_fname>_marsh_search_crossover.tif` |
| `write_marsh_scarp_seed_mask_raster` | bool | `false` | `<write_fname>_marsh_scarp_seed_mask.tif` |
| `write_marsh_scarp_order_raster` | bool | `false` | `<write_fname>_marsh_scarp_order.tif` |
| `write_marsh_scarp_mask_raster` | bool | `true` | `<write_fname>_marsh_scarp_mask.tif` |
| `write_marsh_platform_seed_mask_raster` | bool | `false` | `<write_fname>_marsh_platform_seed_mask.tif` |
| `write_marsh_platform_provisional` | bool | `false` | `<write_fname>_marsh_platform_provisional.tif` |
| `write_marsh_platform_mask_raster` | bool | `true` | `<write_fname>_marsh_platform_mask.tif` |
| `write_marsh_platform_labels_raster` | bool | `true` | `<write_fname>_marsh_platform_labels.tif` |
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
| `write_tidal_channel_mask_raster` | bool | `true` | `<write_fname>_tidal_channel_mask.tif` |
| `write_tidal_upl_raster` | bool | `true` | `<write_fname>_tidal_upl.tif` |
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

write_marsh_search_space_raster:  true
write_marsh_scarp_mask_raster:    true
write_marsh_platform_mask_raster: true
write_marsh_platform_labels_raster: true

# --- Tidal channels ---
tidal_median_window_radii_metres:   5,10,20
tidal_residual_thresholds_metres:   0.1
tidal_scale_mode:                   any_scale_pass
tidal_extract_skeleton:             true
tidal_compute_upl:                  true

write_tidal_channel_mask_raster: true
write_tidal_upl_raster:          true
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

### -

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_marsh_platform` | `bool` | `true` | `switch` | - | Run marsh-platform extraction. | If true, lsdtt_marsh_platform_cpu is run on the DEM. Outputs include search-space, scarp, and marsh-platform rasters. |
| `run_tidal_channels` | `bool` | `true` | `switch` | - | Run tidal-channel extraction. | If true, lsdtt_tidal_channel_cpu is run on the DEM. Outputs include channel mask, optional UPL, and optional skeleton. |
| `run_tidal_skeleton_metrics` | `bool` | `true` | `switch` | - | Run tidal-skeleton metrics. | If true, lsdtt_tidal_skeleton_metrics_cpu is run using a channel mask and a skeleton. These can come from the tidal-channel stage or from supplied rasters. |
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 = use all available hardware threads. |
| `write_output_format` | `string` | `cog` | `io` | - | Raster output format. | Allowed values: cog, gtiff, same_as_input. |
| `marsh_use_supplied_slope` | `bool` | `false` | `param` | - | Use an externally supplied slope raster for marsh-platform extraction. | If false, slope is computed internally inside lsdtt_marsh_platform_cpu. |
| `read_slope_fname` | `string` | `NULL` | `param` | - | Optional slope raster filename. | Used only when marsh_use_supplied_slope is true. |
| `read_channel_mask_fname` | `string` | `NULL` | `param` | - | Optional tidal-channel mask raster filename. | Used when run_tidal_channels is false but tidal-skeleton metrics still need a channel mask. |
| `read_skeleton_fname` | `string` | `NULL` | `param` | - | Optional skeleton raster filename. | Used when run_tidal_channels is false, or when an externally prepared skeleton should be used instead of the extracted one. |
| `read_outlet_mask_fname` | `string` | `NULL` | `param` | - | Optional outlet mask raster filename. | Required only when tidal_outlet_mode is set to outlet_mask. |
| `marsh_histogram_n_bins` | `int` | `100` | `param` | - | Number of bins used in marsh search/cutoff histograms. | Used in the search-space thresholding and final cutoff-elevation identification stages. |
| `marsh_search_space_histogram_derivative_threshold` | `float` | `-2.0` | `param` | - | Histogram-derivative threshold used to define the marsh search space. | This is the legacy LSDMarshPlatform control on the DEM-relative-relief times slope-relative-relief crossover raster. |
| `marsh_max_scarp_growth_iterations` | `int` | `50` | `param` | - | Maximum number of iterative scarp-growth passes. | Controls continuation of marsh scarp lines from seed maxima. |
| `marsh_scarp_clean_relief_fraction` | `float` | `0.85` | `param` | - | Local-relief fraction used when cleaning marsh scarps. | Higher values keep only scarps in stronger local relief settings. |
| `marsh_scarp_clean_min_neighbour_pixels` | `int` | `8` | `param` | - | Minimum number of nearby scarp pixels required after relief cleaning. | Used to remove isolated scarp artefacts. |
| `marsh_max_platform_growth_iterations` | `int` | `100` | `param` | - | Maximum number of iterative marsh-platform growth passes. | Controls outward growth from scarp-seeded platform pixels. |
| `marsh_cutoff_histogram_run_length_bins` | `int` | `8` | `param` | - | Run length below the mean histogram frequency needed to define the cutoff elevation. | This is the legacy marsh-platform elevation cutoff control. |
| `marsh_preserve_scarps_in_final_platform` | `bool` | `true` | `param` | - | Reinsert cleaned scarps into the final marsh-platform mask. | This matches the legacy LSDMarshPlatform behaviour. |
| `marsh_fill_high_gaps` | `bool` | `true` | `param` | - | Apply the legacy high-gap fill stage. | This fills isolated high-elevation gaps inside the platform domain. |
| `marsh_fill_pools` | `bool` | `true` | `param` | - | Apply the legacy pool-filling stage. | This fills enclosed low areas inside the evolving platform mask. |
| `marsh_fill_small_holes` | `bool` | `true` | `param` | - | Fill isolated one-pixel holes in the final marsh-platform mask. | Useful for minor raster cleanup. |
| `marsh_apply_scarp_closing` | `bool` | `false` | `param` | - | Apply binary closing to the cleaned marsh scarp mask. | Useful for bridging tiny gaps in otherwise continuous scarps before platform seeding. |
| `marsh_scarp_closing_radius_metres` | `float` | `0.0` | `param` | - | Binary closing radius for the marsh scarp mask, in metres. | Used only when marsh_apply_scarp_closing is true. |
| `marsh_label_use_scarp_barriers` | `bool` | `false` | `param` | - | Treat scarps as barriers when labelling marsh-platform components. | Scarp pixels are removed from the raster used for final platform component labelling, so scarps separate adjacent platform surfaces. |
| `marsh_label_use_d4_connectivity` | `bool` | `false` | `param` | - | Use d4 connectivity when labelling marsh-platform components. | Useful when 1-pixel scarps should block diagonal leakage between neighbouring platform surfaces. |
| `marsh_min_platform_component_pixels` | `int` | `1` | `param` | - | Minimum component size for the final marsh-platform mask. | Set &gt; 1 to remove isolated tiny platform fragments. |
| `tidal_use_flood_mask` | `bool` | `false` | `param` | - | Use simple elevation-threshold channel mask instead of multi-scale median. | When true, any pixel whose elevation is at or below tidal_flood_elevation_metres is classified as a channel. All cleanup, UPL, and skeleton steps still apply. |
| `tidal_flood_elevation_metres` | `float` | `0.0` | `param` | - | Elevation threshold for the flood-mask channel method (metres). | Used only when tidal_use_flood_mask is true. Pixels at or below this value are marked as channel. |
| `tidal_median_window_radii_metres` | `string` | `5,10,20` | `param` | - | Comma-delimited median-window radii for tidal-channel extraction. | One residual raster is generated for each scale using median minus DEM. |
| `tidal_residual_thresholds_metres` | `string` | `0.1` | `param` | - | Comma-delimited residual thresholds for tidal-channel extraction. | If one value is supplied it is used for all scales. Otherwise there must be one threshold per scale. |
| `tidal_scale_mode` | `string` | `any_scale_pass` | `param` | - | How per-scale tidal-channel masks are combined. | Allowed values: any_scale_pass, all_scales_pass. |
| `tidal_apply_binary_closing` | `bool` | `false` | `param` | - | Apply binary closing to the tidal-channel mask. | Useful for sealing narrow gaps between neighbouring channel pixels. |
| `tidal_binary_closing_radius_metres` | `float` | `0.0` | `param` | - | Closing radius for the tidal-channel mask. | Used only when tidal_apply_binary_closing is true. |
| `tidal_fill_holes` | `bool` | `false` | `param` | - | Fill holes in the tidal-channel mask. | Useful when the thresholded channel mask contains enclosed voids. |
| `tidal_min_channel_component_pixels` | `int` | `1` | `param` | - | Minimum connected-component size for the channel mask. | Set &gt; 1 to remove small isolated channel fragments. |
| `tidal_extract_skeleton` | `bool` | `true` | `param` | - | Extract a 1-pixel skeleton from the tidal-channel mask. | Uses lsdtt3 raster skeletonisation tools. |
| `tidal_skeleton_spur_max_pixels` | `int` | `10` | `param` | - | Maximum spur length removed from the tidal skeleton. | Used only when tidal_extract_skeleton is true. |
| `tidal_skeleton_terminal_loop_max_pixels` | `int` | `10` | `param` | - | Maximum small terminal loop length removed from the tidal skeleton. | Used only when tidal_extract_skeleton is true. |
| `tidal_skeleton_min_component_pixels` | `int` | `1` | `param` | - | Minimum connected-component size for the tidal skeleton. | Set &gt; 1 to remove isolated skeleton fragments. |
| `tidal_connect_skeleton_segments` | `bool` | `false` | `param` | - | Attempt to reconnect nearby tidal skeleton endpoints. | Uses the existing lsdtt3 skeleton graph bridge-routing tools. |
| `tidal_connect_tangent_fit_length_pixels` | `float` | `6.0` | `param` | - | Tangent support length for tidal skeleton endpoint reconnection. | Used only when tidal_connect_skeleton_segments is true. |
| `tidal_connect_tangent_min_support_pixels` | `int` | `2` | `param` | - | Minimum support pixels for endpoint tangent estimation. | Used only when tidal_connect_skeleton_segments is true. |
| `tidal_connect_max_gap_pixels` | `float` | `20.0` | `param` | - | Maximum gap distance considered when reconnecting skeleton endpoints. | Used only when tidal_connect_skeleton_segments is true. |
| `tidal_connect_max_facing_angle_degrees` | `float` | `35.0` | `param` | - | Maximum endpoint-facing angle for skeleton reconnection candidates. | Used only when tidal_connect_skeleton_segments is true. |
| `tidal_connect_max_heading_misalignment_degrees` | `float` | `45.0` | `param` | - | Maximum tangent-heading mismatch for skeleton reconnection candidates. | Used only when tidal_connect_skeleton_segments is true. |
| `tidal_connect_max_candidate_score` | `float` | `1.0` | `param` | - | Maximum accepted candidate score for tidal skeleton reconnection. | Used only when tidal_connect_skeleton_segments is true. |
| `tidal_connect_route_corridor_half_width_pixels` | `float` | `12.0` | `param` | - | Half-width of the least-cost routing corridor for accepted skeleton bridges. | Used only when tidal_connect_skeleton_segments is true. |
| `tidal_compute_upl` | `bool` | `true` | `param` | - | Compute unchannelled path length (UPL). | UPL is the Euclidean distance to the nearest channel pixel. |
| `tidal_outlet_mode` | `string` | `edge_nodes` | `param` | - | How outlet nodes are chosen for tidal skeleton metrics. | Allowed values: edge_nodes, endpoint_nodes, outlet_mask, lowest_endpoint. |
| `tidal_station_spacing_metres` | `float` | `5.0` | `param` | - | Station spacing along tidal skeleton sections. | Stations are sampled at this spacing after section orientation. |
| `tidal_compute_cross_sections` | `bool` | `true` | `param` | - | Compute cross-sections at tidal skeleton stations. | Width, depth, and area metrics are estimated at each station. |
| `tidal_tangent_support_points` | `int` | `2` | `param` | - | Number of neighbouring stations used to estimate local section tangent. | Used in cross-section construction. |
| `tidal_cross_section_half_length_metres` | `float` | `25.0` | `param` | - | Half-length of tidal cross-sections. | Each cross-section spans twice this half-length. |
| `tidal_cross_section_sample_spacing_metres` | `float` | `1.0` | `param` | - | Sampling spacing along tidal cross-sections. | Controls width/depth/area estimation. |
| `tidal_max_width_ratio` | `float` | `3.0` | `param` | - | Maximum accepted section-length to local-width ratio. | Sections exceeding this are rejected as implausible. |
| `tidal_compute_station_strip_metrics` | `bool` | `true` | `param` | - | Assign each DEM pixel to its nearest tidal skeleton station and compute strip metrics. | This produces local and upstream strip areas and optional strip-mean raster values. |
| `tidal_strip_include_dem` | `bool` | `true` | `param` | - | Include DEM elevation in tidal station strip summaries. | If true, local and upstream strip mean elevation are computed. |
| `tidal_strip_include_upl` | `bool` | `false` | `param` | - | Include UPL in tidal station strip summaries. | Only available when the tidal-channel stage computed UPL. |
| `tidal_strip_include_marsh_mask` | `bool` | `false` | `param` | - | Include marsh-platform mask values in tidal station strip summaries. | Only available when the marsh-platform stage was run. |
| `marsh_points_write_vector_format` | `string` | `gpkg` | `param` | - | Output format for point vector outputs. | Allowed values: csv, flatgeobuf, fgb, gpkg. |
| `write_marsh_search_space_raster` | `bool` | `true` | `switch` | - | Write marsh search-space raster. | Binary raster defining the candidate marsh search domain. |
| `write_marsh_search_crossover_raster` | `bool` | `false` | `switch` | - | Write marsh search crossover raster. | Continuous DEM-relative-relief times slope-relative-relief raster. |
| `write_marsh_scarp_seed_mask_raster` | `bool` | `false` | `switch` | - | Write marsh scarp seed mask. | Binary raster of initial slope maxima used to seed scarps. |
| `write_marsh_scarp_order_raster` | `bool` | `false` | `switch` | - | Write marsh scarp order raster. | Stores iterative scarp-growth order values. |
| `write_marsh_scarp_mask_raster` | `bool` | `true` | `switch` | - | Write cleaned marsh scarp mask. | Binary raster of retained marsh scarps. |
| `write_marsh_platform_seed_mask_raster` | `bool` | `false` | `switch` | - | Write marsh platform seed mask. | Binary raster of scarp-seeded initial platform pixels. |
| `write_marsh_platform_provisional` | `bool` | `false` | `switch` | - | Write provisional marsh-platform mask before final cleanup. | Useful for diagnosing the cutoff-elevation stage. |
| `write_marsh_platform_mask_raster` | `bool` | `true` | `switch` | - | Write final marsh-platform mask. | Binary raster of the extracted marsh platform. |
| `write_marsh_platform_labels_raster` | `bool` | `true` | `switch` | - | Write labelled marsh-platform raster. | Connected marsh-platform components are labelled 1,2,3,.... |
| `write_marsh_polygons` | `bool` | `false` | `switch` | - | Vectorise labelled marsh-platform components to polygons. | Uses the existing raster polygonisation tools. |
| `write_tidal_median_rasters` | `bool` | `false` | `switch` | - | Write per-scale tidal moving-window median rasters. | One raster is written for each requested tidal scale. |
| `write_tidal_residual_rasters` | `bool` | `false` | `switch` | - | Write per-scale tidal residual rasters. | Residual = local median minus DEM, so positive values indicate local lows. |
| `write_tidal_threshold_masks` | `bool` | `false` | `switch` | - | Write per-scale thresholded tidal masks. | These are the binary residual-threshold masks prior to multi-scale combination. |
| `write_tidal_channel_mask_raster` | `bool` | `true` | `switch` | - | Write final tidal-channel mask. | Binary raster of detected tidal channels. |
| `write_tidal_upl_raster` | `bool` | `true` | `switch` | - | Write UPL raster. | UPL is the Euclidean distance to the nearest channel pixel. |
| `write_tidal_channel_skeleton` | `bool` | `true` | `switch` | - | Write raw tidal channel skeleton. | Written only when tidal_extract_skeleton is true. |
| `write_tidal_connected_skeleton` | `bool` | `true` | `switch` | - | Write connected tidal channel skeleton. | If no skeleton reconnection is performed this may be empty or match the raw skeleton. |
| `write_tidal_channel_polygons` | `bool` | `false` | `switch` | - | Vectorise tidal channel components to polygons. | Uses binary component labelling plus polygon tracing. |
| `polygonise_min_component_pixels` | `int` | `1` | `param` | - | Minimum component size used when polygonising binary masks. | Set &gt; 1 to remove tiny components before vectorisation. |
| `write_tidal_node_id_raster` | `bool` | `true` | `switch` | - | Write tidal skeleton node-id raster. | Node pixels store node ids; other pixels are nodata. |
| `write_tidal_section_id_raster` | `bool` | `true` | `switch` | - | Write tidal skeleton section-id raster. | Skeleton pixels store section ids. |
| `write_tidal_section_strahler_raster` | `bool` | `true` | `switch` | - | Write tidal skeleton section Strahler-order raster. | Skeleton pixels store section Strahler order. |
| `write_tidal_section_hack_raster` | `bool` | `true` | `switch` | - | Write tidal skeleton section Hack-order raster. | Skeleton pixels store section Hack order. |
| `write_tidal_nearest_station_id_raster` | `bool` | `true` | `switch` | - | Write nearest-station id raster. | Each valid DEM pixel stores the id of its nearest tidal station. |
| `write_tidal_station_downstream_distance_raster` | `bool` | `true` | `switch` | - | Write station downstream-distance raster. | Each valid DEM pixel stores the downstream distance of its nearest station. |
| `write_tidal_sections_vector` | `bool` | `true` | `switch` | - | Write tidal skeleton sections as lines. | Each oriented tidal section becomes one line feature. |
| `write_tidal_stations_vector` | `bool` | `true` | `switch` | - | Write tidal skeleton stations as points. | Each station becomes one point feature with station-scale metrics. |
| `write_tidal_stations_csv` | `bool` | `true` | `switch` | - | Write tidal skeleton stations to CSV. | This is a table export of station-scale metrics. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). CSV is available only through selectors for point/table outputs. Per-output overrides exist only where a container genuinely differs. Raster outputs are governed separately by write_output_format. |

<!-- END GENERATED CLI PARAMETERS -->
