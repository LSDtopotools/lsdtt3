# CLI reference: lsdtt-basin-stats

Computes hydromorphic statistics for basins defined by user-supplied outlet points. Requires a CSV file with outlet coordinates in WGS84; the program builds D8 flow routing, delineates each basin, and writes a summary CSV with elevation, slope, flow-distance, and optionally channel-length statistics.

## Usage

```
lsdtt-basin-stats  <path> <parameter_file>
lsdtt-basin-stats  <parameter_file>        # path defaults to ./
lsdtt-basin-stats  -h                      # print help and exit
lsdtt-basin-stats  -v                      # print version and exit
lsdtt-basin-stats  -c                      # print citation and exit
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
| 1. Sea masking | Pixels at or below `minimum_elevation` are set to nodata |
| 2. Pit filling | Optional (`fill_before_routing`) |
| 3. D8 flow routing | Build the D8 flow-info object |
| 4. Outlet snapping | Optional (`basin_outlet_snapping_method`) |
| 5. Basin delineation | Extract upslope area for each outlet |
| 6. Statistics | Compute elevation, slope, flow-distance stats per basin |
| 7. Channel stats | Optional Strahler-order channel-length calculation |
| 8. Write outputs | `_basin_stats.csv` and optional per-basin PDF/CDF CSVs |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input DEM filename. Any GDAL-supported format. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `write_output_format` | string | `cog` | Raster output format. `cog`, `gtiff`, or `same_as_input`. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` = all available hardware threads. |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `minimum_elevation` | float | `0.0` | Pixels at or below this elevation (m) are masked to nodata. The default removes sea and ocean pixels. |
| `fill_before_routing` | bool | `true` | Fill DEM pits before building flow routing. Recommended. |
| `fill_min_slope` | float | `0.0001` | Minimum slope (m/m) imposed during pit filling. |

---

## Basin outlet input

The outlet CSV **must** contain columns named `latitude` and `longitude` (decimal degrees, EPSG:4326). Optional columns `gauge_id` and `gauge_name` are passed through to the output.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `basin_outlet_csv` | string | `NULL` | **Required.** Path to the outlet CSV file. |

---

## Outlet snapping

| Parameter | Type | Default | Description |
|---|---|---|---|
| `basin_outlet_snapping_method` | string | `none` | Snap outlets to pixels exceeding a contributing-pixel threshold. `none`: use the mapped raster pixel directly. `expanding_square`: search ring-by-ring from the input pixel. `downstream`: follow the D8 receiver chain until a qualifying pixel is found. |
| `basin_outlet_contributing_pixels_threshold` | int | `1000` | Minimum contributing-pixel count for snapping. |
| `basin_outlet_min_search_radius_pixels` | int | `0` | Minimum search radius for `expanding_square` (pixels). `0` tests the input pixel first. |
| `basin_outlet_max_search_radius_pixels` | int | `25` | Maximum search radius for `expanding_square` (pixels). |

---

## Basin quality control

| Parameter | Type | Default | Description |
|---|---|---|---|
| `reject_nodata_edge_basins` | bool | `true` | Skip basins whose upslope area touches any nodata or raster-edge pixel. Prevents incorrect statistics for incomplete basins. |

---

## Hillshade (optional QC output)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_hillshade_raster` | bool | `false` | Write a hillshade raster after sea masking. |
| `hillshade_azimuth` | float | `315` | Sun azimuth (degrees clockwise from north). |
| `hillshade_altitude` | float | `45` | Sun altitude above horizon (degrees). |
| `hillshade_z_factor` | float | `1.0` | Vertical exaggeration factor. |

---

## Channel statistics

When enabled, builds a Strahler-order channel network and computes total channel length for orders 1–6.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `calculate_channel_stats` | bool | `true` | Build a channel/junction network and compute channel length by Strahler order. |
| `threshold_contributing_pixels` | int | `500` | Contributing-pixel threshold defining channel sources for channel statistics. |

---

## Per-basin elevation distribution

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_basin_pdf_cdf_csv` | bool | `false` | Write one elevation PDF/CDF CSV per basin. |
| `n_histogram_bins` | int | `100` | Number of histogram bins for the elevation PDF/CDF. |

Per-basin output: `<write_fname>_basin_<id>_pdf_cdf.csv`

Columns: `basin_id`, `bin_index`, `bin_lower_edge`, `bin_upper_edge`, `bin_centre`, `bin_width`, `count`, `probability`, `pdf`, `cdf`.

---

## Output

Main output: `<write_fname>_basin_stats.csv`

| Column | Description |
|---|---|
| `basin_id` | Sequential integer assigned by the program |
| `gauge_id`, `gauge_name` | Pass-through from the input CSV if present |
| `latitude`, `longitude` | Requested outlet coordinates (EPSG:4326) |
| `outlet_x`, `outlet_y` | Outlet in projected CRS |
| `outlet_row`, `outlet_col`, `outlet_flat` | Raster cell indices |
| `n_pixels` | Number of pixels in the basin |
| `drainage_area_km2` | Basin drainage area (km²) |
| `nodata_on_edge` | 1 if any upslope pixel touches nodata or raster edge |
| `centroid_x`, `centroid_y` | Area-weighted centroid in projected CRS |
| `centroid_latitude`, `centroid_longitude` | Centroid in WGS84 |
| `elevation_min`, `elevation_max`, `elevation_range` | Elevation statistics (m) |
| `elevation_mean`, `elevation_25th`, `elevation_50th`, `elevation_75th` | |
| `elevation_iqr`, `elevation_stddev`, `elevation_mad` | |
| `slope_25th`, `slope_50th`, `slope_75th` | Gradient-magnitude percentiles (m/m) |
| `flow_distance_25th`, `flow_distance_50th`, `flow_distance_75th` | Flow distance from outlet (m) |
| `channel_length_O1` … `channel_length_O6` | Channel length (m) by Strahler order (when `calculate_channel_stats: true`) |

---

## Minimal example parameter file

```
# lsdtt-basin-stats parameter file
read_path:         /data/dems/
read_fname:        catchment_dem.tif

basin_outlet_csv:                       /data/outlets.csv
basin_outlet_snapping_method:           downstream
basin_outlet_contributing_pixels_threshold: 2000

reject_nodata_edge_basins:  true
calculate_channel_stats:    true
threshold_contributing_pixels: 500
```
