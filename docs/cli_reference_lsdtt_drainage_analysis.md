# CLI reference: lsdtt-drainage-analysis

Computes drainage-related rasters from a DEM. Supports D8, D-infinity (Tarboton 1997), and MFD (Quinn et al. 1991 / Freeman 1991) drainage area; Topographic Wetness Index (TWI); sediment connectivity index (Cavalli et al. 2013); and height above channel.

At least one `write_*_raster` flag must be set to `true` or the CLI exits without writing anything.

The DEM must be in a projected CRS (e.g. UTM).

## Usage

```
lsdtt-drainage-analysis <path> <parameter_file>
lsdtt-drainage-analysis <parameter_file>         # path defaults to ./
lsdtt-drainage-analysis -h                       # print help and exit
lsdtt-drainage-analysis -v                       # print version and exit
lsdtt-drainage-analysis -c                       # print citation and exit
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
| `write_output_format` | string | `cog` | Output format for all raster outputs. Options: `cog`, `gtiff`, `same_as_input`. |

---

## Threading and routing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` uses all available hardware threads. Applies to pit filling, flow routing, and polyfit. |
| `boundary_conditions` | string | `oooo` | Boundary conditions for flow routing. Four-character string (N,E,S,W): `o` = open, `n` = no-flux, `p` = periodic, `b` = forced outlet. |

---

## DEM preparation

| Parameter | Type | Default | Description |
|---|---|---|---|
| `raster_is_filled` | bool | `false` | If true, skip pit filling and use the input raster directly. |
| `fill_min_slope` | float | `0.0001` | Minimum gradient imposed during pit filling. Used only when `raster_is_filled: false`. |

---

## Drainage area outputs

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_d8_drainage_area_raster` | bool | `false` | `<write_fname>_d8_drainage_area.tif` — D8 single-receiver drainage area (m²). |
| `write_dinf_drainage_area_raster` | bool | `false` | `<write_fname>_dinf_drainage_area.tif` — D-infinity drainage area (m²). |
| `write_mfd_drainage_area_raster` | bool | `false` | `<write_fname>_mfd_drainage_area.tif` — MFD drainage area (m²). |
| `mfd_exponent` | float | `1.0` | Slope-weighting exponent for MFD partitioning. Quinn et al. (1991) use 1.0; Freeman (1991) suggests 1.1. Used when MFD drainage area or MFD-based TWI is requested. |

---

## Topographic Wetness Index

TWI = ln( (A / pixel_area) / slope ), where A is upstream drainage area (m²), pixel_area = cell_size², and slope is from a local quadratic polyfit. Only valid where area > 0 and slope > 0.

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_twi_raster` | bool | `false` | `<write_fname>_twi.tif` |
| `twi_drainage_area_method` | string | `Dinf` | Flow routing method for TWI drainage area. Options: `Dinf`, `MFD`. |
| `twi_window_radius` | float | `100.0` | Polyfit window radius (m) for slope computation used in TWI. |

---

## Connectivity index

Sediment connectivity index IC = log₁₀(D_up / D_down) following Cavalli et al. (2013). D8 routing is always used for path calculations. Writing connectivity also writes `_connectivity_d_up.tif` and `_connectivity_d_down.tif` for diagnostics.

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_connectivity_raster` | bool | `false` | `<write_fname>_connectivity_index.tif` |
| `connectivity_drainage_area_method` | string | `d8` | Flow routing method for drainage area in D_up. Options: `d8`, `dinf`, `mfd`. |
| `connectivity_area_threshold_m2` | float | `1e6` | Drainage area threshold (m²) defining channel / sink pixels. D_down accumulation stops when the flow path reaches a channel. |
| `connectivity_window_radius` | float | `100.0` | Moving-window radius (m) for slope and roughness computation. Square window always used. |
| `connectivity_slope_method` | string | `polyfit_direct` | Method for computing slope in IC. Options: `polyfit_direct` (fast, square window), `wavelet` (FFT-based, requires FFTW). |
| `connectivity_roughness_method` | string | `madk2` | Roughness metric for the W weighting factor. Options: `madk2` (Median Absolute Deviation of k2 differences, Trevisani & Guth 2025), `iqr_residual` (IQR of polyfit residuals). W = max(1 − R/R_max, 0.001). |

---

## Height above channel

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_height_above_channel_raster` | bool | `false` | `<write_fname>_height_above_channel.tif` — elevation difference from each pixel to the nearest downstream channel pixel along the D8 flow path. |
| `channel_threshold_contributing_pixels` | int | `1000` | Minimum D8 contributing-pixel count to classify a pixel as a channel. Tune to DEM resolution: at 5 m, 1000 pixels ≈ 0.025 km²; at 30 m ≈ 0.9 km². |

---

## Example parameter file: drainage area and TWI

```
read_path:  /data/dems/
read_fname: my_dem
write_path: /data/output/
write_fname: my_dem_drainage

raster_is_filled:    false
fill_min_slope:      0.0001
boundary_conditions: oooo
n_threads:           0
write_output_format: cog

write_d8_drainage_area_raster:   true
write_dinf_drainage_area_raster: true

write_twi_raster:               true
twi_drainage_area_method:       Dinf
twi_window_radius:              100.0
```

## Example parameter file: connectivity index

```
read_path:  /data/dems/
read_fname: my_dem
write_path: /data/output/
write_fname: my_dem_connectivity

raster_is_filled:    true

write_connectivity_raster:             true
connectivity_drainage_area_method:     dinf
connectivity_area_threshold_m2:        500000.0
connectivity_window_radius:            50.0
connectivity_slope_method:             polyfit_direct
connectivity_roughness_method:         madk2
```
