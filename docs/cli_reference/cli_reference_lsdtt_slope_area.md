# CLI reference: lsdtt-slope-area

`lsdtt-slope-area` calculates slope–area statistics from a DEM for one or more
direct-polynomial-fit smoothing lengths. It can route flow with D-infinity or
D8, optionally clip the analysis with a same-grid basin mask, write per-lambda
slope rasters, and write quantile statistics in logarithmic drainage-area
bins.

## Usage

```text
lsdtt-slope-area <path> <parameter_file>
lsdtt-slope-area <parameter_file>
lsdtt-slope-area -h
lsdtt-slope-area -v
lsdtt-slope-area -c
```

## Parameter file format

```text
# Lines beginning with # are comments.
key: value
```

## Pipeline overview

1. Parse the parameter file and load the DEM.
2. Optionally fill pits using the configured minimum slope and boundaries.
3. Build D-infinity or D8 flow information and drainage area.
4. Optionally write flow-information and drainage-area products.
5. For each requested polynomial-fit length, calculate slope.
6. Apply the optional basin mask and logarithmic drainage-area bins.
7. Write per-lambda statistics, slope rasters, and summary results.

## Input/output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | parameter-file directory | Directory containing the DEM. |
| `read_fname` | string | — | Input DEM filename. |
| `write_path` | string | `read_path` | Output directory. |
| `write_fname` | string | DEM stem | Output filename prefix. |
| `write_output_format` | string | `cog` | Raster output format where applicable. |
| `write_flowinfo` | bool | `true` | Write flow-information products. |
| `write_drainage_area` | bool | `true` | Write the drainage-area raster. |
| `write_slope_rasters` | bool | `true` | Write one slope raster per lambda. |
| `summary_csv_fname` | string | automatic | Summary CSV filename. |

## Routing and pit filling parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fill_before_routing` | bool | `true` | Fill pits before routing. |
| `fill_min_slope` | float | `0.0001` | Minimum slope used during pit filling. |
| `flow_algorithm` | string | `dinf` | `dinf` or `d8`. |
| `boundary_conditions` | string | `oooo` | Four-character boundary code, in `N,E,S,W` order. |

## Polynomial-fit and binning parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `polyfit_lambda_metres_list` | list | — | Comma-separated smoothing lengths in metres. |
| `n_log_bins` | int | — | Number of logarithmic drainage-area bins. |

## Basin clipping parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `basin_clip_mode` | string | `none` | `none` or `mask`. |
| `basin_mask_fname` | string | — | Same-grid raster mask when mode is `mask`. |

The basin mask must have exactly the same rows, columns, cell size, origin, and
coordinate reference system as the DEM. Nonzero values are inside the basin.
This CLI does not rasterize vector polygons; use `lsdtt-basin-extraction` first.

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of worker threads; `0` uses the implementation default. |

## Outputs

For each lambda, the tool writes:

```text
<write_fname>_lambda_XXm.csv
<write_fname>_lambda_XXm_slope.tif
```

The per-lambda CSV columns are:

```text
drainage_area,d05_slope,d25_slope,d50_slope,d75_slope,d95_slope,number_of_pixels
```

The summary CSV columns are:

```text
polyfit_lambda,drainage_area_with_max_d50_slope
```

Decimal lambda values use underscores in filenames, for example
`lambda_12_5m`.

When `write_slope_rasters: true`, an existing expected slope raster is opened
and reused instead of recalculating polynomial-fit slope. The raster is
validated against the DEM grid, after which CSV statistics are regenerated with
the current bin count and basin mask.

## Example parameter file

```text
read_path: /data/Moreno
read_fname: Moreno.tif
write_path: /results/Moreno/slope_area
write_fname: Moreno_slope_area
write_output_format: cog
n_threads: 4
boundary_conditions: opop
fill_before_routing: true
fill_min_slope: 0.0001
flow_algorithm: dinf
write_flowinfo: true
write_drainage_area: true
polyfit_lambda_metres_list: 3,15,21
n_log_bins: 200
write_slope_rasters: true
basin_clip_mode: none
```
