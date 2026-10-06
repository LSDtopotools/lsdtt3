# CLI reference: lsdtt-slope-area

Calculates slope-area statistics from a DEM for one or more direct-polyfit
smoothing lengths. The tool can route flow with D-infinity or D8, optionally
clip the analysis with a same-grid basin mask, write per-lambda slope rasters,
and write quantile statistics in logarithmic drainage-area bins.

## Usage

```text
lsdtt-slope-area <path> <parameter_file>
lsdtt-slope-area <parameter_file>       # path defaults to ./
lsdtt-slope-area -h                     # print help and exit
lsdtt-slope-area -v                     # print version and exit
lsdtt-slope-area -c                     # print citation and exit
```

## Parameter file format

```text
# Lines beginning with # are comments.
key: value
```

## Example

```text
read_path: /data/Moreno
read_fname: Moreno.tif
write_path: /results/Moreno/slope_area
write_fname: Moreno_slope_area
n_threads: 4
fill_before_routing: true
fill_min_slope: 0.0001
flow_algorithm: dinf
write_flowinfo: true
write_drainage_area: true
polyfit_lambda_metres_list: 3,15,21
n_log_bins: 200
write_slope_rasters: true
basin_clip_mode: mask
basin_mask_fname: /results/Moreno/slope_area/Moreno_basins.tif
summary_csv_fname: Moreno_slope_area_maximum_d50_slope.csv
```

## Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | parameter-file directory | Directory containing the DEM. |
| `read_fname` | string | — | Input DEM filename. |
| `write_path` | string | `read_path` | Output directory. |
| `write_fname` | string | DEM stem | Output filename prefix. |
| `n_threads` | int | `0` | Number of worker threads. |
| `fill_before_routing` | bool | `true` | Fill pits before flow routing. |
| `fill_min_slope` | float | `0.0001` | Minimum slope used during pit filling. |
| `flow_algorithm` | string | `dinf` | Flow routing algorithm: `dinf` or `d8`. |
| `write_flowinfo` | bool | `true` | Write flow-information products. |
| `write_drainage_area` | bool | `true` | Write the drainage-area raster. |
| `polyfit_lambda_metres_list` | list | — | Comma-separated smoothing lengths in metres. |
| `n_log_bins` | int | — | Number of logarithmic drainage-area bins. |
| `write_slope_rasters` | bool | `true` | Write one slope raster for each lambda. |
| `basin_clip_mode` | string | `none` | Use `none` or `mask`. |
| `basin_mask_fname` | string | — | Raster mask used when `basin_clip_mode: mask`. |
| `summary_csv_fname` | string | automatic | Maximum-d50 summary filename. |

## Basin masks

The basin mask must have exactly the same grid as the DEM: rows, columns, cell
size, origin, and coordinate reference system. Nonzero mask values are treated
as inside the basin. Vector polygon rasterization is not performed by this
CLI; use `lsdtt-basin-extraction` to create a raster basin mask first.

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

## Reusing slope rasters

When `write_slope_rasters: true`, an existing expected slope raster is opened
and reused instead of recalculating polyfit slope. The existing raster is
validated against the DEM grid. CSV statistics are then regenerated using the
current bin count and basin mask.
