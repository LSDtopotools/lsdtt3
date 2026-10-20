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
| `write_flowinfo_rasters` | bool | `true` | Write the drainage-area (and, for d8, flow-direction) rasters. |
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
write_flowinfo_rasters: true
polyfit_lambda_metres_list: 3,15,21
n_log_bins: 200
write_slope_rasters: true
basin_clip_mode: none
```

<!-- BEGIN GENERATED CLI PARAMETERS -->
## Parameters

This section is generated from the executable's `-h json` registry. Do not edit it manually.

**Status: stable.**

### File I/O

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `read_path` | `string` | `./` | `io` | - | Directory from which input files are read. | Defaults to the directory holding the parameter file. |
| `read_fname` | `string` | - | `io` | - | Input filename, with or without extension (e.g. my_dem.tif). | Read from read_path. lsdtt3 uses GDAL, so any GDAL-readable raster format works. If empty, defaults to the parameter filename stem. |
| `write_path` | `string` | - | `io` | - | Directory to which output files are written. | Defaults to read_path when not set. |
| `write_fname` | `string` | - | `io` | - | Output filename prefix (no path, no extension). | All outputs are named &lt;write_fname&gt;_&lt;descriptor&gt;.&lt;ext&gt;. Defaults to the stem of read_fname when not set. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_before_routing` | `bool` | `true` | `switch` | - | Choose whether to pit-fill the DEM before building the flow network. | ASKS: should the DEM be filled at all? Filling is OPTIONAL here and this is a genuine analytical choice. True (the default) routes over a pit-free surface. False builds the flow network on the UNFILLED surface, so flow paths terminate at sinks instead of being carried through them. Do not confuse this with raster_is_filled, which is an assertion that the input is already filled. |

### DEM preparation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_min_slope` | `float` | `0.0001` | `param` | `fill_before_routing` | Minimum slope imposed during pit filling. | Used when fill_before_routing is true. |

### Flow information

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `flow_algorithm` | `string` | `dinf` | `param` | - | Flow-routing algorithm. | Allowed values are dinf or d8. |
| `write_flowinfo_rasters` | `bool` | `true` | `switch` | - | Write the drainage-area and flow-direction rasters. | The drainage-area raster is used by the slope-area analysis. |
| `write_log10_drainage_area_raster` | `bool` | `false` | `switch` | - | Write the base-10 logarithm of the drainage-area raster. | Values are log10 of each positive drainage-area pixel; non-positive and nodata pixels remain nodata. |
| `boundary_conditions` | `string` | `oooo` | `param` | - | Boundary conditions for flow routing. | Four characters in north, east, south, west order. |

### Threading

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads. | 0 uses the runtime default. |

### Slope-area analysis

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `polyfit_lambda_metres_list` | `string` | `15` | `param` | - | Comma-separated direct-polyfit lambdas in metres. | Example: 3,5,7,15,30. Decimal filename labels use _ instead of . |
| `n_log_bins` | `int` | `20` | `param` | - | Number of logarithmic drainage-area bins. | Only non-empty bins are written to each lambda CSV. |
| `write_slope_rasters` | `bool` | `true` | `switch` | - | Write each slope raster, after optional clipping. | Files are named &lt;write_fname&gt;_lambda_XXm_slope.tif. |

### Basin clipping

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `basin_clip_mode` | `string` | `none` | `param` | - | Optional basin clipping mode. | Allowed values are none, mask, or polygon. |
| `basin_mask_fname` | `string` | - | `param` | `basin_clip_mode=mask` | Raster basin mask filename. | The mask must have exactly the DEM grid settings. Cells with value 0 or nodata are excluded; any other value is kept, so a labelled basins raster works. A relative path is read from read_path. |
| `basin_polygon_fname` | `string` | - | `param` | `basin_clip_mode=polygon` | Basin polygon filename. | The polygon must be supplied as a same-grid rasterised polygon mask; vector rasterisation is not performed by this CLI. It is read exactly like basin_mask_fname. A relative path is read from read_path. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `summary_csv_fname` | `string` | - | `param` | - | Summary CSV filename. | Empty uses &lt;write_fname&gt;_maximum_d50_slope_area.csv. |

<!-- END GENERATED CLI PARAMETERS -->
