# CLI reference: lsdtt-linear-stats

Samples a raster along each line feature in a vector file and writes
descriptive statistics of the sampled values for every line. Points are taken
every `interpolating_step_length` metres along each line; samples that fall
outside the raster are skipped, and a line wholly outside it gets no
statistics.

## Usage

```
lsdtt-linear-stats  <path> <parameter_file>
lsdtt-linear-stats  <parameter_file>        # path defaults to ./
lsdtt-linear-stats  -h                      # print help and exit
lsdtt-linear-stats  -v                      # print version and exit
lsdtt-linear-stats  -c                      # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

Keys use `snake_case`. All keys and their defaults are listed under
[Parameters](#parameters) below.

---

!!! note "Multipart lines"
    `vector_fname` accepts multipart geometry (`MultiLineString` and friends). Each feature is split into its parts on read, with the attributes copied onto each, and the split is reported in the log.

    Shapefiles are the commonest source of multipart data, and **versions before August 2026 discarded it silently** — see [Vector data](../what_is_lsdtt3/vector_data.md).

## Example

```
read_path: /home/alice/data/
read_fname: my_dem.tif
vector_fname: /home/alice/data/faults.fgb
interpolating_step_length: 5
```

`vector_fname` is opened as given, **not** relative to `read_path`, so give
the full path (or a path relative to the directory you run the tool from).
Only the first layer of the vector file is read, and the lines must be in the
raster's projected CRS; they are not reprojected. The raster must be in a
projected CRS with metre units.

## Output

`<write_fname>_linear_stats.fgb` (or `.gpkg` with `write_vector_format: gpkg`):
the input lines (lines wholly outside the raster are left out), each with
these attributes of the raster values sampled along it. The input file's own
attributes are not carried over.

`length` (m), `n_values`, `min`, `max`, `range`, `sum`, `mean`, `median`,
`variance_population`, `variance_sample`, `stddev_population`,
`stddev_sample`, `mean_absolute_deviation` (about the mean),
`median_absolute_deviation` (about the median), `mad_normal_scaled` (the
median absolute deviation × 1.4826), `q25`, `q75`, `iqr`.

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

### param

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `vector_fname` | `string` | - | - | - | File name for the dataset containing linear features. | Input |
| `interpolating_step_length` | `float` | `1.0` | - | - | Stride length for sample points when interpolating along linear features. | Process parameters |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vector_format` | `string` | `flatgeobuf` | `param` | - | Format for the output file. | Allowed values: gpkg, flatgeobuf/fgb. |

<!-- END GENERATED CLI PARAMETERS -->
