# CLI reference: lsdtt-swath-profile

A swath profile along an ordered baseline polyline. Every valid pixel within
`half_width_metres` of the baseline is projected onto it and binned by
distance along the baseline. Pixels beyond either end of the baseline (behind
its start or past its end) are not alongside it and are left out. Each bin reports the pixel count, the mean and a
21-level percentile distribution (p0, p5, …, p100) of a value raster, which
is the DEM's own elevation unless you give `value_raster_fname`. The tool also
writes rasters of the along-axis and perpendicular distance of every pixel.

## Usage

```
lsdtt-swath-profile  <path> <parameter_file>
lsdtt-swath-profile  <parameter_file>        # path defaults to ./
lsdtt-swath-profile  -h                      # print help and exit
lsdtt-swath-profile  -v                      # print version and exit
lsdtt-swath-profile  -c                      # print citation and exit
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
    `baseline_vector_fname` accepts multipart geometry (`MultiLineString` and friends). Each feature is split into its parts on read, with the attributes copied onto each, and the split is reported in the log.

    Shapefiles are the commonest source of multipart data, and **versions before August 2026 discarded it silently** — see [Vector data](../what_is_lsdtt3/vector_data.md).

## Example

```
read_fname: my_dem.tif
baseline_vector_fname: valley_axis.fgb
half_width_metres: 2000
bin_width_metres: 250
```

The first LineString feature in `baseline_vector_fname` is the baseline, in
the order its vertices are stored. `half_width_metres: 0` (the default)
includes every valid pixel alongside the baseline.

!!! note "Pixels beyond the ends of the baseline"
    A pixel is in the swath only if its nearest point on the baseline is not
    one of the baseline's ends with the pixel lying beyond it. Versions up to
    v0.9.5 included the pixels within `half_width_metres` of each end point,
    in a half-disc behind the start and past the end, and counted them in the
    first and last bins: a 30 m wide cross section in 2 m bins had about 350
    extra pixels in each end bin against 60 in the others. Those versions
    could also add an empty bin past the end when the baseline's length was a
    whole number of bins plus a rounding error.

A pixel beyond the ends has nodata in the distance rasters, like a pixel
farther than `half_width_metres` from the baseline.

## Outputs

| File | Written when | Contents |
|---|---|---|
| `<write_fname>_swath_profile.csv` | always | One row per bin: `along_axis_centre` (m along the baseline), `n_pixels`, `mean`, then `p0`, `p5`, …, `p100`. |
| `<write_fname>_swath_along_axis.tif` | `write_distance_rasters: true` (default) | Distance along the baseline of each pixel's projection (m), for the pixels in the swath. |
| `<write_fname>_swath_perp_distance.tif` | as above | Perpendicular distance from the baseline (m). |
| `<write_fname>_swath_signed_perp.tif` | as above | Perpendicular distance, signed by which side of the baseline the pixel is on (m). |

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

### Baseline / value input

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `baseline_vector_fname` | `string` | - | `param` | - | Vector file (FlatGeobuf/GeoPackage/Shapefile) with the swath baseline. | The first LineString feature is used as the ordered baseline polyline. Coordinates must be in the raster's projected CRS (metres). A relative path is read from read_path. |
| `value_raster_fname` | `string` | `NULL` | `param` | - | Optional raster whose values are binned along the swath. | If NULL (default), the input DEM's own values (elevation) are binned. Must share the input DEM's geometry. A relative path is read from read_path. |

### Swath geometry

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `half_width_metres` | `float` | `0.0` | `param` | - | Only pixels within this distance of the baseline are included. | Set to 0 to include all valid pixels alongside the baseline (full-extent swath). Pixels beyond either end of the baseline are never included. |
| `bin_width_metres` | `float` | `1000.0` | `param` | - | Along-axis bin width in metres. Must be &gt; 0. | - |
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads (0 = all available). | - |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_distance_rasters` | `bool` | `true` | `switch` | - | Write the along-axis, perpendicular, and signed-perpendicular rasters. | Output: &lt;prefix&gt;_swath_along_axis.tif, _swath_perp_distance.tif, _swath_signed_perp.tif. |
| `write_output_format` | `string` | `cog` | `io` | - | Raster output format for the distance rasters. | Allowed values: cog, gtiff, same_as_input. |

<!-- END GENERATED CLI PARAMETERS -->
