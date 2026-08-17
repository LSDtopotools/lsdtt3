# CLI reference: lsdtt-linear-stats

Extracts statistics from linear features (channels, ridges, faults) in raster or vector form. Computes bearing analysis, sinuosity, length distributions, and optional rose diagrams.

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

Keys use `snake_case`. All keys and their defaults are listed below.

---

## Pipeline overview

Steps run in this order:

| Step | Description |
|---|---|
| 1. Input parsing | Load raster or vector linear features |
| 2. Skeletonisation | Extract 1-pixel-wide centrelines (if raster input) |
| 3. Bearing calculation | Compute local bearing angle along features |
| 4. Statistics | Compute sinuosity, length, orientation distribution |
| 5. Outputs | Write CSV statistics and optional rose diagram |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing input data. |
| `read_fname` | string | parameter file stem | Input feature filename (raster or vector). |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |

---

## Input source

| Parameter | Type | Default | Description |
|---|---|---|---|
| `input_type` | string | `raster` | Input type: `raster` (binary mask) or `vector` (polyline layer). |
| `vector_format` | string | `gpkg` | Vector input format (when `input_type: vector`): `gpkg`, `flatgeobuf`, `shapefile`. |

---

## Raster processing

Used only when `input_type: raster`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `skeletonise_raster` | bool | `true` | Extract 1-pixel centreline via skeletonisation. |
| `skeleton_prune_length_pixels` | int | `3` | Prune skeleton spurs shorter than this (pixels). |

---

## Bearing analysis

| Parameter | Type | Default | Description |
|---|---|---|---|
| `bearing_window_pixels` | int | `11` | Window size (pixels) for local bearing estimation. Odd integers only. |
| `bearing_bin_width_degrees` | float | `10.0` | Bin width (degrees) for bearing histogram. |

---

## Output: Statistics CSV

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_segment_statistics` | bool | `true` | Write per-segment statistics (length, sinuosity, mean bearing). |
| `write_bearing_histogram` | bool | `true` | Write bearing histogram CSV. |
| `write_sinuosity_statistics` | bool | `true` | Write sinuosity distribution statistics. |

Output files:
- `<write_fname>_segment_stats.csv`
- `<write_fname>_bearing_histogram.csv`
- `<write_fname>_sinuosity_stats.csv`

---

## Output: Rose diagram

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_rose_diagram_raster` | bool | `false` | Write bearing distribution as a polar raster (rose diagram). |
| `rose_diagram_size_pixels` | int | `512` | Size of rose diagram raster (pixels × pixels). |

Output file: `<write_fname>_rose_diagram.tif`

---

## Output format

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `cog` | Raster output format (for rose diagram): `cog`, `gtiff`, `same_as_input`. |
