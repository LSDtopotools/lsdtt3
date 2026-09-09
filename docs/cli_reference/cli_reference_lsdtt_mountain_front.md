# CLI reference: lsdtt-mountain-front

Identifies mountain masses in a DEM via adaptive local-relief thresholding. The pipeline classifies each pixel as either mountain or flat/piedmont, then optionally extracts the points where drainage crosses the mountain front boundary.

## Usage

```
lsdtt-mountain-front <path> <parameter_file>
lsdtt-mountain-front <parameter_file>         # path defaults to ./
lsdtt-mountain-front -h                       # print help and exit
lsdtt-mountain-front -v                       # print version and exit
lsdtt-mountain-front -c                       # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Pipeline overview

The pipeline runs all stages unconditionally. Output writing is controlled by the `write_*` parameters.

| Stage | Description |
|---|---|
| 1 | Compute local relief (moving square window) |
| 2 | Build flat mask: adaptive CDF divergence or fixed threshold |
| 3 | Remove small flat components |
| 4 | Morphological closing of flat mask (seal narrow gaps) |
| 5 | Morphological opening of flat mask (sever thin peninsulas) |
| 6 | Invert flat mask → initial mountain mask |
| 7 | Fill holes in mountain mask (enclosed non-mountain pixels become mountain) |
| 8 | Remove small mountain components |
| 9 | Valley infilling via morphological closing (optional) |
| 10 | Identify mountain front outlet points (optional) |

The `_relief_cdf.csv` file is always written when `adaptive_relief_threshold` is true. It contains the 500-bin empirical CDF and the best-fit parametric CDF, and is the primary tool for choosing a manual threshold if the adaptive result is not satisfactory.

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory containing the input DEM. |
| `read_fname` | string | parameter file stem | Input raster filename. |
| `write_path` | string | same as `read_path` | Directory for output files. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |
| `write_output_format` | string | `cog` | Output raster format. `cog`, `gtiff`, or `same_as_input`. |

---

## Threading

| Parameter | Type | Default | Description |
|---|---|---|---|
| `n_threads` | int | `0` | Number of OpenMP threads. `0` uses all available hardware threads. |

---

## Stage 1: Local relief

| Parameter | Type | Default | Description |
|---|---|---|---|
| `relief_window_side_pixels` | int | `5` | Full window **side length** in pixels for local relief (max − min). Not a radius: must be odd and ≥ 3, so `5` is a 5×5 neighbourhood. `lsdtt-fluvial-surfaces` spells its equivalent `fallback_relief_window_radius_pixels`, which is a **radius**, so the same number means a different window there. The old name `relief_window_pixels` is accepted as a deprecated alias. |

---

## Stage 2: Flat mask threshold

| Parameter | Type | Default | Description |
|---|---|---|---|
| `adaptive_relief_threshold` | bool | `true` | If true, automatically select the relief threshold using CDF divergence. Computes a 500-bin empirical CDF, fits the best parametric distribution, and picks the x-value where the empirical CDF most exceeds the fitted CDF. Set to false to supply a fixed threshold instead. |
| `relief_threshold_metres` | float | `1.0` | Fixed local-relief threshold in metres below which a pixel is classed flat. Only used when `adaptive_relief_threshold` is false. Read a good value directly from the `_relief_cdf.csv` output by running the tool once with adaptive thresholding enabled. |

---

## Stages 3–8: Flat and mountain mask filtering

| Parameter | Type | Default | Description |
|---|---|---|---|
| `min_flat_component_pixels` | int | `500` | Minimum 8-connected component size (pixels) to retain in the flat mask (Stage 3). Removes small isolated flat patches before the mask is inverted. |
| `flat_closing_radius_metres` | float | `500` | Radius (m) for morphological closing of the flat mask (Stage 4). Seals narrow background gaps between adjacent flat patches, such as riser pixels between floodplain segments. Set to `0` to skip. |
| `erosion_radius_metres` | float | `500` | Erosion radius (m) for morphological opening of the flat mask (Stage 5). Severs thin flat peninsulas narrower than twice this radius; broad flat bodies are restored to their exact original extent. Set to `0` to skip. |
| `min_component_pixels` | int | `1000` | Minimum 8-connected component size (pixels) to retain in the mountain mask after inverting and hole-filling (Stage 8). |

---

## Stage 9: Valley infilling

| Parameter | Type | Default | Description |
|---|---|---|---|
| `infill_valley_radius_metres` | float | `0` | Radius (m) of the structuring element for morphological closing of the mountain mask. Seals valley notches cut by incised rivers into the mountain front. Set to approximately half the valley width to be bridged. `0` = disabled. |

---

## Stage 10: Mountain front outlets

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_mountain_front_outlets` | bool | `false` | Compute and write mountain front outlet points. Mountain mask edge pixels whose upstream contributing area meets `min_contributing_pixels` are identified as outlets. Upstream duplicates (outlets that have another outlet downstream on the same flow path) are pruned. |
| `min_contributing_pixels` | int | `1000` | Minimum upstream contributing pixels for a mountain edge pixel to qualify as an outlet. |
| `contributing_pixels_fname` | string | _(empty)_ | Filename of a pre-computed `n_contributing_pixels` raster. If empty, the CLI fills pits and builds D8 flow routing from the DEM automatically. Providing a pre-computed raster skips the fill and routing steps. |

---

## Output parameters

| Parameter | Type | Default | Output file |
|---|---|---|---|
| `write_relief_raster` | bool | `false` | `<write_fname>_relief.tif` — local relief raster. Useful for inspecting the threshold. |
| `write_flat_mask_raster` | bool | `false` | `<write_fname>_flat_mask_closing.tif` and `<write_fname>_flat_mask_opening.tif` — flat mask after Stage 4 and Stage 5 respectively. |
| `write_mountain_mask_raster` | bool | `true` | `<write_fname>_mountain_mask.tif` — final binary mountain mask. 1 = mountain, 0 = flat / piedmont. |
| `write_mountain_labels_raster` | bool | `false` | `<write_fname>_mountain_labels.tif` — each 8-connected mountain component labelled with a unique integer (UInt32). |

When `write_mountain_front_outlets` is true, two files are always written:

- `<write_fname>_mountain_front_outlets.fgb` — outlet points as FlatGeobuf with attributes `row`, `col`, `contributing_pixels`
- `<write_fname>_mountain_front_outlets.csv` — same outlets in WGS84 latitude/longitude

---

## Example parameter file

```
read_path:  /data/dems/
read_fname: my_dem
write_path: /data/output/
write_fname: my_dem_mountain

# Relief window
relief_window_side_pixels: 5

# Adaptive threshold (inspect _relief_cdf.csv to tune)
adaptive_relief_threshold: true

# Flat mask filtering
min_flat_component_pixels:    500
flat_closing_radius_metres:   500
erosion_radius_metres:        500
min_component_pixels:         1000

# Seal valley notches 500 m wide
infill_valley_radius_metres:  250

# Outputs
write_mountain_mask_raster:       true
write_relief_raster:              true
write_mountain_front_outlets:     true
min_contributing_pixels:          1000
```
