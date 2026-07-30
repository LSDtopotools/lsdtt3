# CLI reference: lsdtt-drainage-analysis

Computes drainage-related rasters from a DEM. Supports D8, D-infinity (Tarboton 1997), and MFD (Quinn et al. 1991 / Freeman 1991) drainage area; Topographic Wetness Index (TWI); Soil Topographic Index (STI); sediment connectivity index (Cavalli et al. 2013); and height above channel.

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

## Soil Topographic Index (STI)

STI reuses the same TWI formulation and then applies a denominator term:

- `STI = TWI - ln(Ks*D)` (ksd mode)
- `STI_proxy = TWI - term_spr` (spr mode)
  - `term_spr = ln(1 - spr)` when `sti_spr_transform = ln_1_minus_spr` (default)
  - `term_spr = ln(spr)` when `sti_spr_transform = ln_spr`

Both workflows are designed for Scotland. The required datasets are:
- The Land Cover Scotland 1988 dataset, which can be downloaded from:
  https://www.data.gov.uk/dataset/c373c59e-7b4b-4852-ab1d-553324655917/land-cover-scotland-lcs-1988
- The James Hutton Insitute 1981 National Soil Map at 1:250k scale, which can be downloaded from https://www.hutton.ac.uk/soil-maps/

The national soil map has polygon units, which have soil information in the attribute table. The primary soil series is labelled with a series code `SERCDE1`. There is additional data available (licensed, and needs to be requested) that gives additional information per soil series:
- `SSKIB_C_MED_KSAT.csv`: gives per soil series code a soil profile with saturated hydraulic conductivity values, in cm/day. There are different values depending on the Land Cover Class (this corresponds to the Land Cover Scotland 1988 dataset). This data is used to evaluate the ksd term in ksd mode.
- `SERIES_SPR.csv`: gives for each series code a value that represents the percentage of rainfall that is expected to become surface runoff in a given storm event. This data is used to evaluate the term_spr in spr mode.

Spatial workflow for STI:

1. Read DEM.
2. Clip both STI polygon vectors (soil + land use) to DEM bbox.
3. Rasterize `SERIES_CODE` and `LAND_USE` onto DEM grid.
4. Assemble denominator term per pixel using table joins.
5. Compute and write STI outputs.
6. Optional: threshold STI outputs to keep only pixels at or above the global P95 value.

Expected vector fields:

- Soil layer: `SERIES_CODE` (fallback accepted: `SERCDE1`)
- Land-use layer: `LAND_USE` with values interpreted as `CULT`, `SEMI`, `OTHER`

Expected CSV fields:

- `SSKIB_C_MED_KSAT.csv`: `SERIES_CODE`, `SSKIB_PROFILE_ID`, `LUSE`, `HORZ_TOP`, `HORZ_BOTTOM`, `COND_0_HYPRES`
- `SERIES_SPR.csv`: `SERCDE`, `SPR_prc`

Unit warning (critical): before interpreting absolute STI magnitudes, verify physical units of `COND_0_HYPRES` and horizon depths in your SSKIB tables. Current implementation logs explicit runtime warnings and treats values as proxy/relative if units are unverified.

### Pre-processing of Land Cover and Soil data
The Land Cover Scotland 1988 dataset and the JHI 1981 National Soil Map come as shapefiles with many attributes, but the STI-routine expects flatgeobufs extracted from these datasets.

If you downloaded the source datasets into:

- `/data/scotland_soils/Hutton_Soils_250K_v1.4`
- `/data/scotland_soils/SG_LandCoverScotland_1988`

you can build the two STI-ready FlatGeobuf inputs with:

```bash
# 1) Soil series polygons for STI (`SERIES_CODE` field expected by STI)
ogr2ogr -f FlatGeobuf -nlt PROMOTE_TO_MULTI \
  -nln hutton_series_code_min \
  /data/scotland_soils/hutton_series_code_min.fgb \
  /data/scotland_soils/Hutton_Soils_250K_v1.4/qmsoils_UCSS_v1_3.shp \
  -dialect SQLITE \
  -sql "SELECT *, CAST(SERCDE1 AS TEXT) AS SERIES_CODE FROM qmsoils_UCSS_v1_3 WHERE SERCDE1 IS NOT NULL"

# 2) Land-use polygons with 3-class STI mapping (`LAND_USE` expected by STI)
ogr2ogr -f FlatGeobuf -nlt PROMOTE_TO_MULTI \
  -nln sg_land_use_min \
  /data/scotland_soils/sg_land_use_min.fgb \
  /data/scotland_soils/SG_LandCoverScotland_1988/SG_LandCoverScotland_1988.shp \
  -dialect SQLITE \
  -sql "SELECT *, CASE WHEN LCS_LABEL IN ('arable','improved grassland','recently ploughed land') THEN 'CULT' WHEN LCS_LABEL IN ('blanket bog & peatlands','bracken','broadleaved woodland','coarse grassland','coniferous plantation','duneland','heather moor','low scrub','maritime grasslands & heaths','mixed woodland','montane vegetation','open canopy young plantation','salt marsh','smooth grassland','wetland','woodland recently felled') THEN 'SEMI' ELSE 'OTHER' END AS LAND_USE FROM SG_LandCoverScotland_1988"
```

These outputs match the default STI paths and layer names in this CLI.


### Parameters used by both STI modes

| Parameter | Type | Default | Description |
|---|---|---|---|
| `sti_term_mode` | string | `ksd` | Primary STI mode flag. Options: `ksd`, `spr`. |
| `sti_resample_to_1m_before_calculation` | bool | `false` | If true, GDAL-warps the input DEM to 1 m grid spacing using cubic interpolation before pit filling, flow routing, and STI assembly for this run. The resampled DEM is written to `write_path` and used for all subsequent outputs. The DAERA runoff risk tool resamples their 0.5m resolution topographic data to 1m for sti calculation, see https://doi.org/10.1016/j.jenvman.2024.122110. |
| `write_sti_intermediate_rasters` | bool | `false` | Write debug rasters (`_sti_twi`, `_sti_series_codes`, `_sti_landuse_codes`, term rasters, valid masks). |
| `write_sti_95th_percentile_rasters` | bool | `false` | After writing each requested STI raster, compute its 95th percentile over valid pixels and write a companion raster that keeps only pixels `>= P95` (`_sti_p95.tif` and/or `_sti_proxy_p95.tif`); all lower values become nodata. |
| `write_sti_p95_jenks_rasters` | bool | `false` | Classify each written P95 raster into Jenks natural-break classes and write class rasters (`_sti_p95_jenks5.tif` and/or `_sti_proxy_p95_jenks5.tif` by default, here the number 5 is the default value of `sti_p95_jenks_n_classes`). Requires `write_sti_95th_percentile_rasters: true`. The DEARA runoff risk tool uses these natural breaks classes to translate STI data into a decision support tool that can be iterpreted by stakeholders, see https://doi.org/10.1016/j.jenvman.2024.122110.|
| `sti_p95_jenks_n_classes` | int | `5` | Requested number of Jenks classes for P95 categorisation. |
| `sti_p95_jenks_max_histogram_bins` | int | `2048` | Maximum histogram bins used by Jenks fitting (histogram-weighted implementation). |
| `write_sti_lookup_csv` | bool | `false` | Write `<write_fname>_sti_lookup.csv` with join/debug diagnostics. |
| `sti_soil_vector_path` | string | `data/scotland_soils/hutton_series_code_min.fgb` | Soil polygon vector path. |
| `sti_soil_layer_name` | string | `hutton_series_code_min` | Soil layer name. |
| `sti_landuse_vector_path` | string | `data/scotland_soils/sg_land_use_min.fgb` | Land-use polygon vector path. |
| `sti_landuse_layer_name` | string | `sg_land_use_min` | Land-use layer name. |
| `sti_drainage_area_method` | string | `Dinf` | Drainage area method used by STI TWI stage. Options: `Dinf`, `MFD`. |
| `sti_window_radius` | float | `100.0` | Polyfit window radius (m) for STI slope computation. |

### Parameters for `ksd` mode (`sti_term_mode: ksd`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_sti_raster` | bool | `false` | Write `<write_fname>_sti.tif` using `ln(Ks*D)` term from horizon integration. |
| `sti_ksat_csv_path` | string | `data/scotland_soils/SSKIB_C_MED_KSAT.csv` | Horizon/profile table used for `ln(Ks*D)` term. |
| `sti_xlsx_path` | string | `data/scotland_soils/SSKIB_OCT18_V1.2.xlsx` | Optional consistency source (first pass checks file presence). |
| `sti_depth_cap_cm` | float | `100.0` | Depth cap used when integrating horizons for `Ks*D`. |
| `sti_min_valid_depth_cm` | float | `30.0` | Minimum integrated depth required to keep a profile. |
| `sti_ks_unit_scale_to_m_per_day` | float | `1.0` | Explicit scale factor applied to `COND_0_HYPRES` before horizon integration. |

### Parameters for `spr` mode (`sti_term_mode: spr`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_sti_proxy_raster` | bool | `false` | Write `<write_fname>_sti_proxy.tif` using SPR-based denominator term. |
| `sti_spr_csv_path` | string | `data/scotland_soils/SERIES_SPR.csv` | SPR lookup table by soil series code. |
| `sti_spr_transform` | string | `ln_1_minus_spr` | SPR transform. Options: `ln_1_minus_spr`, `ln_spr`. |

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

## Example parameter file: STI (Scotland workflow)

```
read_path:  /path/to/dem/
read_fname: my_dem
write_path: /path/to/out/
write_fname: my_dem

raster_is_filled: false
fill_min_slope: 0.0001
boundary_conditions: oooo
n_threads: 0
write_output_format: cog

write_sti_raster: true
write_sti_proxy_raster: true
write_sti_95th_percentile_rasters: false
write_sti_intermediate_rasters: false
write_sti_lookup_csv: true

sti_term_mode: ksd
sti_drainage_area_method: Dinf
sti_window_radius: 100.0

sti_soil_vector_path: /mnt/c/LSDTopoTools/data/scotland_soils/hutton_series_code_min.fgb
sti_soil_layer_name: hutton_series_code_min
sti_landuse_vector_path: /mnt/c/LSDTopoTools/data/scotland_soils/sg_land_use_min.fgb
sti_landuse_layer_name: sg_land_use_min

sti_ksat_csv_path: /mnt/c/LSDTopoTools/data/scotland_soils/SSKIB_C_MED_KSAT.csv
sti_spr_csv_path: /mnt/c/LSDTopoTools/data/scotland_soils/SERIES_SPR.csv
sti_xlsx_path: /mnt/c/LSDTopoTools/data/scotland_soils/SSKIB_OCT18_V1.2.xlsx

sti_depth_cap_cm: 100.0
sti_min_valid_depth_cm: 30.0
sti_ks_unit_scale_to_m_per_day: 1.0
sti_spr_transform: ln_1_minus_spr
```
