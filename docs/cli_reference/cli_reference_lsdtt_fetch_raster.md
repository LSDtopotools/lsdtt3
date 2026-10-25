# `lsdtt-fetch-raster` — CLI reference

Downloads a DEM for an area and writes it as an analysis-ready projected
GeoTIFF. No Python environment, and by default no account or API key.

Design notes and rationale: `docs/developer_notes/testing_and_benchmarking.md`.

```sh
lsdtt-fetch-raster <path> <parameter_file>
lsdtt-fetch-raster <parameter_file>      # path assumed to be ./
lsdtt-fetch-raster -h                    # parameter help
lsdtt-fetch-raster -h <parameter>        # detailed help for one parameter
lsdtt-fetch-raster -v                    # version
lsdtt-fetch-raster -c                    # citation
```

---

## Quick start

Two ways to say where you want data. Give **one** of them, never both.

### By coordinates

```
write_path: /home/me/mydata
write_fname: ortler

west:  10.34
south: 46.46
east:  10.46
north: 46.54
```

### By a polygon or two points drawn in a GIS

```
write_path: /home/me/mydata
write_fname: ortler

bounds_fname: my_study_area.gpkg
```

### By outlet points (auto bounds, grouping, final verification)

```text
write_path: /home/me/mydata
write_fname: sulawasi

bounds_mode: outlets
outlet_points_fname: sulawasi_outlets.csv

# defaults are keyless-first
dem_source: cop30_aws

# optional: ask for OpenTopography but fall back keylessly if auth fails
allow_auth_fallback_to_keyless: true
auth_fallback_source: cop30_aws

# split distant/large basin groups into multiple DEMs
auto_split_disconnected_basins: true
split_cluster_gap_km: 25
max_final_dem_area_km2: 2500

# final-resolution bounded recheck
verify_final_resolution_basins: true
final_basin_edge_buffer_pixels: 3
final_basin_nodata_buffer_pixels: 1
final_verify_expand_factor: 1.25
final_verify_max_iters: 3
```

Smoke-test profile (network): keep bounds tiny and force two groups with
`max_basins_per_dem: 1`. Validate these outputs:

- `<write_fname>_group_0001_DEM.tif`
- `<write_fname>_group_0002_DEM.tif`
- `<write_fname>_dem_groups.csv`
- `<write_fname>_outlet_bounds_status.csv`

Then run it and carry straight on into an analysis:

```sh
lsdtt-fetch-raster . fetch.param                # writes ortler_DEM.tif
lsdtt-raster-preprocessing . preprocess.param   # fill, hillshade, ...
```

---

## Output

One GeoTIFF named `<write_fname>_DEM.tif` in `write_path`, plus the usual
`<write_fname>_README.csv` / `.html` parameter documentation every lsdtt3 CLI
writes.

When `bounds_mode: outlets`, output is grouped:

- `<write_fname>_group_0001_DEM.tif`, `<write_fname>_group_0002_DEM.tif`, ...
- `<write_fname>_outlet_bounds_status.csv`
- `<write_fname>_dem_groups.csv` (manifest with `source_used`, verification pass,
  verification iterations, and minimum edge/nodata margins)

The raster is:

- reprojected to a **metric** CRS (UTM by default, derived from the area centre)
- on a fixed grid spacing (30 m by default)
- given an explicit NoData value of −9999

The reprojection happens here rather than downstream, so the output can be fed
directly to any other lsdtt3 tool with `enforce_projection: false`.

---

## Parameters

### Bounding box

| Parameter | Type | Default | Description |
|---|---|---|---|
| `west` | string | — | Western edge, WGS84 decimal degrees |
| `east` | string | — | Eastern edge, WGS84 decimal degrees |
| `south` | string | — | Southern edge, WGS84 decimal degrees |
| `north` | string | — | Northern edge, WGS84 decimal degrees |

All four are required if you use coordinates. A missing edge is an error rather
than a silent zero - a zeroed edge would quietly fetch the Gulf of Guinea.

`bounds_mode` controls how area selection works:

| Parameter | Type | Default | Description |
|---|---|---|---|
| `bounds_mode` | string | `manual` | `manual` (existing behavior), `vector` (require `bounds_fname`), or `outlets` (auto-discover from outlet points). |

Rules:

- `bounds_mode: manual` -> use either explicit bbox or `bounds_fname`
- `bounds_mode: vector` -> `bounds_fname` required; explicit bbox forbidden
- `bounds_mode: outlets` -> `outlet_points_fname` required; explicit bbox and `bounds_fname` forbidden

They are declared as strings, not floats, deliberately: the parameter parser
carries `float`, which quantises to about 1e-7 degrees, and that is metres on
the ground — enough to shift the fetched window.

### Bounds from a vector file

| Parameter | Type | Default | Description |
|---|---|---|---|
| `bounds_fname` | string | *(empty)* | Vector file defining the area |
| `bounds_path` | string | *(empty)* | Directory for `bounds_fname`; falls back to `read_path` |
| `bounds_layer` | string | *(empty)* | Layer name; first layer if empty |

FlatGeobuf (`.fgb`), GeoPackage (`.gpkg`), or anything else OGR reads. The layer
must contain **either**:

- **one or more polygons** — the envelope of all of them is used, or
- **exactly two points** — opposite corners of the area.

**Either diagonal works for two points.** Lower-left + upper-right and
lower-right + upper-left give the same box, because it is just the min and max
of the two x and y values. There is nothing to get the wrong way round.

A layer in a projected CRS (UTM, British National Grid, …) is reprojected to
WGS84 automatically, so a box drawn in a GIS works directly. A file with no
declared CRS is assumed to be WGS84 degrees, with a warning.

Rejected, with an explanation: one point, three or more points, a mix of
polygons and points, and other geometry types.

Setting `bounds_fname` **and** any of `west`/`south`/`east`/`north` is an error.
Silently preferring one would leave you unsure which area you actually got.

### Outlet-based bounds discovery

| Parameter | Type | Default | Description |
|---|---|---|---|
| `outlet_points_fname` | string | *(empty)* | Outlet point file for `bounds_mode: outlets` (CSV/FGB/GPKG/SHP and other GDAL vector point formats). |
| `outlet_points_path` | string | *(empty)* | Optional directory for `outlet_points_fname`; falls back to `read_path`. |
| `prefetch_source` | string | `cop30_aws` | Source used for coarse first-guess basin discovery (`cop30_aws`, `planetary_computer`, or `opentopography`). |
| `prefetch_grid_spacing` | float | `90` | Grid spacing (m) used in coarse prefetch routing. |
| `prefetch_initial_radius_km` | float | `100` | Initial radius for first-guess outlet bounds. |
| `prefetch_expand_factor` | float | `1.8` | Expansion factor when retrying bounds. Must be >1. |
| `prefetch_max_iters` | int | `6` | Max first-guess expansion iterations. |
| `prefetch_edge_buffer_pixels` | int | `3` | Coarse-stage edge buffer target. |
| `basin_outlet_snapping_method` | string | `downstream` | Outlet snapping for basin delineation: `none`, `expanding_square`, `downstream`, or `largest_area`. `largest_area` takes the pixel with the most contributing pixels in the square of half-width `basin_outlet_snap_max_search_radius_pixels`, if it reaches the threshold (a wide square can jump onto a bigger neighbouring river). |
| `basin_outlet_contributing_pixels_threshold` | int | `1000` | Minimum contributing-pixel threshold used during outlet snapping. |
| `basin_outlet_snap_min_search_radius_pixels` | int | `0` | Minimum search radius for `expanding_square` snapping. |
| `basin_outlet_snap_max_search_radius_pixels` | int | `25` | Maximum search radius for `expanding_square` snapping, and the square's half-width for `largest_area`. |
| `tight_basin_bbox_refetch` | bool | `true` | After basin delineation, refetch final group DEMs using tight basin bounds. |
| `tight_basin_bbox_buffer_pixels` | int | `3` | Pixel buffer around basin envelope when forming tight refetch bounds. |
| `basin_nested_enable` | bool | `true` | Enable staged outlet-mode nested fetching pipeline. |
| `basin_nested_stop_level` | string | `30m` | Stop level for staged pipeline: `coarse`, `90m`, or `30m`. |
| `coarse_prefetch_demtype` | string | `GEBCOIceTopo` | OT demtype for the coarse prefetch stage when OT prefetch is active. |
| `coastal_mask_mode` | string | `clip_bbox_edges` | `off` or `clip_bbox_edges` pre-download coastal bbox contraction. |
| `coastal_mask_dataset` | string | `natural_earth_110m_land` | Default managed land-mask dataset id (used when no file override is set). |
| `coastal_mask_auto_download` | bool | `false` | If default mask is missing, download Natural Earth `ne_110m_land.zip` into `read_path/coastal_masks/`. |
| `coastal_mask_fname` | string | *(empty)* | Optional land-mask polygon file override. |
| `coastal_mask_path` | string | *(empty)* | Optional directory for `coastal_mask_fname`; falls back to `read_path`. |
| `coastal_mask_layer` | string | *(empty)* | Optional layer name for `coastal_mask_fname`; first layer if empty. |
| `coastal_mask_edge_buffer_km` | float | `2` | Safety margin after coastal edge contraction. Must be >=0. |
| `discard_intermediate_rasters` | bool | `true` | Remove temporary intermediate rasters from staged outlet-mode runs. |
| `write_basin_outlines_fgb` | bool | `true` | Write basin outlines as FlatGeobuf in outlet mode. |
| `write_fetch_extents_manifest` | bool | `true` | Write manifest of exact fetch extents used per stage/group. |

The tool writes per-outlet status to `<write_fname>_outlet_bounds_status.csv`.

When `coastal_mask_mode: clip_bbox_edges`, the tool attempts to contract each
outlet's radius-first bbox against land polygons **before** DEM download. If the
mask is unavailable, it warns and falls back to radius-only bounds.

Default land-mask lookup (`coastal_mask_dataset: natural_earth_110m_land`) is:

- `read_path/ne_110m_land.*`
- `read_path/coastal_masks/ne_110m_land.*`
- `read_path/masks/ne_110m_land.*`

with `*.fgb`, `*.gpkg`, `*.shp`, or `*.zip`.

If none of those exist and `coastal_mask_auto_download: true`, the tool tries
to download:

- `https://naturalearth.s3.amazonaws.com/110m_physical/ne_110m_land.zip`

into `read_path/coastal_masks/ne_110m_land.zip` and then uses it as the mask.

Auto-download example (first run on a new machine):

```text
read_path: ./input
write_path: ./output
write_fname: coast_prefetch

bounds_mode: outlets
outlet_points_fname: outlets.fgb

coastal_mask_mode: clip_bbox_edges
coastal_mask_dataset: natural_earth_110m_land
coastal_mask_auto_download: true
coastal_mask_edge_buffer_km: 2

dem_source: cop30_aws
target_epsg: 32632
grid_spacing: 30
```

Coastal clipping examples:

```text
# Example A: use default Natural Earth lookup in read_path
bounds_mode: outlets
outlet_points_fname: outlets.csv

coastal_mask_mode: clip_bbox_edges
coastal_mask_dataset: natural_earth_110m_land
coastal_mask_auto_download: true
coastal_mask_edge_buffer_km: 2
```

```text
# Example B: override with your own land polygon mask
bounds_mode: outlets
outlet_points_fname: outlets.gpkg

coastal_mask_mode: clip_bbox_edges
coastal_mask_fname: my_land_mask.gpkg
coastal_mask_path: /home/me/masks
coastal_mask_layer: land
coastal_mask_edge_buffer_km: 1
```

#### Outlet points as FlatGeobuf (`.fgb`)

You can provide outlets directly as FlatGeobuf point features.

```text
read_path: ./input
write_path: ./output
write_fname: ortler_outlets

bounds_mode: outlets
outlet_points_fname: outlets_ortler.fgb

# final DEM settings
dem_source: cop30_aws
target_epsg: 32632
grid_spacing: 30
```

Requirements for outlet vectors:

- Geometry must be points.
- Any CRS is accepted if it is declared; points are transformed to WGS84
  internally.
- For `.fgb`, the first (and usually only) layer is used.

#### Basin workflow when you do not know the basin extent

Use `bounds_mode: outlets` and let the tool discover a safe first envelope around
each outlet, then expand until routed basin pixels meet edge-buffer checks.

Important: outlet-mode basin delineation now snaps each provided outlet to a
channel cell by default (`basin_outlet_snapping_method: downstream`) using
`basin_outlet_contributing_pixels_threshold` (default `1000`). This avoids
single-pixel micro-basins caused by unsnapped hillslope points.

Recommended sequence:

1. Start coarse for discovery (`prefetch_grid_spacing: 90` is the default).
2. If you have an OpenTopography key, set `prefetch_source: opentopography` and
   use `coarse_prefetch_demtype` for the coarse stage (default:
   `GEBCOIceTopo`).
3. Keep final output on your preferred source/resolution (`dem_source` and
   `grid_spacing`), with `verify_final_resolution_basins: true`.
4. If any group still fails edge/nodata buffers, raise
   `prefetch_initial_radius_km` and/or `prefetch_max_iters` first.

When nested prefetch is enabled, the coarse-stage fetch writes
`<write_fname>_coarse_prefetch_DEM.tif` before final grouping/refinement. If
`basin_nested_stop_level: coarse`, the workflow exits after that coarse output.
If `discard_intermediate_rasters: true`, this coarse raster is removed when
continuing to later stages.

For grouped refinement:

- `basin_nested_stop_level: 90m` writes group outputs at 90 m and exits.
- `basin_nested_stop_level: 30m` runs a 90 m intermediate stage first, then
  final group fetches at `grid_spacing`.
- Intermediate 90 m stage rasters use
  `<write_fname>_group_####_90m_DEM.tif` and are removed when
  `discard_intermediate_rasters: true`.
- With `tight_basin_bbox_refetch: true` (default), once a basin is delineated,
  final group fetches are reissued on a tight basin envelope with
  `tight_basin_bbox_buffer_pixels` around it (default `3` pixels).
  For a one-pixel basin this gives a `7 x 7` raster footprint before any
  projection/resampling effects.

Example profile (unknown basin extent, OT key available):

```text
bounds_mode: outlets
outlet_points_fname: outlets.csv

# coarse basin discovery pass
prefetch_source: opentopography
coarse_prefetch_demtype: GEBCOIceTopo
prefetch_grid_spacing: 90

# final DEM pass
dem_source: cop30_aws
grid_spacing: 30
verify_final_resolution_basins: true

# snap outlets to channels and keep final fetch tight to basin shape
basin_outlet_snapping_method: downstream
basin_outlet_contributing_pixels_threshold: 1000
tight_basin_bbox_refetch: true
tight_basin_bbox_buffer_pixels: 3
```

### Group splitting and manifests (`bounds_mode: outlets`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `auto_split_disconnected_basins` | bool | `true` | Split distant basin groups into separate DEMs. |
| `split_cluster_gap_km` | float | `25` | Distance threshold for clustering outlets/basins. |
| `max_final_dem_area_km2` | float | `2500` | Maximum approximate area per final DEM group. |
| `split_overlap_km` | float | `2` | Overlap buffer applied around each group bbox. |
| `max_basins_per_dem` | int | `0` | Optional cap per group (`0` = unlimited). |

How grouping decides to split into 2+ rasters:

1. **Disconnected outlets first (optional):** if
   `auto_split_disconnected_basins: true`, outlets are clustered by
   outlet-to-outlet distance. Any pair farther apart than `split_cluster_gap_km`
   ends up in different clusters, so separate rasters are produced.
2. **Size/cap checks per cluster:** each cluster's merged bbox is checked
   against:
   - `max_final_dem_area_km2`
   - `max_basins_per_dem` (if > 0)
3. **Recursive split when too large:** if either limit is exceeded, the cluster
   is split into two parts by median outlet longitude, then each part is checked
   again recursively.
4. **Overlap padding:** once final groups are accepted, each group bbox is
   expanded by `split_overlap_km` to avoid seam-edge clipping.

So, in practice, you get 2 rasters when either:

- two outlet sets are farther apart than `split_cluster_gap_km` (with
  `auto_split_disconnected_basins: true`), or
- one cluster exceeds area/basin-count limits and the recursive split produces
  two valid subgroups.

If you want fewer splits, increase `split_cluster_gap_km`, increase
`max_final_dem_area_km2`, or set `max_basins_per_dem: 0`.

Example tuned to encourage splitting into multiple rasters:

```text
bounds_mode: outlets
outlet_points_fname: outlets_ortler.fgb

auto_split_disconnected_basins: true
split_cluster_gap_km: 15
max_final_dem_area_km2: 800
max_basins_per_dem: 2
split_overlap_km: 1

dem_source: cop30_aws
target_epsg: 32632
grid_spacing: 30
```

Example tuned to keep more basins in one raster:

```text
bounds_mode: outlets
outlet_points_fname: outlets_ortler.fgb

auto_split_disconnected_basins: true
split_cluster_gap_km: 80
max_final_dem_area_km2: 6000
max_basins_per_dem: 0
split_overlap_km: 2

dem_source: cop30_aws
target_epsg: 32632
grid_spacing: 30
```

Manifest output `<write_fname>_dem_groups.csv` records one row per final group
including requested/used source, bbox, area, and final verification metrics.

If `write_fetch_extents_manifest: true`,
`<write_fname>_fetch_extents_manifest.csv` records each staged fetch attempt
(coarse prefetch, group 90 m, group final) with stage label, attempt number,
source requested/used, spacing, exact bbox, area, and status.

If `write_basin_outlines_fgb: true`, basin outlines are written per group as
`<write_fname>_group_####_basin_outlines.fgb` at the retained stage
(`90m` stop-level outputs or final-stage outputs).

### Final-resolution basin safety verification

| Parameter | Type | Default | Description |
|---|---|---|---|
| `verify_final_resolution_basins` | bool | `true` | Re-delineate/verify on the final-resolution DEM before accepting output. |
| `final_basin_edge_buffer_pixels` | int | `3` | Required margin between any selected basin pixel and raster edge. |
| `final_basin_nodata_buffer_pixels` | int | `1` | Required margin between any selected basin pixel and nodata. |
| `final_verify_expand_factor` | float | `1.25` | Bbox growth factor for failed verification retries. |
| `final_verify_max_iters` | int | `3` | Max verification refetch retries per group. |

Verification is basin-pixel-wide (not outlet-only): if any basin pixel violates
edge/nodata buffers, that group is refetched with expanded bounds until pass or
retry cap.

### Source

| Parameter | Type | Default | Description |
|---|---|---|---|
| `dem_source` | string | `cop30_aws` | Which provider to fetch from |
| `ot_dem_type` | string | `COP30` | OpenTopography dataset (that source only); see the table below |
| `pc_collection` | string | `cop-dem-glo-30` | Planetary Computer STAC collection (that source only) |
| `pc_asset` | string | `data` | Asset key within each STAC item |

| `dem_source` | Key needed? | Notes |
|---|---|---|
| `cop30_aws` | **no** | Copernicus GLO-30 on AWS Open Data. Global 30 m. The default; start here. |
| `planetary_computer` | **no** | Microsoft Planetary Computer. Asset URLs are signed automatically. Opens up many other STAC collections. |
| `opentopography` | **yes** | OpenTopography REST API. Free key; needed for its other datasets. |

Anything else is an error, never a silent substitution to a different provider.

`cop30_aws` and `planetary_computer` serve the same underlying Copernicus data;
verified to produce byte-identical output. Prefer `planetary_computer` when you
want a collection other than GLO-30 (`nasadem`, `alos-dem`, `3dep-seamless`, …;
browse <https://planetarycomputer.microsoft.com/catalog>).

#### OpenTopography datasets (`ot_dem_type`)

These are the `demtype` values the OpenTopography
[`/globaldem`](https://portal.opentopography.org/apidocs/#/Public/getGlobalDem)
endpoint serves. All of them need the free API key. A value not in this list is
rejected before the request is sent, with the list printed — the API itself
answers a bad `demtype` with a bare 400.

Run `lsdtt-fetch-raster -h ot_dem_type` for the same table with citations.

| `ot_dem_type` | Dataset | Band type (`GTiff`) | Grid spacing | Coverage | Notes |
|---|---|---|---|---|---|
| `COP30` | Copernicus GLO-30 | `Float32` | 1 arc-second (~30 m) | global | TanDEM-X **DSM**, edited water bodies and coastlines. The default, and the same data `cop30_aws` serves keylessly. |
| `COP90` | Copernicus GLO-90 | `Float32` | 3 arc-second (~90 m) | global | The 90 m Copernicus DSM. Use it when the area exceeds the 30 m request limit. |
| `SRTMGL1` | SRTM GL1 | `Int16` | 1 arc-second (~30 m) | 60 N–56 S | 1 arc-second SRTM, EGM96 heights. The long-standing comparison baseline. |
| `SRTMGL1_E` | SRTM GL1 Ellipsoidal | `Float32` | 1 arc-second (~30 m) | 60 N–56 S | As `SRTMGL1` but WGS84 ellipsoidal heights — for co-registering with GNSS or lidar. |
| `SRTMGL3` | SRTM GL3 | `Int16` | 3 arc-second (~90 m) | 60 N–56 S | 3 arc-second SRTM. Larger request limit than GL1. |
| `NASADEM` | NASADEM | `Int16` | 1 arc-second (~30 m) | 60 N–56 S | SRTM reprocessed: better geolocation, voids filled from ASTER GDEM and ICESat. Generally preferable to raw `SRTMGL1`. |
| `AW3D30` | ALOS World 3D | `Int16` | 1 arc-second (~30 m) | global | JAXA optical-stereo DSM. An *independent* global 30 m surface — useful as a cross-check on the radar products. |
| `AW3D30_E` | ALOS World 3D Ellipsoidal | `Float32` | 1 arc-second (~30 m) | global | `AW3D30` with ellipsoidal rather than geoid heights. |
| `GEDTM30` | Global Ensemble DTM | `Float32` | 1 arc-second (~30 m) | global | Machine-learning ensemble **DTM** — bare earth, vegetation and buildings removed. |
| `EU_DTM` | Continental Europe DTM | `Float32` | 30 m (EPSG:3035 output grid) | Europe | Bare-earth DTM from national lidar and photogrammetric surveys. |
| `ANADEM` | ANADEM | `Float32` | ~30 m (EPSG:4674, API sample grid `0.000269494585` deg) | South America | Bare-earth DTM; vegetation bias removed from Copernicus GLO-30. |
| `CA_MRDEM_DTM` | Canada MRDEM (terrain) | `Float32` | 30 m (EPSG:3979 output grid) | Canada | Natural Resources Canada medium-resolution bare-earth model. |
| `CA_MRDEM_DSM` | Canada MRDEM (surface) | `Float32` | 30 m (EPSG:3979 output grid) | Canada | The surface-model companion to `CA_MRDEM_DTM`. |
| `SRTM15Plus` | SRTM15+ V2.1 | `Float32` | 15 arc-second (~500 m) | global, land and sea | Merged topography and bathymetry. Use it when the domain continues offshore. |
| `GEBCOIceTopo` | GEBCO 2023 (ice surface) | `Int16` | 15 arc-second (~500 m) | global, land and sea | Bathymetry and topography, with the *ice-sheet surface* over Greenland and Antarctica. |
| `GEBCOSubIceTopo` | GEBCO 2023 (bed under ice) | `Int16` | 15 arc-second (~500 m) | global, land and sea | As above, but the ice sheets are replaced by the underlying bed. |
| `GEDI_L3` | GEDI L3 gridded metrics | `Float32` | nominal 1000 m (API sample in EPSG:4326: `0.010573500689` deg) | 52 N–52 S | Gridded GEDI lidar land-surface metrics. A metrics product, **not** a DEM for flow routing. |

Band types and sample grid spacing values above were verified from direct
`/globaldem` `GTiff` responses (tiny bounding boxes, 2026-07-30).

**DSM or DTM matters.** `COP30`, `AW3D30` and the `CA_MRDEM_DSM` product are
*surface* models: they include canopy and buildings, so hillslope curvature and
roughness over forest are vegetation, not ground. `GEDTM30`, `EU_DTM`, `ANADEM`
and `CA_MRDEM_DTM` are bare-earth *terrain* models over their coverage.

**Per-request area limits** are set by the API, not by this tool: 450 000 km² for
the 30 m products, 4 050 000 km² for `COP90` and `SRTMGL3`, 125 000 000 km² for
the 500 m bathymetry products, and 500 000 000 km² for `GEDI_L3`. Exceeding one
comes back as an API error naming the limit.

#### Citation is printed for you

Every OpenTopography dataset is a DOI-minted distribution whose provider
requires attribution. After a successful `dem_source: opentopography` run the
tool prints the citation for the dataset you actually fetched, together with the
OpenTopography NSF acknowledgement — so it lands in your terminal and your log
rather than needing to be looked up months later:

```
--- Please cite this data ---
SRTM GL1 (SRTMGL1)
NASA Shuttle Radar Topography Mission (SRTM) (2013). Shuttle Radar Topography
Mission (SRTM) Global. Distributed by OpenTopography.
https://doi.org/10.5069/G9445JDF

Also acknowledge OpenTopography, which is supported by the US National Science
Foundation (NSF awards 2410799, 2410800, 2410801).
```

The citation is printed on a cached reuse too: the obligation follows the data,
not the download.

### Output grid

| Parameter | Type | Default | Description |
|---|---|---|---|
| `target_epsg` | int | `0` | Output EPSG. `0` derives the UTM zone from the area centre |
| `grid_spacing` | float | `30.0` | Grid spacing in projection units (metres for UTM) |
| `resampling_method` | string | `cubic` | GDAL warp method. Cubic is the default for continuous fetched DEM elevations, not categorical rasters. Options: `near`, `bilinear`, `cubic`, `cubicspline`, `lanczos`, `average`, `mode`, `max`, `min`, `med`, `q1`, `q3`. |
| `size_pixels` | int | `0` | Crop an exact square of this many pixels from the centre. `0` keeps the full extent |
| `no_data_value` | float | `-9999.0` | NoData value written into the output |

Set `target_epsg` explicitly when fetching several adjacent areas that must
share a grid — otherwise each derives its own UTM zone and neighbouring areas
can land in different projections.

`size_pixels` makes output dimensions deterministic, which matters for test
fixtures and for comparing runs. The bounding box must be large enough to
contain the crop, or you get an error telling you to widen it.

`no_data_value` matters more than it looks: cells outside the source footprint
are set to it rather than to zero. Zero would read as sea level and behave as a
large spurious depression in flow routing.

### Network

| Parameter | Type | Default | Description |
|---|---|---|---|
| `http_timeout_seconds` | int | `60` | Maximum seconds for one HTTP request |
| `http_connect_timeout_seconds` | int | `15` | Maximum seconds to establish a connection |
| `http_max_retry` | int | `3` | Retries on a transient failure |
| `cache_dir` | string | *(empty)* | Optional cache directory; empty means no caching |

Keep the connect timeout short so an offline machine fails promptly rather than
hanging. Cloud storage occasionally returns 5xx responses; the retries absorb
those.

---

## API keys

**Only `opentopography` needs one.** The default source and Planetary Computer
need no account at all.

Get a free key from <https://opentopography.org>, then supply it in one of these
ways — checked in this order:

```sh
# 1. environment variable
export OT_API_KEY="your_key"

# 2. or a key file, readable only by you
echo "your_key" > ~/.ot_api_key
chmod 600 ~/.ot_api_key
export OT_API_KEY_FILE="$HOME/.ot_api_key"

# 3. or let it prompt you (input is hidden). Only when running interactively.
```

`OT_API_KEY` is the same variable the `python_helpers/` scripts and
lsdviztools3 use, so an existing setup already works.

### The key is never a parameter — on purpose

> **Do not put an API key in a parameter file.** `lsdtt-fetch-raster` will not
> read one from there, and this is deliberate rather than an oversight.
> Parameter files are routinely committed alongside analyses, and every lsdtt3
> CLI echoes all its parameters at the end of a successful run. A key held as a
> parameter would end up in your repository and in your logs.

For the same reasons the key is never accepted as a command-line argument
(`ps` shows other users your argv, and it lands in shell history) and is never
compiled into the binary.

Keys are masked wherever they might be printed. This is not theoretical: the
OpenTopography API **echoes the submitted key back inside its own error
message**, so a mistyped key would otherwise be printed to your terminal and
into any redirected log. You will see `***` instead.

A key file readable by other users produces a warning suggesting `chmod 600`.

### Quick OpenTopography smoke test (after `source`-ing your key)

If you have just exported/sourced your key, validate the auth path with a tiny
request before a larger run:

```text
write_path: ./dem
write_fname: ot_smoke

west: 10.3373594105
south: 46.4568810636
east: 10.3626405895
north: 46.4818810636

dem_source: opentopography
ot_dem_type: SRTMGL1
target_epsg: 32632
grid_spacing: 30
size_pixels: 128

http_connect_timeout_seconds: 10
http_timeout_seconds: 30
```

```sh
lsdtt-fetch-raster . ot_smoke.param
```

Expected result: `ot_smoke_DEM.tif` is written, the fetch banner shows
`Source      : opentopography`, and the dataset citation block is printed.

### Copernicus Data Space (CDSE) is not supported

Its only distinctive product is EU-DEM at 10 m, which a general public account
does not unlock. Supporting it would introduce the tool's only
username/password credential path in exchange for coverage already available
keylessly from `cop30_aws`.

---

## Worked examples

### A study area from a polygon drawn in QGIS

```
read_path: /home/me/gis
write_path: /home/me/dem
write_fname: catchment

bounds_fname: catchment_outline.gpkg
grid_spacing: 30.0
```

The polygon can be in any CRS; it is reprojected for you.

### A fixed-size tile for a reproducible comparison

```
write_path: ./fixtures
write_fname: alpine

west:  10.3373594105
south: 46.4568810636
east:  10.4626405895
north: 46.5431189364

target_epsg: 32632
grid_spacing: 30.0
size_pixels: 256
```

Pinning `target_epsg` and `size_pixels` makes the output deterministic for a
given GDAL build.

### NASADEM via Planetary Computer, still keyless

```
write_path: ./dem
write_fname: study

bounds_fname: area.fgb

dem_source: planetary_computer
pc_collection: nasadem
```

### SRTM 1-arcsecond via OpenTopography

```
write_path: ./dem
write_fname: study

west:  10.34
south: 46.46
east:  10.46
north: 46.54

dem_source: opentopography
ot_dem_type: SRTMGL1
```

with `OT_API_KEY` set in the environment — not in this file.

---

## Troubleshooting

**`COULD NOT FETCH: ... could not open source tile(s)`**
No network, a proxy or firewall blocking HTTPS to the data host, or a bounding
box outside the source's coverage. The message names the tiles it looked for.

**`FATAL: no area given`**
Set the four edges, or `bounds_fname`.

**`FATAL: bounds_mode=outlets requires outlet_points_fname`**
Set an outlet file when using outlet-based discovery.

**`FATAL: group <id> failed final-resolution verification ...`**
Final DEM did not meet edge/nodata basin-pixel buffers within retry cap.
Increase initial radius or verification retries, or reduce strictness carefully.

**`FATAL: give either bounds_fname or west/south/east/north, not both`**
Exactly what it says — remove one.

**`FATAL: ... exactly two are needed, marking opposite corners`**
The points layer has one point, or three or more. Either diagonal is fine, but
it must be two.

**`FATAL: warped raster is NxM, smaller than the requested ... centre crop`**
`size_pixels` is larger than the area covers. Widen the bounding box or reduce
`size_pixels`.

**`FATAL: no credential found`**
`dem_source: opentopography` without a key. Set `OT_API_KEY`, or use
`cop30_aws`, which needs no key and gives the same global 30 m coverage.

**`FATAL: ... unknown ot_dem_type '...'`**
A typo, or a dataset that is not on the `/globaldem` endpoint. The message lists
every valid value; `-h ot_dem_type` describes them. This is caught locally, so
no request is sent.

**`the OpenTopography API returned an error: ... API Key: ***`**
The key was rejected. The `***` is the masking working as intended — your key
is not being printed.

**`could not open '...' as a vector dataset`**
`bounds_fname` is not a readable vector file. Check the path (it resolves
against `bounds_path`, then `read_path`) and that the format is one GDAL builds
support.

---

## See also

- `docs/developer_notes/testing_and_benchmarking.md` — design rationale, credential rules, findings
- `docs/cli_reference/cli_reference_lsdtt_raster_preprocessing.md` — the usual next step
- `python_helpers/download_dem_tiles.py` — Python helper covering other providers
- <https://portal.opentopography.org/apidocs/#/Public/getGlobalDem> — the
  upstream `demtype` list this tool's catalogue mirrors
- <https://planetarycomputer.microsoft.com/catalog> — collections for
  `pc_collection`

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

### Bounding box

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `west` | `string` | - | `param` | - | Western edge of the bounding box, WGS84 decimal degrees. | Required. Must be less than east. Range -180 to 180. |
| `east` | `string` | - | `param` | - | Eastern edge of the bounding box, WGS84 decimal degrees. | Required. Must be greater than west. Range -180 to 180. |
| `south` | `string` | - | `param` | - | Southern edge of the bounding box, WGS84 decimal degrees. | Required. Must be less than north. Range -90 to 90. |
| `north` | `string` | - | `param` | - | Northern edge of the bounding box, WGS84 decimal degrees. | Required. Must be greater than south. Range -90 to 90. |
| `bounds_mode` | `string` | `manual` | `param` | - | How the fetch area is specified. | manual: existing behavior (west/south/east/north or bounds_fname). vector: require bounds_fname and forbid explicit edges. outlets: E1 first-guess outlet mode using radius-based bounds from outlet points, with status CSV output. |
| `bounds_fname` | `string` | - | `param` | - | Vector file (or raster) defining the area, instead of west/south/east/north. | FlatGeobuf (.fgb), GeoPackage (.gpkg) or any other OGR-readable format. The layer must contain either one or more polygons (the envelope of all of them is used), or exactly two points marking opposite corners. For two points either diagonal works -- lower-left/upper-right and lower-right/upper-left give the same box, so the orientation is worked out automatically. If the layer is in a projected CRS such as UTM it is reprojected to WGS84 for you, so a box drawn in a GIS can be used directly. A raster with a georeference (e.g. an existing DEM GeoTIFF) works too: its footprint is reprojected to WGS84. Set this or the four edge parameters, not both. |
| `bounds_path` | `string` | - | `param` | - | Optional directory for bounds_fname. | If empty, read_path is used. Ignored when bounds_fname is an absolute path. |
| `bounds_layer` | `string` | - | `param` | - | Layer name within bounds_fname. | If empty, the first layer is used (with a warning when the file has more than one). |
| `outlet_points_fname` | `string` | - | `param` | - | Outlet point file used when bounds_mode is outlets. | CSV/FGB/GPKG. Required for bounds_mode: outlets. |
| `outlet_points_path` | `string` | - | `param` | - | Optional directory for outlet_points_fname. | If empty, read_path is used. Ignored when outlet_points_fname is an absolute path. |
| `prefetch_initial_radius_km` | `float` | `100.0` | `param` | - | Initial prefetch radius around each outlet (km). | Expanded iteratively when discovered basins touch edges. |
| `basin_nested_enable` | `bool` | `true` | `param` | - | Enable staged outlet-mode nested fetching. | When true, outlet mode follows the staged coarse-&gt;90 m-&gt;30 m workflow (subject to basin_nested_stop_level). |
| `basin_nested_stop_level` | `string` | `30m` | `param` | - | Where outlet-mode nested fetching stops. | One of: coarse, 90m, 30m. |
| `coastal_mask_mode` | `string` | `clip_bbox_edges` | `param` | - | How outlet mode applies land/coast masking before fetch. | One of: off, clip_bbox_edges. |
| `coastal_mask_dataset` | `string` | `natural_earth_110m_land` | `param` | - | Default built-in/managed land-mask dataset identifier. | Used when coastal_mask_fname is empty. |
| `coastal_mask_auto_download` | `bool` | `false` | `param` | - | Auto-download default Natural Earth land mask when missing. | Only applies to coastal_mask_dataset=natural_earth_110m_land; writes to read_path/coastal_masks/. |
| `coastal_mask_fname` | `string` | - | `param` | - | Optional land-mask vector file overriding coastal_mask_dataset. | Expected polygonal land areas in WGS84-compatible coordinates. |
| `coastal_mask_path` | `string` | - | `param` | - | Optional directory for coastal_mask_fname. | If empty, read_path is used. Ignored for absolute paths. |
| `coastal_mask_layer` | `string` | - | `param` | - | Layer name inside coastal_mask_fname. | If empty, the first layer is used. |
| `coastal_mask_edge_buffer_km` | `float` | `2.0` | `param` | - | Safety margin added after coastal bbox edge contraction (km). | Must be &gt;= 0. |
| `prefetch_expand_factor` | `float` | `1.8` | `param` | - | Expansion multiplier for iterative coarse prefetch. | Must be &gt; 1.0. |
| `prefetch_max_iters` | `int` | `6` | `param` | - | Maximum prefetch expansion iterations per outlet/group. | Must be &gt;= 1. |
| `prefetch_edge_buffer_pixels` | `int` | `3` | `param` | - | Minimum coarse-stage basin margin from raster edge (pixels). | Must be &gt;= 0. |
| `basin_outlet_snapping_method` | `string` | `downstream` | `param` | - | Snapping method for outlet-mode basin delineation. | One of: none, expanding_square, downstream, largest_area. largest_area takes the pixel with the most contributing pixels within basin_outlet_snap_max_search_radius_pixels. |
| `basin_outlet_contributing_pixels_threshold` | `int` | `1000` | `param` | - | Contributing-pixel threshold used for outlet snapping. | Higher values snap outlets onto better-defined channels. |
| `basin_outlet_snap_min_search_radius_pixels` | `int` | `0` | `param` | - | Minimum expanding-square snapping radius in pixels. | Used when basin_outlet_snapping_method=expanding_square. |
| `basin_outlet_snap_max_search_radius_pixels` | `int` | `25` | `param` | - | Maximum expanding-square snapping radius in pixels. | Used when basin_outlet_snapping_method=expanding_square or largest_area (the half-width of the square searched). |
| `tight_basin_bbox_buffer_pixels` | `int` | `3` | `param` | - | Pixel buffer around delineated basin for tight refetch bounds. | A value of 3 gives a 7x7 fetch for a 1-pixel basin. |
| `tight_basin_bbox_refetch` | `bool` | `true` | `param` | - | Refetch group DEMs using tight basin bounds after delineation. | Applied in bounds_mode=outlets. |
| `auto_split_disconnected_basins` | `bool` | `true` | `param` | - | Auto-split distant/disconnected basins into multiple final DEMs. | Used in bounds_mode: outlets. |
| `split_cluster_gap_km` | `float` | `25.0` | `param` | - | Distance threshold for basin grouping (km). | Must be &gt; 0 when auto_split_disconnected_basins is true. |
| `max_final_dem_area_km2` | `float` | `2500.0` | `param` | - | Maximum target area per final DEM group (km^2). | Must be &gt; 0. |
| `split_overlap_km` | `float` | `2.0` | `param` | - | Overlap margin between adjacent split groups (km). | Must be &gt;= 0. |
| `max_basins_per_dem` | `int` | `0` | `param` | - | Operational cap on number of basins per final DEM group. | 0 means unlimited. Must be &gt;= 0. |

### Linear feature frame

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fetch_from_linear_feature` | `bool` | `false` | `switch` | - | Fetch along a feature (a line, two points or a polygon) in fetch_feature_fname. | The area is a rectangle aligned with the feature: for a line or two points, the geodesic between the end points (only the end points of a curved or multipart line are used) lengthened by the end buffer and linear_feature_fetch_width wide; for polygons, their minimum-area enclosing rectangle, whose longest side is the axis. With rotate_after_fetch the DEM is written in an oblique Mercator frame with the feature along +x. Cannot be combined with west/south/east/north, bounds_fname, bounds_mode vector or outlets, or size_pixels. |
| `fetch_feature_fname` | `string` | - | `io` | `fetch_from_linear_feature` | Vector file holding the line, two points or polygon(s). | Any OGR-readable format (GeoPackage, FlatGeobuf, shapefile, ...) with a CRS; it is reprojected to WGS84. Refused before any download: no CRS, a mix of geometry kinds, one point or more than two, a closed line (use a polygon), or end points less than three cells apart. |
| `fetch_feature_path` | `string` | - | `io` | `fetch_from_linear_feature` | Optional directory for fetch_feature_fname. | If empty, read_path is used. Ignored when fetch_feature_fname is an absolute path. |
| `fetch_feature_layer` | `string` | - | `io` | `fetch_from_linear_feature` | Layer name within fetch_feature_fname. | If empty, the first layer is used (with a warning when the file has more than one). |
| `linear_feature_fetch_width` | `float` | `0.0` | `param` | `fetch_from_linear_feature` | Full width (m) of the area across a line or two points, centred on it. | Required (&gt; 0) for a line or two points. Ignored, with a warning, for polygons, whose own short side is the width. |
| `linear_feature_fetch_end_buffer` | `float` | `-1.0` | `param` | `fetch_from_linear_feature` | Buffer (m) added beyond each end of the axis. | -1 means 10 grid cells. For polygons it is added on all four sides. |
| `rotate_after_fetch` | `bool` | `true` | `param` | `fetch_from_linear_feature` | Write the DEM in the feature's oblique frame (feature along +x). | true: the DEM is warped straight from the source tiles into the frame (one resampling) and written as &lt;write_fname&gt;_oblique_DEM.tif with &lt;write_fname&gt;_oblique_DEM_rotation_frame.txt, which records the frame and the lattice to rotate products back onto; target_epsg (or the analysis CRS of the frame centre) is that rotate-back target. false: an ordinary fetch in target_epsg (or UTM) covering the frame's footprint, with no frame file. Read only when fetch_from_linear_feature is true. |

### Source

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `dem_source` | `string` | `cop30_aws` | `param` | - | Which provider to fetch from. | cop30_aws: Copernicus GLO-30 on AWS Open Data. Public and keyless, global 30 m coverage -- start here. opentopography: the OpenTopography REST API, which offers other datasets but needs a free API key (see ot_dem_type). planetary_computer: Microsoft Planetary Computer. Also keyless -- its asset URLs are signed automatically for you -- and gives access to many other STAC collections (see pc_collection). Any other value is an error rather than a silent substitution. |
| `pc_collection` | `string` | `cop-dem-glo-30` | `param` | - | Planetary Computer STAC collection (only used when dem_source is planetary_computer). | Examples: cop-dem-glo-30, cop-dem-glo-90, nasadem, alos-dem, 3dep-seamless. Browse the full catalogue at https://planetarycomputer.microsoft.com/catalog. No account or key is needed: asset URLs are signed automatically. |
| `pc_asset` | `string` | `data` | `param` | - | Asset key within each Planetary Computer STAC item. | Most elevation collections expose the raster as 'data'. Change this only if the collection you are using names it differently. |
| `ot_dem_type` | `string` | `COP30` | `param` | - | OpenTopography dataset identifier (only used when dem_source is opentopography). | One of: COP30, COP90, SRTMGL1, SRTMGL1_E, SRTMGL3, NASADEM, AW3D30, AW3D30_E, GEDTM30, EU_DTM, ANADEM, CA_MRDEM_DTM, CA_MRDEM_DSM, SRTM15Plus, GEBCOIceTopo, GEBCOSubIceTopo, GEDI_L3. Run -h ot_dem_type for a table giving data type, grid spacing, coverage and purpose, or see https://portal.opentopography.org/apidocs/#/Public/getGlobalDem. Anything else is rejected before the request is sent. Every one of these datasets carries a citation requirement; the citation for the dataset you used is printed at the end of a successful run. THE API KEY IS NOT A PARAMETER: set the OT_API_KEY environment variable, or set OT_API_KEY_FILE to a file containing it. Parameter files are routinely committed and every parameter is echoed at the end of a run, so a key must never be written here. |
| `prefetch_source` | `string` | `cop30_aws` | `param` | - | Source for coarse prefetch in bounds_mode: outlets. | Defaults to keyless cop30_aws. Can be cop30_aws, planetary_computer, or opentopography. |
| `coarse_prefetch_demtype` | `string` | `GEBCOIceTopo` | `param` | - | OpenTopography demtype used for coarse nested prefetch. | Used in outlet mode when prefetch_source is opentopography and credentials are available. |
| `allow_auth_fallback_to_keyless` | `bool` | `true` | `param` | - | Allow fallback to a keyless source on OpenTopography auth failure. | When true, auth failures in opentopography can fall back to auth_fallback_source. |
| `auth_fallback_source` | `string` | `cop30_aws` | `param` | - | Fallback source used after OpenTopography auth failure. | Must be keyless: cop30_aws or planetary_computer. |

### Output grid

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `prefetch_grid_spacing` | `float` | `90.0` | `param` | - | Grid spacing (m) for coarse prefetch in bounds_mode: outlets. | Used for the low-resolution basin-discovery stage. |
| `discard_intermediate_rasters` | `bool` | `true` | `param` | - | Delete temporary intermediate rasters from nested outlet stages. | Final outputs and manifests are retained. |
| `write_basin_outlines_fgb` | `bool` | `true` | `param` | - | Write basin outlines as FlatGeobuf outputs in outlet mode. | Applies to retained stages/final outputs. |
| `write_fetch_extents_manifest` | `bool` | `true` | `param` | - | Write manifest of exact fetch extents used in outlet mode. | Includes stage/source/spacing and any bbox contractions. |
| `write_discovered_bounds` | `bool` | `false` | `param` | - | Write discovered outlet-mode bounds for QC. | Used in bounds_mode: outlets. |
| `write_discovered_basin_polygons` | `bool` | `false` | `param` | - | Write discovered basin polygons for QC. | Used in bounds_mode: outlets. |
| `discovered_bounds_fname` | `string` | - | `param` | - | Optional filename override for discovered-bounds output. | Used only when write_discovered_bounds is true. |
| `final_basin_edge_buffer_pixels` | `int` | `3` | `param` | - | Required basin margin from raster edges in final DEMs (pixels). | Must be &gt;= 0. |
| `final_basin_nodata_buffer_pixels` | `int` | `1` | `param` | - | Required basin margin from nodata in final DEMs (pixels). | Must be &gt;= 0. |
| `verify_final_resolution_basins` | `bool` | `true` | `param` | - | Recheck basin extraction on final-resolution DEM. | Guards against coarse-routing mismatch from prefetch. |
| `final_verify_expand_factor` | `float` | `1.25` | `param` | - | Expansion factor for final-stage refetch after failed verification. | Must be &gt; 1.0. |
| `final_verify_max_iters` | `int` | `3` | `param` | - | Maximum final-stage refetch retries after failed verification. | Must be &gt;= 1 when verify_final_resolution_basins is true. |
| `target_epsg` | `int` | `0` | `param` | - | EPSG code of the output projection. | 0 derives the UTM zone from the bounding box centre, which is normally what you want. Set explicitly to pin a projection, for example when fetching several adjacent areas that must share a grid. |
| `grid_spacing` | `float` | `30.0` | `param` | - | Output grid spacing in projection units (metres for UTM). | The source is nominally 30 m. Coarser values resample down; finer values interpolate and add no real detail. |
| `resampling_method` | `string` | `cubic` | `param` | - | GDAL resampling method used when reprojecting. | Cubic is the default because fetched DEM elevations are continuous; it is not a categorical-raster default. Allowed values follow GDAL warp methods: near, bilinear, cubic, cubicspline, lanczos, average, mode, max, min, med, q1, q3. |
| `size_pixels` | `int` | `0` | `param` | - | Crop an exact square of this many pixels from the centre. | 0 keeps the full reprojected extent. Use a positive value when the output dimensions must be deterministic, for example for test fixtures. The bounding box must be large enough to contain the crop. |
| `no_data_value` | `float` | `-9999.0` | `param` | - | NoData value written into the output raster. | Matches the lsdtt3 default. Cells outside the source footprint are set to this rather than zero, which would otherwise read as sea level and behave as a large spurious depression in flow routing. |

### Network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `http_timeout_seconds` | `int` | `60` | `param` | - | Maximum seconds for an individual HTTP request. | Bounds how long an unresponsive server can stall the fetch. |
| `http_connect_timeout_seconds` | `int` | `15` | `param` | - | Maximum seconds to establish a connection. | Keep this short so an offline machine fails promptly. |
| `http_max_retry` | `int` | `3` | `param` | - | Number of retries on a transient HTTP failure. | Cloud storage occasionally returns 5xx responses; retrying absorbs those without failing the run. |
| `cache_dir` | `string` | - | `param` | - | Optional directory for caching fetched rasters. | Empty (the default) means no caching: the fetch always runs. When set, a previously fetched raster covering the same request is reused instead of downloading again. |

<!-- END GENERATED CLI PARAMETERS -->
