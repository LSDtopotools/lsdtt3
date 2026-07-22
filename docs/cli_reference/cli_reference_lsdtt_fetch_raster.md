# `lsdtt-fetch-raster` — CLI reference

Downloads a DEM for an area and writes it as an analysis-ready projected
GeoTIFF. No Python environment, and by default no account or API key.

Design notes and rationale: `docs/testing_and_benchmarking.md`.

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
than a silent zero — a zeroed edge would quietly fetch the Gulf of Guinea.

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

| `ot_dem_type` | Dataset | Res. | Coverage | Notes |
|---|---|---|---|---|
| `COP30` | Copernicus GLO-30 | 30 m | global | TanDEM-X **DSM**, edited water bodies and coastlines. The default, and the same data `cop30_aws` serves keylessly. |
| `COP90` | Copernicus GLO-90 | 90 m | global | The 90 m Copernicus DSM. Use it when the area exceeds the 30 m request limit. |
| `SRTMGL1` | SRTM GL1 | 30 m | 60 N–56 S | 1 arc-second SRTM, EGM96 heights. The long-standing comparison baseline. |
| `SRTMGL1_E` | SRTM GL1 Ellipsoidal | 30 m | 60 N–56 S | As `SRTMGL1` but WGS84 ellipsoidal heights — for co-registering with GNSS or lidar. |
| `SRTMGL3` | SRTM GL3 | 90 m | 60 N–56 S | 3 arc-second SRTM. Larger request limit than GL1. |
| `NASADEM` | NASADEM | 30 m | 60 N–56 S | SRTM reprocessed: better geolocation, voids filled from ASTER GDEM and ICESat. Generally preferable to raw `SRTMGL1`. |
| `AW3D30` | ALOS World 3D | 30 m | global | JAXA optical-stereo DSM. An *independent* global 30 m surface — useful as a cross-check on the radar products. |
| `AW3D30_E` | ALOS World 3D Ellipsoidal | 30 m | global | `AW3D30` with ellipsoidal rather than geoid heights. |
| `GEDTM30` | Global Ensemble DTM | 30 m | global | Machine-learning ensemble **DTM** — bare earth, vegetation and buildings removed. |
| `EU_DTM` | Continental Europe DTM | 30 m | Europe | Bare-earth DTM from national lidar and photogrammetric surveys. |
| `ANADEM` | ANADEM | 30 m | South America | Bare-earth DTM; vegetation bias removed from Copernicus GLO-30. |
| `CA_MRDEM_DTM` | Canada MRDEM (terrain) | 30 m | Canada | Natural Resources Canada medium-resolution bare-earth model. |
| `CA_MRDEM_DSM` | Canada MRDEM (surface) | 30 m | Canada | The surface-model companion to `CA_MRDEM_DTM`. |
| `SRTM15Plus` | SRTM15+ V2.1 | 500 m | global, land and sea | Merged topography and bathymetry. Use it when the domain continues offshore. |
| `GEBCOIceTopo` | GEBCO 2023 (ice surface) | 500 m | global, land and sea | Bathymetry and topography, with the *ice-sheet surface* over Greenland and Antarctica. |
| `GEBCOSubIceTopo` | GEBCO 2023 (bed under ice) | 500 m | global, land and sea | As above, but the ice sheets are replaced by the underlying bed. |
| `GEDI_L3` | GEDI L3 gridded metrics | 1000 m | 52 N–52 S | Gridded GEDI lidar land-surface metrics. A metrics product, **not** a DEM for flow routing. |

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
| `resampling_method` | string | `cubic` | GDAL warp method: `near`, `bilinear`, `cubic`, `cubicspline`, `lanczos`, `average`, `mode`, `max`, `min`, `med`, `q1`, `q3` |
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

- `docs/testing_and_benchmarking.md` — design rationale, credential rules, findings
- `docs/cli_reference/cli_reference_lsdtt_raster_preprocessing.md` — the usual next step
- `python_helpers/download_dem_tiles.py` — Python helper covering other providers
- <https://portal.opentopography.org/apidocs/#/Public/getGlobalDem> — the
  upstream `demtype` list this tool's catalogue mirrors
- <https://planetarycomputer.microsoft.com/catalog> — collections for
  `pc_collection`
