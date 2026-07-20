# `lsdtt-fetch-raster` — CLI reference

Downloads a DEM for an area and writes it as an analysis-ready projected
GeoTIFF. No Python environment, and by default no account or API key.

Design notes and rationale: `docs/fetch_raster_cli_plan.md`.

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
| `ot_dem_type` | string | `COP30` | OpenTopography dataset (that source only) |
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

Useful `ot_dem_type` values: `COP30`, `COP90`, `SRTMGL1`, `SRTMGL3`, `NASADEM`,
`AW3D30`, `SRTM15Plus`, `GEBCOIceTopo`.

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

**`the OpenTopography API returned an error: ... API Key: ***`**
The key was rejected. The `***` is the masking working as intended — your key
is not being printed.

**`could not open '...' as a vector dataset`**
`bounds_fname` is not a readable vector file. Check the path (it resolves
against `bounds_path`, then `read_path`) and that the format is one GDAL builds
support.

---

## See also

- `docs/fetch_raster_cli_plan.md` — design rationale, credential rules, findings
- `docs/cli_reference_lsdtt_raster_preprocessing.md` — the usual next step
- `python_helpers/download_dem_tiles.py` — Python helper covering other providers
