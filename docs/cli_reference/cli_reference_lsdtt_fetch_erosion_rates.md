# CLI reference: lsdtt-fetch-erosion-rates

Downloads basin-averaged cosmogenic erosion rates from
[OCTOPUS](https://octopusdata.org) (the CRN Global compilation) for a bounding
box. It writes the sampled basins as an erosion-rate sample CSV (the outlet
latitude and longitude, the rate, its 1-sigma uncertainty, the units and the
reported drainage area), the basin polygons with every OCTOPUS attribute, a
suggested DEM box around the basins, and the attribution.

OCTOPUS data are licensed CC BY 4.0. Every run writes an attribution file naming
OCTOPUS and each original study; cite them in work that uses the rates.

## Usage

```
lsdtt-fetch-erosion-rates <path> <parameter_file>
lsdtt-fetch-erosion-rates <parameter_file>         # path defaults to ./
lsdtt-fetch-erosion-rates -h                       # print help and exit
```

## What it does

1. Asks the OCTOPUS WFS (`http://geoserver.octopusdata.org/geoserver/wfs`,
   keyless) how many basins of `be10-denude:crn_int_basins` intersect the box
   (`resultType=hits`), and stops if there are more than `max_features`.
2. Reads the basins, then the outlets (`be10-denude:crn_int_outlets`), in pages
   of `octopus_page_size` features sorted by `OBSID1`, and joins them on
   `OBSID1`.
3. Nulls the no-data values in every numeric field and prints how many were
   nulled per field. OCTOPUS marks missing values with -9999 or -9999.99, and
   with the same all-nines pattern truncated to the field's width (-99.99,
   -99.999, -9.999, -9.99999).
4. Writes one CSV row per basin whose outlet lies inside the box and that has a
   positive rate under `octopus_rate_scheme`. Basins with no usable rate are
   listed on stdout and left out of the CSV.

The WFS host answers over HTTP only. The data are public and no credential is
sent, but their integrity in transit is not protected. The box is sent as
`bbox=S,W,N,E,urn:ogc:def:crs:EPSG::4326` (latitude first): the server silently
returns no features for other forms. A box query returns every basin whose
polygon intersects the box, so large basins whose outlets lie far outside it
also come back. They are always in the polygon file, but have a CSV row only
with `include_outlets_outside_bbox: true`.

## Rate schemes

All rates are in mm/kyr (= m/Myr, 1e-6 m/yr).

| `octopus_rate_scheme` | Rate | Uncertainty (`erosion_rate_uncertainty_kind`) | `erosion_rate_internal_uncertainty` |
|---|---|---|---|
| `st` (default) | `EBE_MMKYR`: recalculated by OCTOPUS with CAIRN for every study | `EBE_ERR` (`total_1sigma`) | `ERRBE_AMS` × `EBE_MMKYR` / `EBE_GCMYR` |
| `lsdn` | `EBEL` | `EBEL_EXT` (`external_1sigma`) | `EBEL_INT` |
| `lsdn_no_topo_shielding` | `EBELNS` | `EBELNS_EXT` (`external_1sigma`) | `EBELNS_INT` |
| `published` | `BE10EP`: the source paper's own value | `BE10EP_ERR` (`published_1sigma`) | blank |

`octopus_nuclide: al26` uses the Al-26 fields (`EAL_MMKYR` / `EAL_ERR`, `EALL*`,
`EALLNS*`, `AL26EP*`).

## Outputs

| File | Description |
|---|---|
| `<write_fname>_erosion_rates.csv` | One row per sample: `sample_id` (`OBSID1`), `latitude`, `longitude` (the outlet, WGS84), `erosion_rate`, `erosion_rate_uncertainty` (1 sigma, blank when OCTOPUS has none), `erosion_rate_units` (`mm_per_kyr`), `erosion_rate_uncertainty_kind`, `erosion_rate_internal_uncertainty`, `reported_area_km2` (`AREA`), then `rate_scheme`, `rate_field`, `nuclide`, `obsid2`, `study_id`, `author`, `pub_year`, `ref_doi`, `mean_elevation_m` (`ELEV_AVE`) and `outlet_in_bbox`. With no basin in the box it holds the header only. |
| `<write_fname>_octopus_basins.fgb` (or `.gpkg`) | Layer `octopus_basins`: every basin that intersects the box, in the UTM zone of the box centre (or `output_epsg`), with the WFS feature id (`id`), every OCTOPUS attribute (no-data values null), `outlet_in_bbox`, `outlet_longitude` and `outlet_latitude`. |
| `<write_fname>_octopus_dem_bounds.fgb` | Layer `dem_bounds`: one WGS84 rectangle, the envelope of the basins that have their outlet in the box and `AREA` <= `dem_bounds_max_basin_area_km2`, padded by `dem_bounds_pad_degrees`. It can be given as `bounds_fname` to `lsdtt-fetch-raster` and `lsdtt-fetch-lithology`. Not written when no basin qualifies. |
| `<write_fname>_ATTRIBUTION.txt` | The source, layers, access date and database version (`DBVERNO`), the box and rate scheme, the citations, and one line per original study (`STUDYID, AUTHOR, PUBYEAR, REFDOI, n samples`). |

Stdout gives the number of basins and outlets read, the HTTP requests and bytes,
the no-data values nulled per field, the samples written and left out, the rate
range, and the suggested DEM box.

## Attribution

Cite OCTOPUS and each original study listed in `_ATTRIBUTION.txt`:

- Codilean, A. T., et al. (2018). Earth System Science Data, 10, 2123-2139.
  doi:10.5194/essd-10-2123-2018
- Codilean, A. T., et al. (2022). Earth System Science Data, 14, 3695-3713.
  doi:10.5194/essd-14-3695-2022
- Codilean, A. T. and Munack, H. (2025). Geochronology, 7, 113-122.
  doi:10.5194/gchron-7-113-2025
- The dataset (CRN Global v2.3): doi:10.71747/uow-r3gk326m.28216865.v1

For the `st` rates, also cite CAIRN: Mudd, S. M., et al. (2016). Earth Surface
Dynamics, 4, 655-674. doi:10.5194/esurf-4-655-2016.

## Example: West Virginia

```
# wvoct_large_fetch.param
write_fname: wvoct_large
west: -80.5
east: -78.5
south: 38.3
north: 39.7
octopus_rate_scheme: st
```

On 2026-10-19 this returned 28 basins, 23 of them with the outlet in the box
(studies S046 and S186), with `EBE_MMKYR` from 3.74 to 24.08 mm/kyr, in 4 HTTP
requests (1.03 MB). The suggested DEM box covered 21 basins. Counts grow as
OCTOPUS adds studies.

`bounds_fname` also accepts a raster, e.g. an existing DEM GeoTIFF: its
footprint, reprojected to WGS84, is the box.

## Working offline

`octopus_response_dir` reads `crn_int_basins.json` and `crn_int_outlets.json`
(the GeoJSON bodies of a previous fetch) from a directory instead of the WFS, and
`hits_basins.xml` (a `resultType=hits` response) when present. No request is
made. The committed test responses are in `tests/data/octopus_erosion_rates`.

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

### Bounding box

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `west` | `string` | - | `param` | - | Western edge of the bounding box, WGS84 decimal degrees. | Required. Must be less than east. Range -180 to 180. |
| `east` | `string` | - | `param` | - | Eastern edge of the bounding box, WGS84 decimal degrees. | Required. Must be greater than west. Range -180 to 180. |
| `south` | `string` | - | `param` | - | Southern edge of the bounding box, WGS84 decimal degrees. | Required. Must be less than north. Range -90 to 90. |
| `north` | `string` | - | `param` | - | Northern edge of the bounding box, WGS84 decimal degrees. | Required. Must be greater than south. Range -90 to 90. |
| `bounds_fname` | `string` | - | `param` | - | Vector file (or raster) defining the area, instead of west/south/east/north. | FlatGeobuf (.fgb), GeoPackage (.gpkg) or any other OGR-readable format. The layer must contain either one or more polygons (the envelope of all of them is used), or exactly two points marking opposite corners. For two points either diagonal works -- lower-left/upper-right and lower-right/upper-left give the same box, so the orientation is worked out automatically. If the layer is in a projected CRS such as UTM it is reprojected to WGS84 for you, so a box drawn in a GIS can be used directly. A raster with a georeference (e.g. an existing DEM GeoTIFF) works too: its footprint is reprojected to WGS84. Set this or the four edge parameters, not both. |
| `bounds_path` | `string` | - | `param` | - | Optional directory for bounds_fname. | If empty, read_path is used. Ignored when bounds_fname is an absolute path. |
| `bounds_layer` | `string` | - | `param` | - | Layer name within bounds_fname. | If empty, the first layer is used (with a warning when the file has more than one). |

### OCTOPUS

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `octopus_rate_scheme` | `string` | `st` | `param` | - | Which OCTOPUS erosion rate to write: st, lsdn, lsdn_no_topo_shielding or published. | st: EBE_MMKYR / EBE_ERR, recalculated by OCTOPUS with CAIRN (one calculator for every study; uncertainty is the total 1 sigma). lsdn: EBEL / EBEL_EXT (external 1 sigma; internal EBEL_INT). lsdn_no_topo_shielding: EBELNS / EBELNS_EXT. published: BE10EP / BE10EP_ERR, the source papers' own values. All in mm/kyr. |
| `octopus_nuclide` | `string` | `be10` | `param` | - | Nuclide of the rate: be10 or al26. | al26 uses EAL_MMKYR / EAL_ERR, EALL*, EALLNS* or AL26EP*; most basins have no Al-26 rate and are then left out of the CSV. |
| `octopus_page_size` | `int` | `100` | `param` | - | Features per WFS GetFeature page. | Basins are counted first (resultType=hits) and then read in pages of this size. |
| `max_features` | `int` | `2000` | `param` | - | Refuse a box for which OCTOPUS returns more basins than this. | A guard against an accidental continental request (the full CRN Global layer is about 5400 basins and several hundred MB). Raise it deliberately for large areas. |
| `include_outlets_outside_bbox` | `bool` | `false` | `param` | - | Also write CSV rows for basins whose outlet lies outside the box. | The box query returns every basin whose polygon intersects it, so large basins with outlets far away come back too. Their outlets are fetched by OBSID1 (one more request). The polygon file always holds every intersecting basin. |
| `dem_bounds_pad_degrees` | `float` | `0.02` | `param` | - | Padding of the suggested DEM box around the basins, in degrees. | The DEM box is the envelope of the basins with an outlet in the box, so a DEM fetched with it contains each basin whole. |
| `dem_bounds_max_basin_area_km2` | `float` | `1000` | `param` | - | Basins larger than this (OCTOPUS AREA, km2) do not enlarge the suggested DEM box. | - |
| `request_delay_milliseconds` | `int` | `200` | `param` | - | Pause after every request to the OCTOPUS WFS. | OCTOPUS publishes no rate limit; this keeps a fetch polite. |
| `octopus_response_dir` | `string` | - | `param` | - | Read the WFS responses from a local directory instead of downloading them. | The directory holds crn_int_basins.json and crn_int_outlets.json (the GeoJSON bodies of a previous fetch) and optionally hits_basins.xml (a resultType=hits response). No HTTP request is made. |
| `output_epsg` | `int` | `0` | `param` | - | EPSG code of the basin polygons. 0 (default) = the WGS 84 UTM zone of the box centre. | Web Mercator (EPSG:3857) is refused by lsdtt because it distorts distances. The CSV and the DEM box are always in WGS84 degrees. |

### Network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `http_timeout_seconds` | `int` | `60` | `param` | - | Maximum seconds for an individual HTTP request. | Bounds how long an unresponsive server can stall the fetch. |
| `http_connect_timeout_seconds` | `int` | `15` | `param` | - | Maximum seconds to establish a connection. | Keep this short so an offline machine fails promptly. |
| `http_max_retry` | `int` | `3` | `param` | - | Number of retries on a transient HTTP failure. | Cloud storage occasionally returns 5xx responses; retrying absorbs those without failing the run. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). CSV is available only through selectors for point/table outputs. Per-output overrides exist only where a container genuinely differs. Raster outputs are governed separately by write_output_format. |

<!-- END GENERATED CLI PARAMETERS -->
