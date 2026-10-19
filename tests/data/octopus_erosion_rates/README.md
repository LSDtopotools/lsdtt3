# OCTOPUS erosion-rate fixtures

Responses of the OCTOPUS WFS (`http://geoserver.octopusdata.org/geoserver/wfs`,
keyless, WFS 2.0.0) in the layout `lsdtt-fetch-erosion-rates` reads with
`octopus_response_dir`: `<dir>/crn_int_basins.json`, `<dir>/crn_int_outlets.json`
and, optionally, `<dir>/hits_basins.xml`. Used by `unit_fetch_erosion_rates` and
the offline `smoke_` / `failure_fetch_erosion_rates_*` tests.

| Directory | Contents |
|---|---|
| `wv_large_box/` | The box west -80.5, east -78.5, south 38.3, north 39.7 (West Virginia / Virginia / Maryland). 28 basins, 23 outlets. |
| `wv_small_box/` | The box west -79.80, east -79.05, south 38.85, north 39.15 (the Allegheny Front "wvbig" box). 5 basins, 1 outlet (S186WTS046). |
| `empty_box/` | A hand-written empty FeatureCollection (`numberMatched` 0). |
| `exception_report/` | `crn_int_basins.json` is a hand-written OWS 1.1 ExceptionReport (an unknown feature type), the XML GeoServer returns for a bad request. |

## Provenance

Fetched 2026-10-18 (basins and outlets) and 2026-10-19 (`hits_basins.xml`) with

```
service=WFS&version=2.0.0&request=GetFeature
  &typeNames=be10-denude:crn_int_basins      (or be10-denude:crn_int_outlets)
  &outputFormat=application/json&srsName=EPSG:4326
  &bbox=<south>,<west>,<north>,<east>,urn:ogc:def:crs:EPSG::4326
```

and, for `hits_basins.xml`, the same basins request with `resultType=hits` and no
`outputFormat`. The bbox must be latitude first; other forms return no
features.

- `crn_int_outlets.json` and `hits_basins.xml` are the responses verbatim.
- `crn_int_basins.json` keeps every feature, its `id`, and all 93 attributes with
  their values unchanged (re-serialised as compact JSON), but **each geometry is
  replaced by the bounding rectangle of the basin** (as a one-ring MultiPolygon)
  to keep the files small. They are not the real basin polygons.

The attributes include OCTOPUS's no-data values (-9999, -9999.99 and the
truncated all-nines forms -99.99, -99.999, -9.999, -9.99999), which the tests
expect to be nulled. All basins are Be-10 only, so every Al-26 field is no-data.

## Licence and citation

OCTOPUS data are licensed CC BY 4.0 (https://octopusdata.org). Cite:

- Codilean, A. T., et al. (2018). Earth System Science Data, 10, 2123-2139.
  doi:10.5194/essd-10-2123-2018
- Codilean, A. T., et al. (2022). Earth System Science Data, 14, 3695-3713.
  doi:10.5194/essd-14-3695-2022
- Codilean, A. T. and Munack, H. (2025). Geochronology, 7, 113-122.
  doi:10.5194/gchron-7-113-2025
- Dataset (CRN Global v2.3): doi:10.71747/uow-r3gk326m.28216865.v1

and the original studies of these samples:

- S186: Portenga (2019), doi:10.1130/b31840.1 (17 basins in
  `wv_large_box`, 5 in `wv_small_box`);
- S046: Duxburry (2009), REFID `phdthesis`, no DOI in OCTOPUS (11 basins in
  `wv_large_box`).
