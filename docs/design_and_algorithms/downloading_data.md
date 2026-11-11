# Downloading data

Four lsdtt3 tools download the data an analysis starts from:

| Tool | Downloads | Source |
|---|---|---|
| `lsdtt-fetch-raster` | Digital elevation models | Copernicus, national and polar DEMs, OpenTopography, Planetary Computer |
| `lsdtt-fetch-lithology` | Geological map units | Macrostrat |
| `lsdtt-fetch-river-lines` | Mapped river centrelines | OS Open Rivers, or a line file you already have |
| `lsdtt-fetch-erosion-rates` | Basin-averaged cosmogenic erosion rates | OCTOPUS |

All four take the area the same way, and write data ready for the analysis
tools: projected in metres (by default the WGS 84 UTM zone of the area's
centre), in the formats the other tools read. This page describes what each
tool fetches and how they fit together. Every parameter is in the tools'
command-line reference pages, linked in each section.

**The area.** Give it as a box in WGS 84 decimal degrees:

```
west:  10.34
south: 46.46
east:  10.46
north: 46.54
```

or as `bounds_fname`: a vector file whose polygons' envelope is the area (a
polygon you drew in a GIS, for example), two points marking opposite corners,
or a raster such as an existing DEM, whose footprint is the area.

**Attribution.** The data belong to their publishers, and most licences
require you to cite them. Each tool prints or writes the citation for what it
fetched: `lsdtt-fetch-raster` prints it for OpenTopography datasets, and the
other three write a `<write_fname>_ATTRIBUTION.txt` next to their outputs.

---

## Elevation: `lsdtt-fetch-raster`

Downloads a DEM for the area and writes it as a projected GeoTIFF on square
cells of `grid_spacing` metres (30 by default, or the source's native
resolution for the high-resolution sources), resampled with
`resampling_method` (cubic by default). Pick the provider with `dem_source`.
Full details: [lsdtt-fetch-raster reference](../cli_reference/cli_reference_lsdtt_fetch_raster.md).

**Global, no account needed:**

| `dem_source` | Data |
|---|---|
| `cop30_aws` (default) | Copernicus GLO-30, the global 30 m DEM, from AWS Open Data. Start here. |
| `planetary_computer` | Microsoft Planetary Computer: Copernicus GLO-30 and the other elevation collections it hosts, such as NASADEM, ALOS and USGS 3DEP (`pc_collection`). |

**Global, with a free OpenTopography API key:** `dem_source: opentopography`,
with `ot_dem_type` choosing the dataset: Copernicus 30 and 90 m, SRTM (1 and 3
arc-second, and ellipsoidal heights), NASADEM, ALOS World 3D, the global
bare-earth ensemble DTM (GEDTM30), the European and Canadian DTMs, ANADEM for
South America, and the 500 m land-and-sea products (SRTM15+, GEBCO). The key
is read from the `OT_API_KEY` environment variable or a file named by
`OT_API_KEY_FILE`, never from the parameter file.

**National and polar high-resolution DEMs, no account needed:**

| `dem_source` | Product | Resolution (m) | Coverage |
|---|---|---|---|
| `us_3dep` | USGS 3DEP | 10, 30 | United States |
| `arctic_arcticdem` | ArcticDEM (surface model) | 2, 10, 32 | Land north of 60° N |
| `antarctic_rema` | REMA (surface model) | 2, 10, 32 | Antarctica |
| `no_dtm` | Kartverket height model | 1 | Mainland Norway |
| `nl_ahn` | AHN | 0.5 | The Netherlands |
| `es_mdt` | IGN MDT05 | 5 | Mainland Spain and the Balearics |
| `pl_nmt` | GUGiK NMT | 1 | Poland |
| `ca_hrdem` | NRCan HRDEM | 1, 2 | Canada, where lidar exists |
| `ch_swissalti3d` | swissALTI3D | 0.5, 2 | Switzerland and Liechtenstein |
| `nz_lidar` | LINZ lidar | 1 | New Zealand |
| `gb_sct_lidar` | Lidar for Scotland | 0.5 or 1 | Parts of Scotland |
| `jp_gsi` | GSI elevation tiles | 1, 5, 10 | Japan |

`lsdtt-fetch-raster -h dem_source` prints the full list, with each source's
datum and coverage, without a network connection.

**Surface or terrain model.** Copernicus, ALOS, ArcticDEM and REMA are
*surface* models: over forest and towns they measure the top of the canopy and
the roofs, not the ground. Hillslope curvature and roughness there describe the
vegetation. Bare-earth terrain models (GEDTM30, the national lidar DTMs) avoid
this where they cover your area; `surface_model_choice` chooses between the
two where a source offers both.

**Other ways to choose the area:**

- `bounds_mode: outlets` starts from basin outlet points and finds the area
  itself: it fetches a coarse DEM, delineates the basins draining to the
  outlets, then fetches the final DEM around them, split into separate DEMs
  where the basins are far apart;
- `fetch_from_linear_feature` fetches a strip along a line, such as a fault or
  a river, and can write it straight into an oblique frame with the feature
  along the x axis; see
  [Analysing along a feature](oblique_frames.md).

---

## Geology: `lsdtt-fetch-lithology`

Downloads geological map units from [Macrostrat](https://macrostrat.org), a
global compilation of geological maps, and writes them as a polygon layer
(`<write_fname>_macrostrat_units.fgb`). Each unit carries Macrostrat's
lithology text and `lith_primary`, its first lithology, which the lithology
tools read by default.

It reads Macrostrat's vector map tiles, which blend source maps of different
scales and serve more detailed ones at higher zoom (`macrostrat_tile_zoom`).
It then looks up each unit's linked lithologies and source maps in the
Macrostrat API; if the API cannot be reached the run continues with the
tiles' text and a warning. It also writes a vocabulary of the lithologies
found and a template crosswalk table, for grouping many map units into the
few classes an analysis needs.

The output is the `lithology_fname` input of `lsdtt-lithology-erodibility`
(an erodibility raster for the landscape evolution model) and of
`lsdtt-chi-analysis` (channel steepness by rock type; see
[Lithology and channel steepness statistics](lithology_channel_statistics.md)).
Macrostrat data are CC BY 4.0.
Full details: [lsdtt-fetch-lithology reference](../cli_reference/cli_reference_lsdtt_fetch_lithology.md).

---

## Rivers: `lsdtt-fetch-river-lines`

**Experimental.** Reads mapped river centrelines and prepares them for
`lsdtt-stream-burning`, which burns them into a DEM (see
[Burning mapped rivers into a DEM](burning_river_lines.md)). Choose the source
with `river_line_source`:

- `gb_os_open_rivers`: Ordnance Survey Open Rivers, for Great Britain. The
  national file is downloaded once and reused from a cache (`river_cache_dir`);
- `local_file`: any line file you already have, such as OS MasterMap Water
  Network or an OpenStreetMap extract. Fields naming each line's type, flow
  direction and end nodes can be mapped in, where the file has them.

The lines are clipped to the area plus a buffer, split into single lines,
joined at junctions within `node_snap_tolerance_metres`, oriented downstream
and written in the DEM's CRS (`<write_fname>_river_lines.fgb`, with the nodes
and a topology table). The tool reports the network's problems, such as
loops, bifurcations and lines pointing against the network, but repairs only
line direction, and optionally gaps (`gap_bridge_algorithm_choice`).
Full details: [lsdtt-fetch-river-lines reference](../cli_reference/cli_reference_lsdtt_fetch_river_lines.md).

---

## Erosion rates: `lsdtt-fetch-erosion-rates`

Downloads basin-averaged cosmogenic erosion rates (<sup>10</sup>Be by default,
or <sup>26</sup>Al with `octopus_nuclide: al26`) from
[OCTOPUS](https://octopusdata.org), the global compilation of published
detrital cosmogenic samples, for every basin that intersects the box. It
writes:

- `<write_fname>_erosion_rates.csv`: one row per basin with its outlet in the
  box, with the outlet's latitude and longitude, the rate, its 1-sigma
  uncertainty, the units and the reported drainage area;
- the basin polygons, with every OCTOPUS attribute;
- `<write_fname>_octopus_dem_bounds.fgb`: a box around the basins, to fetch a
  DEM that contains them;
- the attribution, with one line per original study.

`octopus_rate_scheme` chooses which rate is used: `st` (the default) is the
rate OCTOPUS recalculated with one consistent method (CAIRN) for every study;
`lsdn` and `lsdn_no_topo_shielding` use the LSDn production scaling;
`published` is each source paper's own value, calculated with that study's
own methods.

The erosion-rate CSV is the `erosion_rate_fname` input of `lsdtt-chi-analysis`,
which calibrates the erodibility *K* from channel steepness and the measured
rates. OCTOPUS data are CC BY 4.0.
Full details: [lsdtt-fetch-erosion-rates reference](../cli_reference/cli_reference_lsdtt_fetch_erosion_rates.md).

---

## Putting them together

A typical start for an analysis of channel steepness, rock type and erosion
rate:

1. `lsdtt-fetch-erosion-rates` for a box, to find the sampled basins;
2. `lsdtt-fetch-raster` with `bounds_fname: <write_fname>_octopus_dem_bounds.fgb`,
   so the DEM contains whole basins (an incomplete basin is dropped from the
   chi analysis);
3. `lsdtt-fetch-lithology` with the DEM as `bounds_fname`, so the map units
   cover the same area;
4. `lsdtt-chi-analysis` with `lithology_fname` and `erosion_rate_fname` set.

Each tool runs on its own, so any of these steps can use data you already
have instead.

**Working offline.** The DEM tool reuses earlier downloads if you set
`cache_dir`; Open Rivers is cached after its first download; the lithology and
erosion-rate tools have their own offline routes, described in their reference
pages.
