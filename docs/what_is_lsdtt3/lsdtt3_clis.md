# The lsdtt3 command-line tools

lsdtt3 ships **40 command-line programs**. This page lists all of them, grouped
by what they are for, so you can find the one you need without reading forty
reference pages.

Every tool is driven by a parameter file rather than command-line flags — if
you have not met one yet, start with [Your first analysis](first_analysis.md).
Every tool also documents itself: `<tool> --help` prints its full parameter
list, and the same content is on this site under
[Command line reference](../cli_reference/cli_reference_lsdtt_raster_preprocessing.md).

!!! tip "If you only read one line of this page"

    Run [`lsdtt-raster-preprocessing`](../cli_reference/cli_reference_lsdtt_raster_preprocessing.md)
    on your DEM before anything else. It enforces the projection, fills or
    breaches pits, and gives you a hillshade to check the data loaded the way
    you expected.

---

## Getting and preparing data

| Tool | What it does |
|---|---|
| `lsdtt-fetch-raster` | Fetch a DEM for a bounding box and write it as a raster. Keyless COP30 and Planetary Computer sources; OpenTopography needs a key. |
| `lsdtt-fetch-erosion-rates` | Fetch basin-averaged cosmogenic erosion rates from OCTOPUS (keyless) for a bounding box: an erosion-rate sample CSV, the basin polygons, a suggested DEM box and the CC BY 4.0 attribution. |
| `lsdtt-raster-preprocessing` | The preparation pipeline: projection enforcement and reprojection, masking, nodata filling, breaching, pit filling, hillshade. |
| `lsdtt-raster-calculator` | Two-raster algebra. |
| `lsdtt-stream-burning` | Conditions a DEM using a supplied vector channel network. |

## Surfaces, slopes and curvature

| Tool | What it does |
|---|---|
| `lsdtt-surface-metrics` | Slope, aspect, curvature and related derived rasters from polynomial surface fitting, plus the SMAP shape classification. |
| `lsdtt-spectral-analysis` | Perron-style spectral analysis of a DEM. Requires FFTW3. |
| `lsdtt-feature-finder` | Experimental terrain-feature detection for LiDAR-derived DEMs. |
| `lsdtt-mountain-front` | Identifies mountain masses via adaptive relief thresholding. |
| `lsdtt-linear-stats` | Statistics of a raster sampled along linear vector features. |
| `lsdtt-swath-profile` | Swath profile along an ordered polyline baseline. |
| `lsdtt-cross-sections` | Channel network from supplied source points, with extracted cross sections. |

## Flow routing and drainage

| Tool | What it does |
|---|---|
| `lsdtt-drainage-analysis` | Static topographic indices: drainage area, TWI, STI, and the Cavalli connectivity index. No time dimension. |
| `lsdtt-channel-extraction` | Extract a channel network from a DEM. |
| `lsdtt-channel-network-analysis` | Build a channel and junction network, then run one or more network analyses on it. |
| `lsdtt-internal-drainage` | Detect, review, route, tag and analyse internally drained basins — moulins, playas, quarries. |
| `lsdtt-outlet-matching` | Match reference-DEM outlets to target-DEM cells by comparing D8 contributing areas. |

## Basins

| Tool | What it does |
|---|---|
| `lsdtt-basin-extraction` | Extract a channel network and select basins from it. |
| `lsdtt-basin-stats` | Basins from a CSV of outlet points, with per-basin statistics. |
| `lsdtt-basin-surface-metrics` | Basins from outlet coordinates, with multiscale surface metrics per basin. |

## Channel profiles and chi

| Tool | What it does |
|---|---|
| `lsdtt-chi-analysis` | End-to-end chi workflow. Requires Eigen3. |
| `lsdtt-sediment-routing` | 1-D reach-scale bedload sediment routing on a D8 channel network. Also where the GIUH lives. |

## Hillslopes and divides

| Tool | What it does |
|---|---|
| `lsdtt-ridgeline-analysis` | Drainage divides and their elevation profiles. |
| `lsdtt-divide-migration` | Divide-migration metrics from supplied outlets. |
| `lsdtt-hillslope-channel-coupling` | Hilltop curvature, hillslope length and relief, and the non-dimensional E\*/R\* pair, traced from ridgeline to channel. |

## Floodplains, terraces and valleys

| Tool | What it does |
|---|---|
| `lsdtt-flat-surface-finder` | Identifies broadly flat surfaces — floodplains, alluvial fans, fluvial and marine terraces. |
| `lsdtt-fluvial-surfaces` | Stratifies flat ground into fluvial levels: floodplain, terrace treads, fan surfaces. |
| `lsdtt-valley-centreline` | Per-channel valley identification and a medial-axis valley centreline. |
| `lsdtt-water-centreline` | Centreline from each water-area polygon, for wide river reaches. |
| `lsdtt-water-line-repair` | Reconnects stream polylines truncated where roads cross. |

## Hydrology and hydraulics

| Tool | What it does |
|---|---|
| `lsdtt-distributed-hydrology` | A raster TOPMODEL: time-evolving, rainfall-driven runoff simulation with distributed stores. |
| `lsdtt-graphflood` | Stationary graph-based hydraulic routing. Solves for steady-state flow depths given a rainfall-rate field. |

See [Hydrological modelling](../design_and_algorithms/hydrological_modelling.md)
for which of these answers which question — they are complementary, not
alternatives, and neither is a flood-timing model.

## Coastal and tidal

| Tool | What it does |
|---|---|
| `lsdtt-marsh-topography` | Marsh and tidal-channel analysis workflow. |

## Cosmogenic nuclides

| Tool | What it does |
|---|---|
| `lsdtt-cosmo-tool` | Erosion-rate inversion from measured nuclide concentrations at sample points. |
| `lsdtt-cosmo-prediction` | Predicts nuclide concentrations across a DEM: steady state, step change, or accumulated downstream. |
| `lsdtt-cosmo-column-tool` | 1-D column simulation. No DEM required. |
| `lsdtt-cosmo-landslides` | Triggered landslides with an imposed size distribution, excavated and written out for the cosmogenic and sediment-routing tools. |

The physics, the scaling schemes and the CRONUS/OCTOPUS comparisons are in
[Cosmogenic nuclides](../design_and_algorithms/cosmogenic_nuclides.md).

## Landscape evolution and synthetic landscapes

| Tool | What it does |
|---|---|
| `lsdtt-lem` | Landscape evolution model descended from MuddPILE and DONKmodel. Fluvial incision plus a choice of hillslope solvers. |
| `lsdtt-synthetic-surfaces` | Synthetic rasters for testing analysis routines. |
| `lsdtt-synthetic-networks` | Synthetic Optimal Channel Networks with tunable branching. |

[Choosing a synthetic landscape tool](../design_and_algorithms/synthetic_cli_cluster.md)
explains which of the three to reach for.

---

## Other ways to drive lsdtt3

The command-line tools are not the only interface:

- **[Python bindings](../python/index.md)** expose much of the library
  directly, which is often easier for exploratory work and for stitching steps
  together.
- **[The QGIS plugin](../qgis.md)** runs a subset of the tools from inside QGIS.
- **[Snakemake workflows](workflows.md)** chain tools into reproducible
  multi-step pipelines.
