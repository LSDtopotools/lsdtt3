# Overview

`lsdtt3` is a redesign of [LSDTopoTools2](https://github.com/LSDtopotools/LSDTopoTools2) that retains many of the same scientific algorithms while rebuilding the codebase for modern C++17, better memory efficiency, and easier ongoing development. It also includes a number of new analyses and some replacement algorithms where there is a substantial improvement in performance or accuracy. `lsdtt3` also includes testing suites, as well as a QGIS plugin. 

---

## What it does

`lsdtt3` provides a collection of **45 command-line programs** for analysing digital elevation models (DEMs). Each program reads a plain-text parameter file and writes raster or vector outputs. Much of the same functionality is also available through **Python bindings**, and a subset through a **QGIS plugin**. The tools cover:

- **DEM preparation** — reprojection, depression filling, breaching, nodata handling.
- **Surface morphology** — curvature, roughness, hillshade, local relief, cloth-drape upper envelope, Frangi vessel filter.
- **Drainage and flow** — D8, D-infinity, and multiple-flow-direction routing; topographic wetness index.
- **Channel networks** — extraction from a contributing-area threshold, junction and link construction, Strahler ordering, export to geospatial formats.
- **Basin analysis** — delineation, statistics, chi analysis, k_sn and m/n estimation.
- **Cosmogenic nuclides** — basin-averaged erosion rates, concentration prediction, transient 1-D column modelling.
- **Landscape evolution modelling** — stream-power incision, hillslope diffusion, tectonic uplift.

And many other features. 

---

## Key improvements over LSDTopoTools2

| Feature | LSDTopoTools2 | lsdtt3 |
|---|---|---|
| Raster I/O | BIL/ASC | Any GDAL format (GeoTIFF, COG, VRT, BIL, ASC, …) |
| Vector output | CSV only | GeoPackage, FlatGeobuf, CSV with correct CRS |
| Projection handling | Manual | Automatic UTM / polar stereographic / Lambert LAEA |
| Pixel indexing | 32-bit | 64-bit — handles DEMs > 40 000 × 40 000 px |
| Flow routing | D8 | D8, D-infinity, MFD |
| Memory | Higher | Reduced — see [Memory use](#memory-use) below |
| Parallelism | Limited OpenMP | OpenMP across most bottleneck routines |
| Output format | ENVI bil | Cloud-Optimised GeoTIFF (COG) by default, GeoTIFF on small rasters |
| Testing | No test suite | Multiple test suites to compare against alternative analyses (e.g., TopoToolbox3, RichDEM, Cronus Calculator) 

---

## Memory use

Memory, not speed, is usually what stops a large DEM being processed at all. The
object that dominates is the flow-routing structure both codebases build before
almost any drainage analysis: `lsdtt_flowinfo` in lsdtt3, `LSDFlowInfo` in
LSDTopoTools2. Neither stores the elevation raster itself — this is the cost of
the topology on top of the DEM.

For a **10 000 × 10 000 float32 DEM (10^8 cells, all valid)**:

| | LSDTopoTools2 `LSDFlowInfo` | lsdtt3 `lsdtt_flowinfo` |
|---|---|---|
| Resident after construction | ~52 bytes/cell → **5.2 GB** | 43 bytes/cell → **4.3 GB** |
| …allowing for un-reserved `push_back` growth | ~56 bytes/cell → ~5.6 GB | unchanged (arrays are sized once) |
| Without the O(1) upslope index | not separable | 19 bytes/cell → **1.9 GB** |
| Peak during construction | ~76 bytes/cell → ~7.6 GB | lower; transients are scoped, not held to end of function |

**About 20 % less memory for the same D8 topology, and ~65 % less if you do not
need the upslope index.** The gap comes from storing less, not from cleverer
containers — LSDTopoTools2's layout is already flat, contiguous arrays with no
pointer-chasing:

- lsdtt3 keeps **8 raster-scaled arrays; LSDTopoTools2 keeps 13.** Row and
  column indices (8 bytes/cell) are derived from the flat index rather than
  stored, and several cached traversal products are not retained.
- lsdtt3 uses **1-byte types where the values are small** — flow direction,
  donor count and the valid mask. LSDTopoTools2 stores flow direction
  (values −1…7) and flow-length code (values 0…2) as 4-byte integers.
- Construction transients are scoped so they are released before the object is
  returned, rather than staying alive until the constructor exits.

**Where lsdtt3 deliberately pays more:** every node index is 64-bit, against
32-bit in LSDTopoTools2. That doubles the cost of the index arrays, and lsdtt3
is still smaller overall. It is also what makes the next point possible.

### Very large DEMs

LSDTopoTools2's 32-bit indexing caps it at about 2.1 × 10^9 cells, and its
underlying array class computes the allocation size in a 32-bit integer, so
allocation overflows above roughly **46 000 × 46 000**. A 40 000 × 40 000 DEM is
already near that ceiling.

lsdtt3 is 64-bit throughout, so there is no comparable structural limit — the
constraint is the RAM in the machine, not the arithmetic. Combined with the
lower per-cell cost, that is what makes very large lidar tiles and national
mosaics tractable. If memory is still the binding constraint, the upslope index
can be switched off to trade O(1) upslope-set queries for a 2.3× smaller flow
object.

> **How these numbers were obtained.** The lsdtt3 figures are **measured** —
> resident-set growth across construction, 43.0 bytes/cell at both 1.6 × 10⁷ and
> 3.6 × 10⁷ cells, matching the analytic model exactly and linear in cell count,
> extrapolated to 10⁸. The LSDTopoTools2 figures are **derived from its source**
> (member declarations and constructor allocations), not measured on a running
> job, so they are a lower bound on real usage. Both exclude the DEM itself
> (4 bytes/cell) and any per-analysis output rasters.

---

## How it is structured

lsdtt3 is organised as a set of free functions grouped by task, compiled into a shared library and linked into individual command-line executables. Because each analysis file is independent, editing one routine does not trigger a full rebuild.

The build system is CMake. Optional dependencies (FFTW3, Eigen3, NetCDF) are detected automatically; tools that need a missing library are skipped rather than failing the build.

See [Architecture](../developer_notes/lsdtt3_architecture.md) for the full layer design and naming conventions.

---

## Who it is for

lsdtt3 is aimed at geomorphologists, physical geographers, hydrologists and earth scientists who need to:

- process large or many DEMs in batch workflows,
- run reproducible analyses from parameter files,
- integrate topographic analysis into scripted pipelines.

[The lsdtt3 command-line tools](lsdtt3_clis.md) lists all 45 and what each is
for. For reproducible multi-stage pipelines with Snakemake + Pixi, see
[Workflows](workflows.md).

### Three ways to drive it

The command line is the primary interface, but not the only one:

- **Command line.** Every tool, every parameter. Driven by a parameter file
  — see [Your first analysis](first_analysis.md).
- **[Python bindings](../python/index.md).** Much of the library is exposed
  directly to Python, which suits exploratory work and stitching steps
  together without writing parameter files. There are guides for loading and
  saving DEMs, fetching DEMs, flow and basins, chi and k<sub>sn</sub>,
  channel-head extraction, cosmogenics and GraphFlood, plus a full
  [API reference](../python_api/index.md). Note that the Python package is for
  **driving** lsdtt3; the numerical work stays in the C++ library.
- **[QGIS plugin](../qgis.md).** A reduced set of core analyses, run from
  inside QGIS.

It is not a GIS — it has no interactive GUI. Inputs and outputs are files. However we have provided a lightweight QGIS plugin (which is under development). The QGIS plugin has a reduced feature set to the full `lsdtt3` command line interface collection. But for many users, the "core" analyses tools that are in the QGIS plugin should be sufficient for everyday use. 

The visualisation environment for the developers is the sister package `lsdviztools3`, but files can be read and visualised in any platform that takes geospatial data, for example QGIS, ArcMap, Python, or R. `lsdviztools3` is for **visualisation**: computation belongs in lsdtt3 itself.

---

## Compatibility

lsdtt3 is **not backward-compatible** with LSDTopoTools2. Parameter files, output file names, and some algorithm defaults differ. If you have existing LSDTopoTools2 workflows, consult the [command-line reference](../cli_reference/cli_reference_lsdtt_channel_extraction.md) for the updated parameter names.

lsdtt3 is best run on Linux, and pre-built Linux packages are the main way to
install it. They are self-contained: no compiler, no library installation, no
`sudo`.

- **Windows** — run the Linux package inside WSL2. This is the recommended
  route and gives the full tool set; in our testing it is also faster than a
  native Windows build on the same machine. A native Windows package is
  provided as well, but it is a reduced build intended as the QGIS plugin
  backend. 
- **macOS** — no macOS package yet. Run the Linux package inside a Multipass
  Ubuntu VM.
- **Google Colab** — the Linux package works there with nothing installed
  locally at all, and so does the Python package (the `cp313` wheel).

See [Installation](install_lsdtt3.md) for all of these.
