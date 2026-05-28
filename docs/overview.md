# Overview

lsdtt3 is a redesign of [LSDTopoTools2](https://github.com/LSDtopotools/LSDTopoTools2) that retains the same scientific algorithms while rebuilding the codebase for modern C++17, better memory efficiency, and easier ongoing development.

---

## What it does

lsdtt3 provides a collection of command-line programs for analysing digital elevation models (DEMs). Each program reads a plain-text parameter file and writes raster or vector outputs. The tools cover:

- **DEM preparation** — reprojection, depression filling, breaching, nodata handling.
- **Surface morphology** — curvature, roughness, hillshade, local relief, cloth-drape upper envelope, Frangi vessel filter.
- **Drainage and flow** — D8, D-infinity, and multiple-flow-direction routing; topographic wetness index.
- **Channel networks** — extraction from a contributing-area threshold, junction and link construction, Strahler ordering, export to geospatial formats.
- **Basin analysis** — delineation, statistics, chi analysis, k_sn and m/n estimation.
- **Cosmogenic nuclides** — basin-averaged erosion rates, concentration prediction, transient 1-D column modelling.
- **Landscape evolution modelling** — stream-power incision, hillslope diffusion, tectonic uplift.

---

## Key improvements over LSDTopoTools2

| Feature | LSDTopoTools2 | lsdtt3 |
|---|---|---|
| Raster I/O | BIL/ASC | Any GDAL format (GeoTIFF, COG, VRT, BIL, ASC, …) |
| Vector output | CSV only | GeoPackage, FlatGeobuf, CSV with correct CRS |
| Projection handling | Manual | Automatic UTM / polar stereographic / Lambert LAEA |
| Pixel indexing | 32-bit | 64-bit — handles DEMs > 40 000 × 40 000 px |
| Flow routing | D8 | D8, D-infinity, MFD |
| Memory | Higher | Reduced (8/32-bit types used where appropriate) |
| Parallelism | Limited OpenMP | OpenMP across most bottleneck routines |
| Output format | Standard GeoTIFF | Cloud-Optimised GeoTIFF (COG) by default |

---

## How it is structured

lsdtt3 is organised as a set of free functions grouped by task, compiled into a shared library and linked into individual command-line executables. Because each analysis file is independent, editing one routine does not trigger a full rebuild.

The build system is CMake. Optional dependencies (FFTW3, Eigen3, OpenCV, NetCDF) are detected automatically; tools that need a missing library are skipped rather than failing the build.

See [Architecture](lsdtt3_architecture.md) for the full layer design and naming conventions.

---

## Who it is for

lsdtt3 is aimed at geomorphologists, physical geographers, and earth scientists who need to:

- process large or many DEMs in batch workflows,
- run reproducible analyses from parameter files,
- integrate topographic analysis into scripted pipelines.

It is not a GIS — it has no interactive GUI. Inputs and outputs are files; visualisation is done in QGIS, Python, or R.

---

## Compatibility

lsdtt3 is **not backward-compatible** with LSDTopoTools2. Parameter files, output file names, and some algorithm defaults differ. If you have existing LSDTopoTools2 workflows, consult the [command-line reference](cli_reference_lsdtt_channel_extraction.md) for the updated parameter names.

lsdtt3 runs on Linux. Windows users run it inside WSL2; macOS users inside a Multipass Ubuntu VM. See [Installation](installing_locally.md) for details.
