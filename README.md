# lsdtt3

![LSDTopoTools3 logo](docs/LSD-logo.png)

A redesign of [LSDTopoTools2](https://github.com/LSDtopotools/LSDTopoTools2) with a modern C++ architecture.

**Copyright:** Simon M. Mudd  
**Licence:** GPL-3.0  
**Additional developers:** Rutger Vervoordeldonk, Boris Gailleton

**Additional authors of ported code:** Fiona Clubb, Stuart Grieve, David Milodowski, Martin Hurst, Declan Valters 

---

## Overview

lsdtt3 is a C++ library and set of command-line tools for topographic analysis, and various geomorphic and hydrologic analyses. It follows on from [LSDTopoTools2](https://github.com/LSDtopotools/LSDTopoTools2) but has been rewritten from the ground up, while retaining the core scientific functionality of LSDTopoTools2. `lsdtt3` has updated architecture more suited to modern C++17, optional GPU support, and easier incremental development.

If you use this software for research, please also cite the original LSDTopoTools2:

> Mudd, S. M., Clubb, F. J., Grieve, S. W. D., Milodowski, D. T., Gailleton,
> B., Hurst, M. D., Valters, D. V., Wickert, A. D., & Hutton, E. W. H. (2023).
> LSDtopotools/LSDTopoTools2: LSDTopoTools2 v0.9. Zenodo.
> https://doi.org/10.5281/zenodo.8076231

**lsdtt3 is not backward-compatible with LSDTopoTools2.**

Key improvements over LSDTopoTools2:

- Reads any GDAL-supported format (GeoTIFF, COG, VRT, BIL, ASC, …). VRT is especially useful for large tiled datasets.
- A command line interface for downloading data from OpenTopography.org. No additional sofware is needed to get datasets for analysis. 
- Writes vector outputs (GeoPackage, FlatGeobuf, CSV) with correct georeferencing via GDAL/OGR.
- Automatic projection handling: If input data is not projected, converts to UTM below ±75°, polar stereographic at high latitude, optional Lambert LAEA for wide-extent rasters spanning more than two UTM zones.
- If input raster is projected, keeps the projection. So regional projections (e.g., British National Grid) can be retained for both input and output. Can also ingest projected coordinates from other planets (e.g., Mars or Titan)
- Cloud Optimised GeoTIFF (COG) output by default. Calculates internally time savings of COG format so small rasters are output as standard GeoTIFF.
- 64-bit pixel indexing — can process DEMs larger than 40,000 × 40,000 pixels. You will need a very high memory compute node for this! But the authors wanted to retain the ability to analyse the largest drainage basins on Earth at high resolution without resorting to tiling. Elswhere data types are used to reduce memory consumption (8 bit, and 32 bit data used where approriate).
- D8 and DAG (D-infinity / MFD) flow routing. D8 retains the relationships between pixels and their upslope neigbours through an object called the s-vector for O(1) upslope-set queries. 
- CPU and GPU backends behind a shared abstract interface (GPU functionality is significantly reduced compared to CPU, however).
- Analysis routines as stateless free functions as opposed to objects. The reason this was done was to reduce the amount of spaghetti code in the package, and in addition reduce the chances of needing a full rebuild any time a small modification is made, problems both present in LSDTopoTools2. Analyses are also clustered by task so (hopefully) it is easier to find the functions you want (although a side effect of this is that there are many more command line interfaces than previously). In other words, we have not thrown the kitchen sink at something called lsdtt-basic-metrics this time (it no longer exists). The names of the CLI tools should give users a reasonable idea of what they are meant to do. 
- OpenMP multithreading for selected routines; many operations are substantially faster than the equivalent LSDTopoTools2 code. We have a testing and benchmarking document that has some concrete numbers attached to this statement, see it in the `/docs` directory. 
- Many algorithms made more robust and/or faster (e.g., Spectral analysis, proximity tools, metrics based on moving windows). 
- Overall, users should have more flexibility/ease for file I/O, should see less memory consumption (and as a result be able to process bigger datasets and have fewer segmentation faults), and should see speedup of a number of key algorithms.

## LLM Usage:

Developers of `lsdtt3` have used various LLMs (GPT Codex 5.3, GPT 5.6 Sol, Terra and Luna, Claude Haiku, Sonnet, and Opus) and to port legacy code and write new routines, with testing against baselines (e.g., lsdtt2 outputs, previous iterations of algorithms) and regular code review. 

---

## Installation

`lsdtt3` is currently in **private beta**.

Everything below installs a pre-built release package: no compiler, no GDAL
install, no `sudo`. **Full instructions, including Windows, macOS, Google Colab
and troubleshooting, are in [`docs/install_lsdtt3.md`](docs/install_lsdtt3.md).**

Packages are on the [Releases page](../../releases/latest) and at
<https://www.geos.ed.ac.uk/~smudd/lsdtt_packages/>.

> **Check the version number before you copy the commands below.** The version
> in these examples is whatever was current when the page was last edited, and
> it goes stale. The authoritative current version is line 2 of
> [`CMakeLists.txt`](CMakeLists.txt) — `project(lsdtt3 VERSION ...)` — which is
> what release tags and package filenames are built from. Substitute that
> number wherever you see one here.

### Linux (and Windows via WSL2, and macOS via a Multipass VM)

```bash
wget https://www.geos.ed.ac.uk/~smudd/lsdtt_packages/lsdtt3-backend-linux-x86_64-core-v0.4.0.tar.gz
tar -xzf lsdtt3-backend-linux-x86_64-core-v0.4.0.tar.gz
source lsdtt3-backend-linux-x86_64-core-v0.4.0/setenv.sh
lsdtt-channel-extraction --help
```

`setenv.sh` must be **`source`d**, not run with `bash` — it only exports
environment variables, which would vanish with a child shell. It applies to the
current terminal; append the same line (with the full path) to `~/.bashrc` to
make it permanent.

On **Windows**, first `wsl --install` from an administrator PowerShell, install
**Ubuntu 22.04 LTS** from the Microsoft Store, then run the above inside it.
Your Windows drives are visible under `/mnt/c/`. A native Windows `.zip` also
exists but ships a reduced tool set — see the install guide.

On **macOS** there is no package yet; run the Linux one inside a Multipass
Ubuntu VM (`brew install --cask multipass`).

On **Google Colab**, `!source setenv.sh` does not work — set the environment
from Python instead. The install guide has the cell to paste.

### Want the latest code rather than a release?

Then you need to compile it:
[`docs/developer_notes/building_from_source.md`](docs/developer_notes/building_from_source.md).

### Reporting bugs

**WARNING** This is software built by academics for academic use within a research group. It is not supported commercial software. The developers are happy to help if you have questions, comments, bug reports etc, but keep in mind the developers may be working on other things (teaching, doing research, etc) and adjust expectations of response timescales accordingly. 

Please include:
- the exact command you ran
- the full terminal output (copy and paste, not a screenshot)
- your Linux version: `cat /etc/os-release | grep VERSION=`
- the lsdtt3 version: `lsdtt-channel-extraction -v`

---

## Dependencies

**You do not need any of this to use a release package** — the packages bundle
every library they need. This section is for building from source.

`lsdtt3` has more dependencies than `lsdtopotools2`. However all dependencies are stabile and easily installable on Linux (and by extension linux subsystems for Windows and MacOS). We also provide an apptainer container (although not frequently updated) that allows running of the softare without having to manage installation. 

| Dependency | Required | Enables |
|---|---|---|
| CMake ≥ 3.15 | Yes | Build system |
| C++17 compiler (GCC ≥ 9 or Clang ≥ 10) | Yes | All |
| GDAL ≥ 3.0 | Yes | All raster and vector I/O |
| OpenMP | No | Multi-core parallelism across most routines |
| FFTW3 (double precision) | No | `lsdtt-spectral-analysis`; wavelet quadratic polyfit in other tools |
| Eigen3 | No | `lsdtt-chi-analysis`; chi k_sn / u* inversion |
| OpenCV (core + imgproc) | No | Frangi filter, Otsu segmentation, Hough line extraction, adaptive thresholding in `lsdtt-surface-metrics` |
| NetCDF | No | ERA40 atmospheric / palaeomagnetic corrections in `lsdtt-cosmo-tool` and `lsdtt-cosmo-prediction` |
| CUDA toolkit | No | GPU compute backend (not yet implemented) |

### Ubuntu / Debian

```sh
# Required
sudo apt update
sudo apt install cmake build-essential libgdal-dev

# Optional — install whichever you need
sudo apt install libfftw3-dev        # spectral analysis, wavelet polyfit
sudo apt install libeigen3-dev       # chi analysis and k_sn inversion
sudo apt install libopencv-dev       # image analysis in surface metrics
sudo apt install libnetcdf-dev       # CRN palaeomagnetic corrections
```

OpenMP is included with GCC and needs no separate install.

### Windows (WSL2)

Build inside a WSL2 Ubuntu environment and use the Ubuntu instructions above.
`docs/developer_notes/lsdtt3_in_windows.md` has the full walkthrough. A native
Windows build is possible with Visual Studio and vcpkg but is a deliberately
reduced profile (no Eigen, no OpenCV, no OpenMP) intended as the QGIS plugin
backend — see `docs/developer_notes/building_from_source.md`.

### macOS

lsdtt3 does not build with Apple Clang. Build inside an Ubuntu VM managed by
Multipass — see `docs/developer_notes/lsdtt3_in_macos.md`.

### Apptainer container

A pre-configured Apptainer definition is provided that includes all dependencies (GDAL, FFTW3, Eigen3, OpenCV, OpenMP) on Ubuntu 22.04. See `apptainer_def/apptainer_setup.md` for setup on Linux, macOS, and Windows (WSL2).

### Debian packaging

See `docs/developer_notes/debian_packaging.md` for building a `.deb` package for distribution on Debian-based systems.

---

## Building

```sh
git clone https://github.com/LSDtopotools/lsdtt3.git
cd lsdtt3
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

Compiled executables are written to `build/`. To use them from any directory:

```sh
export PATH=/path/to/lsdtt3/build:$PATH
```

Full instructions for every platform, including the CMake option reference:
[`docs/developer_notes/building_from_source.md`](docs/developer_notes/building_from_source.md).

### CMake options

| Option | Default | Description |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Release`, `Debug`, or `RelWithDebInfo` |
| `LSDTT_GPU` | `OFF` | Enable GPU compute backend (requires CUDA toolkit) |
| `LSDTT_OPENCV` | `ON` | Build OpenCV-dependent routines (auto-disabled if OpenCV not found) |
| `LSDTT_EIGEN` | `ON` | Build Eigen3-dependent routines (auto-disabled if Eigen3 not found) |
| `LSDTT_ENABLE_BENCHMARKS` | `OFF` | Register the RichDEM comparison benchmarks (needs network + pixi) |
| `LSDTT_ENABLE_NETWORK_TESTS` | `OFF` | Register tests that download data |
| `LSDTT_ENABLE_PERF_TESTS` | `OFF` | Register flow-routing performance benchmarks |

---

## Testing

Tests are driven by CTest. The default suite is self-contained: it needs no
network, no downloaded data, and no Python packages beyond a plain interpreter.

Short version on Python/pixi:

- Running `lsdtt3` CLIs does **not** need Python.
- Most tests (smoke, failure, unit, rigorous) do **not** need pixi.
- Some opt-in benchmark tests do need pixi because their references are Python
  tools (RichDEM / TopoToolbox and DAGGER adapters).

```sh
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Tests are grouped by label, so you can run just the part you care about:

| Label | What it covers |
|---|---|
| `smoke` | Each CLI runs end to end on a tiny fixture and produces its expected output |
| `rigorous` | Exact checks on outputs — raster checksums, sizes, feature counts, file hashes |
| `failure` | Bad parameters are rejected with the expected error message |
| `unit` | C++ unit tests: parameter parser, cosmogenic physics, flow-routing properties |
| `cosmo` | Cosmogenic-nuclide physics against CRONUS and pinned numerics sentinels |
| `flow` | Flow routing and depression handling: properties, and RichDEM comparison |
| `fetch` | Native DEM fetching (`lsdtt_raster_fetch`); offline logic by default |
| `performance` | In-process timings of the flow routines (opt-in; reported, not asserted) |
| `benchmark` | Opt-in comparison against external references (RichDEM/TopoToolbox flow, DAGGER hydraulic) on fetched real DEMs |

```sh
ctest --test-dir build -L smoke --output-on-failure
ctest --test-dir build -L unit  --output-on-failure
```

`benchmark` tests are **off by default** because they need network access and an
optional RichDEM install; enable with `-DLSDTT_ENABLE_BENCHMARKS=ON`. When their
dependencies are missing they report as skipped rather than failed, so an
optional tool can never turn the build red.

If you enable `LSDTT_ENABLE_BENCHMARKS`, install pixi first so those benchmark
dependencies are available in the benchmark environment.

### Timing lsdtt3 on your own machine

The timing tiers are slow by design — they exist to measure, not to check — so
keep them in a **separate Release build** and leave `build/` fast:

```sh
cmake -S . -B build-perf -DCMAKE_BUILD_TYPE=Release -DLSDTT_ENABLE_PERF_TESTS=ON
cmake --build build-perf -j
OMP_NUM_THREADS=8 ctest --test-dir build-perf -L performance --output-on-failure
```

Results land as CSV in `build-perf/tests/perf_output/` (`flow_performance.csv`,
`lsdn_performance.csv`), each row carrying a git SHA so runs are comparable
across commits. Or drive them directly, which is the useful mode on a server:

```sh
build-perf/tests/lsdtt-test-flow-performance out_dir --full \
    --raster tests/data/dem_cache/alpine_small.tif
```

Three things matter for numbers you can trust: build `Release` (otherwise you
are timing the compiler), pin `OMP_NUM_THREADS` (otherwise a busy node is
unrepeatable), and don't pass `-j` to CTest for this label (parallel tests and
in-process timing measure each other). These tests are reported, never asserted,
so a loaded machine cannot turn the suite red.

Full methodology, the reference numbers to compare against, and the
thread-scaling caveats are in `docs/testing_and_benchmarking.md` §6.

**No DEM is committed to this repository.** Tests that need real topography fetch
it on demand from public keyless sources and cache it under
`tests/data/dem_cache/` (gitignored), using `lsdtt-fetch-raster` itself.

See `tests/README.md` for what each test does, and
`docs/testing_and_benchmarking.md` for methodology, benchmark results, and the
places where lsdtt3 and RichDEM deliberately differ.

---

## Command-line programs

All programs take a plain-text parameter file as their argument. Parameter files use `key: value` syntax, one entry per line; lines beginning with `#` are comments. Pass `-h` to any program for a full parameter reference.


```sh
lsdtt-raster-preprocessing  my_params.driver
lsdtt-drainage-analysis     my_params.driver
```

**Detailed parameter references for each program are in `docs/`.**

### DEM preparation

| Program | Description | Reference |
|---|---|---|
| `lsdtt-fetch-raster` | Downloads a DEM for a WGS84 bounding box and writes it reprojected to UTM, ready to analyse. The default source is Copernicus GLO-30 on AWS Open Data — public and keyless, so no account or API key is needed. The area can be given as `west`/`south`/`east`/`north`, or as a vector file (`.fgb`, `.gpkg`) containing a polygon or two corner points — reprojected from UTM automatically, so a box drawn in a GIS works directly. Bounding boxes spanning several source tiles are mosaicked automatically. OpenTopography is also supported for its other datasets; that needs a free API key, supplied via the `OT_API_KEY` environment variable (or `OT_API_KEY_FILE`) — **never** in a parameter file, since parameter files are commonly committed and every parameter is echoed at the end of a run. | `docs/cli_reference/cli_reference_lsdtt_fetch_raster.md` |
| `lsdtt-raster-preprocessing` | Reprojects geographic/Web Mercator DEMs to UTM/polar, applies elevation threshold masking, fills interior nodata voids, clips to maximum valid rectangle, breaches depressions (Lindsay 2016/2020), fills pits (Barnes 2014), writes hillshade. | `docs/cli_reference/cli_reference_lsdtt_raster_preprocessing.md` |
| `lsdtt-raster-synthetic-tests` | Generates synthetic test surfaces (random, fractal, diamond-square, mountain-range geometries) for algorithm validation. | — |

### Surface analysis

| Program | Description | Reference |
|---|---|---|
| `lsdtt-surface-metrics` | Hillshade, local relief, polynomial surface fit (direct and wavelet), curvature metrics, MADk2 roughness, Frangi filter, Otsu segmentation, Hough line extraction, curvature skeletons, cloth-drape upper envelope and incision depth. | `docs/cli_reference/cli_reference_lsdtt_surface_metrics.md` |
| `lsdtt-flat-surface-finder` | Identifies low-texture surfaces (alluvial fans, terraces, piedmonts) using MADk2, residual IQR, and normal-vector dispersion (NVD). | — |
| `lsdtt-mountain-front` | Identifies mountain fronts from multi-scale relief, piedmont surface fitting, and residual thresholding. Supports synthetic DEM validation mode. | `docs/cli_reference/cli_reference_lsdtt_mountain_front.md` |
| `lsdtt-feature-finder` | Feature-detection workflows (terraces, piedmonts, hillslope morphology). | — |
| `lsdtt-spectral-analysis` *(requires FFTW3)* | Spectral surface analysis using the FFTW double API. | — |

### Flow routing and drainage

| Program | Description | Reference |
|---|---|---|
| `lsdtt-drainage-analysis` | Drainage area via D8, D-infinity, or MFD routing. Topographic Wetness Index (TWI). | `docs/cli_reference/cli_reference_lsdtt_drainage_analysis.md` |
| `lsdtt-channel-extraction` | D8 flow routing, contributing-pixel threshold source extraction, junction/channel network construction, network export to CSV / GeoPackage / FlatGeobuf. | `docs/cli_reference/cli_reference_lsdtt_channel_extraction.md` |
| `lsdtt-basin-extraction` | D8 flow routing, channel/junction network, junction-aware basin extraction. Writes basin outlines, raster labels, and statistics. | `docs/cli_reference/cli_reference_lsdtt_basin_extraction.md` |
| `lsdtt-basin-stats` | Hydromorphic statistics for outlet-point-defined basins: elevation, slope, flow-distance percentiles; centroid; drainage area; channel length by Strahler order. | `docs/cli_reference/cli_reference_lsdtt_basin_stats.md` |
| `lsdtt-graphflood` | Graph-based hydraulic routing (Gailleton et al. 2024) on a fixed D-infinity or MFD flow graph. Writes water depth, water surface, and flux rasters. | — |
| `lsdtt-stream-burning` | Loads a polyline channel network, seeds D8 routing from its source tips, and writes the DEM-extracted network for comparison before burning. | — |

### Channel network analysis

| Program | Description | Reference |
|---|---|---|
| `lsdtt-channel-network-analysis` | Channel-link bearing fitting, junction-angle calculation (donor opening, bending angle), Strahler/Tokunaga self-similarity analysis. | `docs/cli_reference/cli_reference_lsdtt_channel_network_analysis.md` |
| `lsdtt-chi-analysis` *(requires Eigen3)* | End-to-end chi workflow: chi extraction, piecewise-linear network fitting, disorder metrics with jackknife/bootstrap uncertainty, spatial k_sn / u* inversion with m/n model selection. | `docs/cli_reference/cli_reference_lsdtt_chi_analysis.md` |
| `lsdtt-cross-sections` | Cross-section extraction along channel networks or user-defined transects. | `docs/cli_reference/cli_reference_lsdtt_cross_sections.md` |

### Tidal / marsh

| Program | Description | Reference |
|---|---|---|
| `lsdtt-marsh-topography` | Marsh platform identification, tidal channel extraction, skeleton-based channel metrics. | `docs/cli_reference/cli_reference_lsdtt_marsh_topography.md` |

### Cosmogenic nuclides

| Program | Description | Reference |
|---|---|---|
| `lsdtt-cosmo-tool` | Basin-averaged CRN erosion-rate inversion. Reads measured nuclide concentrations, computes Stone 2000 production, snaps points to D8 network, inverts via Newton-Raphson + bisection. Optional ERA40 / palaeomagnetic correction (requires NetCDF). | `docs/cli_reference/cli_reference_lsdtt_cosmo_tool.md` |
| `lsdtt-cosmo-prediction` | Predicts CRN concentrations across a DEM: steady-state, step-change transient, and erosion-weighted accumulated concentration. | `docs/cli_reference/cli_reference_lsdtt_cosmo_prediction.md` |
| `lsdtt-cosmo-column-tool` | 1-D transient cosmogenic nuclide column simulation. Steady-state, step-change, and full-transient (user-supplied erosion history) modes. No DEM required. | `docs/cli_reference/cli_reference_lsdtt_cosmo_column_tool.md` |

### Landscape evolution

| Program | Description |
|---|---|
| `lsdtt-lem` | Landscape evolution model: D-infinity flow routing, implicit stream-power fluvial incision (Anand et al. 2020 / Braun & Willett 2013), Cidre linear/nonlinear hillslope transport (Carretier et al. 2016), uniform uplift. Runs from a synthetic flat or noise surface, or an input DEM. |

### Utilities

| Program | Description |
|---|---|
| `lsdtt-linear-stats` | Statistical summaries for linear (directional/bearing) data. |

---

## Other Documentation

The `docs/` tree is organised by audience:

| Directory | Audience |
|---|---|
| `docs/` (top level) | **Users** — installation and overview. Also the source for the documentation site (`mkdocs.yml`). |
| `docs/cli_reference/` | **Users** — one parameter reference per command-line tool. |
| `docs/developer_notes/` | **Developers** — architecture, building, packaging. |
| `docs/plan_docs/` | Working design documents. Not published, not committed. |
| `docs/tex_technical/` | LaTeX sources for method write-ups with full equations. |

### For users

| Document | Contents |
|---|---|
| `docs/install_lsdtt3.md` | **Start here.** Installing a release package on Linux, Windows (WSL2 or native), macOS, or Google Colab. |
| `docs/overview.md` | What lsdtt3 is, what it does, and how it differs from LSDTopoTools2. |
| `docs/cli_reference/` | Parameter reference for every tool — see the table above for which file. |

### For developers

| Document | Contents |
|---|---|
| `docs/developer_notes/lsdtt3_architecture.md` | Full architecture reference: layer design, naming conventions, all module APIs. Start here for development. |
| `docs/developer_notes/developer_onboarding.md` | Getting-started guide for new contributors. |
| `docs/developer_notes/git_notes.md` | Setting up an SSH key and cloning, for people new to git. |
| `docs/developer_notes/building_from_source.md` | Compiling on Linux, Windows (WSL2 and native) and macOS; CMake options. |
| `docs/developer_notes/making_release_packages.md` | How the release packages are built and published. |
| `docs/developer_notes/building_the_docs.md` | Building the documentation site with MkDocs, and publishing it. |
| `docs/developer_notes/backend_abstraction_rules.md` | When to use abstract bases + factories vs. concrete CPU free functions. |
| `docs/developer_notes/using_lsdtt_flowinfo.md` | How to use the D8 flow routing object in analysis code. |
| `docs/developer_notes/lsdtt_utils_guide.md` | Guide to the `src/utils/` layer: what belongs there and the valid-data contract. |
| `docs/developer_notes/gpu_compute.md` | GPU backend design and the compute-dispatch mechanism. |
| `docs/developer_notes/beta_testing_notes.md` | Getting bundles to beta testers, the release pipeline, and the future conda-forge plan. |
| `docs/developer_notes/lsdtt3_in_windows.md` | Windows (WSL2) development setup. |
| `docs/developer_notes/lsdtt3_in_macos.md` | macOS development setup via a Multipass Ubuntu VM. |
| `docs/developer_notes/debian_packaging.md` | Building a Debian package. |
| `docs/testing_and_benchmarking.md` | How lsdtt3 is tested and benchmarked, and where it deliberately differs from RichDEM. Read before comparing lsdtt3 output against another tool. |
| `tools/README.md` | Building redistributable bundles, and the git/gh commands for cutting a tagged release. |
| `apptainer_def/apptainer_setup.md` | Apptainer container setup. |

---

## Python helpers

We also provide various python helpers (in `python_helpers/`) for simple visualisation and downloading of DEMs. Maintaining python environments is very time consuming and annoying, which is why the core tools are written in c++. We provide helper functions for users familiar with maintaining a geospatial python environment, but if you don't know how to do this just use a GIS and the internet to get and look at `lsdtt3` data inputs and outputs.  

| Script | Purpose |
|---|---|
| `download_dem_tiles.py` | Download DEMs from Planetary Computer (keyless), OpenTopography, a WCS endpoint, or Copernicus Data Space. Note that we also have a command line interface, `lsdtt-fetch-rasters` that also does this. |
| `tif_to_png.py` | Convert a GeoTIFF to a PNG for quick inspection |
| `osm_download_features.py`, `osm_water_download.py` | Download OpenStreetMap features / water layers |
| `download_gem_faults.py` | Download the GEM active-fault database |
| `download_macrostrat_mvt.py` | Download Macrostrat geology vector tiles |

Each has a `--help`; `dem_downloader_instructions.md` and
`tif_to_png_instructions.md` give worked examples.

### Python environments with pixi

None of the above is needed to build or run lsdtt3 — the C++ tools have no
Python dependency. For downloading DEMs in particular, `lsdtt-fetch-raster`
does the common case (a bounding box from Copernicus GLO-30) with no Python at
all; `download_dem_tiles.py` remains for the other providers. If you do want the helpers or the test tooling, we use
[pixi](https://pixi.sh) so you don't have to hand-maintain a geospatial Python
environment. [lsdviztools3](https://github.com/LSDtopotools/lsdviztools3), our
visualisation package, uses pixi on the same conda-forge channel, so the two
environments compose the same way.

```sh
curl -fsSL https://pixi.sh/install.sh | bash     # once
cd lsdtt3
pixi install                                     # solve the default environment
pixi run python python_helpers/download_dem_tiles.py --help
pixi shell                                       # or drop into a shell
```

The manifest is `pixi.toml`. It covers `python_helpers/` and `tests/scripts/`;
heavier or optional packages sit in separate environments so the default solve
stays quick:

| Environment | Adds | For |
|---|---|---|
| default | GDAL, numpy, rasterio, geopandas, shapely, requests, STAC clients, pytest | helpers and test tooling |
| `osm` | osmnx | the OpenStreetMap helpers |
| `benchmark` | RichDEM, scipy | reference implementations used to validate lsdtt3 |

```sh
pixi run -e osm python python_helpers/osm_water_download.py --help
```

Install RichDEM from conda-forge via pixi rather than pip: `pip install richdem`
installs a version that pins `numpy<2` and will downgrade numpy in whatever
environment it lands in.

## References

Mudd, S. M., et al. (2023). LSDtopotools/LSDTopoTools2 v0.9. Zenodo. https://doi.org/10.5281/zenodo.8076231

Barnes, R., Lehman, C., Mulla, D. (2014) Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. *Computers & Geosciences*, 62, 117–127.

Gailleton, B., et al. (2024) Graph-based hydrodynamic modelling. *(see lsdtt-graphflood)*

Lindsay, J.B. (2016) Efficient hybrid breaching-filling sink removal methods for flow path enforcement in digital elevation models. *Hydrological Processes*, 30(6), 846–857.

Lindsay, J.B. (2020) Pit-centric depression removal methods. *Geomorphometry 2020 Conference Proceedings*.

Stone, J.O. (2000) Air pressure and cosmogenic isotope production. *Journal of Geophysical Research*, 105(B10), 23753–23759.

Tarboton, D.G. (1997) A new method for the determination of flow directions and upslope areas in grid digital elevation models. *Water Resources Research*, 33(2), 309–319.

Zhang, W., et al. (2016) An easy-to-use airborne LiDAR data filtering method based on cloth simulation. *Remote Sensing*, 8(6), 501.

Zevenbergen, L.W. and Thorne, C.R. (1987) Quantitative analysis of land surface topography. *Earth Surface Processes and Landforms*, 12, 47–56.
