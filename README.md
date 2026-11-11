# lsdtt3

![Drainage divides of third- and fourth-order basins over shaded relief, drawn with lsdviztools3 from lsdtt3 output](docs/images/readme_banner.png)

[![Release](https://img.shields.io/github/v/release/LSDtopotools/lsdtt3)](https://github.com/LSDtopotools/lsdtt3/releases/latest)
[![Release build](https://img.shields.io/github/actions/workflow/status/LSDtopotools/lsdtt3/release.yml?label=release%20build)](https://github.com/LSDtopotools/lsdtt3/actions/workflows/release.yml)
[![Docs](https://img.shields.io/badge/docs-lsdtopotools.github.io-blue)](https://lsdtopotools.github.io/lsdtt3/)
[![Licence: AGPL-3.0-or-later](https://img.shields.io/github/license/LSDtopotools/lsdtt3)](LICENSE)
[![Docs build](https://img.shields.io/github/actions/workflow/status/LSDtopotools/lsdtt3/docs.yml?label=docs%20build)](https://github.com/LSDtopotools/lsdtt3/actions/workflows/docs.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24110423.svg)](https://doi.org/10.5281/zenodo.24110423)

A ground-up redesign of [LSDTopoTools2](https://github.com/LSDtopotools/LSDTopoTools2). Like LSDTopoTools2 but better in almost every way. lsdtt3 was built with LLM assistance so you will no longer get amazing lsdtt2 messages like "This index is impossible, how did I get here?". We apologise if the bizarre messages were the primary reason you used lsdtt2. If your goal is topographic analysis, or geomorphic research, you will find lsdtt3 is better. 

**Copyright:** Simon M. Mudd  
**Licence:** AGPL-3.0-or-later  
**Additional developers:** Rutger Vervoordeldonk, Boris Gailleton, Martha Hawksworth

**Additional authors of ported lsdtt2 code:** Fiona Clubb, Stuart Grieve, David Milodowski, Martin Hurst, Declan Valters 

**Documentation: <https://lsdtopotools.github.io/lsdtt3/>** — installation,
guides, the design and algorithms, and a reference page for every
command-line tool. Each tool also documents itself: `<tool> -h` lists every
parameter it accepts, and `<tool> -h <parameter_name>` explains one.

---

## Overview

lsdtt3 is a C++ library, a set of command-line tools and a Python package for topographic analysis, and various geomorphic and hydrologic analyses. It follows on from [LSDTopoTools2](https://github.com/LSDtopotools/LSDTopoTools2) but has been rewritten from the ground up, while retaining the core scientific functionality of LSDTopoTools2. `lsdtt3` has updated architecture more suited to modern C++17, optional GPU support (well it will when we write it).

If you use this software for research, please also cite the original LSDTopoTools2:

> Mudd, S. M., Clubb, F. J., Grieve, S. W. D., Milodowski, D. T., Gailleton,
> B., Hurst, M. D., Valters, D. V., Wickert, A. D., & Hutton, E. W. H. (2023).
> LSDtopotools/LSDTopoTools2: LSDTopoTools2 v0.9. Zenodo.
> https://doi.org/10.5281/zenodo.8076231

**lsdtt3 is not backward-compatible with LSDTopoTools2.**

Key improvements over LSDTopoTools2:

- Reads any GDAL-supported format (GeoTIFF, COG, VRT, BIL, ASC, …). VRT is especially useful for large tiled datasets.
- A command line interface for downloading data from OpenTopography.org. No additional sofware is needed to get datasets for analysis. Theer are also a number of additional download tools for lithology, erosion rates, and various vector files. 
- Writes vector outputs (GeoPackage, FlatGeobuf, CSV) with correct georeferencing via GDAL/OGR.
- Automatic projection handling: If input data is not projected, converts to UTM below ±75°, polar stereographic at high latitude, optional Lambert LAEA for wide-extent rasters spanning more than two UTM zones.
- If input raster is projected, keeps the projection. So regional projections (e.g., British National Grid) can be retained for both input and output. Can also ingest projected coordinates from other planets (e.g., Mars or Titan)
- Cloud Optimised GeoTIFF (COG) output by default. Calculates internally time savings of COG format so small rasters are output as standard GeoTIFF.
- 64-bit pixel indexing — can process DEMs larger than 40,000 × 40,000 pixels. You will need a very high memory compute node for this! But the authors wanted to retain the ability to analyse the largest drainage basins on Earth at high resolution without resorting to tiling. Elswhere data types are used to reduce memory consumption (8 bit, and 32 bit data used where approriate).
- D8 and DAG (D-infinity / MFD) flow routing. D8 retains the relationships between pixels and their upslope neigbours through an object called the s-vector for O(1) upslope-set queries. 
- CPU and GPU backends behind a shared abstract interface (GPU functionality is significantly reduced compared to CPU, however).
- Analysis routines as stateless free functions as opposed to objects. The reason this was done was to reduce the amount of spaghetti code in the package, and in addition reduce the chances of needing a full rebuild any time a small modification is made, problems both present in LSDTopoTools2. Analyses are also clustered by task so (hopefully) it is easier to find the functions you want (although a side effect of this is that there are many more command line interfaces than previously). In other words, we have not thrown the kitchen sink at something called lsdtt-basic-metrics this time (it no longer exists). The names of the CLI tools should give users a reasonable idea of what they are meant to do. 
- OpenMP multithreading for selected routines; many operations are substantially faster than the equivalent LSDTopoTools2 code. We have a testing and benchmarking document that has some concrete numbers attached to this statement, see it in the `/docs` directory. 
- A Python package (`import lsdtt3`) that wraps the same C++ routines and keeps rasters and results in memory, so a script or notebook can chain analyses without parameter files. Where a Python function and a CLI share an option, they call the same C++ code and give the same result.
- Many algorithms made more robust and/or faster (e.g., Spectral analysis, proximity tools, metrics based on moving windows, landscape evolution models). 
- Overall, users should have more flexibility/ease for file I/O, should see less memory consumption (and as a result be able to process bigger datasets and have fewer segmentation faults), and should see speedup of a number of key algorithms.

## LLM Usage:

Developers of `lsdtt3` have used various LLMs (GPT Codex 5.3, GPT 5.6 Sol, Terra and Luna, Claude Haiku, Sonnet, and Opus) and to port legacy code and write new routines, with testing against baselines (e.g., lsdtt2 outputs, previous iterations of algorithms) and regular code review. If you have taken the position that all LLM use is evil and all software should be hand crafted, we recommend the previous version, LSDTopoTools2. We also recommend this version if you like unexplained segmentation faults and/or waiting a long time for things to finish. 

---

## Installation

Before you read this section, if you don't like installing software, just go to [our notebooks repositories](https://github.com/LSDtopotools/lsdtt3_notebooks) to see examples where you can just run some cells in google colab and everything should just work. 

Still here? Great. Everything below installs a pre-built release package: no compiler, no GDAL
install, no `sudo`. **Full instructions, including Windows, macOS, Google Colab
and troubleshooting, are in [the installation guide](https://lsdtopotools.github.io/lsdtt3/what_is_lsdtt3/install_lsdtt3.html).**

**Get the packages from [GitHub Releases](../../releases/latest).** The
commands below always fetch the latest release, so there is no version number
to look up. (A copy of some packages on Simon Mudd's university web space is
for when GitHub is unreachable only: it is updated by hand and usually well out
of date.)

### Linux (and Windows via WSL2, and macOS via a Multipass VM)

```bash
wget https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-linux-x86_64-core.tar.gz
tar -xzf lsdtt3-backend-linux-x86_64-core.tar.gz
source lsdtt3-backend-linux-x86_64-core-v*/setenv.sh    # the folder carries the version
lsdtt-channel-extraction --help
```

`setenv.sh` must be **`source`d** (from bash or zsh), not run with `bash` — it
only exports environment variables, which would vanish with a child shell. It
applies to the current terminal; append the same line (with the full path) to
`~/.bashrc` to make it permanent.

On **Windows**, first `wsl --install` from an administrator PowerShell, install
**Ubuntu 22.04 LTS** from the Microsoft Store, then run the above inside it.
Your Windows drives are visible under `/mnt/c/`. A native Windows `.zip` also
exists but ships a reduced tool set — see the install guide.

On **macOS**, run the Linux package inside a Multipass Ubuntu VM
(`brew install --cask multipass`) for the full tool set. Releases from 0.9.1
also include native macOS packages (Apple Silicon and Intel, macOS 15 or later)
with the reduced QGIS tool set, as the backend for the QGIS plugin, and the
Python package runs natively on macOS. The install guide has the steps for
each, including clearing macOS's quarantine flag on the native package.

On **Google Colab**, `!source setenv.sh` does not work — set the environment
from Python instead. The install guide has the cell to paste.

The **QGIS plugin** is [`lsdtt3-qgis-plugin.zip`](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-qgis-plugin.zip)
from the same release, installed with *Install from ZIP*, plus a backend
package; the install guide has the steps. The backend packages, the QGIS
plugin and the Python package share one version number: files from the same
release go together.

### Python package

`lsdtt3` is not on PyPI yet. Wheels for Linux, macOS and Windows are attached
to every GitHub release from 0.9.2 (earlier releases have them only where the
Python API was supported); each one is for one Python version (`cp312` is
Python 3.12, `cp313` is Python 3.13). Keep the wheel's original file name and
install it into a virtual, conda or pixi environment:

```bash
python -m pip install lsdtt3-0.7.1-cp313-cp313-manylinux_2_28_x86_64.whl
```

Google Colab runs Python 3.13, so it needs the `cp313` wheel. The wheel bundles
GDAL and PROJ. Every function is documented in the package (`help(lsdtt3.FlowInfo)`),
and `print(lsdtt3.helper_for_llm("all"))` prints a compact guide written for you
or an AI assistant such as Gemini in Colab. Worked examples are in the
[Python guide](https://lsdtopotools.github.io/lsdtt3/python/index.html) and the
[API reference](https://lsdtopotools.github.io/lsdtt3/python_api/index.html); building the wheels yourself is in
[Building the Python wheels](https://lsdtopotools.github.io/lsdtt3/developer_notes/building_python_wheels.html).

### Want the latest code rather than a release?

Then you need to compile it:
[Building lsdtt3 from source](https://lsdtopotools.github.io/lsdtt3/developer_notes/building_from_source.html).

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
| NetCDF | No | ERA40 atmospheric / palaeomagnetic corrections in `lsdtt-cosmo-tool` and `lsdtt-cosmo-prediction` |
| CUDA toolkit | No | GPU compute backend (not yet implemented) |

### Ubuntu / Debian

```sh
# Required
sudo apt update
sudo apt install cmake build-essential pkg-config libgdal-dev

# Optional — install whichever you need
sudo apt install libfftw3-dev        # spectral analysis, wavelet polyfit
sudo apt install libeigen3-dev       # chi analysis and k_sn inversion
sudo apt install libnetcdf-dev       # CRN palaeomagnetic corrections
```

OpenMP is included with GCC and needs no separate install.

### Windows (WSL2)

Build inside a WSL2 Ubuntu environment and use the Ubuntu instructions above.
A native Windows build is possible with Visual Studio and vcpkg but is a
deliberately reduced profile (no Eigen, no OpenMP) intended as the
QGIS plugin backend. Both are in [Building on Windows](https://lsdtopotools.github.io/lsdtt3/developer_notes/building_on_windows.html).

### macOS

Build natively with Homebrew and the Xcode Command Line Tools, or inside an
Ubuntu VM managed by Multipass — see [Building on macOS](https://lsdtopotools.github.io/lsdtt3/developer_notes/building_on_macos.html).

### Apptainer container

A pre-configured Apptainer definition is provided that includes all dependencies (GDAL, FFTW3, Eigen3, NetCDF, OpenMP) on Ubuntu 22.04. See `apptainer_def/apptainer_setup.md` for setup on Linux, macOS, and Windows (WSL2).

### Debian packaging

See [Building a Debian package for lsdtt3](https://lsdtopotools.github.io/lsdtt3/developer_notes/debian_packaging.html) for building a `.deb` package for distribution on Debian-based systems.

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
[Building lsdtt3 from source](https://lsdtopotools.github.io/lsdtt3/developer_notes/building_from_source.html).

### Gotcha: more than one GDAL or PROJ on your machine

If you have QGIS, conda, `rasterio` or `geopandas` installed — and most people
doing serious analysis do — then you have several GDALs and several PROJs. Two
things go wrong, and neither announces itself clearly.

**At build time, with a conda environment active**, CMake can find GDAL's
*headers* in `/usr/include/gdal` and GDAL's *library* in the conda environment.
Everything compiles and the link fails with undefined references to functions
whose signatures changed between versions. lsdtt3 detects this and repairs it by
preferring the library next to the headers, but if you see it, the fix is to
point both at one installation:

```sh
cmake -S . -B build \
  -DGDAL_INCLUDE_DIR=<prefix>/include \
  -DGDAL_LIBRARY=<prefix>/lib/libgdal.so
```

The simplest cure is to build with the conda environment **deactivated**, unless
you are deliberately building against conda's GDAL.

**At run time**, PROJ finds its coordinate database through `PROJ_DATA` or
`PROJ_LIB`. If one of those is exported and points at a `proj.db` from a
different PROJ generation, reprojection can return coordinates that are
plausible and wrong rather than failing. Every lsdtt3 tool reports what it is
actually using:

```sh
lsdtt-raster-preprocessing -v
```

```
GDAL 3.8.4, PROJ 9.3.1
  proj data: PROJ built-in search
  GDAL data: GDAL built-in search
```

A release package says `bundled:` there and ignores your environment entirely. A
source build uses whatever your system provides, which is normally correct —
but if coordinates look subtly off, this is the first thing to check, and the
first thing to put in a bug report.

### CMake options

| Option | Default | Description |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Release`, `Debug`, or `RelWithDebInfo` |
| `LSDTT_GPU` | `OFF` | Enable GPU compute backend (requires CUDA toolkit) |
| `LSDTT_EIGEN` | `ON` | Build Eigen3-dependent routines (auto-disabled if Eigen3 not found) |
| `LSDTT_ENABLE_BENCHMARKS` | `OFF` | Register the RichDEM comparison benchmarks (needs network + pixi) |
| `LSDTT_ENABLE_PYTHON_TESTS` | `OFF` | Register the Python-driven tests, most CLI coverage (Python ≥ 3.10, GDAL's command-line programs) |
| `LSDTT_ENABLE_PYTHON_GEO_TESTS` | `OFF` | Also the Python tests needing geospatial packages |
| `LSDTT_ENABLE_NETWORK_TESTS` | `OFF` | Register tests that download data, including `cosmo_realworld` |
| `LSDTT_ENABLE_PERF_TESTS` | `OFF` | Register flow-routing performance benchmarks |

---

## Testing

Tests are driven by CTest, in tiers by what they need to run:

- **Native C++ tests** — registered by default, need only the C++ build. The
  offline, reproducible set:
  `ctest --test-dir build -L '^hermetic$' --output-on-failure`.
- **Python-driven tests** — most CLI coverage; opt in with
  `-DLSDTT_ENABLE_PYTHON_TESTS=ON` (Python ≥ 3.10 and GDAL's command-line
  programs), plus `-DLSDTT_ENABLE_PYTHON_GEO_TESTS=ON` for the few that need
  geospatial packages (easiest from the Pixi `test` environment).
- **Network and benchmark tests** — opt in with
  `-DLSDTT_ENABLE_NETWORK_TESTS=ON` / `-DLSDTT_ENABLE_BENCHMARKS=ON`.

Running `lsdtt3` CLIs never needs Python.

```sh
cmake -S . -B build -DLSDTT_ENABLE_PYTHON_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Configure prints which tiers it registered and what each OFF option left out.
If tests fail, `cmake --build build --target test_report` writes one file to
attach to a bug report. See
[Building from source](https://lsdtopotools.github.io/lsdtt3/developer_notes/building_from_source.html#run-the-tests)
and `tests/README.md` for the details, including running the Python tests
from Pixi.

Tests are grouped by label, so you can run just the part you care about:

| Label | What it covers |
|---|---|
| `smoke` | The program runs to completion and writes its expected outputs |
| `negative` | Invalid input or configuration is refused with the expected message |
| `exact_solution`, `property`, `metamorphic`, `differential`, `characterisation`, `conformance` | How the test knows the right answer; see [Testing and benchmarking](https://lsdtopotools.github.io/lsdtt3/developer_notes/testing_and_benchmarking.html#evidence-and-level-labels) |
| `unit` | Calls lsdtt3 routines in process (C++ drivers), or checks repository files |
| `system` | Judges outputs written by a shipped program run end to end (a CLI, the configure, the Python package) |
| `cosmo` | Cosmogenic-nuclide physics against CRONUS and pinned numerics sentinels |
| `flow` | Flow routing and depression handling: properties, and RichDEM comparison |
| `fetch` | Native DEM fetching (`lsdtt_raster_fetch`); offline logic by default |
| `performance` | In-process timings of the flow routines (opt-in; reported, not asserted) |
| `benchmark` | Opt-in comparison against external references (RichDEM/TopoToolbox flow, DAGGER hydraulic) on fetched real DEMs |

```sh
ctest --test-dir build -L '^smoke$' --output-on-failure
ctest --test-dir build -L '^unit$'  --output-on-failure
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
OMP_NUM_THREADS=8 ctest --test-dir build-perf -L '^performance$' --output-on-failure
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
thread-scaling caveats are in [Benchmarking and testing in lsdtt3](https://lsdtopotools.github.io/lsdtt3/developer_notes/testing_and_benchmarking.html) §6.

**No DEM is committed to this repository.** Tests that need real topography fetch
it on demand from public keyless sources and cache it under
`tests/data/dem_cache/` (gitignored), using `lsdtt-fetch-raster` itself.

See `tests/README.md` for what each test does, and
[Benchmarking and testing in lsdtt3](https://lsdtopotools.github.io/lsdtt3/developer_notes/testing_and_benchmarking.html) for methodology, benchmark results, and the
places where lsdtt3 and RichDEM deliberately differ.

---

## Command-line programs

All programs take a plain-text parameter file as their argument. Parameter files use `key: value` syntax, one entry per line; lines beginning with `#` are comments. Pass `-h` to any program for a full parameter reference.


```sh
lsdtt-raster-preprocessing  my_params.driver
lsdtt-drainage-analysis     my_params.driver
```

**Detailed parameter references for each program are in `docs/`.**

> **Command rename:** `lsdtt-synthetic-surfaces` replaces
> `lsdtt-raster-synthetic-tests`. The old executable is not retained; update
> scripts and reinstall lsdtt3 when upgrading.

### DEM preparation

| Program | Description | Reference |
|---|---|---|
| `lsdtt-fetch-raster` | Downloads a DEM for a WGS84 bounding box and writes it reprojected to UTM, ready to analyse. The default source is Copernicus GLO-30 on AWS Open Data — public and keyless, so no account or API key is needed. The area can be given as `west`/`south`/`east`/`north`, or as a vector file (`.fgb`, `.gpkg`) containing a polygon or two corner points — reprojected from UTM automatically, so a box drawn in a GIS works directly. Bounding boxes spanning several source tiles are mosaicked automatically. OpenTopography is also supported for its other datasets; that needs a free API key, supplied via the `OT_API_KEY` environment variable (or `OT_API_KEY_FILE`) — **never** in a parameter file, since parameter files are commonly committed and every parameter is echoed at the end of a run. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_fetch_raster.html) |
| `lsdtt-raster-preprocessing` | Reprojects geographic/Web Mercator DEMs to UTM/polar, applies elevation threshold masking, fills interior nodata voids, clips to maximum valid rectangle, breaches depressions (Lindsay 2016/2020), fills pits (Barnes 2014), writes hillshade. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_raster_preprocessing.html) |
| `lsdtt-synthetic-surfaces` | Generates one-shot analytical, stochastic, and ground-truth fixture surfaces for algorithm validation. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_synthetic_surfaces.html) |
| `lsdtt-synthetic-networks` | Optimal channel networks with tunable branching, for testing network metrics against known structure. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_synthetic_networks.html) |
| `lsdtt-fetch-lithology` | Downloads geological map units for a bounding box from Macrostrat. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_fetch_lithology.html) |
| `lsdtt-fetch-erosion-rates` | Downloads basin-averaged cosmogenic erosion rates for a bounding box from OCTOPUS. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_fetch_erosion_rates.html) |
| `lsdtt-lithology-erodibility` | Turns a lithology map and a property table into an erodibility (K) raster. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_lithology_erodibility.html) |

### Surface analysis

| Program | Description | Reference |
|---|---|---|
| `lsdtt-surface-metrics` | Hillshade, local relief, polynomial surface fit (direct and wavelet), curvature metrics, MADk2 roughness, Frangi filter, Otsu segmentation, Hough line extraction, curvature skeletons, cloth-drape upper envelope and incision depth. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_surface_metrics.html) |
| `lsdtt-flat-surface-finder` | Identifies low-texture surfaces (alluvial fans, terraces, piedmonts) using MADk2, residual IQR, and normal-vector dispersion (NVD). | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_flat_surface_finder.html) |
| `lsdtt-mountain-front` | Identifies mountain fronts from multi-scale relief, piedmont surface fitting, and residual thresholding. Supports synthetic DEM validation mode. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_mountain_front.html) |
| `lsdtt-feature-finder` | Feature-detection workflows (terraces, piedmonts, hillslope morphology). | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_feature_finder.html) |
| `lsdtt-spectral-analysis` *(requires FFTW3)* | Spectral surface analysis using the FFTW double API. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_spectral_analysis.html) |
| `lsdtt-basin-surface-metrics` | Multiscale basin divide and slope-area morphometrics. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_basin_surface_metrics.html) |
| `lsdtt-slope-area` | Polyfit slope against logarithmic drainage area. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_slope_area.html) |
| `lsdtt-swath-profile` | Percentiles of a raster binned along a baseline polyline. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_swath_profile.html) |
| `lsdtt-ridgeline-analysis` | Drainage divides and their ridgelines. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_ridgeline_analysis.html) |
| `lsdtt-fluvial-surfaces` | Sorts flat ground into fluvial levels (floodplains, terrace treads, fan lobes). | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_fluvial_surfaces.html) |

### Flow routing and drainage

| Program | Description | Reference |
|---|---|---|
| `lsdtt-drainage-analysis` | Drainage area via D8, D-infinity, or MFD routing. Topographic Wetness Index (TWI). | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_drainage_analysis.html) |
| `lsdtt-channel-extraction` | D8 flow routing, contributing-pixel threshold source extraction, junction/channel network construction, network export to CSV / GeoPackage / FlatGeobuf. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_channel_extraction.html) |
| `lsdtt-basin-extraction` | D8 flow routing, channel/junction network, junction-aware basin extraction. Writes basin outlines, raster labels, and statistics. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_basin_extraction.html) |
| `lsdtt-outlet-matching` | Matches outlet coordinates from a reference DEM to a target/resampled DEM using spatial overlap of complete D8 upstream masks; independent of channelisation thresholds. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_outlet_matching.html) |
| `lsdtt-basin-stats` | Hydromorphic statistics for outlet-point-defined basins: elevation, slope, flow-distance percentiles; centroid; drainage area; channel length by Strahler order. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_basin_stats.html) |
| `lsdtt-graphflood` | Graph-based hydraulic routing (Gailleton et al. 2024) on a fixed D-infinity or MFD flow graph. Writes water depth, water surface, and flux rasters. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_graphflood.html) |
| `lsdtt-stream-burning` | Loads a polyline channel network, seeds D8 routing from its source tips, and writes the DEM-extracted network for comparison before burning. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_stream_burning.html) |
| `lsdtt-internal-drainage` | Detects and routes closed (internally drained) basins. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_internal_drainage.html) |
| `lsdtt-distributed-hydrology` | A raster TOPMODEL: rainfall and PET drive stores on every pixel, and runoff is routed down the flow network. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_distributed_hydrology.html) |
| `lsdtt-water-centreline` | Centrelines for water-area polygons and their tributaries. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_water_centreline.html) |
| `lsdtt-water-line-repair` | Reconnects water lines truncated by roads. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_water_line_repair.html) |

### Channel network analysis

| Program | Description | Reference |
|---|---|---|
| `lsdtt-channel-network-analysis` | Channel-link bearing fitting, junction-angle calculation (donor opening, bending angle), Strahler/Tokunaga self-similarity analysis. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_channel_network_analysis.html) |
| `lsdtt-chi-analysis` *(requires Eigen3)* | End-to-end chi workflow: chi extraction, piecewise-linear network fitting, disorder metrics with jackknife/bootstrap uncertainty, spatial k_sn / u* inversion with m/n model selection. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_chi_analysis.html) |
| `lsdtt-cross-sections` | Cross-section extraction along channel networks or user-defined transects. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_cross_sections.html) |
| `lsdtt-valley-centreline` | Each channel's valley and its centreline. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_valley_centreline.html) |
| `lsdtt-bankfull-width` | Manning's-equation bankfull channel width, solved per cross section from a high-resolution DEM. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_bankfull_width.html) |
| `lsdtt-divide-migration` | Drainage-divide migration metrics, such as chi-coordinate offsets across divides. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_divide_migration.html) |
| `lsdtt-hillslope-channel-coupling` | Links hilltop curvature to the channels below it. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_hillslope_channel_coupling.html) |

### Tidal / marsh

| Program | Description | Reference |
|---|---|---|
| `lsdtt-marsh-topography` | Marsh platform identification, tidal channel extraction, skeleton-based channel metrics. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_marsh_topography.html) |

### Cosmogenic nuclides

| Program | Description | Reference |
|---|---|---|
| `lsdtt-cosmo-tool` | Basin-averaged CRN erosion-rate inversion. Reads measured nuclide concentrations, computes Stone 2000 production, snaps points to D8 network, inverts via Newton-Raphson + bisection. Optional ERA40 / palaeomagnetic correction (requires NetCDF). | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_cosmo_tool.html) |
| `lsdtt-cosmo-prediction` | Predicts CRN concentrations across a DEM: steady-state, step-change transient, and erosion-weighted accumulated concentration. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_cosmo_prediction.html) |
| `lsdtt-cosmo-column-tool` | 1-D transient cosmogenic nuclide column simulation. Steady-state, step-change, and full-transient (user-supplied erosion history) modes. No DEM required. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_cosmo_column_tool.html) |
| `lsdtt-cosmo-landslides` | Cosmogenic concentrations under triggered landslides with an imposed size distribution. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_cosmo_landslides.html) |

### Landscape evolution

| Program | Description | Reference |
|---|---|---|
| `lsdtt-lem` | Landscape evolution model: D-infinity flow routing, implicit stream-power fluvial incision (Anand et al. 2020 / Braun & Willett 2013), Cidre linear/nonlinear hillslope transport (Carretier et al. 2016), uniform uplift. Runs from a synthetic flat or noise surface, or an input DEM. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_lem.html) |
| `lsdtt-advection-diffusion` | A minimal transient landscape-evolution model. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_advection_diffusion.html) |
| `lsdtt-sediment-routing` | 1-D reach-scale bedload sediment routing, built from a port of d-CASCADE (Battista et al. 2022). | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_sediment_routing.html) |

### Utilities

| Program | Description | Reference |
|---|---|---|
| `lsdtt-linear-stats` | Statistical summaries for linear (directional/bearing) data. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_linear_stats.html) |
| `lsdtt-raster-calculator` | Two-raster algebra and overlap masking, for example DEMs of difference. | [reference](https://lsdtopotools.github.io/lsdtt3/cli_reference/cli_reference_lsdtt_raster_calculator.html) |

---

## Contributing

**Until the code base is stable, as decided by Simon M. Mudd, pull requests are
accepted only from repository collaborators.** Bug reports and feature requests
are welcome as [issues](../../issues). See [`CONTRIBUTING.md`](CONTRIBUTING.md).

## Other Documentation

The documentation is published at <https://lsdtopotools.github.io/lsdtt3/>.
Its source is the `docs/` tree in this repository, organised by audience:

| Directory | Audience |
|---|---|
| `docs/what_is_lsdtt3/` | **Users** — overview, installation, data, flow routing, benchmarks, and workflows. |
| `docs/design_and_algorithms/` | **Users and developers** — scientific design, algorithms, and numerical choices. |
| `docs/cli_reference/` | **Users** — one parameter reference per command-line tool. |
| `docs/python/`, `docs/python_api/` | **Users** — the Python package: worked guides, and the reference for every function. |
| `docs/developer_notes/` | **Developers** — architecture, building, packaging. |
| `docs/plan_docs/` | Working design documents. Not published, not committed. |
| `docs/tex_technical/` | LaTeX sources for method write-ups with full equations. |

### For users

| Document | Contents |
|---|---|
| [Installing lsdtt3](https://lsdtopotools.github.io/lsdtt3/what_is_lsdtt3/install_lsdtt3.html) | **Start here.** Installing a release package on Linux, Windows (WSL2 or native), macOS, or Google Colab. |
| [Overview](https://lsdtopotools.github.io/lsdtt3/what_is_lsdtt3/overview.html) | What lsdtt3 is, what it does, and how it differs from LSDTopoTools2. |
| [Your first analysis](https://lsdtopotools.github.io/lsdtt3/what_is_lsdtt3/first_analysis.html) | From a DEM on disk to a hillshade and a channel network. |
| [The lsdtt3 command-line tools](https://lsdtopotools.github.io/lsdtt3/what_is_lsdtt3/lsdtt3_clis.html) | All the command-line tools, grouped by what they are for, and how their parameter files work. |
| [Workflows](https://lsdtopotools.github.io/lsdtt3/what_is_lsdtt3/workflows.html) | Multi-stage, reproducible analysis chains run with Snakemake (`workflows/`). |
| [Raster Data](https://lsdtopotools.github.io/lsdtt3/what_is_lsdtt3/raster_data.html), `vector_data.md` | Input and output data formats. |
| [Running benchmarks on your own system](https://lsdtopotools.github.io/lsdtt3/what_is_lsdtt3/running_benchmarks.html) | Timing lsdtt3 against other tools on your own machine. |
| [Python bindings](https://lsdtopotools.github.io/lsdtt3/python/index.html) | The Python package: getting started and worked guides. |
| `docs/design_and_algorithms/` | Scientific design and algorithm documentation. |
| `docs/cli_reference/` | Parameter reference for every tool — see the table above for which file. |

### For developers

| Document | Contents |
|---|---|
| [lsdtt3 Architecture Document](https://lsdtopotools.github.io/lsdtt3/developer_notes/lsdtt3_architecture.html) | Full architecture reference: layer design, naming conventions, all module APIs. Start here for development. |
| [lsdtt3 Developer Onboarding](https://lsdtopotools.github.io/lsdtt3/developer_notes/developer_onboarding.html) | Getting-started guide for new contributors. |
| [Notes on git from the command line for new users](https://lsdtopotools.github.io/lsdtt3/developer_notes/git_notes.html) | Setting up an SSH key and cloning, for people new to git. |
| [Building lsdtt3 from source](https://lsdtopotools.github.io/lsdtt3/developer_notes/building_from_source.html) | Start here for compiling: the reference Linux build, CMake options, and which page covers your platform. |
| [Building on Windows](https://lsdtopotools.github.io/lsdtt3/developer_notes/building_on_windows.html) | Windows: WSL2 (full build) and native Visual Studio + vcpkg (QGIS backend). |
| [Building on macOS](https://lsdtopotools.github.io/lsdtt3/developer_notes/building_on_macos.html) | macOS: native Homebrew build, and a Multipass Ubuntu VM. |
| [QGIS plugin development](https://lsdtopotools.github.io/lsdtt3/developer_notes/qgis_plugin_development.html) | The QGIS plugin and its backend bundle: version checks and the development loop. |
| [Making release packages](https://lsdtopotools.github.io/lsdtt3/developer_notes/making_release_packages.html) | How the release packages are built and published. |
| [Building and publishing the documentation site](https://lsdtopotools.github.io/lsdtt3/developer_notes/building_the_docs.html) | Building the documentation site with MkDocs, and publishing it. |
| [CLI Parameter Naming Conventions](https://lsdtopotools.github.io/lsdtt3/developer_notes/cli_parameter_conventions.html) | Naming rules for CLI parameters. |
| [Basin Object Architecture](https://lsdtopotools.github.io/lsdtt3/developer_notes/basin_object_architecture.html) | Design boundaries of the basin object hierarchy. |
| [Python Bindings Architecture And API Contracts](https://lsdtopotools.github.io/lsdtt3/developer_notes/python_bindings_architecture.html) | The Python bindings: architecture and the public-API contracts (manifest, generated reference). |
| [Building the Python wheels](https://lsdtopotools.github.io/lsdtt3/developer_notes/building_python_wheels.html) | Building the manylinux and other Python wheels. |
| [Publishing lsdtt3 to PyPI](https://lsdtopotools.github.io/lsdtt3/developer_notes/publishing_to_pypi.html) | Publishing the wheels to PyPI, by hand or from GitHub Actions (not yet done). |
| [Backend Abstraction Rules for lsdtt3](https://lsdtopotools.github.io/lsdtt3/developer_notes/backend_abstraction_rules.html) | When to use abstract bases + factories vs. concrete CPU free functions. |
| [Using lsdtt_flowinfo](https://lsdtopotools.github.io/lsdtt3/developer_notes/using_lsdtt_flowinfo.html) | How to use the D8 flow routing object in analysis code. |
| [lsdtt3 utils contents](https://lsdtopotools.github.io/lsdtt3/developer_notes/lsdtt_utils_guide.html) | Guide to the `src/utils/` layer: what belongs there and the valid-data contract. |
| [GPU compute backend](https://lsdtopotools.github.io/lsdtt3/developer_notes/gpu_compute.html) | GPU backend design and the compute-dispatch mechanism. |
| [Beta testing and distribution notes](https://lsdtopotools.github.io/lsdtt3/developer_notes/beta_testing_notes.html) | Getting bundles to beta testers, the release pipeline, and the future conda-forge plan. |
| [Building a Debian package for lsdtt3](https://lsdtopotools.github.io/lsdtt3/developer_notes/debian_packaging.html) | Building a Debian package. |
| [Benchmarking and testing in lsdtt3](https://lsdtopotools.github.io/lsdtt3/developer_notes/testing_and_benchmarking.html) | How lsdtt3 is tested and benchmarked, and where it deliberately differs from RichDEM. Read before comparing lsdtt3 output against another tool. |
| `tools/README.md` | Building redistributable bundles, and the git/gh commands for cutting a tagged release. |
| `apptainer_def/apptainer_setup.md` | Apptainer container setup. |

---

## Python helper scripts

Separate from the `lsdtt3` Python package, we also provide various python helper scripts (in `python_helpers/`) for simple visualisation and downloading of DEMs. Maintaining python environments is very time consuming and annoying, which is why the core tools are written in c++. We provide helper functions for users familiar with maintaining a geospatial python environment, but if you don't know how to do this just use a GIS and the internet to get and look at `lsdtt3` data inputs and outputs.  

| Script | Purpose |
|---|---|
| `download_dem_tiles.py` | Download DEMs from Planetary Computer (keyless), OpenTopography, a WCS endpoint, or Copernicus Data Space. Note that we also have a command line interface, `lsdtt-fetch-raster`, that also does this. |
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
environment. [lsdviztools3](https://lsdtopotools.github.io/lsdviztools3/), our
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
