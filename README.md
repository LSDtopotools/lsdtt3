# lsdtt3

A redesign of [LSDTopoTools2](https://github.com/LSDtopotools/LSDTopoTools2) with a modern C++ architecture.

**Copyright:** Simon M. Mudd
**Licence:** GPL-3.0

---

## Overview

lsdtt3 is a C++ library and set of command line tools for topographic analysis
of digital elevation models (DEMs). It is a ground-up redesign of
[LSDTopoTools2](https://github.com/LSDtopotools/LSDTopoTools2), retaining the
core algorithms and scientific functionality while rebuilding the architecture
for modern C++17, cleaner GPU support, and easier incremental development.

If you are using this software for research, please also cite the original
LSDTopoTools2:

> Mudd, S. M., Clubb, F. J., Grieve, S. W. D., Milodowski, D. T., Gailleton,
> B., Hurst, M. D., Valters, D. V., Wickert, A. D., & Hutton, E. W. H. (2023).
> LSDtopotools/LSDTopoTools2: LSDTopoTools2 v0.9 (v0.9). Zenodo.
> https://doi.org/10.5281/zenodo.8076231

lsdtt3 has the following design goals compared to LSDTopoTools2:

- Clean separation between file I/O, data storage, and computation
- CPU and GPU compute backends behind a shared abstract interface
- Analysis routines as free functions — changing one analysis file only
  recompiles that file
- D8 and DAG (D-infinity / MFD) flow routing behind a common interface
- Cloud Optimised GeoTIFF (COG) output by default

---

## Dependencies

| Dependency | Required | Notes |
|---|---|---|
| CMake | Yes | >= 3.15 |
| C++ compiler | Yes | GCC >= 9 or Clang >= 10, C++17 |
| GDAL | Yes | >= 3.0 recommended |
| FFTW3 | No | Double precision. Enables wavelet quadratic fit. |
| OpenMP | No | Enables multi-core support |
| CUDA toolkit | No | Required only for GPU backend |

### Installing dependencies on Ubuntu / Debian

```sh
sudo apt update
sudo apt install cmake build-essential libgdal-dev libfftw3-dev
```

OpenMP is included with GCC on Linux — no separate install needed.

### Installing with the Apptainer container (recommended)

A pre-configured Apptainer container definition is provided at `lsdtt3.def`.
It includes all dependencies (GDAL, FFTW3, Eigen3, OpenCV, OpenMP) on Ubuntu
22.04. See `docs/apptainer_setup.md` for setup instructions on Linux, macOS,
and Windows (WSL2).

### macOS notes

Apple's Clang does not ship with OpenMP. To enable multi-core support on macOS:

```sh
brew install libomp
cmake .. \
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include" \
  -DOpenMP_CXX_LIB_NAMES="omp" \
  -DOpenMP_omp_LIBRARY=$(brew --prefix libomp)/lib/libomp.dylib
```

---

## Building

```sh
git clone https://github.com/LSDtopotools/lsdtt3.git
cd lsdtt3
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

Compiled executables are placed in `bin/`.

### CMake options

| Option | Default | Description |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Release`, `Debug`, or `RelWithDebInfo` |
| `LSDTT_GPU` | `OFF` | Enable GPU compute backend (requires CUDA) |

---

## Command line programs

All programs take a parameter file as their argument. Parameter files are plain
text with `key: value` pairs, one per line. Lines beginning with `#` are
ignored. Pass `-h` to any program for parameter documentation.

```sh
./bin/lsdtt-raster-preprocessing my_params.driver
./bin/lsdtt-test-flowinfo my_params.driver
```

**`lsdtt-raster-preprocessing`** — conditions a DEM before analysis. Operations
include threshold masking, interior nodata filling, rectangular clipping,
depression breaching (Lindsay 2016, Lindsay 2020), pit filling (Barnes 2014
Priority-Flood), and hillshade.

**`lsdtt-test-flowinfo`** — builds D8 flow routing from a DEM and writes
diagnostic rasters: n_contributing_pixels (proportional to drainage area),
flow direction, donor count, and outlet mask.

**`lsdtt-raster-synthetic-tests`** — generates synthetic rasters and tests
the wavelet quadratic fit analysis.

---

## Code structure

### Raster I/O

All file reading and writing goes through `lsdtt_data_loader` in `src/core/`.
This is the only object that touches GDAL — compute objects never do file I/O
directly. Output is Cloud Optimised GeoTIFF by default, with plain GeoTIFF
and same-as-input format also available.

Raster data is held in `lsdtt_raster_cpu` in `src/cpu/`, which stores elevation
data as `lsdtt_array2d<float>` and a validity mask as `lsdtt_array2d<uint8_t>`.
These are simple wrappers over `std::vector` with 2D row/column indexing.

### Analysis free functions

All analysis algorithms live in `src/cpu/` as stateless free functions. They
take a `lsdtt_raster_cpu` reference and return a new `lsdtt_raster_cpu`. The
input is never modified. Functions are grouped by category:

- `lsdtt_raster_surface_metrics_cpu` — hillshade, moving window statistics
- `lsdtt_raster_flow_preprocessing_cpu` — pit filling, depression breaching, threshold masking
- `lsdtt_raster_nodata_filling_cpu` — interior void filling (iterative average or Laplacian)
- `lsdtt_raster_spectral_cpu` — wavelet quadratic fit (requires FFTW3)
- `lsdtt_raster_synthetic_surfaces_cpu` — synthetic test surface generators

### Flow routing

Flow routing is built around two complementary objects, both in `src/cpu/`.

**`lsdtt_flowinfo`** handles D8 routing — one receiver per cell. The flow
graph is a forest of trees, which allows the upslope set of any pixel to be
stored as a contiguous interval in a flat array (the s_vector). This makes
repeated upslope queries O(1): chi analysis, steepness, knickpoint detection,
and basin statistics all use this. See `docs/using_lsdtt_flowinfo.md` for the
full interface.

**`lsdtt_flowgraph_cpu`** (stub, pending implementation) handles D-infinity
and multiple-flow-direction (MFD) routing — fractional weighted flow to
multiple receivers. Use this for distributed-flow analyses such as the
topographic wetness index and sediment flux.

Both objects share the same pixel index system: a flat row-major index
(`flat_index_type`, currently `int32_t`) that converts trivially to and from
row/column coordinates.

### CPU / GPU architecture

Each major object has an abstract base class in `src/core/` and a CPU
implementation in `src/cpu/`:

```
src/core/lsdtt_raster_base   ←  abstract raster interface
src/cpu/lsdtt_raster_cpu     ←  CPU implementation

src/core/lsdtt_flow_base     ←  abstract D8 flow interface
src/cpu/lsdtt_flowinfo_cpu   ←  CPU implementation
```

Analysis code includes the CPU headers directly and works with the CPU types.
CLI programs use the factory functions (`make_flow_base()`) which return the
correct backend — currently CPU, and GPU once implemented — without the CLI
needing to change. A GPU backend for flow routing is added by implementing
`lsdtt_flow_base` and adding a branch to the factory; all existing CLIs and
analysis code work unchanged.

---

## Documentation

- `docs/lsdtt3_architecture.md` — full architecture reference for contributors
- `docs/using_lsdtt_flowinfo.md` — how to use the flow routing object in analysis code
- `docs/apptainer_setup.md` — container setup on Linux, macOS, and Windows

---

## References

Mudd, S. M., Clubb, F. J., Grieve, S. W. D., Milodowski, D. T., Gailleton,
B., Hurst, M. D., Valters, D. V., Wickert, A. D., & Hutton, E. W. H. (2023).
LSDtopotools/LSDTopoTools2: LSDTopoTools2 v0.9 (v0.9). Zenodo.
https://doi.org/10.5281/zenodo.8076231

Barnes, R., Lehman, C., Mulla, D. (2014) Priority-flood: An optimal
depression-filling and watershed-labeling algorithm for digital elevation
models. *Computers & Geosciences*, 62, 117–127.

Lindsay, J.B. (2016) Efficient hybrid breaching-filling sink removal methods
for flow path enforcement in digital elevation models. *Hydrological
Processes*, 30(6), 846–857.

Zevenbergen, L.W. and Thorne, C.R. (1987) Quantitative analysis of land
surface topography. *Earth Surface Processes and Landforms*, 12, 47–56.
