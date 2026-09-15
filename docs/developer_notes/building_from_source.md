# Building lsdtt3 from source

This is the developer route. If you only want to *use* lsdtt3, install a
release package instead — see [Installation](../what_is_lsdtt3/install_lsdtt3.md).
Build from source when you need code that has not been released yet, or when
you are modifying lsdtt3.

lsdtt3 is a Linux-native C++17 codebase. **Linux is the reference platform**:
every tool builds there, CI builds there, and the test suite is run there. This
page is the reference build; the platform pages only cover what is different on
their platform, and link back here for everything else.

| Your machine | Route | Page |
|---|---|---|
| Linux | Native build — full tool set | this page, [Linux](#linux) |
| Windows | WSL2 (Ubuntu) — full tool set, **recommended** | [Building on Windows](building_on_windows.md#wsl2-recommended) |
| Windows, native `.exe` needed | Visual Studio + vcpkg — reduced, QGIS backend only | [Building on Windows](building_on_windows.md#native-windows-qgis-backend-only) |
| macOS | Native (Homebrew) or a Multipass Ubuntu VM | [Building on macOS](building_on_macos.md) |

Once you have a build, turning it into something other people can install is
covered in [Making release packages](making_release_packages.md), and the QGIS
plugin workflow in [QGIS plugin development](qgis_plugin_development.md).

---

## Linux

### Install dependencies

```bash
sudo apt update
sudo apt install -y \
    build-essential cmake pkg-config \
    libgdal-dev \
    libfftw3-dev \
    libeigen3-dev \
    libnetcdf-dev
```

| Package | Required? | Role |
|---|---|---|
| `build-essential` | **yes** | GCC C++17 compiler and standard build tools |
| `cmake` (≥ 3.15) | **yes** | Build system |
| `pkg-config` | **yes** | Lets CMake locate FFTW3 and others |
| `libgdal-dev` (≥ 3.0) | **yes** | Raster and vector I/O, CRS handling |
| `libfftw3-dev` | no | `lsdtt-spectral-analysis`, wavelet quadratic polyfit |
| `libeigen3-dev` | no | k_sn / u\* inversion in `lsdtt-chi-analysis` |
| `libnetcdf-dev` | no | ERA40 / palaeomagnetic corrections in the cosmo tools |

OpenMP ships with GCC and needs no separate install.

Optional dependencies are detected automatically. A tool that needs a missing
library is **skipped**, and the build still succeeds — so a partial dependency
set gives you a partial tool set rather than an error.

### Clone and build

If you have not cloned from GitHub before, [Git from the command
line](git_notes.md) walks through setting up an SSH key.

```bash
git clone git@github.com:lsdtopotools/lsdtt3.git
cd lsdtt3
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

CMake prints a summary of which optional libraries it found:

```
=== lsdtt3 build configuration ===
  Build type  : Release
  C++ std     : 17
  GDAL        : <version>
  OpenMP      : TRUE
  FFTW3       : TRUE
  GPU (CUDA)  : FALSE
  Eigen3      : TRUE
===================================
```

Read it — it is the quickest way to discover that, say, Eigen was missing and
chi inversion was silently dropped.

### Put the tools on your PATH

```bash
echo 'export PATH=/path/to/lsdtt3/build:$PATH' >> ~/.bashrc
source ~/.bashrc
lsdtt-channel-extraction --help
```

Replace `/path/to/lsdtt3` with wherever you cloned the repository.

### Run the tests

```bash
ctest --test-dir build --output-on-failure
```

The default suite needs no network, no downloaded data and no Python packages.
While developing, run one domain rather than the whole suite, e.g.
`ctest --test-dir build -L <domain> -LE heavy`. See `tests/README.md` in the
repository for what each test covers, and
[Testing and benchmarking](testing_and_benchmarking.md) for methodology.

---

## CMake options

These apply on every platform.

| Option | Default | Effect |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Debug` for development, `RelWithDebInfo` for profiling |
| `LSDTT_EIGEN` | `ON` | Eigen3-dependent routines (auto-off if Eigen3 is absent) |
| `LSDTT_OPENMP` | `ON` | OpenMP parallelism. **Must be `OFF` on MSVC** — see [Building on Windows](building_on_windows.md#why-openmp-must-be-off-on-msvc) |
| `LSDTT_GPU` | `OFF` | GPU compute backend (needs the CUDA toolkit) |
| `BUILD_TESTING` | `ON` | Build the CTest suite |
| `LSDTT_ENABLE_BENCHMARKS` | `OFF` | RichDEM comparison benchmarks (needs network + pixi) |
| `LSDTT_ENABLE_NETWORK_TESTS` | `OFF` | Tests that download data |
| `LSDTT_ENABLE_PERF_TESTS` | `OFF` | Flow-routing performance benchmarks |

## Build directory conventions

Keep one build directory per configuration. The packaging scripts look in these
by default, and mixing configurations in one tree silently ships the wrong
binaries:

| Build | Directory |
|---|---|
| Everyday Linux / WSL2 development | `build/` |
| Linux release (made by `tools/make_release.sh`) | `build/release/` |
| Native Windows QGIS backend | `build/windows-qgis/` |
| Native macOS | `build/macos-core/` |

---

## See also

- [Developer onboarding](developer_onboarding.md) — start here as a new contributor
- [Architecture](lsdtt3_architecture.md) — layer design, naming conventions, module APIs
- [Making release packages](making_release_packages.md) — turning a build into a distributable package
- [QGIS plugin development](qgis_plugin_development.md) — the plugin and its backend bundle
- [Debian packaging](debian_packaging.md) — building a `.deb`
