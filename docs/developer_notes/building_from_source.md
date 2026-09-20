# Building lsdtt3 from source

You really only need to do this if you are a developer. If you only want to *use* lsdtt3, install a
release package instead — see [Installation](../what_is_lsdtt3/install_lsdtt3.md).
Build from source when you need code that has not been released yet, or when
you are modifying `lsdtt3`.

`lsdtt3` is a Linux-native C++17 codebase. **Linux is the default platform**:
every tool builds there, the automated builds are run there, and the test suite is run there. This
page is the reference build; the platform pages only cover what is different on
their platform, and link back here for everything else.

| Your machine | Route | Page |
|---|---|---|
| Linux | Native build — full tool set | this page, [Linux](#linux) |
| Windows | WSL2 (Ubuntu) — full tool set, if you have windows we suggest you use this | [Building on Windows](building_on_windows.md#wsl2-recommended) |
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
  OpenMP      : ON
  FFTW3       : ON
  GPU (CUDA)  : OFF
  Eigen3      : ON
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

If you don't want to modify your `.bashrc` file you can also just put the path to the binaries in your `$PATH`, but you need to do this every time you open a new shell:

```bash
export PATH=/path/to/lsdtt3/build:$PATH
lsdtt-channel-extraction --help
```

### Run the tests

Pick the row for what you want. **Configure** means add the options to the
`cmake -S . -B build` command above, then rebuild (`cmake --build build -j`);
**run** is the command afterwards. Options stay set in `build/` until you
change them.

| I want to… | Configure with | Also install | Run |
|---|---|---|---|
| Check the C++ build with nothing extra installed | *(defaults)* | — | `ctest --test-dir build -L '^hermetic$' --output-on-failure` |
| Run all the C++ tests | *(defaults)* | — | `ctest --test-dir build --output-on-failure` |
| Test the command-line tools as well (most CLI coverage) | `-DLSDTT_ENABLE_PYTHON_TESTS=ON` | `python3` (≥ 3.10), `gdal-bin` | `ctest --test-dir build --output-on-failure` |
| Run every offline test, including the geospatial Python ones | `-DLSDTT_ENABLE_PYTHON_TESTS=ON -DLSDTT_ENABLE_PYTHON_GEO_TESTS=ON` plus the Pixi interpreter ([below](#python-tests-from-pixi)) | Pixi `test` environment | `pixi run -e test ctest --test-dir build --output-on-failure` |
| Skip the slow tests | *(any)* | — | add `-LE '^heavy$'` to the run command |
| Test only the landscape evolution model | `-DLSDTT_ENABLE_PYTHON_TESTS=ON` (without it, only the LEM's C++ unit tests are registered) | `python3`, `gdal-bin` | `ctest --test-dir build -L '^lem_component$' --output-on-failure` |
| Test everything except the landscape evolution model | *(any)* | — | `ctest --test-dir build -LE '^lem_component$' --output-on-failure` |
| Test one part of the code, e.g. cosmogenic nuclides | `-DLSDTT_ENABLE_PYTHON_TESTS=ON` | `python3`, `gdal-bin` | `ctest --test-dir build -L '^cosmo$' -LE '^heavy$' --output-on-failure` |
| Run the tests that download data (e.g. the `cosmo_realworld` erosion-rate validations) | `-DLSDTT_ENABLE_PYTHON_TESTS=ON -DLSDTT_ENABLE_NETWORK_TESTS=ON` | `python3`, `gdal-bin`; a network connection | `ctest --test-dir build -L '^needs_network$' --output-on-failure` |
| Benchmark against RichDEM, TopoToolbox and Whitebox | `-DLSDTT_ENABLE_BENCHMARKS=ON` | Pixi, then `pixi install -e benchmark`; a network connection | `ctest --test-dir build -L '^benchmark$' --output-on-failure` |
| Time lsdtt3 on this machine | in a separate directory: `cmake -S . -B build-perf -DCMAKE_BUILD_TYPE=Release -DLSDTT_ENABLE_PERF_TESTS=ON` | — | `ctest --test-dir build-perf -L '^performance$' --output-on-failure` |
| Keep a full log, and send a report when tests fail | *(any)* | — | `cmake -P build/run_tests.cmake -LE '^heavy$' -j4`, then `cmake --build build --target test_report` |
| Build without any tests | `-DBUILD_TESTING=OFF` | — | — |

Labels are regular expressions, so keep the `^` and `$`: `-L network` would
also match `synthetic_networks`. The domain labels (`cosmo`, `network`,
`hillslope`, `raster`, …) are listed in `tests/README.md` in the repository.

#### What configure tells you

At the end of configure, CMake prints what it registered and why anything was
left out:

```
=== lsdtt3 test configuration (<n> tests) ===
  tier_native        : <n> tests (hermetic <n>)
  tier_python_stdlib : disabled (LSDTT_ENABLE_PYTHON_TESTS=OFF; <n> not registered)
  tier_python_geo    : disabled (LSDTT_ENABLE_PYTHON_GEO_TESTS=OFF; <n> not registered)
  ...
```

A default configure registers only the C++ tests. They do **not** exercise
most of the command-line tools, which are tested by Python scripts; a green
C++-only run says the libraries are sound, not that every CLI is. If an option
needs something you have not installed — a Python version, GDAL's
command-line programs, a Python package — configure stops with an error that
names it and says how to get it.

#### Python tests from Pixi

The four geospatial Python tests need the GDAL Python bindings, NumPy,
Rasterio, GeoPandas and Shapely. The simplest source is the repository's Pixi
`test` environment. **Configure outside Pixi** — so the C++ build keeps using
the system GDAL — and hand CMake Pixi's interpreter; then run the tests inside
the environment:

```bash
pixi install -e test
cmake -S . -B build -DLSDTT_ENABLE_PYTHON_TESTS=ON -DLSDTT_ENABLE_PYTHON_GEO_TESTS=ON \
      -DPython3_EXECUTABLE="$(pixi run -e test python -c 'import sys; print(sys.executable)')"
cmake --build build -j
pixi run -e test ctest --test-dir build --output-on-failure
```

`pixi install` on its own does not change the interpreter an existing build
tree was configured with; re-run `cmake` with `-DPython3_EXECUTABLE` to switch.
Running the tests inside Pixi is safe: the test launcher keeps Pixi's GDAL and
PROJ data away from lsdtt3 programs built against the system GDAL.

#### If tests fail

Run the suite through the logging wrapper, which keeps a complete log in
`build/test-logs/` even if the run is interrupted, then build a report:

```bash
cmake -P build/run_tests.cmake -LE '^heavy$' -j4
cmake --build build --target test_report
```

The report is one plain-text file in `build/test-logs/` describing your
machine, toolchain, dependency versions and build options, with the failing
tests' output. Your home directory appears as `~` and no environment variables
beyond a short fixed list are included — but read it before sending it with
your bug report.

See `tests/README.md` for what each test covers, and
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
| `BUILD_TESTING` | `ON` | Build the CTest suite (the native C++ tests; add the options below for more) |
| `LSDTT_ENABLE_PYTHON_TESTS` | `OFF` | Python-driven tests, most CLI coverage (needs Python ≥ 3.10 and GDAL's command-line programs) |
| `LSDTT_ENABLE_PYTHON_GEO_TESTS` | `OFF` | The Python tests that need geospatial packages (needs `LSDTT_ENABLE_PYTHON_TESTS`) |
| `LSDTT_ENABLE_BENCHMARKS` | `OFF` | External comparison benchmarks (needs network + pixi) |
| `LSDTT_ENABLE_NETWORK_TESTS` | `OFF` | Tests that download data, including `cosmo_realworld` |
| `LSDTT_ENABLE_PERF_TESTS` | `OFF` | Performance benchmarks |

### Examples

An everyday development build — optimised, with the native C++ tests:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

A debugging build in its own tree, so the release tree is left intact:

```bash
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug -j
```

A build with no optional dependencies — no Eigen3, no OpenMP. This is the
configuration to use on MSVC, and the one to try when a dependency problem
needs to be ruled out:

```bash
cmake -S . -B build -DLSDTT_EIGEN=OFF -DLSDTT_OPENMP=OFF
cmake --build build -j
```

The full test configuration, adding the Python and geospatial tiers to the
native suite (see [Run the tests](#run-the-tests) for the interpreter
requirements):

```bash
cmake -S . -B build -DLSDTT_ENABLE_PYTHON_TESTS=ON -DLSDTT_ENABLE_PYTHON_GEO_TESTS=ON \
      -DLSDTT_ENABLE_NETWORK_TESTS=ON
cmake --build build -j
ctest --test-dir build -LE '^heavy$' --output-on-failure
```

Check the test configuration block CMake prints at the end of configure (see
[What configure tells you](#what-configure-tells-you)) to confirm which tiers
were registered before running the suite.

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
