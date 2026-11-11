# Building the Python wheels

How the `lsdtt3` Python wheels are built, tested and released. For how the
bindings themselves are structured, see
[`python_bindings_architecture.md`](python_bindings_architecture.md).

---

## What a wheel contains

A wheel bundles the compiled extension module together with the shared
libraries it needs — GDAL, PROJ, FFTW, netCDF and the OpenMP runtime — so the
user needs none of them installed. That closure is produced by `auditwheel` on
Linux, `delocate` on macOS and `delvewheel` on Windows, all driven by
`cibuildwheel`.

Wheels are built with `-DLSDTT_PYTHON_RELEASE_WHEEL=ON`. That option exists to
make the release configuration explicit and repeatable rather than something
assembled by hand.

---

## Building locally

For development you do not want a wheel at all — an editable install is faster
and avoids the repair step:

```bash
python -m pip install -e python_bindings
```

Use a virtual environment or a conda/pixi environment for this; never install
into the global interpreter.

To build a real wheel on your own machine, the dependencies must be present
first. On Ubuntu:

```bash
sudo apt install cmake g++ pkg-config libgdal-dev libfftw3-dev \
                 libnetcdf-dev libeigen3-dev
python -m pip install --upgrade pip cibuildwheel
CIBW_CONTAINER_ENGINE=podman CIBW_MANYLINUX_X86_64_IMAGE=manylinux_2_28 \
    python -m cibuildwheel python_bindings --output-dir wheelhouse
```

On Linux `cibuildwheel` builds inside the manylinux container, and it uses
Docker unless told otherwise. This project does not use Docker: set
`CIBW_CONTAINER_ENGINE=podman` (`sudo apt install podman`).
`tools/make_all_packages.sh --only wheel` does all of this for you, from a tag.

If you have Apptainer but no podman (and no sudo to install it), see
[Building the manylinux wheels with Apptainer](#building-the-manylinux-wheels-with-apptainer).
`make_all_packages.sh` picks podman first, then Apptainer, then falls back to
a host wheel; `--wheel-engine podman|apptainer|host` forces one.

`cibuildwheel` will build one wheel per supported Python version, which takes
a while. To build a single wheel for the interpreter you already have, bypass
`cibuildwheel` and use `build` directly:

```bash
python -m pip install build
CMAKE_ARGS="-DLSDTT_PYTHON_RELEASE_WHEEL=ON" python -m build --wheel python_bindings
```

That wheel is **not** self-contained — it links whatever GDAL your machine has.
Fine for a local check, not for distribution.

---

## Building the manylinux wheels with Apptainer

cibuildwheel can only drive Docker or podman. On Linux, though, all it does is
run three steps inside the manylinux_2_28 image, and that image already holds
every CPython the wheels target and `auditwheel`. So
`tools/build_manylinux_wheel_apptainer.sh` runs the same three steps under
Apptainer:

1. **build**: `pip wheel python_bindings`, with `CMAKE_ARGS` from
   `[tool.cibuildwheel] environment`;
2. **repair**: `auditwheel repair`, which vendors GDAL, PROJ, FFTW and netCDF
   into `lsdtt3.libs/`;
3. **test**: a clean venv holding only the repaired wheel, then the
   `[tool.cibuildwheel] test-command` (the three checks in
   [What is tested before a wheel ships](#what-is-tested-before-a-wheel-ships)).

The Python versions, the system dependencies (`before-all`), the environment
and the test command are all read from `python_bindings/pyproject.toml`, so
this route and the cibuildwheel route (CI, or podman) build the same wheels.
The result is a real manylinux wheel, fit to upload.

The image is defined in `apptainer_def/manylinux_wheel.def`: the
`quay.io/pypa/manylinux_2_28_x86_64` image with `before-all` already run. The
script builds it once and reuses it. It is named after a hash of the `.def`
file, the `before-all` command and the script that command runs, so it is
rebuilt only when one of those changes, not at every release. The first build
takes 10 to 20 minutes, because `before-all` compiles PROJ, GDAL, GEOS, SQLite and curl
(see [The Linux wheel's GDAL and PROJ](#the-linux-wheels-gdal-and-proj)).

The podman route is unchanged: on a machine with podman,
`make_all_packages.sh` still uses cibuildwheel exactly as before.

### What you need on the server

- **Apptainer 1.1 or newer** (`apptainer --version`). From 1.1 an ordinary
  user can run `apptainer build --fakeroot`, which is how the image is built.
  **No sudo is needed.** If your account has `/etc/subuid` and `/etc/subgid`
  entries Apptainer uses them; if not, it falls back to its unprivileged
  fakeroot mode, which needs user namespaces enabled on the server and the
  `fakeroot` program available. `apptainer build --fakeroot` reports which
  of these is missing.
- **Network access** to `quay.io` (the base image, first build only), to the
  package mirrors the image's `dnf` uses and to `www.sqlite.org`, `curl.se`, `download.osgeo.org`
  and `github.com` for the PROJ, GDAL, GEOS, SQLite and curl sources (first build
  only), and to PyPI
  (`scikit-build-core`, `nanobind` and `pytest`, at every build).
- **Disk space**: a few GB for the image (default
  `~/.cache/lsdtt3/manylinux_2_28-lsdtt3-<hash>.sif`, or `--sif <path>` /
  `$LSDTT_WHEEL_SIF`), the same again transiently in Apptainer's cache
  (`~/.apptainer/cache`; set `APPTAINER_CACHEDIR` to put it elsewhere, and
  `apptainer cache clean` once the image exists), and a staged copy of the
  source tree under `$TMPDIR` during the build.
- **git**, **bash**, `awk` and `sha256sum`. No Python, pixi or conda
  environment is needed on the host for this route: everything Python runs
  inside the image. The one host Python call is `make_all_packages.sh`'s
  version stamp, which works with any `python3`.

Everything runs as you. The build uses `apptainer exec --cleanenv`, so an
activated pixi or conda shell's `GDAL_DATA`, `PROJ_LIB` or `PYTHONPATH` cannot
leak into it.

### The steps

```bash
git clone <lsdtt3 repo> && cd lsdtt3      # or git fetch --tags in an existing clone

# 1. Check: apptainer found, image path, Python versions.
tools/make_all_packages.sh --clean --tag v0.7.1 --only wheel \
    --wheel-engine apptainer --wheel-python all --check-only

# 2. Build the image once (optional: step 3 does it if it is missing).
#    Runs the image's %test section afterwards.
tools/build_manylinux_wheel_apptainer.sh --image-only
apptainer cache clean                      # reclaim the cache once the .sif exists

# 3. Build and repair the abi3 wheel, and test it with every Python, from the tag.
tools/make_all_packages.sh --clean --tag v0.7.1 --only wheel \
    --wheel-engine apptainer --wheel-python all
```

The wheel lands in `dist/v0.7.1/`: one `cp312-abi3` wheel, with a
`manylinux` platform tag in its name (`auditwheel` picks the exact tag). The
first Python builds it; the others install and test the same file. Upload them as
described in [`publishing_to_pypi.md`](publishing_to_pypi.md).

`--clean --tag` builds exactly what the tag contains, but the wheel script
and image definition are taken from your checkout, so the route also works
for tags that predate it. While iterating, `--wheel-python 3.12` builds a
single wheel, and the script can be run on its own against the working tree:

```bash
tools/build_manylinux_wheel_apptainer.sh --python 3.12     # wheels in ./wheelhouse
tools/build_manylinux_wheel_apptainer.sh --help
```

### If the image build fails

- `apptainer build --fakeroot` refused: your Apptainer is older than 1.1, or
  the administrator has disabled user namespaces. Ask for Apptainer >= 1.1, or
  build the `.sif` on another machine and copy it over with `--sif <path>`
  (the script's error message gives the exact command).
- `faked: ... GLIBC_2.33 not found` (followed by `fakeroot: error while
  starting the faked daemon`): your account has no `/etc/subuid` entry, so
  Apptainer runs `%post` under the host's fakeroot helper. That helper needs a
  newer glibc than the manylinux_2_28 image has. The script now detects the
  missing entry and builds with `--ignore-fakeroot-command`, and the `.def`
  runs `before-all` under the image's own fakeroot. If you still see this
  error, your copy of the script is older than that fix: pull and rerun.
- `no space left on device` during `%post`: the build unpacks the image and
  installs the packages in `APPTAINER_TMPDIR` (default `/tmp`), which needs
  about 3 GB. The script warns if there is less than 4 GB free. Set
  `APPTAINER_TMPDIR` and `APPTAINER_CACHEDIR` to directories on local disk.
  NFS often fails for fakeroot builds.
- `cpio: chown` or `unpacking of archive failed` from `dnf`: the build ran
  without `/etc/subuid` entries and the image's own `fakeroot` could not be
  installed. Ask for subuid/subgid entries, or build the `.sif` on a machine
  that has them and pass it with `--sif`.
- Once the image exists, `pip` fails to reach PyPI: the build runs with
  `--cleanenv`. The script passes the proxy variables, `PIP_INDEX_URL`,
  `PIP_EXTRA_INDEX_URL`, `PIP_TRUSTED_HOST` and the certificate variables
  through, so set those in the shell you run it from.
- `c++: fatal error: Killed signal terminated program cc1plus`: the compiler
  ran out of memory. Lower `--jobs`.
- A failure in `build_manylinux_native_deps.sh` (a download, a checksum, a
  compile, or its closing driver check): that script is `before-all`, shared
  with CI. Fix it there and every route picks it up. The image hash changes,
  so the next run rebuilds the image.
- A build or test failure after the image exists is a code or packaging
  failure, the same as it would be in CI: see
  [When a wheel build fails](#when-a-wheel-build-fails).

---

## The build matrix

Configured in `python_bindings/pyproject.toml` under `[tool.cibuildwheel]`:

- **Python versions:** one wheel per platform, `cp312-abi3`, for CPython 3.12
  and later. The extension is built against Python's stable ABI (nanobind
  `STABLE_ABI` in `python_bindings/bindings/CMakeLists.txt`, and
  `wheel.py-api = "cp312"` under `[tool.scikit-build]`). cibuildwheel builds it
  once, on 3.12, then installs and tests that same wheel on 3.13 and 3.14
  (`build = "cp312-* cp313-* cp314-*"`). Binding code may call only
  limited-API C functions; anything else fails to compile. Free-threaded
  CPython (`cp313t`, `cp314t`) cannot load abi3 wheels and is not built.
- **Skipped:** PyPy, and musllinux.

The four platform legs, from `.github/workflows/python-wheels.yml`. A run builds
either the Linux leg alone or all four (the workflow's `platforms` input; see
[Releasing](#releasing)):

| Leg | Runner | Arch | Notes |
|---|---|---|---|
| Linux x86_64 | `ubuntu-24.04` | `x86_64` | `manylinux_2_28` image |
| macOS ARM64 | `macos-15` | `arm64` | `MACOSX_DEPLOYMENT_TARGET=15.0` |
| macOS Intel | `macos-15-intel` | `x86_64` | `MACOSX_DEPLOYMENT_TARGET=15.0` |
| Windows x86_64 | `windows-2025` | `AMD64` | vcpkg toolchain; **OpenMP off** |

System dependencies come from `before-all`: inside the manylinux image,
`python_bindings/tools/build_manylinux_native_deps.sh`; Homebrew on macOS; and
vcpkg on Windows (with the installed tree cached between runs).

### The Linux wheel's GDAL and PROJ

The Linux wheel bundles whatever GDAL and PROJ the manylinux image has. In
AlmaLinux 8 (the base of manylinux_2_28), EPEL 8 has only GDAL 3.0.4 and
PROJ 6.3. PROJ before 7.2 cannot convert a CRS that uses a datum ensemble,
which is how current pyproj and GDAL describe WGS 84 and most CRSs based on
it. With those versions, reprojection fails with `DatumEnsemble can only be
exported to WKT2:2019`. So `build_manylinux_native_deps.sh` compiles the
following from pinned, checksummed source tarballs, using the image's own
compiler (gcc-toolset):

- PROJ 9.5.1
- GDAL 3.10.3
- GEOS 3.13.1. GDAL needs 3.8 or newer.
- SQLite 3.47.2. GDAL needs 3.31 or newer, and AlmaLinux 8 has 3.26.
- curl 8.11.1. GDAL needs curl 7.68 or newer for `/vsicurl`, and AlmaLinux 8
  has 7.61.

Using the image's compiler keeps the result manylinux-compliant. Libraries
from conda-forge are not, because they need a newer libstdc++.

GDAL is built with `GDAL_USE_EXTERNAL_LIBS=OFF`, so it uses only the
libraries the script names (those above, zlib, and glibc's iconv) plus GDAL's
own internal copies of libtiff, libgeotiff, libpng, libjpeg, json-c and a few
others. Every GDAL driver that needs nothing more is built. Drivers that need
another library are not, for example netCDF, HDF5, JPEG 2000, WebP and Parquet.
At the end, the script checks that the formats lsdtt3 reads and writes are
present: GTiff/COG, MEM, VRT, PNG, JPEG, ENVI/EHdr, AAIGrid and HFA for
rasters, and Shapefile, GeoJSON, GPKG, FlatGeobuf, MVT and CSV for vectors.
It fails if any of them is missing. If lsdtt3 starts depending on another
format, add it to that check.

The script also has a local-test mode. Run it as an ordinary user with
`LSDTT_NATIVE_DEPS_PREFIX` set to a writable directory. It then skips `dnf`
and `ldconfig` and runs the downloads, the builds and the closing checks on
the host. That is a quick way to check a version bump before rebuilding the
image.

This adds 10 to 20 minutes to every run of the CI Linux leg, and to the first
Apptainer image build. The Apptainer image is then reused. The macOS and
Windows legs are unchanged.

The bundled curl only knows AlmaLinux's CA-certificate path. On Debian,
Ubuntu (including Colab) and other distributions, `lsdtt3/__init__` sets
`CURL_CA_BUNDLE` to the system certificate bundle at import, unless
`CURL_CA_BUNDLE` or `SSL_CERT_FILE` is already set. Without that, every
`https://` read through `/vsicurl` would fail certificate verification.

!!! warning "Windows wheels are single-threaded"

    The Windows leg builds with `-DLSDTT_PYTHON_OPENMP=OFF`. MSVC's OpenMP
    support does not cover the constructs lsdtt3 uses — see the C3016/C4849
    discussion in [`building_on_windows.md`](building_on_windows.md). Windows
    users therefore get correct results at single-thread speed. Do not quote a
    Windows wheel timing as an lsdtt3 benchmark.

---

## What is tested before a wheel ships

`cibuildwheel` runs `test-command` against every wheel it builds, in a fresh
environment with only the wheel installed. Three checks:

1. `pytest -q python_bindings/tests` — the binding test suite.
2. `python_bindings/tools/check_api_examples.py --run` — executes the examples
   in the API documentation, so a doc example that no longer works fails the
   release.
3. `python_bindings/tools/check_installed_wheel.py --manifest
   python_bindings/api/public_api.json --smoke-console-scripts` — checks the
   installed wheel actually exposes the declared public API, and that the
   console scripts start.

The third is the one that catches packaging faults rather than code faults: a
wheel can import cleanly and still be missing a symbol or a data file.

---

## Releasing

**Every `vX.Y.Z` tag builds the wheels automatically.** `release.yml` calls
`python-wheels.yml` (as a reusable workflow) once the Linux backend has passed
and created the release: all four platform legs, one abi3 wheel each,
tested as above on CPython 3.12, 3.13 and 3.14, then attached to the release. It
never fires on an ordinary push.

Run it by hand for a `vX.Y.Z.postN` wheel-only rebuild (which `release.yml`
ignores), or to rebuild one tag's wheels:

1. Tag first (`vX.Y.Z.postN`). The workflow builds an existing tag; it never
   creates one.
2. Go to Actions → *Build Python wheels* → Run workflow, give the tag, and
   choose **platforms** (see below).
3. Tick *upload* to attach the wheels to the `vX.Y.Z` GitHub release (which
   `release.yml` created); otherwise they stay as workflow artefacts.

### Compiler cache (ccache) in CI

When each platform built five wheels, one per CPython version, about 250 of
the roughly 300 C++ files (the lsdtt3 runtime and the bundled CLIs) compiled
identically for every Python version, so `python-wheels.yml` runs the compiler
through ccache and wheels 2 to 5 reused those objects. Each platform now
builds a single abi3 wheel, so within a run there is nothing to reuse and
ccache saves little; it stays until CI timings show whether to remove it. Three settings
make that work, all in the workflow, not in `pyproject.toml`:

- `-DCMAKE_C_COMPILER_LAUNCHER` / `-DCMAKE_CXX_COMPILER_LAUNCHER`, passed in
  `CMAKE_ARGS` (only scikit-build-core reads it). Not the `CMAKE_<LANG>_COMPILER_LAUNCHER`
  environment variables: `CIBW_ENVIRONMENT` also reaches the Linux dependency
  build in before-all, which runs before ccache is installed;
- `SKBUILD_BUILD_DIR=<project>/python_bindings/build/cibw-{wheel_tag}`, a fixed
  per-wheel build folder instead of scikit-build-core's random temporary one;
- `CCACHE_BASEDIR=<project>`, so ccache compares paths relative to the project
  and the per-wheel folders look identical to it.

`python_bindings/tools/ci_ccache.sh` installs ccache before each build (EPEL in
the manylinux container, Homebrew on macOS) and prints its hit statistics, so
the log shows whether it is working. If ccache cannot be installed it puts a
pass-through stand-in on `PATH`, so the build never fails because of it.

**Windows** builds with Ninja instead of cibuildwheel's default Visual Studio
generator, which ignores compiler launchers and compiles the files of one
project one at a time. Ninja needs the compiler on `PATH`, so the job first
sets up an x64 Visual Studio developer environment (`ilammy/msvc-dev-cmd`,
as the Windows backend build does). The "Install ccache (Windows)" step
downloads a pinned ccache release and sets `LSDTT_WINDOWS_LAUNCHER` to the
launcher flags, which the Windows `CMAKE_ARGS` expands; if the download fails
the variable is empty and the wheels build uncached. ccache handles MSVC only
without `/Zi` debug information, which a Release build does not use. The
statistics are printed in a step of their own after the build.

The cache lives only for the length of the job; nothing is saved to the
Actions cache. Local builds and the Apptainer route do not use ccache.

### A Linux-only wheel build

The `platforms` input of a manual run defaults to `linux`, which builds only
the Linux x86_64 leg: one wheel, on one `ubuntu-24.04` runner, in a fraction
of the time of the full matrix. Useful for a quick `.postN` fix that only
Linux users need.

From the Actions tab: *Build Python wheels* → **Run workflow** → tag
`v0.7.1`, platforms **linux**, upload ticked or not. Or from a terminal with an
authenticated `gh`:

```bash
gh workflow run python-wheels.yml -f tag=v0.7.1 -f platforms=linux -f upload=true
gh run watch                       # follow it; or: gh run list --workflow python-wheels.yml
```

`gh workflow run` dispatches the workflow file as it is on the default branch,
and that version checks out and builds the tag you name.

**What a Linux-only run means for users.** There is no source
distribution, so `pip install lsdtt3` on macOS or native Windows finds no
wheel and fails with "no matching distribution". Linux users, including WSL2
and Google Colab, are covered. With *upload* ticked, a later `all` run for the
same tag replaces the Linux wheels an earlier run attached.

The local alternative needs no CI at all: `tools/make_all_packages.sh
--clean --tag v0.7.1 --only wheel` builds one Linux manylinux wheel, for the
interpreter it runs under, inside `podman` (see [Building locally](#building-locally)).

Before building anything the workflow checks the tag against the tagged tree:
the `X.Y.Z` stem must equal `CMakeLists.txt`, the tag minus its `v` must equal
`pyproject.toml`, and every version copy must agree. A mismatched tag fails in
seconds rather than after the matrix.

!!! note "The Python package shares the lsdtt3 version number"

    `python_bindings/pyproject.toml` carries the same `X.Y.Z` as
    `project(lsdtt3 VERSION ...)` in the top-level `CMakeLists.txt`
    (`tools/set_version.py` sets both). Release tags build wheels in CI
    (the first to attach them was v0.9.2); before that they were built by hand for some releases only.
    `tools/make_release.sh` also builds the bindings once on Linux and runs the
    whole `python_bindings/tests` suite in the first release job, so a C++
    change that breaks them stops the release before the wheel matrix starts
    (also available as the opt-in CTest test `python_bindings_smoke`).

---

## When a wheel build fails

**A dependency is missing on one platform only.** Check `before-all` for that
platform in `pyproject.toml`. On macOS, FFTW is found through pkg-config (so
`pkgconf` must be installed) and OpenMP through
`cmake/lsdtt_homebrew_openmp.cmake`, because Homebrew's `libomp` is keg-only;
a mac wheel built without either is single-threaded or lacks the spectral
routines, and the configure log says so. The three package managers do not agree on
package names, and a dependency added for Linux is easy to forget on the other
two.

**The repair step cannot find a library.** On Windows this is nearly always the
`--add-path` in `CIBW_REPAIR_WHEEL_COMMAND_WINDOWS`, which points at the vcpkg
bin directory. On Linux and macOS the repair tools follow the linker, so a
missing library usually means it was linked from outside the expected prefix.

**`check_installed_wheel.py` fails but tests pass.** The code is fine and the
packaging is not: something in `public_api.json` is not reaching the installed
wheel. Check the manifest matches what the bindings actually export.

**It works locally but not in CI.** The manylinux image is much older than your
machine. A local build linking your system GDAL proves nothing about the
manylinux leg.

---

## See also

- [`python_bindings_architecture.md`](python_bindings_architecture.md) — how the
  bindings are laid out and what the checker scripts do
- [`publishing_to_pypi.md`](publishing_to_pypi.md) — getting the wheels to
  users, by hand with twine or from CI
- [`making_release_packages.md`](making_release_packages.md) — the CLI backend
  tarballs, which are a separate release path
- [`building_on_windows.md`](building_on_windows.md) — the MSVC OpenMP situation
