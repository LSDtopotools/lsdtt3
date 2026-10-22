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
file and the `before-all` command, so it is rebuilt only when one of those
changes, not at every release.

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
  package mirrors the image's `dnf` uses (first build only), and to PyPI
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

# 3. Build, repair and test one wheel per Python version, from the tag.
tools/make_all_packages.sh --clean --tag v0.7.1 --only wheel \
    --wheel-engine apptainer --wheel-python all
```

The wheels land in `dist/v0.7.1/`, five of them for CPython 3.10 to 3.14,
each with a `manylinux` platform tag in its name (`auditwheel` picks the
exact tag). Upload them as
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
- `dnf` cannot find `gdal-devel`: the package set is `before-all` in
  `pyproject.toml`, shared with CI; fix it there and both routes pick it up
  (the image hash changes, so the next run rebuilds it).
- A build or test failure after the image exists is a code or packaging
  failure, the same as it would be in CI: see
  [When a wheel build fails](#when-a-wheel-build-fails).

---

## The build matrix

Configured in `python_bindings/pyproject.toml` under `[tool.cibuildwheel]`:

- **Python versions:** CPython 3.10 through 3.14 (`cp310-*` … `cp314-*`).
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

System dependencies come from `before-all`: `dnf`/`yum` inside the manylinux
image, Homebrew on macOS, and vcpkg on Windows (with the installed tree cached
between runs).

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

The release is **manual by design** — the workflow is `workflow_dispatch` only
and never fires on a push.

1. Tag the release first (`vX.Y.Z`, or `vX.Y.Z.postN` for a wheel-only
   rebuild). The workflow builds an existing tag; it never creates one.
2. Go to Actions → *Build Python wheels* → Run workflow, give the tag, and
   choose **platforms** (see below).
3. Tick *upload* to attach the wheels to the `vX.Y.Z` GitHub release (which
   `release.yml` created); otherwise they stay as workflow artefacts.

### A Linux-only wheel build (the cheap one)

While the repository is private, Actions minutes are metered, and runners are
billed at a multiple of the Linux rate: **macOS at 10×, Windows at 2×**. Two of
the four legs are macOS, so a full run costs many times a Linux-only one. The
`platforms` input therefore defaults to `linux`, which builds only the Linux
x86_64 leg: five wheels (CPython 3.10–3.14), each tested as above, on one
`ubuntu-24.04` runner.

From the Actions tab: *Build Python wheels* → **Run workflow** → tag
`v0.7.1`, platforms **linux**, upload ticked or not. Or from a terminal with an
authenticated `gh`:

```bash
gh workflow run python-wheels.yml -f tag=v0.7.1 -f platforms=linux -f upload=true
gh run watch                       # follow it; or: gh run list --workflow python-wheels.yml
```

`gh workflow run` dispatches the workflow file as it is on the default branch,
and that version checks out and builds the tag you name.

**What a Linux-only release means for users.** There is no source
distribution, so `pip install lsdtt3` on macOS or native Windows finds no
wheel and fails with "no matching distribution". Linux users, including WSL2
and Google Colab, are covered. Choose `all` when a release should reach
macOS and Windows Python users. You can start with `linux` and run `all` for
the same tag later; with *upload* ticked, the second run replaces the Linux
wheels the first one attached.

The local alternative costs no minutes at all: `tools/make_all_packages.sh
--clean --tag v0.7.1 --only wheel` builds one Linux manylinux wheel, for the
interpreter it runs under, inside `podman` (see [Building locally](#building-locally)).

Before building anything the workflow checks the tag against the tagged tree:
the `X.Y.Z` stem must equal `CMakeLists.txt`, the tag minus its `v` must equal
`pyproject.toml`, and every version copy must agree. A mismatched tag fails in
seconds rather than after the matrix.

!!! note "The Python package shares the lsdtt3 version number"

    `python_bindings/pyproject.toml` carries the same `X.Y.Z` as
    `project(lsdtt3 VERSION ...)` in the top-level `CMakeLists.txt`
    (`tools/set_version.py` sets both). Wheels are built only for releases
    where the Python API is supported, so not every lsdtt3 version is on PyPI.
    Between wheel builds, `tools/make_release.sh` builds the bindings once on
    Linux and runs the whole `python_bindings/tests` suite at every release
    (also available as the opt-in CTest test `python_bindings_smoke`).

---

## When a wheel build fails

**A dependency is missing on one platform only.** Check `before-all` for that
platform in `pyproject.toml`. The three package managers do not agree on
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
