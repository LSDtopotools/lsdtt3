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
`tools/make_all_packages.sh --only wheel` does all of this for you, for one
Python version, from a tag.

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

## The build matrix

Configured in `python_bindings/pyproject.toml` under `[tool.cibuildwheel]`:

- **Python versions:** CPython 3.10 through 3.14 (`cp310-*` … `cp314-*`).
- **Skipped:** PyPy, and musllinux.

The four platform legs, from `.github/workflows/python-wheels.yml`:

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
2. Go to Actions → *Build Python wheels* → Run workflow, and give the tag.
3. Tick *upload* to attach the wheels to the `vX.Y.Z` GitHub release (which
   `release.yml` created); otherwise they stay as workflow artefacts.

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
