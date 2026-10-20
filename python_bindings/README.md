# lsdtt3 Python bindings

`lsdtt3` is a nanobind-based Python interface to LSDTopoTools3.

## Install

```bash
python -m pip install lsdtt3
```

Wheels bundle GDAL and PROJ, and also install the lsdtt3 command-line tools as
commands. Install into a virtual, conda or pixi environment rather than your
system Python. See "Versioning" below for why some lsdtt3 versions have no
wheel.

To install from a clone of the repository instead (needs a C++17 compiler,
CMake and the GDAL development files), from the repository root:

```bash
python -m pip install ./python_bindings
```

## Quickstart

```python
import numpy as np
import lsdtt3

dem = lsdtt3.Raster.read("input.tif")
resampled = lsdtt3.reproject(dem, epsg=32612, cell_size=30.0)
filled = lsdtt3.fill_pits(resampled)

np.save("filled_dem.npy", filled.to_numpy())
filled.write("filled_dem.tif", format="gtiff")
```

## Structure

```text
python_bindings/
├── CMakeLists.txt          standalone runtime and package build
├── pyproject.toml          Python package metadata
├── bindings/               nanobind module and feature-specific C++ bindings
├── src/lsdtt3/             Python import package
├── examples/               small plotting workflows to check the framework
└── tests/                  synthetic binding tests
```

## Design

The binding exposes minimal data-owning `RasterMeta` and `Raster` classes for
grid metadata, values, and nodata masks. Analysis is deliberately functional:
operations accept rasters and return new rasters rather than wrapping CLI
parameter files or mutating inputs.

See [Python bindings architecture and API contracts](https://github.com/LSDtopotools/lsdtt3/blob/main/docs/developer_notes/python_bindings_architecture.md)
for the backend/API distinction, contract manifests, extension workflow, and
verification commands.

## Versioning

The Python package shares its version number with lsdtt3: `lsdtt3` 0.7.0 on
PyPI is built from the lsdtt3 `v0.7.0` release, the same source as the 0.7.0
command-line tools and backend bundles.

Wheels are not published for every lsdtt3 release, because the Python API is
not supported for every version of the C++ core. A version missing from PyPI is
expected; it is not a withdrawn release. `pip install lsdtt3` installs the
newest version that has wheels.

A version ending in `.postN` (for example `0.7.0.post1`) is the same source as
`0.7.0`, rebuilt or repackaged.

## Authors

Boris Gailleton — boris.gailleton@univ-rennes.fr
Simon Mudd
