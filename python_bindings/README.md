# lsdtt3 Python bindings

`lsdtt3` is a nanobind-based Python interface to LSDTopoTools3.

**Documentation:** guides at
<https://lsdtopotools.github.io/lsdtt3/python/index.html> and the reference
for every function at
<https://lsdtopotools.github.io/lsdtt3/python_api/index.html>. In Python,
`help(lsdtt3.<function>)` shows its docstring.

## Install

```bash
python -m pip install lsdtt3
```

Python 3.12 or later, on Linux (x86_64), macOS (Apple Silicon and Intel,
macOS 15 or later) and Windows (x64). Install into a virtual, conda or pixi
environment rather than your system Python. The wheels bundle GDAL and PROJ,
so nothing else needs installing, and they also install the lsdtt3
command-line tools (`lsdtt-*`) as commands. In Google Colab, run
`!pip install lsdtt3` in a cell.

There is one wheel per platform, `cp312-abi3`: it is built against Python's
stable ABI, so the same file installs on 3.12, 3.13, 3.14 and later. The
same wheels are attached to every
[GitHub release](https://github.com/LSDtopotools/lsdtt3/releases/latest),
together with wheels for older lsdtt3 versions that are not on PyPI; the
[installation guide](https://lsdtopotools.github.io/lsdtt3/what_is_lsdtt3/install_lsdtt3.html#python-package)
explains how to install one of those.

To install from a clone of the repository instead (needs a C++17 compiler,
CMake and the GDAL development files), from the repository root:

```bash
python -m pip install ./python_bindings
```

## Documentation in the package

Everything is documented in the installed package, so no repository access
is needed: `help(lsdtt3.fill_pits)`, `help(lsdtt3.FlowInfo)` and
`inspect.signature(...)` show parameters, units and return values. For a
compact overview written for you or for an AI assistant (Gemini in Colab,
for example), print the LLM helper:

```python
import lsdtt3
print(lsdtt3.helper_for_llm("all"))
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

flow = lsdtt3.FlowInfo(filled)
scope = flow.channel_scope(threshold_contributing_pixels=1000)
network = flow.channel_network(scope)
network.lines_geodataframe().to_file("channels.gpkg")
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

See [Python bindings architecture and API contracts](https://lsdtopotools.github.io/lsdtt3/developer_notes/python_bindings_architecture.html)
for the backend/API distinction, contract manifests, extension workflow, and
verification commands.

## Versioning

The Python package shares its version number with lsdtt3: the `lsdtt3` 0.7.1
wheels are built from the lsdtt3 `v0.7.1` release, the same source as the 0.7.1
command-line tools and backend bundles.

Versions up to 0.9.3 are not on PyPI; their wheels (one per Python version)
are attached to their GitHub releases. Some early releases have no wheels at
all, because the Python API was not supported for every version of the C++
core; that is expected, not a withdrawn release.

A version ending in `.postN` (for example `0.7.1.post1`) is the same source as
`0.7.1`, rebuilt or repackaged.

## Authors

Boris Gailleton — boris.gailleton@univ-rennes.fr
Simon Mudd
