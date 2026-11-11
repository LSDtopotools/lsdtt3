# Python bindings

`lsdtt3` is the Python interface to LSDTopoTools3. It keeps raster data and
analysis results in memory, so a script can compose terrain, flow, channel,
basin, chi, cosmogenic, and hydraulic operations without temporary parameter
files.

## Installation

Install it from PyPI into a virtual, conda or pixi environment, not your
system Python:

```bash
python -m pip install lsdtt3
python -c "import lsdtt3; print(lsdtt3.__version__)"
```

It needs Python 3.12 or later, on Linux (x86_64), macOS (Apple Silicon and
Intel) or Windows (x64). In Google Colab, run `!pip install lsdtt3` in a cell.
The same wheels, and those of older versions that are not on PyPI, are
attached to each [GitHub release](https://github.com/LSDtopotools/lsdtt3/releases/latest);
the [installation guide](../what_is_lsdtt3/install_lsdtt3.md#python-package)
shows how to install one of those.

The wheel is self-contained: it bundles GDAL and PROJ and also installs the
command-line tools as commands, so you do not need a backend package. The
`python_binding_examples` notebook in the lsdtt3_notebooks repository
installs the package and runs a complete basin, slope-area and chi workflow.
For development from a clone, see
[Building the Python wheels](../developer_notes/building_python_wheels.md).

## Start here

1. [Load, inspect, and save a DEM](getting_started.md).
2. [Fetch a DEM](fetching.md) when you do not already have one.
3. Build a [flow, channel, and basin workflow](flow_networks.md).
4. Continue with [chi and ksn](chi_ksn.md), channel-head extraction,
   [floodplains, terraces and valleys](fluvial_surfaces.md),
   [ridgeline profiles](ridgelines.md), or the cosmogenic
   and GraphFlood tools as needed.

## Core ideas

- `Raster` stores values, nodata mask, georeferencing, and CRS together.
- `FlowInfo` is the D8 foundation for channels, basins, and chi.
- `FlowGraph` provides distributed D-infinity or MFD routing products.
- Result objects retain related rasters and tables rather than writing files
  implicitly.
- The wheels carry their own GDAL and PROJ, with their data files, and
  `lsdtt3` points itself at them internally. It does not set `GDAL_DATA`,
  `PROJ_DATA` or `PROJ_LIB`, so other geospatial libraries in the same Python
  session (pyproj, rasterio, GMT) keep using their own data. (The 0.7.1 wheels
  did set them; if pyproj then fails to open its database, remove the three
  variables after `import lsdtt3`.)

The [API reference](../python_api/index.md) is generated from the installed
Python API. Guides explain a workflow; API pages specify the exact callable.
