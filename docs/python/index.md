# Python bindings

`lsdtt3` is the Python interface to LSDTopoTools3. It keeps raster data and
analysis results in memory, so a script can compose terrain, flow, channel,
basin, chi, cosmogenic, and hydraulic operations without temporary parameter
files.

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

## Installation

`lsdtt3` shares its version number with the lsdtt3 command-line tools, but
wheels are built only for releases where the Python API is supported, so not
every version has one. It is **not on PyPI yet** while the repository is in
private beta: when a release has wheels they are attached to its GitHub
release, and you install the file for your Python and platform into a virtual,
conda or pixi environment:

```bash
python -m pip install lsdtt3-0.7.1-cp312-cp312-manylinux_2_28_x86_64.whl
```

Each wheel is for one Python version (`cp312` is Python 3.12) and keeps its
original file name, which pip reads. In Google Colab, which runs Python 3.13 on
Linux, the `python_binding_examples` notebook in the lsdtt3_notebooks
repository installs the `cp313` wheel and runs a complete basin, slope-area and
chi workflow.

The [installation guide](../what_is_lsdtt3/install_lsdtt3.md#python-package)
has the details. For development from this checkout, see
[Building the Python wheels](../developer_notes/building_python_wheels.md).
