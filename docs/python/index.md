# Python bindings

`lsdtt3` is the Python interface to LSDTopoTools3. It keeps raster data and
analysis results in memory, so a script can compose terrain, flow, channel,
basin, chi, cosmogenic, and hydraulic operations without temporary parameter
files.

## Start here

1. [Load, inspect, and save a DEM](getting_started.md).
2. [Fetch a DEM](fetching.md) when you do not already have one.
3. Build a [flow, channel, and basin workflow](flow_networks.md).
4. Continue with [chi and ksn](chi_ksn.md), channel-head extraction, or the
   cosmogenic and GraphFlood tools as needed.

## Core ideas

- `Raster` stores values, nodata mask, georeferencing, and CRS together.
- `FlowInfo` is the D8 foundation for channels, basins, and chi.
- `FlowGraph` provides distributed D-infinity or MFD routing products.
- Result objects retain related rasters and tables rather than writing files
  implicitly.

The [API reference](../python_api/index.md) is generated from the installed
Python API. Guides explain a workflow; API pages specify the exact callable.

## Installation

Install a released wheel with:

```bash
python -m pip install lsdtt3
```

For development from this checkout, build/install `python_bindings/` in an
environment that has the native dependencies required by lsdtt3.
