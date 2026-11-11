# Python bindings

`lsdtt3` is the Python interface to LSDTopoTools3. It keeps raster data and
analysis results in memory, so a script can compose terrain, flow, channel,
basin, chi, cosmogenic, and hydraulic operations without temporary parameter
files.

## Installation

The package is a wheel attached to each
[GitHub release](https://github.com/LSDtopotools/lsdtt3/releases/latest), for
Linux (x86_64), macOS (Apple Silicon and Intel) and Windows (x64), and Python
3.10 to 3.14. It is **not on PyPI yet**, so `pip install lsdtt3` does not work.
Install it into a virtual, conda or pixi environment, not your system Python.

**1. Find the wheel for your Python and system.** Run this with the Python you
will use; it prints the address of the matching wheel in the latest release:

```bash
python - <<'EOF'
import json, platform, sys, urllib.request
rel = json.load(urllib.request.urlopen(
    "https://api.github.com/repos/LSDtopotools/lsdtt3/releases/latest"))
py = f"-cp{sys.version_info.major}{sys.version_info.minor}-"
plat = {"linux": "manylinux", "darwin": "macosx", "win32": "win_"}[sys.platform]
arch = platform.machine().lower()      # x86_64, arm64 or amd64
for a in rel["assets"]:
    n = a["name"]
    if n.endswith(".whl") and py in n and plat in n and arch in n.lower():
        print(a["browser_download_url"])
EOF
```

(On Windows, save the lines between `<<'EOF'` and `EOF` as `find_wheel.py`
and run `python find_wheel.py`.)

**2. Install it** from the address it printed, for example:

```bash
python -m pip install https://github.com/LSDtopotools/lsdtt3/releases/download/v0.9.3/lsdtt3-0.9.3-cp313-cp313-manylinux_2_28_x86_64.whl
```

**3. Check it works:**

```bash
python -c "import lsdtt3; print(lsdtt3.__version__)"
```

If step 1 prints nothing, the latest release has no wheel for your Python
version or system; the [list of releases](https://github.com/LSDtopotools/lsdtt3/releases)
shows what each release has.

The wheel is self-contained: it bundles GDAL and PROJ and also installs the
command-line tools as commands, so you do not need a backend package. **In
Google Colab** the lookup and the install fit in one cell; it is in the
[installation guide](../what_is_lsdtt3/install_lsdtt3.md#python-package), and
the `python_binding_examples` notebook in the lsdtt3_notebooks repository
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
