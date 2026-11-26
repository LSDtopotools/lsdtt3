# Command-line tools

**The command-line tools are lsdtt3.** Every analysis lsdtt3 can do is done by
one of its `lsdtt-*` programs. The [Python bindings](../python/index.md) and
the [QGIS plugin](../qgis.md) will feel more familiar if you already work in
Python or QGIS, but each contains only a subset of what the command-line tools
do. To install the tools, see
[Installing lsdtt3](../what_is_lsdtt3/install_lsdtt3.md): a release package,
or `pip install lsdtt3`, which installs every tool as well as the Python
bindings.

## How a tool is run

Each tool reads one plain-text **parameter file** and is run with that file as
its only argument:

```bash
lsdtt-raster-preprocessing my_parameters.param
```

A parameter file has one `key: value` per line; lines starting with `#` are
comments. Keys you leave out take their defaults. Most files name where to
read and write:

| Key | Meaning |
|---|---|
| `read_path`, `read_fname` | folder and file name of the input (give the extension, e.g. `my_dem.tif`) |
| `write_path`, `write_fname` | folder for the outputs, and the prefix of every output file name |

Then come the switches and settings of the analysis, listed on each tool's
page. Every run also writes `<write_fname>_README.html` / `.csv`, documenting
every parameter and the value used, and `<write_fname>_ingested_params.param`,
a parameter file that repeats the run exactly. `lsdtt-<tool> -h` lists the
parameters and `lsdtt-<tool> -v` prints the version.

## Example: a DEM of Glen Coe and its hillshade

**1. Fetch the DEM.** Save this as `fetch_glencoe.param`. The area is a
longitude and latitude box; the default source, Copernicus GLO-30, needs no
account or key, and the DEM is written reprojected to UTM at 30 m:

```
# Copernicus GLO-30 DEM of Glen Coe
write_path: ./
write_fname: glencoe

west:  -5.15
south: 56.62
east:  -4.90
north: 56.72
```

```bash
lsdtt-fetch-raster fetch_glencoe.param
```

This writes `glencoe_DEM.tif`.

**2. Make a hillshade.** Save this as `hillshade_glencoe.param`:

```
# Hillshade of the DEM from step 1
read_path: ./
read_fname: glencoe_DEM.tif
write_path: ./
write_fname: glencoe

write_hillshade_raster: true
```

```bash
lsdtt-raster-preprocessing hillshade_glencoe.param
```

This writes `glencoe_hillshade.tif`, along with the default outputs of the
preprocessing steps. Open it in QGIS or any GIS.

## The tools

The pages in this section document each tool: what it does, its parameters
and defaults, its outputs and worked examples. They are grouped in the menu by
task: DEM preparation, surface analysis, flow routing and drainage, channel
networks, and so on. The table of all tools, with a one-line description of
each, is in the
[README](https://github.com/LSDtopotools/lsdtt3#command-line-programs).
