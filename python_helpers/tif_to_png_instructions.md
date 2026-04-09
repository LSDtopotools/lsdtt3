# tif_to_png.py — instructions

Converts GeoTIFF rasters to PNG images for quick visual inspection.  Uses GDAL for raster I/O and matplotlib for rendering.

---

## Dependencies

```
pip install gdal numpy matplotlib
```

`gdal` must match the GDAL system library version.  On a Debian/Ubuntu system with the UbuntuGIS PPA:

```bash
sudo apt install python3-gdal
```

Or in a conda environment:

```bash
conda install -c conda-forge gdal numpy matplotlib
```

---

## Usage

### Convert all tifs in a directory

```bash
python tif_to_png.py /path/to/outputs
```

Finds every `*.tif` in the directory and writes a matching `.png` alongside it.

### Convert a single named file

```bash
python tif_to_png.py /path/to/outputs --file my_dem_mountain_mask.tif
```

The `--file` argument can be a filename relative to the directory argument, or an absolute path:

```bash
python tif_to_png.py --file /absolute/path/to/my_raster.tif
```

### Run in the current directory

```bash
python tif_to_png.py
```

---

## Options

| Flag | Default | Description |
|---|---|---|
| `directory` | `.` | Directory to search for `.tif` files |
| `--file`, `-f` | — | Process only this single `.tif` file |
| `--colormap`, `-c` | `terrain` | Matplotlib colormap for single-band rasters |
| `--width` | `5.0` | Output image width in inches |
| `--dpi` | `300` | Output resolution in dots per inch |

### Useful colormaps for lsdtt3 outputs

| Raster type | Suggested colormap |
|---|---|
| DEM / elevation | `terrain`, `gist_earth` |
| Residual (positive + negative) | `RdBu_r`, `bwr` |
| Binary mask (0/1) | `grey` |
| Local relief | `hot_r`, `plasma` |
| Slope / curvature | `viridis`, `magma` |

---

## Output

Each `.tif` produces a `.png` with the same stem in the same directory.  Output images are **5 inches wide at 300 DPI** (1 500 pixels wide) by default, with height set to preserve the raster's aspect ratio.

Single-band rasters are stretched with a **2 %–98 % percentile linear stretch**.  Nodata pixels are rendered transparent (white when viewed on a white background).  Three-band rasters are treated as RGB.

---

## Example: inspect lsdtt-mountain-front outputs

After running `lsdtt-mountain-front` with default outputs in `/data/tabernas/`:

```bash
python tif_to_png.py /data/tabernas/
```

This produces:

```
tabernas_synthetic_dem.png
tabernas_synthetic_front.png
tabernas_low_ground_mask.png
tabernas_fitted_surface.png
tabernas_residual.png
tabernas_mountain_mask.png
```

To inspect just the residual with a diverging colormap:

```bash
python tif_to_png.py /data/tabernas/ --file tabernas_residual.tif --colormap RdBu_r
```
