# tif_to_png.py — instructions

Converts GeoTIFF rasters to PNG images for quick visual inspection.
Calls `gdalinfo` and `gdal_translate` from the GDAL command-line tools.
No Python GDAL bindings or numpy required — only the standard library plus
a working GDAL installation.

---

## Dependencies

Only the GDAL command-line tools are needed:

```bash
# Debian / Ubuntu
sudo apt install gdal-bin

# conda
conda install -c conda-forge gdal
```

Verify with:

```bash
gdal_translate --version
gdalinfo --version
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

`--file` can be a filename relative to the directory argument, or an absolute path:

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
| `--width` | `5.0` | Output image width in inches |
| `--dpi` | `300` | Pixels per inch (output pixel width = width × dpi) |

Default output is **1500 pixels wide** (5 in × 300 DPI).
Height is calculated to preserve the raster's aspect ratio.

---

## Rendering

`gdal_translate` is called with:
- `-of PNG` — PNG output
- `-outsize W H` — pixel dimensions calculated from width and DPI
- `-scale` — auto-stretches values to the full data range (min → 0, max → 255)
- `-ot Byte` — output as 8-bit

This gives a usable image for any raster type (DEM, mask, residual, etc.)
without requiring knowledge of the data range in advance.

---

## Example: inspect lsdtt-mountain-front outputs

After running `lsdtt-mountain-front` in `/data/tabernas/`:

```bash
python tif_to_png.py /data/tabernas/
```

Produces:

```
tabernas_synthetic_dem.png
tabernas_synthetic_front.png
tabernas_low_ground_mask.png
tabernas_fitted_surface.png
tabernas_residual.png
tabernas_mountain_mask.png
```

To convert just one file at a lower resolution:

```bash
python tif_to_png.py /data/tabernas/ --file tabernas_residual.tif --dpi 150
```
