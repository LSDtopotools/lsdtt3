# Download DEM Tiles

`download_dem_tiles.py` supports four download methods:

| Method | Flag | API key needed? |
|---|---|---|
| Microsoft Planetary Computer | *(default)* | No |
| OpenTopography REST API | `--ot-dataset` | Yes (free) |
| WCS endpoint | `--wcs-url` | Depends on service |
| Copernicus Data Space EU-DEM | `--cdse-eu-dem` | Yes for `product` downloads |

Coordinates are always given as quoted `"lat, lon"` pairs, copied directly from Google Maps.

> **Run location note:** If you are in a data directory (not `python_helpers`), call the script with its full path and set `-d .` to save outputs in the current data folder.
>
> ```bash
> python3 /mnt/c/LSDTopoTools/lsdtt3/python_helpers/download_dem_tiles.py \
>   "55.95, -3.21" "55.96, -3.20" -d .
> ```

---

## Quick credential setup with `export`

Use environment variables so you do not have to pass credentials on every call:

```bash
export OT_API_KEY="your_ot_key"
export WCS_API_KEY="your_wcs_key"
export CDSE_USERNAME="your_cdse_username"
export CDSE_PASSWORD="your_cdse_password"
```

Persist them by adding the same lines to `~/.bashrc` or `~/.zshrc`.

---

## Method 1: Microsoft Planetary Computer (default)

No account or API key required. Downloads full source tiles intersecting a point or bounding box. The default dataset is the Copernicus GLO-30 (30 m global DEM).

### Requirements

```bash
pip install requests pystac-client planetary-computer
```

### Single point — download the tile containing that point

```bash
python3 download_dem_tiles.py "55.95, -3.21"
```

### Two corners — download all tiles intersecting the bounding box

```bash
python3 download_dem_tiles.py "55.95, -3.21" "55.96, -3.20"
```

### Custom output prefix and directory

```bash
python3 download_dem_tiles.py \
  "55.95, -3.21" "55.96, -3.20" \
  --output-prefix edinburgh_dem \
  --output-dir dem_tiles
```

Output files are named `edinburgh_dem.tif` (one tile) or `edinburgh_dem_001.tif`, `edinburgh_dem_002.tif`, … (multiple tiles).

### Use the 90 m DEM instead

```bash
python3 download_dem_tiles.py "55.95, -3.21" --collection cop-dem-glo-90
```

### Key optional arguments

| Argument | Short | Default | Description |
|---|---|---|---|
| `--output-prefix` | `-o` | `my_DEM` | Prefix for output filenames |
| `--output-dir` | `-d` | `.` | Directory to save files |
| `--collection` | `-c` | `cop-dem-glo-30` | Planetary Computer collection |
| `--asset` | `-a` | `data` | Asset key within each STAC item |
| `--limit` | | `200` | Max tiles to fetch |

---

## Method 2: OpenTopography REST API

Downloads an exact clip of your bounding box (no tile seams). Supports a wider range of global and US datasets including lidar. Requires a **free** API key from [opentopography.org](https://opentopography.org).

### Requirements

```bash
pip install requests
```

### Available datasets

| Dataset | Resolution | Notes |
|---|---|---|
| `SRTMGL1` | 30 m | SRTM global |
| `SRTMGL3` | 90 m | SRTM global |
| `AW3D30` | 30 m | ALOS global |
| `COP30` | 30 m | Copernicus global |
| `COP90` | 90 m | Copernicus global |
| `SRTM15Plus` | 500 m | Ocean bathymetry included |
| `USGS30m` | 30 m | USGS (USA only) |
| `USGS10m` | 10 m | USGS (USA only) |
| `USGS1m` | 1 m | USGS (USA only, academic account required) |

### Providing your API key

The script looks for the key in this order — first found wins:

**1. Environment variable (recommended)**

```bash
export OT_API_KEY="your_key_here"
```

Set this in your shell profile (`~/.bashrc` or `~/.zshrc`) so you never have to think about it again.

**2. Key file**

```bash
python3 download_dem_tiles.py ... --ot-api-key-file ~/ot_key.txt
```

The file should contain just the key on the first line. Protect it and keep it out of any repository:

```bash
chmod 600 ~/ot_key.txt
echo "ot_key.txt" >> ~/.gitignore
```

**3. Interactive prompt**

If neither of the above is set the script will prompt you:

```
OpenTopography API key:
```

The key is not echoed to the terminal.

### Example calls

```bash
# Copernicus 30 m, key from environment variable
export OT_API_KEY="your_key_here"
python3 download_dem_tiles.py \
  "37.036, -2.408" "37.063, -2.374" \
  --ot-dataset COP30 \
  --output-prefix tabernas

# SRTM 30 m, key from file
python3 download_dem_tiles.py \
  "55.95, -3.21" "55.96, -3.20" \
  --ot-dataset SRTMGL1 \
  --ot-api-key-file ~/ot_key.txt \
  --output-prefix edinburgh_dem \
  --output-dir dem_tiles

# USGS 10 m (USA only), interactive key prompt
python3 download_dem_tiles.py \
  "37.20, -112.05" "37.25, -111.98" \
  --ot-dataset USGS10m \
  --output-prefix bryce_canyon
```

Output filename includes the dataset name, e.g. `tabernas_COP30.tif`.

---

## Method 3: WCS endpoint

Downloads from any OGC Web Coverage Service. Useful for national mapping agency datasets (e.g. Spanish IGN, Copernicus EU DEM via EEA). The bounding box is passed directly in the requested CRS — since lsdtt3 handles reprojection, the default CRS is `EPSG:4326` (WGS84 lat/lon).

The downloader now does a `GetCapabilities` preflight before requesting coverage, so endpoint/version issues fail early with a clear error message.

### Requirements

```bash
pip install requests owslib
```

### API keys for WCS

Many public WCS services require no key. If one is needed, the same priority order applies:

1. Environment variable: `export WCS_API_KEY="your_key_here"`
2. Key file: `--api-key-file /path/to/keyfile.txt`
3. Interactive prompt

The key is appended to the URL as a query parameter. The parameter name defaults to `key`; override with `--api-key-param`.

### Example calls

```bash
# Spanish IGN 5 m DEM (no key required), default EPSG:4326
python3 download_dem_tiles.py \
  "37.036, -2.408" "37.063, -2.374" \
  --wcs-url "https://servicios.idee.es/wcs-inspire/mdt" \
  --wcs-coverage "Elevacion25830_5" \
  --wcs-res 0.00005 \
  --output-prefix tabernas_ign

# NRW DGM WCS (no key required, NRW extent only)
python3 download_dem_tiles.py \
  "51.45, 7.00" "51.55, 7.20" \
  --wcs-url "https://www.wcs.nrw.de/geobasis/wcs_nw_dgm" \
  --wcs-coverage "nw_dgm" \
  --wcs-version 1.0.0 \
  --output-prefix nrw_dgm

# Service requiring an API key passed as ?token=...
python3 download_dem_tiles.py \
  "55.95, -3.21" "55.96, -3.20" \
  --wcs-url "https://example.service/wcs" \
  --wcs-coverage "MyDEM" \
  --api-key-file ~/wcs_key.txt \
  --api-key-param token \
  --output-prefix my_dem
```

### Key WCS optional arguments

| Argument | Default | Description |
|---|---|---|
| `--wcs-url` | | WCS endpoint URL (triggers WCS mode) |
| `--wcs-coverage` | | Coverage ID (required) |
| `--wcs-version` | `1.0.0` | WCS protocol version |
| `--wcs-crs` | `EPSG:4326` | CRS for the request |
| `--wcs-res` | | Cell size in CRS units |
| `--wcs-format` | `image/tiff` | Requested output format |
| `--api-key-file` | | Path to WCS key file |
| `--api-key-param` | `key` | URL parameter name for the key |

### Endpoint notes

- `https://www.wcs.nrw.de/geobasis/wcs_nw_dgm` responds correctly for WCS `1.0.0`, but coverage extent is Nordrhein-Westfalen only (does not include Stuttgart/Baden-Wuerttemberg).
- The previously used EEA EU-DEM WCS endpoint has recently returned HTTP 400 for `GetCapabilities`; use CDSE mode (`--cdse-eu-dem`) for EU-DEM when that WCS endpoint is unavailable.

---

## Method 4: Copernicus Data Space EU-DEM (STAC)

Downloads DEM tiles from the Copernicus Data Space Ecosystem STAC catalogue.

- Default CDSE collection: `cop-dem-eea-10-laea-tif` (EU-DEM 10 m)
- Default asset: `product` (zip download via CDSE OData endpoint)
- You can request other CDSE DEM collections with `--cdse-collection`.

### Requirements

```bash
pip install requests pystac-client
```

### Credentials for CDSE mode

For `--cdse-asset product`, authentication is required. Credential resolution order:

1. `--cdse-username` and `--cdse-password`
2. `CDSE_USERNAME` and `CDSE_PASSWORD` environment variables
3. Interactive prompt

Recommended:

```bash
export CDSE_USERNAME="your_cdse_username"
export CDSE_PASSWORD="your_cdse_password"
```

### Example calls

```bash
# EU-DEM 10 m by bbox (default CDSE collection)
python3 download_dem_tiles.py \
  "57.0, -4.5" "58.0, -3.0" \
  --cdse-eu-dem \
  --output-prefix feshie_eudem

# EU-DEM 10 m point query, explicit credentials
python3 download_dem_tiles.py \
  "55.95, -3.21" \
  --cdse-eu-dem \
  --cdse-username "your_user" \
  --cdse-password "your_password" \
  --output-prefix edinburgh_eudem

# Global Copernicus DEM 30 m via CDSE STAC
python3 download_dem_tiles.py \
  "55.95, -3.21" "55.96, -3.20" \
  --cdse-eu-dem \
  --cdse-collection cop-dem-glo-30-dged-cog \
  --output-prefix copdem30_cdse
```

### Key CDSE optional arguments

| Argument | Default | Description |
|---|---|---|
| `--cdse-eu-dem` | off | Enable CDSE STAC mode |
| `--cdse-collection` | `cop-dem-eea-10-laea-tif` | CDSE STAC collection |
| `--cdse-asset` | `product` | Asset key to download (`product` or `data`) |
| `--cdse-username` | | CDSE username (or `CDSE_USERNAME`) |
| `--cdse-password` | | CDSE password (or `CDSE_PASSWORD`) |

### Recommended presets

Use these copy/paste presets for common DEM choices:

```bash
# Preset A: EU-DEM 10 m (Europe)
# Collection: cop-dem-eea-10-laea-tif
python3 download_dem_tiles.py \
  "57.0, -4.5" "58.0, -3.0" \
  --cdse-eu-dem \
  --cdse-collection cop-dem-eea-10-laea-tif \
  --cdse-asset product \
  --output-prefix eudem10

# Preset B: Copernicus DEM Global 30 m
# Collection: cop-dem-glo-30-dged-cog
python3 download_dem_tiles.py \
  "57.0, -4.5" "58.0, -3.0" \
  --cdse-eu-dem \
  --cdse-collection cop-dem-glo-30-dged-cog \
  --cdse-asset product \
  --output-prefix copdem30

# Preset C: Copernicus DEM Global 90 m
# Collection: cop-dem-glo-90-dged-cog
python3 download_dem_tiles.py \
  "57.0, -4.5" "58.0, -3.0" \
  --cdse-eu-dem \
  --cdse-collection cop-dem-glo-90-dged-cog \
  --cdse-asset product \
  --output-prefix copdem90
```

---

## Running from a data folder

The script can be called from any directory. Use `--output-dir` (or `-d`) to control where files are saved, regardless of where you are when you run the command.

```bash
# cd into your data folder first, then save there with -d .
cd /my/data/folder
python3 /mnt/c/LSDTopoTools/lsdtt3/python_helpers/download_dem_tiles.py \
  "55.95, -3.21" -d .

# Or stay anywhere and point -d at the target folder
python3 /mnt/c/LSDTopoTools/lsdtt3/python_helpers/download_dem_tiles.py \
  "55.95, -3.21" -d /my/data/folder
```

### Avoiding the full path

**Option 1 — shell alias (recommended)**

Add to your `~/.bashrc` (or `~/.zshrc`):

```bash
alias download_dem_tiles='python3 /mnt/c/LSDTopoTools/lsdtt3/python_helpers/download_dem_tiles.py'
```

Reload with `source ~/.bashrc`, then call it directly from any folder:

```bash
download_dem_tiles "55.95, -3.21" -d .
```

**Option 2 — add the directory to PATH**

Make the script executable once, then add its directory to PATH:

```bash
chmod +x /mnt/c/LSDTopoTools/lsdtt3/python_helpers/download_dem_tiles.py
```

Add to `~/.bashrc`:

```bash
export PATH="$PATH:/mnt/c/LSDTopoTools/lsdtt3/python_helpers"
```

Then call it without `python3`:

```bash
download_dem_tiles.py "55.95, -3.21" -d .
```

**Option 3 — symlink into `~/.local/bin`**

A symlink keeps your PATH clean while pointing back to the live script file:

```bash
mkdir -p ~/.local/bin
ln -s /mnt/c/LSDTopoTools/lsdtt3/python_helpers/download_dem_tiles.py \
      ~/.local/bin/download_dem_tiles
chmod +x /mnt/c/LSDTopoTools/lsdtt3/python_helpers/download_dem_tiles.py
```

**Note (WSL2):** if the script was edited on Windows it may have `\r\n` line endings, causing the error `/usr/bin/env: 'python3\r': No such file or directory`. Fix once with:

```bash
sed -i 's/\r//' /mnt/c/LSDTopoTools/lsdtt3/python_helpers/download_dem_tiles.py
```

`~/.local/bin` is usually already on PATH in modern Linux/WSL2. Then:

```bash
download_dem_tiles "55.95, -3.21" -d .
```

---

## Notes

- Coordinates must be quoted so the shell treats each `lat, lon` pair as one argument.
- Planetary Computer downloads **full source tiles** and does not clip or merge them.
- OpenTopography and WCS both return an **exact clip** of the bounding box.
- lsdtt3 handles reprojection, so downloading in `EPSG:4326` is fine for all methods.
- In CDSE mode, multiple versions of the same DEM grid tile can exist; the helper keeps the latest item per `grid:code`.
