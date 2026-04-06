# Download DEM Tiles from Microsoft Planetary Computer


This script downloads DEM tiles from Microsoft Planetary Computer using coordinates copied from Google Maps. We use Planetary Computer because it does not require an API token. 

Its functionality is somewhat limitied, you can only get the Copernicus DEM and you get the data in tiles, so you will need to do postprocesing for clipping. If you want greater flexibility you can try the OpenTopography API for different formats, exact clips, and lidar data but this requires an API key so we keep this lighter weight version if you don't have that or don't want to bother with an API key. 

---

Google Maps gives coordinates in this format:

```text
55.95244789111162, -3.2193971081922665
```

That is:

```text
latitude, longitude
```

This script accepts either:

- **one coordinate string**: download the full DEM tile containing that point
- **two coordinate strings**: download all full DEM tiles intersecting the bounding box defined by those two corners

---

## Requirements

Install the required Python packages:

```bash
pip install requests pystac-client planetary-computer
```

---

## Usage

### 1. Download the full tile containing a single point

```bash
python3 download_dem_tiles.py "55.95244789111162, -3.2193971081922665"
```

This will usually save the output as:

```text
my_DEM.tif
```

---

### 2. Download all tiles intersecting a box

Pass two opposite corners:

```bash
python3 download_dem_tiles.py \
  "55.95244789111162, -3.2193971081922665" \
  "55.96244789111162, -3.2093971081922665"
```

If multiple tiles are found, they will be saved as:

```text
my_DEM_001.tif
my_DEM_002.tif
my_DEM_003.tif
```

You can paste the two corners in any order.

---

## Choose the output filename prefix

Use `--output-prefix`:

```bash
python3 download_dem_tiles.py \
  "55.95244789111162, -3.2193971081922665" \
  --output-prefix edinburgh_dem
```

Single-tile output:

```text
edinburgh_dem.tif
```

Multiple-tile output:

```text
edinburgh_dem_001.tif
edinburgh_dem_002.tif
```

---

## Save files to a specific folder

Use `--output-dir`:

```bash
python3 download_dem_tiles.py \
  "55.95244789111162, -3.2193971081922665" \
  --output-prefix my_DEM \
  --output-dir dem_tiles
```

---

## Choose DEM collection

By default the script uses:

```text
cop-dem-glo-30
```

To use the 90 m DEM instead:

```bash
python3 download_dem_tiles.py \
  "55.95244789111162, -3.2193971081922665" \
  --collection cop-dem-glo-90
```

---

## Optional arguments

### Output prefix

```bash
-o my_DEM
--output-prefix my_DEM
```

### Output directory

```bash
-d dem_tiles
--output-dir dem_tiles
```

### Collection

```bash
-c cop-dem-glo-30
--collection cop-dem-glo-30
```

### Asset key

Default:

```text
data
```

You can override it if needed:

```bash
--asset data
```

### Search limit

Default:

```text
200
```

Override with:

```bash
--limit 200
```

---

## Notes

- Coordinates must be quoted so the shell treats each `lat, lon` pair as one argument.
- The script downloads **full source tiles**.
- It does **not** clip, crop, or merge rasters.
- With one coordinate, the script searches for tile(s) containing that point.
- With two coordinates, the script searches for tile(s) intersecting the bounding box formed by those corners.

---

## Make the script executable

```bash
chmod +x download_dem_tiles.py
```

Then run it as:

```bash
./download_dem_tiles.py "55.95244789111162, -3.2193971081922665"
```

or:

```bash
./download_dem_tiles.py \
  "55.95244789111162, -3.2193971081922665" \
  "55.96244789111162, -3.2093971081922665"
```

---

## Examples

### Single point, default output name

```bash
python3 download_dem_tiles.py "55.95244789111162, -3.2193971081922665"
```

### Single point, custom output name

```bash
python3 download_dem_tiles.py \
  "55.95244789111162, -3.2193971081922665" \
  --output-prefix edinburgh_dem
```

### Two corners, custom output folder

```bash
python3 download_dem_tiles.py \
  "55.95244789111162, -3.2193971081922665" \
  "55.96244789111162, -3.2093971081922665" \
  --output-prefix edinburgh_dem \
  --output-dir dem_tiles
```

