# A note on datums and PROJ grids

**Short version.** Changing from one projection to another on the same datum
(for example WGS 84 longitude and latitude to WGS 84 / UTM) is exact. Changing
**datum** is not: the most accurate methods use a correction grid. lsdtt3
carries a small set of grids. If a conversion needs one it does not carry,
the tool prints a warning naming the grid and still runs, using a less
accurate method that is usually off by about a metre, sometimes several, and
by tens of metres for some old datums. If you need metre-scale accuracy, fetch
the grid it names (see [Getting a missing grid](#getting-a-missing-grid)).

## Why datums matter

A coordinate system has two parts: a **datum** (the shape and position of the
Earth model the coordinates refer to) and a **projection** (how that model is
flattened onto a map). Two datums can differ by tens to hundreds of metres, and
that difference is not a simple shift: it bends across a country. National
mapping agencies publish the correction as a grid, and PROJ (the library GDAL,
QGIS and lsdtt3 use for coordinates) reads those grids.

Without the grid, PROJ falls back to a formula fitted over a whole country
(a Helmert transformation) or, for some datums, to a rough "ballpark"
conversion. It does not stop: the coordinates look right, in the right place,
and are wrong by a metre or more. At 1 m lidar resolution that is a shift of
one to several cells.

Vertical datums are the same story. Heights "above sea level" are measured
above a geoid, and global DEMs use different ones; converting between them
needs a geoid grid.

## What lsdtt3 carries

Every release package and every pip wheel includes these grids (about 6 MB):

| Grid | What it is for |
|---|---|
| `uk_os_OSTN15_NTv2_OSGBtoETRS.tif` | British National Grid (OSGB36, EPSG:27700) to ETRS89 / WGS 84, from the Ordnance Survey |
| `us_nga_egm96_15.tif` | the EGM96 geoid (heights of SRTM, NASADEM, AW3D30) |
| `us_noaa_conus.tif` | NAD27 to NAD83 in the contiguous US (NADCON) |
| `us_noaa_nadcon5_nad27_nad83_1986_conus.tif` | NAD27 to NAD83(1986) in the contiguous US (NADCON5) |

The full PROJ grid collection is about 760 MB, too large to bundle. In
particular **the EGM2008 geoid grid (`us_nga_egm08_25.tif`, 77 MB) is not
included.**

To see what your installation has, run any tool with `-v`:

```
  proj grids: 4 of 4 default; grid check on
```

## Where this matters

**No grid needed (exact):**

- WGS 84 longitude/latitude to UTM, and between UTM zones on WGS 84. This is
  what `lsdtt-fetch-raster` does with Copernicus, SRTM, NASADEM and the other
  global DEMs.
- Polar stereographic and Lambert azimuthal equal-area projections on WGS 84
  (EPSG:3413, 3031, 3995, 3976, 6931 and similar): ArcticDEM and REMA.
- ETRS89-based European grids (EPSG:258xx) to WGS 84: PROJ treats the two as the
  same, which is good to about half a metre to a metre (the plates have moved
  since 1989). A grid would not change that.

**Covered by the included grids:**

- **British National Grid (EPSG:27700) to UTM or WGS 84**: Scottish and
  English lidar, Ordnance Survey data, and inventories drawn in BNG. Without
  OSTN15 the result is off by 0.5 to 4 m (measured: 2 m at Glen Coe and Arran,
  3.7 m in Cornwall).
- **NAD27 to NAD83 or WGS 84 in the contiguous US**: old USGS maps and
  datasets. Without the grids, errors are tens of metres.

**Not covered: you will see a warning, or for geoids an error:**

- **Converting heights to or from EGM2008** with `lsdtt-raster-calculator`
  (`first_raster_vertical_datum` / `second_raster_vertical_datum` set to
  `egm2008`). This is how you remove the geoid difference between Copernicus
  GLO-30 (EGM2008) and SRTM, NASADEM or AW3D30 (EGM96) before differencing
  them. The difference is typically decimetres and reaches metres in high
  mountains. The run stops before reading the rasters, naming
  `us_nga_egm08_25.tif` and the command below. If the datums are not declared,
  nothing is converted, and the calculator prints a note that the difference
  contains that offset.
- **US data on NAD83 (most 3DEP lidar, state plane coordinates) to WGS 84 or
  UTM.** PROJ's most accurate method uses the state's NOAA HARN grid
  (`us_noaa_<state>hpgn.tif`, for example `us_noaa_cohpgn.tif` for Colorado).
  Without it the fallback is good to about 4 m instead of 2 m. Data that stays
  in NAD83 / UTM throughout needs no conversion at all.
- **Other national grids**: New Zealand (NZGD49), Australia (AGD66/84), older
  European national datums, and others. Any conversion whose best method uses a
  grid lsdtt3 does not carry produces the warning.

Data that stays in one coordinate system from start to finish never needs a
grid. Two DEMs converted the same way share the same error, so their
difference is unaffected; the error matters when you combine lsdtt3 output with
data converted elsewhere (for example in QGIS, which usually has the full
grid set), or when absolute position matters (field points, GNSS, inventories).

## The warning

When a tool reprojects a raster or vector layer, it checks the conversion over
the data's area. If the most accurate method needs a grid that is not
installed, it prints one warning per pair of coordinate systems, for example:

```
WARNING: reprojecting the raster: the most accurate transformation (Inverse of
British National Grid + OSGB36 to WGS 84 (9) + UTM zone 30N, accuracy 1 m) needs
the PROJ grid uk_os_OSTN15_NTv2_OSGBtoETRS.tif, which is not installed. Using
Inverse of British National Grid + OSGB36 to WGS 84 (6) + UTM zone 30N
(accuracy 2 m) instead. To install it: python3 fetch_proj_grids.py fetch
uk_os_OSTN15_NTv2_OSGBtoETRS.tif
```

The run continues. The accuracies are those published for each method. A
source build without the PROJ development headers has no check, and `-v` says
`grid check off`.

## Getting a missing grid

Grids are installed into lsdtt3's own `share/proj` folder, which the tools
always use (grids installed elsewhere, for example by `projsync`, conda or
Homebrew, are not seen). The script `tools/fetch_proj_grids.py` does this. It
needs only Python 3; if you do not have the lsdtt3 repository, download it on
its own:

```bash
curl -O https://raw.githubusercontent.com/LSDtopotools/lsdtt3/main/tools/fetch_proj_grids.py
```

Run it in the shell or environment where you use lsdtt3 (after
`source setenv.sh` for a release package, or in the Python environment where
you ran `pip install lsdtt3`), so that it finds the right folder:

```bash
python3 fetch_proj_grids.py where                       # which share/proj it will use
python3 fetch_proj_grids.py list                        # what is installed there
python3 fetch_proj_grids.py fetch us_nga_egm08_25.tif   # add a grid (the name from the warning)
```

If it cannot find the folder, give it: `--dir <lsdtt3 folder>/share/proj`.
Grids come from the PROJ content delivery network, <https://cdn.proj.org>,
where every available grid is listed. The default grids are checked against
known checksums.

A new release package is a new folder with the default set, so fetch your
extra grids again after upgrading it. A pip upgrade (`pip install -U lsdtt3`)
replaces only the files the wheel installed, so grids you added stay; check
with `list`.

**Without installing anything**, you can also let PROJ read grids from the CDN
as it needs them, by setting the environment variable `PROJ_NETWORK=ON` before
running a tool. That needs a network connection on every run.

## See also

- [Raster data](raster_data.md#coordinate-reference-systems) and
  [Vector data](vector_data.md#coordinate-systems)
- [If you already have GDAL or PROJ installed](install_lsdtt3.md#if-you-already-have-gdal-or-proj-installed)
- PROJ's own documentation of grids: <https://proj.org/en/latest/resource_files.html>
