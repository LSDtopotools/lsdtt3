# Fetch a DEM

`fetch_dem` retrieves, mosaics, projects, and returns a `Raster`; it does not
write an output file unless your script asks it to.

```python
import lsdtt3

dem = lsdtt3.fetch_dem(
    west=25.35,
    south=47.00,
    east=25.65,
    north=47.20,
    grid_spacing=30.0,
    cache_dir="dem_cache",
)

dem.write("carpathians_dem.tif")
```

Bounds are WGS84 longitude/latitude. The default `source`, `"cop30_aws"`, is
the Copernicus GLO-30 DEM on AWS, which needs no account or API key. The
returned DEM is reprojected to `target_epsg`; the default, `0`, means the UTM
zone at the centre of the box (`lsdtt3.utm_epsg_for_point(lon, lat)` gives the
code), so the pixels are square and in metres. The Copernicus DEM stores the
sea as 0 m rather than nodata: `lsdtt3.mask_to_nodata(dem, 0.0)` removes it. Use `bounds_from_vector` for a vector
study area, or `fetch_dem_for_outlets` when the task starts from outlet points
and needs complete upstream basins.

## High-resolution national and polar DEMs

`fetch_dem` reaches the same keyless national and polar sources as
`lsdtt-fetch-raster`, with the same parameter names. `lsdtt3.dem_sources()`
lists them, with each one's resolutions, surface models (DSM or DTM),
vertical datum, coverage and the attribution to cite.

```python
dem = lsdtt3.fetch_dem(
    west=-105.65, south=40.00, east=-105.55, north=40.06,
    source="us_3dep",
)
```

When `grid_spacing` is not given, a national or polar source is returned on
its native resolution (10 m for `us_3dep`), and the global sources on 30 m.
The other source options:

- `source_resolution_metres` picks one of a source's native resolutions,
  for example `2` for `arctic_arcticdem`.
- `surface_model_choice` picks `"dtm"` or `"dsm"` where a source offers
  both, for example `nl_ahn`.
- `source_collection_choice` picks one Scottish lidar phase for
  `gb_sct_lidar`. The default, `"latest"`, uses the newest phase with tiles
  in the area. Phases are never mixed.
- `max_output_pixels`, `max_download_megabytes` and `max_tiles` limit the size
  of a request, and `request_delay_milliseconds` paces the requests to tile
  services. The size limits are checked before any download;
  `max_download_megabytes` also stops a transfer that runs past it, where
  lsdtt3 does the HTTP itself.

Coverage gaps, overviews read at a coarse `grid_spacing`, and other notes
from the fetch come back as Python `UserWarning`s. Asking for an area
outside a source's coverage, or over a size limit, raises `RuntimeError`.
`nz_lidar` needs a GDAL built with LERC; it warns, and fails on a GDAL
without it. The same data are available from OpenTopography with an API
key. `fetch_dem_for_outlets` accepts only the global sources, as
`lsdtt-fetch-raster` does in outlets mode.

The [lsdtt-fetch-raster reference](../cli_reference/cli_reference_lsdtt_fetch_raster.md)
describes each source: DSM or DTM, vertical datums (ArcticDEM and REMA
heights are ellipsoidal), and the attribution each provider asks for.

See [Fetching and vector API](../python_api/fetch.md) for the return objects
and outlet-fetch options.
