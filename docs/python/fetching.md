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

Bounds are WGS84 longitude/latitude. The returned DEM is projected to a local
working CRS suitable for its extent. Use `bounds_from_vector` for a vector
study area, or `fetch_dem_for_outlets` when the task starts from outlet points
and needs complete upstream basins.

See [Fetching and vector API](../python_api/fetch.md) for the return objects
and outlet-fetch options.
