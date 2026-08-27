"""Fetch a small projected Copernicus DEM in the Romanian Carpathians."""

import lsdtt3


dem = lsdtt3.fetch_dem(
    west=25.35,
    south=47.001,
    east=25.65,
    north=47.20,
    grid_spacing=30.0,
    cache_dir="dem_cache",
)

dem.write("carpathians_dem.tif")
print(dem.shape)
