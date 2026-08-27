"""Fetch a projected Copernicus DEM from a study-area vector."""

import lsdtt3


bounds = lsdtt3.bounds_from_vector("study_area.gpkg")
if not bounds.ok:
    raise RuntimeError(bounds.message)

dem = lsdtt3.fetch_dem(
    *bounds.bbox,
    grid_spacing=30.0,
    cache_dir="dem_cache",
)

dem.write("study_area_dem.tif")
print(dem.shape, dem.metadata.epsg_code)
