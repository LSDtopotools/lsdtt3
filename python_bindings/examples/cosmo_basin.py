"""Invert basin-average 10Be erosion from point samples with native CAIRN."""

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point

import lsdtt3


# Fetch a small 30 m Copernicus DEM of the Green River canyon beside
# Canyonlands National Park, Utah (projected to UTM zone 12N). The tile is
# cached in dem_cache/ so later runs do not download it again.
dem = lsdtt3.fetch_dem(
    west=-110.10, south=38.36, east=-109.96, north=38.47,
    grid_spacing=30.0, cache_dir="dem_cache",
)
dem = lsdtt3.fill_pits(dem)
flow = lsdtt3.FlowInfo(dem)
production = lsdtt3.crn_production_scaling(dem, atmosphere="standard")
meta = dem.metadata

# In real work these are your measured outlet samples. They may be in another
# CRS: the binding reprojects them to the flow-grid CRS before snapping.
x = meta.x_min + (meta.n_cols // 2 + 0.5) * meta.cell_size
y = meta.y_max - (meta.n_rows // 2 + 0.5) * meta.cell_size
samples = gpd.GeoDataFrame(
    {"sample": ["example"], "concentration": [200_000.0]},
    geometry=[Point(x, y)], crs=meta.crs_wkt,
)

result = lsdtt3.crn_invert_basin_erosion(
    flow, samples, production,
    channel_threshold_pixels=1,  # use a science-appropriate threshold for real samples
    snap_radius_pixels=10,
)
print(result[["sample", "crn_status", "erosion_rate_g_cm2_yr", "erosion_rate_mm_kyr"]])

terrain = np.where(dem.mask_numpy(), dem.to_numpy(), np.nan)
fig, ax = plt.subplots(figsize=(8, 6))
ax.imshow(terrain, cmap="terrain", extent=(meta.x_min, meta.x_max, meta.y_min, meta.y_max))
ax.scatter(result.input_x, result.input_y, marker="x", color="white", s=70, label="input sample")
ax.scatter(result.snapped_x, result.snapped_y, marker="o", color="crimson", s=35, label="snapped D8 outlet")
ax.legend()
ax.set_title(r"CAIRN basin inversion: $^{10}$Be outlet sample")
plt.show()
