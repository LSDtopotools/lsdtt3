import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon

import lsdtt3


def as_nan(raster):
    return np.where(raster.mask_numpy(), raster.to_numpy(), np.nan)


# Fetch a small 30 m Copernicus DEM of the Green River canyon beside
# Canyonlands National Park, Utah (projected to UTM zone 12N). The tile is
# cached in dem_cache/ so later runs do not download it again.
dem = lsdtt3.fetch_dem(
    west=-110.10, south=38.36, east=-109.96, north=38.47,
    grid_spacing=30.0, cache_dir="dem_cache",
)
threshold = np.nanpercentile(as_nan(dem), 5)
masked = lsdtt3.mask_to_nodata(dem, threshold, below=True)
scaled = lsdtt3.scale(lsdtt3.normalize_unit_interval(masked), 100.0)

rows, cols = scaled.shape
half_size = max(1, min(rows, cols) // 20)
with_hole = lsdtt3.set_nodata_rectangle(
    scaled,
    rows // 2 - half_size,
    rows // 2 + half_size,
    cols // 2 - half_size,
    cols // 2 + half_size,
)
repaired = lsdtt3.fill_interior_nodata(with_hole)
filled = lsdtt3.fill_pits_boundary_aware(
    repaired, boundary_conditions="nnon", min_slope=0.0001
)

meta = filled.metadata
x0 = meta.x_min + 0.2 * (meta.x_max - meta.x_min)
x1 = meta.x_min + 0.8 * (meta.x_max - meta.x_min)
y0 = meta.y_min + 0.2 * (meta.y_max - meta.y_min)
y1 = meta.y_min + 0.8 * (meta.y_max - meta.y_min)
polygon = gpd.GeoDataFrame(
    geometry=[Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])],
    crs=meta.crs_wkt,
)
clipped = lsdtt3.clip_to_polygons(repaired, polygon)
slope = lsdtt3.slope(clipped)
shade = lsdtt3.multidirectional_hillshade(clipped)

fig, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)

axes[0, 0].imshow(as_nan(dem), cmap="terrain")
axes[0, 0].set_title("Original DEM")

axes[0, 1].imshow(as_nan(repaired), cmap="terrain")
axes[0, 1].set_title("Thresholded, scaled, and repaired")

fill_difference = as_nan(filled) - as_nan(repaired)
axes[1, 0].imshow(fill_difference, cmap="magma")
axes[1, 0].set_title("Boundary-aware pit-fill difference")

slope_values = as_nan(slope)
vmin, vmax = np.nanpercentile(slope_values, (2, 98))
image = axes[1, 1].imshow(slope_values, cmap="magma", vmin=vmin, vmax=vmax)
axes[1, 1].imshow(as_nan(shade), cmap="gray", alpha=0.35)
axes[1, 1].set_title("Polygon-clipped slope with hillshade")
fig.colorbar(image, ax=axes[1, 1], shrink=0.8)

for axis in axes.flat:
    axis.set_axis_off()

plt.show()
