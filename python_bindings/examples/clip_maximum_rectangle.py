import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


def display_array(raster):
    return np.ma.masked_where(~raster.mask_numpy(), raster.to_numpy())


# Fetch a small 30 m Copernicus DEM of the Green River canyon beside
# Canyonlands National Park, Utah (projected to UTM zone 12N). The tile is
# cached in dem_cache/ so later runs do not download it again.
dem = lsdtt3.fetch_dem(
    west=-110.10, south=38.36, east=-109.96, north=38.47,
    grid_spacing=30.0, cache_dir="dem_cache",
)
rows, cols = dem.shape
with_edge_nodata = lsdtt3.set_nodata_rectangle(dem, 0, rows, 0, cols // 10)
clipped = lsdtt3.clip_to_maximum_rectangle(with_edge_nodata)

fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
axes[0].imshow(display_array(with_edge_nodata), cmap="terrain")
axes[0].set_title("Edge nodata")
axes[1].imshow(display_array(clipped), cmap="terrain")
axes[1].set_title("Maximum valid rectangle")
for axis in axes:
    axis.set_axis_off()

plt.show()
