import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


# Fetch a small 30 m Copernicus DEM of the Green River canyon beside
# Canyonlands National Park, Utah (projected to UTM zone 12N). The tile is
# cached in dem_cache/ so later runs do not download it again.
dem = lsdtt3.fetch_dem(
    west=-110.10, south=38.36, east=-109.96, north=38.47,
    grid_spacing=30.0, cache_dir="dem_cache",
)
front = lsdtt3.identify_mountain_front(
    dem,
    relief_scales=[[500.0, 50.0], [2000.0, 200.0]],
    return_diagnostics=True,
)

fig, axes = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)
axes[0].imshow(dem.to_numpy(), cmap="terrain")
axes[0].set_title("DEM")
axes[1].imshow(np.where(front.mountain_mask.mask_numpy(), front.mountain_mask.to_numpy(), np.nan), cmap="gray_r")
axes[1].set_title("mountain mask")
for axis in axes:
    axis.set_axis_off()
plt.show()
