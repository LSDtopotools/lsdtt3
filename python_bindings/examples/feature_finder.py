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
channels = lsdtt3.find_channel_candidates(dem, lambda_min_metres=60, lambda_max_metres=600)
valleys = lsdtt3.find_valleys(dem, smoothing_window_metres=100)

fig, axes = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)
axes[0].imshow(np.where(channels.score.mask_numpy(), channels.score.to_numpy(), np.nan), cmap="magma")
axes[0].set_title("channel candidate score")
axes[1].imshow(np.where(valleys.valley_mask.mask_numpy(), valleys.valley_mask.to_numpy(), np.nan), cmap="Blues")
axes[1].set_title("broad valley mask")
for axis in axes:
    axis.set_axis_off()
plt.show()
