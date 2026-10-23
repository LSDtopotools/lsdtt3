"""Compare MFD and D-infinity drainage area on Green River."""

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
dem = lsdtt3.fill_pits(dem)
mfd = lsdtt3.FlowGraph(dem, method="mfd")
dinf = lsdtt3.FlowGraph(dem, method="dinf")

valid = dem.mask_numpy()
mfd_area = np.where(valid, np.log10(mfd.drainage_area().to_numpy()), np.nan)
dinf_area = np.where(valid, np.log10(dinf.drainage_area().to_numpy()), np.nan)

fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
axes[0].imshow(mfd_area, cmap="viridis")
axes[0].set_title("MFD log drainage area")
axes[1].imshow(dinf_area, cmap="viridis")
axes[1].set_title("D∞ log drainage area")
axes[2].imshow(dinf.direction().to_numpy(), cmap="twilight", vmin=0, vmax=360)
axes[2].set_title("D∞ direction (degrees)")
for axis in axes:
    axis.set_axis_off()
plt.show()
