"""Minimal D8 routing, channel, and HAND example on Green River."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


dem = lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif")
filled = lsdtt3.fill_pits(dem)
threshold_contributing_pixels = max(100, filled.valid_pixel_count // 250)

flow = lsdtt3.FlowInfo(filled)
channels = flow.channel_mask(threshold_contributing_pixels)
hand = flow.height_above_channel(filled, channels)


def as_nan(raster):
    return np.where(raster.mask_numpy(), raster.to_numpy(), np.nan)

fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
axes[0].imshow(as_nan(filled), cmap="terrain")
axes[0].set_title("DEM")
axes[1].imshow(np.log10(as_nan(flow.drainage_area())), cmap="viridis")
axes[1].contour(as_nan(channels), levels=[0.5], colors="white", linewidths=0.7)
axes[1].set_title("Drainage area + channels")
axes[2].imshow(as_nan(hand), cmap="magma")
axes[2].set_title("HAND")
plt.show()
