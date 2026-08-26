import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


dem = lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif")
channels = lsdtt3.find_channel_candidates(dem, lambda_min_metres=10, lambda_max_metres=100)
valleys = lsdtt3.find_valleys(dem, smoothing_window_metres=100)

fig, axes = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)
axes[0].imshow(np.where(channels.score.mask_numpy(), channels.score.to_numpy(), np.nan), cmap="magma")
axes[0].set_title("channel candidate score")
axes[1].imshow(np.where(valleys.valley_mask.mask_numpy(), valleys.valley_mask.to_numpy(), np.nan), cmap="Blues")
axes[1].set_title("broad valley mask")
for axis in axes:
    axis.set_axis_off()
plt.show()
