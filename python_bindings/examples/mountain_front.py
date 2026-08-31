import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


dem = lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif")
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
