import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


def display_array(raster):
    return np.ma.masked_where(~raster.mask_numpy(), raster.to_numpy())


dem = lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif")
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
