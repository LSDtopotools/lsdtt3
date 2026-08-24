import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


def display_array(raster):
    return np.ma.masked_where(~raster.mask_numpy(), raster.to_numpy())


dem = lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif")
rows, cols = dem.shape
with_nodata = lsdtt3.set_nodata_rectangle(
    dem, rows // 2 - 10, rows // 2 + 10, cols // 2 - 10, cols // 2 + 10
)
filled = lsdtt3.fill_interior_nodata(with_nodata)

fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
axes[0].imshow(display_array(with_nodata), cmap="terrain")
axes[0].set_title("Interior nodata")
axes[1].imshow(filled.to_numpy() - dem.to_numpy(), cmap="coolwarm")
axes[1].set_title("Filled − original")
for axis in axes:
    axis.set_axis_off()

plt.show()
