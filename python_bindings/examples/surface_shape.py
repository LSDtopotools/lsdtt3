import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


def as_nan(raster):
    return np.where(raster.mask_numpy(), raster.to_numpy(), np.nan)


dem = lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif")
shape = lsdtt3.classify_surface_shape(dem, lambda_metres=90.0, kt=1e-4)

plt.imshow(as_nan(shape["shape_class"]), cmap="RdBu", vmin=-4, vmax=3)
plt.colorbar(label="surface-shape class")
plt.axis("off")
plt.show()
