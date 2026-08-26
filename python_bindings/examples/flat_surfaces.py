import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


dem = lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif")
surfaces = lsdtt3.find_flat_surfaces(
    dem,
    madk2_window_radius_metres=100.0,
    use_residual_iqr=False,
    min_component_pixels=100,
)

fig, axes = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)
axes[0].imshow(np.where(dem.mask_numpy(), dem.to_numpy(), np.nan), cmap="terrain")
axes[0].set_title("DEM")
axes[1].imshow(
    np.where(surfaces.labelled_flat_mask.mask_numpy(),
             surfaces.labelled_flat_mask.to_numpy(), np.nan),
    cmap="tab20",
)
axes[1].set_title("labelled flat surfaces")
for axis in axes:
    axis.set_axis_off()
plt.show()
