import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


def as_nan(raster):
    return np.where(raster.mask_numpy(), raster.to_numpy(), np.nan)


dem = lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif")
metrics = lsdtt3.polyfit_metrics(
    dem,
    ["slope_degrees", "mean_curvature", "plan_curvature"],
    lambda_metres=90.0,
)

fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
for axis, (name, raster) in zip(axes, metrics.items()):
    image = axis.imshow(as_nan(raster), cmap="terrain")
    axis.set_title(name.replace("_", " "))
    axis.set_axis_off()
    fig.colorbar(image, ax=axis, shrink=0.8)

plt.show()
