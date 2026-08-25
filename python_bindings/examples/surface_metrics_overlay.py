import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


def as_nan(raster):
    return np.where(raster.mask_numpy(), raster.to_numpy(), np.nan)


dem = lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif")
threshold = np.nanpercentile(as_nan(dem), 5)
masked = lsdtt3.mask_to_nodata(dem, threshold, below=True)
scaled = lsdtt3.scale(lsdtt3.normalize_unit_interval(masked), 100.0)

shade = as_nan(lsdtt3.multidirectional_hillshade(scaled))
metrics = (
    (as_nan(lsdtt3.slope(scaled)), "Slope"),
    (as_nan(lsdtt3.plan_curvature(scaled)), "Plan curvature"),
)

fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
for axis, (metric, title) in zip(axes, metrics):
    image = axis.imshow(metric, cmap="viridis")
    axis.imshow(shade, cmap="gray", alpha=0.35)
    axis.set_title(title)
    axis.set_axis_off()
    fig.colorbar(image, ax=axis, shrink=0.8)

plt.show()
