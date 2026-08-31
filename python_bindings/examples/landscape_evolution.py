"""Coupled fluvial and hillslope landscape evolution on a small synthetic DEM."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3

meta = lsdtt3.RasterMeta.from_shape(80, 100, cell_size=20.0)
rows, cols = np.indices((80, 100), dtype=np.float32)
z = 800.0 - 2.0 * rows - 0.012 * (cols - 50.0) ** 2
dem = lsdtt3.Raster.from_numpy(z, np.ones(z.shape, dtype=bool), meta)

model = lsdtt3.LandscapeModel(
    dem,
    timestep_years=100.0,
    erodibility=2.0e-5,
    uplift_rate=2.0e-4,
    hillslope_model="linear",
    hillslope_diffusivity=0.02,
    incision_threshold_area=5_000.0,
    coupling_regime="channels_as_boundaries",
    boundary_conditions="nnon",  # fixed south edge outlet
)
result = model.run(25, compute_erosion_rate=True)

fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
for ax, raster, title in zip(
    axes[:2], (dem, result.elevation), ("Initial elevation", "After 2,500 years")
):
    image = ax.imshow(raster.to_numpy(), cmap="terrain")
    ax.contour(model.channel_mask().to_numpy(), levels=[0.5], colors="k", linewidths=0.5)
    ax.set_title(title)
    fig.colorbar(image, ax=ax, label="elevation (m)")

rate = result.erosion_rate.to_numpy()
image = axes[2].imshow(rate, cmap="RdBu_r")
axes[2].set_title("Final signed erosion rate")
fig.colorbar(image, ax=axes[2], label="m yr$^{-1}$")
plt.show()
