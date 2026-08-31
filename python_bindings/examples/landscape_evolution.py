"""Erode an uplifted random landscape for one million years."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3

# A gently south-sloping, spatially correlated random surface. The imposed
# slope provides a domain outlet; pit filling makes the initial routing surface
# drainage-ready without changing the model's subsequent in-memory workflow.
shape = (256, 256)
cell_size = 30.0
rng = np.random.default_rng(42)
noise = rng.standard_normal(shape).astype(np.float32)
for _ in range(12):
    noise = (
        noise
        + np.roll(noise, 1, axis=0)
        + np.roll(noise, -1, axis=0)
        + np.roll(noise, 1, axis=1)
        + np.roll(noise, -1, axis=1)
    ) / 5.0

rows, _ = np.indices(shape, dtype=np.float32)
elevation = 1200.0 - 1.5 * rows + 80.0 * noise / noise.std()
meta = lsdtt3.RasterMeta.from_shape(*shape, cell_size=cell_size)
initial = lsdtt3.Raster.from_numpy(
    elevation, np.ones(shape, dtype=bool), meta
)
initial = lsdtt3.fill_pits(initial)

model = lsdtt3.LandscapeModel(
    initial,
    timestep_years=1_000.0,
    erodibility=2.0e-4,
    m_exponent=0.45,
    n_exponent=1.0,
    uplift_rate=2.0e-3,
    hillslope_model="roering_implicit",
    hillslope_diffusivity=0.02,
    incision_threshold_area=1.0,
    coupling_regime="coupled",
    boundary_conditions="nnon",  # south edge is the fixed outlet
)
result = model.run(1_000, compute_erosion_rate=True)
final = result.elevation

fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
for ax, raster, title in zip(
    axes[:2], (initial, final), ("Initial random topography", "After 1 Myr")
):
    ax.imshow(raster.to_numpy(), cmap="terrain")
    ax.imshow(lsdtt3.multidirectional_hillshade(raster).to_numpy(),
              cmap="gray", alpha=0.35)
    ax.contour(model.channel_mask().to_numpy(), levels=[0.5],
               colors="dodgerblue", linewidths=0.4)
    ax.set_title(title)
    ax.set_axis_off()

erosion_rate = result.erosion_rate.to_numpy()
image = axes[2].imshow(erosion_rate, cmap="RdBu_r")
axes[2].set_title("Final signed erosion rate")
axes[2].set_axis_off()
fig.colorbar(image, ax=axes[2], label="m yr$^{-1}$")

print(f"completed {result.time_years:,.0f} years in {result.iteration} steps")
print(f"mean erosion rate: {result.mean_erosion_rate:.3e} m yr-1")
plt.show()
