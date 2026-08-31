"""Solve a channel-bounded hillslope surface and inspect one transport step."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


meta = lsdtt3.RasterMeta.from_shape(31, 31, cell_size=5.0)
rows, cols = np.indices((31, 31), dtype=np.float32)
initial = 120.0 + 0.03 * (cols - 15.0) ** 2 + 0.02 * rows
dem = lsdtt3.Raster.from_numpy(initial, np.ones_like(initial, dtype=bool), meta)

# The two edge rows are fixed channel/base-level elevations.
channel_data = np.zeros_like(initial)
channel_data[[0, -1], :] = 1.0
channels = lsdtt3.Raster.from_numpy(
    channel_data, np.ones_like(channel_data, dtype=bool), meta
)

steady = lsdtt3.hillslope_steady_state(
    dem, channels, diffusivity=0.05, uplift_rate=4e-4,
    model="roering", nonlinear=False, boundary_conditions="nnnn", tolerance=1e-4,
    sor_omega=1.7,
)
evolved = lsdtt3.hillslope_transport_implicit(
    steady.elevation, channels, duration_years=100.0,
    diffusivity=0.05, uplift_rate=4e-4, boundary_conditions="nnnn",
)

fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
extent = [meta.x_min, meta.x_max, meta.y_min, meta.y_max]
for axis, values, title in zip(
    axes,
    [dem.to_numpy(), steady.elevation.to_numpy(), evolved.elevation.to_numpy()],
    ["Initial surface", "Linear steady state", "After 100 years"],
):
    image = axis.imshow(values, cmap="terrain", extent=extent, origin="upper")
    axis.contour(channel_data, levels=[0.5], colors="deepskyblue", linewidths=1.0,
                 extent=extent, origin="upper")
    axis.set(title=title, xlabel="x (m)", ylabel="y (m)")
    fig.colorbar(image, ax=axis, label="elevation (m)")

plt.show()
