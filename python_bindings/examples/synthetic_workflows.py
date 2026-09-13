"""Create synthetic terrain, a hillslope test case, and a small OCN.

Authors
-------
B. Gailleton
"""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


meta = lsdtt3.RasterMeta.from_shape(128, 128, cell_size=30.0, epsg_code=32631)

# Reproducible correlated DEM, suitable for algorithm experiments.
noise = lsdtt3.synthetic_noise(meta, kind="diamond_square", seed=7, stddev=250.0)

# A terrain-plus-ground-truth-mask pair for mountain-front validation.
range_case = lsdtt3.synthetic_mountain_range(
    meta, geometry="arc", edge="north", front_fraction=0.35, seed=7
)

# A compact OCN experiment.  Use more iterations for quantitative branching work.
ocn = lsdtt3.synthetic_network(meta, n_iterations=5_000, relief_metres=800.0, seed=7)

fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
axes[0].imshow(noise.to_numpy(), cmap="terrain")
axes[0].set_title("Diamond-square terrain")
axes[0].set_axis_off()
axes[1].imshow(range_case["elevation"].to_numpy(), cmap="terrain")
axes[1].contour(range_case["mountain_mask"].to_numpy(), levels=[0.5], colors="white")
axes[1].set_title("Arc range + known mask")
axes[1].set_axis_off()
axes[2].imshow(np.log10(ocn.drainage_area_cells.to_numpy()), cmap="viridis")
axes[2].set_title(f"OCN area (energy {ocn.final_energy:.2g})")
axes[2].set_axis_off()
plt.show()
