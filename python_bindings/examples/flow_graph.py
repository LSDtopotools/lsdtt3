"""Compare MFD and D-infinity drainage area on Green River."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


dem = lsdtt3.fill_pits(lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif"))
mfd = lsdtt3.FlowGraph(dem, method="mfd")
dinf = lsdtt3.FlowGraph(dem, method="dinf")

valid = dem.mask_numpy()
mfd_area = np.where(valid, np.log10(mfd.drainage_area().to_numpy()), np.nan)
dinf_area = np.where(valid, np.log10(dinf.drainage_area().to_numpy()), np.nan)

fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
axes[0].imshow(mfd_area, cmap="viridis")
axes[0].set_title("MFD log drainage area")
axes[1].imshow(dinf_area, cmap="viridis")
axes[1].set_title("D∞ log drainage area")
axes[2].imshow(dinf.direction().to_numpy(), cmap="twilight", vmin=0, vmax=360)
axes[2].set_title("D∞ direction (degrees)")
for axis in axes:
    axis.set_axis_off()
plt.show()
