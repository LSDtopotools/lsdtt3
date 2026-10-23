"""Steady-state 10Be prediction and D8 concentration accumulation."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


def as_nan(raster):
    return np.where(raster.mask_numpy(), raster.to_numpy(), np.nan)


# Fetch a small 30 m Copernicus DEM of the Green River canyon beside
# Canyonlands National Park, Utah (projected to UTM zone 12N). The tile is
# cached in dem_cache/ so later runs do not download it again.
dem = lsdtt3.fetch_dem(
    west=-110.10, south=38.36, east=-109.96, north=38.47,
    grid_spacing=30.0, cache_dir="dem_cache",
)
dem = lsdtt3.fill_pits(dem)
flow = lsdtt3.FlowInfo(dem)

# ``standard`` keeps this small example independent of crn_parameters.nc.
# Use a supplied crn_parameters.nc path for ERA40/NCEP, paleomagnetism, or LSDn.
production = lsdtt3.crn_production_scaling(dem, atmosphere="standard")
erosion = lsdtt3.crn_uniform_raster(dem, 0.01)  # effective g cm^-2 yr^-1
be10 = lsdtt3.crn_steady_state(erosion, production, nuclide="Be10")
be10_at_outlets = lsdtt3.crn_accumulated_concentration(be10, erosion, flow)

fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
for ax, raster, title, cmap in zip(
    axes,
    [production, be10, be10_at_outlets],
    ["production scaling", r"local $^{10}$Be", r"D8 accumulated $^{10}$Be"],
    ["viridis", "magma", "magma"],
):
    image = ax.imshow(as_nan(raster), cmap=cmap)
    ax.set_title(title)
    fig.colorbar(image, ax=ax, shrink=0.8)

plt.show()
