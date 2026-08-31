"""Route mixed bedload over a Green River channel network."""

import matplotlib.pyplot as plt
import numpy as np

import lsdtt3


def as_nan(raster):
    return np.where(raster.mask_numpy(), raster.to_numpy(), np.nan)


dem = lsdtt3.fill_pits(lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif"))
flow = lsdtt3.FlowInfo(dem)
scope = flow.channel_scope(max(100, dem.valid_pixel_count // 250))
network = flow.channel_network(scope)

result = lsdtt3.route_sediment(
    network,
    dem,
    formula="wilcock_crowe",
    n_timesteps=30,
    timestep_seconds=86_400.0,
    grain_sizes_mm=[4.0, 8.0, 16.0, 32.0, 64.0],
    grain_fractions=[0.10, 0.25, 0.35, 0.20, 0.10],
)

fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
axes[0].imshow(as_nan(dem), cmap="terrain")
axes[0].contour(network.channel_mask.to_numpy(), levels=[0.5],
                colors="dodgerblue", linewidths=0.4)
axes[0].set_title("Channel network")

deposition = axes[1].imshow(as_nan(result.deposition), cmap="magma")
axes[1].set_title("Cumulative deposition")
fig.colorbar(deposition, ax=axes[1], label="m³ per reach")

d50 = axes[2].imshow(as_nan(result.active_layer_d50) * 1_000.0, cmap="viridis")
axes[2].set_title("Final active-layer D50")
fig.colorbar(d50, ax=axes[2], label="mm")

for ax in axes:
    ax.set_axis_off()

print(f"outlet flux: {result.total_flux_out_m3:.3e} m³")
print(result.reaches[["junction_id", "slope", "deposited_m3", "mobilised_m3"]].head())
plt.show()
