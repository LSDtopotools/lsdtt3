"""Extract and measure a synthetic tidal-channel network.

Authors
-------
B. Gailleton
"""

import numpy as np
import matplotlib.pyplot as plt

import lsdtt3


# Replace this synthetic case with a projected lidar DEM in real analyses.  The
# embayment has a sinuous tidal trunk plus branching creeks, but stays fully
# self-contained for a quick API demonstration.
rows, cols, cell_size = 240, 300, 2.0
meta = lsdtt3.RasterMeta.from_shape(rows, cols, cell_size=cell_size, epsg_code=32631)
y, x = np.mgrid[0:rows, 0:cols]
z = 1.7 + 0.003 * (rows - y) + 0.06 * np.sin(x / 22) * np.sin(y / 31)
channel = np.zeros((rows, cols), dtype=bool)


def burn_polyline(points, half_width_pixels):
    """Rasterise a smooth line into the synthetic channel mask."""
    for start, stop in zip(points[:-1], points[1:]):
        n = int(np.hypot(stop[0] - start[0], stop[1] - start[1])) * 3
        rr = np.linspace(start[0], stop[0], max(n, 2)).astype(int)
        cc = np.linspace(start[1], stop[1], max(n, 2)).astype(int)
        for row, col in zip(rr, cc):
            r0, r1 = max(0, row - half_width_pixels), min(rows, row + half_width_pixels + 1)
            c0, c1 = max(0, col - half_width_pixels), min(cols, col + half_width_pixels + 1)
            local_rows, local_cols = np.ogrid[r0:r1, c0:c1]
            channel[r0:r1, c0:c1] |= (
                (local_rows - row) ** 2 + (local_cols - col) ** 2
                <= half_width_pixels**2
            )


# Main tidal trunk reaches the seaward boundary at the bottom of the map.
trunk_rows = np.linspace(rows - 1, 28, 160)
trunk_cols = cols / 2 + 30 * np.sin((trunk_rows - 40) / 48)
burn_polyline(list(zip(trunk_rows, trunk_cols)), half_width_pixels=4)

# Tributary creeks join it from both marsh flanks.
for join_row in (55, 85, 120, 155, 188):
    join_col = cols / 2 + 30 * np.sin((join_row - 40) / 48)
    for direction in (-1, 1):
        t = np.linspace(0, 1, 50)
        creek_rows = join_row - 30 * (1 - t) + 12 * np.sin(np.pi * t)
        creek_cols = join_col + direction * (
            95 * (1 - t) + 14 * np.sin(np.pi * t)
        )
        burn_polyline(list(zip(creek_rows, creek_cols)), half_width_pixels=2)

# Channels are low relative to the salt-marsh platform.
z[channel] = 0.15 + 0.03 * np.sin(y[channel] / 13)
dem = lsdtt3.Raster.from_numpy(
    z.astype(np.float32), np.ones_like(z, dtype=bool), meta
)

channels = lsdtt3.find_tidal_channels(
    dem,
    use_flood_mask=True,
    flood_elevation_metres=0.5,
    extract_skeleton=True,
    compute_upl=True,
)
metrics = lsdtt3.analyse_tidal_skeleton(
    dem,
    channels.channel_mask,
    channels.skeleton,
    outlet_mode="edge_nodes",
    station_spacing_metres=10.0,
)

extent = (0, cols * cell_size, 0, rows * cell_size)
fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)

terrain = axes[0].imshow(dem.to_numpy(), cmap="terrain", extent=extent, origin="upper")
channel_overlay = np.ma.masked_where(channels.channel_mask.to_numpy() < 0.5, channels.channel_mask.to_numpy())
axes[0].imshow(channel_overlay, cmap="Blues", alpha=0.65, extent=extent, origin="upper")
axes[0].contour(channels.skeleton.to_numpy(), [0.5], colors="white", linewidths=0.7, extent=extent, origin="upper")
axes[0].set(title="Synthetic tidal embayment", xlabel="Easting (m)", ylabel="Northing (m)")
fig.colorbar(terrain, ax=axes[0], label="Elevation (m)", shrink=0.8)

upl = axes[1].imshow(channels.upl.to_numpy(), cmap="magma_r", extent=extent, origin="upper")
metrics.sections.plot(ax=axes[1], color="deepskyblue", linewidth=1.8)
metrics.stations.plot(ax=axes[1], color="white", edgecolor="black", markersize=10, zorder=3)
axes[1].set(title="Tidal sections, order, and stations", xlabel="Easting (m)", ylabel="Northing (m)")
fig.colorbar(upl, ax=axes[1], label="Unchannelled path length (m)", shrink=0.8)
plt.show()
