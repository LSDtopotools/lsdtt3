"""Compare threshold and Wiener channel-head extraction on Green River."""

import matplotlib.pyplot as plt
import lsdtt3

dem = lsdtt3.fill_pits(lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif"))
flow = lsdtt3.FlowInfo(dem)
hillshade = lsdtt3.multidirectional_hillshade(dem).to_numpy(nodata_to_nan=True)

methods = [
    ("Threshold: 100 pixels", dict(method="threshold", threshold_pixels=100)),
    ("Threshold: 500 pixels", dict(method="threshold", threshold_pixels=500)),
    ("Wiener: 500 m²", dict(method="wiener", pruning_drainage_area_m2=500.0,
                              connected_components_threshold_pixels=10)),
    ("Wiener: 2,000 m²", dict(method="wiener", pruning_drainage_area_m2=2000.0,
                                connected_components_threshold_pixels=10)),
]

fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
for ax, (title, options) in zip(axes.flat, methods):
    heads = lsdtt3.extract_channel_heads(dem, flow, **options)
    network = flow.channel_network(flow.scope_from_heads(heads))
    channels = network.points_geodataframe()
    source_table = heads.table

    ax.imshow(hillshade, cmap="gray", origin="upper")
    ax.scatter(channels.col, channels.row, s=1, c="dodgerblue", alpha=0.7)
    ax.scatter(source_table.col, source_table.row, marker="x", s=42,
               c="crimson", linewidths=1.4)
    ax.set(title=f"{title} ({heads.n_heads} heads)", xticks=[], yticks=[])

plt.show()
