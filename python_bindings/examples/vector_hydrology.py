"""Prepare mapped hydrography before building a flow network.

Authors
-------
B. Gailleton
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj import CRS
from shapely.geometry import LineString, Polygon

import lsdtt3


crs = "EPSG:32633"
areas = gpd.GeoDataFrame(geometry=[Polygon([(40, 30), (160, 30), (160, 70), (40, 70)])], crs=crs)
tributaries = gpd.GeoDataFrame(
    geometry=[LineString([(0, 50), (40, 50)]), LineString([(160, 50), (200, 50)])], crs=crs
)
water = lsdtt3.water_centrelines(areas, tributaries, cell_size_metres=2.0).lines

# Repair a short mapped gap, then combine it with the polygon-derived channels.
broken = gpd.GeoDataFrame(geometry=[LineString([(0, 10), (30, 10)]), LineString([(38, 10), (80, 10)])], crs=crs)
roads = gpd.GeoDataFrame(geometry=[LineString([(34, 0), (34, 20)])], crs=crs)
repaired = lsdtt3.repair_water_lines(broken, roads, max_gap_length_metres=12.0).lines
channels = gpd.GeoDataFrame(pd.concat([water, repaired], ignore_index=True), crs=crs)

meta = lsdtt3.RasterMeta.from_shape(20, 20, cell_size=10.0, crs_wkt=CRS.from_user_input(crs).to_wkt(), epsg_code=32633)
z = np.tile(np.linspace(250, 50, 20, dtype=np.float32), (20, 1))
dem = lsdtt3.Raster.from_numpy(z, np.ones_like(z, dtype=bool), meta)
burned = lsdtt3.burn_streams(dem, channels)

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
extent = [meta.x_min, meta.x_max, meta.y_min, meta.y_max]
for axis, raster, title in zip(axes, (dem, burned.raster), ("Original DEM", "Conditioned DEM")):
    axis.imshow(raster.to_numpy(), extent=extent, origin="upper", cmap="terrain")
    channels.plot(ax=axis, color="deepskyblue", linewidth=2)
    axis.set(title=title, xlabel="x (m)", ylabel="y (m)")
plt.tight_layout()
plt.show()
