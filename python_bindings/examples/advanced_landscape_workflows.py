"""Cellular hillslopes, depression topology, and CRN-history products."""

import numpy as np

import lsdtt3


dem = lsdtt3.Raster.read("perso_tests/dem/green_river_1.tif")
channels = lsdtt3.crn_uniform_raster(dem, 0.0)

# Conservative cellular transport can also run independently of fluvial incision.
cellular = lsdtt3.hillslope_transport_cellular(
    dem, channels, duration_years=1_000.0, diffusivity=1.0e-2
)

# The stateful LEM accepts the same cellular process model.
model = lsdtt3.LandscapeModel(
    dem, timestep_years=1_000.0, hillslope_model="cellular",
    fluvial_method="dinf", coupling_regime="coupled", uplift_rate=0.0,
)
final_surface = model.run(10).elevation
print(model.coupling_regime, cellular.n_face_events)

# Inspect all nested depressions rather than only candidate sinks.
hierarchy = lsdtt3.depression_hierarchy(final_surface)
print(hierarchy.nodes[["depression_id", "parent_id", "volume_m3"]].head())

# One present-to-past erosion history drives Be-10, Al-26, and C-14 together.
column = lsdtt3.CRNColumn()
column.set_erosion_history(
    [20_000.0, 5_000.0, -9_999.0],
    [2.0e-4, 8.0e-4, 2.0e-4],
    [0.0, 0.0, 0.0],
)
production = lsdtt3.crn_uniform_raster(final_surface, 1.0)
concentrations = lsdtt3.crn_transient_history(column, production)
be10 = concentrations["be10"]
