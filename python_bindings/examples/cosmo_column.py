"""Minimal one-dimensional cosmogenic-nuclide column workflow."""

import matplotlib.pyplot as plt
import numpy as np
import lsdtt3

# Point forward model: 10Be, 26Al, and in-situ 14C at the surface.
column = lsdtt3.CRNColumn(rock_density_kg_m3=2650.0)
steady = column.steady_state(
    effective_erosion_rate_g_cm2_yr=0.01,
    production_scaling=1.0,
)
print("surface 10Be:", steady.concentration_be10, "atoms/g")

# The history is youngest -> oldest. The final -9999 row supplies the
# long-term background erosion rate used to initialise the column.
column.set_erosion_history(
    time_yr_bp=[0.0, 10_000.0, -9999.0],
    erosion_rate_m_yr=[0.0002, 0.00005, 0.00005],
    removal_metres=[0.0, 0.0, 0.0],
)
profile = column.transient(np.linspace(0.0, 2.0, 101))

depth = np.array([sample.depth_m for sample in profile])
be10 = np.array([sample.concentration_be10 for sample in profile])
plt.plot(be10, depth)
plt.gca().invert_yaxis()
plt.xlabel(r"$^{10}$Be concentration (atoms g$^{-1}$)")
plt.ylabel("depth below present surface (m)")
plt.show()

# A one-nuclide apparent steady-state erosion-rate inversion.
print(lsdtt3.crn_apparent_erosion_rate("Be10", steady.concentration_be10))
