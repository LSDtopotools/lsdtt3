# SPDX-License-Identifier: AGPL-3.0-or-later
"""Authors of the public lsdtt3 Python API, for the generated docstrings.

Every public export is credited to B. Gailleton, who wrote the bindings,
except those S. M. Mudd wrote himself, which list him first and B. Gailleton
second. Exports S. M. Mudd changed list him second. The lists were made on
2026-10-22 from ``git blame -w -M -C`` over each export's binding code: its
``.def``/``nb::class_`` statement, the C++ function it binds, and, for
exports defined in Python, their ``def``. An export is "changed" if any
non-blank line of that code is his (new functions and parameters, renamed
keywords, changed defaults, fixes such as the ``gdal_init()`` calls), and
"written" if every line is.

When you change an export, add yourself here. Removing a name needs a reason
recorded in the commit message, since it takes away a credit.

Authors
-------
S. M. Mudd
"""

DEFAULT_AUTHORS = ("B. Gailleton",)

# Bindings S. M. Mudd wrote himself: he is first author, B. Gailleton second.
WRITTEN_BY_MUDD = (
    "absolute",
    "coherence_diffusion",
    "coregister_dem",
    "gdal_data_description",
    "gdal_version",
    "perona_malik_smoothing",
    "proj_data_description",
    "proj_version",
    "stencil_difference_magnitude",
    "vertical_units_to_metres",
)

# Bindings by B. Gailleton that S. M. Mudd changed: he is second author.
CHANGED_BY_MUDD = (
    "CRNColumn",
    "FlowInfo",
    "LandscapeModel",
    "RasterMeta",
    "bounds_from_vector",
    "breach_least_cost",
    "crn_transient_history",
    "divide",
    "extract_channel_heads",
    "extract_valley_centrelines",
    "fetch_dem",
    "fetch_dem_for_outlets",
    "find_flat_surfaces",
    "find_internal_drainage",
    "find_road_flow_discrimination",
    "frangi_filter",
    "helper_for_llm",
    "least_cost_carve_depth",
    "mountain_front_outlets",
    "route_sediment",
    "run_graphflood",
    "run_graphflood_induced",
    "run_graphflood_multigrid",
    "run_graphflood_optimised",
    "run_graphflood_optimised_v2",
    "run_landslide_series",
    "spectral_analysis",
    "synthetic_network",
    "wiener_filter",
)

AUTHORS = {name: ("S. M. Mudd", "B. Gailleton") for name in WRITTEN_BY_MUDD}
AUTHORS.update({name: ("B. Gailleton", "S. M. Mudd") for name in CHANGED_BY_MUDD})


def authors_for(name):
    """Return the authors credited for the public export ``name``.

    Parameters
    ----------
    name : str
        Name of the export in ``lsdtt3.__all__``.

    Returns
    -------
    tuple of str
        The authors in order: ``("B. Gailleton",)`` unless the export is in
        :data:`AUTHORS`.

    Authors
    -------
    S. M. Mudd
    """
    return AUTHORS.get(name, DEFAULT_AUTHORS)
