# Cosmogenic erosion tools

The cosmogenic bindings cover three related scales:

- one-dimensional columns (`CRNColumn` and paired-nuclide diagnostics);
- raster production/concentration prediction;
- basin-scale inversion and outlet predictions.

For a basin workflow, first construct the same D8 flow/network context used by
the basin and chi tools, then use `crn_invert_basin_erosion` with a selected
basin and a measured concentration. The result retains inversion diagnostics
instead of writing them implicitly.

The runnable examples in `python_bindings/examples/` are a better starting
point than a synthetic one-line call:

- `cosmo_column.py`
- `cosmo_raster.py`
- `cosmo_basin.py`
- `cosmo_basin_carpathians.py`

See [Cosmogenic API](../python_api/cosmogenic.md) for the complete callable
reference and units.
