# Python API reference

These pages are rendered from a Python-only documentation façade for the
public `lsdtt3` API. It gives native nanobind callables readable Python
signatures and NumPy-style descriptions without changing their C++
implementation or runtime behaviour. Use the [Python guides](../python/index.md)
for workflow context and complete scripts.

## Groups

- [Raster and terrain](raster.md)
- [Flow, channels, and basins](flow.md)
- [Chi and ksn](chi.md)
- [Fetching and vectors](fetch.md)
- [Cosmogenic nuclides](cosmogenic.md)
- [GraphFlood](graphflood.md)

Each documented entry identifies its public arguments, result type, notes, and
author. Native overloads remain authoritative at runtime; the façade describes
their supported public calling patterns.
