# Python API reference

These pages are generated directly from the built runtime API manifest, including
nanobind overloads, defaults, methods, and properties. Use the
[Python guides](../python/index.md) for workflow context and complete scripts.

## Groups

- [Raster and terrain](raster.md)
- [Flow, channels, and basins](flow.md)
- [Chi and ksn](chi.md)
- [Fetching and vectors](fetch.md)
- [Cosmogenic nuclides](cosmogenic.md)
- [GraphFlood](graphflood.md)
- [Landscape and hillslopes](landscape.md)
- [Drainage and depressions](drainage.md)
- [Sediment, landslides, and landforms](sediment.md)
- [Synthetic workflows](synthetic.md)
- [Fluvial, valley, marsh, and vector hydrology](fluvial.md)
- [Statistics and utilities](statistics.md)

Every public export has exactly one primary entry. Runtime nanobind and Python
signatures are authoritative; regeneration fails if the committed pages drift.

Developers changing bindings or signatures should read
[Python bindings architecture and API contracts](../developer_notes/python_bindings_architecture.md).
