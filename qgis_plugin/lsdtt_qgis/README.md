# lsdtt_qgis

QGIS Processing provider for selected lsdtt3 command-line workflows.

Current status:

- plugin scaffold and provider registration
- backend path manager and health check helpers
- MVP wrappers for:
  - `lsdtt-raster-preprocessing`
  - `lsdtt-channel-extraction`
  - `lsdtt-chi-analysis`
  - `lsdtt-basin-extraction`

Backend requirements:

- set backend path in plugin settings to a folder containing `lsdtt-*` executables
- for Windows backend bundles, use the packaged `bin/` and `share/` layout
