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

- set backend path to a versioned bundle with `backend_version >= 0.5.0` in
  `manifest.json`; older and unversioned backends are rejected because the
  wrappers emit the canonical `threshold_contributing_pixels` key
- for Windows backend bundles, use the packaged `bin/` and `share/` layout
