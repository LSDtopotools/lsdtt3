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

- set backend path to a versioned bundle with `backend_version >= 0.6.0` in
  `manifest.json`; older and unversioned backends are rejected because the
  wrappers emit the normalised parameter names (older backends do not
  recognise them) and the least-cost breach keys `breach_lc_edge_policy` and
  `breach_lc_minimize_dist`
- for Windows backend bundles, use the packaged `bin/` and `share/` layout
