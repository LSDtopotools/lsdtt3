# lsdtt_qgis

QGIS Processing provider for selected lsdtt3 command-line workflows.

Contents:

- `provider.py`: registers the Processing algorithms
- `algorithms/`: one module per algorithm (the five drainage algorithms
  share `drainage_analysis.py`), plus shared helpers in `support.py`
- `backend/`: backend path, manifest version checks, and `health.py`, whose
  `REQUIRED_CLIS` lists the eight CLIs the wrappers run
- `naming.py`: QGIS-free path and format helpers, unit-tested without QGIS

The algorithms and what they need are listed in `qgis_plugin/README.md`.

Backend requirements:

- set backend path to a versioned bundle with `backend_version >= 0.6.0` in
  `manifest.json`; older and unversioned backends are rejected because the
  wrappers emit the normalised parameter names (older backends do not
  recognise them) and the least-cost breach keys `breach_lc_edge_policy` and
  `breach_lc_minimize_dist`
- for Windows backend bundles, use the packaged `bin/` and `share/` layout
