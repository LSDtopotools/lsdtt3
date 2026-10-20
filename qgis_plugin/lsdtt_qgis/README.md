# lsdtt_qgis

QGIS Processing provider for selected lsdtt3 command-line workflows.

Contents:

- `provider.py`: registers the Processing algorithms
- `algorithms/`: one module per algorithm (the five drainage algorithms
  share `drainage_analysis.py`), plus shared helpers in `support.py`
- `backend/`: backend path, manifest version checks, and `health.py`, whose
  `REQUIRED_CLIS` lists the nine CLIs the wrappers run
- `naming.py`: QGIS-free path and format helpers, unit-tested without QGIS

The algorithms and what they need are listed in `qgis_plugin/README.md`.

Backend requirements:

- set backend path to a versioned bundle whose `backend_version` in
  `manifest.json` has the plugin's major.minor version (plugin `X.Y.*` accepts
  backend `X.Y.*`); other series and unversioned backends are rejected,
  because the wrappers emit parameter names the backend must recognise
- for Windows backend bundles, use the packaged `bin/` and `share/` layout
