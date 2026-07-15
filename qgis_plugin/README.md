# qgis_plugin

Utilities and source for the `lsdtt_qgis` QGIS plugin.

## Package plugin zip

From the repository root:

```bash
python3 qgis_plugin/package_plugin.py
```

This writes an installable zip to:

- `qgis_plugin/dist/lsdtt3-<version>.zip`

where `<version>` is read from `qgis_plugin/lsdtt_qgis/metadata.txt`.

## Install in QGIS

In QGIS: `Plugins -> Manage and Install Plugins... -> Install from ZIP` and
select the generated zip.

## Backend command-line tools (required)

The plugin calls external `lsdtt-*` executables. Install or unpack a backend
bundle before running algorithms.

### Recommended (Windows): use packaged backend zip

1. Build/package backend (developer side):

```powershell
.\tools\package_windows_backend.ps1 -Version 0.2.0
```

2. Share and unzip `dist/lsdtt3-backend-win64-minimal-v0.2.0.zip`.
3. In a PowerShell session, load runtime env vars:

```powershell
cd C:\path\to\lsdtt3-backend-win64-minimal-v0.2.0
. .\setenv.ps1
.\bin\lsdtt-raster-preprocessing.exe --version
```

If `--version` works, the backend is ready.

### Point plugin to backend folder

In QGIS, open `lsdtt3 -> LSDTT Backend Settings` and select the unzipped
backend root folder (the one containing `bin/`, `share/`, and `manifest.json`).
