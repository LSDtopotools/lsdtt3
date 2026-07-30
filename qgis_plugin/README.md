# qgis_plugin

Utilities and source for the `lsdtt_qgis` QGIS plugin.

## Finding the latest backend version number

Use one of these from the repository root:

```bash
# latest git tag in this checkout (for example v0.2.2)
git describe --tags --abbrev=0

# latest GitHub release tag (requires gh auth)
gh release list --limit 1
```

Use the numeric part (for example `0.2.2`) as the `-Version` argument for
packaging scripts.

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
.\tools\package_windows_backend.ps1 -Version <version>
```

2. Share and unzip `dist/lsdtt3-backend-win64-qgis-v<version>.zip`.
3. In a PowerShell session, load runtime env vars:

```powershell
cd C:\path\to\lsdtt3-backend-win64-qgis-v<version>
. .\setenv.ps1
.\bin\lsdtt-raster-preprocessing.exe --version
```

If `--version` works, the backend is ready.

### Point plugin to backend folder

In QGIS, open `lsdtt3 -> LSDTT Backend Settings` and select the unzipped
backend root folder (the one containing `bin/`, `share/`, and `manifest.json`).
