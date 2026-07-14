# tools

Utility scripts for building and packaging developer artifacts.

## `package_windows_backend.ps1`

Creates a redistributable Windows backend bundle for the QGIS-plugin-oriented
minimal profile.

### What it does

- copies selected `lsdtt-*.exe` tools from the Windows build directory
- resolves and copies dependent runtime DLLs from vcpkg
- copies GDAL and PROJ data folders (`share/gdal`, `share/proj`)
- writes `manifest.json` with SHA256 checksums
- writes `setenv.ps1` to configure runtime environment variables
- builds `dist/lsdtt3-backend-win64-minimal-v<version>.zip`

### Requirements

- run from repository root
- x64 Visual Studio developer shell (for `dumpbin.exe`)
- completed Windows build (default: `build/minimal-windows`)
- vcpkg install with GDAL/PROJ runtime files (default: `%USERPROFILE%/tools/vcpkg`)

### Basic usage

```powershell
.\tools\package_windows_backend.ps1 -Version 0.2.0
```

### Common options

```powershell
.\tools\package_windows_backend.ps1 `
  -Version 0.2.0 `
  -BuildDir build/minimal-windows `
  -OutDir dist
```

### Outputs

- zip artifact: `dist/lsdtt3-backend-win64-minimal-v<version>.zip`
- staging folder: `dist/lsdtt3-backend-win64-minimal-v<version>/`

## End-user install (Windows)

After you share the zip with a user, they can run the tools without installing
Visual Studio, vcpkg, or GDAL separately.

1. Unzip `lsdtt3-backend-win64-minimal-v<version>.zip` to a folder, for example:
   `C:\tools\lsdtt3-backend-win64-minimal-v<version>`
2. Open PowerShell in that unzipped folder.
3. Load runtime environment variables for the current session:

```powershell
. .\setenv.ps1
```

4. Run a CLI from `bin`:

```powershell
.\bin\lsdtt-raster-preprocessing.exe --version
```

Notes:

- The leading dot-space in `. .\setenv.ps1` is required (dot-sourcing).
- Users must run `setenv.ps1` in each new PowerShell session unless they add
  those environment variables to their PowerShell profile.
