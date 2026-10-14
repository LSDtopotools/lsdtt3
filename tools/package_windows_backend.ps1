# Packages a minimal, redistributable Windows backend bundle for QGIS plugin use.
#
# The script stages selected lsdtt executables, discovers and copies dependent
# DLLs from vcpkg, copies GDAL/PROJ data folders, writes metadata files, and
# creates a final zip artifact.
param(
    [Parameter(Mandatory = $true)]
    [string]$Version,

    [string]$BuildDir = "build/windows-qgis",
    [string]$OutDir = "dist",
    [string]$VcpkgRoot = "$env:USERPROFILE/tools/vcpkg",
    [string]$Triplet = "x64-windows",
    [string]$PluginVersionRange = ">=0.1.0 <0.5.0",

    # CLI set to package: qgis / core (from tools/cli_profiles.txt) or all
    # (every lsdtt-*.exe in the build directory, excluding lsdtt-test-*).
    # Windows bundles feed the QGIS plugin, so the default is qgis.
    [ValidateSet("qgis", "core", "all")]
    [string]$CliProfile = "qgis",

    # Explicit list; overrides -CliProfile (bundle is labelled 'custom').
    [string[]]$CliNames = @(),

    # Staged when present in the build, skipped with a warning otherwise
    # (lsdtt-spectral-analysis is only built when FFTW3 is available).
    [string[]]$OptionalCliNames = @(
        "lsdtt-spectral-analysis"
    )
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# ---------------------------------------------------------------------------
# Resolve the CLI list: explicit -CliNames wins; otherwise -CliProfile
# decides. qgis/core come from tools/cli_profiles.txt; all = every
# lsdtt-*.exe in the build directory (excluding lsdtt-test-*).
# ---------------------------------------------------------------------------
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRootEarly = (Resolve-Path (Join-Path $scriptDir "..")).Path
$cmakeLists = Join-Path $repoRootEarly "CMakeLists.txt"

# The version is typed by hand here, unlike the Linux path where
# make_release.sh derives and checks it. Mislabelling a bundle is silent and
# permanent once published, so check it against the single source of truth.
$cmakeVersion = (Select-String -Path $cmakeLists -Pattern '^project\(lsdtt3 VERSION ([0-9.]+)' |
    Select-Object -First 1).Matches.Groups[1].Value
if (-not $cmakeVersion) {
    throw "Could not read project version from $cmakeLists"
}
$versionStem = ($Version -split '-')[0]
if ($Version -eq $versionStem) {
    if ($Version -ne $cmakeVersion) {
        throw ("Packaging version $Version but CMakeLists.txt says $cmakeVersion. " +
               "Bump project(lsdtt3 VERSION ...) or pass the matching -Version.")
    }
}
elseif ($versionStem -ne $cmakeVersion) {
    Write-Warning "Pre-release stem $versionStem != CMake version $cmakeVersion (continuing)"
}

# The CLIs that currently exist are declared in CMakeLists.txt. A long-lived
# build tree also holds binaries for tools since renamed or removed; they look
# identical by name, so packaging one ships a tool with no source behind it.
$knownClis = @(Select-String -Path $cmakeLists -Pattern '^\s*lsdtt_cli\((lsdtt-[A-Za-z0-9_-]+)' |
    ForEach-Object { $_.Matches.Groups[1].Value })

if ($CliNames.Count -gt 0) {
    $CliProfile = "custom"
}
elseif ($CliProfile -eq "all") {
    if ($knownClis.Count -eq 0) {
        throw "Could not read any lsdtt_cli() targets from $cmakeLists"
    }
    $exeDirs = @($BuildDir, (Join-Path $BuildDir "bin")) | Where-Object { Test-Path $_ }
    $discovered = $exeDirs |
        ForEach-Object { Get-ChildItem -Path $_ -Filter "lsdtt-*.exe" -File -ErrorAction SilentlyContinue } |
        Where-Object { $_.BaseName -notlike "lsdtt-test-*" } |
        ForEach-Object { $_.BaseName } |
        Sort-Object -Unique
    $stale = @($discovered | Where-Object { $knownClis -notcontains $_ })
    $CliNames = @($discovered | Where-Object { $knownClis -contains $_ })
    if ($stale.Count -gt 0) {
        Write-Warning ("Skipping $($stale.Count) stale binaries in $BuildDir no longer " +
                       "lsdtt_cli() targets in CMakeLists.txt: " + ($stale -join ", ") +
                       ". Delete them, or build into a clean directory, to silence this.")
    }
    if ($CliNames.Count -eq 0) {
        throw "CliProfile 'all' found no current lsdtt-*.exe in $BuildDir"
    }
    Write-Host "Profile 'all': packaging $($CliNames.Count) CLIs discovered in the build directory"
}
else {
    $profilesFile = Join-Path $scriptDir "cli_profiles.txt"
    if (-not (Test-Path $profilesFile)) {
        throw "profiles file not found: $profilesFile"
    }
    $CliNames = Get-Content $profilesFile |
        Where-Object { $_ -notmatch '^\s*(#|$)' } |
        ForEach-Object { $parts = $_ -split '\s+'; if ($parts[0] -eq $CliProfile) { $parts[1] } } |
        Where-Object { $_ }
    if ($CliNames.Count -eq 0) {
        throw "profile '$CliProfile' has no entries in $profilesFile"
    }
}

# Resolve a CLI executable from either:
# - <build_dir>/<cli>.exe
# - <build_dir>/bin/<cli>.exe
function Resolve-ExePath {
    param(
        [string]$BaseDir,
        [string]$CliName
    )

    $candidates = @(
        (Join-Path $BaseDir ("{0}.exe" -f $CliName)),
        (Join-Path $BaseDir (Join-Path "bin" ("{0}.exe" -f $CliName)))
    )

    foreach ($candidate in $candidates) {
        if (Test-Path $candidate) {
            return (Resolve-Path $candidate).Path
        }
    }

    throw "Could not find $CliName.exe in '$BaseDir' (checked root and bin/)."
}

# Uses dumpbin to list direct DLL dependencies of a PE binary (.exe or .dll).
# This intentionally requires a VS developer shell where dumpbin.exe is present.
function Get-PeDependents {
    param([string]$BinaryPath)

    $dumpbin = Get-Command dumpbin.exe -ErrorAction SilentlyContinue
    if (-not $dumpbin) {
        throw "dumpbin.exe not found. Run this script from a Visual Studio developer shell."
    }

    $output = & $dumpbin.Path /nologo /dependents $BinaryPath
    $deps = @()
    foreach ($line in $output) {
        if ($line -match '^\s+([A-Za-z0-9._-]+\.dll)\s*$') {
            $deps += $matches[1].ToLowerInvariant()
        }
    }
    return ($deps | Sort-Object -Unique)
}

# True for DLLs the OS itself provides on every supported Windows (10+):
# core system DLLs, API sets, and the universal CRT. These must NOT be
# bundled; everything else a staged binary imports must ship in bin/.
function Test-IsSystemDll {
    param([string]$DllName)

    if ($DllName -like "api-ms-win-*") { return $true }
    if ($DllName -like "ext-ms-*") { return $true }

    $systemDlls = @(
        "ucrtbase.dll",
        "kernel32.dll", "kernelbase.dll", "user32.dll", "gdi32.dll",
        "advapi32.dll", "shell32.dll", "shlwapi.dll", "ole32.dll",
        "oleaut32.dll", "comdlg32.dll", "comctl32.dll", "uxtheme.dll",
        "ws2_32.dll", "wsock32.dll", "iphlpapi.dll", "winhttp.dll",
        "wininet.dll", "urlmon.dll", "normaliz.dll", "dnsapi.dll",
        "crypt32.dll", "bcrypt.dll", "ncrypt.dll", "secur32.dll",
        "sspicli.dll", "wldap32.dll", "version.dll", "winmm.dll",
        "setupapi.dll", "cfgmgr32.dll", "userenv.dll", "psapi.dll",
        "dbghelp.dll", "imm32.dll", "msimg32.dll", "netapi32.dll",
        "ntdll.dll", "rpcrt4.dll", "msvcrt.dll"
    )
    return ($systemDlls -contains $DllName)
}

# Copies transitive runtime DLL dependencies for the selected executables.
#
# Strategy:
# 1) Seed queue with direct exe dependencies.
# 2) For each DLL found in vcpkg bin, copy it into stage/bin.
# 3) Parse that DLL's dependencies and continue until closure.
#
# Only DLLs that exist in vcpkg's triplet bin directory are copied here;
# the MSVC CRT is staged separately and self-containment is verified after.
function Copy-RuntimeDllClosure {
    param(
        [string[]]$ExePaths,
        [string]$VcpkgBin,
        [string]$TargetBin
    )

    $queue = New-Object System.Collections.Generic.Queue[string]
    $seen = New-Object System.Collections.Generic.HashSet[string]

    foreach ($exe in $ExePaths) {
        foreach ($dll in (Get-PeDependents -BinaryPath $exe)) {
            if ($seen.Add($dll)) {
                $queue.Enqueue($dll)
            }
        }
    }

    while ($queue.Count -gt 0) {
        $dllName = $queue.Dequeue()
        $source = Join-Path $VcpkgBin $dllName
        if (-not (Test-Path $source)) {
            continue
        }

        $dest = Join-Path $TargetBin $dllName
        Copy-Item -Path $source -Destination $dest -Force

        foreach ($dep in (Get-PeDependents -BinaryPath $source)) {
            if ($seen.Add($dep)) {
                $queue.Enqueue($dep)
            }
        }
    }
}

# Stages the MSVC C/C++ runtime (vcruntime140, msvcp140, ...) next to the
# executables. Users must not need the VC++ Redistributable preinstalled:
# without these DLLs the tools fail on a clean machine with
# "VCRUNTIME140.dll was not found". The redist folder is located via the
# VS developer-shell environment the script already requires for dumpbin.
function Copy-MsvcRuntime {
    param([string]$TargetBin)

    $redistRoot = $env:VCToolsRedistDir
    if (-not $redistRoot) {
        throw "VCToolsRedistDir is not set. Run this script from a Visual Studio developer shell so the MSVC runtime can be bundled."
    }

    $crtDir = Get-ChildItem -Path (Join-Path $redistRoot "x64") -Directory -Filter "Microsoft.VC*.CRT" |
        Select-Object -First 1
    if (-not $crtDir) {
        throw "MSVC CRT redist folder not found under $redistRoot\x64"
    }

    foreach ($dll in (Get-ChildItem $crtDir.FullName -Filter *.dll)) {
        Copy-Item -Path $dll.FullName -Destination $TargetBin -Force
    }
}

# Verifies the staged bundle is self-contained: every import of every staged
# binary must resolve inside stage/bin or be an OS-provided DLL. Fails the
# packaging rather than ship a bundle that breaks on user machines.
function Assert-SelfContained {
    param([string]$TargetBin)

    $missing = @()
    foreach ($binary in (Get-ChildItem $TargetBin -File | Where-Object { $_.Extension -in ".exe", ".dll" })) {
        foreach ($dep in (Get-PeDependents -BinaryPath $binary.FullName)) {
            if (Test-IsSystemDll -DllName $dep) { continue }
            if (-not (Test-Path (Join-Path $TargetBin $dep))) {
                $missing += "$($binary.Name) needs '$dep' but it is not staged in bin/"
            }
        }
    }

    if ($missing.Count -gt 0) {
        $missing | Sort-Object -Unique | ForEach-Object { Write-Error -ErrorAction Continue $_ }
        throw "Bundle is not self-contained; aborting."
    }
}

# Resolve key input/output locations from script parameters.
# Use repository root (parent of this tools directory) so the script works even
# when launched from inside tools/ or any other working directory.
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$buildAbs = [System.IO.Path]::GetFullPath((Join-Path $repoRoot $BuildDir))
$outAbs = [System.IO.Path]::GetFullPath((Join-Path $repoRoot $OutDir))
$vcpkgAbs = [System.IO.Path]::GetFullPath($VcpkgRoot)

$vcpkgBin = Join-Path $vcpkgAbs (Join-Path "installed" (Join-Path $Triplet "bin"))
$gdalData = Join-Path $vcpkgAbs (Join-Path "installed" (Join-Path $Triplet "share/gdal"))
$projData = Join-Path $vcpkgAbs (Join-Path "installed" (Join-Path $Triplet "share/proj"))

# Preflight checks fail early with actionable messages.
if (-not (Test-Path $buildAbs)) {
    throw "Build directory not found: $buildAbs"
}
if (-not (Test-Path $vcpkgBin)) {
    throw "vcpkg runtime bin not found: $vcpkgBin"
}
if (-not (Test-Path $gdalData)) {
    throw "GDAL data directory not found: $gdalData"
}
if (-not (Test-Path $projData)) {
    throw "PROJ data directory not found: $projData"
}

$bundleName = "lsdtt3-backend-win64-$CliProfile-v$Version"
$stage = Join-Path $outAbs $bundleName
$stageBin = Join-Path $stage "bin"
$stageShare = Join-Path $stage "share"

# Recreate clean staging tree.
if (Test-Path $stage) {
    Remove-Item -Recurse -Force $stage
}

New-Item -ItemType Directory -Path $stageBin -Force | Out-Null
New-Item -ItemType Directory -Path $stageShare -Force | Out-Null

# Copy selected CLI executables into stage/bin.
$exePaths = @()
$stagedCliNames = @()
foreach ($cli in $CliNames) {
    $exe = Resolve-ExePath -BaseDir $buildAbs -CliName $cli
    Copy-Item -Path $exe -Destination $stageBin -Force
    $exePaths += $exe
    $stagedCliNames += $cli
}
foreach ($cli in $OptionalCliNames) {
    if ($stagedCliNames -contains $cli) { continue }
    try {
        $exe = Resolve-ExePath -BaseDir $buildAbs -CliName $cli
    } catch {
        Write-Warning "Optional CLI $cli.exe not found in build - skipping (build with FFTW3 to include it)."
        continue
    }
    Copy-Item -Path $exe -Destination $stageBin -Force
    $exePaths += $exe
    $stagedCliNames += $cli
}

# Populate runtime closure and required data directories.
Copy-RuntimeDllClosure -ExePaths $exePaths -VcpkgBin $vcpkgBin -TargetBin $stageBin
Copy-MsvcRuntime -TargetBin $stageBin
Assert-SelfContained -TargetBin $stageBin
Copy-Item -Path $gdalData -Destination (Join-Path $stageShare "gdal") -Recurse -Force
Copy-Item -Path $projData -Destination (Join-Path $stageShare "proj") -Recurse -Force

# lsdtt3 bundled runtime data. The cosmo tools resolve crn_parameters.nc as
# <exe>\..\share\lsdtt3, the same bin/../share contract as the GDAL and PROJ
# data above, so it must be staged here or those tools cannot run without a
# download. Missing data is an error when a cosmo CLI is in the bundle; for a
# profile with none (the qgis default) it is only a warning.
$crnSource = Join-Path $repoRoot "data/crn_parameters.nc"
$cosmoStaged = @($CliNames | Where-Object { $_ -like "lsdtt-cosmo-*" })
if (Test-Path $crnSource) {
    $stageLsdtt = Join-Path $stageShare "lsdtt3"
    New-Item -ItemType Directory -Path $stageLsdtt -Force | Out-Null
    Copy-Item -Path $crnSource -Destination $stageLsdtt -Force
    Write-Host "Staged crn_parameters.nc into share/lsdtt3"
}
elseif ($cosmoStaged.Count -gt 0) {
    throw ("data/crn_parameters.nc not found, but this bundle ships " +
           ($cosmoStaged -join ", ") +
           " which need it at bin/../share/lsdtt3")
}
else {
    Write-Warning "data/crn_parameters.nc not found; bundle ships no cosmo CLIs, continuing"
}

$licensePath = Join-Path $repoRoot "LICENSE"
if (Test-Path $licensePath) {
    Copy-Item -Path $licensePath -Destination (Join-Path $stage "LICENSE") -Force
}

# Record lightweight bundle contents for traceability.
$thirdParty = @(
    "This backend bundle includes runtime components from vcpkg package builds",
    "and the Microsoft Visual C++ redistributable runtime DLLs.",
    "Triplet: $Triplet",
    "Generated: $(Get-Date -Format o)",
    "",
    "Included executables:"
)
$thirdParty += ($stagedCliNames | ForEach-Object { "- $_.exe" })
$thirdParty += ""
$thirdParty += "Included runtime DLLs:"
$thirdParty += (Get-ChildItem $stageBin -Filter *.dll | Sort-Object Name | ForEach-Object { "- $($_.Name)" })
Set-Content -Path (Join-Path $stage "THIRD_PARTY_NOTICES.txt") -Value $thirdParty -Encoding utf8

# setenv.ps1 makes the bundle runnable from an arbitrary folder.
$envScript = @(
    "`$ErrorActionPreference = 'Stop'",
    "`$root = Split-Path -Parent `$MyInvocation.MyCommand.Path",
    "`$env:PATH = `"`$(Join-Path `$root 'bin');`$env:PATH`"",
    "`$env:GDAL_DATA = Join-Path `$root 'share/gdal'",
    "`$env:PROJ_DATA = Join-Path `$root 'share/proj'",
    "`$env:PROJ_LIB = `$env:PROJ_DATA"
)
Set-Content -Path (Join-Path $stage "setenv.ps1") -Value $envScript -Encoding utf8

# End-user README: audience is often non-technical, so cover the full path
# from download to a working command. Single-quoted here-string, then
# placeholder replacement, so $env: examples survive verbatim.
$readme = @'
# lsdtt3 backend for Windows (v@VERSION@)

Self-contained LSDTopoTools command-line tools for 64-bit Windows. All
required DLLs (including the Microsoft Visual C++ runtime) and the GDAL/PROJ
support data are included -- there is nothing to install and no admin rights
are needed. Just unzip and run.

## Install and run

1. Download `@BUNDLE_NAME@.zip`.
2. Right-click the zip file and choose **Extract All...** (any folder is
   fine, for example `C:\lsdtt3`).
3. Open PowerShell and go to the extracted folder:

   ```powershell
   cd C:\lsdtt3\@BUNDLE_NAME@
   ```

4. Activate the tools for the current PowerShell session:

   ```powershell
   . .\setenv.ps1
   ```

   The leading dot-space matters. If Windows refuses to run the script
   ("running scripts is disabled"), first run:

   ```powershell
   Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
   ```

5. Check that it works:

   ```powershell
   lsdtt-raster-preprocessing --help
   ```

You need to repeat step 4 in every new PowerShell window.

## How the tools work: parameter files

Every `lsdtt-*` tool reads its options from a plain-text **parameter file**
that you write. The tools do NOT take analysis options as command-line flags
(something like `--input my_dem.tif` is rejected). The only flags are `-h`
(help), `-v` (version) and `-c` (citation).

```powershell
lsdtt-surface-metrics                                     # no arguments: prints how to run it and the file format
lsdtt-surface-metrics -h                                  # list every parameter, with type and default
lsdtt-surface-metrics -h write_hillshade_raster           # full guidance on one parameter
lsdtt-surface-metrics my_params.driver                    # run; parameter file in the current folder
lsdtt-surface-metrics C:\data\my_project my_params.driver # run; parameter file in C:\data\my_project
```

The `-h` output is the complete reference for each tool: the parameter names
it lists are exactly the keys you put in the parameter file. `-h html` writes
the same reference as HTML and CSV files, and `--help-json` prints it as JSON.

### Parameter file format

- One `key: value` per line. The key is a parameter name from `-h`.
- Keys are case-sensitive. Any parameter you leave out takes its default.
- A value must not contain spaces. Lists are comma-separated with no
  spaces (`a,b,c`), and paths must not contain spaces.
- Booleans are `true` or `false`.
- `#` starts a comment; blank lines are ignored.
- File I/O keys, accepted by every tool:
  - `read_path`: folder holding the input files (default: the folder
    containing the parameter file)
  - `read_fname`: input raster, e.g. `my_dem.tif` (any GDAL-readable format)
  - other input files (any key ending in `_fname`) are read from `read_path`
    too. A relative path such as `inputs\outlets.csv` is relative to
    `read_path`; an absolute path is used as given.
  - `write_path`: folder for outputs (default: `read_path`)
  - `read_path` and `write_path` are optional and usually left out, so that
    everything is read from and written to the parameter file's folder. If
    you do set one, the folder must already exist, or the tool stops with
    an error.
  - `write_fname`: prefix for output files (default: the stem of `read_fname`)
- Outputs are written as `<write_fname>_<descriptor>.<ext>` in `write_path`,
  together with `<write_fname>_ingested_params.param` (every parameter value
  actually used, including defaults) and `<write_fname>_README.html`.
- Keys a tool does not recognise are reported as warnings and ignored, so
  check the run output for typos.

### Example

Write a hillshade of `C:\data\my_project\my_dem.tif`. Put the parameter file
in the same folder as the DEM and leave out `read_path` and `write_path`. Save
this as `C:\data\my_project\my_params.driver` (in Notepad, choose "All files"
so it is not saved as `.txt`):

```text
# my_params.driver -- hillshade a DEM
read_fname: my_dem.tif
write_fname: my_dem
write_hillshade_raster: true
```

Then, in a PowerShell window where `setenv.ps1` has been run:

```powershell
cd C:\data\my_project
lsdtt-surface-metrics my_params.driver
```

This writes `C:\data\my_project\my_dem_hillshade.tif`. From another folder,
run `lsdtt-surface-metrics C:\data\my_project my_params.driver` instead. To
send outputs elsewhere, add a `write_path:` line naming an existing folder. For other
analyses, run the relevant tool with `-h` to find its parameters, and add
them to the parameter file.

## What is in this bundle

- `bin\` -- the lsdtt command-line tools and their DLLs:
@TOOL_LIST@
- `share\gdal`, `share\proj` -- support data needed at runtime
- `share\lsdtt3` -- bundled lsdtt3 data (cosmogenic nuclide parameters)
- `setenv.ps1` -- sets up PATH and data paths for the current session
- `manifest.json` -- version info and SHA256 checksums
'@

$toolListMd = ($stagedCliNames | ForEach-Object { "  - ``$_.exe``" }) -join "`n"
$readme = $readme.Replace('@VERSION@', $Version)
$readme = $readme.Replace('@BUNDLE_NAME@', $bundleName)
$readme = $readme.Replace('@TOOL_LIST@', $toolListMd)
[System.IO.File]::WriteAllText(
    (Join-Path $stage "README.md"),
    $readme + "`n",
    (New-Object System.Text.UTF8Encoding $false))

# Build manifest with metadata and per-file SHA256 checksums.
$manifest = [ordered]@{
    backend_version = $Version
    profile = "$CliProfile-windows"
    architecture = "win64"
    plugin_version_range = $PluginVersionRange
    generated_utc = (Get-Date).ToUniversalTime().ToString("o")
    tools = @()
    checksums_sha256 = [ordered]@{}
}

$filesForHash = @()
$filesForHash += Get-Item (Join-Path $stage "README.md")
$filesForHash += Get-ChildItem $stageBin -File
$filesForHash += Get-ChildItem (Join-Path $stageShare "gdal") -Recurse -File
$filesForHash += Get-ChildItem (Join-Path $stageShare "proj") -Recurse -File
if (Test-Path (Join-Path $stageShare "lsdtt3")) {
    $filesForHash += Get-ChildItem (Join-Path $stageShare "lsdtt3") -Recurse -File
}
foreach ($f in $filesForHash) {
    $rel = $f.FullName.Substring($stage.Length + 1).Replace('\', '/')
    $hash = (Get-FileHash -Path $f.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    $manifest.checksums_sha256[$rel] = $hash
}

foreach ($cli in $stagedCliNames) {
    $manifest.tools += "$cli.exe"
}

# Write without a BOM: PowerShell 5.1's utf8 encoding prepends one, which
# strict JSON parsers (e.g. Python json.load) reject.
$manifestJson = $manifest | ConvertTo-Json -Depth 8
[System.IO.File]::WriteAllText(
    (Join-Path $stage "manifest.json"),
    $manifestJson + "`n",
    (New-Object System.Text.UTF8Encoding $false))

# Zip staged contents into the distributable backend artifact.
if (-not (Test-Path $outAbs)) {
    New-Item -ItemType Directory -Path $outAbs -Force | Out-Null
}

$zipPath = Join-Path $outAbs ("$bundleName.zip")
if (Test-Path $zipPath) {
    Remove-Item -Force $zipPath
}
# Use .NET ZipFile rather than Compress-Archive: it includes the top-level
# bundle folder (matching the Linux tar.gz layout, so unzipping never
# scatters files into the current directory) and writes forward-slash
# entry names, which non-Windows unzip tools require.
Add-Type -AssemblyName System.IO.Compression.FileSystem
[System.IO.Compression.ZipFile]::CreateFromDirectory(
    $stage, $zipPath,
    [System.IO.Compression.CompressionLevel]::Optimal,
    $true)

Write-Host "Created backend bundle: $zipPath"
Write-Host "Staging directory: $stage"
