# Packages a minimal, redistributable Windows backend bundle for QGIS plugin use.
#
# The script stages selected lsdtt executables, discovers and copies dependent
# DLLs from vcpkg, copies GDAL/PROJ data folders, writes metadata files, and
# creates a final zip artifact.
param(
    [Parameter(Mandatory = $true)]
    [string]$Version,

    [string]$BuildDir = "build/minimal-windows",
    [string]$OutDir = "dist",
    [string]$VcpkgRoot = "$env:USERPROFILE/tools/vcpkg",
    [string]$Triplet = "x64-windows",
    [string]$PluginVersionRange = ">=0.1.0 <0.3.0",

    [string[]]$CliNames = @(
        "lsdtt-raster-preprocessing",
        "lsdtt-channel-extraction",
        "lsdtt-chi-analysis",
        "lsdtt-basin-extraction"
    ),

    # Staged when present in the build, skipped with a warning otherwise
    # (lsdtt-spectral-analysis is only built when FFTW3 is available).
    [string[]]$OptionalCliNames = @(
        "lsdtt-spectral-analysis"
    )
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

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

# Copies transitive runtime DLL dependencies for the selected executables.
#
# Strategy:
# 1) Seed queue with direct exe dependencies.
# 2) For each DLL found in vcpkg bin, copy it into stage/bin.
# 3) Parse that DLL's dependencies and continue until closure.
#
# Only DLLs that exist in vcpkg's triplet bin directory are copied.
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

$bundleName = "lsdtt3-backend-win64-minimal-v$Version"
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
Copy-Item -Path $gdalData -Destination (Join-Path $stageShare "gdal") -Recurse -Force
Copy-Item -Path $projData -Destination (Join-Path $stageShare "proj") -Recurse -Force

$licensePath = Join-Path $repoRoot "LICENSE"
if (Test-Path $licensePath) {
    Copy-Item -Path $licensePath -Destination (Join-Path $stage "LICENSE") -Force
}

# Record lightweight bundle contents for traceability.
$thirdParty = @(
    "This backend bundle includes runtime components from vcpkg package builds.",
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

# Build manifest with metadata and per-file SHA256 checksums.
$manifest = [ordered]@{
    backend_version = $Version
    profile = "minimal-windows"
    architecture = "win64"
    plugin_version_range = $PluginVersionRange
    generated_utc = (Get-Date).ToUniversalTime().ToString("o")
    tools = @()
    checksums_sha256 = [ordered]@{}
}

$filesForHash = @()
$filesForHash += Get-ChildItem $stageBin -File
$filesForHash += Get-ChildItem (Join-Path $stageShare "gdal") -Recurse -File
$filesForHash += Get-ChildItem (Join-Path $stageShare "proj") -Recurse -File
foreach ($f in $filesForHash) {
    $rel = $f.FullName.Substring($stage.Length + 1).Replace('\\', '/')
    $hash = (Get-FileHash -Path $f.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    $manifest.checksums_sha256[$rel] = $hash
}

foreach ($cli in $stagedCliNames) {
    $manifest.tools += "$cli.exe"
}

$manifestJson = $manifest | ConvertTo-Json -Depth 8
Set-Content -Path (Join-Path $stage "manifest.json") -Value $manifestJson -Encoding utf8

# Zip staged contents into the distributable backend artifact.
if (-not (Test-Path $outAbs)) {
    New-Item -ItemType Directory -Path $outAbs -Force | Out-Null
}

$zipPath = Join-Path $outAbs ("$bundleName.zip")
if (Test-Path $zipPath) {
    Remove-Item -Force $zipPath
}
Compress-Archive -Path (Join-Path $stage "*") -DestinationPath $zipPath -CompressionLevel Optimal

Write-Host "Created backend bundle: $zipPath"
Write-Host "Staging directory: $stage"
