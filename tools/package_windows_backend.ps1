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
    )
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

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

$repoRoot = (Get-Location).Path
$buildAbs = [System.IO.Path]::GetFullPath((Join-Path $repoRoot $BuildDir))
$outAbs = [System.IO.Path]::GetFullPath((Join-Path $repoRoot $OutDir))
$vcpkgAbs = [System.IO.Path]::GetFullPath($VcpkgRoot)

$vcpkgBin = Join-Path $vcpkgAbs (Join-Path "installed" (Join-Path $Triplet "bin"))
$gdalData = Join-Path $vcpkgAbs (Join-Path "installed" (Join-Path $Triplet "share/gdal"))
$projData = Join-Path $vcpkgAbs (Join-Path "installed" (Join-Path $Triplet "share/proj"))

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

if (Test-Path $stage) {
    Remove-Item -Recurse -Force $stage
}

New-Item -ItemType Directory -Path $stageBin -Force | Out-Null
New-Item -ItemType Directory -Path $stageShare -Force | Out-Null

$exePaths = @()
foreach ($cli in $CliNames) {
    $exe = Resolve-ExePath -BaseDir $buildAbs -CliName $cli
    Copy-Item -Path $exe -Destination $stageBin -Force
    $exePaths += $exe
}

Copy-RuntimeDllClosure -ExePaths $exePaths -VcpkgBin $vcpkgBin -TargetBin $stageBin
Copy-Item -Path $gdalData -Destination (Join-Path $stageShare "gdal") -Recurse -Force
Copy-Item -Path $projData -Destination (Join-Path $stageShare "proj") -Recurse -Force

$licensePath = Join-Path $repoRoot "LICENSE"
if (Test-Path $licensePath) {
    Copy-Item -Path $licensePath -Destination (Join-Path $stage "LICENSE") -Force
}

$thirdParty = @(
    "This backend bundle includes runtime components from vcpkg package builds.",
    "Triplet: $Triplet",
    "Generated: $(Get-Date -Format o)",
    "",
    "Included executables:"
)
$thirdParty += ($CliNames | ForEach-Object { "- $_.exe" })
$thirdParty += ""
$thirdParty += "Included runtime DLLs:"
$thirdParty += (Get-ChildItem $stageBin -Filter *.dll | Sort-Object Name | ForEach-Object { "- $($_.Name)" })
Set-Content -Path (Join-Path $stage "THIRD_PARTY_NOTICES.txt") -Value $thirdParty -Encoding utf8

$envScript = @(
    "`$ErrorActionPreference = 'Stop'",
    "`$root = Split-Path -Parent `$MyInvocation.MyCommand.Path",
    "`$env:PATH = `"`$(Join-Path `$root 'bin');`$env:PATH`"",
    "`$env:GDAL_DATA = Join-Path `$root 'share/gdal'",
    "`$env:PROJ_DATA = Join-Path `$root 'share/proj'",
    "`$env:PROJ_LIB = `$env:PROJ_DATA"
)
Set-Content -Path (Join-Path $stage "setenv.ps1") -Value $envScript -Encoding utf8

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

foreach ($cli in $CliNames) {
    $manifest.tools += "$cli.exe"
}

$manifestJson = $manifest | ConvertTo-Json -Depth 8
Set-Content -Path (Join-Path $stage "manifest.json") -Value $manifestJson -Encoding utf8

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
