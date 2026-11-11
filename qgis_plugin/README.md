# lsdtt3 QGIS plugin

The plugin adds an **lsdtt3** group to the QGIS Processing Toolbox. Each
algorithm runs one of the lsdtt3 command-line tools, so you install two
things: the **plugin** and a **backend** (the command-line tools for your
system).

## Installation

**1. Download the plugin and the backend for your system.** These links always
fetch the latest release:

| What | Download |
|---|---|
| Plugin (all systems) | [`lsdtt3-qgis-plugin.zip`](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-qgis-plugin.zip) |
| Backend, Windows | [`lsdtt3-backend-win64-qgis.zip`](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-win64-qgis.zip) |
| Backend, macOS Apple Silicon (M1 and later) | [`lsdtt3-backend-macos-arm64-qgis.tar.gz`](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-macos-arm64-qgis.tar.gz) |
| Backend, macOS Intel | [`lsdtt3-backend-macos-x86_64-qgis.tar.gz`](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-macos-x86_64-qgis.tar.gz) |
| Backend, Linux | [`lsdtt3-backend-linux-x86_64-core.tar.gz`](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-linux-x86_64-core.tar.gz) |

Download both from GitHub, not from another copy you come across: the plugin
and the backend must come from the same `X.Y` series (plugin 0.9.x works with
any 0.9.x backend; see [Versioning](#versioning-plugin-and-backend)). Older
releases are on the [list of releases](https://github.com/LSDtopotools/lsdtt3/releases).

**2. Unpack the backend** into a folder of its own, for example in your
Documents folder. Do not unpack the plugin zip.

- Windows: right-click the zip, **Extract All…**
- macOS: double-click the `.tar.gz`, then, in Terminal, clear the quarantine
  flag the browser set, or macOS blocks the tools:
  `xattr -dr com.apple.quarantine <unpacked folder>`. The backend's own
  `README.md` has the details.
- Linux: `tar -xzf lsdtt3-backend-linux-x86_64-core.tar.gz`

**3. Install the plugin in QGIS:** **Plugins → Manage and Install Plugins… →
Install from ZIP**, choose `lsdtt3-qgis-plugin.zip`, and click **Install
Plugin**.

**4. Point the plugin at the backend:** open **lsdtt3 → LSDTT Backend
Settings** and select the unpacked backend's root folder, the one containing
`bin/`, `share/` and `manifest.json` (not `bin/` itself). Saving checks the
folder for every tool the plugin uses.

**5. Check it works:** run **Raster Preprocessing** from the **lsdtt3** group
of the Processing Toolbox on a small DEM, with hillshade switched on.

To update later, repeat steps 1–4: **Install from ZIP** replaces the plugin in
place, and the new backend goes in a new folder, which you then select in the
settings. More on that, and on fixing errors, is in
[Updating to a new version](#updating-to-a-new-version) and
[Backend command-line tools](#backend-command-line-tools-required) below.

Developers: how the plugin, its backend bundle and their version checks fit
together, and the loop for changing them, is in
[QGIS plugin development](https://lsdtopotools.github.io/lsdtt3/developer_notes/qgis_plugin_development.html)
on the documentation site. Packaging commands are in
`tools/README.md`.

## Algorithms

The plugin adds an **lsdtt3** group to the Processing Toolbox. Every algorithm
writes a parameter file next to the input DEM and runs one backend CLI, so any
run can be repeated from the command line (the parameter file's header gives
the command).

| Algorithm | Backend CLI | What it does |
|---|---|---|
| Raster Preprocessing | `lsdtt-raster-preprocessing` | Reproject, mask, fill nodata, edge-preserving smoothing, pixel replacement, breach or fill depressions (with flat resolution and sink candidates), vertical-unit checks, hillshade |
| Clip Raster to Polygon | `lsdtt-raster-preprocessing` | Crop a raster to a polygon, optionally masking outside it |
| Clip Vector to Raster Extent | `lsdtt-raster-preprocessing` | Clip a vector layer to a raster's bounding box (any CRS) |
| Mosaic Rasters | `lsdtt-raster-preprocessing` | Merge pixel-aligned rasters on the same grid into one |
| Sample Rasters at Points | `lsdtt-raster-preprocessing` | Sample a list of rasters at each point of a point layer |
| Channel Extraction | `lsdtt-channel-extraction` | Channel sources (threshold or Wiener), network lines, points and junctions |
| Chi Analysis | `lsdtt-chi-analysis` | Chi extraction, chi-network fitting with an m/n sweep, k_sn maps |
| Basin Extraction | `lsdtt-basin-extraction` | Basins by junction, outlet, nesting, order or crossing line, with filters |
| Side Tributary Basins | `lsdtt-basin-extraction` | Basins draining into the side of one main channel |
| Basin Statistics | `lsdtt-basin-stats` | Per-basin relief, hypsometry, shape and channel statistics |
| Surface Metrics | `lsdtt-surface-metrics` | Hillshade, relief, polyfit slope, aspect and curvature, roughness; Frangi, Otsu, Hough, curvature skeletons and adaptive VEC thresholding |
| Channel Network Analysis | `lsdtt-channel-network-analysis` | Link bearings, junction angles, Strahler and Tokunaga layers, area ratios |
| Swath Profile | `lsdtt-swath-profile` | Percentiles of a raster binned along a baseline line |
| Drainage Area | `lsdtt-drainage-analysis` | D8, D-infinity and/or MFD drainage area |
| Topographic Wetness Index | `lsdtt-drainage-analysis` | ln(a / tan b) from D-infinity or MFD routing |
| Sediment Connectivity Index | `lsdtt-drainage-analysis` | Cavalli et al. (2013) index, with its D_up and D_down rasters |
| Height Above Channel | `lsdtt-drainage-analysis` | Height of each pixel above the channel its flow path reaches |
| Soil Topographic Index | `lsdtt-drainage-analysis` | Wetness index less a soil term from soil and land-use data |

Each algorithm's form shows the parameter help from the backend itself, and
the full parameter lists are in the CLI reference pages.

### Which CLIs the backend needs

The plugin uses the nine CLIs in the table. The `qgis` backend bundle ships
exactly these (the `qgis` profile in `tools/cli_profiles.txt`).

- **When you save a backend folder** in `lsdtt3 -> LSDTT Backend Settings`,
  the plugin checks for all nine and warns, naming every one that is missing.
  The folder is still saved.
- **When you run an algorithm**, only that algorithm's CLI must be present.
  An older bundle still runs the algorithms it has; the others stop before
  anything runs, with a message naming the missing CLI.

Bundles made before this release have no `lsdtt-swath-profile` or
`lsdtt-drainage-analysis`. Install a new bundle to use Swath Profile and the
five drainage algorithms.

## Versioning (plugin and backend)

The plugin and the backend share one version number, the lsdtt3 version. It is
set in `project(lsdtt3 VERSION ...)` in the root `CMakeLists.txt`;
`qgis_plugin/lsdtt_qgis/metadata.txt` carries a copy, and `tools/set_version.py`
changes both together.

- **Plugin version** comes from `metadata.txt` (`[general] version=...`) and is
  used in the plugin zip name generated by `qgis_plugin/package_plugin.py`.
- **Backend version** comes from the version passed to the packaging script and
  is written into the backend bundle name and `manifest.json`
  (`backend_version`).

**The rule: a plugin `X.Y.*` accepts a backend `X.Y.*`.** Major and minor must
match; the patch level is free, so a plugin-only or backend-only patch release
still works with the other. A minor version bump marks a change the other side
depends on, such as a renamed parameter the wrappers emit. A backend from
another series is refused before anything runs, with a message saying whether
to install the matching backend or upgrade the plugin.

The bundle's `manifest.json` also carries `plugin_version_range`. The packaging
scripts set it to the bundle's own series (`>=X.Y.0 <X.(Y+1).0`) unless
`--plugin-version-range` / `-PluginVersionRange` overrides it, and the plugin
refuses a bundle whose range excludes it.

Plugin 0.4.0 and older do not know this rule. With a 0.7 backend they report
that the backend "supports plugin versions '>=0.7.0 <0.8.0'": upgrade the
plugin.

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

- `qgis_plugin/dist/lsdtt3-qgis-plugin-<version>.zip`

where `<version>` is read from `qgis_plugin/lsdtt_qgis/metadata.txt`.

## Installing a plugin zip you built

Install a zip from `qgis_plugin/dist/` exactly as the downloaded one
([Installation](#installation), step 3). It must match the backend's `X.Y`
series like any other plugin.

## Basin extraction: select all basins of one stream order

The QGIS `Basin Extraction` tool includes a selection mode:

- `All basins of stream order`

When this mode is chosen, set the companion numeric parameter:

- `Selected stream order`

This maps to CLI parameters `select_all_basins_of_stream_order: true` and
`selected_stream_order: <n>` and extracts every basin whose junction stream
order is exactly `<n>`.

## Swath Profile

Bins every valid pixel within a half-width of a baseline line by its distance
along the line, and writes the pixel count, mean and 21 percentiles (p0, p5,
... p100) of each bin to `<prefix>_swath_profile.csv`. It profiles the DEM's
elevation unless you give another raster on the same grid.

- **Baseline:** the first selected feature of the line layer, or its first
  feature if none is selected. The layer can be in any CRS: the plugin
  reprojects the line to the DEM's CRS and saves it as
  `<prefix>_qgis_baseline.geojson`, which the backend reads (the backend
  itself does not reproject vectors). A feature with more than one part is
  refused; split or merge it first.
- **Half-width 0** includes every valid pixel of the DEM.
- **Distance rasters** (optional): each pixel's distance along the baseline,
  its distance from it, and that distance signed by side.

## Drainage algorithms

Drainage Area, Topographic Wetness Index, Sediment Connectivity Index, Height
Above Channel and Soil Topographic Index all run `lsdtt-drainage-analysis`,
one product each. They share the DEM, pit-filling and output options; the
advanced options include how flats are made drainable and the boundary
conditions. Each run writes only the outputs its own form asks for.

- **Drainage Area** writes any combination of D8, D-infinity and MFD
  drainage area, in square metres.
- **Sediment Connectivity Index** always writes the index and its D_up and
  D_down components. Channels, where downslope paths stop, are pixels whose
  drainage area reaches the threshold in square metres.
- **Height Above Channel** defines channels by contributing pixels instead.
- **Soil Topographic Index** takes soil-series and land-use polygon layers and
  two tables: saturated conductivity profiles (needed for the STI raster and
  the lookup CSV) and standard percentage runoff (needed for the proxy
  raster, the lookup CSV, and spr mode). The tables follow the Scottish soil
  data (Hutton series codes, SSKIB conductivities); see the
  `lsdtt-drainage-analysis` reference for their columns. The plugin says
  which table is missing before anything runs.

## Re-running an analysis

Every algorithm names its outputs `<prefix>_<suffix>` in the output folder, so
running it again with the same folder and prefix writes the same files. The
**Overwrite existing outputs with the same prefix** option (on by default)
makes that work:

- Before the run, the plugin closes any layer in the project that is loaded
  from one of those files. QGIS keeps a loaded layer's file open, and on
  Windows an open file cannot be replaced, so without this a re-run fails or
  keeps the old file. The log names each layer it closed; the new results are
  loaded in their place.
- With the option off, a run whose folder already holds `<prefix>_...` files
  stops before anything runs, so earlier results are never replaced by
  accident. Use a new prefix to keep both.
- An input that is itself one of this run's outputs (for example re-running
  on `<prefix>_fill.tif` with the same prefix) is refused either way: choose
  another prefix.
- After the run, every output is checked to have been written by this run. A
  file that another program still holds open, and so was not replaced, is
  reported rather than loaded as if it were new.

The run log (`<prefix>_<tool>.log`, replaced each run) and the numbered
parameter file next to the input DEM never count as outputs.

## Updating to a new version

You usually do **not** need to manually uninstall first.

- **Plugin update:** download the new plugin zip from the latest release (or
  build one) and install it again via `Install from ZIP`. QGIS replaces the
  existing plugin files.
- **Backend update:** unzip the new backend bundle to a new folder (recommended
  to keep versioned folder names), then point the plugin to that folder in
  `lsdtt3 -> LSDTT Backend Settings`.

Recommended update flow:

1. Close QGIS.
2. Install the new plugin zip (or reinstall over the old one).
3. Unzip the new backend bundle.
4. Open QGIS and update backend path in settings.
5. Run a quick check (`lsdtt-raster-preprocessing --version` from that backend,
   or run a small tool in QGIS).

If QGIS still shows old behavior, restart QGIS once more; if needed, uninstall
the plugin and reinstall from the new zip.

## Backend command-line tools (required)

The plugin calls external `lsdtt-*` executables, so a backend must be
unpacked and selected before any algorithm runs
([Installation](#installation), steps 1, 2 and 4).

### Checking a Windows backend outside QGIS

Developers build the Windows backend with
`.\tools\package_windows_backend.ps1 -Version <version>` (see
`tools/README.md`). To check an unpacked one works on its own, load its
runtime environment in PowerShell and ask a tool for its version:

```powershell
cd C:\path\to\lsdtt3-backend-win64-qgis-v<version>
. .\setenv.ps1
.\bin\lsdtt-raster-preprocessing.exe --version
```

If `--version` works, the backend is ready. On Linux and macOS, run
`bin/lsdtt-raster-preprocessing --version` from the unpacked folder.

### Troubleshooting PROJ / `proj.db` errors

If you see errors like:

- `proj_create_from_database ... proj.db contains DATABASE.LAYOUT.VERSION...`
- `input raster has no CRS defined` immediately after a PROJ error

then QGIS' own PROJ data is being used instead of the backend bundle's PROJ
data. Usually this means the backend path is wrong.

Fix:

1. In `lsdtt3 -> LSDTT Backend Settings`, select the **bundle root**
   (contains `bin/`, `share/`, `manifest.json`), not `bin/` itself.
2. Save the settings again (the plugin then checks the bundle for every CLI
   it uses), and run a small preprocessing run.

The plugin now auto-corrects a `.../bin` selection to its parent bundle root
when possible, but selecting the root folder directly is still recommended.

### Troubleshooting GDAL plugin DLL warnings on Windows

If processing succeeds but logs show messages like:

- `Can't load requested DLL ... gdal_GEOR.dll ... The specified procedure could not be found`
- similar lines for `ogr_FileGDB.dll` or `ogr_OCI.dll`

these usually come from QGIS GDAL plugin-driver paths leaking into the backend
subprocess environment. They are typically noisy warnings (not fatal) when the
output is still produced.

Recent plugin builds isolate backend GDAL by disabling external GDAL plugin
autoload for subprocesses. If you still see this, update the plugin zip to the
latest build and reinstall it.
