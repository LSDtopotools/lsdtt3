# QGIS plugin development

How the lsdtt3 QGIS plugin is put together, how its versions fit together, and
the loop for changing it.

> **Where the commands live.** The exact packaging and release commands are in
> [`tools/README.md`](https://github.com/LSDtopotools/lsdtt3/blob/main/tools/README.md)
> in the repository, and the end-user install and troubleshooting guide is the
> [QGIS plugin](../qgis.md) page (generated from `qgis_plugin/README.md`). This
> page explains how the pieces relate; it does not repeat those command
> sequences, because copies drift.

---

## Two artefacts, built separately

The plugin does no analysis itself. Its Processing algorithms write a parameter
file and run the `lsdtt-*` command-line tools as subprocesses. A working
install therefore needs **two** things:

| Artefact | Source | Built by | Output |
|---|---|---|---|
| **Plugin zip** — pure Python | `qgis_plugin/lsdtt_qgis/` | `python3 qgis_plugin/package_plugin.py` (any platform) | `qgis_plugin/dist/lsdtt3-<plugin version>.zip` |
| **Backend bundle** — the compiled CLIs plus their libraries and GDAL/PROJ data | `cli/`, `src/` | a per-platform packaging script, **on that platform** | `dist/lsdtt3-backend-<platform>-<profile>-v<backend version>.{zip,tar.gz}` |

In QGIS, the plugin is installed from the zip, and pointed at the unpacked
backend in **lsdtt3 → LSDTT Backend Settings**. Select the bundle root (the
folder containing `bin/`, `share/` and `manifest.json`), not `bin/`.

### Which CLIs the backend needs

The `qgis` profile in `tools/cli_profiles.txt` is the list of CLIs the plugin
wraps. It is the single source of truth for every packaging script. Each CLI
has a wrapper in `qgis_plugin/lsdtt_qgis/algorithms/`.

### Backend per platform

| QGIS on | Build | Package |
|---|---|---|
| Windows | [Native Windows](building_on_windows.md#native-windows-qgis-backend-only) in `build\windows-qgis` | `tools\package_windows_backend.ps1` (default profile `qgis`) |
| Linux | [Linux](building_from_source.md#linux) | `tools/package_linux_backend.sh --cli-profile qgis` |
| macOS | [Native](building_on_macos.md#native-build) in `build/macos-core` | `tools/package_macos_backend.sh --cli-profile qgis` (untested on real hardware) |

A WSL2 build cannot serve as the backend for a Windows QGIS: its binaries are
Linux executables.

---

## Versions: three numbers, two checks

| Number | Where it is set | What reads it |
|---|---|---|
| **Plugin version** | `qgis_plugin/lsdtt_qgis/metadata.txt` (`[general] version=`) | the zip filename; the plugin's range check |
| **Backend version** | the `-Version` / `--version` you pass to the packaging script | the bundle filename; `manifest.json` → `backend_version` |
| **CLI version** | `k_version` in each `cli/lsdtt-*.cpp`, printed in the banner | people — below 1.0.0 means the CLI is still experimental. The plugin does not check it |

When QGIS loads a backend, `qgis_plugin/lsdtt_qgis/backend/manifest.py` rejects
it unless **both** of these hold:

1. **`backend_version` ≥ `MINIMUM_BACKEND_VERSION`** (a constant in
   `manifest.py`). The minimum is raised when the wrappers start depending on
   a CLI behaviour older backends lack.
2. **The plugin version satisfies the bundle's `plugin_version_range`.** The
   packaging scripts write this into `manifest.json`; override the default with
   `-PluginVersionRange` (PowerShell) or `--plugin-version-range` (bash).

> [!WARNING]
> **`backend_version` must be plain `X.Y.Z`.** The plugin's version parser
> accepts digits and dots only. A suffixed version such as `0.5.0-rc1` or
> `0.5.0-local-dirty` packages without complaint and is then refused by the
> plugin as an invalid `backend_version`. For a local test bundle, use a plain
> version number.

Rules of thumb:

- Change the plugin code or UI → bump the plugin version.
- Change what the bundle contains → bump the backend version.
- Make an old backend unusable (a wrapper now emits a parameter older CLIs do
  not understand) → raise `MINIMUM_BACKEND_VERSION`, and update the pinned
  value in `qgis_plugin/lsdtt_qgis/tests/test_backend_health.py` and
  `test_package_content_contract.py`.
- Keep `plugin_version_range` wide enough to cover the plugin versions you
  intend to support.
- Check before packaging: `CMakeLists.txt` line 2 has the release version, and
  a backend packaged at that version must still clear `MINIMUM_BACKEND_VERSION`.

---

## Development loop

### Changing the plugin (Python only)

1. Edit `qgis_plugin/lsdtt_qgis/`.
2. Run the plugin unit tests. They do not need QGIS:
   ```bash
   cd qgis_plugin
   pixi run python -m pytest lsdtt_qgis/tests
   ```
3. Bump `metadata.txt` if this is going to anyone else.
4. `python3 qgis_plugin/package_plugin.py`, then in QGIS **Plugins → Manage and
   Install Plugins… → Install from ZIP**. Installing over the old version
   replaces it; restart QGIS if the old behaviour persists.

### Changing a CLI the plugin wraps

1. Change and test the CLI on Linux as usual.
2. If a parameter was added or renamed, update the wrapper in
   `qgis_plugin/lsdtt_qgis/algorithms/` and its tests.
3. Rebuild the backend on the platform QGIS runs on, and package it (see the
   table above, and `tools/README.md` for the commands).
4. Unzip to a **new**, versioned folder, and point **LSDTT Backend Settings**
   at it. Keeping old bundles makes it easy to switch back.

### Adding a CLI to the plugin

1. Write the wrapper in `qgis_plugin/lsdtt_qgis/algorithms/` and register it
   in `qgis_plugin/lsdtt_qgis/provider.py`.
2. Add the CLI to the `qgis` profile in `tools/cli_profiles.txt`, so every
   packaging script ships it.
3. Add its target to the native Windows build command in
   [Building on Windows](building_on_windows.md#3-configure-and-build).

---

## See also

- [QGIS plugin](../qgis.md) — end-user install, update and troubleshooting
- [`tools/README.md`](https://github.com/LSDtopotools/lsdtt3/blob/main/tools/README.md) — packaging and release commands
- [Making release packages](making_release_packages.md) — what the bundles are and how releases are published
- [Building on Windows](building_on_windows.md) — the native Windows backend build
