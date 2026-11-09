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
| **Plugin zip** — pure Python | `qgis_plugin/lsdtt_qgis/` | `python3 qgis_plugin/package_plugin.py` (any platform) | `qgis_plugin/dist/lsdtt3-qgis-plugin-<plugin version>.zip` |
| **Backend bundle** — the compiled CLIs plus their libraries and GDAL/PROJ data | `cli/`, `src/` | a per-platform packaging script, **on that platform** | `dist/lsdtt3-backend-<platform>-<profile>-v<backend version>.{zip,tar.gz}` |

In QGIS, the plugin is installed from the zip, and pointed at the unpacked
backend in **lsdtt3 → LSDTT Backend Settings**. Select the bundle root (the
folder containing `bin/`, `share/` and `manifest.json`), not `bin/`.

### Which CLIs the backend needs

The `qgis` profile in `tools/cli_profiles.txt` is the list of CLIs the plugin
wraps. It is the single source of truth for every packaging script. Each CLI
has a wrapper in `qgis_plugin/lsdtt_qgis/algorithms/`.

The plugin keeps the same list as `REQUIRED_CLIS` in
`qgis_plugin/lsdtt_qgis/backend/health.py`, and checks that the configured
backend has them. Saving a backend in **LSDTT Backend Settings** checks all
of them and warns, naming every missing CLI; the path is still saved. Running
an algorithm checks only that algorithm's CLI, so an older bundle still runs
the algorithms it has, and the others stop with an error naming the missing
CLI before anything runs. `test_backend_health.py` fails if `REQUIRED_CLIS` differs from the
`qgis` profile or from the CLIs named in `algorithms/*.py`.

### Backend per platform

| QGIS on | Build | Package |
|---|---|---|
| Windows | [Native Windows](building_on_windows.md#native-windows-qgis-backend-only) in `build\windows-qgis` | `tools\package_windows_backend.ps1` (default profile `qgis`) |
| Linux | [Linux](building_from_source.md#linux) | `tools/package_linux_backend.sh --cli-profile qgis` |
| macOS | [Native](building_on_macos.md#native-build) in `build/macos-qgis` | `tools/package_macos_backend.sh` (default profile `qgis`; built by `release.yml` on every tag) |

A WSL2 build cannot serve as the backend for a Windows QGIS: its binaries are
Linux executables.

---

## Versions: one number, one rule

lsdtt3 has one version number. The plugin, the backend bundles, the C++ package
and the Python package all carry the `X.Y.Z` in `project(lsdtt3 VERSION ...)`
on line 2 of `CMakeLists.txt`.

| Copy | Where it is set | What reads it |
|---|---|---|
| **Plugin version** | `qgis_plugin/lsdtt_qgis/metadata.txt` (`[general] version=`), a checked copy | the zip filename; the compatibility check |
| **Backend version** | the `-Version` / `--version` you pass to the packaging script, which must equal `CMakeLists.txt` | the bundle filename; `manifest.json` → `backend_version` |

`tools/set_version.py X.Y.Z` changes every copy together, and the CTest test
`unit_version_consistency` fails if they differ. CLIs have no version of their
own: each carries a stability label, `k_status` (`experimental` or `stable`),
which the banner, `-v` and `--help-json` report. The plugin does not read it.

When QGIS loads a backend, `qgis_plugin/lsdtt_qgis/backend/manifest.py` rejects
it unless **both** of these hold:

1. **The backend has the plugin's major.minor version.** A plugin `X.Y.*`
   accepts a backend `X.Y.*`; the patch level is free, so a plugin-only or
   backend-only patch release still works with the other. A backend from an
   older series gets "install the lsdtt3 X.Y backend bundle"; a newer one gets
   "upgrade the plugin".
2. **The plugin version satisfies the bundle's `plugin_version_range`.** The
   packaging scripts write this into `manifest.json`. The default is the
   bundle's own series, `>=X.Y.0 <X.(Y+1).0`; override it with
   `-PluginVersionRange` (PowerShell) or `--plugin-version-range` (bash).

> [!WARNING]
> **`backend_version` must be plain `X.Y.Z`.** The plugin's version parser
> accepts digits and dots only. A suffixed version such as `0.5.0-rc1` or
> `0.5.0-local-dirty` packages without complaint and is then refused by the
> plugin as an invalid `backend_version`. For a local test bundle, use a plain
> version number.

Rules of thumb:

- A wrapper change that needs a new CLI behaviour (a renamed parameter, a new
  key) is a **minor** version bump. That is what keeps an older backend from
  being accepted: the bump moves the plugin into a series the old backend is
  not in.
- A fix that needs no change on the other side is a **patch** bump.
- Plugins up to 0.4.0 predate this rule. Faced with a 0.7 backend they report
  that it "supports plugin versions '>=0.7.0 <0.8.0'"; the fix is to upgrade
  the plugin.

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
   in `qgis_plugin/lsdtt_qgis/provider.py`. A CLI that computes several
   independent products can be several algorithms, one short form each:
   `drainage_analysis.py` does this for `lsdtt-drainage-analysis` with a
   shared base class, and writes every `write_*` switch so a run produces
   only what its form asked for.
2. List the algorithm in the table in `qgis_plugin/README.md` (the user page
   on the documentation site).
3. Add the CLI to the `qgis` profile in `tools/cli_profiles.txt`, so every
   packaging script ships it, and to `REQUIRED_CLIS` in
   `qgis_plugin/lsdtt_qgis/backend/health.py`. The plugin tests fail until
   all three lists agree.
4. Add its target to the native Windows build command in
   [Building on Windows](building_on_windows.md#3-configure-and-build).

---

## See also

- [QGIS plugin](../qgis.md) — end-user install, update and troubleshooting
- [`tools/README.md`](https://github.com/LSDtopotools/lsdtt3/blob/main/tools/README.md) — packaging and release commands
- [Making release packages](making_release_packages.md) — what the bundles are and how releases are published
- [Building on Windows](building_on_windows.md) — the native Windows backend build
