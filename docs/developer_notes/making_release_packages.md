# Making release packages

How the downloadable packages that [Installation](../what_is_lsdtt3/install_lsdtt3.md)
tells users to install are actually produced.

> **The operational cheat sheet is `tools/README.md` in the repository.**
> It has the exact commands, the two release routes, the `gh` invocations, and
> how to undo a mistake. This page give details: what exists, what runs where,
> and why. If you are reading this in code and thinking of making edits: don't change the command sequences or record them here, that stays in `tools/README.md`

---

## The basics

One bundle per platform, each built **on that platform**. There is no
cross-compilation path: every packaging script interrogates the native loader
(`ldd`/`objdump` on Linux, `dumpbin` on Windows, `otool` on macOS) to work out
what to bundle.

| Platform | Profile | Built by | Script | Artefact |
|---|---|---|---|---|
| Linux x86_64 | `core` | `release.yml`, first job, on a version tag | `tools/package_linux_backend.sh` | `lsdtt3-backend-linux-x86_64-core-v<ver>.tar.gz` |
| Windows x64 | `qgis` | `release.yml`, after Linux | `tools/package_windows_backend.ps1` | `lsdtt3-backend-win64-qgis-v<ver>.zip` |
| macOS arm64, x86_64 | `qgis` | `release.yml`, after Linux | `tools/package_macos_backend.sh` | `lsdtt3-backend-macos-<arch>-qgis-v<ver>.tar.gz` |
| QGIS plugin | — | `release.yml`, Linux job | `qgis_plugin/package_plugin.py` | `lsdtt3-qgis-plugin-<ver>.zip` |
| Python wheels | — | `python-wheels.yml`, called by `release.yml` after Linux | cibuildwheel | `lsdtt3-<ver>-cp3XX-*.whl` (20) |

**Windows and macOS ship the `qgis` profile only.** Those bundles exist to
back the QGIS plugin. Anyone who wants the full CLI set on Windows uses the
Linux bundle under WSL2, and on macOS the Linux bundle in a VM; both are the
reference build, and faster.

The same packaging scripts run locally, so any artefact can still be rebuilt
by hand on its platform (`tools/README.md`).

`tools/make_release.sh` is the single driver: resolve version → check it
against `project(lsdtt3 VERSION ...)` → configure → **check that no optional
dependency was silently dropped** → build → smoke-check every CLI with `-v` →
package → optionally publish. CI is a thin wrapper that apt-installs
dependencies and calls it, so a local release runs the identical code path.

The release build configures with `-DBUILD_TESTING=OFF`. A release packages
nothing from the CTest suite, so compiling it is billed CI minutes spent on
artefacts that are thrown away. Run the tests before releasing, not during.

**lsdtt3 has one version number, and it lives on line 2 of `CMakeLists.txt`.**
The C++ package, the backend bundles, the Python package and the QGIS plugin
all carry it, and every release is the tag `vX.Y.Z`. Tags, bundle filenames,
the manifest, and the version the binaries report are derived from
`project(lsdtt3 VERSION ...)`. Five files hold checked copies, because their
ecosystems read them as they stand: `python_bindings/pyproject.toml`,
`qgis_plugin/lsdtt_qgis/metadata.txt`, `conda-recipe/meta.yaml`, `pixi.toml`
and the top entry of `debian/changelog`. `tools/set_version.py` rewrites all
six together; the CTest test `unit_version_consistency` fails if any differs,
and `make_release.sh` and `make_all_packages.sh` refuse to build. Only the
ecosystem rebuild suffixes may differ: `.postN` in `pyproject.toml` (a
rebuilt wheel) and `+pkgN` in the Debian changelog.

For a new X.Y.Z, `tools/set_version.py` then regenerates the test inventory
(`docs/developer_notes/test_inventory.md`), which names the version it was
generated from. It configures `build/inventory_tests` with every test tier
registered (nothing is compiled; about a minute the first time, then
seconds), which needs what the full test suite needs: Python 3 with
`osgeo.gdal`, numpy, rasterio, geopandas, shapely and pytest, and `gdal-bin`.
`--inventory-build-dir DIR` reuses a build directory you have already
configured that way; `--skip-inventory` skips the step, and
`unit_test_inventory_is_current` then fails until the inventory is
regenerated. Commit the regenerated page with the version change, before
tagging: the documentation site is built from the tag. The version numbers quoted
in the install documentation (`README.md`,
[Installing lsdtt3](../what_is_lsdtt3/install_lsdtt3.md)) are examples and have no authority.

CLIs have no version of their own. Each carries a stability label, `k_status`
(`experimental` or `stable`), shown in its banner, its `-v` output and its
`--help-json`, and on its reference page. Promoting a CLI to `stable` is a
review decision.

**Choosing the number.** Bump the **minor** version when any component changes
in a way another component or a user's parameter files depend on: a parameter
rename, a removed output, a CLI behaviour the plugin wrappers rely on, a Python
API break. Bump the **patch** version for fixes that need no change elsewhere.
Components with the same major.minor work together; the QGIS plugin enforces
that against the backend.

**Every tag gets every artefact.** A tag is one consistent snapshot, so each
artefact built from it contains everything on `main` at that commit. Releases
up to 0.8.0 got wheels only by hand; from then on `release.yml` builds them on
every tag. `make_release.sh` also builds the bindings once on Linux and runs
their whole test suite in the first job, so a C++ change that breaks them stops
the release before any other platform starts. Releases are tagged on `main`; a `release/X.Y`
branch is cut only when a fix must ship while `main` holds work that should
not. The commands, the worked example and the rebuild-without-a-version-change
table are in `tools/README.md` ("Release cheat sheet").

Which CLIs go in a bundle is a **profile** (`tools/cli_profiles.txt`, the single
source of truth for all four scripts): `qgis` (the plugin subset; Windows and
macOS), `core` (Linux), `all` (local releases only). The profile name appears
in the filename.

---

## How the release workflow is ordered

`release.yml` runs only on `v*` tag pushes (and manual dispatch), never on
ordinary pushes or pull requests. Its first job builds the Linux bundle and the
plugin zip, runs the bindings tests, checks the CLI reference and the docs
site, and creates the GitHub release. **Nothing else starts unless it
succeeds**, so a broken tag fails once, on the cheapest runner. The Windows and
macOS backends then build in parallel and attach themselves to that release;
the wheels follow the Windows backend, because the two share one pinned vcpkg
binary cache and a cold cache would otherwise build GDAL from source twice.

The **glibc of the build host is the Linux bundle's compatibility floor.**
Building on `ubuntu-22.04` gives a deliberately old, reproducible floor that a
personal machine's rolling distro would not. The macOS bundles have the same
property for macOS: Homebrew bottles target the runner's release, so a bundle
built on `macos-15` needs macOS 15 or later.

**macOS is the least exercised path.** No maintainer builds on a Mac, so the
first real run of `package_macos_backend.sh` is in CI. It was reviewed and
partly tested off-Mac for that reason: the whole tree passes a clang + libc++
syntax check (including with macOS's `int64_t` = `long long`, and against the
GDAL headers of the version Homebrew ships, which is newer than Ubuntu's or
vcpkg's: GDAL 3.13 changed `GDALGetMetadata` to return `CSLConstList`, and
that broke the first macOS release build), the dependency
resolver and the dyld load-log parser are tested against a fake `otool`, and
the script fails hard, naming the library, if the bundle is not
self-contained. A missing optional dependency (OpenMP, FFTW, Eigen, netCDF)
is reported as a warning annotation on the run, not a failure.

**The Windows vcpkg cache.** A cold cache builds GDAL from source: about 95
minutes. A GitHub Actions cache is visible only to the branch or tag that saved
it, plus the default branch, so a cache saved by one release tag cannot be used
by the next: every release would start cold. `vcpkg-cache.yml` therefore keeps
the cache on `main`, where every tag can restore it. It runs when the package
list changes, every five days (GitHub deletes caches unused for 7 days), and
by hand. All three workflows pin the same `VCPKG_TAG` and package list
(`.github/vcpkg-packages.txt`); change them together, and run *Windows vcpkg
cache* by hand on `main` before tagging a release that uses a new `VCPKG_TAG`.

---

## Publishing

**GitHub Releases is the one place users are sent.** The first job of the tag
build creates the release and every other CI artefact attaches itself to it.
The Debian and conda packages are built locally and attached with
`gh release upload v<version> <file>`.

Each backend bundle and the plugin zip is attached twice: under its versioned
name, and under an unversioned one (`lsdtt3-backend-linux-x86_64-core.tar.gz`,
`lsdtt3-qgis-plugin.zip`, ...). The unversioned names make permanent links,
`https://github.com/LSDtopotools/lsdtt3/releases/latest/download/<name>`,
which always serve the newest release; the install documentation uses them, so
it never needs a version number. Wheels cannot be renamed (pip reads the
version and Python tag from the file name), so the documentation looks up the
latest release's wheels through the GitHub API instead. Each bundle's own
README links its exact versioned release.

The release is created as a **draft** and published only by the last job,
`publish-release`, once every build has finished (whether or not the Windows
and macOS jobs succeeded, as long as Linux did). Until then the `latest` links
and the API's latest release stay on the previous version, so notebooks and
the install documentation never pick up a half-filled release. The same job
deletes the unversioned copies from every older release: `latest` never
resolves to them, so they would only be about 260 MB of dead weight per
release. Their versioned files are kept. The names it deletes are listed
explicitly in the job (`LATEST_COPIES`); add any new unversioned artefact
there.

A copy of some packages on the personal web space
(<https://www.geos.ed.ac.uk/~smudd/lsdtt_packages/>) is mentioned in the
install documentation as a last resort only, for networks that cannot reach
GitHub, with a warning that it is usually out of date. It is not part of the
release procedure and nothing links to it as the primary route.

---

## Every script ends with a self-containment check

Each packaging script verifies that every dependency of every staged binary
either resolves inside the bundle or is OS-provided, and **fails the packaging
run** otherwise. This is not decoration: before it existed, the scripts silently
skipped libraries they could not resolve, which is exactly how v0.2.0 shipped a
bundle that died with `libgdal.so.30: cannot open shared object file` on users'
machines.

Related invariant: **no symlinks in bundles.** Each library is staged as a
regular file named by its SONAME (Linux) or install-name basename (macOS), so
the bundle survives any filesystem or archiver.

---

## The bundle README is the only documentation shipped

Each packaging script writes the bundle's `README.md` from a template held in
the script itself (a heredoc in the `.sh` scripts, a here-string in the `.ps1`).
A bundle contains no other documentation, and its reader may be an LLM with
nothing else to go on. So besides install steps, the README has a section
"How the tools work: parameter files": options are set in a parameter file, not
with flags; `-h` lists the parameters; the file format; the file I/O keys; and
a worked hillshade example. The binaries print the same guide themselves, when
run with no arguments and at the top of `-h` (see `lsdtt_param_parser` in the
architecture document). If the parameter-file rules change, update both, in
all three scripts.

---

## Bundled runtime data

A bundle stages more than shared libraries. GDAL and PROJ both read **data
files** at run time, and lsdtt3 ships one of its own, and a bundle that has the
libraries but not the data is broken in a way the self-containment check above
cannot see — `ldd` is satisfied and the tools still fail, or worse, silently use
someone else's data.

| Staged path | Comes from | Why |
|---|---|---|
| `share/proj/` | `proj.db` and the transformation grids | every CRS lookup and datum shift |
| `share/gdal/` | `gdal-config --datadir` | driver support files (`gdalvrt.xsd` and friends) |
| `share/lsdtt3/` | `data/crn_parameters.nc` in the repository | the cosmo tools' production-rate and palaeomagnetic parameters |

`share/lsdtt3/` is the same contract as the other two — the cosmo tools resolve
it as `<exe>/../share/lsdtt3` — and it is just as invisible when it is missing:
the tools start, report their version, and then cannot run without a download.

### The binaries find these themselves

Since the `lsdtt_gdal_init()` change, every tool pins these directories at
startup from its own location — `<exe>/../share/proj` and
`<exe>/../share/gdal`, which is exactly the layout the packaging scripts stage.
**A bundled `proj.db` takes precedence over `PROJ_DATA`, `PROJ_LIB` and
`GDAL_DATA` in the user's environment.**

Two consequences for packaging:

1. **`setenv.sh` is no longer load-bearing for data.** It still sets `PATH` and
   `LD_LIBRARY_PATH`, which the binaries cannot pin for themselves, but its
   `GDAL_DATA` / `PROJ_DATA` exports are now belt-and-braces. The practical gain
   is that the tools work when it has *not* been sourced — the Google Colab case,
   where each `!` line is a separate shell and the exports never persist.
2. **The staged layout is now a contract, not a convention.** If a future
   packaging change moves `share/proj` relative to `bin/`, the binaries will stop
   finding it and fall back to the user's environment — silently. Keep
   `bin/../share/proj`, and the same for `bin/../share/lsdtt3`.

### Verifying it worked

The self-containment check does not cover data files, so check by hand on the
staged bundle:

```bash
./bin/lsdtt-raster-preprocessing -v
```

```
GDAL 3.8.4, PROJ 9.3.1
  proj data: bundled: /path/to/bundle/share/proj
  GDAL data: bundled: /path/to/bundle/share/gdal
```

**Both lines must say `bundled:`.** If either says `PROJ_DATA=…`,
`PROJ built-in search` or `GDAL built-in search`, the data did not stage where
the binaries look, and the bundle will use whatever the user happens to have.
Run this before publishing, with `PROJ_DATA` and `GDAL_DATA` deliberately set to
junk in your shell — that is the condition a user's machine simulates.

### Why this matters more than it sounds

A mismatched `proj.db` does not usually error. It resolves a datum differently
and returns coordinates that are plausible and wrong — right shape, right part
of the world, tens of metres out. That is why `-v` reports the data directories
at all, and why it is the first thing to ask a user for in a projection bug
report.

---

## After adding a dependency

When a new library becomes a build dependency in `CMakeLists.txt`:

1. Add the `-dev` package to the Linux `Install dependencies` step of
   `.github/workflows/release.yml`.
2. Add the Homebrew formula to the macOS `Install dependencies` step of
   `release.yml`, to `before-all` under `[tool.cibuildwheel.macos]` in
   `python_bindings/pyproject.toml`, and to the `brew install` line in
   `building_on_macos.md`.
3. Add the vcpkg port to `.github/vcpkg-packages.txt` (Windows backend and
   Windows wheel; editing it changes the binary-cache key).
4. Add it to `apptainer_def/lsdtt3.def`.

There is **no** separate runtime package list to maintain — the packaging
scripts resolve the runtime closure themselves and fail if anything is missing.

---

## Alternative distribution routes

- [`debian_packaging.md`](debian_packaging.md) — building a `.deb`, and hosting
  a minimal APT repository. The decision that shapes the whole package is
  vendored-versus-system libraries: the `.deb` links the system GDAL, PROJ,
  FFTW and netCDF through `Depends:` so the package stays small and tracks
  distribution security updates, whereas the tarball backends vendor their
  closure so they run on a machine with no matching system libraries at all.
- **conda-forge** — the intended route for the public repository;
  `conda-recipe/meta.yaml` is the draft recipe. It is also the right answer for
  QGIS integration, because QGIS and lsdtt3 would then share one conda-managed
  GDAL. Details in [`beta_testing_notes.md`](beta_testing_notes.md).
- **Apptainer** — `apptainer_def/lsdtt3.def` builds a container with all
  dependencies on Ubuntu 22.04; see `apptainer_def/apptainer_setup.md`.

---

## See also

- `tools/README.md` (in the repository) — the commands, in order
- [`building_from_source.md`](building_from_source.md) — the builds these package
- [`qgis_plugin_development.md`](qgis_plugin_development.md) — the QGIS plugin, and the version checks a backend bundle must pass
- [`beta_testing_notes.md`](beta_testing_notes.md) — getting bundles to beta testers, conda-forge plan

