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

| Platform | Built by | Script | Artefact |
|---|---|---|---|
| Linux x86_64 | GitHub Actions on a version tag | `tools/package_linux_backend.sh` | `lsdtt3-backend-linux-x86_64-<profile>-v<ver>.tar.gz` |
| Windows x64 | by hand, VS developer shell | `tools/package_windows_backend.ps1` | `lsdtt3-backend-win64-<profile>-v<ver>.zip` |
| macOS | by hand, on a Mac — **not yet in CI** | `tools/package_macos_backend.sh` | `lsdtt3-backend-macos-<arch>-<profile>-v<ver>.tar.gz` |

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
rebuilt wheel) and `+pkgN` in the Debian changelog. The version numbers quoted
in the install documentation (`README.md`,
`docs/what_is_lsdtt3/install_lsdtt3.md`) are examples and have no authority.

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

**Every tag gets every artefact, except the Python wheel.** A tag is one
consistent snapshot, so each artefact built from it contains everything on
`main` at that commit. The wheel is built by hand, and only for releases where
the Python API is supported; between wheels, `make_release.sh` builds the
bindings once on Linux and runs their whole test suite, so a C++ change that
breaks them stops the release. Releases are tagged on `main`; a `release/X.Y`
branch is cut only when a fix must ship while `main` holds work that should
not. The commands, the worked example and the rebuild-without-a-version-change
table are in `tools/README.md` ("Release cheat sheet").

Which CLIs go in a bundle is a **profile** (`tools/cli_profiles.txt`, the single
source of truth for all four scripts): `qgis` (the plugin subset, Windows
default), `core` (Linux/macOS/CI default), `all` (local releases only). The
profile name appears in the filename.

---

## Why Linux is the CI build and the others are not

The **glibc of the build host is the bundle's compatibility floor.** Building
Linux bundles on `ubuntu-22.04` in Actions gives a deliberately old,
reproducible floor that a personal machine's rolling distro would not — that is
the reason for CI, not convenience.

Windows in CI would mean building GDAL through vcpkg on every run: hours on a
cold cache, and Windows runners bill at 2× minutes. macOS runners bill at 10×.
Both platforms are available locally, so both are built by hand and attached to
the release the Linux build already created.

The repository is private, so Actions minutes are metered. CI therefore fires
**only** on `v*` tag pushes and manual dispatch — never on ordinary pushes or
pull requests.

---

## Publishing

Two destinations, same files:

1. **GitHub Releases** — created automatically by the tag build; the Windows
   (and later macOS) bundles are attached afterwards with
   `gh release upload v<version> dist/<bundle>`.
2. **<https://www.geos.ed.ac.uk/~smudd/lsdtt_packages/>** — a plain directory
   on the personal web space, copied up by hand. This is the route that works
   for people who are not repository collaborators, and it is the URL baked
   into the bundle README's Colab instructions. **Keep it in step with the
   GitHub release**: a user following the Colab instructions downloads from
   here, so a missing or stale file there is a broken install path even when
   GitHub is fine.

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

1. Add the `-dev` package to the `Install dependencies` step of
   `.github/workflows/release.yml`.
2. Add it to `apptainer_def/lsdtt3.def`.

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
- **conda-forge** — the intended route once the repository is public;
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

