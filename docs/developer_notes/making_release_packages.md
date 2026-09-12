# Making release packages

How the downloadable packages that [Installation](../what_is_lsdtt3/install_lsdtt3.md)
tells users to install are actually produced.

> **The operational cheat sheet is `tools/README.md` in the repository.**
> It has the exact commands, the two release routes, the `gh` invocations, and
> how to undo a mistake. This page is the map: what exists, what runs where,
> and why. Do not duplicate command sequences here — they drift.

---

## The shape of it

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
against `project(lsdtt3 VERSION ...)` → configure → build → smoke-check every
CLI with `-v` → package → optionally publish. CI is a thin wrapper that
apt-installs dependencies and calls it, so a local release runs the identical
code path.

**The version lives in exactly one place: line 2 of `CMakeLists.txt`.** Tags,
bundle filenames, the manifest, and the version the binaries report are all
derived from `project(lsdtt3 VERSION ...)`, and `make_release.sh` aborts if the
version it is asked to release disagrees with it. Bumping is therefore a
one-line edit plus the version numbers quoted in the install documentation
(`README.md`, `docs/what_is_lsdtt3/install_lsdtt3.md`), which are examples and have no
authority — always check `CMakeLists.txt` rather than trusting a version you
read in prose.

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
  a minimal APT repository. See also the `.deb` design section in
  `docs/plan_docs/release_ci_plan.md` for the vendored-vs-system-libraries
  trade-off, which is the decision that shapes the whole package.
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
- [`beta_testing_notes.md`](beta_testing_notes.md) — getting bundles to beta testers, conda-forge plan
- `docs/plan_docs/release_ci_plan.md` — design rationale and open questions
