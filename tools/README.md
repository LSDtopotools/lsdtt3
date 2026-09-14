# Making lsdtt3 releases

This directory contains utility scripts for building and packaging
redistributable backend bundles. This readme also contains a cheat sheet for
spinning up a new `lsdtt3` release.

## Finding the latest version number

**`CMakeLists.txt` line 2 is the single source of truth.** Every other version
string in the project — release tags, package filenames, the version the
binaries print — is derived from it:

```bash
sed -n 2p CMakeLists.txt      # project(lsdtt3 VERSION 0.3.0 LANGUAGES CXX)
```

Read it before quoting a version anywhere. Documentation that spells out a
version in an example command is a snapshot and drifts; this line does not.

To see what has actually been *published* (which can lag the version in
`CMakeLists.txt`, if it has been bumped ahead of a release):

```bash
# latest tag available in your local checkout
git describe --tags --abbrev=0

# latest published GitHub release (requires gh auth)
gh release list --limit 1
```

For script arguments, use the numeric part without `v` (for example `0.3.0`).
If you are cutting a new release, that value must match
`project(lsdtt3 VERSION ...)` in `CMakeLists.txt` — `make_release.sh` reads
that line and refuses to run on a mismatch.

For local Windows packaging, `package_windows_backend.ps1` does **not** validate
the version string against git tags or `CMakeLists.txt`: it simply writes your
`-Version` into the zip name, bundle README, and `manifest.json`.

## Release cheat sheet

> **The scripts never create, move, or push a git tag.** `make_release.sh`
> only *reads* a tag, and only when you leave `--version` off. Tagging is
> yours to do — with one exception: `--publish` calls `gh release create`,
> and if the tag does not already exist on GitHub, **gh creates it
> server-side** from the tip of the default branch. See route 2.

Releasing `0.3.0`, start to finish:

**1. Bump the version and get it onto main.** 
*You must first of all increment the version in the CMakeLists.txt file.* This version number will need to match your tag.
`make_release.sh` refuses to run
if the version it is releasing does not match `project(lsdtt3 VERSION ...)` in
`CMakeLists.txt`, so this is step one, not an afterthought. The tag must be
exactly `v` + that version (`v0.3.0` for `VERSION 0.3.0`); `CMakeLists.txt` uses
the same rule to decide whether a build reports a plain release version or a
`-dev` one, so a mismatch silently produces `-dev` binaries in what looks like
a release.

```bash
# CMakeLists.txt line 2:  project(lsdtt3 VERSION 0.3.0 LANGUAGES CXX)
git add CMakeLists.txt
git commit -m "Bump version to 0.3.0"
git push
git status                       # release from a clean tree: see step 6
```

**2. Rehearse locally** (optional but cheap — it is the same script CI runs,
and it catches a version mismatch or a CLI that fails `-v` before a tag exists):

```bash
tools/make_release.sh --version 0.3.0        # builds + packages, publishes nothing
```

**3. Tag and push — this is what triggers the release.**

```bash
# annotated (-a with -m), not lightweight: it carries an author, date and
# message, and `git describe` prefers it
git tag -a v0.3.0 -m "lsdtt3 0.3.0"
git push origin v0.3.0           # pushing the branch does NOT push tags
```

**4. Wait for Actions**, which builds and publishes the Linux `core` bundle
(see [route 1](#route-1--tag-and-push-let-actions-build-it) for what it
actually runs):

```bash
gh run list --workflow "Build and release" --limit 3
gh release view v0.3.0
```

**5. Attach the Windows and macOS bundles**, built by hand on those machines
(see [Attaching the Windows and macOS bundles](#attaching-the-windows-and-macos-bundles)).

**6. Sanity checks** — before step 3, and after step 4:

```bash
git describe --tags --dirty --always   # what CMakeLists.txt will embed; -dirty = uncommitted changes in the build
git log --oneline $(git describe --tags --abbrev=0)..HEAD   # what is going into this release
gh release view v0.3.0                 # what ended up attached, and the notes
```

Went wrong? See [Undoing a mistake](#undoing-a-mistake) — and do not move a
tag anyone may have pulled.

## Cutting a release from the command line

Two routes. Both build the identical bundle, because both end up running
`tools/make_release.sh`. They differ in **who runs it** and, more importantly,
**who creates the tag and the GitHub release**.

|  | Route 1 — push a tag, Actions builds | Route 2 — build and publish locally |
|---|---|---|
| Runs `make_release.sh` | `ubuntu-22.04` runner | your machine |
| Version comes from | the tag name, passed by the workflow as `--version` | your `--version` flag |
| CLI profile | `core` (the script default; the workflow passes no `--profile`) | whatever you pass |
| Tag | **must already exist** — pushing it is the trigger | created by `gh` if absent, at default-branch tip |
| Creates the GitHub release | `softprops/action-gh-release` step in the workflow | `gh release create`, inside the script |
| Release notes | generated **plus a commit-author roll-call** (`git shortlog` since the previous tag) | `--generate-notes` only — no roll-call |
| Actions minutes | one Linux build | none |

Route 1 is the one to use for an official release: it produces the fuller
notes and it is the path that stays exercised. Route 2 is for prereleases and
one-off bundles where you do not want to spend Actions minutes.

### Route 1 — tag and push, let Actions build it

```bash
git tag -a v0.3.0 -m "lsdtt3 0.3.0"
git push origin v0.3.0
```

That matches the workflow's `on: push: tags: - "v*"` trigger. The workflow then
installs the build dependencies, runs

```bash
tools/make_release.sh --version 0.3.0      # VERSION="${GITHUB_REF_NAME#v}"
```

— i.e. **exactly the command you would run locally**, without `--publish` — and
a separate `softprops/action-gh-release` step creates the release and attaches
`dist/lsdtt3-backend-linux-*.tar.gz`. So `make_release.sh` builds and packages
in both routes; the difference is that in route 1 the *publishing* is done by
the workflow instead of by the script's `--publish`.

The workflow can also be run manually from the Actions tab (`workflow_dispatch`).
With no tag it builds at the current CMakeLists.txt version and uploads the
bundle as a 14-day workflow artifact — it does not create a release.

To push every local tag at once (rarely what you want — it publishes anything
you were still thinking about): `git push origin --tags`.

### Route 2 — build and publish locally, zero Actions minutes

**Do not pre-create a local tag for this route.** `gh release create` (which
`--publish` calls) says it plainly: *"If a matching git tag does not yet exist,
one will automatically get created from the latest state of the default
branch."* So:

```bash
git push                                   # make sure the default branch tip IS the commit you want
tools/make_release.sh --version 0.3.0-rc1 --profile all --publish --prerelease
git fetch --tags origin                    # pull down the tag gh just made
```

`--prerelease` keeps it from becoming "latest" for users. Re-running is safe:
if the release already exists the script uploads with `--clobber`.

> **Why not tag locally first?** If you create a local *annotated* tag and let
> gh create the remote one, the two are different objects, and pushing yours
> afterwards is **rejected** — verified: `! [rejected] v0.3.0 -> v0.3.0
> (already exists)` — even when both point at the same commit. You are then
> left with a lightweight remote tag and an unpushable local one. If you want a
> proper annotated tag on the remote, push the tag first and let route 1 do the
> build.

### Attaching the Windows and macOS bundles

Neither is built by CI: you build and package them by hand on that machine,
then attach the artefact to the release the Linux build already created.

**Each packaging script looks in a specific build directory**, and they do not
all use the same one. Keep to these conventions and the commands below work
unedited:

| Platform | Build directory | Script default | Profile | Artefact |
|---|---|---|---|---|
| Linux | `build/release` | `build` — but `make_release.sh` passes `--build-dir build/release` for you | `core` | `dist/lsdtt3-backend-linux-x86_64-core-v<version>.tar.gz` |
| Windows | `build/windows-qgis` | `build/windows-qgis` ✔ nothing to pass | `qgis` | `dist/lsdtt3-backend-win64-qgis-v<version>.zip` |
| macOS | `build/macos-core` | `build/macos-core` ✔ nothing to pass | `core` | `dist/lsdtt3-backend-macos-<arch>-core-v<version>.tar.gz` |

The build directory is separate per platform on purpose: the Windows tree is
configured with a different, deliberately minimal option set (no OpenMP, no
Eigen, no OpenCV), so reusing one `build/` for two profiles would silently ship
the wrong binaries.

#### Windows

In an **x64 Visual Studio Developer PowerShell**, from the repository root. Full
walkthrough (vcpkg, GDAL, the mandatory `-DLSDTT_OPENMP=OFF`, troubleshooting)
is in
[`docs/developer_notes/building_on_windows.md`](../docs/developer_notes/building_on_windows.md#native-windows-qgis-backend-only):

```powershell
# configure once per checkout
cmake -S . -B build\windows-qgis -G "Ninja" `
  -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_TOOLCHAIN_FILE="$env:USERPROFILE\tools\vcpkg\scripts\buildsystems\vcpkg.cmake" `
  -DVCPKG_TARGET_TRIPLET=x64-windows `
  -DBUILD_TESTING=OFF -DLSDTT_OPENMP=OFF -DLSDTT_EIGEN=OFF `
  -DLSDTT_OPENCV=OFF -DLSDTT_GPU=OFF `
  -DCMAKE_DISABLE_FIND_PACKAGE_PkgConfig=ON

# build just the CLIs the qgis profile ships
cmake --build build\windows-qgis --target `
  lsdtt-raster-preprocessing lsdtt-channel-extraction lsdtt-chi-analysis `
  lsdtt-basin-extraction lsdtt-basin-stats lsdtt-surface-metrics `
  lsdtt-channel-network-analysis

# package (BuildDir defaults to build/windows-qgis, CliProfile to qgis)
.\tools\package_windows_backend.ps1 -Version 0.3.0

gh release upload v0.3.0 dist\lsdtt3-backend-win64-qgis-v0.3.0.zip
```

`-DLSDTT_OPENMP=OFF` is not optional on MSVC — see the building-from-source note.
The script also expects vcpkg at `%USERPROFILE%/tools/vcpkg`; override with
`-VcpkgRoot`.

##### Packaging current code when your tree is dirty

If you just want a backend zip from your **current checkout** (including
uncommitted changes), package it with the CMake project version:

```powershell
# read project version from CMakeLists.txt: project(lsdtt3 VERSION X.Y.Z ...)
$ver = (Select-String -Path .\CMakeLists.txt -Pattern 'project\(lsdtt3 VERSION ([0-9.]+)' ).Matches[0].Groups[1].Value

.\tools\package_windows_backend.ps1 -Version $ver
```

This produces `dist\lsdtt3-backend-win64-qgis-v<version>.zip`, with
`backend_version: <version>` in `manifest.json`.

**Do not add a suffix** such as `-local-dirty` or `-rc1` to a bundle meant for
QGIS. The script accepts it, but the plugin only accepts a plain `X.Y.Z`
`backend_version` and refuses the bundle. The version must also be at least the
plugin's `MINIMUM_BACKEND_VERSION`. See
[`docs/developer_notes/qgis_plugin_development.md`](../docs/developer_notes/qgis_plugin_development.md)
for both checks. The binaries themselves still report a `-dirty` git
description with `-v`, so a dirty build remains identifiable.

#### macOS

Native macOS building and bundling are supported but **not yet integrated into CI** —
each release requires a manual native build on a Mac. See
[`docs/developer_notes/building_on_macos.md`](../docs/developer_notes/building_on_macos.md)
for the full walkthrough (Xcode CLT, Homebrew, CMake configuration, and bundling).

The `package_macos_backend.sh` script was written on a Linux machine and has not yet
run on real macOS hardware. The self-containment check fails hard, so the first run
will report exactly what is missing rather than shipping a broken bundle — but it may
find issues the Linux analogue does not anticipate.

```bash
cmake -S . -B build/macos-core -DCMAKE_BUILD_TYPE=Release
cmake --build build/macos-core -j"$(sysctl -n hw.ncpu)"

# BuildDir defaults to build/macos-core, cli-profile to core
tools/package_macos_backend.sh --version 0.3.0

gh release upload v0.3.0 dist/lsdtt3-backend-macos-arm64-core-v0.3.0.tar.gz
```

Bundles are architecture-specific: package on Apple Silicon for `arm64` users,
on Intel for `x86_64` users. If you need both architectures, build on both
machines separately and upload both tarballs with `gh release upload --clobber`.

#### Replacing a bundle already attached

```bash
gh release upload v0.3.0 dist/<bundle> --clobber
```

### Managing the release itself

```bash
gh release view   v0.3.0                  # what is attached, and the notes
gh release list                           # all releases
gh release create v0.3.0 --generate-notes # notes from merged PRs/commits
```

> `gh release edit` — for promoting a draft or an rc in place — **does not exist
> in gh 2.4.0**, which is what `apt install gh` gives on Ubuntu 22.04 (verified
> on this machine). `view`, `list`, `create`, `delete` and `upload --clobber` all
> work there. If you want `edit`, install a current gh from
> <https://cli.github.com>; otherwise promote a prerelease through the web UI, or
> delete and recreate it before anyone has fetched the tag.

### Undoing a mistake

Order matters: delete the release before the tag, or you leave a release
pointing at a tag that no longer exists.

```bash
gh release delete v0.3.0 --yes    # remove the GitHub release
git tag -d v0.3.0                 # delete locally
git push origin --delete v0.3.0   # delete on the remote
```

> **Do not move a tag that others may have pulled.** Re-pointing a published tag
> (`git tag -f` + `git push --force`) means two people can have different code
> under the same version, and anyone who already fetched it keeps the old commit
> silently. Cut `v0.3.1` instead. Deleting is only safe in the minutes after a
> mistake, before anyone has fetched.

### Checking what you are about to release

```bash
git describe --tags --dirty --always   # what CMakeLists.txt will embed
git tag -l                             # tags that already exist
git log --oneline $(git describe --tags --abbrev=0)..HEAD   # since last tag
```

A `-dirty` suffix means uncommitted changes are in the build. Release from a
clean tree.

## `make_release.sh` (the release driver)

Single entry point for producing a Linux release, used identically by you and
by GitHub Actions. In order it:

1. **Resolves the version** — from `--version`, or, if you omit it, from an
   exact `v*` tag on HEAD (`git describe --tags --exact-match`). With
   `--version` it does not consult git for the version at all.
2. **Checks it against `project(lsdtt3 VERSION ...)`** in CMakeLists.txt and
   aborts on a mismatch. A pre-release (`0.3.0-rc1`) compares only the `0.3.0`
   stem and warns instead of failing.
3. **Configures and builds** `build/release` (Release, OpenCV + Eigen on),
   unless `--skip-build`.
4. **Smoke-checks every CLI in the profile** by running it with `-v`; a
   non-zero exit fails the release.
5. **Packages** via `package_linux_backend.sh`.
6. **Publishes** with `gh`, only if you pass `--publish`.

What it does **not** do: create a tag, push a tag, move a tag, bump the
version in CMakeLists.txt, or commit anything. The only git operations are the
version lookup in step 1 and the dirty-tree check (a warning, not an error).

```bash
# local test build only (artifact in dist/), default core profile
tools/make_release.sh --version 0.3.0

# on an already-tagged commit: version taken from the tag
tools/make_release.sh

# local prerelease with every CLI, published to GitHub (zero Actions minutes)
tools/make_release.sh --version 0.3.0-rc1 --profile all --publish --prerelease

# repack an existing build without rebuilding
tools/make_release.sh --version 0.3.0 --skip-build --build-dir build
```

`tools/make_release.sh --help` lists every flag.

## CLI profiles (`cli_profiles.txt`)

Which CLIs go into a bundle is controlled by a profile; the qgis and core
lists live in `tools/cli_profiles.txt` (single source of truth for all four
scripts).

| Profile | Contents | Default for |
|---|---|---|
| `qgis` | small subset wired into the QGIS plugin | Windows bundles |
| `core` | tested tools for official releases | Linux and macOS bundles, CI |
| `all` | every `lsdtt-*` executable in the build dir (excl. `lsdtt-test-*`) | local releases only |

Select with `--profile` on `make_release.sh`, `--cli-profile` on the bash
packaging scripts, or `-CliProfile` on the PowerShell script. An explicit
`--cli` / `-CliNames` list overrides the profile (bundle labelled `custom`).
The profile appears in the bundle filename, e.g.
`lsdtt3-backend-linux-x86_64-core-v0.3.0.tar.gz`.

Windows and macOS bundles are built by hand on the respective machines
with the per-platform scripts below, then attached with
`gh release upload v<version> dist/<bundle>`.

There is one packaging script per platform, all producing the same bundle
layout (`bin/`, `lib/` on Unix, `share/gdal`, `share/proj`, `setenv.*`,
`README.md`, `manifest.json`, `THIRD_PARTY_NOTICES.txt`, `LICENSE`):

| Script | Platform | Artifact |
|---|---|---|
| `package_linux_backend.sh` | Linux x86_64 | `dist/lsdtt3-backend-linux-<arch>-<profile>-v<version>.tar.gz` |
| `package_windows_backend.ps1` | Windows x64 | `dist/lsdtt3-backend-win64-<profile>-v<version>.zip` |
| `package_macos_backend.sh` | macOS (arm64/x86_64) | `dist/lsdtt3-backend-macos-<arch>-<profile>-v<version>.tar.gz` *UNTESTED*|

Bundles are aimed at (often non-technical) end users: download, unpack, run
one activation command, and the tools work. Each bundle contains a generated
`README.md` with platform-specific install instructions (the Linux one also
covers Google Colab). Every script ends with a **self-containment check**
that fails the packaging if any staged binary has a dependency that is
neither bundled nor OS-provided — a broken bundle cannot ship silently.

## `package_linux_backend.sh`

### What it does

- copies selected `lsdtt-*` CLIs from the Linux build directory
- resolves the transitive shared-library closure with `ldd` and stages each
  library as a **regular file named by its SONAME** (the name the dynamic
  loader searches for, i.e. what `ldd` prints left of `=>`). No symlinks are
  used, so the bundle survives any filesystem or archiver. Core glibc
  libraries (`libc`, `libm`, `ld-linux`, ...) are deliberately excluded and
  come from the user's system
- sets rpath `$ORIGIN/../lib` on the executables when `patchelf` is available
- verifies every `DT_NEEDED` of every staged binary/library resolves inside
  `lib/` (via `objdump`); aborts otherwise
- copies GDAL and PROJ data (`share/gdal`, `share/proj`)
- writes `setenv.sh` (must be `source`d), end-user `README.md`,
  `manifest.json` with SHA256 checksums, `THIRD_PARTY_NOTICES.txt`
- builds the tar.gz (which includes the top-level bundle folder)

### Requirements

- Linux build host (a plain server or WSL both work; the script has no
  host-specific behaviour)
- completed build (default: `build`)
- `ldd`, `objdump`, `gdal-config`, `tar`, `python3`; `patchelf` optional but
  recommended

### Usage

```bash
tools/package_linux_backend.sh --version 0.2.0
tools/package_linux_backend.sh --version 0.2.0 --build-dir build/minimal-linux --cli lsdtt-chi-analysis
```

### End-user install (also in the bundle's README.md)

```bash
tar -xzf lsdtt3-backend-linux-x86_64-<profile>-v<version>.tar.gz
source lsdtt3-backend-linux-x86_64-<profile>-v<version>/setenv.sh   # per shell session
lsdtt-raster-preprocessing --help
```

`setenv.sh` must be **`source`d**, not run with `bash`: it only exports
environment variables (`PATH`, `LD_LIBRARY_PATH`, `GDAL_DATA`, `PROJ_DATA`),
which vanish with the child shell if executed.

### Google Colab

`!source setenv.sh` does not work in a notebook — every `!` line runs in its
own throwaway shell. Set the environment from Python once per session
instead:

```python
!wget https://www.geos.ed.ac.uk/~smudd/lsdtt_packages/lsdtt3-backend-linux-x86_64-<profile>-v<version>.tar.gz
!tar -xzf lsdtt3-backend-linux-x86_64-<profile>-v<version>.tar.gz
```

```python
import os

root = "/content/lsdtt3-backend-linux-x86_64-<profile>-v<version>"
os.environ["PATH"] = f"{root}/bin:" + os.environ["PATH"]
os.environ["LD_LIBRARY_PATH"] = f"{root}/lib:" + os.environ.get("LD_LIBRARY_PATH", "")
os.environ["GDAL_DATA"] = f"{root}/share/gdal"
os.environ["PROJ_DATA"] = f"{root}/share/proj"
os.environ["PROJ_LIB"]  = os.environ["PROJ_DATA"]
```

After that cell, `!lsdtt-raster-preprocessing --help` (etc.) works in any
later cell.

## `package_windows_backend.ps1`

### What it does

- copies selected `lsdtt-*.exe` tools from the Windows build directory
- resolves and copies dependent runtime DLLs from vcpkg (via `dumpbin`)
- bundles the **MSVC C/C++ runtime** (`vcruntime140.dll`, `msvcp140.dll`,
  ...) from the Visual Studio redist directory, so users do not need the
  VC++ Redistributable preinstalled
- verifies every import of every staged binary resolves inside `bin/` or is
  an OS-provided DLL (system DLLs, `api-ms-win-*` API sets, `ucrtbase`);
  aborts otherwise
- copies GDAL and PROJ data folders (`share/gdal`, `share/proj`)
- writes `setenv.ps1`, end-user `README.md`, `manifest.json` (UTF-8 without
  BOM, forward-slash paths), `THIRD_PARTY_NOTICES.txt`
- builds the zip via .NET `ZipFile` so the archive contains the top-level
  bundle folder and forward-slash entry names (PowerShell 5.1's
  `Compress-Archive` provides neither)

### Requirements

- x64 Visual Studio developer shell (for `dumpbin.exe` and
  `$env:VCToolsRedistDir`)
- completed Windows build (default: `build/windows-qgis`) — see
  `docs/developer_notes/building_on_windows.md` for the configure/build commands, and note
  `-DLSDTT_OPENMP=OFF` is required on MSVC
- vcpkg install with GDAL/PROJ runtime files (default:
  `%USERPROFILE%/tools/vcpkg`)

### Usage

```powershell
.\tools\package_windows_backend.ps1 -Version 0.2.0
.\tools\package_windows_backend.ps1 -Version 0.2.0 -BuildDir build/windows-qgis -OutDir dist
```

### End-user install (also in the bundle's README.md)

1. Unzip (right-click → Extract All...).
2. In PowerShell, `cd` into the extracted bundle folder.
3. `. .\setenv.ps1` (dot-source; repeat per PowerShell session). If script
   execution is blocked:
   `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass`
4. `lsdtt-raster-preprocessing --help`

## `package_macos_backend.sh`

**Status: UNTESTED** — written by analogy with the Linux script on a machine
without macOS. The self-containment check fails hard, so the first run on a
real Mac reports exactly what (if anything) is still wrong rather than
shipping a broken bundle.

### What it does

- copies selected CLIs from the macOS build directory
- resolves the transitive dylib closure with `otool -L` and stages each
  library as a regular file named by its install-name basename (no
  symlinks); `/usr/lib` and `/System` libraries are excluded (they live in
  the dyld shared cache and must not be copied)
- rewrites install names with `install_name_tool`: library references become
  `@rpath/<name>`, executables get rpath `@executable_path/../lib`, staged
  dylibs get id `@rpath/<name>` — so no `DYLD_*` environment variables are
  needed at runtime (SIP strips them anyway)
- re-signs every modified binary ad-hoc (`codesign --force --sign -`),
  required on Apple Silicon because `install_name_tool` invalidates
  signatures
- verifies every reference of every staged binary is either OS-provided or
  an `@rpath` reference satisfied in `lib/`; aborts otherwise
- copies GDAL/PROJ data; writes `setenv.sh` (PATH + data dirs only),
  end-user `README.md` (including the Gatekeeper quarantine step:
  `xattr -dr com.apple.quarantine <bundle>`), `manifest.json`,
  `THIRD_PARTY_NOTICES.txt`; builds the tar.gz

### Requirements

- macOS build host with Xcode command-line tools (`otool`,
  `install_name_tool`, `codesign`), `gdal-config`, `python3`
- completed macOS build in `build/macos-core` (the default; mirrors
  `build/windows-qgis` on Windows so the two platform trees never collide —
  override with `--build-dir`)
- bundles are architecture-specific: package on an Apple Silicon Mac for
  `arm64` users, on an Intel Mac (or with a cross-built tree and
  `--arch x86_64`) for `x86_64` users

### Usage

```bash
tools/package_macos_backend.sh --version 0.2.0
```

## Design notes (why it is done this way)

- **No symlinks in bundles.** Linux sonames and macOS install names are
  satisfied by staging the real file under the exact name the loader
  searches for. This keeps one workflow across build hosts and avoids the
  v0.2.0 failure where soname links were lost and users hit
  `libgdal.so.30: cannot open shared object file`.
- **Fail-hard verification.** Each script ends by checking the staged
  closure is complete. Silent gaps (the pre-fix scripts silently skipped
  unresolvable DLLs/libraries) are exactly how broken bundles reach users.
- **Bundle the MSVC runtime on Windows.** Without it, every exe fails with
  "VCRUNTIME140.dll was not found" on machines lacking the VC++
  Redistributable — i.e. most non-developer machines.
- **Archives contain the top-level bundle folder** on all platforms, so
  unpacking never scatters files into the user's current directory.
- **`manifest.json` is strict JSON** (UTF-8 without BOM, forward-slash
  paths) so the QGIS plugin and any Python tooling can parse it.
