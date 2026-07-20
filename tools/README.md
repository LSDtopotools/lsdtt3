# tools

Utility scripts for building and packaging redistributable backend bundles.

## `make_release.sh` (the release driver)

Single entry point for producing a Linux release: resolves the version
(from `--version` or an exact `v*` tag on HEAD), checks it against
`project(lsdtt3 VERSION ...)` in CMakeLists.txt, configures and builds,
smoke-checks every release CLI (`-v`), packages via
`package_linux_backend.sh`, and optionally publishes with `gh`.

GitHub Actions (`.github/workflows/release.yml`) calls this same script on
`v*` tag pushes, so a local run exercises exactly the CI path.

```bash
# local test build only (artifact in dist/), default core profile
tools/make_release.sh --version 0.3.0

# local prerelease with every CLI, published to GitHub (zero Actions minutes)
tools/make_release.sh --version 0.3.0-rc1 --profile all --publish --prerelease

# repack an existing build without rebuilding
tools/make_release.sh --version 0.3.0 --skip-build --build-dir build
```

## Cutting a release from the command line

Two routes. Both produce the same artefacts, because `.github/workflows/release.yml`
just calls `tools/make_release.sh` — the difference is only who runs it.

### Before anything else: bump the version in CMakeLists.txt

`make_release.sh` **refuses to run** if the version it is releasing does not
match `project(lsdtt3 VERSION ...)` — so this is step one, not an afterthought:

```bash
# CMakeLists.txt line 2:  project(lsdtt3 VERSION 0.3.0 LANGUAGES CXX)
git add CMakeLists.txt
git commit -m "Bump version to 0.3.0"
git push
```

The tag must be exactly `v` + that version (`v0.3.0` for `VERSION 0.3.0`).
CMakeLists.txt uses the same rule to decide whether a build reports a plain
release version or a `-dev` one, so a mismatch silently produces `-dev` binaries
in what looks like a release.

### Route 1 — tag and push, let Actions build it

```bash
# annotated tag (-a with -m), not lightweight: it carries an author, date and
# message, and `git describe` prefers it
git tag -a v0.3.0 -m "lsdtt3 0.3.0"

# push the tag on its own; pushing the branch does NOT push tags
git push origin v0.3.0
```

That matches the workflow's `on: push: tags: - "v*"` trigger, which runs
`make_release.sh --version 0.3.0` and publishes the Linux `core` bundle.

To push every local tag at once (rarely what you want — it publishes anything
you were still thinking about): `git push origin --tags`.

### Route 2 — build and publish locally, zero Actions minutes

```bash
git tag -a v0.3.0 -m "lsdtt3 0.3.0"
tools/make_release.sh --version 0.3.0 --publish
git push origin v0.3.0
```

`--publish` creates the GitHub release and uploads the bundle via `gh`
(re-running is safe — it uploads with `--clobber` if the release exists). Push
the tag afterwards so the release and the tag agree; if you push first, Actions
will also build it and you get two runs producing the same artefact.

For a prerelease, which does not become "latest" for users:

```bash
git tag -a v0.3.0-rc1 -m "lsdtt3 0.3.0-rc1"
tools/make_release.sh --version 0.3.0-rc1 --profile all --publish --prerelease
git push origin v0.3.0-rc1
```

### Attaching the Windows and macOS bundles

Those are built by hand on the respective machines (see the per-platform scripts
below), then attached to the existing release:

```bash
gh release upload v0.3.0 dist/lsdtt3-backend-win64-qgis-v0.3.0.zip
gh release upload v0.3.0 dist/lsdtt3-backend-macos-arm64-core-v0.3.0.tar.gz

# replacing a bundle already attached
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

Full release workflow and rationale: `docs/release_ci_plan.md`.
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
- completed Windows build (default: `build/minimal-windows`)
- vcpkg install with GDAL/PROJ runtime files (default:
  `%USERPROFILE%/tools/vcpkg`)

### Usage

```powershell
.\tools\package_windows_backend.ps1 -Version 0.2.0
.\tools\package_windows_backend.ps1 -Version 0.2.0 -BuildDir build/minimal-windows -OutDir dist
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
- completed macOS build (default: `build`)
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
