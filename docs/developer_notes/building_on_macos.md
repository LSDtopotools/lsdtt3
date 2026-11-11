# Building on macOS

There are two ways to build lsdtt3 on a Mac:

| Route | Use it for |
|---|---|
| [Native build](#native-build) (Homebrew + Xcode Command Line Tools) | Development at full speed, and producing a macOS backend bundle |
| [Multipass Ubuntu VM](#multipass-ubuntu-vm) | Avoiding macOS-specific problems; matches the Linux reference environment exactly |

Linux is the reference platform and the native macOS build is less exercised
(see [Known gaps](#known-gaps)). If the native build gives you trouble, the VM
route always works.

> **Status (November 2026).** No maintainer builds on a Mac. The native route
> below is what the `build-macos` job of `.github/workflows/release.yml` runs on
> every release (on `macos-15` for Apple Silicon and `macos-15-intel`), so CI
> is where it is exercised: it has built and packaged the `qgis` bundles for
> both architectures since v0.9.1, and the macOS Python wheels since v0.9.2.
> The CMake handles Homebrew's keg-only `libomp` and pkg-config FFTW
> explicitly.

---

## Native build

### Prerequisites

- Apple Silicon or Intel Mac running macOS 11 (Big Sur) or later
- Xcode Command Line Tools (Apple Clang, Git)
- Homebrew

The build targets your machine's native architecture (`uname -m` reports
`arm64` on Apple Silicon, `x86_64` on Intel). Cross-compiling is not covered
here.

### 1. Install Xcode Command Line Tools

```bash
xcode-select --install
clang --version     # verify
```

### 2. Install Homebrew and the dependencies

Install Homebrew from [brew.sh](https://brew.sh) if you don't have it, then:

```bash
brew install cmake gdal fftw libomp netcdf eigen pkgconf
```

| Package | Required? | Role |
|---|---|---|
| `cmake` (≥ 3.15) | **yes** | Build system |
| `gdal` | **yes** | Raster and vector I/O, CRS handling |
| `libomp` | for multithreading | OpenMP runtime (see below) |
| `fftw` | no | Spectral analysis, wavelet polyfit |
| `pkgconf` | for `fftw` | CMake finds Homebrew's FFTW through pkg-config |
| `netcdf` | no | ERA40 / palaeomagnetic corrections in cosmo tools |
| `eigen` | no | k_sn / u\* inversion in chi analysis |

These are the Homebrew names for the same libraries listed in
[Building from source](building_from_source.md#install-dependencies).

**About `libomp`:** Apple Clang ships **no OpenMP runtime or headers at all**,
so `-fopenmp` is unusable out of the box. Homebrew's `libomp` supplies the
runtime, but it is keg-only (not linked into `/opt/homebrew`), so plain
`find_package(OpenMP)` cannot see it. `cmake/lsdtt_homebrew_openmp.cmake`
points CMake at the keg (from `$LIBOMP_PREFIX`, `brew --prefix libomp`, or the
standard prefixes) and prints `Apple Clang: using Homebrew libomp from ...`.
If you see `OpenMP not found - CPU operations will be single-threaded`
instead, install `libomp` or set `LIBOMP_PREFIX`.

### 3. Clone, configure and build

If you have not cloned from GitHub before, see [Git from the command
line](git_notes.md).

```bash
git clone git@github.com:lsdtopotools/lsdtt3.git
cd lsdtt3
cmake -S . -B build/macos -DCMAKE_BUILD_TYPE=Release
cmake --build build/macos -j"$(sysctl -n hw.ncpu)"
```

Read the configuration summary CMake prints (`OpenMP`, `FFTW3`, `Eigen3`): a
missing optional library means the corresponding tool is skipped, and the build
still succeeds. To make a bundle, see
[Packaging a macOS bundle](#packaging-a-macos-bundle), which uses its own
`build/macos-qgis` tree.

### 4. Put the tools on your PATH and test

```bash
echo 'export PATH=/path/to/lsdtt3/build/macos:$PATH' >> ~/.zshrc
source ~/.zshrc
lsdtt-raster-preprocessing --version

ctest --test-dir build/macos -L <domain> -LE heavy --output-on-failure
```

(`~/.zshrc` for the default zsh; use `~/.bashrc` if you use bash.)

### Known gaps

- **Only the `qgis` profile is packaged.** The macOS bundle exists to back the
  QGIS plugin. For the full CLI set on a Mac, use the
  [Multipass Ubuntu VM](#multipass-ubuntu-vm) and the Linux bundle.
- **Not run on a maintainer's Mac.** The build and packaging run in CI on
  every release; a local run that fails differently from CI is worth a bug
  report.

### Packaging a macOS bundle

The command sequence is in `tools/README.md` in the repository (see
[Making release packages](making_release_packages.md)); in short,
`tools/package_macos_backend.sh --version <version>` packages the `qgis` profile
from `build/macos-qgis` into
`dist/lsdtt3-backend-macos-<arch>-qgis-v<version>.tar.gz`. Add
`--tag v<version>` to build the release tag itself: the script exports it with
`git archive` into a scratch folder, configures and builds it there, and
packages that build. Bundles are
architecture-specific: package on Apple Silicon for `arm64` users and on Intel
for `x86_64` users. A bundle runs on the macOS release it was built on and
newer, because Homebrew's libraries target the build machine's release.

The script needs `otool`, `install_name_tool` and `codesign` (Xcode CLT),
`gdal-config` (Homebrew `gdal`) and `python3`. It rewrites library references
to `@rpath/...` and re-signs each modified binary ad hoc.

It deletes every rpath except `@executable_path/../lib` / `@loader_path`, so
the bundle never prefers a Homebrew library on the user's machine, runs each
CLI under `DYLD_PRINT_LIBRARIES=1` to prove nothing loads from outside the
bundle, and ends with a self-containment check. It **aborts** if a dependency
is not staged:

```
ERROR: <binary> needs '@loader_path/...' but it is not staged in lib/
ERROR: bundle is not self-contained; aborting
```

Do not work around it — this check exists because v0.2.0 shipped a bundle with
a missing library that failed on users' machines. To inspect a bundle by hand:

```bash
cd dist/lsdtt3-backend-macos-<arch>-qgis-v<version>
otool -L bin/lsdtt-raster-preprocessing   # non-system libs should be @rpath/...
source setenv.sh
lsdtt-raster-preprocessing --help
```

System libraries in `/usr/lib` and `/System/Library` do not appear as files —
they live in the dyld shared cache.

**Gatekeeper quarantine:** a bundle downloaded through a browser is marked
`com.apple.quarantine`, and unsigned binaries carrying it may be blocked. Clear
it with:

```bash
xattr -dr com.apple.quarantine lsdtt3-backend-macos-<arch>-qgis-v<version>
```

The bundle's own README tells end users the same. To use the bundle from QGIS,
see [QGIS plugin development](qgis_plugin_development.md).

---

## Multipass Ubuntu VM

[Multipass](https://multipass.run) is Canonical's lightweight Ubuntu VM manager.
Inside the VM you are on Linux, so the build is the reference one.

### 1. Install Multipass and launch an instance

```bash
brew install --cask multipass
multipass launch --name lsdtt --cpus 4 --memory 8G --disk 40G 22.04
multipass shell lsdtt
```

Adjust `--cpus`, `--memory` and `--disk` to your machine; 8 GB of RAM and 40 GB
of disk are comfortable. You are now at an Ubuntu 22.04 prompt. `exit` leaves
the shell; `multipass stop lsdtt` / `multipass start lsdtt` stop and restart
the VM.

### 2. Build

Inside the VM, follow [Linux](building_from_source.md#linux) in *Building from
source*. For cloning, [Git from the command line](git_notes.md) covers creating
an SSH key inside the VM and adding it to GitHub.

### 3. Share files between macOS and the VM

Multipass mounts your macOS home directory in the VM at `/home/ubuntu/Home` by
default. To mount another folder, run this **from macOS**, not inside the VM:

```bash
multipass mount /Users/yourname/lsdtt_data lsdtt:/lsdtt_data
multipass umount lsdtt:/lsdtt_data      # to unmount
```

That keeps your DEMs on the Mac while the tools run in the VM.

---

## See also

- [Building from source](building_from_source.md) — the reference Linux build and CMake options
- [Making release packages](making_release_packages.md) — how macOS bundles fit into a release
- [QGIS plugin development](qgis_plugin_development.md) — using a backend bundle from QGIS
