# Building on macOS

There are two ways to build lsdtt3 on a Mac:

| Route | Use it for |
|---|---|
| [Native build](#native-build) (Homebrew + Xcode Command Line Tools) | Development at full speed, and producing a macOS backend bundle |
| [Multipass Ubuntu VM](#multipass-ubuntu-vm) | Avoiding macOS-specific problems; matches the Linux reference environment exactly |

Linux is the reference platform and the native macOS build is less exercised
(see [Known gaps](#known-gaps)). If the native build gives you trouble, the VM
route always works.

**WARNING** None of this has been tested as of September 2026

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
brew install cmake gdal fftw libomp netcdf eigen
```

| Package | Required? | Role |
|---|---|---|
| `cmake` (≥ 3.15) | **yes** | Build system |
| `gdal` | **yes** | Raster and vector I/O, CRS handling |
| `libomp` | for multithreading | OpenMP runtime (see below) |
| `fftw` | no | Spectral analysis, wavelet polyfit |
| `netcdf` | no | ERA40 / palaeomagnetic corrections in cosmo tools |
| `eigen` | no | k_sn / u\* inversion in chi analysis |

These are the Homebrew names for the same libraries listed in
[Building from source](building_from_source.md#install-dependencies).

**About `libomp`:** Apple Clang ships **no OpenMP runtime or headers at all**,
so `-fopenmp` is unusable out of the box, `find_package(OpenMP)` fails, and the
build falls back to single-threaded (CMake prints
`OpenMP not found - CPU operations will be single-threaded`). Homebrew's
`libomp` supplies the runtime — but see [Known gaps](#known-gaps): it is
keg-only, and CMake probably will not find it without help.

### 3. Clone, configure and build

If you have not cloned from GitHub before, see [Git from the command
line](git_notes.md).

```bash
git clone git@github.com:lsdtopotools/lsdtt3.git
cd lsdtt3
cmake -S . -B build/macos-core -DCMAKE_BUILD_TYPE=Release
cmake --build build/macos-core -j"$(sysctl -n hw.ncpu)"
```

`build/macos-core` is the directory `tools/package_macos_backend.sh` reads by
default. Read the configuration summary CMake prints: a missing optional
library means the corresponding tool is skipped, and the build still succeeds.

### 4. Put the tools on your PATH and test

```bash
echo 'export PATH=/path/to/lsdtt3/build/macos-core:$PATH' >> ~/.zshrc
source ~/.zshrc
lsdtt-raster-preprocessing --version

ctest --test-dir build/macos-core --output-on-failure
```

(`~/.zshrc` for the default zsh; use `~/.bashrc` if you use bash.)

### Known gaps

- **OpenMP discovery is unverified.** Homebrew's `libomp` is keg-only, so it is
  not linked into the default search paths, and `CMakeLists.txt` has no
  `APPLE` branch pointing at it. Expect to pass the libomp paths explicitly, or
  configure with `-DLSDTT_OPENMP=OFF` for a working single-threaded build.
- **Packaging is untested on real hardware.** `tools/package_macos_backend.sh`
  was written by analogy with the Linux script. See [Packaging a macOS
  bundle](#packaging-a-macos-bundle).

### Packaging a macOS bundle

The command sequence is in `tools/README.md` in the repository (see
[Making release packages](making_release_packages.md)); in short,
`tools/package_macos_backend.sh --version <version>` produces
`dist/lsdtt3-backend-macos-<arch>-core-v<version>.tar.gz`. Bundles are
architecture-specific: package on Apple Silicon for `arm64` users and on Intel
for `x86_64` users.

The script needs `otool`, `install_name_tool` and `codesign` (Xcode CLT),
`gdal-config` (Homebrew `gdal`) and `python3`. It rewrites library references
to `@rpath/...` and re-signs each modified binary ad hoc.

It ends with a self-containment check and **aborts** if a dependency is not
staged:

```
ERROR: <binary> needs '@loader_path/...' but it is not staged in lib/
ERROR: bundle is not self-contained; aborting
```

Do not work around it — this check exists because v0.2.0 shipped a bundle with
a missing library that failed on users' machines. To inspect a bundle by hand:

```bash
cd dist/lsdtt3-backend-macos-<arch>-core-v<version>
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
xattr -dr com.apple.quarantine lsdtt3-backend-macos-<arch>-core-v<version>
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
