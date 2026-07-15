# Installing lsdtt3

lsdtt3 runs on Linux. This page describes the three ways to get it running and
points to the detailed instructions for each platform.

---

## Which route is right for you?

| I want to… | Route |
|---|---|
| Try lsdtt3 without compiling anything | [Pre-built install (beta)](#pre-built-install-beta) |
| Run lsdtt3 on a Windows machine | [Windows (WSL2)](#windows-wsl2) |
| Build a limited native Windows backend for QGIS plugin development | [Windows native (developers only, limited)](#windows-native-developers-only-limited) |
| Run lsdtt3 on a Mac | [macOS (Multipass)](#macos-multipass) |
| Run lsdtt3 on a Linux machine | [Build from source on Linux](#build-from-source-on-linux) |
| Develop or modify lsdtt3 | [Build from source on Linux](#build-from-source-on-linux) |
| Create a redistributable Linux bundle | [Build from source on Linux](#build-from-source-on-linux) |

---

## Pre-built install (beta)

Pre-built executables for Ubuntu 22.04 are attached to each
[GitHub Release](https://github.com/LSDtopotools/lsdtt3/releases/latest).
You do not need to install a compiler or manage any build dependencies.

Download `lsdtt3-ubuntu-22.04.tar.gz` from the Releases page and run:

```bash
bash install.sh /path/to/lsdtt3-ubuntu-22.04.tar.gz
```

The script installs the required runtime libraries via `apt`, unpacks all
command-line tools to `~/.local/lsdtt3/bin`, and adds them to your `PATH`.
Open a new terminal afterwards and run `lsdtt-channel-extraction --help` to
verify.

Full details, troubleshooting, and instructions for reporting bugs are in
[`docs/beta_testing_notes.md`](beta_testing_notes.md).

---

## Windows (WSL2)

lsdtt3 runs inside **Windows Subsystem for Linux (WSL2)**, which gives you a
full Ubuntu 22.04 environment without leaving Windows.

**For beta testers (no compilation needed):**  
Install WSL2 and Ubuntu, then follow the [pre-built install](#pre-built-install-beta)
steps above. The `install.sh` script works identically inside WSL2.

**For developers who want to build from source:**  
See [`docs/lsdtt3_in_windows.md`](lsdtt3_in_windows.md) for the complete
walkthrough:

1. Install WSL2 and Ubuntu 22.04
2. Install build tools and libraries (`apt`)
3. Set up an SSH key and clone the repository
4. Build with CMake
5. Add the executables to your PATH
6. Access Windows files from the Ubuntu environment

---

## Windows native (developers only, limited)

This path is intended for plugin/backend developers who need a native Windows build for a subset of tools used by the QGIS plugin. It is not the primary development route.

Important limitations:

- Full `lsdtt3` functionality is available on Linux (including WSL2).
- Native Windows builds are currently treated as a reduced-feature profile.
- Optional components such as Eigen/OpenCV may be disabled in this mode;
  FFTW3 is included (via vcpkg) so spectral analysis and Wiener channel
  extraction work on Windows.
- Testing suggests that the Linux version run within a windows computer actually runs faster than versions compiled for Windows and run on Windows. So these builds are really only for use as the backend for the QGIS pluging. 

### What you need

- Visual Studio Community 2022 with C++ build tools
- Developer shell from Visual Studio (`Developer PowerShell for VS 2022`)
- CMake and Ninja available in that shell
- `vcpkg` in your user folder to install dependencies without admin rights

### Open Developer PowerShell

Open Start menu and launch `Developer PowerShell for VS 2022`.

> [!WARNING]
> Some Visual Studio shell entries start in an `x86` toolchain context. If this
> happens, CMake uses a 32-bit compiler and will reject `gdal:x64-windows` as
> incompatible.
>
> Always verify `cl` reports `for x64` before configuring.

If Start menu shortcuts are unavailable, open PowerShell and run:

```powershell
& "$env:ProgramW6432\Microsoft Visual Studio\2022\Community\Common7\Tools\Launch-VsDevShell.ps1" -Arch amd64 -HostArch amd64
```

If `cl` still reports `for x86`, switch to x64 explicitly:

```powershell
& "$env:ProgramW6432\Microsoft Visual Studio\2022\Community\Common7\Tools\Launch-VsDevShell.ps1" -Arch amd64 -HostArch amd64
cl
```

You can also launch `x64 Native Tools Command Prompt for VS 2022`, then run
PowerShell commands from there.

Verify toolchain setup:

```powershell
cmake --version
ninja --version
cl
```

Expected compiler line:

```text
Microsoft (R) C/C++ Optimizing Compiler Version ... for x64
```

### Install GDAL via vcpkg (no admin)


This installs required dependencies for `lsdtt3` on windows. 


> [!WARNING]
> This process takes a very long time: on a standard laptop in 2026 it took the author of this documentation 2 hours.
 
To install `vcpkg` in your user profile and use it to install GDAL:

```powershell
git clone https://github.com/microsoft/vcpkg "$env:USERPROFILE\tools\vcpkg"
& "$env:USERPROFILE\tools\vcpkg\bootstrap-vcpkg.bat"
& "$env:USERPROFILE\tools\vcpkg\vcpkg.exe" install gdal:x64-windows fftw3:x64-windows
```

Notes:

- GDAL is required.
- FFTW3 is included in the Windows profile: it enables `lsdtt-spectral-analysis`
  and the Wiener source-extraction option in `lsdtt-channel-extraction`.
  CMake finds the vcpkg fftw3 port via its CMake config package
  (`FFTW3::fftw3`), so pkg-config is not needed on Windows.
- If Visual Studio sets a conflicting `VCPKG_ROOT`, clear it in-session with
  `Remove-Item Env:VCPKG_ROOT`.

### Configure and build minimal Windows profile

From the repository root:

```powershell
cmake -S . -B build\minimal-windows -G "Ninja" `
  -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_TOOLCHAIN_FILE="$env:USERPROFILE\tools\vcpkg\scripts\buildsystems\vcpkg.cmake" `
  -DVCPKG_TARGET_TRIPLET=x64-windows `
  -DLSDTT_EIGEN=OFF `
  -DLSDTT_OPENCV=OFF `
  -DLSDTT_GPU=OFF `
  -DCMAKE_DISABLE_FIND_PACKAGE_PkgConfig=ON `
  -DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=ON

cmake --build build\minimal-windows
```

Note: disabling PkgConfig does not disable FFTW3 - on Windows FFTW3 is found
through the vcpkg CMake config package. The configure summary should report
`FFTW3       : ON`.

### Troubleshooting (Windows native)

- If configure fails with GDAL "configuration files were considered but not
  accepted" and mentions `64bit`, your shell is likely using an x86 toolchain.
  Re-open an x64 Visual Studio shell and confirm `cl` reports `for x64`.
- If build fails with `cl : Command line error D8021 : invalid numeric argument
  '/Wextra'`, update to a revision that includes MSVC-specific warning flags in
  `CMakeLists.txt` (`/W4` on MSVC instead of GCC-style `-Wextra`).

Optional verification:

```powershell
Get-ChildItem build\minimal-windows\bin\lsdtt-*.exe
.\build\minimal-windows\bin\lsdtt-raster-preprocessing.exe --version
```

### Package a redistributable Windows backend zip

After a successful build, create a runtime bundle that users can run without
installing Visual Studio, vcpkg, or GDAL themselves.

From the repository root (x64 Visual Studio developer shell):

```powershell
.\tools\package_windows_backend.ps1 -Version 0.2.0
```

This creates:

- `dist/lsdtt3-backend-win64-minimal-v0.2.0.zip` (the file to distribute)
- `dist/lsdtt3-backend-win64-minimal-v0.2.0/` (staging folder used to build the zip)

What the bundle includes:

- `bin/` with the MVP CLIs (`lsdtt-raster-preprocessing`,
  `lsdtt-channel-extraction`, `lsdtt-chi-analysis`,
  `lsdtt-basin-extraction`) plus `lsdtt-spectral-analysis` when the build
  has FFTW support
- required runtime DLLs discovered from those executables (including
  `fftw3.dll`)
- `share/gdal` and `share/proj` data folders (including `proj.db`)
- `manifest.json`, `LICENSE`, `THIRD_PARTY_NOTICES.txt`, and `setenv.ps1`

How recipients run it:

```powershell
cd <unzipped-backend-folder>
. .\setenv.ps1
.\bin\lsdtt-raster-preprocessing.exe --version
```

Advanced options:

```powershell
.\tools\package_windows_backend.ps1 `
  -Version 0.2.0 `
  -BuildDir build/minimal-windows `
  -OutDir dist
```

This profile is the recommended starting point for native Windows work related
to QGIS plugin packaging/testing. For complete capability and general lsdtt3
development, use Linux (native, WSL2, or Multipass).

---

## macOS (Multipass)

lsdtt3 does not build natively on macOS (Apple Clang lacks several required
features). The recommended approach is to run it inside an Ubuntu 22.04 virtual
machine managed by **Multipass**, a lightweight tool from Canonical.

See [`docs/lsdtt3_in_macos.md`](lsdtt3_in_macos.md) for the complete
walkthrough:

1. Install Homebrew and Multipass
2. Launch an Ubuntu 22.04 instance
3. Install build tools and libraries inside the VM
4. Set up an SSH key and clone the repository
5. Build with CMake
6. Share your macOS data folders into the VM

---

## Build from source on Linux

If you are on a native Linux system (or already inside WSL2 or Multipass), the
build process is straightforward.

### Install dependencies

```bash
sudo apt update
sudo apt install -y \
    build-essential cmake pkg-config \
    libgdal-dev \
    libfftw3-dev \
    libopencv-dev \
    libeigen3-dev \
    libnetcdf-dev
```

| Package | Role |
|---|---|
| `build-essential` | GCC C++ compiler and standard build tools |
| `cmake` | Build system (version 3.15+ required) |
| `pkg-config` | Lets CMake locate FFTW3 and other libraries |
| `libgdal-dev` | **Required** — raster I/O and CRS handling |
| `libfftw3-dev` | Optional — spectral analysis tools |
| `libopencv-dev` | Optional — image-analysis routines in surface metrics |
| `libeigen3-dev` | Optional — chi/k_sn inversion tools |
| `libnetcdf-dev` | Optional — ERA40 atmosphere correction in cosmo tools |

OpenMP is included with GCC and needs no separate install.

### Clone the repository

```bash
git clone git@github.com:lsdtopotools/lsdtt3.git
cd lsdtt3
```

### Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

CMake prints a summary of which optional libraries were found. Tools that
depend on a missing library are skipped — the build does not fail.

### Add the executables to your PATH

```bash
echo 'export PATH=/path/to/lsdtt3/build:$PATH' >> ~/.bashrc
source ~/.bashrc
```

Replace `/path/to/lsdtt3` with the actual path to the repository. Verify with:

```bash
lsdtt-channel-extraction --help
```

### Package a redistributable Linux backend tarball

If you want to distribute lsdtt3 on Linux without asking users to install all
runtime dependencies, create a bundle that includes executables, required
shared libraries (`.so` files), and GDAL/PROJ data.

From the repository root:

```bash
chmod +x tools/package_linux_backend.sh
./tools/package_linux_backend.sh --version 0.2.0
```

This creates:

- `dist/lsdtt3-backend-linux-x86_64-minimal-v0.2.0.tar.gz` (file to distribute)
- `dist/lsdtt3-backend-linux-x86_64-minimal-v0.2.0/` (staging folder)

Bundle contents:

- `bin/` with selected `lsdtt-*` CLIs
- `lib/` with transitive runtime `.so` dependencies (except glibc loader/core)
- `share/gdal` and `share/proj` (including `proj.db`)
- `setenv.sh`, `manifest.json`, `LICENSE`, `THIRD_PARTY_NOTICES.txt`

How recipients run it:

```bash
cd <unzipped-backend-folder>
source ./setenv.sh
./bin/lsdtt-raster-preprocessing --version
```

Compatibility note:

- Linux uses `.so` files as the shared-library equivalent of Windows `.dll`.
- Build on an older baseline distro (for example Ubuntu 22.04) to maximise
  compatibility with newer Linux systems.
- This bundle does not package glibc itself; kernel/glibc compatibility still
  applies.

### CMake options

| Option | Default | Effect |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | Use `Debug` for development or `RelWithDebInfo` for profiling |
| `LSDTT_OPENCV` | `ON` | Build OpenCV-dependent routines (auto-disabled if OpenCV not found) |
| `LSDTT_EIGEN` | `ON` | Build Eigen3-dependent routines (auto-disabled if Eigen3 not found) |
| `LSDTT_GPU` | `OFF` | Enable GPU backend (requires CUDA toolkit) |
