# Building lsdtt3 from source

This is the developer route. If you only want to *use* lsdtt3, install a
release package instead — see [`docs/install_lsdtt3.md`](../install_lsdtt3.md).
Build from source when you need code that has not been released yet, or when
you are modifying lsdtt3.

lsdtt3 is a Linux-native C++17 codebase. **Linux is the reference platform**
and is the only one where every tool builds. The other sections below get you
to a Linux environment, or explain the deliberately reduced native Windows
build.

| Your machine | Route |
|---|---|
| Linux | [Linux](#linux) |
| Windows | [Windows via WSL2](#windows-via-wsl2) — full build |
| Windows, native `.exe` needed | [Windows native (QGIS backend only)](#windows-native-qgis-backend-only) — reduced |
| macOS | [macOS via Multipass](#macos-via-multipass) |

---

## Linux

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

| Package | Required? | Role |
|---|---|---|
| `build-essential` | **yes** | GCC C++17 compiler and standard build tools |
| `cmake` (≥ 3.15) | **yes** | Build system |
| `pkg-config` | **yes** | Lets CMake locate FFTW3 and others |
| `libgdal-dev` (≥ 3.0) | **yes** | Raster and vector I/O, CRS handling |
| `libfftw3-dev` | no | `lsdtt-spectral-analysis`, wavelet quadratic polyfit |
| `libopencv-dev` | no | Frangi filter, Otsu segmentation, Hough lines in `lsdtt-surface-metrics` |
| `libeigen3-dev` | no | k_sn / u\* inversion in `lsdtt-chi-analysis` |
| `libnetcdf-dev` | no | ERA40 / palaeomagnetic corrections in the cosmo tools |

OpenMP ships with GCC and needs no separate install.

Optional dependencies are detected automatically. A tool that needs a missing
library is **skipped**, and the build still succeeds — so a partial dependency
set gives you a partial tool set rather than an error.

### Clone and build

```bash
git clone git@github.com:lsdtopotools/lsdtt3.git
cd lsdtt3
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

CMake prints a summary of which optional libraries it found. Read it — it is
the quickest way to discover that, say, Eigen was missing and chi inversion
was silently dropped.

### Put the tools on your PATH

```bash
echo 'export PATH=/path/to/lsdtt3/build:$PATH' >> ~/.bashrc
source ~/.bashrc
lsdtt-channel-extraction --help
```

### CMake options

| Option | Default | Effect |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Debug` for development, `RelWithDebInfo` for profiling |
| `LSDTT_OPENCV` | `ON` | OpenCV-dependent routines (auto-off if OpenCV is absent) |
| `LSDTT_EIGEN` | `ON` | Eigen3-dependent routines (auto-off if Eigen3 is absent) |
| `LSDTT_OPENMP` | `ON` | OpenMP parallelism. **Must be `OFF` on MSVC** — see below |
| `LSDTT_GPU` | `OFF` | GPU compute backend (needs the CUDA toolkit) |
| `BUILD_TESTING` | `ON` | Build the CTest suite |
| `LSDTT_ENABLE_BENCHMARKS` | `OFF` | RichDEM comparison benchmarks (needs network + pixi) |
| `LSDTT_ENABLE_NETWORK_TESTS` | `OFF` | Tests that download data |
| `LSDTT_ENABLE_PERF_TESTS` | `OFF` | Flow-routing performance benchmarks |

### Run the tests

```bash
ctest --test-dir build --output-on-failure
```

The default suite needs no network, no downloaded data and no Python packages.
See `tests/README.md` in the repository for what each test covers, and
[`docs/testing_and_benchmarking.md`](../testing_and_benchmarking.md) for
methodology.

---

## Windows via WSL2

WSL2 gives you a full Ubuntu environment, so you build exactly as in
[Linux](#linux) above. This is the recommended way to develop lsdtt3 on a
Windows machine — every tool builds, and the resulting binaries run faster
than a native Windows build on the same hardware.

[`lsdtt3_in_windows.md`](lsdtt3_in_windows.md) has the full walkthrough:
installing WSL2 and Ubuntu 22.04, apt dependencies, SSH keys and cloning,
building with CMake, PATH setup, and reaching your Windows files from Ubuntu.

---

## Windows native (QGIS backend only)

This path produces native `.exe` files for the subset of tools the QGIS plugin
drives. It is **not** the primary development route.

Limitations, up front:

- Full lsdtt3 functionality exists only on Linux (including WSL2).
- Native Windows is a reduced-feature profile: Eigen and OpenCV are off.
  FFTW3 *is* included via vcpkg, so spectral analysis and Wiener channel
  extraction work.
- Benchmarking suggests the Linux build run under WSL2 on a Windows machine is
  faster than a native Windows build on the same machine. These builds exist
  for the QGIS plugin backend, not for performance.

### What you need

- Visual Studio Community 2022 with the C++ build tools
- `Developer PowerShell for VS 2022`
- CMake and Ninja on the path in that shell
- `vcpkg` in your user folder, so dependencies install without admin rights

### Open an x64 developer shell

Launch `Developer PowerShell for VS 2022` from the Start menu.

> [!WARNING]
> Some Visual Studio shell entries start in an **x86** toolchain context. CMake
> then picks a 32-bit compiler and rejects `gdal:x64-windows` as incompatible.
> Always confirm `cl` reports `for x64` before configuring.

If the Start menu shortcut is unavailable, or `cl` reports x86:

```powershell
& "$env:ProgramW6432\Microsoft Visual Studio\2022\Community\Common7\Tools\Launch-VsDevShell.ps1" -Arch amd64 -HostArch amd64
cl
```

`x64 Native Tools Command Prompt for VS 2022` works too. Verify the toolchain:

```powershell
cmake --version
ninja --version
cl        # expect: ... Optimizing Compiler Version ... for x64
```

### Install GDAL and FFTW3 via vcpkg

> [!WARNING]
> This takes a long time — around 2 hours on a standard 2026 laptop.

```powershell
git clone https://github.com/microsoft/vcpkg "$env:USERPROFILE\tools\vcpkg"
& "$env:USERPROFILE\tools\vcpkg\bootstrap-vcpkg.bat"
& "$env:USERPROFILE\tools\vcpkg\vcpkg.exe" install gdal:x64-windows fftw3:x64-windows
```

Notes:

- GDAL is required; FFTW3 is optional but cheap (one extra DLL).
- CMake finds the vcpkg FFTW3 port through its CMake config package
  (`FFTW3::fftw3`), so pkg-config is not needed on Windows.
- If Visual Studio has set a conflicting `VCPKG_ROOT`, clear it for the
  session with `Remove-Item Env:VCPKG_ROOT`.

### Configure and build

From the repository root:

```powershell
cmake -S . -B build\windows-qgis -G "Ninja" `
  -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_TOOLCHAIN_FILE="$env:USERPROFILE\tools\vcpkg\scripts\buildsystems\vcpkg.cmake" `
  -DVCPKG_TARGET_TRIPLET=x64-windows `
  -DBUILD_TESTING=OFF `
  -DLSDTT_OPENMP=OFF `
  -DLSDTT_EIGEN=OFF `
  -DLSDTT_OPENCV=OFF `
  -DLSDTT_GPU=OFF `
  -DCMAKE_DISABLE_FIND_PACKAGE_PkgConfig=ON

cmake --build build\windows-qgis
```

Or build only the CLIs the QGIS bundle needs, rather than the whole tree:

```powershell
cmake --build build\windows-qgis --target `
  lsdtt-raster-preprocessing lsdtt-channel-extraction lsdtt-chi-analysis `
  lsdtt-basin-extraction lsdtt-basin-stats lsdtt-surface-metrics `
  lsdtt-channel-network-analysis
```

**`-DLSDTT_OPENMP=OFF` is not optional on MSVC.** Visual Studio's `/openmp`
selects OpenMP **2.0** (2002), which rejects unsigned loop indices, silently
ignores `collapse`, and has no min/max reductions — all of which this codebase
uses. Leaving it on fails with `C3016`, `C4849` and `C7660` part-way through
`lsdtt3_cpu`. CMake warns about this at configure time.

The loss is smaller than it sounds: `C4849` means `collapse` was being ignored
anyway, so MSVC was only ever parallelising outer loops. `/openmp:llvm` would
give real OpenMP 3.1+, but it switches the runtime from `vcomp140.dll` to
`libomp140.x86_64.dll`, which `package_windows_backend.ps1` does not bundle —
so that is a deliberate change, not a flag flip.

**`-DBUILD_TESTING=OFF`** drops test targets the bundle does not need and which
pull in code only exercised on Linux.

Disabling PkgConfig does **not** disable FFTW3; the configure summary should
still report `FFTW3       : ON`.

### Verify

```powershell
Get-ChildItem build\windows-qgis\lsdtt-*.exe
.\build\windows-qgis\lsdtt-raster-preprocessing.exe --version
```

### Troubleshooting (Windows native)

- GDAL configure error mentioning "configuration files were considered but not
  accepted" and `64bit`: your shell is on an x86 toolchain. Re-open an x64
  developer shell and confirm `cl` reports `for x64`.
- `cl : Command line error D8021 : invalid numeric argument '/Wextra'`: update
  to a revision whose `CMakeLists.txt` uses MSVC warning flags (`/W4`) rather
  than GCC-style `-Wextra`.

---

## macOS via Multipass

lsdtt3 does not build with Apple Clang — it lacks several features the code
relies on. Build inside an Ubuntu 22.04 VM managed by
[Multipass](https://multipass.run), then follow [Linux](#linux) inside it.

[`lsdtt3_in_macos.md`](lsdtt3_in_macos.md) has the full walkthrough: installing
Homebrew and Multipass, launching the instance, dependencies, SSH keys and
cloning, building, and sharing macOS folders into the VM.

---

## See also

- [`developer_onboarding.md`](developer_onboarding.md) — start here as a new contributor
- [`lsdtt3_architecture.md`](lsdtt3_architecture.md) — layer design, naming conventions, module APIs
- [`making_release_packages.md`](making_release_packages.md) — turning a build into a distributable package
- [`debian_packaging.md`](debian_packaging.md) — building a `.deb`
