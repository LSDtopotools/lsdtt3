# Building on Windows

There are two ways to build lsdtt3 on a Windows machine:

| Route | Tool set | Use it for |
|---|---|---|
| [WSL2](#wsl2-recommended) | **Full** — identical to Linux | Development, analysis, everything |
| [Native Windows](#native-windows-qgis-backend-only) (Visual Studio + vcpkg) | **Reduced** — no OpenMP or Eigen | Producing `.exe` files for the QGIS plugin backend |

Use WSL2 unless you specifically need native `.exe` files. Benchmarking
suggests the Linux build under WSL2 is also *faster* than a native Windows
build on the same machine.

---

## WSL2 (recommended)

WSL2 gives you a full Ubuntu environment inside Windows. Once it is installed
you build exactly as on Linux.

### 1. Install WSL2 and Ubuntu

1. Open PowerShell with administrator privileges. On a managed machine you may
   need to ask IT to do this.
2. Run:
   ```
   wsl --install
   ```
3. Open the Microsoft Store and install **Ubuntu 22.04**.
4. When Ubuntu first starts, create a username and password — you will need
   this password for `sudo` commands.
5. You can launch your Ubuntu environment from PowerShell at any time with:
   ```
   wsl -d Ubuntu-22.04
   ```

### 2. Build

From the Ubuntu terminal, follow [Linux](building_from_source.md#linux) in
*Building from source*: install the apt dependencies, clone, build, and put the
tools on your PATH. Nothing about the build is WSL-specific.

### Accessing Windows files from WSL2

Your Windows drives are mounted under `/mnt/`. `C:\MyData\topography` is
`/mnt/c/MyData/topography` from Ubuntu.

Building *in* `/mnt/c/...` works, but file access across that boundary is much
slower than inside the Linux filesystem (`~/`). For the fastest builds, clone
into your Linux home directory and keep only data on the Windows side.

### WSL2 quirks

**No internet from WSL2** (package downloads fail): `/etc/resolv.conf` is
usually wrong. Two files are involved, and the settings are not
interchangeable. First stop WSL2 regenerating the file, by adding to
`/etc/wsl.conf`:

```
[network]
generateResolvConf = false
```

Then set a nameserver in `/etc/resolv.conf`:

```
nameserver 8.8.8.8
```

Save each with `Ctrl+O`, then `Ctrl+X`. Restart WSL2 from a Windows prompt
(`wsl --shutdown`) for the `wsl.conf` change to take effect.

**Executable bits and git:** on a clone under `/mnt/c/`, `chmod +x` is not
visible to git (NTFS has no exec bit). To mark a script executable in the
repository, use `git update-index --chmod=+x <file>`.

---

## Native Windows (QGIS backend only)

This path produces native `.exe` files for the subset of tools the QGIS plugin
drives. It is **not** a development route. If you are working on Windows it is still probably better for you to do a WSL2 window and then run lsdtt3 from Ubuntu. So this section is really just notes for the developers who are compiling the backend for QGIS. 

Limitations, up front:

- Full lsdtt3 functionality exists only on Linux (including WSL2).
- Native Windows is a reduced-feature profile: OpenMP and Eigen are
  off. FFTW3 *is* included via vcpkg, so spectral analysis and Wiener channel
  extraction work.
- These builds exist for the QGIS plugin backend, not for performance.

Packaging the result into a bundle, and installing it in QGIS, is covered in
[QGIS plugin development](qgis_plugin_development.md).

### What you need

- Visual Studio Community 2022 with the C++ build tools
- `Developer PowerShell for VS 2022`
- CMake and Ninja on the path in that shell
- `vcpkg` in your user folder, so dependencies install without admin rights

### 1. Open an x64 developer shell

Launch `Developer PowerShell for VS 2022` from the Start menu.

> [!WARNING]
> Some Visual Studio shell entries start in an **x86** toolchain context. CMake
> then picks a 32-bit compiler and rejects `gdal:x64-windows` as incompatible.
> Always confirm `cl` reports `for x64` before configuring.

If the Start menu shortcut is unavailable, or `cl` reports x86:

```powershell
& "$env:ProgramW6432\Microsoft Visual Studio\2022\Community\Common7\Tools\Launch-VsDevShell.ps1" -Arch amd64 -HostArch amd64
```

`x64 Native Tools Command Prompt for VS 2022` works too. Verify the toolchain:

```powershell
cmake --version
ninja --version
cl        # expect: ... Optimizing Compiler Version ... for x64
```

### 2. Install GDAL and FFTW3 via vcpkg

> [!WARNING]
> This takes a long time — around 2 hours on a standard 2026 laptop. It is a
> one-off: later builds reuse the installed packages.

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
- `tools/package_windows_backend.ps1` expects vcpkg at
  `%USERPROFILE%\tools\vcpkg`. If you put it elsewhere, pass `-VcpkgRoot`.

### 3. Configure and build

From the repository root:

```powershell
cmake -S . -B build\windows-qgis -G "Ninja" `
  -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_TOOLCHAIN_FILE="$env:USERPROFILE\tools\vcpkg\scripts\buildsystems\vcpkg.cmake" `
  -DVCPKG_TARGET_TRIPLET=x64-windows `
  -DBUILD_TESTING=OFF `
  -DLSDTT_OPENMP=OFF `
  -DLSDTT_EIGEN=OFF `
  -DLSDTT_GPU=OFF `
  -DCMAKE_DISABLE_FIND_PACKAGE_PkgConfig=ON
```

Configuring is once per checkout. Then build only the CLIs the QGIS bundle
ships (the `qgis` profile in `tools/cli_profiles.txt`):

```powershell
cmake --build build\windows-qgis --target `
  lsdtt-raster-preprocessing lsdtt-channel-extraction lsdtt-chi-analysis `
  lsdtt-basin-extraction lsdtt-basin-stats lsdtt-surface-metrics `
  lsdtt-channel-network-analysis lsdtt-swath-profile `
  lsdtt-drainage-analysis
```

or everything with `cmake --build build\windows-qgis`.

`-DBUILD_TESTING=OFF` drops test targets the bundle does not need and which
pull in code only exercised on Linux. Disabling PkgConfig does **not** disable
FFTW3; the configure summary should still report `FFTW3       : ON`.

### Why OpenMP must be off on MSVC

`-DLSDTT_OPENMP=OFF` is not optional. Visual Studio's `/openmp` selects OpenMP
**2.0** (2002), which rejects unsigned loop indices, silently ignores
`collapse`, and has no min/max reductions — all of which this codebase uses.
Leaving it on fails with `C3016`, `C4849` and `C7660` part-way through
`lsdtt3_cpu`. CMake warns about this at configure time.

The loss is smaller than it sounds: `C4849` means `collapse` was being ignored
anyway, so MSVC was only ever parallelising outer loops. `/openmp:llvm` would
give real OpenMP 3.1+, but it switches the runtime from `vcomp140.dll` to
`libomp140.x86_64.dll`, which `package_windows_backend.ps1` does not bundle —
so that is a deliberate change, not a flag flip.

To build a **release tag** rather than your working tree, let the packaging
script export, configure and build it in a scratch folder (this configure line,
the `qgis` profile's targets, then packaging):

```powershell
.\tools\package_windows_backend.ps1 -Version 0.7.1 -Tag v0.7.1
```

### 4. Verify

```powershell
Get-ChildItem build\windows-qgis\lsdtt-*.exe
.\build\windows-qgis\lsdtt-raster-preprocessing.exe --version
```

### Troubleshooting

- GDAL configure error mentioning "configuration files were considered but not
  accepted" and `64bit`: your shell is on an x86 toolchain. Re-open an x64
  developer shell and confirm `cl` reports `for x64`.
- `cl : Command line error D8021 : invalid numeric argument '/Wextra'`: update
  to a revision whose `CMakeLists.txt` uses MSVC warning flags (`/W4`) rather
  than GCC-style `-Wextra`.

---

## See also

- [Building from source](building_from_source.md) — the reference Linux build and CMake options
- [QGIS plugin development](qgis_plugin_development.md) — packaging the native build for QGIS
- [Making release packages](making_release_packages.md) — how Windows bundles fit into a release
