# Installing lsdtt3

lsdtt3 runs on Linux. This page describes the three ways to get it running and
points to the detailed instructions for each platform.

---

## Which route is right for you?

| I want to… | Route |
|---|---|
| Try lsdtt3 without compiling anything | [Pre-built install (beta)](#pre-built-install-beta) |
| Run lsdtt3 on a Windows machine | [Windows (WSL2)](#windows-wsl2) |
| Run lsdtt3 on a Mac | [macOS (Multipass)](#macos-multipass) |
| Run lsdtt3 on a Linux machine | [Build from source on Linux](#build-from-source-on-linux) |
| Develop or modify lsdtt3 | [Build from source on Linux](#build-from-source-on-linux) |

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

### CMake options

| Option | Default | Effect |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | Use `Debug` for development or `RelWithDebInfo` for profiling |
| `LSDTT_OPENCV` | `ON` | Build OpenCV-dependent routines (auto-disabled if OpenCV not found) |
| `LSDTT_EIGEN` | `ON` | Build Eigen3-dependent routines (auto-disabled if Eigen3 not found) |
| `LSDTT_GPU` | `OFF` | Enable GPU backend (requires CUDA toolkit) |
