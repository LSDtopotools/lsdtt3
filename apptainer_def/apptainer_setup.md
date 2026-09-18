# Building lsdtt3 with Apptainer

Apptainer (formerly Singularity) runs Linux containers. On Linux it runs
natively. On macOS and Windows it requires a Linux virtual machine first.
This document covers all three platforms.

---

## Linux

Apptainer is available in most distribution package managers.

```bash
# Ubuntu / Debian
sudo apt-get install -y apptainer

# Fedora / RHEL / Rocky
sudo dnf install -y apptainer
```

If your distribution does not package it, follow the official installation
guide at https://apptainer.org/docs/admin/latest/installation.html.

Once installed, skip to [Building the container](#building-the-container).

---

## macOS

Apptainer does not run directly on macOS. The recommended approach is to run
it inside a Linux virtual machine using Lima, which is a lightweight VM manager
for macOS that works on both Intel and Apple Silicon.

### Install Lima

```bash
# Install Homebrew if you do not have it
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install Lima
brew install lima
```

### Create a Lima VM with Apptainer

Lima provides a ready-made Apptainer template:

```bash
limactl start template://apptainer
```

This downloads an Ubuntu VM image and installs Apptainer inside it. It takes
a few minutes on first run.

### Open a shell in the VM

```bash
limactl shell apptainer
```

You are now inside the Linux VM. Your macOS home directory is mounted read-only
at `/Users/<your-username>` and read-write under `~/`. The lsdtt3 source tree
is accessible from inside the VM.

Skip to [Building the container](#building-the-container).

### Stopping and starting the VM

```bash
limactl stop apptainer    # stop the VM
limactl start apptainer   # start it again later
limactl delete apptainer  # remove the VM entirely
```

---

## Windows (WSL2)

Apptainer runs inside WSL2 (Windows Subsystem for Linux 2). If you do not
have WSL2 set up, do that first.

### Enable WSL2

Open PowerShell as Administrator and run:

```powershell
wsl --install
```

Restart your machine when prompted. This installs WSL2 with Ubuntu by default.

If you already have WSL1, upgrade to WSL2:

```powershell
wsl --set-default-version 2
wsl --set-version Ubuntu 2
```

### Open a WSL2 terminal

From the Start menu open **Ubuntu** (or whichever distribution you installed).
All subsequent commands run inside this Linux terminal.

### Install Apptainer inside WSL2

```bash
sudo apt-get update
sudo apt-get install -y apptainer
```

If your Ubuntu version does not have Apptainer in its repositories, install
it from the official release:

```bash
# Download the latest .deb from GitHub releases
APPTAINER_VERSION=1.3.3
wget https://github.com/apptainer/apptainer/releases/download/v${APPTAINER_VERSION}/apptainer_${APPTAINER_VERSION}_amd64.deb
sudo dpkg -i apptainer_${APPTAINER_VERSION}_amd64.deb
```

Check https://github.com/apptainer/apptainer/releases for the current version
number and update `APPTAINER_VERSION` accordingly.

### Access your Windows files from WSL2

Your Windows drives are mounted under `/mnt/`. For example:

```
C:\Users\you\lsdtt3  →  /mnt/c/Users/you/lsdtt3
```

---

## Building the container

These steps are the same on all platforms (Linux, Lima VM on macOS, WSL2 on
Windows). Run them from inside your Linux environment.

Navigate to the directory containing `lsdtt3.def`:

```bash
cd /path/to/lsdtt3
```

Build the container image. This requires either root access or `--fakeroot`
(available on most systems without root):

```bash
# Without root (recommended on shared systems)
apptainer build --fakeroot lsdtt3.sif lsdtt3.def

# With root
sudo apptainer build lsdtt3.sif lsdtt3.def
```

The build downloads Ubuntu 22.04 and installs all dependencies. It takes
5–15 minutes depending on your connection. The output is a single file
`lsdtt3.sif` that contains the complete build environment.

Verify the container built correctly:

```bash
apptainer test lsdtt3.sif
```

This runs the `%test` section of the definition file and reports the versions
of GCC, CMake, GDAL, FFTW3, Eigen3, and NetCDF.

---

## Compiling lsdtt3

### Configure and build

```bash
# From the lsdtt3 source directory
apptainer exec lsdtt3.sif bash -c "
    mkdir -p build &&
    cd build &&
    cmake .. -DCMAKE_BUILD_TYPE=Release &&
    make -j\$(nproc)
"
```

Compiled binaries are written to the build directory, not the source tree:

| What | Location |
|---|---|
| CLI programs (`lsdtt-*`) | `build/` — e.g. `build/lsdtt-raster-preprocessing` |
| C++ test drivers (`lsdtt-test-*`) | `build/tests/` |

To collect the CLIs somewhere else, install them. `cmake --install` copies each
CLI into `<prefix>/bin`:

```bash
apptainer exec lsdtt3.sif cmake --install build --prefix "$HOME/lsdtt3_install"
# CLIs are then in $HOME/lsdtt3_install/bin
```

### Build options

| Option | Default | Effect |
|---|---|---|
| `-DCMAKE_BUILD_TYPE=Release` | Release | Optimised build |
| `-DCMAKE_BUILD_TYPE=Debug` | — | Debug symbols, no optimisation |
| `-DLSDTT_GPU=ON` | OFF | Enable CUDA GPU backend (requires CUDA in container) |

### Parallel compilation

`$(nproc)` uses all available CPU cores. On a machine with many cores you may
want to limit this to avoid running out of memory:

```bash
make -j4
```

---

## Running lsdtt3 programs

All CLI programs are in `build/` after compilation (or `<prefix>/bin` if you
ran `cmake --install`, above). Run them through the container so they have
access to the correct GDAL and library versions:

```bash
apptainer exec lsdtt3.sif ./build/lsdtt-raster-preprocessing my_params.driver
apptainer exec lsdtt3.sif ./build/lsdtt-channel-extraction my_params.driver
```

Or open an interactive shell in the container:

```bash
apptainer shell lsdtt3.sif
```

From inside the shell you can run programs directly without the
`apptainer exec` prefix.

### Binding data directories

By default Apptainer makes your home directory and the current directory
available inside the container. If your data lives elsewhere, bind it
explicitly:

```bash
apptainer exec --bind /data/dems:/data/dems lsdtt3.sif \
    ./build/lsdtt-raster-preprocessing my_params.driver
```

---

## Rebuilding after code changes

The container only needs to be rebuilt if the dependencies change (i.e. the
`.def` file changes). For code changes, just recompile:

```bash
apptainer exec lsdtt3.sif bash -c "cd build && make -j\$(nproc)"
```

Or to do a clean rebuild:

```bash
apptainer exec lsdtt3.sif bash -c "
    rm -rf build &&
    mkdir build &&
    cd build &&
    cmake .. -DCMAKE_BUILD_TYPE=Release &&
    make -j\$(nproc)
"
```

---

## Troubleshooting

**`--fakeroot` not available**
Some older systems or restricted HPC environments disable fakeroot. Try
building with `sudo` instead, or ask your system administrator to enable
user namespaces (`kernel.unprivileged_userns_clone=1`).

**GDAL version mismatch at runtime**
If you get GDAL errors when running on a system that also has GDAL installed
natively, the system GDAL may be interfering. Force the container's libraries:

```bash
apptainer exec --no-home lsdtt3.sif ./build/lsdtt-raster-preprocessing my_params.driver
```

**Lima VM runs out of disk space**
The default Lima VM has a small disk. Increase it by editing the Lima config
before starting:

```bash
limactl edit apptainer
# increase disk field, e.g. disk: "50GiB"
```

**WSL2 slow file I/O on Windows drives**
File access across the Windows/WSL2 boundary (`/mnt/c/...`) is slow. For best
performance, keep the lsdtt3 source and data files inside the WSL2 filesystem
(`~/lsdtt3`) rather than on a Windows drive.
