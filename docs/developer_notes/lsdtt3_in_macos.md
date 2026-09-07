# Running lsdtt3 on macOS

lsdtt3 can run on macOS in two ways:

1. **Native build** (recommended for developers): compile directly on macOS using Homebrew
   and Xcode Command Line Tools. The resulting tools run natively at full speed.
2. **Inside a virtual machine** (alternative): run it in an Ubuntu 22.04 VM managed by
   **Multipass**, a lightweight tool made by Canonical. This is the safer option for
   end users, since lsdtt3 has always been tested primarily on Linux.

Choose native build if you are developing or testing the QGIS plugin backend. Choose
Multipass if you want to avoid platform-specific issues or prefer a sandboxed setup.

---

## Native macOS Build (Recommended for Developers)

### Prerequisites

- **Apple Silicon or Intel Mac** running macOS 11 (Big Sur) or later
- **Xcode Command Line Tools** (includes Apple Clang, Git)
- **Homebrew** (the macOS package manager)
- **CMake** 3.15 or later
- **Build dependencies**: GDAL, FFTW3, optional libraries

The build defaults to your machine's native architecture (`uname -m` reports `arm64` on
Apple Silicon, `x86_64` on Intel). Cross-compilation to a different architecture is
possible but not covered here.

### Install Xcode Command Line Tools

Open Terminal and run:

```bash
xcode-select --install
```

A dialog will appear. Click "Install" and wait for it to complete (20–30 minutes on a
fresh Mac). Verify the installation:

```bash
clang --version
git --version
```

### Install Homebrew and build dependencies

From [brew.sh](https://brew.sh), install Homebrew if you don't have it:

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Then install the required and optional lsdtt3 dependencies:

```bash
brew install cmake gdal fftw libomp netcdf opencv eigen
```

What each package provides:

| Package | Role | Required? |
|---|---|---|
| `cmake` (≥ 3.15) | Build system | **yes** |
| `gdal` | Raster and vector I/O, CRS handling | **yes** |
| `libomp` | OpenMP parallelism (Apple Clang's default OpenMP is outdated) | **yes** |
| `fftw` | Spectral analysis, wavelet polyfit | no |
| `netcdf` | ERA40 / palaeomagnetic corrections in cosmo tools | no |
| `opencv` | Image analysis in surface metrics | no |
| `eigen` | k_sn / u\* inversion in chi analysis | no |

Notes:

- **`libomp` must be installed via Homebrew.** The reason is not that Apple's OpenMP is
  an old *version* — it is that Apple Clang ships **no OpenMP runtime or headers at
  all**. `-fopenmp` is not usable out of the box, so `find_package(OpenMP)` fails and the
  build falls back to single-threaded (CMakeLists.txt prints
  `OpenMP not found - CPU operations will be single-threaded`). Homebrew's `libomp`
  supplies the missing runtime.

  **Not verified on macOS:** whether `find_package(OpenMP)` locates Homebrew's libomp
  without help. It very likely does **not** — libomp is keg-only, so Homebrew does not
  link it into the default search paths, and CMakeLists.txt currently has no `APPLE`
  branch to point at it (see "Known gaps" below). Expect to pass the paths explicitly, or
  to configure with `-DLSDTT_OPENMP=OFF` for a working single-threaded build while this is
  resolved.
- **Architecture-specific builds**: Homebrew bottles are pre-compiled for your machine's
  native architecture. If you need to build for a different architecture (e.g. Intel on
  Apple Silicon), use `--build-from-source` or see the Homebrew documentation for
  cross-compilation. This is uncommon; skip it unless you have a specific reason.

### Clone the repository

Use SSH (recommended) or HTTPS:

```bash
git clone git@github.com:lsdtopotools/lsdtt3.git
cd lsdtt3
```

### Configure and build

```bash
mkdir build/macos-core
cd build/macos-core
cmake -DCMAKE_BUILD_TYPE=Release ../..
cmake --build . -j$(sysctl -n hw.ncpu)
```

The `-j$(sysctl -n hw.ncpu)` flag parallelizes the build across all CPU cores. On Apple
Silicon or modern Intel Macs, the full build takes 5–15 minutes depending on which
optional libraries are present.

CMake will print a summary of which optional libraries it found:

```
=== lsdtt3 build configuration ===
  Build type  : Release
  C++ std     : 17
  GDAL        : 3.x
  OpenMP      : TRUE
  FFTW3       : TRUE
  OpenCV      : TRUE
  GPU (CUDA)  : FALSE
  Eigen3      : TRUE
===================================
```

Read this summary — it tells you which tools were built. A missing library means the
corresponding tool is skipped, and the build still succeeds.

### Put the tools on your PATH

After a successful build, the executables are in `build/macos-core/`. To use them from
any directory:

```bash
export PATH=/path/to/lsdtt3/build/macos-core:$PATH
```

Test one:

```bash
lsdtt-raster-preprocessing --version
```

To make this permanent, add the export line to your shell profile:

```bash
# For zsh (default on macOS 10.15+)
echo 'export PATH=/path/to/lsdtt3/build/macos-core:$PATH' >> ~/.zshrc
source ~/.zshrc

# Or for bash (if you use it)
echo 'export PATH=/path/to/lsdtt3/build/macos-core:$PATH' >> ~/.bashrc
source ~/.bashrc
```

### Run the test suite (optional)

```bash
cd /path/to/lsdtt3/build/macos-core
ctest --output-on-failure
```

See `tests/README.md` in the repository for what each test covers.

---

## Bundling for the QGIS Plugin (Native macOS Build)

If you are building the QGIS plugin backend, you must package the tools and libraries
into a self-contained bundle. The script `tools/package_macos_backend.sh` automates this.

### Prerequisites

- A complete native macOS build (see above)
- `otool`, `install_name_tool`, `codesign` (included with Xcode Command Line Tools)
- `gdal-config` (from Homebrew's gdal package)
- `python3` (included in Homebrew or as part of Xcode CLT)

### Create the bundle

From the repository root:

```bash
tools/package_macos_backend.sh --version 0.4.0
```

This creates:

- Staging directory: `dist/lsdtt3-backend-macos-<arch>-core-v0.4.0/` (where `<arch>` is
  `arm64` on Apple Silicon or `x86_64` on Intel)
- Tarball: `dist/lsdtt3-backend-macos-<arch>-core-v0.4.0.tar.gz`

The bundle contains:

```
bin/                              # lsdtt-* executables
lib/                              # bundled shared libraries (libc++, libomp, GDAL, etc.)
share/gdal, share/proj            # GDAL/PROJ data directories
setenv.sh                         # environment setup script
README.md                         # end-user instructions
manifest.json                     # checksums and version info
THIRD_PARTY_NOTICES.txt          # library attribution
LICENSE                          # license file
```

### Verify the bundle is self-contained (macOS-specific)

The script runs a self-containment check at the end and **aborts if any dependency is
missing**. If you see an error like:

```
ERROR: <binary> needs '@loader_path/...' but it is not staged in lib/
ERROR: bundle is not self-contained; aborting
```

Do not ignore it. This check catches the class of errors that broke v0.2.0 (shipped a
bundle missing `libgdal.so.30`, then silently failed for users).

**If the self-containment check passes**, the bundle is safe to distribute. Verify
manually if needed:

```bash
cd dist/lsdtt3-backend-macos-<arch>-core-v0.4.0
otool -L bin/lsdtt-raster-preprocessing
```

Every non-system library should show an `@rpath/` reference:

```
@rpath/libc++.1.dylib
@rpath/libomp.dylib
@rpath/libgdal.dylib
...
```

(System libraries in `/usr/lib`, `/System/Library` are not shown by `otool -L` — they
live in the dyld shared cache.)

### Test the bundle before distribution

```bash
cd dist/lsdtt3-backend-macos-<arch>-core-v0.4.0
source setenv.sh
lsdtt-raster-preprocessing --help
```

If `--help` works, the bundle is ready. If you see a `dyld` error about a missing
library, the self-containment check may have missed something — report it.

### Known macOS issues

**Downloaded bundles and Gatekeeper quarantine:** When you download the bundle via a
web browser, macOS marks it with the `com.apple.quarantine` attribute. Older macOS
versions (pre-Monterey) block execution of unsigned binaries with this attribute. The
bundle README includes instructions for users, but during development you can clear it:

```bash
# In the distribution directory
xattr -dr com.apple.quarantine lsdtt3-backend-macos-<arch>-core-v0.4.0
```

After that, the bundle runs freely in the current Terminal session.

### For QGIS plugin integration

Point the QGIS plugin to the unpacked bundle folder (containing `bin/`, `lib/`, `share/`,
`manifest.json`). See `qgis_plugin/README.md` for detailed plugin setup.

---

## Alternative: Multipass Virtual Machine (No Native Build)

If you prefer not to build natively on macOS, or if you encounter platform-specific issues,
you can build and run lsdtt3 inside an Ubuntu virtual machine. This avoids any macOS-specific
complications and matches the environment used in CI.

### 1. Install Homebrew and Multipass

If you do not already have Homebrew (the standard macOS package manager), install it by
following the instructions at [brew.sh](https://brew.sh). Then install Multipass:

```bash
brew install --cask multipass
```

Alternatively, download the `.pkg` installer directly from the [Multipass website](https://multipass.run).

### 2. Launch an Ubuntu instance

```bash
multipass launch --name lsdtt --cpus 4 --memory 8G --disk 40G 22.04
```

Adjust `--cpus`, `--memory`, and `--disk` to match your machine. 8 GB of RAM and 40 GB of disk are comfortable for building and running lsdtt3. The `22.04` argument requests Ubuntu 22.04 LTS.

Once the instance is running, open a shell inside it:

```bash
multipass shell lsdtt
```

You are now inside a Ubuntu 22.04 environment. All subsequent steps run here unless stated otherwise.

To leave the Ubuntu shell at any time, type `exit`. To stop and restart the instance later:

```bash
multipass stop lsdtt
multipass start lsdtt
```

### 3. Install required libraries

From inside the Ubuntu shell, install the build tools and libraries that lsdtt3 needs:

```bash
sudo apt update
sudo apt install -y \
    build-essential \
    cmake \
    pkg-config \
    libgdal-dev \
    libfftw3-dev \
    libopencv-dev \
    libeigen3-dev \
    libnetcdf-dev
```

What each package provides:

| Package | Role |
|---|---|
| `build-essential` | GCC C++ compiler and standard build tools |
| `cmake` | Build system generator (version 3.15+ required) |
| `pkg-config` | Lets CMake locate libraries like FFTW3 |
| `libgdal-dev` | **Required** — raster I/O and CRS handling |
| `libfftw3-dev` | Optional — enables spectral analysis tools |
| `libopencv-dev` | Optional — enables image-analysis routines |
| `libeigen3-dev` | Optional — enables chi/ksn inversion tools |
| `libnetcdf-dev` | Optional — enables ERA40 atmosphere correction in cosmogenic tools |

### 4. Set up an SSH key and clone the repository

GitHub no longer accepts password authentication for `git clone`. The recommended approach is an SSH key.

### Generate a key

Run this inside the Ubuntu shell:

```bash
ssh-keygen -t ed25519 -C "your_email@example.com"
```

Press Enter to accept the default file location (`~/.ssh/id_ed25519`). Add a passphrase if you want extra security.

### Add the key to the SSH agent

```bash
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
```

### Add the public key to GitHub

Print your public key:

```bash
cat ~/.ssh/id_ed25519.pub
```

Copy the entire output. Then in GitHub:

1. Go to **Settings → SSH and GPG keys**.
2. Click **New SSH key**.
3. Paste the key and give it a name (e.g. `multipass-ubuntu`).
4. Click **Add SSH key**.

### Verify the connection

```bash
ssh -T git@github.com
```

You should see a message like `Hi <username>! You've successfully authenticated`.

### Clone lsdtt3

```bash
git clone git@github.com:lsdtopotools/lsdtt3.git
```

This creates a directory called `lsdtt3` in your current working directory inside the VM.

### 5. Build the code

lsdtt3 uses CMake. The standard out-of-source build pattern is:

```bash
cd lsdtt3
mkdir build
cd build
cmake ..
make -j$(nproc)
```

CMake will print a summary at the end showing which optional libraries it found:

```
=== lsdtt3 build configuration ===
  Build type  : Release
  C++ std     : 17
  GDAL        : <version>
  OpenMP      : TRUE
  FFTW3       : TRUE
  OpenCV      : TRUE
  GPU (CUDA)  : FALSE
  Eigen3      : TRUE
===================================
```

Tools that depend on a missing optional library are simply skipped — the build does not fail.

### Add the binaries to your PATH

After a successful build, the executables are in `build/`. To make them available from any directory:

```bash
export PATH=/path/to/lsdtt3/build:$PATH
```

To make this permanent, add that line to your `~/.bashrc`:

```bash
echo 'export PATH=/path/to/lsdtt3/build:$PATH' >> ~/.bashrc
source ~/.bashrc
```

Replace `/path/to/lsdtt3` with the actual path where you cloned the repository inside the VM.

### 6. Sharing files between macOS and the Ubuntu VM

Multipass mounts your macOS home directory inside the VM at `/home/ubuntu/Home` by default. You can also mount any macOS directory explicitly:

```bash
multipass mount /Users/yourname/lsdtt_data lsdtt:/lsdtt_data
```

This makes `/Users/yourname/lsdtt_data` on your Mac visible as `/lsdtt_data` inside the `lsdtt` instance. Run this command from your macOS terminal (not from inside the VM shell). To unmount:

```bash
multipass umount lsdtt:/lsdtt_data
```

This lets you keep your DEM data on your Mac while running lsdtt3 tools inside the Ubuntu VM.
