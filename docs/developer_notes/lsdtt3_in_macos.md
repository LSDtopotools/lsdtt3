# Running lsdtt3 on macOS

lsdtt3 is a Linux-native codebase. On macOS the simplest path is to run it inside an Ubuntu virtual machine managed by **Multipass**, a lightweight tool made by Canonical (the Ubuntu developers) that spins up Ubuntu instances with a single command.

## 1. Install Homebrew and Multipass

If you do not already have Homebrew (the standard macOS package manager), install it first by following the instructions at [brew.sh](https://brew.sh). Then install Multipass:

```bash
brew install --cask multipass
```

Alternatively, download the `.pkg` installer directly from the Multipass website.

## 2. Launch an Ubuntu instance

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

---

## 3. Install required libraries

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

---

## 4. Set up an SSH key and clone the repository

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

---

## 5. Build the code

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

---

## 6. Sharing files between macOS and the Ubuntu VM

Multipass mounts your macOS home directory inside the VM at `/home/ubuntu/Home` by default. You can also mount any macOS directory explicitly:

```bash
multipass mount /Users/yourname/lsdtt_data lsdtt:/lsdtt_data
```

This makes `/Users/yourname/lsdtt_data` on your Mac visible as `/lsdtt_data` inside the `lsdtt` instance. Run this command from your macOS terminal (not from inside the VM shell). To unmount:

```bash
multipass umount lsdtt:/lsdtt_data
```

This lets you keep your DEM data on your Mac while running lsdtt3 tools inside the Ubuntu VM.
