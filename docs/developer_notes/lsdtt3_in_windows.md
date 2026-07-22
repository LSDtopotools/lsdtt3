# Running lsdtt3 on Windows

lsdtt3 is a Linux-native codebase. On Windows you run it inside **Windows Subsystem for Linux (WSL2)**, which gives you a full Ubuntu environment without leaving your Windows machine.

## 1. Install WSL2 and Ubuntu

1. Open PowerShell with administrator privileges. On a managed machine you may need to ask IT to do this.
2. Run:
   ```
   wsl --install
   ```
3. Open the Microsoft Store and install **Ubuntu 22.04**.
4. When Ubuntu first starts, create a username and password — you will need this password for `sudo` commands.
5. You can launch your Ubuntu environment from PowerShell at any time with:
   ```
   wsl -d Ubuntu-22.04
   ```

### WSL2 internet troubleshooting

If WSL2 cannot reach the internet (package downloads fail), the `/etc/resolv.conf` file is likely wrong. Fix it:

```bash
sudo nano /etc/resolv.conf
```

Make the file contain:

```
generateResolvConf = false
nameserver 8.8.8.8
```

Save with `Ctrl+O`, then `Ctrl+X`.

---

## 2. Install required libraries

From inside your Ubuntu terminal, install the build tools and libraries that lsdtt3 needs:

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

## 3. Set up an SSH key and clone the repository

GitHub no longer accepts password authentication for `git clone`. The recommended approach is an SSH key.

### Generate a key

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
3. Paste the key and give it a name (e.g. `wsl2-ubuntu`).
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

This creates a directory called `lsdtt3` in your current working directory.

---

## 4. Build the code

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

After a successful build, the executables are in `build/` (or `build/bin/` if you ran `cmake --install`). To make them available from any directory, add the build folder to your PATH:

```bash
export PATH=/path/to/lsdtt3/build:$PATH
```

To make this permanent, add that line to your `~/.bashrc`:

```bash
echo 'export PATH=/path/to/lsdtt3/build:$PATH' >> ~/.bashrc
source ~/.bashrc
```

Replace `/path/to/lsdtt3` with the actual path where you cloned the repository.

---

## 5. Accessing Windows files from WSL2

Your Windows drives are mounted under `/mnt/`. For example, your `C:\` drive is at `/mnt/c/`. If your data lives in `C:\MyData\topography`, you can reach it from Ubuntu at `/mnt/c/MyData/topography`.
