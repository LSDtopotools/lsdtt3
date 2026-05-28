# Beta testing and distribution notes

This document covers three things:

1. **For beta testers without GitHub access** — how to install lsdtt3 on a Windows machine using WSL2 from a tarball shared directly.
2. **For beta testers with GitHub access** — how to authenticate with GitHub and download the release tarball themselves.
3. **For the developer** — how the release pipeline works and how to create a release.

A fourth section covers the future conda-forge route for when the repository goes public.

---

## For beta testers: installing on Windows via WSL2

lsdtt3 runs on Linux. On Windows the easiest way to get a Linux environment is
Windows Subsystem for Linux (WSL2), which is built into all modern Windows
installations.

### Step 1 — Install WSL2 and Ubuntu

1. Open PowerShell **as administrator** (right-click the Start menu → Windows
   PowerShell (Admin) or Terminal (Admin)).
2. Run:
   ```
   wsl --install
   ```
3. When prompted, restart your computer.
4. After restarting, open the Microsoft Store, search for **Ubuntu 22.04 LTS**,
   and install it.
5. Launch Ubuntu from the Start menu. The first time you open it you will be
   asked to create a username and password — choose something you will remember.
   This password is used whenever you run `sudo`.

#### WSL2 internet troubleshooting

If Ubuntu cannot reach the internet (package downloads fail), edit the DNS
configuration:

```bash
sudo nano /etc/resolv.conf
```

Make the file contain exactly:
```
nameserver 8.8.8.8
```

Save with Ctrl-O, then Ctrl-X. If the file reverts after a reboot, also create
`/etc/wsl.conf` with:
```
[network]
generateResolvConf = false
```

### Step 2 — Install lsdtt3

**If you have been given a `.tar.gz` file directly** (the repo is private so the
download link requires a GitHub account):

```bash
bash install.sh /path/to/lsdtt3-ubuntu-22.04.tar.gz
```

Replace `/path/to/` with the actual location. If the file is on your Windows
desktop, the path in WSL2 is:
```bash
bash install.sh /mnt/c/Users/YourWindowsUsername/Desktop/lsdtt3-ubuntu-22.04.tar.gz
```

**If the repo is public and you have a GitHub account with access:**

```bash
curl -fsSL https://raw.githubusercontent.com/LSDtopotools/lsdtt3/main/install.sh | bash
```

The installer:
- installs the required system libraries (GDAL, FFTW, OpenCV) via `apt`
- unpacks the pre-built executables to `~/.local/lsdtt3/bin`
- adds that directory to your `PATH` in `~/.bashrc`

### Step 3 — Verify the installation

After the installer finishes, either open a new terminal or run:

```bash
source ~/.bashrc
```

Then check that the tools are visible:

```bash
lsdtt-channel-extraction --help
lsdtt-basin-extraction --help
```

Each command should print a list of parameters. If you see `command not found`,
check that `~/.local/lsdtt3/bin` is in your `PATH`:

```bash
echo $PATH | tr ':' '\n' | grep lsdtt3
```

If nothing is printed, run `source ~/.bashrc` again and retry.

### Step 4 — Accessing your Windows files from WSL2

Your Windows drives are mounted under `/mnt/`. Your `C:` drive is `/mnt/c/`.
So a file at `C:\Users\Alice\data\dem.tif` is visible in WSL2 as:

```
/mnt/c/Users/Alice/data/dem.tif
```

You can run lsdtt3 tools directly on files in your Windows file system:

```bash
lsdtt-channel-extraction --dem /mnt/c/Users/Alice/data/dem.tif --output_directory /mnt/c/Users/Alice/output/
```

### Reporting bugs

When filing a bug report please include:

1. The exact command you ran.
2. The full terminal output (copy and paste it).
3. The Ubuntu version: `cat /etc/os-release | grep VERSION=`
4. The lsdtt3 version: check the release tag you installed, e.g. `v0.1.0`.

Send this to the developer or open an issue on the GitHub repository if you have
access.

---

## For beta testers with GitHub access: downloading the tarball yourself

If you have been added to the GitHub repository you can download the release
tarball directly using the `gh` command-line tool rather than waiting for the
developer to share the file.

### Step 1 — Install the GitHub CLI

Inside your WSL2 Ubuntu terminal:

```bash
sudo apt install gh
```

If the apt version is outdated, install the current release from GitHub's own
package feed:

```bash
curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg \
    | sudo dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] \
    https://cli.github.com/packages stable main" \
    | sudo tee /etc/apt/sources.list.d/github-cli.list
sudo apt update && sudo apt install gh
```

### Step 2 — Authenticate with GitHub

```bash
gh auth login
```

When prompted:

1. **Where do you use GitHub?** → `GitHub.com`
2. **What is your preferred protocol?** → `SSH`
3. **Upload your public SSH key?** → select your existing key (e.g. `~/.ssh/id_rsa.pub`), or press Enter to skip if you only want HTTPS token access
4. **How would you like to authenticate?** → `Paste an authentication token`

To generate the token, go to **GitHub.com → profile picture → Settings →
Developer settings** (bottom of the left sidebar) **→ Personal access tokens →
Tokens (classic) → Generate new token (classic)**. Give it any name, set an
expiry, and check exactly two scopes:

- **`repo`** (full private repo access)
- **`read:org`** (required to verify organisation membership)

Click **Generate token**, copy the `ghp_...` string, and paste it when `gh auth
login` asks.

### Step 3 — Download the release tarball and install

```bash
gh release download v0.1.0 --repo LSDtopotools/lsdtt3 \
    --pattern "lsdtt3-ubuntu-22.04.tar.gz" -D /tmp/
bash install.sh /tmp/lsdtt3-ubuntu-22.04.tar.gz
```

Replace `v0.1.0` with the tag the developer told you to use. Then follow
Steps 3 and 4 from the section above (verify the installation and access your
Windows files).

---

## For the developer: creating a release

### How the build pipeline works

The file `.github/workflows/release.yml` defines a GitHub Actions workflow that
runs automatically whenever a version tag is pushed to the repository. GitHub
Actions works on private repositories.

The workflow:

1. Starts an Ubuntu 22.04 virtual machine.
2. Installs the same apt packages that beta testers will have (`libgdal-dev`,
   `libfftw3-dev`, `libopencv-dev`, `libeigen3-dev`).
3. Runs `cmake` and builds all `lsdtt-*` executables in Release mode.
4. Packs them into `lsdtt3-ubuntu-22.04.tar.gz`.
5. Attaches the tarball to a GitHub Release named after the tag.

### Creating a release

1. Make sure `main` is in the state you want to release and all tests pass.
2. Choose a version number following semantic versioning, e.g. `v0.2.1`.
3. Tag and push:
   ```bash
   git tag v0.2.1
   git push origin v0.2.1
   ```
4. Go to the repository on GitHub → **Actions**. You will see the
   "Build and release" workflow running. It takes roughly 5–10 minutes.
5. When it finishes, go to **Releases**. A new release named `v0.2.1` will have
   been created automatically with `lsdtt3-ubuntu-22.04.tar.gz` attached.
6. Download that file and share it with beta testers (email, shared drive, etc.).

### Running the build manually (without a tag)

You can trigger the workflow at any time from the GitHub Actions tab:

1. Go to **Actions** → **Build and release** → **Run workflow**.
2. Leave the tag input blank and click **Run workflow**.
3. After it completes the tarball is available as a workflow **artifact**
   (visible in the run summary, downloadable for 14 days).

This is useful for producing a test build to share with a beta tester without
creating an official release.

### Updating the install script after a dependency change

If a new optional dependency is added to CMakeLists.txt, update two places:

- **`.github/workflows/release.yml`** — add the `-dev` package to the
  `Install dependencies` step so the CI build picks it up.
- **`install.sh`** — add the corresponding runtime package to `RUNTIME_PACKAGES`
  so beta testers have it installed before running the binary.

---

## Future: conda-forge (when the repository goes public)

[conda-forge](https://conda-forge.org) is a community-maintained collection of
conda packages. Once lsdtt3 is public, publishing there means users can install
with a single command and all dependency management is handled automatically:

```bash
conda install -c conda-forge lsdtt3
```

This is also the right route for QGIS integration, because QGIS and lsdtt3
would share the same conda-managed GDAL installation.

### The recipe

The file `conda-recipe/meta.yaml` in this repository is the draft recipe. It is
not yet submitted to conda-forge because the repository is private.

### Submitting to conda-forge

When the repository goes public:

1. Update the `source.sha256` field in `conda-recipe/meta.yaml`:
   ```bash
   curl -fsSL https://github.com/LSDtopotools/lsdtt3/archive/refs/tags/vX.Y.Z.tar.gz \
       | sha256sum
   ```
2. Fork [conda-forge/staged-recipes](https://github.com/conda-forge/staged-recipes).
3. Copy `conda-recipe/meta.yaml` to `recipes/lsdtt3/meta.yaml` in your fork.
4. Open a pull request. The conda-forge bot will build the package on Linux,
   macOS, and (if not skipped) Windows and leave review comments.
5. After the PR is merged, the package is live. Update the project README to
   show the install command.

### Testing the recipe locally before submitting

```bash
conda install conda-build
conda build conda-recipe/
```

This builds the package on your local machine using your local source tree. The
result is a `.tar.bz2` file you can install with `conda install --use-local lsdtt3`.

### Personal Anaconda.org channel (interim, before conda-forge)

While the repo is still private it is possible to publish packages to a personal
channel on [anaconda.org](https://anaconda.org) without making the source public.
You build the package locally (or in CI) and push the compiled artifact:

```bash
conda install anaconda-client
conda build conda-recipe/
anaconda upload ~/miniconda3/pkgs/lsdtt3-0.2.0-linux-64.tar.bz2
```

Users then install via:
```bash
conda install -c your_anaconda_username lsdtt3
```

This requires building separately on each platform you want to support. It is a
reasonable interim step if beta testers use conda environments rather than bare
Ubuntu.
