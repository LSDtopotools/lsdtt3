# Installing lsdtt3

The main option for installation of `lsdtt3` is the **pre-built release package** which 
installs a series of command line interface programs: you run these in a Linux terminal (or a terminal on your operating system).
This method is the way the developers primarily use this software, and is built for automated, repeatable analysis of large topographic datasets. 

However, you can also use `lsdtt3` in [Google Colab](#google-colab), in our [QGIS plugin](#qgis-plugin), and using our [Python package](#python-package)

> **What lsdtt3 is made of.** The software is its **command-line tools**
> (`lsdtt-*`, 46 of them on Linux), each driven by a parameter file. Everything
> lsdtt3 can do, it does through these. The **Python bindings** (`import
> lsdtt3`) expose part of the same C++ library for use in scripts and
> notebooks; they do not cover everything the tools do. `pip install lsdtt3`
> installs both: the bindings and every command-line tool.

**`pip install lsdtt3` is a second way to get the command-line tools.** The
Python package installs every `lsdtt-*` tool as a command, the same set as the
Linux release package, on Linux, macOS and Windows. If you already work in a
Python environment (conda, pixi or a virtual environment) this is often the
quickest route, and on macOS and Windows it is the only native way to get the
full tool set. See [Release package or pip?](#release-package-or-pip)

For the **pre-built release package** you do not
need a compiler, you do not need to install GDAL or any other library, and on
Linux you do not need `sudo` — the package carries the command-line tools and
every library they depend on.

> **Want the very latest code instead of a release?**
> Then you need to compile it yourself. If you are not a developer you really do not want to do this. That activity is
> documented separately in
> [Building lsdtt3 from source](../developer_notes/building_from_source.md).

---

## Where the packages live

**Download everything from [GitHub Releases](https://github.com/LSDtopotools/lsdtt3/releases/latest).**
Every package is built and attached automatically when a release is made, and
no account is needed. These links always fetch the newest release:

| Package | Permanent link to the latest version |
|---|---|
| Linux | [`lsdtt3-backend-linux-x86_64-core.tar.gz`](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-linux-x86_64-core.tar.gz) |
| Windows (QGIS tool set) | [`lsdtt3-backend-win64-qgis.zip`](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-win64-qgis.zip) |
| macOS, Apple Silicon (QGIS tool set) | [`lsdtt3-backend-macos-arm64-qgis.tar.gz`](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-macos-arm64-qgis.tar.gz) |
| macOS, Intel (QGIS tool set) | [`lsdtt3-backend-macos-x86_64-qgis.tar.gz`](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-macos-x86_64-qgis.tar.gz) |
| QGIS plugin | [`lsdtt3-qgis-plugin.zip`](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-qgis-plugin.zip) |
| Python wheels | see [Python package](#python-package) |

These unversioned names exist from release 0.9.2 on. Each release also carries
the same files with the version in the name, which is what to download for a
particular older release (open it from the
[list of releases](https://github.com/LSDtopotools/lsdtt3/releases)). The
folder inside each archive is always named with its version.

> A copy of some packages is also kept at
> <https://www.geos.ed.ac.uk/~smudd/lsdtt_packages/>. Use it **only as a last
> resort**, if GitHub is unreachable from your network: it is updated by hand
> and is usually several releases out of date.

Versioned file names look like this:

```
lsdtt3-backend-linux-x86_64-core-v0.7.1.tar.gz     Linux
lsdtt3-backend-win64-qgis-v0.7.1.zip               Windows (native, limited to the QGIS tool set)
lsdtt3-backend-macos-arm64-qgis-v0.7.1.tar.gz      macOS, Apple Silicon (QGIS tool set; from 0.9.1)
lsdtt3-backend-macos-x86_64-qgis-v0.7.1.tar.gz     macOS, Intel (QGIS tool set; from 0.9.1)
lsdtt3-qgis-plugin-0.7.1.zip                       QGIS plugin (shares the number from 0.7.0)
```

`core` and `qgis` are **profiles** — which set of command-line tools the
package contains. `core` is the general-purpose set; `qgis` is the smaller
subset the QGIS plugin drives. Native Windows and macOS packages are `qgis`
only: for the full tool set on those systems, `pip install lsdtt3` (see
[Python package](#python-package)), or run the Linux package under
[WSL2](#windows-via-wsl2) or in a [Linux VM](#macos).

> lsdtt3 has **one version number** for everything: the command-line tools,
> the backend packages, the QGIS plugin and the Python package. Files from the
> same release go together.

> **You do not need to know the version number to install.** The links and
> commands on this page fetch the latest release and match the versioned
> folder inside it with `*`. A version number that appears below is an example
> of a versioned name. To see which version you have, look at the folder name,
> or run any tool with `-v`.

---

## Choose your route

**1. Where do you want to run lsdtt3?**

| | |
|---|---|
| On my own computer | continue to question 2 |
| In a browser, with nothing installed | [Google Colab](#google-colab) |
| Inside QGIS | [QGIS plugin](#qgis-plugin) (plus a backend package from question 2) |
| From Python scripts | [Python package](#python-package) |
| From a terminal, and I already use Python environments | [Python package](#python-package): it installs all the command-line tools too |

**2. What is your computer running?**

| | |
|---|---|
| Linux | [Linux package](#linux-package) |
| Windows | [Windows, via WSL2](#windows-via-wsl2) — **recommended**, full functionality |
| Windows, and I want native `.exe` files | [Python package](#python-package) (`pip install lsdtt3`) for the full tool set, or the [Windows native package](#windows-native-package) — reduced tool set, the QGIS plugin's backend |
| macOS | [Python package](#python-package) (`pip install lsdtt3`) for the full tool set natively, or [macOS](#macos) for the native QGIS package and the Linux VM |

### Release package or pip?

Both give you the same command-line tools, built from the same release, each
carrying its own GDAL and PROJ. They differ in what they need and in what they
touch:

| | Release package | `pip install lsdtt3` |
|---|---|---|
| Needs | nothing: download and unpack | Python 3.12 or later, in a conda, pixi or virtual environment |
| Tools on Linux | all 46 (`core`) | the same 46 |
| Tools on macOS and Windows | the 9 the QGIS plugin uses | all 46 |
| Oldest Linux | Ubuntu 22.04 (glibc 2.35) | glibc 2.28 (e.g. RHEL 8, Ubuntu 20.04) |
| macOS | 15 or later, Apple Silicon or Intel | 15 or later, Apple Silicon or Intel |
| Its own GDAL, PROJ and their data | yes | yes |
| Affects other software | sourcing `setenv.sh` puts its `lib/`, `GDAL_DATA` and `PROJ_DATA` into that shell, so other GDAL programs started from it may pick them up | no: only the `lsdtt-*` commands see the wheel's libraries and data |
| Python bindings | no | yes |
| Updating | download the new release | `python -m pip install -U lsdtt3` |
| Removing | delete the folder | `python -m pip uninstall lsdtt3` |

**GDAL and PROJ conflicts.** Neither route depends on a GDAL or PROJ you have
installed, and both make the tools use their own `proj.db` even when
`PROJ_DATA` or `PROJ_LIB` points elsewhere (check with any tool's `-v`, which
should say `bundled:`; see
[If you already have GDAL or PROJ installed](#if-you-already-have-gdal-or-proj-installed)).
The wheel is, if anything, the better isolated of the two: its libraries carry
renamed file names, and it changes nothing in your shell. The release package
is the better choice when there is no Python to install into, on a shared
server or in a container, or when you want one folder you can copy between
machines.

### Having both installed

You can have a release package and a pip install on the same machine, and
their GDAL, PROJ and data stay separate: each tool uses the libraries and
`proj.db` that came with it. What is not separate is **which copy runs**. Both
put `lsdtt-*` commands on your `PATH`, and the shell runs whichever it finds
first:

| Order in this shell | `lsdtt-*` runs from |
|---|---|
| `source setenv.sh`, then activate the environment (`conda activate`, `pixi shell`, `source venv/bin/activate`) | the pip install |
| activate the environment, then `source setenv.sh` | the release package |
| `pixi run lsdtt-…` | the pip install |
| no environment active | the release package (if its `source` line is in `~/.bashrc`) |

So if you have a release package in `~/.bashrc` and later `pip install` a newer
lsdtt3 into an environment, you get the newer tools inside that environment and
the older ones everywhere else, and a script can pick up either. Nothing warns
you. Check before a run that matters:

```bash
type -a lsdtt-raster-preprocessing     # every copy on PATH, the one used first
lsdtt-raster-preprocessing -v          # its version, and whose PROJ data it uses
hash -r                                 # after changing PATH, if bash still runs the old one
```

**Our advice: use one route per machine.** If you move to pip, remove the
`source …/setenv.sh` line from `~/.bashrc` and delete the package folder. Keep
both only on purpose, for example an older release kept for a published
analysis, and then call its tools by full path
(`~/lsdtt3-backend-…-v0.9.3/bin/lsdtt-…`). Do not mix versions within one
analysis: the Python package and the tools it calls share one version number,
and so do their parameter names and outputs.

`setenv.sh` also puts the package's `lib/` on `LD_LIBRARY_PATH`, which affects
other programs started from that shell; see
[If you already have GDAL or PROJ installed](#if-you-already-have-gdal-or-proj-installed).
The pip install's own tools are not affected by it: they load their libraries
from the wheel.

> **macOS, v0.10.0:** both the macOS release packages and the macOS wheel of
> v0.10.0 were built without PROJ data, so coordinate conversions fail. Use a
> later release.

---

## Linux package

This is the main way to install, and in other operating systems you just end up emulating Linux, so it works for all the pathways. 

### 1. Download and unpack

```bash
wget https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-linux-x86_64-core.tar.gz
tar -xzf lsdtt3-backend-linux-x86_64-core.tar.gz
ls -d lsdtt3-backend-linux-x86_64-core-v*      # the folder name shows the version
```

### 2. Activate it

```bash
source lsdtt3-backend-linux-x86_64-core-v*/setenv.sh
```

(If you have unpacked more than one version in the same place, give the full
folder name instead of the `*`.)

`setenv.sh` must be **`source`d**, not run with `bash`. All it does is export
environment variables (`PATH`, `LD_LIBRARY_PATH`, `GDAL_DATA`, `PROJ_DATA`),
and those vanish with the child process if you execute it instead.

### 3. Check it worked

```bash
lsdtt-channel-extraction --help
```

A list of parameters means you are done.

### 4. Make it permanent (optional)

`setenv.sh` applies to the current terminal only. To have the tools available
in every new terminal, append the same line, with the **full** path, to the
end of your `~/.bashrc`:

```bash
# run this in the folder where you unpacked it; it writes the full path
echo "source $(cd lsdtt3-backend-linux-x86_64-core-v* && pwd)/setenv.sh" >> ~/.bashrc
```

After installing a newer release, edit that line in `~/.bashrc` to point at
the new folder.

### What you need

Very little: a 64-bit x86 Linux system. The package is built on Ubuntu
22.04, which sets the oldest glibc it will run against — anything that age or
newer works. It does **not** bundle glibc itself, so a distribution older than
Ubuntu 22.04 may refuse to start it.

### Debian package (experimental, untested)

The repository also contains a Debian package definition under `debian/`, which
builds a `.deb` linking the **system** GDAL, PROJ, FFTW and netCDF through
`Depends:` rather than bundling them. That keeps the package small and lets it
track your distribution's security updates, which is the opposite trade-off
from the tarball above.

!!! warning "This route is not tested"

    The `.deb` build has not been through the testing the tarball packages
    have, and no `.deb` is published. It is included because the definition
    exists and works for the developers, not because it is ready to rely on.
    **Use the tarball above unless you specifically want the system-library
    build and are prepared to debug it.** If you try it and something breaks,
    a bug report is genuinely useful.

Building it is a developer task rather than an install step — the procedure,
the dependency list and the smoke test are in
[`debian_packaging.md`](../developer_notes/debian_packaging.md).

---

## Windows, via WSL2

lsdtt3 runs best on Windows through **WSL2** (Windows Subsystem for Linux),
which gives you a real Ubuntu terminal inside Windows. You then use the
[Linux package](#linux-package) exactly as above.

This is the recommended Windows route: you get every tool, and in our testing
the Linux build running under WSL2 is *faster* than a native Windows build on
the same machine. We don't understand how that can be the case but friends, strange things happen. 

### 1. Install WSL2 and Ubuntu

1. Open PowerShell **as administrator** (right-click the Start menu →
   Terminal (Admin)). You may need to become friendly with your IT team to be able to do this. 
2. Run `wsl --install`.
3. Restart when prompted.
4. Open the Microsoft Store, search for **Ubuntu 22.04 LTS**, install it.
5. Launch Ubuntu from the Start menu. The first launch asks you to create a
   username and password (don't forget your password!) the password is what `sudo` will ask for later.

### 2. Install the Linux package

Follow [Linux package](#linux-package) inside the Ubuntu terminal.

If you downloaded the tarball with a Windows browser it is on your Windows
disk, which WSL2 mounts under `/mnt/`. A file at
`C:\Users\Alice\Downloads\lsdtt3-backend-linux-x86_64-core.tar.gz` is:

```
/mnt/c/Users/Alice/Downloads/lsdtt3-backend-linux-x86_64-core.tar.gz
```

You will need to update this path unless your username is Alice.

### 3. Working with your Windows files

The same `/mnt/` mapping lets the tools read and write DEMs that live on your
Windows drives, so you can keep your data where your GIS expects it:

```bash
lsdtt-raster-preprocessing /mnt/c/Users/Alice/data/ my_analysis.driver
```

Note that reading and writing across `/mnt/c` is noticeably slower than
working inside the Linux filesystem (`~/`). For large DEMs, copy the data into
your Linux home directory first.

(The developers of `lsdtt3` try to avoid using a GIS since it takes forever to start, so they use [`lsdviztools3`](https://lsdtopotools.github.io/lsdviztools3/) instead).

### WSL2 troubleshooting

**No internet inside Ubuntu** (downloads fail, `apt` cannot reach anything).
Set the DNS server manually:

```bash
sudo nano /etc/resolv.conf
```

Make the file contain exactly:

```
nameserver 8.8.8.8
```

Save with `Ctrl-O`, exit with `Ctrl-X`. If the file reverts after a reboot,
also create `/etc/wsl.conf` containing:

```
[network]
generateResolvConf = false
```

---

## Windows native package

A native Windows build exists, distributed as a `.zip`. **It is a reduced
build** and exists mainly as the backend for the QGIS plugin:

- it ships the `qgis` profile, not the full tool set;
- Eigen-dependent routines are disabled, so some analyses are
  missing or skip their optional steps;
- it is generally *slower* than the same machine running the Linux package
  under WSL2.

For every tool natively on Windows, use `pip install lsdtt3` instead: the
Windows wheel carries all the tools, built with Eigen. See
[Python package](#python-package).

Unless you specifically need `.exe` files — for the QGIS plugin, or because
WSL2 is unavailable on a managed machine — use [WSL2](#windows-via-wsl2).

Download [`lsdtt3-backend-win64-qgis.zip`](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-win64-qgis.zip)
(always the latest release) and unzip it (right-click, **Extract All…**). Then
in PowerShell, in the folder you extracted to:

```powershell
cd lsdtt3-backend-win64-qgis-v*
. .\setenv.ps1
.\bin\lsdtt-raster-preprocessing.exe --version
```

Note the leading `. ` (dot, space) — like `source` on Linux, it runs
`setenv.ps1` in your current shell so the environment it sets survives.

---

## macOS

There are three ways to run lsdtt3 on a Mac. Choose by what you need:

| You want | Use |
|---|---|
| The QGIS plugin | The [native macOS package](#native-macos-package-qgis-tool-set) |
| lsdtt3 from Python | The [Python package](#python-package), which runs natively |
| Every command-line tool | The [Linux package in a virtual machine](#the-full-tool-set-in-a-linux-vm) |

### Native macOS package (QGIS tool set)

A native package is built for every release from 0.9.1, for Apple Silicon
(M1 and later) and Intel Macs, and needs macOS 15 or later. **It is a reduced
build:** it ships the `qgis` profile, the tools the
[QGIS plugin](#qgis-plugin) drives, not the full tool set. It is new, so please
report any problem with it.

1. Download the package for your Mac (these links always fetch the latest
   release): [Apple Silicon](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-macos-arm64-qgis.tar.gz) or
   [Intel](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-macos-x86_64-qgis.tar.gz). **Apple menu → About
   This Mac** shows which you have: an "Apple M…" chip is Apple Silicon.
2. Open Terminal in the folder you downloaded it to (in Finder, right-click
   the folder, then **New Terminal at Folder**) and unpack it:

   ```bash
   tar -xzf lsdtt3-backend-macos-arm64-qgis.tar.gz
   ```

3. Clear the quarantine flag the browser put on the download. Without this,
   macOS blocks the tools because they are not from the App Store:

   ```bash
   xattr -drs com.apple.quarantine lsdtt3-backend-macos-arm64-qgis-v*
   ```

   `-s` acts on symbolic links themselves instead of following them. Without
   it, `xattr` stops at the first link it cannot follow. In the v0.10.0 macOS
   packages that is `share/proj`, with the error
   `xattr: No such file: .../share/proj`.

   **v0.10.0 on macOS only:** that `share/proj` is a broken link, so the
   package has no PROJ data, and anything that converts between coordinate
   systems fails. Use a later release. If you must use v0.10.0 and have
   Homebrew PROJ 9.9 (arm64) or 9.8 (Intel), point the link at it:
   `ln -sfn "$(brew --prefix proj)/share/proj" lsdtt3-backend-macos-arm64-qgis-v0.10.0/share/proj`.

4. Activate the tools for this Terminal window, and check they run:

   ```bash
   source lsdtt3-backend-macos-arm64-qgis-v*/setenv.sh
   lsdtt-raster-preprocessing --version
   ```

   Use `source`, not `bash`: running the script in a new shell does not set up
   the current one. Repeat it in each new Terminal window, or add the `source`
   line, with the full path, to your `~/.zshrc`.

For Intel Macs, use `x86_64` in place of `arm64` in these names. The
package's own `README.md` repeats these steps. To use it with QGIS, select the
unpacked folder in the plugin's settings; see [QGIS plugin](#qgis-plugin). For
the tools the package does not include, `pip install lsdtt3` (macOS 15 or later;
see [Python package](#python-package)) gives all of them natively, or use the VM
below.

### The full tool set in a Linux VM

For every command-line tool, run the Linux package inside a lightweight Ubuntu
virtual machine using [Multipass](https://multipass.run), which Canonical
maintain for exactly this.

```bash
brew install --cask multipass
multipass launch --name lsdtt --cpus 4 --memory 8G --disk 40G 22.04
multipass shell lsdtt
```

You are now at an Ubuntu 22.04 prompt — follow
[Linux package](#linux-package) from there.

To share a macOS folder into the VM so the tools can reach your data:

```bash
multipass mount ~/topographic_data lsdtt:/home/ubuntu/data
```

Stop and restart the VM later with `multipass stop lsdtt` /
`multipass start lsdtt`.

---

## Google Colab

**The simplest way is pip.** Colab's Python is 3.12 or later, so one cell
gives you the Python package and every command-line tool, with nothing else to
set:

```python
!pip install lsdtt3
!lsdtt-channel-extraction -v
```

**Or use the Linux release package.** Colab already runs Linux, so the Linux
package works there with nothing installed on your own machine at all. Two
things differ from a normal terminal.

Download and unpack the latest release in a code cell (`!` runs a shell
command):

```python
!wget -q https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-linux-x86_64-core.tar.gz
!tar -xzf lsdtt3-backend-linux-x86_64-core.tar.gz
```

**Do not use `!source setenv.sh`.** Every `!` line runs in its own throwaway
shell, so nothing it sets survives to the next cell. Set the environment from
Python instead, once per session:

```python
import glob, os

root = sorted(glob.glob("/content/lsdtt3-backend-linux-x86_64-core-v*"))[-1]
print("Using", root)
os.environ["PATH"] = f"{root}/bin:" + os.environ["PATH"]
os.environ["LD_LIBRARY_PATH"] = f"{root}/lib:" + os.environ.get("LD_LIBRARY_PATH", "")
os.environ["GDAL_DATA"] = f"{root}/share/gdal"
os.environ["PROJ_DATA"] = f"{root}/share/proj"
os.environ["PROJ_LIB"]  = os.environ["PROJ_DATA"]
```

After that cell the tools work anywhere in the notebook:

```python
!lsdtt-channel-extraction --help
```

Remember that Colab sessions are wiped when they end — you will repeat both
cells next time, and anything you want to keep must be copied to Google Drive.

---

## QGIS plugin

The plugin adds an **lsdtt3** group to the QGIS Processing Toolbox. It runs the
command-line tools, so it needs a backend package as well:

1. Install a backend package for your system: the
   [Windows native package](#windows-native-package) on Windows, the
   [Linux package](#linux-package) on Linux, or the macOS `qgis` package
   ([Apple Silicon](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-macos-arm64-qgis.tar.gz) or
   [Intel](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-backend-macos-x86_64-qgis.tar.gz)) on macOS. The macOS
   package's `README.md` covers unpacking it and clearing the browser's
   quarantine flag.
2. Download the plugin from the same release:
   [`lsdtt3-qgis-plugin.zip`](https://github.com/LSDtopotools/lsdtt3/releases/latest/download/lsdtt3-qgis-plugin.zip) is always the
   latest. Download it from GitHub, not from any other copy you come across:
   an older plugin may refuse a newer backend.
3. In QGIS: **Plugins → Manage and Install Plugins… → Install from ZIP**, and
   choose the zip.
4. Open **lsdtt3 → LSDTT Backend Settings** and select the backend package's
   root folder (the one containing `bin/`, `share/` and `manifest.json`).

**The plugin and the backend must come from the same `X.Y` series**: plugin
0.7.x works with any 0.7.x backend. A mismatch is refused before anything
runs, with a message saying whether to update the plugin or the backend.
Plugins up to 0.4.0 predate this rule; if one reports that the backend
"supports plugin versions …", update the plugin.

To update either, repeat its install step: **Install from ZIP** replaces the
plugin in place, and a new backend goes in a new folder, which you then select
in the settings.

---

## Python package

The `lsdtt3` Python package wraps the same C++ code for use from scripts and
notebooks, keeping rasters and results in memory instead of files. It shares
the lsdtt3 version number, and is built for Linux, macOS and Windows on every
release.

The bindings cover part of what lsdtt3 does; the complete software is the
command-line tools. The package installs all of those as well (see
[The command-line tools from pip](#the-command-line-tools-from-pip)), so with
Python 3.12 or later, `pip install lsdtt3` gets you everything.

Install it from PyPI, into a virtual, conda or pixi environment rather than
your system Python:

```bash
python -m pip install lsdtt3
```

It needs Python 3.12 or later, on Linux (x86_64), macOS (Apple Silicon and
Intel, macOS 15 or later) or Windows (x64). There is one wheel per platform,
`cp312-abi3`: built against Python's stable ABI, the same file installs on
3.12, 3.13, 3.14 and later.

**In Google Colab**, run this in a code cell:

```python
!pip install lsdtt3
```

**From a GitHub release.** The same wheels are attached to each
[GitHub release](https://github.com/LSDtopotools/lsdtt3/releases/latest), as
are the wheels of versions up to 0.9.3, which are not on PyPI (those have a
separate wheel for each Python version, `cp310` to `cp314`). pip reads the
Python and platform from the file name, so keep the name as it is. This
command, run with the Python you will use, prints the address of the matching
wheel in the latest release:

```bash
python - <<'EOF'
import json, platform, sys, urllib.request
if sys.version_info < (3, 12):
    sys.exit("lsdtt3 needs Python 3.12 or later")
rel = json.load(urllib.request.urlopen(
    "https://api.github.com/repos/LSDtopotools/lsdtt3/releases/latest"))
plat = {"linux": "manylinux", "darwin": "macosx", "win32": "win_"}[sys.platform]
arch = platform.machine().lower()      # x86_64, arm64 or amd64
for a in rel["assets"]:
    n = a["name"]
    if n.endswith(".whl") and "-abi3-" in n and plat in n and arch in n.lower():
        print(a["browser_download_url"])
EOF
```

Then install from that address:

```bash
python -m pip install https://github.com/LSDtopotools/lsdtt3/releases/download/vX.Y.Z/lsdtt3-X.Y.Z-cp312-abi3-manylinux_2_28_x86_64.whl
```

If it prints nothing, that release has no wheel for your platform, or your
Python is older than 3.12; the [list of releases](https://github.com/LSDtopotools/lsdtt3/releases)
shows which wheels each one has.

You do not need the backend tarball or the `os.environ` cell from
[Google Colab](#google-colab) for the Python package; the wheel finds its own
GDAL and PROJ data.

The wheel is self-contained: it bundles GDAL and PROJ, and also installs the
command-line tools as commands. Install it into a virtual environment, conda
or pixi environment, not into your system Python.

### The command-line tools from pip

After `python -m pip install lsdtt3`, every `lsdtt-*` tool is a command in
that environment, used exactly as with the release package (parameter files,
`-h`, `-v`); there is no `setenv.sh` to source. For example, with pixi:

```bash
pixi init lsdtt && cd lsdtt
pixi add "python>=3.12" pip
pixi run python -m pip install lsdtt3
pixi run lsdtt-raster-preprocessing -v
pixi shell          # a shell where the lsdtt-* commands are on PATH
```

or with a virtual environment:

```bash
python3.12 -m venv ~/lsdtt3-env
source ~/lsdtt3-env/bin/activate      # Windows: ~\lsdtt3-env\Scripts\activate
python -m pip install lsdtt3
lsdtt-raster-preprocessing -v
```

The `-v` output should show `proj data: bundled:` with a path inside the
environment's `site-packages/lsdtt3`. The commands are only on `PATH` while the
environment is active. See [Release package or pip?](#release-package-or-pip)
for how this compares with the release package.

Check what you have with:

```bash
python -c "import lsdtt3; print(lsdtt3.__version__, lsdtt3._core.lsdtt3_version())"
```

To use the Python package from unreleased code, build it from a clone; see
[Building the Python wheels](../developer_notes/building_python_wheels.md).
The [Python guide](../python/index.md) shows what to do once it is installed.
Every function is documented in the package itself (`help(lsdtt3.FlowInfo)`),
and `print(lsdtt3.helper_for_llm())` prints a compact guide written for you
or for an AI assistant such as Gemini in Colab. It lists topics, and
`helper_for_llm("symbol:<name>")` gives the exact arguments of one function.
To paste the whole guide into a chat, use `print(lsdtt3.helper_for_llm("all"))`.

---

## What is inside a package

| Path | Contents |
|---|---|
| `bin/` | the `lsdtt-*` command-line tools |
| `lib/` | every shared library they need (GDAL, PROJ, …) |
| `share/gdal`, `share/proj` | support data those libraries read at runtime |
| `share/lsdtt3` | lsdtt3's own data (the cosmogenic-nuclide parameter file) |
| `setenv.sh` / `setenv.ps1` | sets `PATH` and the library/data variables |
| `manifest.json` | version and a SHA256 checksum of every file |
| `README.md` | the same instructions, offline |
| `LICENSE`, `THIRD_PARTY_NOTICES.txt` | the AGPL licence, and the licences of the bundled libraries |

Nothing is installed system-wide. To uninstall, delete the folder (and the
`source` line from `~/.bashrc` if you added one).

---

## Troubleshooting

**`command not found` after sourcing `setenv.sh`**
Check the bundle's `bin` is on your `PATH`:

```bash
echo $PATH | tr ':' '\n' | grep lsdtt3
```

If nothing prints, you probably ran `setenv.sh` instead of `source`-ing it, or
you are in a different terminal from the one where you sourced it.

**`error while loading shared libraries`**
The `lib/` folder is not being found. This happens if you moved or renamed the
bundle directory after sourcing `setenv.sh` — re-source it from its new
location.

**`version 'GLIBC_2.35' not found`**
Your Linux distribution is older than the Ubuntu 22.04 the package was built
against. Use a newer distribution, or build from source.

---

### If you already have GDAL or PROJ installed

This is the most common source of confusing behaviour, and it gets more likely
the more serious your GIS setup is. If you have QGIS, conda, `rasterio`,
`geopandas` or a system `libgdal` — and anyone doing large-scale analysis
does — then there is more than one GDAL and more than one PROJ on your machine.

**The good news:** lsdtt3 pins its own data directories at startup, in the
release packages and in the pip wheel alike. A `proj.db` bundled next to the
binary always wins over `PROJ_DATA`, `PROJ_LIB` or `GDAL_DATA` in your
environment. You do not have to unset anything, and you do
not have to source `setenv.sh` for the tools to find their own data — that
matters in places like Google Colab where each `!` line is a separate shell.

**Check what is actually in use** with any tool's `-v` (shown here as
0.7.0 and later print it; the first line names the lsdtt3 version and whether
the tool is `stable` or `experimental`):

```bash
lsdtt-raster-preprocessing -v
```

```
lsdtt-raster-preprocessing (lsdtt3 0.7.1, stable)
GDAL 3.8.4, PROJ 9.3.1
  proj data: bundled: /home/you/lsdtt3/share/proj
  GDAL data: bundled: /home/you/lsdtt3/share/gdal
```

`bundled:` means lsdtt3 is using its own copy and your environment is being
ignored, which is what you want from a package. If it says `PROJ_DATA=…` or
`PROJ built-in search` instead, it is using something else — expected for a
source build, worth knowing about for a package or a pip install.

**The other direction.** `setenv.sh` exports the package's `lib/`,
`GDAL_DATA` and `PROJ_DATA` into the shell you source it in. GDAL, PROJ or QGIS
programs that you then start from that same shell may pick up lsdtt3's copies,
and an older PROJ can fail on a newer `proj.db`. If that happens, run them from
a shell where `setenv.sh` was not sourced, or use the pip route, which changes
nothing outside the `lsdtt-*` commands.

**Why we care so much about this.** A mismatched `proj.db` does not usually
produce an error. It produces *coordinates*. They will be the right shape, in
the right part of the world, and wrong by tens of metres, because the datum
transformation your CRS needs was resolved differently or not at all. That is
much worse than a crash, and it is why the version banner reports the data
directories at all.

**Please include the whole `-v` output in any bug report** that involves
projections, reprojection, or coordinates that look subtly off. It is the first
thing we will ask for.

**Symptoms that point here:**

| What you see | What it usually means |
|---|---|
| `Cannot find proj.db` | `PROJ_DATA` or `PROJ_LIB` points somewhere without one, and nothing is bundled |
| `Cannot find proj.db` *from a conda environment* | conda's PROJ and a system PROJ disagree; activate the environment, or use a release package instead |
| Coordinates off by ~100 m in a national grid | a datum shift was not applied — a wrong `proj.db`, not a bug in the analysis |
| `undefined symbol: GDAL…` when building | headers from one GDAL, library from another — see the build-from-source notes |
| GeoPackage writes fail under WSL2 on `/mnt/c` | SQLite on an NTFS mount; write FlatGeobuf instead |

---

## Reporting bugs

> **This is software built by academics for academic use within a research
> group.** It is not supported commercial software. The developers are happy to
> help with questions, comments and bug reports, but they may be teaching or
> doing research, or sitting in administrative meetings, or filling in tedious paperwork, or reviewing your grant or paper.  Please adjust expectations of response times accordingly.

Please include:

1. The exact command you ran.
2. The full terminal output — copied and pasted as text, not a screenshot.
3. Your Linux version: `cat /etc/os-release | grep VERSION=`
4. The lsdtt3 version: `lsdtt-channel-extraction -v`

Open an issue on the GitHub repository if you have access, or send it to the
developers.

---

## See also

- [Building lsdtt3 from source](../developer_notes/building_from_source.md)
  — compiling lsdtt3 yourself, on Linux, Windows or macOS
- [Making release packages](../developer_notes/making_release_packages.md)
  — how these packages are built and published
- [Command-line reference](../cli_reference/cli_reference_lsdtt_raster_preprocessing.md)
  — what to do once it is installed
