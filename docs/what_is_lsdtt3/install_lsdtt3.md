# Installing lsdtt3

Every route on this page installs a **pre-built release package**. You do not
need a compiler, you do not need to install GDAL or any other library, and on
Linux you do not need `sudo` — the package carries the command-line tools and
every library they depend on.

> **Want the very latest code instead of a release?**
> Then you need to compile it yourself. If you are not a developer do not do this. That activity is
> documented separately in
> [`docs/developer_notes/building_from_source.md`](../developer_notes/building_from_source.md).

---

## Choose your route

**1. Where do you want to run lsdtt3?**

| | |
|---|---|
| On my own computer | continue to question 2 |
| In a browser, with nothing installed | [Google Colab](#google-colab) |

**2. What is your computer running?**

| | |
|---|---|
| Linux | [Linux package](#linux-package) |
| Windows | [Windows, via WSL2](#windows-via-wsl2) — **recommended**, full functionality |
| Windows, and I want native `.exe` files | [Windows native package](#windows-native-package) — reduced tool set |
| macOS | [macOS](#macos) |

---

## Where the packages live

Both of these hold the same files. Use whichever you can reach.

| Source | Notes |
|---|---|
| [GitHub Releases](https://github.com/LSDtopotools/lsdtt3/releases/latest) | The canonical location. While the repository is in private beta you need to be a collaborator to see it. |
| <https://www.geos.ed.ac.uk/~smudd/lsdtt_packages/> | A plain download directory, no account needed. It is annoying for Simon to upload these so they go out of date.  |

File names look like this:

```
lsdtt3-backend-linux-x86_64-core-v0.5.1.tar.gz     Linux
lsdtt3-backend-win64-qgis-v0.5.1.zip               Windows (native)
```

`core` and `qgis` are **profiles** — which set of command-line tools the
package contains. `core` is the general-purpose set; `qgis` is the smaller
subset the QGIS plugin drives.

> **The version numbers on this page are examples and will go stale.** Every
> command below spells out a version because a bare filename would not work,
> but don't assume it is the current one. The authoritative version is line 2
> of `CMakeLists.txt` in the repository — `project(lsdtt3 VERSION ...)` (we do not recommend reading beyond this line or you will get a headache) and
> release tags and package filenames are derived from it directly (`VERSION
> 0.5.1` → tag `v0.5.1` → `...-v0.5.1.tar.gz`). If you cannot see the
> repository, take the newest file from one of the two download locations
> above instead. Substitute that number throughout.

---

## Linux package

This is the main way to install, and in other operating systems you just end up emulating Linux, so it works for all the pathways. 

### 1. Download and unpack

```bash
wget https://www.geos.ed.ac.uk/~smudd/lsdtt_packages/lsdtt3-backend-linux-x86_64-core-v0.5.1.tar.gz
tar -xzf lsdtt3-backend-linux-x86_64-core-v0.5.1.tar.gz
```

### 2. Activate it

```bash
source lsdtt3-backend-linux-x86_64-core-v0.5.1/setenv.sh
```

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
echo "source $HOME/lsdtt3-backend-linux-x86_64-core-v0.5.1/setenv.sh" >> ~/.bashrc
```

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
`C:\Users\Alice\Downloads\lsdtt3-backend-linux-x86_64-core-v0.5.1.tar.gz` is:

```
/mnt/c/Users/Alice/Downloads/lsdtt3-backend-linux-x86_64-core-v0.5.1.tar.gz
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

(The developers of `lsdtt3` try to avoid using a GIS since it takes forever to start, so they use `lsdviztools3` instead).

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

Unless you specifically need `.exe` files — for the QGIS plugin, or because
WSL2 is unavailable on a managed machine — use [WSL2](#windows-via-wsl2).

Unzip it, then in PowerShell:

```powershell
cd lsdtt3-backend-win64-qgis-v0.5.1
. .\setenv.ps1
.\bin\lsdtt-raster-preprocessing.exe --version
```

Note the leading `. ` (dot, space) — like `source` on Linux, it runs
`setenv.ps1` in your current shell so the environment it sets survives.

---

## macOS

**There is no macOS package yet.** One is planned; until then, run the Linux
package inside a lightweight Ubuntu virtual machine using
[Multipass](https://multipass.run), which Canonical maintain for exactly this.

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

Colab already runs Linux, so the Linux package works there with nothing
installed on your own machine at all. Two things differ from a normal terminal.

Download and unpack in a code cell (`!` runs a shell command):

```python
!wget https://www.geos.ed.ac.uk/~smudd/lsdtt_packages/lsdtt3-backend-linux-x86_64-core-v0.5.1.tar.gz
!tar -xzf lsdtt3-backend-linux-x86_64-core-v0.5.1.tar.gz
```

**Or get the URL from https://www.geos.ed.ac.uk/~smudd/lsdtt_packages** 

**Do not use `!source setenv.sh`.** Every `!` line runs in its own throwaway
shell, so nothing it sets survives to the next cell. Set the environment from
Python instead, once per session:

```python
import os

root = "/content/lsdtt3-backend-linux-x86_64-core-v0.5.1"
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

## What is inside a package

| Path | Contents |
|---|---|
| `bin/` | the `lsdtt-*` command-line tools |
| `lib/` | every shared library they need (GDAL, PROJ, …) |
| `share/gdal`, `share/proj` | support data those libraries read at runtime |
| `setenv.sh` / `setenv.ps1` | sets `PATH` and the library/data variables |
| `manifest.json` | version and a SHA256 checksum of every file |
| `README.md` | the same instructions, offline |

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

**The good news:** lsdtt3 pins its own data directories at startup. A `proj.db`
bundled next to the binary always wins over `PROJ_DATA`, `PROJ_LIB` or
`GDAL_DATA` in your environment. You do not have to unset anything, and you do
not have to source `setenv.sh` for the tools to find their own data — that
matters in places like Google Colab where each `!` line is a separate shell.

**Check what is actually in use** with any tool's `-v`:

```bash
lsdtt-raster-preprocessing -v
```

```
lsdtt-raster-preprocessing version 1.0.0 (lsdtt3 0.5.1)
GDAL 3.8.4, PROJ 9.3.1
  proj data: bundled: /home/you/lsdtt3/share/proj
  GDAL data: bundled: /home/you/lsdtt3/share/gdal
```

`bundled:` means lsdtt3 is using its own copy and your environment is being
ignored, which is what you want from a package. If it says `PROJ_DATA=…` or
`PROJ built-in search` instead, it is using something else — expected for a
source build or a conda install, worth knowing about for a package.

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

- [`docs/developer_notes/building_from_source.md`](../developer_notes/building_from_source.md)
  — compiling lsdtt3 yourself, on Linux, Windows or macOS
- [`docs/developer_notes/making_release_packages.md`](../developer_notes/making_release_packages.md)
  — how these packages are built and published
- [Command-line reference](../cli_reference/cli_reference_lsdtt_raster_preprocessing.md)
  — what to do once it is installed
