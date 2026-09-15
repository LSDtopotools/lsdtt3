# Installing lsdtt3

Every route on this page installs a **pre-built release package**. You do not
need a compiler, you do not need to install GDAL or any other library, and on
Linux you do not need `sudo` — the package carries the command-line tools and
every library they depend on.

> **Want the very latest code instead of a release?**
> Then you need to compile it yourself. That is a developer activity and is
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
| <https://www.geos.ed.ac.uk/~smudd/lsdtt_packages/> | A plain download directory, no account needed. |

File names look like this:

```
lsdtt3-backend-linux-x86_64-core-v0.4.0.tar.gz     Linux
lsdtt3-backend-win64-qgis-v0.4.0.zip               Windows (native)
```

`core` and `qgis` are **profiles** — which set of command-line tools the
package contains. `core` is the general-purpose set; `qgis` is the smaller
subset the QGIS plugin drives.

> **The version numbers on this page are examples and will go stale.** Every
> command below spells out a version because a bare filename would not work,
> but do not assume it is the current one. The authoritative version is line 2
> of `CMakeLists.txt` in the repository — `project(lsdtt3 VERSION ...)` — and
> release tags and package filenames are derived from it directly (`VERSION
> 0.4.0` → tag `v0.4.0` → `...-v0.4.0.tar.gz`). If you cannot see the
> repository, take the newest file from one of the two download locations
> above instead. Substitute that number throughout.

---

## Linux package

This is the main route, and the one every other route below ends up using.

### 1. Download and unpack

```bash
wget https://www.geos.ed.ac.uk/~smudd/lsdtt_packages/lsdtt3-backend-linux-x86_64-core-v0.4.0.tar.gz
tar -xzf lsdtt3-backend-linux-x86_64-core-v0.4.0.tar.gz
```

### 2. Activate it

```bash
source lsdtt3-backend-linux-x86_64-core-v0.4.0/setenv.sh
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
in every new terminal, append the same line — with the **full** path — to the
end of your `~/.bashrc`:

```bash
echo "source $HOME/lsdtt3-backend-linux-x86_64-core-v0.4.0/setenv.sh" >> ~/.bashrc
```

### What you need

Almost nothing: a 64-bit x86 Linux system. The package is built on Ubuntu
22.04, which sets the oldest glibc it will run against — anything that age or
newer works. It does **not** bundle glibc itself, so a distribution older than
Ubuntu 22.04 may refuse to start it.

---

## Windows, via WSL2

lsdtt3 runs best on Windows through **WSL2** (Windows Subsystem for Linux),
which gives you a real Ubuntu terminal inside Windows. You then use the
[Linux package](#linux-package) exactly as above.

This is the recommended Windows route: you get every tool, and in our testing
the Linux build running under WSL2 is *faster* than a native Windows build on
the same machine.

### 1. Install WSL2 and Ubuntu

1. Open PowerShell **as administrator** (right-click the Start menu →
   Terminal (Admin)).
2. Run `wsl --install`.
3. Restart when prompted.
4. Open the Microsoft Store, search for **Ubuntu 22.04 LTS**, install it.
5. Launch Ubuntu from the Start menu. The first launch asks you to create a
   username and password — the password is what `sudo` will ask for later.

### 2. Install the Linux package

Follow [Linux package](#linux-package) inside the Ubuntu terminal.

If you downloaded the tarball with a Windows browser it is on your Windows
disk, which WSL2 mounts under `/mnt/`. A file at
`C:\Users\Alice\Downloads\lsdtt3-backend-linux-x86_64-core-v0.4.0.tar.gz` is:

```
/mnt/c/Users/Alice/Downloads/lsdtt3-backend-linux-x86_64-core-v0.4.0.tar.gz
```

### 3. Working with your Windows files

The same `/mnt/` mapping lets the tools read and write DEMs that live on your
Windows drives, so you can keep your data where your GIS expects it:

```bash
lsdtt-raster-preprocessing /mnt/c/Users/Alice/data/ my_analysis.driver
```

Note that reading and writing across `/mnt/c` is noticeably slower than
working inside the Linux filesystem (`~/`). For large DEMs, copy the data into
your Linux home directory first.

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
cd lsdtt3-backend-win64-qgis-v0.4.0
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
!wget https://www.geos.ed.ac.uk/~smudd/lsdtt_packages/lsdtt3-backend-linux-x86_64-core-v0.4.0.tar.gz
!tar -xzf lsdtt3-backend-linux-x86_64-core-v0.4.0.tar.gz
```

**Do not use `!source setenv.sh`.** Every `!` line runs in its own throwaway
shell, so nothing it sets survives to the next cell. Set the environment from
Python instead, once per session:

```python
import os

root = "/content/lsdtt3-backend-linux-x86_64-core-v0.4.0"
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

## Reporting bugs

> **This is software built by academics for academic use within a research
> group.** It is not supported commercial software. The developers are happy to
> help with questions, comments and bug reports, but they may be teaching or
> doing research — please adjust expectations of response times accordingly.

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
