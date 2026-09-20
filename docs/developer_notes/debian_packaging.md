# Building a Debian package for lsdtt3

This document explains how to build a `.deb` package for distributing lsdtt3
binaries. The package installs the compiled CLI programs and the bundled
runtime data, so end users do not need to build from source.

The packaging lives in the repository, in `debian/`. You do not create it —
you check it is current for the release you are packaging, then build.

---

## Prerequisites

Install the Debian packaging tools on your build machine:

```bash
sudo apt-get install -y \
    build-essential \
    devscripts \
    debhelper \
    dpkg-dev \
    fakeroot \
    lintian
```

You also need all of lsdtt3's build dependencies installed — see
[Building from source](building_from_source.md), `apptainer_def/apptainer_setup.md`,
or the main README for the full list. `dpkg-buildpackage` will refuse to start
if anything in `Build-Depends` is missing, so a clean machine tells you what it
needs.

---

## What is in `debian/`

```
lsdtt3/
├── debian/
│   ├── control        — package metadata and dependencies
│   ├── rules          — build instructions (a Makefile, driven by debhelper)
│   ├── changelog      — version history (required by Debian policy)
│   ├── copyright      — DEP-5 copyright and the full AGPL-3 text
│   ├── smoke-test.sh  — post-build and post-install checks (see below)
│   └── source/format  — source package format
├── src/
├── CMakeLists.txt
└── ...
```

Those six files are the whole component, and all six are committed. Everything
else that appears under `debian/` during a build — `debian/lsdtt3/`,
`debian/.debhelper/`, `debian/files`, `*.substvars` — is generated and
gitignored.

There is deliberately **no `debian/compat`**. The compatibility level is set
once, by `debhelper-compat (= 13)` in `Build-Depends`. Specifying it in both
places is a fatal error — debhelper stops with *"Please specify the compatibility
level exactly once"*.

### Two traps on WSL

**You cannot build the package under `/mnt/c`.** The drive is mounted without
the `metadata` option, so every file reports mode 777 and `chmod` does nothing.
`dpkg-deb` refuses a control directory more permissive than 0775 and the build
dies at the very last step, after the full compile:

```
dpkg-deb: error: control directory has bad permissions 777 (must be >=0775 and <=0775)
```

Copy the tree to a Linux-native path and build there. Copy only what git
tracks, or you will also copy `.pixi/` and `build/`:

```bash
mkdir -p ~/lsdtt3-debbuild
git ls-files -z > /tmp/tracked.z
rsync -a --from0 --files-from=/tmp/tracked.z . ~/lsdtt3-debbuild/
cd ~/lsdtt3-debbuild && dpkg-buildpackage -us -uc -b -rfakeroot
```

**`debian/rules` and `debian/smoke-test.sh` must be executable**, and for the
same reason git cannot see the exec bit on this filesystem — `core.filemode` is
`false`, so git knows only what the index records. Any unstage-and-re-add
silently drops it back to `100644`, and neither `git status` nor `ls -l` will
tell you:

```bash
git update-index --chmod=+x debian/rules debian/smoke-test.sh
git ls-tree HEAD debian/        # both must read 100755
```

A `debian/rules` committed non-executable means `dpkg-buildpackage` fails
immediately for everyone who clones the repository.

### debian/control

`Build-Depends` lists **every** optional dependency, not just the required
ones. This matters more here than in a normal package, because this build
degrades *silently* when a library is absent — see
[Optional dependencies degrade silently](#optional-dependencies-degrade-silently)
below.

`Depends` is just `${shlibs:Depends}, ${misc:Depends}`. There is no hand-written
list of library packages: `dpkg-shlibdeps` reads the linked binaries and works
out the correct runtime packages and versions for whichever Ubuntu release you
built on. A hand-written list would have to name `libgdal30` on Ubuntu 22.04,
`libgdal34t64` on 24.04, and so on; letting `shlibdeps` do it removes a whole
class of per-release mistakes.

The CLIs are **not** a `gdal-bin` dependency. They call gdalwarp through the
GDAL C API (`GDALWarpAppOptions` in `src/core/lsdtt_raster_reproject.cpp`), not
by running the command-line programs. GDAL's command-line tools are a *test*
requirement only.

### debian/rules

Five things are worth knowing about it.

**The tests are turned off.** The top-level `CMakeLists.txt` calls
`include(CTest)`, so `BUILD_TESTING` defaults to `ON`. Left alone, every
package build would compile the entire CTest suite and `dh_auto_test` would
then run it — including the `heavy` tests — so a slow or environment-sensitive
test would fail the package build. `override_dh_auto_configure` passes
`-DBUILD_TESTING=OFF`. Test before you package, not while you package.

**Nothing needs moving after install.** Every `install()` rule in
`CMakeLists.txt` uses a destination relative to `CMAKE_INSTALL_PREFIX`, and
debhelper's cmake buildsystem sets that to `/usr`. So `dh_auto_install` already
puts the CLIs in `debian/lsdtt3/usr/bin` and the bundled `crn_parameters.nc` in
`debian/lsdtt3/usr/share/lsdtt3` — which is exactly where the cosmo tools look
for it, since they resolve it as `<exe>/../share/lsdtt3`. Nothing in this build
installs to `/usr/local`; if you find packaging instructions that shuffle files
out of `/usr/local`, they predate this layout, and a `rm -rf` of that tree would
take the bundled cosmo parameter file with it.

**The clean step does not walk your source tree.** `dh_clean` normally finishes
with a `find` over everything, looking for editor leftovers (`*~`, `*.orig`,
`*.rej`). This repository produces none of those, but the walk descends into
`build/` and `.pixi/` — gigabytes — and on network storage or a Windows drive
mounted into WSL that single `find` takes longer than compiling the whole
package. `override_dh_clean` does the removals debhelper actually needs and
skips the sweep: 1.4 seconds instead of minutes. Nothing about the built package
changes; `dh_auto_clean` has already removed the CMake build directory by then.

If you are on a release of debhelper where that override misbehaves, the escape
hatch is `dpkg-buildpackage -nc`, which skips the clean step altogether. The
build is out of tree, so nothing stale is reused.

**Static libraries and headers are dropped.** `CMakeLists.txt` also installs
`liblsdtt3_core.a`, `liblsdtt3_utils.a`, `liblsdtt3_cpu.a` and the public
headers, which are useful for a development package but dead weight in a
binary one. `override_dh_auto_install` removes them. If you ever want an
`lsdtt3-dev` package, that is where to start.

**Smoke checks run in place of the test suite.** With `BUILD_TESTING=OFF`
nothing else in the build confirms that a single binary runs, so
`debian/smoke-test.sh` is called twice:

| Hook | Mode | Checks |
|---|---|---|
| `override_dh_auto_test` | `built <build-dir>` | every `lsdtt-*` starts, initialises GDAL and exits 0 on `--version`, and reports a GDAL version; at least 40 programs built; `lsdtt-spectral-analysis` exists (FFTW3 was found); `lsdtt-cosmo-tool` links `libnetcdf` (the cosmo corrections are in) |
| `override_dh_auto_install` | `staged <destdir>` | ≥40 programs in `usr/bin`; `usr/share/lsdtt3/crn_parameters.nc` present; no static libraries or headers left staged |

These take seconds and catch the packaging-specific failures a unit test
cannot see: a binary that will not start, a tool missing because an optional
dependency was absent at configure time, or the bundled cosmo data not reaching
its installed location. A successful run prints:

```
smoke-test: 40 programs started, reported GDAL/PROJ and exited 0
smoke-test: optional dependencies present: FFTW3 (spectral analysis), netCDF (cosmo corrections)
smoke-test: 40 programs in usr/bin, crn_parameters.nc present, no static libs or headers
```

---

## Building without root

On a shared machine where you cannot `apt-get install` anything, two separate
things block the build, and only one of them is about the compiler.

**`dpkg-checkbuilddeps` reads the dpkg database.** It has no knowledge of conda,
Spack, modules or anything else in your home directory, so a dependency
satisfied that way still aborts the build before CMake runs:

```
dpkg-checkbuilddeps: error: Unmet build dependencies: libeigen3-dev
```

That check is metadata only. Skip it with `-d`.

**A fresh configure finds nothing your cache remembers.** `debian/rules` builds
in a clean `obj-<arch>/` directory, so a `CMakeCache.txt` in your usual `build/`
proves nothing about what the package build will find — a cached `Eigen3_DIR`
keeps reporting `Eigen3 found` even after you deactivate the environment that
provided it. To see what the package build sees, configure into an empty
directory.

`CMakeLists.txt` looks for Eigen in **config mode only**
(`find_package(Eigen3 NO_MODULE)`), so it needs the directory containing
`Eigen3Config.cmake`. In a conda environment that is
`$CONDA_PREFIX/share/eigen3/cmake`. `debian/rules` takes it from
`LSDTT_EIGEN3_DIR`:

```bash
conda activate <env>
LSDTT_EIGEN3_DIR="$CONDA_PREFIX/share/eigen3/cmake" \
    dpkg-buildpackage -us -uc -b -rfakeroot -d
```

Confirm it worked by reading `Eigen3 : ON` in the configure block in the build
log — with `-d` in play, a silently missing Eigen would otherwise produce a
package without chi inversion.

### Why only Eigen, and not the whole prefix

It is tempting to pass `-DCMAKE_PREFIX_PATH=$CONDA_PREFIX` instead. Do not.
That also lets GDAL, FFTW and netCDF resolve inside the environment, and the
resulting package then depends on libraries that exist on no other machine: it
installs cleanly and fails at run time. `dpkg-shlibdeps` usually catches it, but
reports it as an unmappable path rather than as the mistake it is.

Eigen is the exception because it is **header-only** — it is compiled into the
binaries and leaves no entry in `ldd` output, so taking it from a private prefix
cannot affect what the package depends on. That is why it is safe to point at
conda and the others are not.

`debian/smoke-test.sh` enforces this: it fails the build if any binary links a
library from a path containing `conda`, `/envs/` or `virtualenv`, and prints

```
smoke-test: no binary links against a conda or virtualenv prefix
```

when the build is clean.

### The package is still a normal package

Nothing about building this way marks the result. Because only Eigen came from
the private prefix, the `.deb` depends on the same system libraries as any other
build, and runs on machines with no conda at all. The environment is needed to
**build**, not to run.

---

## Optional dependencies degrade silently

This is the failure mode to watch for. None of these stop the build — they
just quietly produce a smaller package:

| Missing build dependency | What you lose | What CMake says |
|---|---|---|
| `libnetcdf-dev` | ERA40 atmosphere and palaeomagnetic corrections in the cosmo tools; they fall back to the standard atmosphere | `STATUS` line only |
| `pkg-config` or `libfftw3-dev` | `lsdtt-spectral-analysis` is not built **at all** | `WARNING` |
| `libeigen3-dev` | chi / k_sn inversion in `lsdtt-chi-analysis` | `WARNING` |

All three are in `Build-Depends`, so a correctly configured build machine gets
them. Confirm it anyway by reading the configuration block CMake prints in the
build log:

```
=== lsdtt3 build configuration ===
  Build type  : Release
  C++ std     : 17
  GDAL        : 3.4.1
  OpenMP      : ON
  FFTW3       : ON
  GPU (CUDA)  : OFF
  Eigen3      : ON
===================================
```

A package built with `FFTW3 : FALSE` is missing a tool and will not tell the
user why.

---

## Versioning

`debian/changelog` is the authoritative version for the package, and it should
track `project(lsdtt3 VERSION ...)` in `CMakeLists.txt`. Bump both together.

The packaging is **native format** (`debian/source/format` is `3.0 (native)`),
because it lives in the upstream tree rather than being applied to a separate
upstream tarball. That means the version has no Debian revision: `0.5.1`, not
`0.5.1-1`. To re-release the same code with changed packaging, append a
packaging suffix — `0.5.1+pkg2`.

Add a new entry with `dch` rather than editing the file by hand; the format is
strict:

```bash
dch --newversion 0.5.2 --distribution jammy "Release 0.5.2."
dch --release ""
```

**The distribution field must name the Ubuntu series you are building for**, not
a Debian suite. Leaving it as `unstable` makes lintian report
`E: bad-distribution-in-changes-file`. Since a separate `.deb` is built on each
target release anyway (see [Building for multiple Ubuntu
versions](#building-for-multiple-ubuntu-versions)), set it to match the machine
you are on — `jammy` on 22.04, `noble` on 24.04:

```bash
dch --distribution "$(lsb_release -cs)" --release ""
```

The committed `debian/changelog` names whichever series it was last edited on,
so **check it against the machine you are building on**. Two quick tells that
you are not on the series the changelog claims: the GDAL version in the CMake
configure block (3.4.1 on jammy, 3.8.4 on noble), and `-D_FORTIFY_SOURCE=3` in
the compile lines, which is 23.10 and later.

The field only appears in the `.changes` file, which matters for an upload to an
archive; it does not affect the `.deb` itself or a repository served with
`dpkg-scanpackages`.

One thing to know about what the binaries report: `CMakeLists.txt` captures
`git describe` at configure time, so `lsdtt-<tool> --version` prints the plain
version only when it is built from a checkout sitting exactly on the matching
`v<version>` tag. Any other build reports `<version>-dev, git <describe>`, and
a build from an exported tarball with no `.git` reports the plain version.
Build releases from a tagged checkout so the two agree.

---

## Building the package

From the root of the lsdtt3 source tree:

```bash
dpkg-buildpackage -us -uc -b -rfakeroot
```

Flags:

- `-us` — do not sign the source package
- `-uc` — do not sign the `.changes` file
- `-b` — build binary package only (no source package)
- `-rfakeroot` — use fakeroot rather than requiring `sudo`

This is a full optimised build of around forty CLI programs; expect it to take
a while. The output files appear one directory above the source tree:

```
../lsdtt3_0.5.1_amd64.deb
../lsdtt3_0.5.1_amd64.buildinfo
../lsdtt3_0.5.1_amd64.changes
```

To build without re-doing everything after a failure, `dpkg-buildpackage -nc`
skips the clean step. To start genuinely fresh, `debian/rules clean` (or
`dh clean`) removes the build tree.

### If the compiler is killed

```
c++: fatal error: Killed signal terminated program cc1plus
```

That is the OOM killer, not an interrupt and not a code problem. debhelper runs
a parallel build sized by `nproc`, and on a many-core machine that can mean
dozens of concurrent `cc1plus` processes. Two translation units are far heavier
than the rest — `lsdtt_chi_ksn_inversion_cpu.cpp` (Eigen template
instantiation) and `lsdtt_raster_coregistration_cpu.cpp` — and Debian's default
flags add LTO with fat objects, which roughly doubles peak memory.

Cap the job count:

```bash
DEB_BUILD_OPTIONS="parallel=4" dpkg-buildpackage -us -uc -b -rfakeroot
```

Halve it until the build completes; `parallel ≈ available_GB / 2` is a
reasonable starting point. Dropping LTO sheds the largest single consumer:

```bash
DEB_BUILD_OPTIONS="parallel=2 optimize=-lto" dpkg-buildpackage -us -uc -b -rfakeroot
```

On a shared login node with a memory cgroup, the better answer is to build on a
compute node through the batch system rather than tuning flags.

### Checking the package with lintian

```bash
lintian ../lsdtt3_0.5.1_amd64.deb
```

Lintian checks the package against Debian policy. A clean run for this package
looks like this:

- **40 × `W: no-manual-page`** — one per CLI, expected. The tools document
  themselves through `--help` and the docs site rather than man pages.
- **Nothing else on the `.deb`.** No `no-copyright-file`, no static-library or
  embedded-library tags, no permission or ownership tags.

Run it on the `.changes` file instead of the `.deb` and you also get the debug
package, which adds **`W: elf-error ... Unable to find program interpreter
name`** once per build-id. That is lintian failing to parse a detached `.debug`
file as a complete ELF executable — it is a known false positive for every
`-dbgsym` package, not a defect in this one.

Anything else should be investigated.

---

## Inspecting the built package

List the files the package will install:

```bash
dpkg -c ../lsdtt3_0.5.1_amd64.deb
```

You should see the `lsdtt-*` programs under `usr/bin/` and
`usr/share/lsdtt3/crn_parameters.nc`, and **no** `usr/lib/*.a` or
`usr/include/`.

Show package metadata, including the dependencies `shlibdeps` worked out:

```bash
dpkg -I ../lsdtt3_0.5.1_amd64.deb
```

Test install locally without a repository:

```bash
sudo dpkg -i ../lsdtt3_0.5.1_amd64.deb
sudo apt-get install -f    # resolves any missing dependencies
lsdtt-channel-extraction --version
```

`--version` prints the GDAL and PROJ versions the binaries actually loaded,
which is the quickest check that the package pulled in a working runtime.

---

## Building for multiple Ubuntu versions

A package built on Ubuntu 22.04 (Jammy) will not install cleanly on 24.04
(Noble) or 20.04 (Focal), because the shared library versions `shlibdeps`
recorded differ. Build a separate `.deb` on each target release.

The easiest way to do that without maintaining several machines is to use the
Apptainer definitions in `apptainer_def/` as a base: copy
`apptainer_def/lsdtt3.def`, change the `From:` line to the target release,
rebuild the container, and build the package inside it.

```bash
apptainer exec lsdtt3.sif bash -c "
    cd /path/to/lsdtt3 &&
    dpkg-buildpackage -us -uc -b -rfakeroot
"
```

Do not use Docker for this.

---

## Hosting the package on your website

You can distribute the `.deb` file directly as a download link, or set up a
minimal APT repository so users can install and update with `apt`.

### Simple APT repository with dpkg-scanpackages

On your web server, create a directory for the packages:

```bash
mkdir -p /var/www/html/lsdtt3/dists/jammy/main/binary-amd64
cp lsdtt3_0.5.1_amd64.deb /var/www/html/lsdtt3/dists/jammy/main/binary-amd64/
```

Generate the `Packages` index:

```bash
cd /var/www/html/lsdtt3
dpkg-scanpackages dists/jammy/main/binary-amd64 /dev/null \
    > dists/jammy/main/binary-amd64/Packages
gzip -k dists/jammy/main/binary-amd64/Packages
```

Users add your repository with:

```bash
echo "deb [trusted=yes] https://yoursite.com/lsdtt3 jammy main" \
    | sudo tee /etc/apt/sources.list.d/lsdtt3.list
sudo apt-get update
sudo apt-get install lsdtt3
```

The `[trusted=yes]` flag bypasses GPG signature checking. For a more robust
setup, sign your repository with a GPG key — see `man apt-secure` for details.
Repeat the `dpkg-scanpackages` step every time you add a new package version.

Serve each Ubuntu release from its own distribution directory (`dists/jammy`,
`dists/noble`, …) built on that release.

---

## See also

- [Building from source](building_from_source.md) — the build this packages
- [Making release packages](making_release_packages.md) — the tarball release route
