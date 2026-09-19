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
│   └── source/format  — source package format
├── src/
├── CMakeLists.txt
└── ...
```

There is deliberately **no `debian/compat`**. The compatibility level is set
once, by `debhelper-compat (= 13)` in `Build-Depends`. Specifying it in both
places is a fatal error — debhelper stops with *"Please specify the compatibility
level exactly once"*.

`debian/rules` must be executable. On WSL over NTFS a `chmod +x` is invisible
to git, so if you have just created or restored the file, tell git about it
explicitly:

```bash
chmod +x debian/rules
git update-index --chmod=+x debian/rules
```

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

Three things are worth knowing about it.

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

**Static libraries and headers are dropped.** `CMakeLists.txt` also installs
`liblsdtt3_core.a`, `liblsdtt3_utils.a`, `liblsdtt3_cpu.a` and the public
headers, which are useful for a development package but dead weight in a
binary one. `override_dh_auto_install` removes them. If you ever want an
`lsdtt3-dev` package, that is where to start.

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
dch --newversion 0.5.2 "Release 0.5.2."
dch --release ""
```

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

### Checking the package with lintian

```bash
lintian ../lsdtt3_0.5.1_amd64.deb
```

Lintian checks the package against Debian policy. Warnings about missing man
pages and a non-official maintainer address are expected for a personal
distribution — errors should be investigated.

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
