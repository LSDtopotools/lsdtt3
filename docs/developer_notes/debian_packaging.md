# Building a Debian package for lsdtt3

This document explains how to build a `.deb` package for distributing lsdtt3
binaries. The package installs the compiled CLI programs and their runtime
dependencies, so end users do not need to build from source.

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
`apptainer_def/apptainer_setup.md`, `building_from_source.md`, or the main README for the full list.

---

## Directory structure

Debian packaging expects a `debian/` directory at the root of the source tree
containing at minimum four files:

```
lsdtt3/
├── debian/
│   ├── control       — package metadata and dependencies
│   ├── rules         — build instructions (essentially a Makefile)
│   ├── changelog     — version history (required by Debian policy)
│   └── compat        — debhelper compatibility level
├── src/
├── CMakeLists.txt
└── ...
```

Create the directory:

```bash
mkdir debian
```

---

## debian/control

This file describes the package, its dependencies, and what it installs.
Adjust the version and description as needed.

```
Source: lsdtt3
Section: science
Priority: optional
Maintainer: Simon M. Mudd <simon.m.mudd@ed.ac.uk>
Build-Depends: debhelper-compat (= 13),
               cmake (>= 3.15),
               libgdal-dev,
               libfftw3-dev,
               libeigen3-dev,
               libgomp1
Standards-Version: 4.6.0
Homepage: https://github.com/LSDtopotools/lsdtt3

Package: lsdtt3
Architecture: any
Depends: ${shlibs:Depends}, ${misc:Depends},
         gdal-bin,
         libgdal30 | libgdal-dev,
         libfftw3-3,
         libgomp1
Description: Topographic analysis tools for digital elevation models
 lsdtt3 is a set of command line tools for topographic analysis of digital
 elevation models (DEMs). It is a redesign of LSDTopoTools2 using modern
 C++17. Tools include raster preprocessing (pit filling, depression
 breaching, nodata repair), D8 flow routing, hillshade, and spectral
 surface analysis.
 .
 See https://github.com/LSDtopotools/lsdtt3 for documentation.
```

**Notes on `Depends`:**

- `${shlibs:Depends}` is automatically expanded by `dpkg-shlibdeps` to the
  correct shared library versions for your build machine.
- `libgdal30` is the runtime library name on Ubuntu 22.04. On Ubuntu 24.04 it
  is `libgdal34`. Check the correct name with `dpkg -l 'libgdal*'` on your
  build machine and adjust accordingly. The `| libgdal-dev` fallback handles
  machines where only the dev package is installed.
- If you build without FFTW support (`LSDTT_HAS_FFTW` not set), remove
  `libfftw3-3` from `Depends`.

---

## debian/rules

The rules file drives the build. This version uses CMake via debhelper's
`dh` command with the `cmake` buildsystem helper.

```makefile
#!/usr/bin/make -f

export DH_VERBOSE = 1

%:
	dh $@ --buildsystem=cmake

override_dh_auto_configure:
	dh_auto_configure -- \
		-DCMAKE_BUILD_TYPE=Release \
		-DLSDTT_GPU=OFF

override_dh_auto_install:
	dh_auto_install
	# The CMakeLists install target puts binaries in bin/ — move them to
	# the standard Debian location /usr/bin
	mkdir -p debian/lsdtt3/usr/bin
	cp -a debian/lsdtt3/usr/local/bin/lsdtt-* debian/lsdtt3/usr/bin/ || true
	rm -rf debian/lsdtt3/usr/local
```

Make the rules file executable:

```bash
chmod +x debian/rules
```

---

## debian/changelog

Debian policy requires a changelog. The format is strict — use `dch` to
create and maintain it rather than editing by hand.

Create the initial entry:

```bash
dch --create --package lsdtt3 --newversion 3.0.0-1 \
    "Initial Debian packaging of lsdtt3."
```

This creates `debian/changelog` with the correct format. Future releases:

```bash
dch --newversion 3.0.1-1 "Release 3.0.1."
dch --release ""
```

The version string has two parts: the upstream version (`3.0.0`) and the
Debian revision (`1`). Increment the Debian revision when you change the
packaging without changing the upstream code; reset it to `1` when you bump
the upstream version.

---

## debian/compat

Set the debhelper compatibility level:

```bash
echo 13 > debian/compat
```

---

## Building the package

From the root of the lsdtt3 source tree:

```bash
dpkg-buildpackage -us -uc -b
```

Flags:
- `-us` — do not sign the source package
- `-uc` — do not sign the `.changes` file
- `-b` — build binary package only (no source package)

The output files appear one directory above the source tree:

```
../lsdtt3_3.0.0-1_amd64.deb
../lsdtt3_3.0.0-1_amd64.buildinfo
../lsdtt3_3.0.0-1_amd64.changes
```

### Building with fakeroot (no sudo required)

```bash
dpkg-buildpackage -us -uc -b -rfakeroot
```

### Checking the package with lintian

```bash
lintian ../lsdtt3_3.0.0-1_amd64.deb
```

Lintian checks the package against Debian policy. Warnings about missing
man pages or a non-official maintainer address are expected for a personal
distribution — errors should be investigated.

---

## Inspecting the built package

List the files the package will install:

```bash
dpkg -c ../lsdtt3_3.0.0-1_amd64.deb
```

Show package metadata:

```bash
dpkg -I ../lsdtt3_3.0.0-1_amd64.deb
```

Test install locally without a repository:

```bash
sudo dpkg -i ../lsdtt3_3.0.0-1_amd64.deb
sudo apt-get install -f    # resolves any missing dependencies
```

---

## Building for multiple Ubuntu versions

The package built on Ubuntu 22.04 (Jammy) will not install cleanly on
Ubuntu 24.04 (Noble) or 20.04 (Focal) because shared library versions differ.
Build a separate `.deb` on each target Ubuntu version.

The easiest way to do this without maintaining multiple machines is to use
the Apptainer container (`lsdtt3.def`) as a base and create variant `.def`
files for each target Ubuntu version — change the `From:` line to
`ubuntu:24.04` etc., rebuild the container, and build the package inside it:

```bash
# Build the package inside the Ubuntu 22.04 container
apptainer exec lsdtt3.sif bash -c "
    cd /path/to/lsdtt3 &&
    dpkg-buildpackage -us -uc -b -rfakeroot
"
```

---

## Hosting the package on your website

You can distribute the `.deb` file directly as a download link, or set up a
minimal APT repository so users can install and update with `apt`.

### Simple APT repository with dpkg-scanpackages

On your web server, create a directory for the packages:

```bash
mkdir -p /var/www/html/lsdtt3/dists/jammy/main/binary-amd64
cp lsdtt3_3.0.0-1_amd64.deb /var/www/html/lsdtt3/dists/jammy/main/binary-amd64/
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
