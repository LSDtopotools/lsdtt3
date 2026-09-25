#!/bin/sh
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Packaging smoke checks for the lsdtt3 Debian package.
#
# This is NOT the test suite. The CTest suite is disabled during packaging
# (debian/rules passes -DBUILD_TESTING=OFF) because a package build is the
# wrong place to discover a test regression. What this script covers instead
# is the class of failure that unit tests cannot see: a package that compiled
# cleanly but produces binaries that will not run on the target machine, or
# that is quietly missing a tool or a data file because an optional build
# dependency was absent.
#
# Two modes:
#   built  <build-dir>   run after dh_auto_build, on the freshly built tree
#   staged <destdir>     run after dh_auto_install, on the staged package tree
#
# Any failure exits non-zero and aborts the package build.

set -eu

fail() {
    echo "smoke-test: FAIL: $*" >&2
    exit 1
}

note() {
    echo "smoke-test: $*"
}

# True for an ELF binary. The build tree also holds generated help files named
# lsdtt-*.csv / lsdtt-*.html, and a source tree on an NTFS mount (WSL) reports
# every file as executable, so a mode test alone does not distinguish them.
is_elf() {
    [ -f "$1" ] || return 1
    magic=$(head -c 4 "$1" 2>/dev/null | od -An -tx1 | tr -d ' \n')
    [ "$magic" = "7f454c46" ]
}

# Minimum number of lsdtt-* programs expected. A real drop below this means a
# group of tools silently failed to build. Kept well under the true count so
# that adding or removing a single tool does not need a packaging change.
MIN_TOOLS=40

check_built() {
    build_dir=$1
    [ -d "$build_dir" ] || fail "build directory '$build_dir' does not exist"

    tools=""
    for candidate in "$build_dir"/lsdtt-*; do
        is_elf "$candidate" || continue
        tools="$tools$candidate
"
    done
    [ -n "$tools" ] || fail "no lsdtt-* programs found in '$build_dir'"

    count=$(echo "$tools" | grep -c .)
    [ "$count" -ge "$MIN_TOOLS" ] ||
        fail "only $count lsdtt-* programs built, expected at least $MIN_TOOLS"

    # Every program must start, load its shared libraries, initialise GDAL and
    # exit cleanly. --version is the cheapest path that does all of that: it
    # reports the GDAL and PROJ actually loaded, so a binary linked against one
    # GDAL and running against another shows up here rather than in a user's
    # analysis.
    for tool in $tools; do
        out=$("$tool" --version 2>&1) ||
            fail "$(basename "$tool") --version exited non-zero:
$out"
        case $out in
            *GDAL*) ;;
            *) fail "$(basename "$tool") --version did not report a GDAL version:
$out" ;;
        esac
    done
    note "$count programs started, reported GDAL/PROJ and exited 0"

    # Optional dependencies degrade silently in this build: a missing library
    # is a CMake STATUS or WARNING line, never an error, so the package builds
    # and ships without the feature. Assert the two that are visible from
    # outside the binary.

    # FFTW3: lsdtt-spectral-analysis is only added as a target when FFTW3 is
    # found, so its absence means the package is a tool short.
    is_elf "$build_dir/lsdtt-spectral-analysis" ||
        fail "lsdtt-spectral-analysis was not built - FFTW3 was not found at
configure time (is pkg-config installed alongside libfftw3-dev?)"

    # netCDF: without it the cosmo tools compile, run, and quietly fall back to
    # the standard atmosphere instead of ERA40 / palaeomagnetic corrections.
    # The library is linked, not dlopened, so ldd is a reliable witness.
    if command -v ldd >/dev/null 2>&1; then
        ldd "$build_dir/lsdtt-cosmo-tool" 2>/dev/null | grep -q libnetcdf ||
            fail "lsdtt-cosmo-tool is not linked against libnetcdf - ERA40 and
palaeomagnetic corrections were compiled out (is libnetcdf-dev installed?)"
    fi
    note "optional dependencies present: FFTW3 (spectral analysis), netCDF (cosmo corrections)"

    # A package must link only libraries that exist on a plain system. On a
    # machine without root the build may run inside conda or another private
    # prefix, and if CMake resolves GDAL, FFTW or netCDF there instead of in
    # /usr, the .deb depends on libraries no other machine has - it installs
    # cleanly and then fails at run time. dpkg-shlibdeps usually catches this,
    # but it reports it as an unmappable path rather than as the mistake it is.
    #
    # Eigen is exempt by nature: header-only, so it never appears here however
    # it was found.
    if command -v ldd >/dev/null 2>&1; then
        for tool in $tools; do
            private=$(ldd "$tool" 2>/dev/null |
                grep -oE '/[^ ]*(conda|/envs/|virtualenv)[^ ]*' || true)
            [ -z "$private" ] && continue
            fail "$(basename "$tool") links libraries from a private prefix:
$private
A package built this way works only where that environment exists. Configure so
that GDAL, FFTW and netCDF resolve in /usr - do not put the whole prefix on
CMAKE_PREFIX_PATH (see the LSDTT_EIGEN3_DIR note in debian/rules)."
        done
        note "no binary links against a conda or virtualenv prefix"
    fi
}

check_staged() {
    destdir=$1
    [ -d "$destdir" ] || fail "staged directory '$destdir' does not exist"

    # The CLIs must be in /usr/bin. Nothing in this build installs to
    # /usr/local; if that ever changes, the bundled data file below moves with
    # it and the cosmo tools stop finding it.
    count=0
    for candidate in "$destdir"/usr/bin/lsdtt-*; do
        is_elf "$candidate" && count=$((count + 1))
    done
    [ "$count" -ge "$MIN_TOOLS" ] ||
        fail "only $count lsdtt-* programs staged in usr/bin, expected at least $MIN_TOOLS"

    # The cosmo tools resolve this as <exe>/../share/lsdtt3/crn_parameters.nc,
    # so /usr/bin + /usr/share/lsdtt3 is the layout they need. Without it they
    # cannot run without a download.
    crn="$destdir/usr/share/lsdtt3/crn_parameters.nc"
    [ -s "$crn" ] ||
        fail "bundled crn_parameters.nc missing from usr/share/lsdtt3 - the
cosmo tools will not find it at <exe>/../share/lsdtt3"

    # Static libraries and headers are of no use in a runtime package and
    # debian/rules removes them. Catch a regression in that removal.
    if find "$destdir/usr" -name 'liblsdtt3_*.a' 2>/dev/null | grep -q .; then
        fail "static libraries are still staged - they belong in a -dev package"
    fi
    [ ! -d "$destdir/usr/include/lsdtt3" ] ||
        fail "headers are still staged - they belong in a -dev package"

    note "$count programs in usr/bin, crn_parameters.nc present, no static libs or headers"
}

mode=${1:-}
case $mode in
    built)  [ $# -eq 2 ] || fail "usage: $0 built <build-dir>";  check_built  "$2" ;;
    staged) [ $# -eq 2 ] || fail "usage: $0 staged <destdir>"; check_staged "$2" ;;
    *) fail "usage: $0 built <build-dir> | staged <destdir>" ;;
esac

note "$mode checks passed"
