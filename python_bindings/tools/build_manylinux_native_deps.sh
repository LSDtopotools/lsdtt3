#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# -----------------------------------------------------------------------------
# build_manylinux_native_deps.sh
#
# Installs the native libraries the Linux wheels link against, inside the
# manylinux_2_28 image. This is the [tool.cibuildwheel.linux] before-all step,
# so the same script runs for CI (python-wheels.yml), for podman via
# tools/make_all_packages.sh, and for the Apptainer route
# (apptainer_def/manylinux_wheel.def).
#
# Why PROJ, GDAL, GEOS, SQLite and curl are built from source: AlmaLinux 8 / EPEL 8 has
# GDAL 3.0.4 and PROJ 6.3. Every Linux wheel bundles whatever is installed here,
# and PROJ < 7.2 cannot export a CRS that uses a datum ensemble, which is how
# current pyproj and GDAL describe WGS 84 and most other CRSs. The bindings'
# reprojection fails on those ("DatumEnsemble can only be exported to
# WKT2:2019"). GDAL 3.10 also needs newer versions than AlmaLinux 8 has of
# libcurl (>= 7.68 for /vsicurl, used by the fetch tools; 8 has 7.61) and
# SQLite (>= 3.31; 8 has 3.26), and GEOS >= 3.8, so those are built too.
#
# They are built with the image's own toolchain (gcc-toolset), which keeps the
# result manylinux-compliant; auditwheel repair then bundles them into each
# wheel. Do not use conda-forge builds here: they need a newer libstdc++ than
# manylinux_2_28 allows.
#
# Everything goes to /usr/local (libraries in /usr/local/lib) and is registered
# with ldconfig, so auditwheel and the lsdtt3 build find it. The GDAL build is
# checked at the end for the drivers lsdtt3 reads and writes, so a missing
# driver fails here rather than in a user's session.
#
# Versions and checksums are pinned. To update one, change the version, run
# `sha256sum` on the new tarball, and update the checksum beside it.
#
# Usage (as root, inside the manylinux_2_28 image):
#   bash python_bindings/tools/build_manylinux_native_deps.sh
# Environment:
#   LSDTT_NATIVE_DEPS_JOBS     parallel build jobs (default: nproc)
#   LSDTT_NATIVE_DEPS_PREFIX   install prefix (default: /usr/local)
#
# Run as an ordinary user (not root) it is a local test: it skips the dnf
# install and ldconfig, uses LD_LIBRARY_PATH instead, and needs
# LSDTT_NATIVE_DEPS_PREFIX set to a writable directory. That checks the
# downloads, the configure options, the builds and the closing checks without
# the image.
# -----------------------------------------------------------------------------

set -euo pipefail

PREFIX="${LSDTT_NATIVE_DEPS_PREFIX:-/usr/local}"
JOBS="${LSDTT_NATIVE_DEPS_JOBS:-$(nproc)}"
AS_ROOT=0
[[ "$(id -u)" -eq 0 ]] && AS_ROOT=1
if (( ! AS_ROOT )) && [[ -z "${LSDTT_NATIVE_DEPS_PREFIX:-}" ]]; then
    echo "ERROR: not root; for a local test set LSDTT_NATIVE_DEPS_PREFIX to a writable directory" >&2
    exit 1
fi

SQLITE_VERSION=3470200   # 3.47.2, in SQLite's tarball numbering
SQLITE_YEAR=2024
SQLITE_SHA256=f1b2ee412c28d7472bc95ba996368d6f0cdcf00362affdadb27ed286c179540b
CURL_VERSION=8.11.1
CURL_SHA256=c7ca7db48b0909743eaef34250da02c19bc61d4f1dcedd6603f109409536ab56
GEOS_VERSION=3.13.1
GEOS_SHA256=df2c50503295f325e7c8d7b783aca8ba4773919cde984193850cf9e361dfd28c
PROJ_VERSION=9.5.1
PROJ_SHA256=a8395f9696338ffd46b0feb603edbb730fad6746fba77753c77f7f997345e3d3
GDAL_VERSION=3.10.3
GDAL_SHA256=e4bf7f104acbcb3e2d16c97fd1af2b92b28d0ba59d17d976e3ef08b794f4153b

say() { printf '\n== %s ==\n' "$*"; }

# libtiff-devel is for PROJ's grid support only; GDAL uses its own internal
# libtiff. No sqlite-devel: SQLite is built below, and leaving the old headers
# out means nothing can pick them up by mistake. fftw, netcdf and eigen are
# lsdtt3's own dependencies.
if (( AS_ROOT )); then
    say "System packages"
    dnf install -y \
        cmake make gcc-c++ pkgconf-pkg-config xz bzip2 perl \
        libtiff-devel openssl-devel zlib-devel \
        fftw-devel netcdf-devel eigen3-devel libgomp
fi
mkdir -p "$PREFIX"

# Use the image's gcc-toolset (the newest present), as cibuildwheel's build
# does. Its PATH entry comes from the image's environment, which an Apptainer
# %post may not have; without it the build would fall back to AlmaLinux's gcc 8,
# which is marginal for GDAL 3.10 and PROJ 9.5 (C++17). The enable scripts
# read unset variables, hence set +u around them.
toolset_enable="$(ls -d /opt/rh/gcc-toolset-*/enable 2>/dev/null | sort -V | tail -1 || true)"
if [[ -n "$toolset_enable" ]]; then
    set +u
    # shellcheck disable=SC1090
    source "$toolset_enable"
    set -u
fi
echo "compiler: $(command -v g++) ($(g++ -dumpfullversion))"

WORK="$(mktemp -d /tmp/lsdtt3-native-deps-XXXXXX)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

# fetch <url> <sha256>: download into $WORK and check the checksum.
fetch() {
    local url="$1" sha="$2" file
    file="$(basename "$url")"
    curl -fsSL --retry 3 -o "$file" "$url"
    echo "${sha}  ${file}" | sha256sum -c -
}

# Download everything first. Once the curl built below is registered with
# ldconfig, the image's own /usr/bin/curl (7.61) loads it and fails ("unknown
# option": it asks for features this build leaves out). Nothing after that step
# may use the system curl, or dnf, which also uses libcurl.
say "Downloading sources"
fetch "https://www.sqlite.org/${SQLITE_YEAR}/sqlite-autoconf-${SQLITE_VERSION}.tar.gz" "$SQLITE_SHA256"
fetch "https://curl.se/download/curl-${CURL_VERSION}.tar.xz" "$CURL_SHA256"
fetch "https://download.osgeo.org/geos/geos-${GEOS_VERSION}.tar.bz2" "$GEOS_SHA256"
fetch "https://download.osgeo.org/proj/proj-${PROJ_VERSION}.tar.gz" "$PROJ_SHA256"
fetch "https://github.com/OSGeo/gdal/releases/download/v${GDAL_VERSION}/gdal-${GDAL_VERSION}.tar.gz" "$GDAL_SHA256"

# tar --no-same-owner below: as root, tar restores each file's owner from
# the archive (uid 1001 in GEOS's). Without /etc/subuid entries Apptainer runs
# this in a namespace that maps only one uid, so that chown fails.

# The libraries all go in $PREFIX/lib (not lib64) so one ldconfig entry covers
# them, and so CMake finds them before any older system copy.
CMAKE_COMMON=(
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX="$PREFIX"
    -DCMAKE_INSTALL_LIBDIR=lib
    -DCMAKE_PREFIX_PATH="$PREFIX"
    -DBUILD_SHARED_LIBS=ON
    -DBUILD_TESTING=OFF
)
register_libs() {
    if (( AS_ROOT )); then
        echo "${PREFIX}/lib" > /etc/ld.so.conf.d/lsdtt3-native-deps.conf
        ldconfig
    else
        export LD_LIBRARY_PATH="${PREFIX}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    fi
}
# Each tool below must use the SQLite built here, not a system copy.
SQLITE_ARGS=(
    -DSQLite3_INCLUDE_DIR="${PREFIX}/include"
    -DSQLite3_LIBRARY="${PREFIX}/lib/libsqlite3.so"
)

say "SQLite ${SQLITE_VERSION}"
# COLUMN_METADATA is required by GDAL; RTREE is used by GeoPackage spatial
# indexes. The sqlite3 program built here is the one PROJ uses to make proj.db.
tar --no-same-owner -xzf "sqlite-autoconf-${SQLITE_VERSION}.tar.gz"
(
    cd "sqlite-autoconf-${SQLITE_VERSION}"
    CFLAGS="-O2 -DSQLITE_ENABLE_COLUMN_METADATA=1 -DSQLITE_ENABLE_RTREE=1" \
        ./configure --prefix="$PREFIX" --libdir="${PREFIX}/lib" --disable-static
    make -j"$JOBS"
    make install
)
register_libs

say "curl ${CURL_VERSION}"
# Autotools rather than CMake for --enable-versioned-symbols: the system
# libnetcdf (from EPEL) also needs libcurl.so.4 and was linked against
# versioned symbols (CURL_OPENSSL_4). One libcurl.so.4 ends up in the wheel, so
# it has to satisfy both. CA certificates: built-in paths are AlmaLinux's; on
# other distributions lsdtt3/__init__ points CURL_CA_BUNDLE at the system bundle.
tar --no-same-owner -xJf "curl-${CURL_VERSION}.tar.xz"
(
    cd "curl-${CURL_VERSION}"
    ./configure --prefix="$PREFIX" --libdir="${PREFIX}/lib" \
        --with-openssl --enable-versioned-symbols --with-ca-fallback \
        --disable-static --disable-manual --disable-docs --disable-ldap \
        --without-libpsl --without-libidn2 --without-nghttp2 \
        --without-brotli --without-zstd --without-libssh2
    make -j"$JOBS"
    make install
)
register_libs

say "GEOS ${GEOS_VERSION}"
tar --no-same-owner -xjf "geos-${GEOS_VERSION}.tar.bz2"
cmake -S "geos-${GEOS_VERSION}" -B build-geos "${CMAKE_COMMON[@]}" \
    -DBUILD_DOCUMENTATION=OFF -DBUILD_GEOSOP=OFF
cmake --build build-geos -j"$JOBS"
cmake --install build-geos
register_libs

say "PROJ ${PROJ_VERSION}"
# No network access (ENABLE_CURL=OFF): the wheel does not download grids.
tar --no-same-owner -xzf "proj-${PROJ_VERSION}.tar.gz"
cmake -S "proj-${PROJ_VERSION}" -B build-proj "${CMAKE_COMMON[@]}" "${SQLITE_ARGS[@]}" \
    -DEXE_SQLITE3="${PREFIX}/bin/sqlite3" \
    -DENABLE_CURL=OFF -DENABLE_TIFF=ON -DBUILD_APPS=OFF
cmake --build build-proj -j"$JOBS"
cmake --install build-proj
register_libs

say "GDAL ${GDAL_VERSION}"
# GDAL_USE_EXTERNAL_LIBS=OFF: only the dependencies named here are used, so a
# library that happens to be in the image cannot creep into the wheel. Every
# driver those dependencies (and GDAL's internal copies of libtiff, libgeotiff,
# libpng, libjpeg, json-c and so on) allow is built, so users can open the
# same wide range of formats as with any GDAL; drivers needing another
# library (netCDF, HDF5, JPEG2000, WebP, Parquet ...) are not.
# Iconv is glibc's (no extra library); Shapefile reading uses it to recode
# attribute text from the code page in the .cpg file.
tar --no-same-owner -xzf "gdal-${GDAL_VERSION}.tar.gz"
cmake -S "gdal-${GDAL_VERSION}" -B build-gdal "${CMAKE_COMMON[@]}" "${SQLITE_ARGS[@]}" \
    -DBUILD_APPS=ON \
    -DBUILD_PYTHON_BINDINGS=OFF -DBUILD_JAVA_BINDINGS=OFF -DBUILD_CSHARP_BINDINGS=OFF \
    -DGDAL_USE_EXTERNAL_LIBS=OFF \
    -DGDAL_USE_GEOS=ON -DGDAL_USE_SQLITE3=ON -DGDAL_USE_CURL=ON -DGDAL_USE_ZLIB=ON \
    -DGDAL_USE_ICONV=ON \
    -DCURL_INCLUDE_DIR="${PREFIX}/include" -DCURL_LIBRARY="${PREFIX}/lib/libcurl.so" \
    -DGDAL_USE_TIFF_INTERNAL=ON -DGDAL_USE_GEOTIFF_INTERNAL=ON \
    -DGDAL_USE_JSONC_INTERNAL=ON -DGDAL_USE_PNG_INTERNAL=ON -DGDAL_USE_JPEG_INTERNAL=ON \
    -DGDAL_BUILD_OPTIONAL_DRIVERS=ON -DOGR_BUILD_OPTIONAL_DRIVERS=ON
cmake --build build-gdal -j"$JOBS"
cmake --install build-gdal
register_libs

say "Checking the result"
fail() { echo "ERROR: $*" >&2; exit 1; }
export PATH="${PREFIX}/bin:${PATH}"
# AlmaLinux pkgconf does not search /usr/local by default.
export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
echo "GDAL $(gdal-config --version), PROJ $(pkg-config --modversion proj), GEOS $(geos-config --version)"
[[ "$(gdal-config --version)" == "$GDAL_VERSION" ]] || fail "gdal-config reports $(gdal-config --version), expected ${GDAL_VERSION}"
[[ -f "${PREFIX}/share/proj/proj.db" ]] || fail "no ${PREFIX}/share/proj/proj.db"
[[ -f "$(gdal-config --datadir)/gdalvrt.xsd" ]] || fail "GDAL data directory incomplete: $(gdal-config --datadir)"

raster_formats="$(gdalinfo --formats)"
for d in GTiff COG MEM VRT PNG JPEG ENVI EHdr AAIGrid HFA; do
    grep -qE "^ *${d} -" <<< "$raster_formats" || fail "GDAL raster driver ${d} was not built"
done
vector_formats="$(ogrinfo --formats)"
for d in "ESRI Shapefile" GeoJSON GPKG FlatGeobuf MVT CSV; do
    grep -qE "^ *${d} -" <<< "$vector_formats" || fail "OGR vector driver ${d} was not built"
done

# Captured first: with pipefail, `ldd | grep -q` can fail when grep exits
# early and ldd gets SIGPIPE.
gdal_ldd="$(ldd "$(readlink -f "${PREFIX}/lib/libgdal.so")")"
proj_ldd="$(ldd "$(readlink -f "${PREFIX}/lib/libproj.so")")"
grep -q "${PREFIX}/lib/libcurl.so" <<< "$gdal_ldd" || fail "libgdal does not link the libcurl built here"
grep -q "${PREFIX}/lib/libproj.so" <<< "$gdal_ldd" || fail "libgdal does not link the PROJ built here"
grep -q "${PREFIX}/lib/libgeos_c" <<< "$gdal_ldd" || fail "libgdal does not link the GEOS built here"
grep -q "${PREFIX}/lib/libsqlite3.so" <<< "$gdal_ldd" || fail "libgdal does not link the SQLite built here"
grep -q "${PREFIX}/lib/libsqlite3.so" <<< "$proj_ldd" || fail "libproj does not link the SQLite built here"
echo "native dependencies ready in ${PREFIX}"
