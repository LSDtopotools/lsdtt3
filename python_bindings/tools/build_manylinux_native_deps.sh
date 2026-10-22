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
# Why PROJ, GDAL, GEOS and curl are built from source: AlmaLinux 8 / EPEL 8 has
# GDAL 3.0.4 and PROJ 6.3. Every Linux wheel bundles whatever is installed here,
# and PROJ < 7.2 cannot export a CRS that uses a datum ensemble, which is how
# current pyproj and GDAL describe WGS 84 and most other CRSs. The bindings'
# reprojection fails on those ("DatumEnsemble can only be exported to
# WKT2:2019"). GDAL 3.10 needs libcurl >= 7.68 for /vsicurl (the fetch tools)
# and AlmaLinux 8 has 7.61, so curl is built too. GEOS from source keeps it
# current with GDAL.
#
# They are built with the image's own toolchain (gcc-toolset), which keeps the
# result manylinux-compliant; auditwheel repair then bundles them into each
# wheel. Do not use conda-forge builds here: they need a newer libstdc++ than
# manylinux_2_28 allows.
#
# Everything goes to /usr/local (libraries in /usr/local/lib) and is registered
# with ldconfig, so auditwheel and the lsdtt3 build find it. The GDAL build is
# reduced to the drivers lsdtt3 reads and writes; the list is checked at the
# end, so a missing driver fails here rather than in a user's session.
#
# Versions and checksums are pinned. To update one, change the version, run
# `sha256sum` on the new tarball, and update the checksum beside it.
#
# Usage (as root, inside the manylinux_2_28 image):
#   bash python_bindings/tools/build_manylinux_native_deps.sh
# Environment:
#   LSDTT_NATIVE_DEPS_JOBS   parallel build jobs (default: nproc)
# -----------------------------------------------------------------------------

set -euo pipefail

PREFIX=/usr/local
JOBS="${LSDTT_NATIVE_DEPS_JOBS:-$(nproc)}"

CURL_VERSION=8.11.1
CURL_SHA256=c7ca7db48b0909743eaef34250da02c19bc61d4f1dcedd6603f109409536ab56
GEOS_VERSION=3.13.1
GEOS_SHA256=df2c50503295f325e7c8d7b783aca8ba4773919cde984193850cf9e361dfd28c
PROJ_VERSION=9.5.1
PROJ_SHA256=a8395f9696338ffd46b0feb603edbb730fad6746fba77753c77f7f997345e3d3
GDAL_VERSION=3.10.3
GDAL_SHA256=e4bf7f104acbcb3e2d16c97fd1af2b92b28d0ba59d17d976e3ef08b794f4153b

say() { printf '\n== %s ==\n' "$*"; }

say "System packages"
# sqlite provides the sqlite3 program, which the PROJ build uses to make
# proj.db. libtiff-devel is for PROJ's grid support only; GDAL uses its own
# internal libtiff. fftw, netcdf and eigen are lsdtt3's own dependencies.
dnf install -y \
    cmake make gcc-c++ pkgconf-pkg-config xz bzip2 perl \
    sqlite sqlite-devel libtiff-devel openssl-devel zlib-devel \
    fftw-devel netcdf-devel eigen3-devel libgomp

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
fetch "https://curl.se/download/curl-${CURL_VERSION}.tar.xz" "$CURL_SHA256"
fetch "https://download.osgeo.org/geos/geos-${GEOS_VERSION}.tar.bz2" "$GEOS_SHA256"
fetch "https://download.osgeo.org/proj/proj-${PROJ_VERSION}.tar.gz" "$PROJ_SHA256"
fetch "https://github.com/OSGeo/gdal/releases/download/v${GDAL_VERSION}/gdal-${GDAL_VERSION}.tar.gz" "$GDAL_SHA256"

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
    echo "${PREFIX}/lib" > /etc/ld.so.conf.d/lsdtt3-native-deps.conf
    ldconfig
}

say "curl ${CURL_VERSION}"
# Autotools rather than CMake for --enable-versioned-symbols: the system
# libnetcdf (from EPEL) also needs libcurl.so.4 and was linked against
# versioned symbols (CURL_OPENSSL_4). One libcurl.so.4 ends up in the wheel, so
# it has to satisfy both. CA certificates: built-in paths are AlmaLinux's; on
# other distributions lsdtt3/__init__ points CURL_CA_BUNDLE at the system bundle.
tar xJf "curl-${CURL_VERSION}.tar.xz"
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
tar xjf "geos-${GEOS_VERSION}.tar.bz2"
cmake -S "geos-${GEOS_VERSION}" -B build-geos "${CMAKE_COMMON[@]}" \
    -DBUILD_DOCUMENTATION=OFF -DBUILD_GEOSOP=OFF
cmake --build build-geos -j"$JOBS"
cmake --install build-geos
register_libs

say "PROJ ${PROJ_VERSION}"
# No network access (ENABLE_CURL=OFF): the wheel does not download grids.
tar xzf "proj-${PROJ_VERSION}.tar.gz"
cmake -S "proj-${PROJ_VERSION}" -B build-proj "${CMAKE_COMMON[@]}" \
    -DENABLE_CURL=OFF -DENABLE_TIFF=ON -DBUILD_APPS=OFF
cmake --build build-proj -j"$JOBS"
cmake --install build-proj
register_libs

say "GDAL ${GDAL_VERSION}"
# GDAL_USE_EXTERNAL_LIBS=OFF: only the dependencies named here are used, so a
# library that happens to be in the image cannot creep into the wheel. The
# optional drivers are off and the ones lsdtt3 uses are switched back on.
tar xzf "gdal-${GDAL_VERSION}.tar.gz"
cmake -S "gdal-${GDAL_VERSION}" -B build-gdal "${CMAKE_COMMON[@]}" \
    -DBUILD_APPS=ON \
    -DBUILD_PYTHON_BINDINGS=OFF -DBUILD_JAVA_BINDINGS=OFF -DBUILD_CSHARP_BINDINGS=OFF \
    -DGDAL_USE_EXTERNAL_LIBS=OFF \
    -DGDAL_USE_GEOS=ON -DGDAL_USE_SQLITE3=ON -DGDAL_USE_CURL=ON -DGDAL_USE_ZLIB=ON \
    -DCURL_INCLUDE_DIR="${PREFIX}/include" -DCURL_LIBRARY="${PREFIX}/lib/libcurl.so" \
    -DGDAL_USE_TIFF_INTERNAL=ON -DGDAL_USE_GEOTIFF_INTERNAL=ON \
    -DGDAL_USE_JSONC_INTERNAL=ON -DGDAL_USE_PNG_INTERNAL=ON -DGDAL_USE_JPEG_INTERNAL=ON \
    -DGDAL_BUILD_OPTIONAL_DRIVERS=OFF -DOGR_BUILD_OPTIONAL_DRIVERS=OFF \
    -DGDAL_ENABLE_DRIVER_GTIFF=ON -DGDAL_ENABLE_DRIVER_VRT=ON -DGDAL_ENABLE_DRIVER_DERIVED=ON \
    -DGDAL_ENABLE_DRIVER_PNG=ON -DGDAL_ENABLE_DRIVER_JPEG=ON -DGDAL_ENABLE_DRIVER_RAW=ON \
    -DGDAL_ENABLE_DRIVER_AAIGRID=ON -DGDAL_ENABLE_DRIVER_HFA=ON \
    -DOGR_ENABLE_DRIVER_SHAPE=ON -DOGR_ENABLE_DRIVER_GEOJSON=ON -DOGR_ENABLE_DRIVER_CSV=ON \
    -DOGR_ENABLE_DRIVER_FLATGEOBUF=ON -DOGR_ENABLE_DRIVER_SQLITE=ON -DOGR_ENABLE_DRIVER_GPKG=ON \
    -DOGR_ENABLE_DRIVER_MVT=ON -DOGR_ENABLE_DRIVER_VRT=ON -DOGR_ENABLE_DRIVER_KML=ON
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

gdal_lib="$(readlink -f "${PREFIX}/lib/libgdal.so")"
ldd "$gdal_lib" | grep -q "${PREFIX}/lib/libcurl.so" || fail "libgdal does not link the libcurl built here"
ldd "$gdal_lib" | grep -q "${PREFIX}/lib/libproj.so" || fail "libgdal does not link the PROJ built here"
ldd "$gdal_lib" | grep -q "libgeos_c" || fail "libgdal was built without GEOS"
echo "native dependencies ready in ${PREFIX}"
