#!/usr/bin/env bash
# conda-build script for lsdtt3: standard CMake build + install into $PREFIX.
set -euo pipefail

# ${CMAKE_ARGS} carries conda-build's toolchain and prefix settings (and the
# cross-compilation setup on conda-forge); it must be passed through.
# BUILD_TESTING=OFF: the CTest suite is not packaged and roughly doubles the
# build time. GPU stays at its default (OFF).
cmake ${CMAKE_ARGS} -S . -B build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DBUILD_TESTING=OFF \
    -DLSDTT_OPENMP=ON \
    -DLSDTT_OPENCV=ON \
    -DLSDTT_EIGEN=ON

cmake --build build -j"${CPU_COUNT:-2}"

cmake --install build

# The install rules also emit the internal static libraries and a partial set
# of headers (core/ and utils/ only, not cpu/). They are not a usable C++ SDK,
# and conda-forge does not accept static libraries in a package, so ship the
# command-line tools and their runtime data only.
rm -f "${PREFIX}"/lib/liblsdtt3_*.a
rm -rf "${PREFIX}/include/lsdtt3"
