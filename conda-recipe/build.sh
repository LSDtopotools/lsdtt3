#!/usr/bin/env bash
# conda-build script for lsdtt3: standard CMake build + install into $PREFIX.
set -euo pipefail

cmake -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DLSDTT_OPENCV=ON \
    -DLSDTT_EIGEN=ON

cmake --build build -j"${CPU_COUNT:-2}"

cmake --install build
