#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# ci_ccache.sh
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Run by cibuildwheel before EACH wheel build in CI (CIBW_BEFORE_BUILD_LINUX /
# CIBW_BEFORE_BUILD_MACOS in .github/workflows/python-wheels.yml). Not used by
# local builds or the Apptainer route.
#
# Why: each platform builds five wheels (CPython 3.10-3.14), and all of the
# lsdtt3 runtime library -- about 250 of the ~300 C++ files -- compiles
# identically for every Python version. With ccache as the compiler launcher,
# wheels 2-5 reuse those objects and recompile only the nanobind and binding
# sources. The workflow sets the launcher and the path settings that make the
# commands identical between builds (CCACHE_BASEDIR, a fixed per-wheel
# SKBUILD_BUILD_DIR).
#
# It must never fail a build. If ccache cannot be installed, a pass-through
# stand-in named ccache is put on PATH, so the configured launcher still works
# and simply runs the compiler.

set -u

if ! command -v ccache >/dev/null 2>&1; then
    case "$(uname -s)" in
        Darwin)
            brew install ccache >/dev/null 2>&1 || true
            ;;
        Linux)
            # manylinux_2_28 is AlmaLinux 8: ccache comes from EPEL.
            { dnf -y -q install epel-release && dnf -y -q install ccache; } >/dev/null 2>&1 || true
            ;;
    esac
fi

if ccache --version >/dev/null 2>&1; then
    echo "ci_ccache: $(ccache --version | head -n 1)"
    ccache --show-stats 2>/dev/null || true
else
    # Pass-through stand-in: `ccache <compiler> <args>` -> `<compiler> <args>`.
    if [ "$(uname -s)" = "Darwin" ] && command -v brew >/dev/null 2>&1; then
        bin_dir="$(brew --prefix)/bin"
    else
        bin_dir=/usr/local/bin
    fi
    printf '#!/bin/sh\nexec "$@"\n' > "${bin_dir}/ccache" && chmod +x "${bin_dir}/ccache"
    echo "ci_ccache: ccache unavailable; installed a pass-through stand-in in ${bin_dir}"
fi

exit 0
