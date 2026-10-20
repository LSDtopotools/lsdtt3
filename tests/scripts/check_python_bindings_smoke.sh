#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# check_python_bindings_smoke.sh - build the Python bindings once on this
# machine and run the whole python_bindings/tests suite against them.
#
# The full wheel matrix (python-wheels.yml: 4 platforms x 5 Pythons, every
# wheel compiled from scratch) is run by hand and only for some releases. This
# is the cheap check in between: one Linux build, one interpreter, CLIs off.
# It catches a C++ change that breaks the bindings (a compile break, or a
# behaviour change a binding test sees) at the commit that made it.
#
# Two callers:
#   - CTest test python_bindings_smoke (labels python, heavy; opt-in);
#   - tools/make_release.sh, as a mandatory step, so every v* tag runs it.
#
# It builds through pip (scikit-build-core), not direct CMake, because:
#   - lsdtt3.__version__ comes from installed distribution metadata, so only an
#     installed package can check _core.lsdtt3_version() against the real
#     version;
#   - SKBUILD_PROJECT_VERSION, the source of the C++ LSDTT3_VERSION in a
#     wheel, exists only under scikit-build-core.
# The private keyword-compatibility probe is built and required, so
# test_keyword_compatibility.py does not skip its probe test.
#
# Needs: a Python >= 3.10 with the venv module, a C++ toolchain, GDAL
# development files, and network access (or a pip cache) for nanobind,
# scikit-build-core, numpy, geopandas, shapely, pyproj and pytest.
#
# Usage:
#   check_python_bindings_smoke.sh --repo-root <dir> --work-dir <dir>
#                                  [--python <exe>] [--reuse]
#
#   --work-dir  holds the venv and the persistent CMake build directory, so a
#               second run recompiles only what changed.
#   --reuse     keep an existing venv instead of recreating it. pip still
#               rebuilds (incrementally) and the tests always run.
#
# Exit codes: 0 pass, 1 fail, 77 skipped (no Python with venv available).

set -euo pipefail

REPO_ROOT=""
WORK_DIR=""
PYTHON="python3"
REUSE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --repo-root) REPO_ROOT="${2:-}"; shift 2 ;;
        --work-dir)  WORK_DIR="${2:-}"; shift 2 ;;
        --python)    PYTHON="${2:-}"; shift 2 ;;
        --reuse)     REUSE=1; shift ;;
        -h|--help)   sed -n '2,/^set -euo/p' "$0" | sed 's/^# \{0,1\}//;/^set -euo/d'; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; exit 1 ;;
    esac
done
[[ -n "$REPO_ROOT" && -n "$WORK_DIR" ]] || {
    echo "ERROR: --repo-root and --work-dir are required" >&2; exit 1; }

REPO_ROOT="$(cd "$REPO_ROOT" && pwd)"
mkdir -p "$WORK_DIR"
WORK_DIR="$(cd "$WORK_DIR" && pwd)"
VENV="$WORK_DIR/venv"

# A pixi or conda shell exports GDAL/PROJ data paths for ITS libraries. The
# bindings link whichever GDAL CMake finds, so inherited paths can point a
# build at another installation's proj.db. Start clean.
unset GDAL_DATA GDAL_DRIVER_PATH PROJ_DATA PROJ_LIB

if ! "$PYTHON" -c "import venv, ensurepip" >/dev/null 2>&1; then
    echo "SKIP: $PYTHON has no working venv module (on Debian/Ubuntu: apt install python3-venv)." >&2
    exit 77
fi

start=$(date +%s)

if [[ $REUSE -eq 0 || ! -x "$VENV/bin/python" ]]; then
    rm -rf "$VENV"
    "$PYTHON" -m venv "$VENV"
fi
"$VENV/bin/python" -m pip install --quiet --upgrade pip

echo "== Building the bindings (CLIs off, compatibility probe on) =="
"$VENV/bin/python" -m pip install --quiet \
    -C cmake.define.LSDTT_PYTHON_CLI=OFF \
    -C cmake.define.LSDTT_PYTHON_COMPATIBILITY_PROBE=ON \
    -C "build-dir=$WORK_DIR/build" \
    "$REPO_ROOT/python_bindings" pytest
built=$(date +%s)

echo "== Running python_bindings/tests =="
# Run from outside the source tree: python_bindings/src/lsdtt3 must not
# shadow the installed package.
cd "$WORK_DIR"
status=0
LSDTT_REQUIRE_COMPATIBILITY_PROBE=1 \
    "$VENV/bin/python" -m pytest -q -p no:cacheprovider \
    "$REPO_ROOT/python_bindings/tests" || status=$?
finished=$(date +%s)

echo "== Python bindings smoke check: build $((built - start)) s, tests $((finished - built)) s =="
if [[ $status -ne 0 ]]; then
    echo "FAIL: python_bindings/tests failed against this build (pytest exit $status)." >&2
    exit 1
fi
echo "PASS: bindings build and the whole python_bindings/tests suite passes."
