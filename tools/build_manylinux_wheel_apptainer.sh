#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# -----------------------------------------------------------------------------
# build_manylinux_wheel_apptainer.sh
#
# Build, repair and test the lsdtt3 manylinux Python wheels with Apptainer,
# without podman, Docker or sudo.
#
# cibuildwheel can only drive Docker or podman. On Linux all it does is run
# three steps inside the manylinux_2_28 image, and that image already holds
# every CPython and auditwheel, so this script runs the same steps under
# Apptainer:
#
#   1. build   pip wheel python_bindings   (CMAKE_ARGS from pyproject.toml)
#   2. repair  auditwheel repair           (vendors GDAL, PROJ, FFTW, netCDF)
#   3. test    a clean venv with only the repaired wheel, then the
#              [tool.cibuildwheel] test-command from pyproject.toml
#
# The build matrix, the system dependencies (before-all), the environment and
# the test command are all read from python_bindings/pyproject.toml, so this
# route and the cibuildwheel route (CI, or podman via make_all_packages.sh)
# build the same thing.
#
# The image is built once from apptainer_def/manylinux_wheel.def and reused.
# It is named after a hash of the .def file and the before-all command, so it
# is rebuilt only when one of those changes. Building it needs
# `apptainer build --fakeroot`, which Apptainer >= 1.1 allows an ordinary user
# (no sudo). It takes a few GB in the image directory and Apptainer's cache;
# `apptainer cache clean` empties the cache once the image exists.
#
# tools/make_all_packages.sh calls this script when podman is absent and
# Apptainer is present. It can also be run directly; see --help.
# -----------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEF_FILE="${REPO_ROOT}/apptainer_def/manylinux_wheel.def"

if [[ -t 1 ]]; then
    C_RED=$'\033[31m'; C_GREEN=$'\033[32m'; C_YELLOW=$'\033[33m'
    C_BOLD=$'\033[1m'; C_RESET=$'\033[0m'
else
    C_RED=""; C_GREEN=""; C_YELLOW=""; C_BOLD=""; C_RESET=""
fi
info() { printf '  %s\n' "$*"; }
ok()   { printf '  %s%s%s\n' "$C_GREEN" "$*" "$C_RESET"; }
warn() { printf '  %s%s%s\n' "$C_YELLOW" "$*" "$C_RESET" >&2; }
fail() { printf '%sERROR:%s %s\n' "$C_RED" "$C_RESET" "$*" >&2; exit 1; }
step() { printf '\n%s== %s ==%s\n' "$C_BOLD" "$*" "$C_RESET"; }

usage() {
    cat <<'EOF'
Usage:
  tools/build_manylinux_wheel_apptainer.sh [options]

Builds manylinux_2_28 x86_64 wheels of lsdtt3 inside Apptainer, repairs them
with auditwheel and tests each one, the same way cibuildwheel does. Needs only
apptainer on the host: no podman, Docker, sudo or host Python.

Options:
  --source <dir>     Source tree to build (must contain python_bindings/).
                     Default: this repository.
  --out-dir <dir>    Where the repaired wheels go. Default: <source>/wheelhouse
  --python <list>    Comma-separated Python versions, e.g. 3.12 or 3.12,3.13,
                     or "all" for every version in [tool.cibuildwheel] build.
                     The first builds the abi3 wheel; the others test it.
                     Default: all.
  --sif <path>       The image to use, built there if missing.
                     Default: $LSDTT_WHEEL_SIF, else
                     ${XDG_CACHE_HOME:-~/.cache}/lsdtt3/manylinux_2_28-lsdtt3-<hash>.sif
  --rebuild-image    Rebuild the image even if it exists.
  --image-only       Build (or check) the image, then stop.
  --no-test          Skip the test step. Do not publish untested wheels.
  --jobs <N>         Parallel compile jobs. Default: nproc.
  --check-only       Check tooling and print the plan, then exit.
  --help             Show this help.

Examples:
  tools/build_manylinux_wheel_apptainer.sh --check-only
  tools/build_manylinux_wheel_apptainer.sh --image-only
  tools/build_manylinux_wheel_apptainer.sh --python 3.12
  tools/build_manylinux_wheel_apptainer.sh            # all five versions

For a release, build from a tag with
  tools/make_all_packages.sh --clean --tag vX.Y.Z --only wheel --wheel-python all
which stamps the version and calls this script when podman is absent.
EOF
}

SOURCE="$REPO_ROOT"
OUT_DIR=""
PYTHONS="all"
SIF="${LSDTT_WHEEL_SIF:-}"
REBUILD_IMAGE=0
IMAGE_ONLY=0
RUN_TESTS=1
JOBS=""
CHECK_ONLY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --source)        SOURCE="${2:-}"; shift 2 ;;
        --out-dir)       OUT_DIR="${2:-}"; shift 2 ;;
        --python)        PYTHONS="${2:-}"; shift 2 ;;
        --sif)           SIF="${2:-}"; shift 2 ;;
        --rebuild-image) REBUILD_IMAGE=1; shift ;;
        --image-only)    IMAGE_ONLY=1; shift ;;
        --no-test)       RUN_TESTS=0; shift ;;
        --jobs)          JOBS="${2:-}"; shift 2 ;;
        --check-only)    CHECK_ONLY=1; shift ;;
        --help|-h)       usage; exit 0 ;;
        *)               fail "unknown option: $1 (try --help)" ;;
    esac
done

JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
SOURCE="$(cd "$SOURCE" 2>/dev/null && pwd)" || fail "--source directory does not exist"
PYPROJECT="${SOURCE}/python_bindings/pyproject.toml"
DEPS_SCRIPT="${SOURCE}/python_bindings/tools/build_manylinux_native_deps.sh"
OUT_DIR="${OUT_DIR:-${SOURCE}/wheelhouse}"

# ----------------------------------------------------------------------------
# Prerequisites
# ----------------------------------------------------------------------------
command -v apptainer >/dev/null 2>&1 || fail "apptainer is not installed (or not on PATH)"
[[ -f "$DEF_FILE" ]]  || fail "image definition not found: ${DEF_FILE}"
[[ -f "$PYPROJECT" ]] || fail "no python_bindings/pyproject.toml under ${SOURCE}"
[[ -f "$DEPS_SCRIPT" ]] || fail "no python_bindings/tools/build_manylinux_native_deps.sh under ${SOURCE}"

# toml_value <section> <key>: the raw right-hand side of a one-line key in a
# pyproject section. Only used for the image hash and the version list, so
# the host needs no Python; inside the image tomllib reads the file properly.
toml_value() {
    awk -v want="[$1]" -v key="$2" '
        /^\[/ { section = $0; next }
        section == want && $0 ~ "^" key "[ \t]*=" {
            sub("^" key "[ \t]*=[ \t]*", ""); print; exit
        }' "$PYPROJECT"
}

BEFORE_ALL="$(toml_value tool.cibuildwheel.linux before-all)"
[[ -n "$BEFORE_ALL" ]] || fail "no before-all in [tool.cibuildwheel.linux] of ${PYPROJECT}"
BUILD_SELECTORS="$(toml_value tool.cibuildwheel build | tr -d '"')"
[[ -n "$BUILD_SELECTORS" ]] || fail "no build list in [tool.cibuildwheel] of ${PYPROJECT}"

IMAGE_HASH="$( { cat "$DEF_FILE"; printf '%s\n' "$BEFORE_ALL"; cat "$DEPS_SCRIPT"; } | sha256sum | cut -c1-12 )"
if [[ -z "$SIF" ]]; then
    SIF="${XDG_CACHE_HOME:-${HOME}/.cache}/lsdtt3/manylinux_2_28-lsdtt3-${IMAGE_HASH}.sif"
fi

# Resolve the Python versions to CPython tags (3.12 -> cp312).
TAGS=()
if [[ "$PYTHONS" == "all" ]]; then
    for sel in $BUILD_SELECTORS; do
        [[ "$sel" =~ ^(cp3[0-9]+)-\*$ ]] && TAGS+=("${BASH_REMATCH[1]}")
    done
else
    IFS=',' read -r -a wanted <<< "$PYTHONS"
    for v in "${wanted[@]}"; do
        [[ "$v" =~ ^3\.([0-9]+)$ ]] || fail "--python expects versions like 3.12, got '$v'"
        tag="cp3${BASH_REMATCH[1]}"
        [[ " $BUILD_SELECTORS " == *" ${tag}-* "* ]] || \
            fail "Python $v is not in the [tool.cibuildwheel] build list (${BUILD_SELECTORS})"
        TAGS+=("$tag")
    done
fi
((${#TAGS[@]})) || fail "no Python versions selected"

step "Plan"
info "$(printf '%-16s %s' "Apptainer:" "$(apptainer --version)")"
info "$(printf '%-16s %s' "Source:" "$SOURCE")"
info "$(printf '%-16s %s' "Output:" "$OUT_DIR")"
info "$(printf '%-16s %s' "Pythons:" "${TAGS[*]}")"
info "$(printf '%-16s %s' "Jobs:" "$JOBS")"
info "$(printf '%-16s %s' "Tests:" "$( ((RUN_TESTS)) && echo yes || echo "SKIPPED (--no-test)")")"
if [[ -f "$SIF" ]] && (( ! REBUILD_IMAGE )); then
    info "$(printf '%-16s %s' "Image:" "${SIF} (exists, reused)")"
else
    info "$(printf '%-16s %s' "Image:" "${SIF} (will be built: a few GB, needs network)")"
fi

if (( CHECK_ONLY )); then
    echo
    ok "check-only: apptainer and pyproject.toml are present; the plan above is what would run."
    exit 0
fi

# ----------------------------------------------------------------------------
# 1. The image
# ----------------------------------------------------------------------------
if [[ ! -f "$SIF" ]] || (( REBUILD_IMAGE )); then
    step "Building the wheel image"
    mkdir -p "$(dirname "$SIF")"
    STAGE="$(mktemp -d "${TMPDIR:-/tmp}/lsdtt3-wheel-image-XXXXXX")"
    trap 'rm -rf "$STAGE"' EXIT
    # %files in the .def copies these relative to the current directory, so
    # build from a directory holding just them.
    cp "$PYPROJECT" "${STAGE}/pyproject.toml"
    cp "$DEPS_SCRIPT" "${STAGE}/build_manylinux_native_deps.sh"
    tmp_sif="${SIF}.partial"
    rm -f "$tmp_sif"
    build_flags=(--fakeroot)
    # Without /etc/subuid entries Apptainer falls back to a root-mapped
    # namespace and runs %post under the host's fakeroot helper. That helper
    # is built against the host glibc, which is newer than the image's (2.28),
    # so it fails to start ("GLIBC_2.33 not found"). Skip it; the .def then
    # runs before-all under the image's own fakeroot instead.
    # Apptainer needs entries in BOTH files to avoid the fallback.
    # (Entries served by LDAP rather than these files would make this flag
    # unnecessary, but it is harmless then: no fakeroot command is injected
    # when the full mapping is available.)
    user_re="^($(id -un)|$(id -u)):"
    if ! grep -qE "$user_re" /etc/subuid 2>/dev/null || ! grep -qE "$user_re" /etc/subgid 2>/dev/null; then
        build_flags+=(--ignore-fakeroot-command)
        info "not in /etc/subuid and /etc/subgid: building with --ignore-fakeroot-command"
    fi
    # The build unpacks the image and installs the -devel packages in
    # APPTAINER_TMPDIR (default /tmp): about 3 GB. A small /tmp fails part
    # way through with "no space left on device", so warn up front.
    build_tmp="${APPTAINER_TMPDIR:-/tmp}"
    free_kb="$(df -Pk "$build_tmp" 2>/dev/null | awk 'NR==2 {print $4}')"
    if [[ -n "$free_kb" ]] && (( free_kb < 4 * 1024 * 1024 )); then
        warn "only $((free_kb / 1024)) MB free in ${build_tmp}; the image build needs about 3 GB there."
        warn "set APPTAINER_TMPDIR to a local directory with more space (not NFS) and rerun."
    fi
    ( cd "$STAGE" && apptainer build "${build_flags[@]}" "$tmp_sif" "$DEF_FILE" ) || fail "apptainer build failed.
  Without sudo this needs Apptainer >= 1.1 (unprivileged --fakeroot).
  'no space left on device': set APPTAINER_TMPDIR (and APPTAINER_CACHEDIR)
  to a local directory with a few GB free. NFS often fails for fakeroot
  builds, so prefer local disk. If it still fails, ask an administrator for
  /etc/subuid and /etc/subgid entries, or build the image on another machine with:
      cd <dir holding copies of python_bindings/pyproject.toml and
          python_bindings/tools/build_manylinux_native_deps.sh>
      apptainer build --fakeroot manylinux_wheel.sif ${DEF_FILE}
  and pass it here with --sif <path>."
    mv "$tmp_sif" "$SIF"
    rm -rf "$STAGE"
    trap - EXIT
    ok "image built: ${SIF}"
fi

if (( IMAGE_ONLY )); then
    # From /, which always exists in the image (avoids the chdir warning).
    ( cd / && apptainer test "$SIF" ) || fail "the image's %test section failed"
    ok "image ready: ${SIF}"
    exit 0
fi

# ----------------------------------------------------------------------------
# 2. Build, repair and test, one Python version at a time
# ----------------------------------------------------------------------------
mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"

# Runs inside the image as the calling user. $1 = cpXY tag, $2 = 1 to test,
# $3 = the file name of an abi3 wheel already in /wheelhouse to test with this
# Python instead of building one (empty: build). /project is the source tree
# and /wheelhouse the output directory.
read -r -d '' INNER <<'INNER_EOF' || true
set -euo pipefail
tag="$1"; run_tests="$2"; reuse="${3:-}"
py="/opt/python/${tag}-${tag}/bin/python"
[[ -x "$py" ]] || { echo "this image has no ${tag} (/opt/python: $(ls /opt/python | tr '\n' ' '))" >&2; exit 3; }

# Read the cibuildwheel settings with a Python that has tomllib.
cfg_py=/opt/python/cp312-cp312/bin/python
eval "$("$cfg_py" - <<'EOF'
import shlex, tomllib
with open("/project/python_bindings/pyproject.toml", "rb") as f:
    cfg = tomllib.load(f)["tool"]["cibuildwheel"]
for k, v in cfg.get("environment", {}).items():
    print(f"export {k}={shlex.quote(v)}")
print(f"TEST_REQUIRES={shlex.quote(' '.join(cfg.get('test-requires', [])))}")
print(f"TEST_COMMAND={shlex.quote(cfg.get('test-command', ''))}")
EOF
)"

work="$(mktemp -d /tmp/lsdtt3-wheel-${tag}-XXXXXX)"
trap 'rm -rf "$work"' EXIT

if [[ -n "$reuse" ]]; then
    wheel="/wheelhouse/${reuse}"
    echo "-- reuse ${reuse} (abi3: built once, tested with ${tag})"
else
    echo "-- build (${tag}, CMAKE_ARGS=${CMAKE_ARGS:-})"
    "$py" -m pip wheel --no-deps --wheel-dir "${work}/raw" /project/python_bindings

    echo "-- repair"
    auditwheel repair --wheel-dir "${work}/repaired" "${work}"/raw/*.whl
    wheel="$(ls "${work}"/repaired/*.whl)"
fi

if [[ "$run_tests" == 1 ]]; then
    echo "-- test $(basename "$wheel")"
    "$py" -m venv "${work}/venv"
    "${work}/venv/bin/python" -m pip install --quiet --upgrade pip
    # shellcheck disable=SC2086  # TEST_REQUIRES is a space-separated list
    "${work}/venv/bin/python" -m pip install --quiet $TEST_REQUIRES "$wheel"
    cmd="${TEST_COMMAND//\{project\}//project}"
    ( export PATH="${work}/venv/bin:${PATH}"; cd /tmp && bash -c "$cmd" )
fi

[[ -n "$reuse" ]] || cp "$wheel" /wheelhouse/
echo "-- done: $(basename "$wheel")"
INNER_EOF

# --cleanenv drops the caller's whole environment, including any proxy and
# pip index settings that pip needs to reach PyPI (build backend, pytest).
# Pass those few through explicitly.
pass_env=()
for v in http_proxy https_proxy no_proxy HTTP_PROXY HTTPS_PROXY NO_PROXY \
         PIP_INDEX_URL PIP_EXTRA_INDEX_URL PIP_TRUSTED_HOST PIP_CERT \
         SSL_CERT_FILE REQUESTS_CA_BUNDLE; do
    [[ -n "${!v:-}" ]] && pass_env+=(--env "${v}=${!v}")
done

BUILT=()
# The wheel is abi3 (cp312-abi3, pyproject.toml wheel.py-api), so the first
# Python builds it and every later one installs and tests that same wheel, as
# cibuildwheel does. ABI3_WHEEL holds its file name once built.
ABI3_WHEEL=""
for tag in "${TAGS[@]}"; do
    step "Wheel for ${tag}"
    t0=$SECONDS
    # --cleanenv: nothing from the caller's environment (a pixi or conda
    # shell's GDAL_DATA, PROJ_LIB, PYTHONPATH ...) leaks into the build.
    # PYTHONNOUSERSITE: the host home directory is mounted, so without it the
    # image's Pythons would import packages from ~/.local/lib/python3.X (an
    # old pip, say) that were installed on the host.
    # --pwd: the inner script uses only absolute paths. Without this Apptainer
    # tries to start in the caller's directory, which is not mounted in the
    # image when it is outside $HOME's mount (a harmless chdir warning).
    apptainer exec --cleanenv --pwd /project \
        --env "CMAKE_BUILD_PARALLEL_LEVEL=${JOBS}" \
        --env "PYTHONNOUSERSITE=1" \
        ${pass_env[@]+"${pass_env[@]}"} \
        --bind "${SOURCE}:/project" \
        --bind "${OUT_DIR}:/wheelhouse" \
        "$SIF" bash -c "$INNER" inner "$tag" "$RUN_TESTS" "$ABI3_WHEEL" \
        || fail "wheel build for ${tag} failed"
    BUILT+=("${tag} ($((SECONDS - t0))s)")
    if [[ -z "$ABI3_WHEEL" ]]; then
        newest="$(ls -t "$OUT_DIR"/*.whl 2>/dev/null | head -n 1)"
        [[ "$newest" == *-abi3-* ]] && ABI3_WHEEL="$(basename "$newest")"
    fi
done

step "Summary"
info "built: ${BUILT[*]}"
info "wheels in ${OUT_DIR}:"
ls -1 "$OUT_DIR"/*.whl | sed 's/^/    /'
(( RUN_TESTS )) || warn "tests were skipped (--no-test): do not publish these wheels."
ok "done"
