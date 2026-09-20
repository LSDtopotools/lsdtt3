#!/usr/bin/env bash
#
# make_all_packages.sh - build every lsdtt3 distributable artefact locally.
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
# Released under the GNU Affero General Public Licence v3 or later
# See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
#
# -----------------------------------------------------------------------------
# WHY THIS EXISTS
#
# The repository is private and lsdtt3 takes a long time to compile, so running
# the release and wheel workflows on GitHub Actions burns the CI minute
# allowance quickly. This script produces the same three artefacts on a local
# machine, costing nothing but wall-clock time:
#
#   a) the backend tarball      - what .github/workflows/release.yml uploads
#   b) the Debian package       - built from debian/ via dpkg-buildpackage
#   c) the Python wheel         - what .github/workflows/python-wheels.yml ships
#
# EACH TARGET COMPILES THE C++ SEPARATELY. There is no shared object cache
# between them: the tarball builds with the release CMake preset, the .deb
# builds under debian/rules with its own flags and its own build directory, and
# the wheel builds the bundled runtime inside python_bindings. Building all
# three is therefore roughly three full compiles. Use --only to build just the
# one you need while iterating.
#
# -----------------------------------------------------------------------------
# THE WHEEL IS THE ONE CAVEAT
#
# A DISTRIBUTABLE Linux wheel must be a manylinux wheel, which means compiling
# inside the manylinux container image so the wheel does not depend on this
# machine's glibc and GDAL. cibuildwheel drives that, and it needs a container
# engine.
#
# This project does not use Docker (see docs/developer_notes/). If podman is
# installed, this script uses it and produces a real manylinux wheel identical
# in kind to the CI artefact. If podman is absent it falls back to building a
# HOST wheel: correct, importable, fine for local testing, and NOT
# distributable, because it links whatever GDAL this machine has. The script
# says which one it built, and names the wheel accordingly.
#
# Docker is never invoked, even if present.
# -----------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ----------------------------------------------------------------------------
# Output helpers
# ----------------------------------------------------------------------------
if [[ -t 1 ]]; then
    C_RESET=$'\033[0m'; C_BOLD=$'\033[1m'; C_RED=$'\033[31m'
    C_GREEN=$'\033[32m'; C_YELLOW=$'\033[33m'; C_BLUE=$'\033[34m'
else
    C_RESET=""; C_BOLD=""; C_RED=""; C_GREEN=""; C_YELLOW=""; C_BLUE=""
fi

info()  { printf '%s\n' "${C_BLUE}==>${C_RESET} $*"; }
step()  { printf '\n%s\n' "${C_BOLD}${C_BLUE}=== $* ===${C_RESET}"; }
warn()  { printf '%s\n' "${C_YELLOW}warning:${C_RESET} $*" >&2; }
ok()    { printf '%s\n' "${C_GREEN}ok:${C_RESET} $*"; }
fail()  { printf '%s\n' "${C_RED}error:${C_RESET} $*" >&2; exit 1; }

# copy_source <src> <dst>
#
# Copy a source tree WITHOUT the heavy directories. This is not an
# optimisation: on a WSL2 checkout under /mnt/c the filesystem is NTFS, where
# copying many small files is slow enough that a blind `cp -a` of the whole
# repository (.git, build/, .pixi/, site/, dist/) takes longer than the
# compile it is preparing for. None of those directories belongs in a package
# build anyway.
copy_source() {
    local src="$1" dst="$2"
    mkdir -p "$dst"
    tar -C "$src" \
        --exclude=./.git \
        --exclude=./.pixi \
        --exclude='./build' \
        --exclude='./build_*' \
        --exclude='./cmake-build-*' \
        --exclude=./dist \
        --exclude=./site \
        --exclude=./wheelhouse \
        --exclude=./benchmark_campaigns \
        --exclude=./.venv \
        --exclude='./*.deb' \
        --exclude='*.o' \
        --exclude='__pycache__' \
        -cf - . | tar -x -C "$dst"
}

# collect <src-file> <dest-path>
#
# Move a finished artefact into the output directory. Deliberately copy +
# verify + remove rather than `mv`: the scratch tree is on the Linux
# filesystem while the output directory is typically on /mnt/c (drvfs), and a
# cross-filesystem `mv` there has been observed to complete the copy and still
# report a non-zero status, which under `set -e` aborts a run whose artefact is
# actually fine. Verifying the byte count is the check that matters anyway.
collect() {
    local src="$1" dest="$2"
    [[ -f "$src" ]] || { warn "collect: no such file: $src"; return 1; }
    # Already where it belongs. make_release.sh writes into dist/, which is
    # also this script's default --out-dir, so the tarball routinely arrives
    # pre-collected. Copying a file onto itself is an error, not a no-op.
    if [[ "$src" -ef "$dest" ]]; then
        return 0
    fi
    local src_size dest_size
    src_size="$(stat -c%s "$src")"
    cp -f "$src" "$dest" || { warn "collect: copy failed: $src -> $dest"; return 1; }
    dest_size="$(stat -c%s "$dest" 2>/dev/null || echo 0)"
    if [[ "$src_size" != "$dest_size" ]]; then
        warn "collect: size mismatch after copy ($src_size vs $dest_size): $dest"
        return 1
    fi
    rm -f "$src"
    return 0
}

usage() {
    cat <<'EOF'
Usage:
  tools/make_all_packages.sh [--clean | --dirty] [options]

Builds the backend tarball, the Debian package and the Python wheel locally,
so a private repository does not spend GitHub Actions minutes on them.

Source selection (choose one; --clean is the default):
  --clean          Build the latest tagged version. The tree is exported with
                   `git archive <tag>` into a scratch directory, so the build
                   sees exactly what the tag contains. Your working tree is
                   not touched and uncommitted changes are NOT included.
  --dirty          Build the current working tree, uncommitted changes and
                   all. Versions get a suffix derived from `git describe`, so
                   a dirty artefact can never be mistaken for a release.

Target selection (default: all three):
  --only <target>  Build only this target. Repeatable.
  --skip <target>  Skip this target. Repeatable.
                   Targets: tarball, deb, wheel

Options:
  --out-dir <path> Where finished artefacts are collected.
                   Default: dist/
  --jobs <N>       Parallel compile jobs. Default: nproc.
  --profile <name> CLI profile for the tarball: qgis / core / all.
                   Default: core (what the official Linux release ships).
  --keep-scratch   Keep the scratch build trees for inspection.
  --check-only     Check tooling and print the plan, then exit without
                   building. Use this first.
  --help           Show this help.

Examples:
  tools/make_all_packages.sh --check-only
  tools/make_all_packages.sh --clean
  tools/make_all_packages.sh --dirty --only wheel
  tools/make_all_packages.sh --dirty --skip deb --jobs 4

Notes:
  * Each target compiles the C++ independently; all three is ~3 full builds.
  * Without podman the wheel is a host wheel, not a distributable manylinux
    wheel. See the header comment. Docker is never used.
EOF
}

# ----------------------------------------------------------------------------
# Argument parsing
# ----------------------------------------------------------------------------
MODE="clean"
OUT_DIR=""
JOBS=""
CLI_PROFILE="core"
KEEP_SCRATCH=0
CHECK_ONLY=0
ONLY=()
SKIP=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --clean)        MODE="clean"; shift ;;
        --dirty)        MODE="dirty"; shift ;;
        --only)         ONLY+=("${2:-}"); shift 2 ;;
        --skip)         SKIP+=("${2:-}"); shift 2 ;;
        --out-dir)      OUT_DIR="${2:-}"; shift 2 ;;
        --jobs)         JOBS="${2:-}"; shift 2 ;;
        --profile)      CLI_PROFILE="${2:-}"; shift 2 ;;
        --keep-scratch) KEEP_SCRATCH=1; shift ;;
        --check-only)   CHECK_ONLY=1; shift ;;
        --help|-h)      usage; exit 0 ;;
        *)              fail "unknown option: $1 (try --help)" ;;
    esac
done

JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
OUT_DIR="${OUT_DIR:-${REPO_ROOT}/dist}"
[[ "$OUT_DIR" = /* ]] || OUT_DIR="${REPO_ROOT}/${OUT_DIR}"

# Resolve which targets to build.
want_tarball=1; want_deb=1; want_wheel=1
if ((${#ONLY[@]})); then
    want_tarball=0; want_deb=0; want_wheel=0
    for t in "${ONLY[@]}"; do
        case "$t" in
            tarball) want_tarball=1 ;;
            deb)     want_deb=1 ;;
            wheel)   want_wheel=1 ;;
            *)       fail "unknown target for --only: $t" ;;
        esac
    done
fi
for t in "${SKIP[@]-}"; do
    case "$t" in
        "")      ;;
        tarball) want_tarball=0 ;;
        deb)     want_deb=0 ;;
        wheel)   want_wheel=0 ;;
        *)       fail "unknown target for --skip: $t" ;;
    esac
done
(( want_tarball || want_deb || want_wheel )) || fail "every target was skipped"

cd "$REPO_ROOT"

# ----------------------------------------------------------------------------
# Version resolution
#
# Three artefacts, three version grammars, one git state:
#
#   tarball  make_release.sh compares the part before the first '-' against
#            project(lsdtt3 VERSION ...) and treats the rest as a pre-release
#            suffix. So a hyphen is exactly what we want here.
#   deb      debian/source/format is '3.0 (native)', so the version must NOT
#            contain a hyphen - a hyphen would be read as separating a Debian
#            revision. Use '+git...' instead.
#   wheel    PEP 440. python_bindings/pyproject.toml carries its own version,
#            independent of the lsdtt3 version, so the dev suffix is appended
#            to THAT as a local version identifier.
# ----------------------------------------------------------------------------
git rev-parse --git-dir >/dev/null 2>&1 || fail "not a git repository: $REPO_ROOT"

LATEST_TAG="$(git describe --tags --abbrev=0 2>/dev/null || true)"
[[ -n "$LATEST_TAG" ]] || fail "no git tags found; cannot determine a version"
BASE_VERSION="${LATEST_TAG#v}"

CMAKE_VERSION="$(sed -n 's/^project(lsdtt3 VERSION \([0-9.]*\).*/\1/p' CMakeLists.txt)"
[[ -n "$CMAKE_VERSION" ]] || fail "could not read project version from CMakeLists.txt"

PY_BASE_VERSION="$(sed -n 's/^version *= *"\([^"]*\)".*/\1/p' python_bindings/pyproject.toml | head -1)"
[[ -n "$PY_BASE_VERSION" ]] || fail "could not read version from python_bindings/pyproject.toml"

if [[ "$MODE" == "clean" ]]; then
    TARBALL_VERSION="$BASE_VERSION"
    DEB_VERSION="$BASE_VERSION"
    WHEEL_VERSION="$PY_BASE_VERSION"
    SOURCE_DESC="tag ${LATEST_TAG}"

    if [[ "$BASE_VERSION" != "$CMAKE_VERSION" ]]; then
        warn "tag ${LATEST_TAG} implies ${BASE_VERSION} but CMakeLists.txt says ${CMAKE_VERSION}."
        warn "In --clean mode the TAGGED CMakeLists.txt is what gets built, so this"
        warn "mismatch reflects the working tree having moved on. That is expected"
        warn "between releases and is not an error."
    fi
else
    DESCRIBE="$(git describe --tags --dirty --always 2>/dev/null || echo "$LATEST_TAG")"
    # e.g. v0.5.1-14-g04460f4-dirty  ->  14-g04460f4-dirty
    SUFFIX="${DESCRIBE#"$LATEST_TAG"}"
    SUFFIX="${SUFFIX#-}"
    if [[ -z "$SUFFIX" ]]; then
        # Exactly on the tag, but the user asked for dirty: still mark it, so
        # the artefact is never confused with the released one.
        SUFFIX="0-g$(git rev-parse --short HEAD)"
    fi
    if [[ -n "$(git status --porcelain 2>/dev/null)" && "$SUFFIX" != *dirty* ]]; then
        SUFFIX="${SUFFIX}-dirty"
    fi

    TARBALL_VERSION="${CMAKE_VERSION}-${SUFFIX}"
    DEB_VERSION="${CMAKE_VERSION}+git$(printf '%s' "$SUFFIX" | tr -c 'A-Za-z0-9.+~' '.')"
    WHEEL_VERSION="${PY_BASE_VERSION}+git$(printf '%s' "$SUFFIX" | tr -c 'A-Za-z0-9' '.')"
    SOURCE_DESC="working tree (${DESCRIBE})"
fi

# ----------------------------------------------------------------------------
# Tooling checks - fail fast, and say exactly what is missing
# ----------------------------------------------------------------------------
have() { command -v "$1" >/dev/null 2>&1; }

MISSING=()
need() { have "$1" || MISSING+=("$1 -- $2"); }

need git "version resolution"
need cmake "compiling lsdtt3"

if (( want_tarball )); then
    [[ -x tools/make_release.sh ]] || MISSING+=("tools/make_release.sh -- not executable")
fi

if (( want_deb )); then
    need dpkg-buildpackage "building the .deb"
    need dh "debhelper, building the .deb"
    [[ -f debian/control ]] || MISSING+=("debian/control -- missing")
fi

# Interpreter selection for the wheel.
#
# `python -m build` creates an isolated venv and installs the backend pinned in
# pyproject.toml INTO THAT VENV, so nothing is installed into the interpreter we
# pick. That needs two things from the interpreter: ensurepip (to make the venv)
# and the `build` module itself.
#
# The system python on Debian/Ubuntu often lacks ensurepip, because it is split
# into a separate python3-venv package. The project's pixi environments do not
# have that problem - they are conda-forge Pythons with ensurepip and build
# already present - so prefer one of those and fall back to the system python.
# This also keeps the wheel build inside the environment tooling the project
# already uses, rather than depending on what happens to be installed globally.
WHEEL_ENGINE="none"
WHEEL_PYTHON=""
if (( want_wheel )); then
    # Probe from a directory that is NOT the repository root. Python puts the
    # current directory on sys.path, and this repo has a build/ DIRECTORY, so
    # `import build` run from the root succeeds as a namespace package and
    # reports a build tool that is not there. Probing from / avoids that.
    usable_python() {
        local py="$1"
        [[ -x "$py" ]] || return 1
        ( cd / && "$py" -c '
import importlib.util as u, sys
need = ["ensurepip", "build"]
missing = [m for m in need if u.find_spec(m) is None]
sys.exit(1 if missing else 0)
' ) >/dev/null 2>&1
    }

    for candidate in \
        "${REPO_ROOT}/.pixi/envs/packaging/bin/python" \
        "${REPO_ROOT}/.pixi/envs/default/bin/python" \
        "${REPO_ROOT}/.pixi/envs/docs/bin/python" \
        "$(command -v python3 || true)"
    do
        if usable_python "$candidate"; then
            WHEEL_PYTHON="$candidate"
            break
        fi
    done

    if [[ -z "$WHEEL_PYTHON" ]]; then
        MISSING+=("a Python with BOTH ensurepip and the 'build' module.
      Checked: .pixi/envs/default, .pixi/envs/docs, and the system python3.
      As things stand neither side has both:
        * the pixi environments have ensurepip but not 'build'
        * the system python has 'build' but not ensurepip (python3-venv)

      Two ways to fix it, in order of preference:

      1. Add 'build' to a pixi environment (keeps Python tooling in pixi,
         which is this project's convention). In pixi.toml:

             [feature.packaging.pypi-dependencies]
             build = \"*\"

             [environments]
             packaging = [\"packaging\"]

         then 'pixi install', and re-run this script.

      2. Fix the system python instead:

             sudo apt install python3-venv

      This script never installs build backends into your interpreter, so an
      isolated build is the only supported route.")
    fi

    if have podman; then
        WHEEL_ENGINE="podman"
    fi
fi

if ((${#MISSING[@]})); then
    printf '%s\n' "${C_RED}Missing prerequisites:${C_RESET}" >&2
    printf '  - %s\n' "${MISSING[@]}" >&2
    cat >&2 <<'EOF'

On Ubuntu/Debian most of these come from:
  sudo apt install cmake g++ pkg-config libgdal-dev libfftw3-dev \
                   libnetcdf-dev libeigen3-dev dpkg-dev debhelper \
                   python3-pip python3-venv
EOF
    exit 1
fi

# ----------------------------------------------------------------------------
# Plan
# ----------------------------------------------------------------------------
step "Plan"
printf '  %-22s %s\n' "Mode:"        "$MODE"
printf '  %-22s %s\n' "Source:"      "$SOURCE_DESC"
printf '  %-22s %s\n' "Latest tag:"  "$LATEST_TAG"
printf '  %-22s %s\n' "CMake version:" "$CMAKE_VERSION"
printf '  %-22s %s\n' "Output:"      "$OUT_DIR"
printf '  %-22s %s\n' "Jobs:"        "$JOBS"
echo
printf '  %-22s %s\n' "Tarball:" "$( ((want_tarball)) && echo "yes  (version ${TARBALL_VERSION}, profile ${CLI_PROFILE})" || echo "skipped" )"
printf '  %-22s %s\n' "Debian package:" "$( ((want_deb)) && echo "yes  (version ${DEB_VERSION})" || echo "skipped" )"
if (( want_wheel )); then
    if [[ "$WHEEL_ENGINE" == "podman" ]]; then
        printf '  %-22s %s\n' "Python wheel:" "yes  (version ${WHEEL_VERSION}, manylinux via podman)"
    else
        printf '  %-22s %s\n' "Python wheel:" "yes  (version ${WHEEL_VERSION}, HOST wheel - not distributable)"
    fi
else
    printf '  %-22s %s\n' "Python wheel:" "skipped"
fi

if (( want_wheel )) && [[ "$WHEEL_ENGINE" != "podman" ]]; then
    echo
    warn "PYTHON WHEEL ONLY -- this does not affect the tarball or the .deb,"
    warn "which build fully and correctly on this machine with no container."
    warn ""
    warn "podman is not installed, so the WHEEL will be built against this"
    warn "machine's libraries. It will import and run here, but it is not a"
    warn "manylinux wheel and must not be published. For a distributable wheel:"
    warn "    sudo apt install podman"
    warn "Docker is deliberately not used by this project, and is never invoked."
fi

n_targets=$(( want_tarball + want_deb + want_wheel ))
if (( n_targets > 1 )); then
    echo
    info "${n_targets} targets selected; each compiles the C++ separately."
fi

if (( CHECK_ONLY )); then
    echo
    ok "check-only: tooling is present and the plan above is what would run."
    exit 0
fi

# ----------------------------------------------------------------------------
# Scratch space
# ----------------------------------------------------------------------------
SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/lsdtt3-pkg-XXXXXX")"
cleanup() {
    if (( KEEP_SCRATCH )); then
        info "scratch kept at ${SCRATCH}"
    else
        rm -rf "$SCRATCH"
    fi
}
trap cleanup EXIT

mkdir -p "$OUT_DIR"

# In clean mode, export the tagged tree and build from there. git archive is
# used rather than a worktree: it cannot pick up anything from the working
# tree, needs no bookkeeping, and leaves no state behind if this script dies.
if [[ "$MODE" == "clean" ]]; then
    SRC_DIR="${SCRATCH}/src"
    mkdir -p "$SRC_DIR"
    info "exporting ${LATEST_TAG} to a scratch tree"
    git archive --format=tar "$LATEST_TAG" | tar -x -C "$SRC_DIR"
    [[ -f "${SRC_DIR}/CMakeLists.txt" ]] || fail "export produced no CMakeLists.txt"
else
    SRC_DIR="$REPO_ROOT"
fi

# ----------------------------------------------------------------------------
# Per-target source checks
#
# A tag can predate a target. v0.5.1, for example, was cut three days before
# debian/ was added, so `--clean --only deb` against that tag has nothing to
# build from. Catch that here with a clear message rather than letting
# dpkg-buildpackage fail obscurely after the checks have already passed.
# ----------------------------------------------------------------------------
missing_in_source=()
if (( want_tarball )) && [[ ! -x "${SRC_DIR}/tools/make_release.sh" ]]; then
    missing_in_source+=("tarball needs tools/make_release.sh")
fi
if (( want_deb )) && [[ ! -f "${SRC_DIR}/debian/control" ]]; then
    missing_in_source+=("deb needs debian/control")
fi
if (( want_wheel )) && [[ ! -f "${SRC_DIR}/python_bindings/pyproject.toml" ]]; then
    missing_in_source+=("wheel needs python_bindings/pyproject.toml")
fi

if ((${#missing_in_source[@]})); then
    printf '%s\n' "${C_RED}Target(s) not buildable from this source:${C_RESET}" >&2
    printf '  - %s\n' "${missing_in_source[@]}" >&2
    if [[ "$MODE" == "clean" ]]; then
        cat >&2 <<EOF

You are in --clean mode, building ${LATEST_TAG}, and that tag does not contain
everything the selected targets need. This happens when a target was added to
the repository after the last release was tagged.

Either:
  * cut a newer tag that includes it, then re-run; or
  * use --dirty to build from the working tree instead; or
  * drop the target, e.g. --skip deb
EOF
    else
        echo >&2
        echo "The working tree is missing these files." >&2
    fi
    exit 1
fi

declare -a RESULTS=()
declare -a TIMINGS=()
OVERALL_START=$SECONDS

record() { RESULTS+=("$1"); }
timing() { TIMINGS+=("$1"); }

# ----------------------------------------------------------------------------
# a) Backend tarball
# ----------------------------------------------------------------------------
if (( want_tarball )); then
    step "a) Backend tarball"
    t0=$SECONDS
    (
        cd "$SRC_DIR"
        tools/make_release.sh \
            --version "$TARBALL_VERSION" \
            --profile "$CLI_PROFILE" \
            --build-dir "${SCRATCH}/build-tarball"
    ) || fail "tarball build failed"

    shopt -s nullglob
    produced=( "${SRC_DIR}"/lsdtt3-backend-*.tar.gz "${SRC_DIR}"/dist/lsdtt3-backend-*.tar.gz )
    shopt -u nullglob
    if ((${#produced[@]})); then
        for f in "${produced[@]}"; do
            collect "$f" "${OUT_DIR}/$(basename "$f")" \
                || fail "could not collect $(basename "$f")"
            record "tarball: $(basename "$f")"
        done
    else
        warn "make_release.sh reported success but no tarball was found."
        warn "Look in ${SRC_DIR} for the bundle and check tools/README.md."
        record "tarball: BUILT, but not collected (see warning)"
    fi
    timing "tarball: $((SECONDS - t0))s"
    ok "tarball stage finished"
fi

# ----------------------------------------------------------------------------
# b) Debian package
#
# dpkg-buildpackage writes its output to the PARENT of the source directory,
# so the source is copied into a scratch subdirectory first. That also keeps
# debian/rules from leaving build artefacts in the repository.
# ----------------------------------------------------------------------------
if (( want_deb )); then
    step "b) Debian package"
    t0=$SECONDS

    DEB_SRC="${SCRATCH}/deb/lsdtt3-${DEB_VERSION}"
    mkdir -p "$(dirname "$DEB_SRC")"
    info "staging source for dpkg-buildpackage"
    copy_source "$SRC_DIR" "$DEB_SRC"

    # Stamp the version into debian/changelog. The existing top entry is
    # rewritten rather than appended to, because this is a throwaway build
    # tree and a real changelog entry belongs in a commit, not in a script.
    CHANGELOG="${DEB_SRC}/debian/changelog"
    [[ -f "$CHANGELOG" ]] || fail "no debian/changelog in the build tree"
    DIST_CODENAME="$(sed -n '1s/.*) \([^;]*\);.*/\1/p' "$CHANGELOG")"
    DIST_CODENAME="${DIST_CODENAME:-jammy}"
    {
        printf 'lsdtt3 (%s) %s; urgency=medium\n\n' "$DEB_VERSION" "$DIST_CODENAME"
        if [[ "$MODE" == "clean" ]]; then
            printf '  * Local build of %s.\n\n' "$LATEST_TAG"
        else
            printf '  * Local development build from %s.\n' "$SOURCE_DESC"
            printf '  * NOT A RELEASE - built by tools/make_all_packages.sh --dirty.\n\n'
        fi
        printf ' -- %s <%s>  %s\n' \
            "$(git config user.name  || echo 'lsdtt3 build')" \
            "$(git config user.email || echo 'noreply@example.invalid')" \
            "$(date -R)"
    } > "${CHANGELOG}.new"
    mv "${CHANGELOG}.new" "$CHANGELOG"

    # debian/rules and the smoke test must be executable. On an NTFS checkout
    # the exec bit does not survive, so set it explicitly rather than trusting
    # the copy (see docs/developer_notes/git_notes.md).
    chmod +x "${DEB_SRC}/debian/rules" 2>/dev/null || true
    [[ -f "${DEB_SRC}/debian/smoke-test.sh" ]] && chmod +x "${DEB_SRC}/debian/smoke-test.sh"

    (
        cd "$DEB_SRC"
        # -b       binary only (a native source package adds nothing here)
        # -uc -us  do not sign; signing is a release step, not a build step
        DEB_BUILD_OPTIONS="parallel=${JOBS}" \
            dpkg-buildpackage -b -uc -us
    ) || fail "dpkg-buildpackage failed"

    shopt -s nullglob
    debs=( "${SCRATCH}/deb"/*.deb )
    shopt -u nullglob
    if ((${#debs[@]})); then
        for f in "${debs[@]}"; do
            collect "$f" "${OUT_DIR}/$(basename "$f")" \
                || fail "could not collect $(basename "$f")"
            record "deb: $(basename "$f")"
        done
    else
        fail "dpkg-buildpackage reported success but produced no .deb"
    fi
    timing "deb: $((SECONDS - t0))s"
    ok "Debian package stage finished"
fi

# ----------------------------------------------------------------------------
# c) Python wheel
# ----------------------------------------------------------------------------
if (( want_wheel )); then
    step "c) Python wheel"
    t0=$SECONDS

    WHEEL_SRC="${SCRATCH}/wheel"
    info "staging source for the wheel build"
    copy_source "$SRC_DIR" "$WHEEL_SRC"

    # Stamp the version into the copied pyproject.toml only.
    PYPROJECT="${WHEEL_SRC}/python_bindings/pyproject.toml"
    [[ -f "$PYPROJECT" ]] || fail "no python_bindings/pyproject.toml in the build tree"
    "${WHEEL_PYTHON:-python3}" - "$PYPROJECT" "$WHEEL_VERSION" <<'PY'
import re, sys
path, version = sys.argv[1], sys.argv[2]
text = open(path).read()
new, n = re.subn(r'(?m)^version *= *"[^"]*"', f'version = "{version}"', text, count=1)
if n != 1:
    sys.exit(f"could not rewrite version in {path}")
open(path, "w").write(new)
PY

    if [[ "$WHEEL_ENGINE" == "podman" ]]; then
        info "building a manylinux wheel with cibuildwheel (podman)"
        have cibuildwheel || "$WHEEL_PYTHON" -c 'import cibuildwheel' >/dev/null 2>&1 || \
            fail "cibuildwheel is not installed, and this script does not install
into your interpreter. Install it in an environment you manage:
    pipx install cibuildwheel
or
    pip install --user cibuildwheel"
        (
            cd "$WHEEL_SRC"
            # Only the running interpreter's version, to keep this to one
            # compile rather than one per supported Python.
            py_tag="cp$("$WHEEL_PYTHON" -c 'import sys; print(f"{sys.version_info.major}{sys.version_info.minor}")')-*"
            CIBW_CONTAINER_ENGINE=podman \
            CIBW_BUILD="$py_tag" \
            CIBW_MANYLINUX_X86_64_IMAGE=manylinux_2_28 \
            CIBW_ENVIRONMENT='CMAKE_ARGS="-DLSDTT_PYTHON_RELEASE_WHEEL=ON"' \
                "$WHEEL_PYTHON" -m cibuildwheel python_bindings --output-dir wheelhouse
        ) || fail "cibuildwheel failed"
        WHEEL_GLOB="${WHEEL_SRC}/wheelhouse/*.whl"
        WHEEL_KIND="manylinux"
    else
        info "building a host wheel (not distributable)"
        info "interpreter: ${WHEEL_PYTHON}"
        # Isolated build: `build` creates a venv and installs the backend
        # (scikit-build-core, nanobind at the pin in pyproject.toml) into it.
        # Nothing is installed into the caller's interpreter, and the pins are
        # honoured exactly rather than resolved against whatever happens to be
        # installed here.
        (
            cd "$WHEEL_SRC"
            CMAKE_ARGS="-DLSDTT_PYTHON_RELEASE_WHEEL=ON" \
            CMAKE_BUILD_PARALLEL_LEVEL="$JOBS" \
                "$WHEEL_PYTHON" -m build --wheel python_bindings --outdir wheelhouse
        ) || fail "wheel build failed"
        WHEEL_GLOB="${WHEEL_SRC}/wheelhouse/*.whl"
        WHEEL_KIND="host"
    fi

    shopt -s nullglob
    wheels=( $WHEEL_GLOB )
    shopt -u nullglob
    ((${#wheels[@]})) || fail "no wheel was produced"
    for f in "${wheels[@]}"; do
        base="$(basename "$f")"
        if [[ "$WHEEL_KIND" == "host" ]]; then
            # Make it obvious on disk that this one must not be published.
            base="LOCAL-ONLY-${base}"
        fi
        collect "$f" "${OUT_DIR}/${base}" || fail "could not collect ${base}"
        record "wheel: ${base}"
    done
    timing "wheel: $((SECONDS - t0))s"
    ok "wheel stage finished"
fi

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
step "Summary"
printf '  %-22s %s\n' "Mode:"   "$MODE"
printf '  %-22s %s\n' "Source:" "$SOURCE_DESC"
echo
echo "  Artefacts in ${OUT_DIR}:"
if ((${#RESULTS[@]})); then
    printf '    %s\n' "${RESULTS[@]}"
else
    echo "    (none)"
fi
echo
echo "  Time:"
printf '    %s\n' "${TIMINGS[@]}"
printf '    %-12s %ss\n' "total:" "$((SECONDS - OVERALL_START))"

if (( want_wheel )) && [[ "$WHEEL_ENGINE" != "podman" ]]; then
    echo
    warn "The wheel is prefixed LOCAL-ONLY- because it links this machine's"
    warn "libraries. Do not upload it. Install podman, or use the GitHub"
    warn "workflow, for a publishable manylinux wheel."
fi

if [[ "$MODE" == "dirty" ]]; then
    echo
    warn "These artefacts were built from an untagged tree and are versioned"
    warn "accordingly (${TARBALL_VERSION}). They are for testing, not release."
fi

echo
ok "done"
