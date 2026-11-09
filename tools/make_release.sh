#!/usr/bin/env bash
#
# make_release.sh - single entry point for building an lsdtt3 release.
#
# Used identically by GitHub Actions (on a v* tag) and by a local machine
# (test/pre-releases between versions), so local releases exercise exactly
# the path CI uses. See tools/README.md for the full release workflow.
#
# Steps:
#   1. Resolve the release version (--version flag, else exact git tag).
#   2. Check it against project(lsdtt3 VERSION ...) in CMakeLists.txt.
#      Pre-release versions (e.g. 0.3.0-rc1) compare only the X.Y.Z stem.
#      Then check that every copy of the version number agrees
#      (tests/scripts/check_version_consistency.py; fix with tools/set_version.py).
#   3. Configure and build (Release; Eigen on; tests off - a release build
#      packages nothing from the CTest suite, so compiling it is wasted CI
#      minutes). Then refuse to continue if an optional dependency was missing:
#      they degrade silently, so this is the only place it can be caught.
#   4. Smoke-check every packaged CLI (runs each with -v and checks the output
#      reports GDAL, not merely that the process exited 0). Then build the
#      Python bindings once and run the whole python_bindings/tests suite
#      against them (tests/scripts/check_python_bindings_smoke.sh). This is
#      mandatory: a release whose C++ breaks the bindings does not ship.
#   5. Package the Linux backend bundle (tools/package_linux_backend.sh).
#   6. Optionally publish to a GitHub release with gh (--publish).

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  tools/make_release.sh [options]

Options:
  --version <X.Y.Z[-suffix]>  Release version. Default: taken from an exact
                              v* git tag on HEAD (git describe --exact-match).
  --profile <name>            CLI set to ship: qgis / core / all
                              (default: core). Defined in
                              tools/cli_profiles.txt; all = every lsdtt-*
                              executable in the build directory. Official
                              Linux releases use core; all is for local
                              releases only.
  --build-dir <path>          Build directory (default: build/release)
  --skip-build                Reuse the build directory without reconfiguring
                              or rebuilding (fast repack iterations)
  --publish                   Create/update GitHub release v<version> and
                              upload the bundle (requires gh, authenticated)
  --prerelease                Mark the published release as a prerelease
  --cli <name>                Override the CLI list entirely (repeatable;
                              bundle is labelled 'custom')
  --smoke-python <exe>        Interpreter for the Python bindings smoke check
                              (default: python3; needs the venv module and
                              network access or a pip cache)
  --help                      Show this help

Examples:
  tools/make_release.sh --version 0.3.0                # local build, core set
  tools/make_release.sh --version 0.3.0 --profile all --publish --prerelease
  tools/make_release.sh                                # on a tagged commit
EOF
}

VERSION=""
CLI_PROFILE="core"
BUILD_DIR="build/release"
SKIP_BUILD=0
PUBLISH=0
PRERELEASE=0
CLI_OVERRIDE=()
SMOKE_PYTHON="python3"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --version)    VERSION="${2:-}"; shift 2 ;;
        --profile)    CLI_PROFILE="${2:-}"; shift 2 ;;
        --build-dir)  BUILD_DIR="${2:-}"; shift 2 ;;
        --skip-build) SKIP_BUILD=1; shift ;;
        --publish)    PUBLISH=1; shift ;;
        --prerelease) PRERELEASE=1; shift ;;
        --cli)        CLI_OVERRIDE+=("${2:-}"); shift 2 ;;
        --smoke-python) SMOKE_PYTHON="${2:-}"; shift 2 ;;
        --help|-h)    usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

case "$CLI_PROFILE" in qgis|core|all) ;; *)
    echo "ERROR: unknown --profile '$CLI_PROFILE' (use qgis, core, or all)" >&2
    exit 1 ;;
esac

info()  { echo "[make_release] $*"; }
fail()  { echo "[make_release] ERROR: $*" >&2; exit 1; }

# Every optional dependency in this tree degrades SILENTLY: a missing library is
# a CMake STATUS or WARNING line, never an error, so the build succeeds and the
# bundle ships without the feature. Nothing downstream notices - the affected
# CLIs still build, still run, and still exit 0 on -v. netCDF is the worst case,
# because the cosmo tools then fall back to the standard atmosphere and produce
# plausible, different numbers. Read the configure log and refuse to release a
# crippled bundle.
check_optional_deps() {
    local log="$1"
    [[ -f "$log" ]] || fail "configure log not found: $log"

    grep -Eq "FFTW3 +: +(ON|TRUE|1)" "$log" || fail \
        "FFTW3 was not found at configure time: lsdtt-spectral-analysis will not
be built and the bundle would ship without it.
Install libfftw3-dev (and pkg-config, which is how it is discovered on Linux)."

    grep -Eq "Eigen3 +: +(ON|TRUE|1)" "$log" || fail \
        "Eigen3 was not found at configure time: chi / k_sn inversion is compiled
out of lsdtt-chi-analysis, which still builds and still runs.
Install libeigen3-dev."

    if grep -q "netCDF-C not found" "$log"; then
        fail "netCDF was not found at configure time: ERA40 atmosphere and
palaeomagnetic corrections are compiled out of the cosmo tools, which still
build and still run - they silently fall back to the standard atmosphere.
Install libnetcdf-dev."
    fi

    info "optional dependencies present: FFTW3, Eigen3, netCDF"
}

# The set of CLIs that currently exist is declared in CMakeLists.txt, not by
# whatever happens to be sitting in the build directory. A long-lived build tree
# accumulates binaries for tools that have since been renamed or removed, and
# they are indistinguishable from current ones by name, mode or mtime. Packaging
# one ships a tool that no longer has source behind it.
current_cli_targets() {
    sed -n 's/^[[:space:]]*lsdtt_cli(\(lsdtt-[A-Za-z0-9_-]*\).*/\1/p' CMakeLists.txt
}

# -----------------------------------------------------------------------------
# 1. Resolve version
# -----------------------------------------------------------------------------
if [[ -z "$VERSION" ]]; then
    tag="$(git describe --tags --exact-match 2>/dev/null || true)"
    [[ -n "$tag" ]] || fail "no --version given and HEAD carries no exact v* tag.
Either tag the commit (git tag vX.Y.Z) or pass --version X.Y.Z."
    [[ "$tag" == v* ]] || fail "tag on HEAD is '$tag'; release tags must look like vX.Y.Z"
    VERSION="${tag#v}"
fi
info "release version: $VERSION"

# -----------------------------------------------------------------------------
# 2. Check against CMakeLists.txt
# -----------------------------------------------------------------------------
cmake_version="$(sed -n 's/^project(lsdtt3 VERSION \([0-9.]*\).*/\1/p' CMakeLists.txt)"
[[ -n "$cmake_version" ]] || fail "could not read project version from CMakeLists.txt"

version_stem="${VERSION%%-*}"
if [[ "$VERSION" == "$version_stem" ]]; then
    # Full release: exact match required.
    [[ "$VERSION" == "$cmake_version" ]] || fail \
        "version mismatch: releasing $VERSION but CMakeLists.txt says $cmake_version.
Bump project(lsdtt3 VERSION ...) and commit before tagging."
else
    # Pre-release (e.g. 0.3.0-rc1): compare only the X.Y.Z stem, warn on drift.
    if [[ "$version_stem" != "$cmake_version" ]]; then
        info "WARNING: pre-release stem $version_stem != CMake version $cmake_version (continuing)"
    fi
fi

python3 tests/scripts/check_version_consistency.py --repo-root "$REPO_ROOT" \
    || fail "the version copies disagree (above). Run tools/set_version.py $cmake_version
and commit before releasing."

if ! git diff --quiet || ! git diff --cached --quiet; then
    info "WARNING: working tree has uncommitted changes; binaries will report -dirty"
fi

# -----------------------------------------------------------------------------
# 3. Configure and build
# -----------------------------------------------------------------------------
if [[ "$SKIP_BUILD" -eq 1 ]]; then
    [[ -d "$BUILD_DIR" ]] || fail "--skip-build given but build dir not found: $BUILD_DIR"
    info "skipping build; using existing $BUILD_DIR"
    if [[ -f "$BUILD_DIR/configure.log" ]]; then
        check_optional_deps "$BUILD_DIR/configure.log"
    else
        info "WARNING: no configure.log in $BUILD_DIR; cannot confirm that FFTW3,
Eigen3 and netCDF were found for this build. Reconfigure without --skip-build to
check, or accept that the bundle may silently lack those features."
    fi
else
    info "configuring $BUILD_DIR"
    # BUILD_TESTING is ON by default (CMakeLists.txt calls include(CTest)), which
    # would compile the whole CTest suite into every release build - minutes that
    # CI is billed for and that package nothing. Test before releasing, not while.
    configure_log="$BUILD_DIR/configure.log"
    mkdir -p "$BUILD_DIR"
    cmake -B "$BUILD_DIR" \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_TESTING=OFF \
        -DLSDTT_EIGEN=ON 2>&1 | tee "$configure_log"
    # set -o pipefail is in force, so a cmake failure aborts here despite tee.
    check_optional_deps "$configure_log"
    n_jobs="$(nproc 2>/dev/null || echo 4)"
    info "building with $n_jobs jobs"
    cmake --build "$BUILD_DIR" -j "$n_jobs"
fi

# -----------------------------------------------------------------------------
# 4. Resolve the CLI list and smoke-check it
# -----------------------------------------------------------------------------
RELEASE_CLIS=()
if [[ ${#CLI_OVERRIDE[@]} -gt 0 ]]; then
    RELEASE_CLIS=("${CLI_OVERRIDE[@]}")
elif [[ "$CLI_PROFILE" == "all" ]]; then
    known="$(current_cli_targets)"
    [[ -n "$known" ]] || fail "could not read any lsdtt_cli() targets from CMakeLists.txt"
    skipped=()
    for exe in "$BUILD_DIR"/lsdtt-* "$BUILD_DIR"/bin/lsdtt-*; do
        [[ -f "$exe" && -x "$exe" ]] || continue
        name="$(basename "$exe")"
        [[ "$name" == *.* ]] && continue
        [[ "$name" == lsdtt-test-* ]] && continue
        [[ " ${RELEASE_CLIS[*]:-} " == *" $name "* ]] && continue
        # Stale binary: present in the build tree, no longer a CMake target.
        if ! grep -qx "$name" <<<"$known"; then
            [[ " ${skipped[*]:-} " == *" $name "* ]] || skipped+=("$name")
            continue
        fi
        RELEASE_CLIS+=("$name")
    done
    if [[ ${#skipped[@]} -gt 0 ]]; then
        info "WARNING: skipping ${#skipped[@]} stale binaries in $BUILD_DIR that are
no longer lsdtt_cli() targets in CMakeLists.txt: ${skipped[*]}
Delete them, or build into a clean directory, to silence this."
    fi
else
    while read -r profile name; do
        [[ "$profile" == "$CLI_PROFILE" ]] && RELEASE_CLIS+=("$name")
    done < <(grep -Ev '^[[:space:]]*(#|$)' tools/cli_profiles.txt)
fi
[[ ${#RELEASE_CLIS[@]} -gt 0 ]] || fail "no CLIs resolved for profile '$CLI_PROFILE'"

for cli in "${RELEASE_CLIS[@]}"; do
    exe="$BUILD_DIR/$cli"
    [[ -x "$exe" ]] || exe="$BUILD_DIR/bin/$cli"
    [[ -x "$exe" ]] || fail "release CLI not found in build dir: $cli"
    # -v both exercises startup (shared libraries, GDAL initialisation) and
    # reports the GDAL/PROJ actually loaded. Checking only the exit status
    # throws away the half that catches a binary linked against one GDAL and
    # running against another, so read the output.
    if ! v_out="$("$exe" -v 2>&1)"; then
        fail "smoke check failed: $cli -v returned nonzero
$v_out"
    fi
    [[ "$v_out" == *GDAL* ]] || fail \
        "smoke check failed: $cli -v did not report a GDAL version
$v_out"
done
info "smoke check passed for ${#RELEASE_CLIS[@]} CLIs (profile: $CLI_PROFILE)"

# The wheel matrix (python-wheels.yml) runs only after this script has passed
# in release.yml's first job, so this is the cheap early check that the
# bindings still compile and pass their tests against this C++. Exit 77 (no
# usable venv) is a failure here, not a skip.
smoke_args=(--repo-root "$REPO_ROOT" --work-dir "$BUILD_DIR/python_smoke"
            --python "$SMOKE_PYTHON")
[[ "$SKIP_BUILD" -eq 1 ]] && smoke_args+=(--reuse)
info "Python bindings smoke check (whole python_bindings/tests suite)"
bash tests/scripts/check_python_bindings_smoke.sh "${smoke_args[@]}" || fail \
    "the Python bindings do not build or their tests fail against this C++ (above).
Fix the bindings before releasing: whoever broke them fixes them."

# -----------------------------------------------------------------------------
# 5. Package
# -----------------------------------------------------------------------------
if [[ ${#CLI_OVERRIDE[@]} -gt 0 ]]; then
    cli_args=()
    for cli in "${CLI_OVERRIDE[@]}"; do
        cli_args+=(--cli "$cli")
    done
    tools/package_linux_backend.sh \
        --version "$VERSION" \
        --build-dir "$BUILD_DIR" \
        "${cli_args[@]}"
else
    tools/package_linux_backend.sh \
        --version "$VERSION" \
        --build-dir "$BUILD_DIR" \
        --cli-profile "$CLI_PROFILE"
fi

bundle_label="$CLI_PROFILE"
[[ ${#CLI_OVERRIDE[@]} -gt 0 ]] && bundle_label="custom"
bundle="$(ls dist/lsdtt3-backend-linux-*-"$bundle_label"-v"$VERSION".tar.gz 2>/dev/null | head -1)"
[[ -n "$bundle" ]] || fail "packaging reported success but no bundle found in dist/"
info "bundle: $bundle"

# -----------------------------------------------------------------------------
# 6. Publish (optional)
# -----------------------------------------------------------------------------
if [[ "$PUBLISH" -eq 1 ]]; then
    command -v gh >/dev/null 2>&1 || fail "--publish requires the gh CLI"
    release_tag="v$VERSION"
    if gh release view "$release_tag" >/dev/null 2>&1; then
        info "release $release_tag exists; uploading bundle (clobber)"
        gh release upload "$release_tag" "$bundle" --clobber
    else
        create_args=(--generate-notes)
        [[ "$PRERELEASE" -eq 1 ]] && create_args+=(--prerelease)
        info "creating release $release_tag"
        gh release create "$release_tag" "${create_args[@]}" "$bundle"
    fi
    info "published: $release_tag"
else
    info "not publishing (pass --publish to create/update the GitHub release)"
fi

info "done"
