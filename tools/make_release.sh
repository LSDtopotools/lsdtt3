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
#   3. Configure and build (Release; OpenCV + Eigen on).
#   4. Smoke-check every packaged CLI (runs each with -v).
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

while [[ $# -gt 0 ]]; do
    case "$1" in
        --version)    VERSION="${2:-}"; shift 2 ;;
        --profile)    CLI_PROFILE="${2:-}"; shift 2 ;;
        --build-dir)  BUILD_DIR="${2:-}"; shift 2 ;;
        --skip-build) SKIP_BUILD=1; shift ;;
        --publish)    PUBLISH=1; shift ;;
        --prerelease) PRERELEASE=1; shift ;;
        --cli)        CLI_OVERRIDE+=("${2:-}"); shift 2 ;;
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

if ! git diff --quiet || ! git diff --cached --quiet; then
    info "WARNING: working tree has uncommitted changes; binaries will report -dirty"
fi

# -----------------------------------------------------------------------------
# 3. Configure and build
# -----------------------------------------------------------------------------
if [[ "$SKIP_BUILD" -eq 1 ]]; then
    [[ -d "$BUILD_DIR" ]] || fail "--skip-build given but build dir not found: $BUILD_DIR"
    info "skipping build; using existing $BUILD_DIR"
else
    info "configuring $BUILD_DIR"
    cmake -B "$BUILD_DIR" \
        -DCMAKE_BUILD_TYPE=Release \
        -DLSDTT_OPENCV=ON \
        -DLSDTT_EIGEN=ON
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
    for exe in "$BUILD_DIR"/lsdtt-* "$BUILD_DIR"/bin/lsdtt-*; do
        [[ -f "$exe" && -x "$exe" ]] || continue
        name="$(basename "$exe")"
        [[ "$name" == *.* ]] && continue
        [[ "$name" == lsdtt-test-* ]] && continue
        [[ " ${RELEASE_CLIS[*]:-} " == *" $name "* ]] && continue
        RELEASE_CLIS+=("$name")
    done
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
    "$exe" -v >/dev/null || fail "smoke check failed: $cli -v returned nonzero"
done
info "smoke check passed for ${#RELEASE_CLIS[@]} CLIs (profile: $CLI_PROFILE)"

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
