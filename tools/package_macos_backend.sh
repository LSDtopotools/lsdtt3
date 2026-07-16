#!/usr/bin/env bash
#
# Package a minimal redistributable macOS backend bundle for lsdtt3.
#
# STATUS: UNTESTED. Written by analogy with package_linux_backend.sh on a
# machine without macOS. The self-containment check at the end fails hard,
# so the first run on a real Mac will report exactly what (if anything) is
# still wrong rather than shipping a broken bundle.
#
# macOS differences from the Linux packaging:
# - Dependencies are discovered with `otool -L` (not ldd).
# - Mach-O binaries reference libraries by full "install name"; every
#   staged binary is rewritten with `install_name_tool` so references
#   become @rpath/<libname>, and executables get an rpath of
#   @executable_path/../lib. No LD_LIBRARY_PATH-style env var is needed
#   at runtime (and SIP strips DYLD_* vars anyway).
# - Modifying a binary invalidates its code signature, which Apple Silicon
#   refuses to run, so every staged file is re-signed ad-hoc (codesign -s -).
# - Downloaded bundles carry the com.apple.quarantine attribute; the README
#   tells users how to clear it.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  tools/package_macos_backend.sh --version <X.Y.Z> [options]

Options:
  --version <ver>              Backend version (required)
  --build-dir <path>           Build directory (default: build)
  --out-dir <path>             Output directory (default: dist)
  --arch <name>                Architecture label (default: uname -m)
  --profile <name>             Profile label (default: minimal-macos)
  --plugin-version-range <rng> Plugin compatibility range (default: >=0.1.0 <0.3.0)
  --cli <name>                 Add CLI name to package (repeatable)
  --help                       Show this help

Examples:
  tools/package_macos_backend.sh --version 0.2.0
  tools/package_macos_backend.sh --version 0.2.0 --build-dir build/minimal-macos
EOF
}

VERSION=""
BUILD_DIR="build"
OUT_DIR="dist"
ARCH="$(uname -m)"
PROFILE="minimal-macos"
PLUGIN_VERSION_RANGE=">=0.1.0 <0.3.0"
CLI_NAMES=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --version)
            VERSION="${2:-}"
            shift 2
            ;;
        --build-dir)
            BUILD_DIR="${2:-}"
            shift 2
            ;;
        --out-dir)
            OUT_DIR="${2:-}"
            shift 2
            ;;
        --arch)
            ARCH="${2:-}"
            shift 2
            ;;
        --profile)
            PROFILE="${2:-}"
            shift 2
            ;;
        --plugin-version-range)
            PLUGIN_VERSION_RANGE="${2:-}"
            shift 2
            ;;
        --cli)
            CLI_NAMES+=("${2:-}")
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$VERSION" ]]; then
    echo "ERROR: --version is required" >&2
    usage
    exit 1
fi

if [[ ${#CLI_NAMES[@]} -eq 0 ]]; then
    CLI_NAMES=(
        "lsdtt-raster-preprocessing"
        "lsdtt-channel-extraction"
        "lsdtt-chi-analysis"
        "lsdtt-basin-extraction"
    )
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# macOS ships no realpath on older systems; use python3 everywhere for it.
realpath_py() {
    python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "$1"
}

if [[ "$BUILD_DIR" = /* ]]; then
    BUILD_ABS="$(realpath_py "$BUILD_DIR")"
else
    BUILD_ABS="$(realpath_py "$REPO_ROOT/$BUILD_DIR")"
fi

if [[ "$OUT_DIR" = /* ]]; then
    OUT_ABS="$(realpath_py "$OUT_DIR")"
else
    OUT_ABS="$(realpath_py "$REPO_ROOT/$OUT_DIR")"
fi

if [[ ! -d "$BUILD_ABS" ]]; then
    echo "ERROR: build directory not found: $BUILD_ABS" >&2
    exit 1
fi

for cmd in otool install_name_tool codesign gdal-config tar python3; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $cmd" >&2
        exit 1
    fi
done

resolve_exe_path() {
    local base_dir="$1"
    local cli_name="$2"
    local c1="$base_dir/$cli_name"
    local c2="$base_dir/bin/$cli_name"
    if [[ -x "$c1" ]]; then
        realpath_py "$c1"
        return
    fi
    if [[ -x "$c2" ]]; then
        realpath_py "$c2"
        return
    fi
    echo "ERROR: could not find executable '$cli_name' in $base_dir (checked root and bin/)" >&2
    exit 1
}

# OS-provided libraries must not be bundled (and cannot be: they live in the
# dyld shared cache, and copying them out breaks them).
is_system_lib() {
    case "$1" in
        /usr/lib/*|/System/*)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

# Emits the raw install-name references of a Mach-O file, one per line,
# excluding the file's own install-name id line that otool -L prints first
# for dylibs.
mach_o_deps() {
    local item="$1"
    local self_id
    self_id="$(otool -D "$item" 2>/dev/null | tail -n +2 | head -n 1 || true)"
    otool -L "$item" 2>/dev/null | tail -n +2 | awk '{print $1}' | while IFS= read -r dep; do
        [[ -z "$dep" ]] && continue
        [[ "$dep" == "$self_id" ]] && continue
        echo "$dep"
    done
}

# Resolves an install-name reference to an on-disk file. Absolute paths are
# used as-is; @rpath/@loader_path references are tried relative to the
# directory of the referencing binary (covers the common Homebrew/vcpkg
# layouts where dependent dylibs sit side by side).
resolve_dep_path() {
    local ref_dir="$1"
    local dep="$2"
    case "$dep" in
        /*)
            echo "$dep"
            ;;
        @rpath/*|@loader_path/*)
            local rel="${dep#@*/}"
            if [[ -e "$ref_dir/$rel" ]]; then
                echo "$ref_dir/$rel"
            else
                echo ""
            fi
            ;;
        *)
            echo ""
            ;;
    esac
}

BUNDLE_NAME="lsdtt3-backend-macos-${ARCH}-minimal-v${VERSION}"
STAGE_DIR="$OUT_ABS/$BUNDLE_NAME"
STAGE_BIN="$STAGE_DIR/bin"
STAGE_LIB="$STAGE_DIR/lib"
STAGE_SHARE="$STAGE_DIR/share"

rm -rf "$STAGE_DIR"
mkdir -p "$STAGE_BIN" "$STAGE_LIB" "$STAGE_SHARE"

EXE_PATHS=()
for cli in "${CLI_NAMES[@]}"; do
    exe_path="$(resolve_exe_path "$BUILD_ABS" "$cli")"
    cp "$exe_path" "$STAGE_BIN/"
    EXE_PATHS+=("$exe_path")
done

# ---- Transitive dependency closure -----------------------------------------
# Queue entries are "referenced_name<TAB>resolved_path". Each library is
# staged as a single regular file under the basename of its install-name
# reference (no symlinks), mirroring the Linux script's SONAME staging.
declare -a QUEUE=()
declare -a SEEN_NAMES=()

seen_name() {
    local n="$1"
    local s
    for s in ${SEEN_NAMES[@]+"${SEEN_NAMES[@]}"}; do
        [[ "$s" == "$n" ]] && return 0
    done
    return 1
}

for exe in "${EXE_PATHS[@]}"; do
    exe_dir="$(dirname "$exe")"
    while IFS= read -r dep; do
        [[ -z "$dep" ]] && continue
        is_system_lib "$dep" && continue
        resolved="$(resolve_dep_path "$exe_dir" "$dep")"
        if [[ -z "$resolved" || ! -e "$resolved" ]]; then
            echo "ERROR: cannot resolve dependency '$dep' of $(basename "$exe")" >&2
            exit 1
        fi
        QUEUE+=("$(basename "$dep")"$'\t'"$resolved")
    done < <(mach_o_deps "$exe")
done

while [[ ${#QUEUE[@]} -gt 0 ]]; do
    entry="${QUEUE[0]}"
    QUEUE=(${QUEUE[@]+"${QUEUE[@]:1}"})

    name="${entry%%$'\t'*}"
    dep_path="${entry#*$'\t'}"

    if seen_name "$name"; then
        continue
    fi
    SEEN_NAMES+=("$name")

    dep_real="$(realpath_py "$dep_path")"
    [[ -e "$dep_real" ]] || continue

    if [[ ! -e "$STAGE_LIB/$name" ]]; then
        cp "$dep_real" "$STAGE_LIB/$name"
        chmod u+w "$STAGE_LIB/$name"
    fi

    dep_dir="$(dirname "$dep_real")"
    while IFS= read -r next_dep; do
        [[ -z "$next_dep" ]] && continue
        is_system_lib "$next_dep" && continue
        next_resolved="$(resolve_dep_path "$dep_dir" "$next_dep")"
        if [[ -z "$next_resolved" || ! -e "$next_resolved" ]]; then
            echo "ERROR: cannot resolve dependency '$next_dep' of $name" >&2
            exit 1
        fi
        QUEUE+=("$(basename "$next_dep")"$'\t'"$next_resolved")
    done < <(mach_o_deps "$dep_real")
done

# ---- Rewrite install names so the bundle is relocatable ---------------------
for exe in "$STAGE_BIN"/lsdtt-*; do
    [[ -f "$exe" ]] || continue
    install_name_tool -add_rpath '@executable_path/../lib' "$exe" 2>/dev/null || true
    while IFS= read -r dep; do
        [[ -z "$dep" ]] && continue
        is_system_lib "$dep" && continue
        install_name_tool -change "$dep" "@rpath/$(basename "$dep")" "$exe"
    done < <(mach_o_deps "$exe")
done

for lib in "$STAGE_LIB"/*; do
    [[ -f "$lib" ]] || continue
    install_name_tool -id "@rpath/$(basename "$lib")" "$lib"
    install_name_tool -add_rpath '@loader_path' "$lib" 2>/dev/null || true
    while IFS= read -r dep; do
        [[ -z "$dep" ]] && continue
        is_system_lib "$dep" && continue
        case "$dep" in
            @rpath/*) continue ;;
        esac
        install_name_tool -change "$dep" "@rpath/$(basename "$dep")" "$lib"
    done < <(mach_o_deps "$lib")
done

# install_name_tool invalidates signatures; Apple Silicon refuses unsigned
# binaries, so re-sign everything ad-hoc.
for f in "$STAGE_BIN"/lsdtt-* "$STAGE_LIB"/*; do
    [[ -f "$f" ]] || continue
    codesign --force --sign - "$f"
done

# ---- Verify the bundle is self-contained ------------------------------------
# Every reference of every staged binary must be a system library or an
# @rpath reference satisfied inside lib/. Fail rather than ship broken.
MISSING=0
for f in "$STAGE_BIN"/lsdtt-* "$STAGE_LIB"/*; do
    [[ -f "$f" ]] || continue
    while IFS= read -r dep; do
        [[ -z "$dep" ]] && continue
        is_system_lib "$dep" && continue
        case "$dep" in
            @rpath/*)
                if [[ ! -e "$STAGE_LIB/$(basename "$dep")" ]]; then
                    echo "ERROR: $(basename "$f") needs '$dep' but it is not staged in lib/" >&2
                    MISSING=1
                fi
                ;;
            *)
                echo "ERROR: $(basename "$f") still references '$dep' (not rewritten to @rpath)" >&2
                MISSING=1
                ;;
        esac
    done < <(mach_o_deps "$f")
done
if [[ $MISSING -ne 0 ]]; then
    echo "ERROR: bundle is not self-contained; aborting" >&2
    exit 1
fi

# ---- GDAL / PROJ data --------------------------------------------------------
find_proj_data_dir() {
    if [[ -n "${PROJ_LIB:-}" && -f "${PROJ_LIB}/proj.db" ]]; then
        echo "$PROJ_LIB"
        return
    fi
    if [[ -n "${PROJ_DATA:-}" && -f "${PROJ_DATA}/proj.db" ]]; then
        echo "$PROJ_DATA"
        return
    fi

    local candidates=(
        "/usr/local/share/proj"
        "/opt/homebrew/share/proj"
    )

    if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists proj; then
        local d
        d="$(pkg-config --variable=datadir proj)/proj"
        candidates+=("$d")
    fi

    for d in "${candidates[@]}"; do
        if [[ -f "$d/proj.db" ]]; then
            echo "$d"
            return
        fi
    done

    echo ""
}

GDAL_DATA_DIR="$(gdal-config --datadir)"
if [[ ! -d "$GDAL_DATA_DIR" ]]; then
    echo "ERROR: GDAL data directory not found: $GDAL_DATA_DIR" >&2
    exit 1
fi

PROJ_DATA_DIR="$(find_proj_data_dir)"
if [[ -z "$PROJ_DATA_DIR" ]]; then
    echo "ERROR: could not locate PROJ data directory (proj.db)" >&2
    exit 1
fi

cp -R "$GDAL_DATA_DIR" "$STAGE_SHARE/gdal"
cp -R "$PROJ_DATA_DIR" "$STAGE_SHARE/proj"

if [[ -f "$REPO_ROOT/LICENSE" ]]; then
    cp "$REPO_ROOT/LICENSE" "$STAGE_DIR/LICENSE"
fi

{
    echo "This backend bundle includes runtime shared libraries from the build host."
    echo "Profile: $PROFILE"
    echo "Architecture: $ARCH"
    echo "Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo
    echo "Included executables:"
    for cli in "${CLI_NAMES[@]}"; do
        echo "- $cli"
    done
    echo
    echo "Included shared libraries:"
    (cd "$STAGE_LIB" && ls -1 | sort | sed 's/^/- /')
} > "$STAGE_DIR/THIRD_PARTY_NOTICES.txt"

# setenv.sh only needs PATH and the data dirs: library lookup is baked into
# the binaries via @rpath (DYLD_* vars are stripped by SIP anyway).
cat > "$STAGE_DIR/setenv.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$ROOT_DIR/bin:${PATH:-}"
export GDAL_DATA="$ROOT_DIR/share/gdal"
export PROJ_DATA="$ROOT_DIR/share/proj"
export PROJ_LIB="$PROJ_DATA"
EOF
chmod +x "$STAGE_DIR/setenv.sh"

# End-user README: audience is often non-technical, so cover download to a
# working command, including the Gatekeeper quarantine step.
cat > "$STAGE_DIR/README.md" <<'EOF'
# lsdtt3 backend for macOS (v@VERSION@, @ARCH@)

Self-contained LSDTopoTools command-line tools for macOS. All required
libraries and the GDAL/PROJ support data are included -- there is nothing to
install and no admin rights are needed. Just unpack and run.

## Install and run (Terminal)

1. Download `@BUNDLE_NAME@.tar.gz` and open Terminal in the folder
   containing it (in Finder: right-click the folder, then
   "New Terminal at Folder").
2. Unpack it:

   ```bash
   tar -xzf @BUNDLE_NAME@.tar.gz
   ```

3. Tell macOS the downloaded tools are OK to run. Files downloaded with a
   browser are "quarantined" and macOS will otherwise block them because
   they are not from the App Store:

   ```bash
   xattr -dr com.apple.quarantine @BUNDLE_NAME@
   ```

4. Activate the tools for the current Terminal session (note the `source` --
   running the script with `bash` will NOT work):

   ```bash
   source @BUNDLE_NAME@/setenv.sh
   ```

5. Check that it works:

   ```bash
   lsdtt-raster-preprocessing --help
   ```

You need to repeat step 4 in every new Terminal window. To make it
permanent, add that `source` line (with the full path to `setenv.sh`) to
your `~/.zshrc`.

## What is in this bundle

- `bin/` -- the lsdtt command-line tools:
@TOOL_LIST@
- `lib/` -- bundled shared libraries (GDAL, PROJ, ...)
- `share/gdal`, `share/proj` -- support data the libraries need at runtime
- `setenv.sh` -- sets up PATH and data paths (must be `source`d)
- `manifest.json` -- version info and SHA256 checksums of every file

Note: this bundle only runs on @ARCH@ Macs (check yours with `uname -m` --
Apple Silicon reports `arm64`, older Intel Macs report `x86_64`).
EOF

TOOL_LIST_MD=""
for cli in "${CLI_NAMES[@]}"; do
    TOOL_LIST_MD+="  - \`$cli\`"$'\n'
done
# BSD sed (macOS) needs -i '' and does not expand \n in replacements the way
# GNU sed does; do the substitution with python3 instead for portability.
TOOL_LIST_MD="$TOOL_LIST_MD" \
python3 - "$STAGE_DIR/README.md" "$VERSION" "$BUNDLE_NAME" "$ARCH" <<'PY'
import os
import sys
path, version, bundle, arch = sys.argv[1:5]
tools = os.environ["TOOL_LIST_MD"].rstrip("\n")
text = open(path).read()
text = (text.replace("@VERSION@", version)
            .replace("@BUNDLE_NAME@", bundle)
            .replace("@ARCH@", arch)
            .replace("@TOOL_LIST@", tools))
open(path, "w").write(text)
PY

python3 - <<PY
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone

stage = Path(r"$STAGE_DIR")
manifest = {
    "backend_version": "$VERSION",
    "profile": "$PROFILE",
    "architecture": "$ARCH",
    "plugin_version_range": "$PLUGIN_VERSION_RANGE",
    "generated_utc": datetime.now(timezone.utc).isoformat(),
    "tools": [],
    "checksums_sha256": {}
}

# Fill tools from staged binaries to avoid shell quoting issues.
manifest["tools"] = sorted(p.name for p in (stage / "bin").glob("lsdtt-*"))

for path in sorted(stage.rglob("*")):
    if path.is_file():
        rel = path.relative_to(stage).as_posix()
        h = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest["checksums_sha256"][rel] = h

(stage / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
PY

mkdir -p "$OUT_ABS"
TARBALL="$OUT_ABS/$BUNDLE_NAME.tar.gz"
rm -f "$TARBALL"
tar -C "$OUT_ABS" -czf "$TARBALL" "$BUNDLE_NAME"

echo "Created backend bundle: $TARBALL"
echo "Staging directory: $STAGE_DIR"
echo
echo "REMINDER: this script is untested on real macOS hardware."
echo "If the self-containment check passed, test at minimum:"
echo "  source $BUNDLE_NAME/setenv.sh && lsdtt-raster-preprocessing --help"
