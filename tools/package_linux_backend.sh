#!/usr/bin/env bash
#
# Package a minimal redistributable Linux backend bundle for lsdtt3.
#
# The bundle includes selected CLIs, transitive runtime shared libraries,
# GDAL/PROJ data directories, metadata, and a tar.gz artifact.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  tools/package_linux_backend.sh --version <X.Y.Z> [options]

Options:
  --version <ver>              Backend version (required)
  --build-dir <path>           Build directory (default: build)
  --out-dir <path>             Output directory (default: dist)
  --arch <name>                Architecture label (default: x86_64)
  --profile <name>             Profile label (default: minimal-linux)
  --plugin-version-range <rng> Plugin compatibility range (default: >=0.1.0 <0.3.0)
  --cli <name>                 Add CLI name to package (repeatable)
  --help                       Show this help

Examples:
  tools/package_linux_backend.sh --version 0.2.0
  tools/package_linux_backend.sh --version 0.2.0 --build-dir build/minimal-linux
EOF
}

VERSION=""
BUILD_DIR="build"
OUT_DIR="dist"
ARCH="x86_64"
PROFILE="minimal-linux"
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

if [[ "$BUILD_DIR" = /* ]]; then
    BUILD_ABS="$(realpath "$BUILD_DIR")"
else
    BUILD_ABS="$(realpath "$REPO_ROOT/$BUILD_DIR")"
fi

if [[ "$OUT_DIR" = /* ]]; then
    OUT_ABS="$(realpath -m "$OUT_DIR")"
else
    OUT_ABS="$(realpath -m "$REPO_ROOT/$OUT_DIR")"
fi

if [[ ! -d "$BUILD_ABS" ]]; then
    echo "ERROR: build directory not found: $BUILD_ABS" >&2
    exit 1
fi

for cmd in ldd gdal-config tar python3; do
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
        realpath "$c1"
        return
    fi
    if [[ -x "$c2" ]]; then
        realpath "$c2"
        return
    fi
    echo "ERROR: could not find executable '$cli_name' in $base_dir (checked root and bin/)" >&2
    exit 1
}

is_core_glibc_lib() {
    local name="$1"
    case "$name" in
        ld-linux*.so*|libc.so.*|libm.so.*|libpthread.so.*|libdl.so.*|librt.so.*)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

copy_lib_with_symlink_chain() {
    local src="$1"
    local dst_dir="$2"

    if [[ ! -e "$src" ]]; then
        return
    fi

    local base
    base="$(basename "$src")"
    if [[ ! -e "$dst_dir/$base" ]]; then
        cp -a "$src" "$dst_dir/"
    fi

    if [[ -L "$src" ]]; then
        local target
        target="$(readlink "$src")"
        if [[ "$target" != /* ]]; then
            target="$(dirname "$src")/$target"
        fi
        target="$(realpath "$target")"
        if [[ -e "$target" ]]; then
            copy_lib_with_symlink_chain "$target" "$dst_dir"
        fi
    fi
}

extract_ldd_paths() {
    local item="$1"
    ldd "$item" 2>/dev/null | while IFS= read -r line; do
        if [[ "$line" == *"=>"* ]]; then
            local p
            p="$(awk -F'=> ' '{print $2}' <<<"$line" | awk '{print $1}')"
            if [[ "$p" == /* ]]; then
                echo "$p"
            fi
        else
            local p
            p="$(awk '{print $1}' <<<"$line")"
            if [[ "$p" == /* ]]; then
                echo "$p"
            fi
        fi
    done
}

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
        "/usr/share/proj"
        "/usr/local/share/proj"
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

BUNDLE_NAME="lsdtt3-backend-linux-${ARCH}-minimal-v${VERSION}"
STAGE_DIR="$OUT_ABS/$BUNDLE_NAME"
STAGE_BIN="$STAGE_DIR/bin"
STAGE_LIB="$STAGE_DIR/lib"
STAGE_SHARE="$STAGE_DIR/share"

rm -rf "$STAGE_DIR"
mkdir -p "$STAGE_BIN" "$STAGE_LIB" "$STAGE_SHARE"

EXE_PATHS=()
for cli in "${CLI_NAMES[@]}"; do
    exe_path="$(resolve_exe_path "$BUILD_ABS" "$cli")"
    cp -a "$exe_path" "$STAGE_BIN/"
    EXE_PATHS+=("$exe_path")
done

declare -A SEEN_LIBS=()
QUEUE=()

for exe in "${EXE_PATHS[@]}"; do
    while IFS= read -r dep; do
        [[ -z "$dep" ]] && continue
        QUEUE+=("$dep")
    done < <(extract_ldd_paths "$exe")
done

while [[ ${#QUEUE[@]} -gt 0 ]]; do
    dep="${QUEUE[0]}"
    QUEUE=("${QUEUE[@]:1}")

    dep="$(realpath "$dep")"
    [[ -e "$dep" ]] || continue

    if [[ -n "${SEEN_LIBS[$dep]:-}" ]]; then
        continue
    fi
    SEEN_LIBS[$dep]=1

    dep_base="$(basename "$dep")"
    if is_core_glibc_lib "$dep_base"; then
        continue
    fi

    copy_lib_with_symlink_chain "$dep" "$STAGE_LIB"

    while IFS= read -r next_dep; do
        [[ -z "$next_dep" ]] && continue
        QUEUE+=("$next_dep")
    done < <(extract_ldd_paths "$dep")
done

if command -v patchelf >/dev/null 2>&1; then
    for exe in "$STAGE_BIN"/lsdtt-*; do
        if [[ -f "$exe" ]]; then
            patchelf --set-rpath '$ORIGIN/../lib' "$exe"
        fi
    done
fi

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

cp -a "$GDAL_DATA_DIR" "$STAGE_SHARE/gdal"
cp -a "$PROJ_DATA_DIR" "$STAGE_SHARE/proj"

if [[ -f "$REPO_ROOT/LICENSE" ]]; then
    cp -a "$REPO_ROOT/LICENSE" "$STAGE_DIR/LICENSE"
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
    find "$STAGE_LIB" -maxdepth 1 -type f -name '*.so*' -printf '%f\n' | sort | sed 's/^/- /'
} > "$STAGE_DIR/THIRD_PARTY_NOTICES.txt"

cat > "$STAGE_DIR/setenv.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$ROOT_DIR/bin:${PATH:-}"
export LD_LIBRARY_PATH="$ROOT_DIR/lib:${LD_LIBRARY_PATH:-}"
export GDAL_DATA="$ROOT_DIR/share/gdal"
export PROJ_DATA="$ROOT_DIR/share/proj"
export PROJ_LIB="$PROJ_DATA"
EOF
chmod +x "$STAGE_DIR/setenv.sh"

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
