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
  --cli-profile <name>         CLI set: qgis / core / all (default: core).
                               Defined in tools/cli_profiles.txt; all = every
                               lsdtt-* executable in the build directory
  --profile <name>             Manifest profile label (default: <cli-profile>-linux)
  --plugin-version-range <rng> Plugin compatibility range (default: >=0.1.0 <0.3.0)
  --cli <name>                 Package exactly this CLI (repeatable; overrides
                               --cli-profile, bundle is labelled 'custom')
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
CLI_PROFILE="core"
PROFILE=""
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
        --cli-profile)
            CLI_PROFILE="${2:-}"
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

# Resolve the CLI list: explicit --cli wins; otherwise the profile decides.
# qgis / core come from cli_profiles.txt; all = every lsdtt-* executable in
# the build directory (resolved after the build dir is validated below).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ ${#CLI_NAMES[@]} -gt 0 ]]; then
    CLI_PROFILE="custom"
elif [[ "$CLI_PROFILE" == "qgis" || "$CLI_PROFILE" == "core" ]]; then
    PROFILES_FILE="$SCRIPT_DIR/cli_profiles.txt"
    if [[ ! -f "$PROFILES_FILE" ]]; then
        echo "ERROR: profiles file not found: $PROFILES_FILE" >&2
        exit 1
    fi
    while read -r profile name; do
        if [[ "$profile" == "$CLI_PROFILE" ]]; then
            CLI_NAMES+=("$name")
        fi
    done < <(grep -Ev '^[[:space:]]*(#|$)' "$PROFILES_FILE")
    if [[ ${#CLI_NAMES[@]} -eq 0 ]]; then
        echo "ERROR: profile '$CLI_PROFILE' has no entries in $PROFILES_FILE" >&2
        exit 1
    fi
elif [[ "$CLI_PROFILE" != "all" ]]; then
    echo "ERROR: unknown --cli-profile '$CLI_PROFILE' (use qgis, core, or all)" >&2
    exit 1
fi

if [[ -z "$PROFILE" ]]; then
    PROFILE="${CLI_PROFILE}-linux"
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

# 'all' profile: discover every built lsdtt-* executable (skip test harnesses).
if [[ "$CLI_PROFILE" == "all" ]]; then
    declare -A seen_cli
    for exe in "$BUILD_ABS"/lsdtt-* "$BUILD_ABS"/bin/lsdtt-*; do
        [[ -f "$exe" && -x "$exe" ]] || continue
        name="$(basename "$exe")"
        [[ "$name" == *.* ]] && continue
        [[ "$name" == lsdtt-test-* ]] && continue
        [[ -n "${seen_cli[$name]:-}" ]] && continue
        seen_cli[$name]=1
        CLI_NAMES+=("$name")
    done
    if [[ ${#CLI_NAMES[@]} -eq 0 ]]; then
        echo "ERROR: --cli-profile all found no lsdtt-* executables in $BUILD_ABS" >&2
        exit 1
    fi
    echo "Profile 'all': packaging ${#CLI_NAMES[@]} CLIs discovered in the build directory"
fi

for cmd in ldd gdal-config tar python3 objdump; do
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

# Emits "requested_name<TAB>absolute_path" per resolved dependency. The
# requested name (left of '=>') is the SONAME the dynamic loader searches
# for; the staged file must carry exactly that name. Staging real files
# under their SONAME instead of copying symlink chains keeps the bundle
# working on filesystems and archivers that cannot represent symlinks.
extract_ldd_deps() {
    local item="$1"
    ldd "$item" 2>/dev/null | while IFS= read -r line; do
        if [[ "$line" == *"=>"* ]]; then
            local name p
            name="$(awk '{print $1}' <<<"$line")"
            p="$(awk -F'=> ' '{print $2}' <<<"$line" | awk '{print $1}')"
            if [[ "$p" == /* ]]; then
                printf '%s\t%s\n' "$name" "$p"
            fi
        else
            local p
            p="$(awk '{print $1}' <<<"$line")"
            if [[ "$p" == /* ]]; then
                printf '%s\t%s\n' "$(basename "$p")" "$p"
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

BUNDLE_NAME="lsdtt3-backend-linux-${ARCH}-${CLI_PROFILE}-v${VERSION}"
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
    while IFS= read -r entry; do
        [[ -z "$entry" ]] && continue
        QUEUE+=("$entry")
    done < <(extract_ldd_deps "$exe")
done

while [[ ${#QUEUE[@]} -gt 0 ]]; do
    entry="${QUEUE[0]}"
    QUEUE=("${QUEUE[@]:1}")

    name="${entry%%$'\t'*}"
    dep="${entry#*$'\t'}"

    if [[ -n "${SEEN_LIBS[$name]:-}" ]]; then
        continue
    fi
    SEEN_LIBS[$name]=1

    if is_core_glibc_lib "$name"; then
        continue
    fi

    dep="$(realpath "$dep")"
    [[ -e "$dep" ]] || continue

    if [[ ! -e "$STAGE_LIB/$name" ]]; then
        cp "$dep" "$STAGE_LIB/$name"
    fi

    while IFS= read -r next_entry; do
        [[ -z "$next_entry" ]] && continue
        QUEUE+=("$next_entry")
    done < <(extract_ldd_deps "$dep")
done

if command -v patchelf >/dev/null 2>&1; then
    for exe in "$STAGE_BIN"/lsdtt-*; do
        if [[ -f "$exe" ]]; then
            patchelf --set-rpath '$ORIGIN/../lib' "$exe"
        fi
    done
fi

# Verify the bundle is self-contained: every DT_NEEDED of every staged
# binary and library must resolve inside the stage (or be core glibc).
# Fail the packaging rather than ship a bundle that breaks for users.
MISSING=0
for f in "$STAGE_BIN"/lsdtt-* "$STAGE_LIB"/*.so*; do
    [[ -f "$f" ]] || continue
    while IFS= read -r needed; do
        [[ -z "$needed" ]] && continue
        if is_core_glibc_lib "$needed"; then
            continue
        fi
        if [[ ! -e "$STAGE_LIB/$needed" ]]; then
            echo "ERROR: $(basename "$f") needs '$needed' but it is not staged in lib/" >&2
            MISSING=1
        fi
    done < <(objdump -p "$f" 2>/dev/null | awk '/NEEDED/{print $2}')
done
if [[ $MISSING -ne 0 ]]; then
    echo "ERROR: bundle is not self-contained; aborting" >&2
    exit 1
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

# End-user README: the audience is often non-technical, so it must cover the
# full path from download to a working command, including Google Colab where
# `!source` does not work (each ! line is a separate shell).
cat > "$STAGE_DIR/README.md" <<'EOF'
# lsdtt3 backend for Linux (v@VERSION@)

Self-contained LSDTopoTools command-line tools for Linux (x86_64). All
required libraries and the GDAL/PROJ support data are included -- there is
nothing to install and no admin rights are needed. Just unpack and run.

## Install and run (Linux terminal)

1. Download `@BUNDLE_NAME@.tar.gz` and open a terminal in the folder
   containing it.
2. Unpack it:

   ```bash
   tar -xzf @BUNDLE_NAME@.tar.gz
   ```

3. Activate the tools for the current terminal session (note the `source` --
   running the script with `bash` will NOT work):

   ```bash
   source @BUNDLE_NAME@/setenv.sh
   ```

4. Check that it works:

   ```bash
   lsdtt-raster-preprocessing --help
   ```

You need to repeat step 3 in every new terminal. To make it permanent, add
that `source` line (with the full path to `setenv.sh`) to your `~/.bashrc`.

## Google Colab

In Colab, download and unpack in a code cell (the leading `!` runs a shell
command):

```python
!wget https://www.geos.ed.ac.uk/~smudd/lsdtt_packages/@BUNDLE_NAME@.tar.gz
!tar -xzf @BUNDLE_NAME@.tar.gz
```

Do NOT use `!source setenv.sh` -- every `!` line runs in its own throwaway
shell, so nothing it sets survives to the next line. Instead, set the
environment from Python once per session:

```python
import os

root = "/content/@BUNDLE_NAME@"
os.environ["PATH"] = f"{root}/bin:" + os.environ["PATH"]
os.environ["LD_LIBRARY_PATH"] = f"{root}/lib:" + os.environ.get("LD_LIBRARY_PATH", "")
os.environ["GDAL_DATA"] = f"{root}/share/gdal"
os.environ["PROJ_DATA"] = f"{root}/share/proj"
os.environ["PROJ_LIB"]  = os.environ["PROJ_DATA"]
```

After that cell, the tools work in any later cell:

```python
!lsdtt-raster-preprocessing --help
```

## What is in this bundle

- `bin/` -- the lsdtt command-line tools:
@TOOL_LIST@
- `lib/` -- bundled shared libraries (GDAL, PROJ, ...)
- `share/gdal`, `share/proj` -- support data the libraries need at runtime
- `setenv.sh` -- sets up PATH and library/data paths (must be `source`d)
- `manifest.json` -- version info and SHA256 checksums of every file
EOF

TOOL_LIST_MD=""
for cli in "${CLI_NAMES[@]}"; do
    TOOL_LIST_MD+="  - \`$cli\`\n"
done
sed -i \
    -e "s/@VERSION@/$VERSION/g" \
    -e "s/@BUNDLE_NAME@/$BUNDLE_NAME/g" \
    -e "s|@TOOL_LIST@|${TOOL_LIST_MD%\\n}|" \
    "$STAGE_DIR/README.md"

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
