#!/usr/bin/env bash
set -Eeuo pipefail

SRC_DIR="${1:-.}"
SRC_DIR="$(cd "$SRC_DIR" && pwd)"
OUT_DIR="$SRC_DIR/for_elm"
MANIFEST="$OUT_DIR/directory_structure.txt"

TMP_DIR="$(mktemp -d)"
REL_LIST="$TMP_DIR/relevant_files.txt"
MAP_LIST="$TMP_DIR/mapping.txt"
SKIP_LIST="$TMP_DIR/skipped.txt"

cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

log() {
    printf '%s\n' "$*" >&2
}

sanitize_path_to_name() {
    local rel="$1"
    local flat
    flat="${rel//\//__}"
    flat="${flat// /_}"
    flat="$(printf '%s' "$flat" | sed 's/[^A-Za-z0-9._-]/_/g')"
    flat="${flat#.}"
    [ -n "$flat" ] || flat="unnamed"
    printf '%s' "$flat"
}

has_source_like_name() {
    local path="$1"
    local base ext lower_ext

    base="$(basename "$path")"
    ext="${base##*.}"
    lower_ext="$(printf '%s' "$ext" | tr '[:upper:]' '[:lower:]')"

    case "$base" in
        CMakeLists.txt|Makefile|makefile|GNUmakefile|README|README.*|LICENSE|LICENSE.*|COPYING|COPYING.*|CHANGELOG|CHANGELOG.*)
            return 0
            ;;
    esac

    case "$lower_ext" in
        c|cc|cpp|cxx|h|hh|hpp|hxx|ipp|inl|tpp|ixx|md|txt|rst|cmake|sh|py|json|yaml|yml|toml|ini|cfg|conf|xml|html)
            return 0
            ;;
    esac

    return 1
}

should_ignore_path() {
    local rel="$1"
    case "$rel" in
        for_elm/*|*/for_elm/*) return 0 ;;
        .git/*|*/.git/*) return 0 ;;
        .svn/*|*/.svn/*) return 0 ;;
        .hg/*|*/.hg/*) return 0 ;;
        .idea/*|*/.idea/*) return 0 ;;
        .vscode/*|*/.vscode/*) return 0 ;;
        .vs/*|*/.vs/*) return 0 ;;
        __pycache__/*|*/__pycache__/*) return 0 ;;
        node_modules/*|*/node_modules/*) return 0 ;;
        dist/*|*/dist/*) return 0 ;;
        out/*|*/out/*) return 0 ;;
        build/*|*/build/*) return 0 ;;
        cmake-build-debug/*|*/cmake-build-debug/*) return 0 ;;
        cmake-build-release/*|*/cmake-build-release/*) return 0 ;;
        .venv/*|*/.venv/*) return 0 ;;
        venv/*|*/venv/*) return 0 ;;
    esac
    return 1
}

is_binary() {
    local f="$1"

    if command -v file >/dev/null 2>&1; then
        local enc
        enc="$(file --brief --mime-encoding "$f" 2>/dev/null || true)"
        enc="$(printf '%s' "$enc" | tr '[:upper:]' '[:lower:]')"

        case "$enc" in
            binary)
                return 0
                ;;
            us-ascii|ascii|utf-8|utf-16|utf-16le|utf-16be|utf-32|utf-32le|utf-32be|iso-8859-1|unknown-8bit|*)
                return 1
                ;;
        esac
    fi

    # Fallback only if file is unavailable:
    # treat as text unless proven otherwise
    return 1
}

detect_bom_encoding() {
    local f="$1"
    local hex
    hex="$(head -c 4 "$f" | od -An -tx1 | tr -d ' \n' | tr '[:upper:]' '[:lower:]')"
    case "$hex" in
        efbbbf*) echo "UTF-8" ;;
        fffe0000*) echo "UTF-32LE" ;;
        0000feff*) echo "UTF-32BE" ;;
        fffe*) echo "UTF-16LE" ;;
        feff*) echo "UTF-16BE" ;;
        *) echo "" ;;
    esac
}

convert_with_iconv() {
    local src="$1"
    local dest="$2"
    local enc="$3"
    iconv -f "$enc" -t UTF-8 "$src" -o "$dest" >/dev/null 2>&1
}

convert_with_python_fallback() {
    local src="$1"
    local dest="$2"

    command -v python3 >/dev/null 2>&1 || return 1

    python3 - "$src" "$dest" <<'PY'
import sys
from pathlib import Path

src = Path(sys.argv[1])
dest = Path(sys.argv[2])
data = src.read_bytes()

encodings = [
    "utf-8",
    "utf-8-sig",
    "utf-16",
    "utf-16-le",
    "utf-16-be",
    "utf-32",
    "utf-32-le",
    "utf-32-be",
    "cp1252",
    "iso-8859-1",
    "iso-8859-15",
    "cp1250",
    "cp1251",
]

for enc in encodings:
    try:
        text = data.decode(enc)
        with open(dest, "w", encoding="utf-8") as f:
            f.write(text)
        sys.exit(0)
    except Exception:
        pass

sys.exit(1)
PY
}

convert_to_utf8() {
    local src="$1"
    local dest="$2"

    if command -v iconv >/dev/null 2>&1; then
        local bom enc
        bom="$(detect_bom_encoding "$src")"

        if [ -n "$bom" ]; then
            if convert_with_iconv "$src" "$dest" "$bom"; then
                return 0
            fi
        fi

        if convert_with_iconv "$src" "$dest" "UTF-8"; then
            return 0
        fi

        if command -v file >/dev/null 2>&1; then
            enc="$(file --brief --mime-encoding "$src" 2>/dev/null || true)"
            enc="$(printf '%s' "$enc" | tr '[:upper:]' '[:lower:]')"

            if [ "$enc" = "us-ascii" ] || [ "$enc" = "ascii" ]; then
                if convert_with_iconv "$src" "$dest" "ASCII"; then
                    return 0
                fi
            fi

            if [ -n "$enc" ] && [ "$enc" != "binary" ] && [ "$enc" != "unknown-8bit" ]; then
                if convert_with_iconv "$src" "$dest" "$enc"; then
                    return 0
                fi
            fi
        fi

        local try_enc
        for try_enc in UTF-8 UTF-16 UTF-16LE UTF-16BE UTF-32 UTF-32LE UTF-32BE WINDOWS-1252 ISO-8859-1 ISO-8859-15 CP1250 CP1251 ASCII; do
            if convert_with_iconv "$src" "$dest" "$try_enc"; then
                return 0
            fi
        done
    fi

    if convert_with_python_fallback "$src" "$dest"; then
        return 0
    fi

    return 1
}

append_skip() {
    printf '%s: %s\n' "$1" "$2" >> "$SKIP_LIST"
}

mkdir -p "$OUT_DIR"
: > "$REL_LIST"
: > "$MAP_LIST"
: > "$SKIP_LIST"

log "Scanning source tree..."

find "$SRC_DIR" -type f | while IFS= read -r file; do
    rel="${file#$SRC_DIR/}"

    if should_ignore_path "$rel"; then
        continue
    fi

    if ! has_source_like_name "$rel"; then
        continue
    fi

    if is_binary "$file"; then
        append_skip "SKIPPED binary" "$rel"
        continue
    fi

    printf '%s\n' "$rel" >> "$REL_LIST"
done

sort -u "$REL_LIST" -o "$REL_LIST"

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

log "Found $(wc -l < "$REL_LIST" | tr -d ' ') relevant files"
log "Copying and converting files..."

while IFS= read -r rel; do
    [ -n "$rel" ] || continue
    src="$SRC_DIR/$rel"

    flat="$(sanitize_path_to_name "$rel")"
    base="$flat"
    n=1
    while [ -e "$OUT_DIR/$flat.txt" ]; do
        flat="${base}__dup${n}"
        n=$((n + 1))
    done

    dest="$OUT_DIR/$flat.txt"

    if convert_to_utf8 "$src" "$dest"; then
        printf '%s\t%s.txt\n' "$rel" "$flat" >> "$MAP_LIST"
    else
        rm -f "$dest"
        append_skip "SKIPPED unconvertible text" "$rel"
    fi
done < "$REL_LIST"

sort -u "$MAP_LIST" -o "$MAP_LIST"
sort -u "$SKIP_LIST" -o "$SKIP_LIST"

{
    echo "ELM upload directory structure"
    echo "=============================="
    echo
    echo "Original source directory:"
    echo "$SRC_DIR"
    echo
    echo "This folder contains flattened UTF-8 copies of source-like text files."
    echo
    echo "Flattening rule:"
    echo "- original relative path separators '/' become '__'"
    echo "- unsafe filename characters are replaced with '_'"
    echo "- '.txt' is appended to the copied filename"
    echo
    echo "Relevant directory tree"
    echo "-----------------------"
    echo "."
} > "$MANIFEST"

if [ -s "$REL_LIST" ]; then
    python3 - "$REL_LIST" >> "$MANIFEST" <<'PY'
import sys
from pathlib import PurePosixPath

paths = []
with open(sys.argv[1], "r", encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n")
        if p:
            paths.append(PurePosixPath(p))

tree = {}
for path in paths:
    node = tree
    for i, part in enumerate(path.parts):
        if i == len(path.parts) - 1:
            node.setdefault("__files__", []).append(part)
        else:
            node = node.setdefault(part, {})

def emit(node, prefix=""):
    dirs = sorted(k for k in node.keys() if k != "__files__")
    files = sorted(node.get("__files__", []))
    items = [("dir", d) for d in dirs] + [("file", f) for f in files]
    for i, (kind, name) in enumerate(items):
        last = i == len(items) - 1
        branch = "└── " if last else "├── "
        if kind == "dir":
            print(prefix + branch + name + "/")
            emit(node[name], prefix + ("    " if last else "│   "))
        else:
            print(prefix + branch + name)

emit(tree)
PY
else
    echo "(no relevant files found)" >> "$MANIFEST"
fi

{
    echo
    echo "File mapping"
    echo "------------"
    echo
} >> "$MANIFEST"

if [ -s "$MAP_LIST" ]; then
    while IFS="$(printf '\t')" read -r rel flat; do
        printf '%s -> %s\n' "$rel" "$flat"
    done < "$MAP_LIST" >> "$MANIFEST"
else
    echo "(no files copied)" >> "$MANIFEST"
fi

{
    echo
    echo "Skipped files"
    echo "-------------"
    echo
} >> "$MANIFEST"

if [ -s "$SKIP_LIST" ]; then
    cat "$SKIP_LIST" >> "$MANIFEST"
else
    echo "(none)" >> "$MANIFEST"
fi

log "Done."
log "Output directory: $OUT_DIR"
log "Manifest: $MANIFEST"
log "Copied files: $(wc -l < "$MAP_LIST" | tr -d ' ')"
log "Skipped files: $(wc -l < "$SKIP_LIST" | tr -d ' ')"
