#!/usr/bin/env bash
set -u

section() {
    printf '\n===== %s =====\n' "$1"
}

line() {
    printf '%s\n' "$*"
}

try_enable_module_cmd() {
    if type module >/dev/null 2>&1; then
        return 0
    fi

    for f in \
        /etc/profile.d/modules.sh \
        /usr/share/Modules/init/bash \
        /etc/profile.d/lmod.sh \
        /usr/share/lmod/lmod/init/bash
    do
        if [ -r "$f" ]; then
            # shellcheck disable=SC1090
            . "$f" >/dev/null 2>&1 || true
            if type module >/dev/null 2>&1; then
                return 0
            fi
        fi
    done

    return 1
}

find_compiler() {
    if [ -n "${CC:-}" ] && command -v "${CC}" >/dev/null 2>&1; then
        printf '%s' "${CC}"
        return
    fi

    for c in cc gcc clang icc icx; do
        if command -v "$c" >/dev/null 2>&1; then
            printf '%s' "$c"
            return
        fi
    done

    printf ''
}

ROOTS=()

add_root() {
    local d="$1"

    [ -n "$d" ] || return 0
    [ -d "$d" ] || return 0

    for r in "${ROOTS[@]:-}"; do
        if [ "$r" = "$d" ]; then
            return 0
        fi
    done

    ROOTS+=("$d")
}

add_path_list() {
    local list="${1:-}"
    local old_ifs="$IFS"

    [ -n "$list" ] || return 0

    IFS=':'
    for p in $list; do
        add_root "$p"
    done
    IFS="$old_ifs"
}

run_cmd_test() {
    local name="$1"
    shift

    local err_file="$WORKDIR/${name}.err"
    local out_file="$WORKDIR/${name}.out"

    line "-- ${name}"
    printf 'command:'
    printf ' %q' "$@"
    printf '\n'

    if "$@" >"$out_file" 2>"$err_file"; then
        line "result: OK"
    else
        line "result: FAIL"
        if [ -s "$err_file" ]; then
            line "stderr (first 20 lines):"
            sed -n '1,20p' "$err_file"
        fi
    fi
    line
}

PC_MODULES=(fftw3 fftw3f fftw3_threads fftw3f_threads)

section "Basic environment"
line "date            : $(date)"
line "hostname        : $(hostname 2>/dev/null || true)"
line "pwd             : $(pwd)"
line "shell           : ${SHELL:-<unset>}"
line "CONDA_PREFIX    : ${CONDA_PREFIX:-<unset>}"
line "EBROOTFFTW      : ${EBROOTFFTW:-<unset>}"
line "FFTW_ROOT       : ${FFTW_ROOT:-<unset>}"
line "FFTW_DIR        : ${FFTW_DIR:-<unset>}"
line "PKG_CONFIG_PATH : ${PKG_CONFIG_PATH:-<unset>}"
line "CPATH           : ${CPATH:-<unset>}"
line "LIBRARY_PATH    : ${LIBRARY_PATH:-<unset>}"
line "LD_LIBRARY_PATH : ${LD_LIBRARY_PATH:-<unset>}"
line "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH:-<unset>}"

section "Compiler"
CCBIN="$(find_compiler)"
if [ -n "$CCBIN" ]; then
    line "compiler: $CCBIN"
    "$CCBIN" --version 2>/dev/null | sed -n '1,3p'
else
    line "compiler: <not found>"
fi

section "Module system"
if try_enable_module_cmd; then
    line "module command: available"
    line
    line "module list:"
    module list 2>&1 || true
else
    line "module command: not available inside this script"
fi

section "pkg-config"
if command -v pkg-config >/dev/null 2>&1; then
    line "pkg-config binary: $(command -v pkg-config)"
    line "pkg-config version: $(pkg-config --version 2>/dev/null || true)"
else
    line "pkg-config not found"
fi

# Seed search roots from environment
add_root "${CONDA_PREFIX:-}"
add_root "${EBROOTFFTW:-}"
add_root "${FFTW_ROOT:-}"
add_root "${FFTW_DIR:-}"
add_path_list "${PKG_CONFIG_PATH:-}"
add_path_list "${CPATH:-}"
add_path_list "${LIBRARY_PATH:-}"
add_path_list "${LD_LIBRARY_PATH:-}"
add_path_list "${CMAKE_PREFIX_PATH:-}"
add_root /usr
add_root /usr/local
add_root /opt

if command -v pkg-config >/dev/null 2>&1; then
    for m in "${PC_MODULES[@]}"; do
        section "pkg-config module: ${m}"

        if pkg-config --exists "$m"; then
            line "status   : FOUND"
            line "version  : $(pkg-config --modversion "$m" 2>/dev/null || true)"
            line "cflags   : $(pkg-config --cflags "$m" 2>/dev/null || true)"
            line "libs     : $(pkg-config --libs "$m" 2>/dev/null || true)"
            line "prefix   : $(pkg-config --variable=prefix "$m" 2>/dev/null || true)"
            line "includedir: $(pkg-config --variable=includedir "$m" 2>/dev/null || true)"
            line "libdir   : $(pkg-config --variable=libdir "$m" 2>/dev/null || true)"
            line "pcfiledir: $(pkg-config --variable=pcfiledir "$m" 2>/dev/null || true)"

            add_root "$(pkg-config --variable=prefix "$m" 2>/dev/null || true)"
            add_root "$(pkg-config --variable=includedir "$m" 2>/dev/null || true)"
            add_root "$(pkg-config --variable=libdir "$m" 2>/dev/null || true)"
            add_root "$(pkg-config --variable=pcfiledir "$m" 2>/dev/null || true)"
        else
            line "status   : NOT FOUND"
        fi
    done
fi

section "Candidate search roots"
if [ "${#ROOTS[@]}" -eq 0 ]; then
    line "(none)"
else
    for r in "${ROOTS[@]}"; do
        line "$r"
    done
fi

section "Discovered FFTW files"
FOUND_ANY=0
for r in "${ROOTS[@]}"; do
    [ -d "$r" ] || continue
    while IFS= read -r f; do
        [ -n "$f" ] || continue
        FOUND_ANY=1
        line "$f"
    done < <(
        find "$r" -maxdepth 4 \
            \( -name 'fftw3.h' \\
            -o -name 'libfftw3*.so' \\
            -o -name 'libfftw3*.so.*' \\
            -o -name 'libfftw3*.a' \\
            -o -name 'fftw3*.pc' \\) \\
            -print 2>/dev/null | sort -u
    )
done

if [ "$FOUND_ANY" -eq 0 ]; then
    line "(no FFTW headers/libraries/pkg-config files found under candidate roots)"
fi

section "ldconfig cache"
if command -v ldconfig >/dev/null 2>&1; then
    ldconfig -p 2>/dev/null | grep -i fftw || true
else
    line "ldconfig not available"
fi

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

cat > "$WORKDIR/header_only.c" <<'EOF'
#include <fftw3.h>
int main(void) { return 0; }
EOF

cat > "$WORKDIR/double_core.c" <<'EOF'
#include <fftw3.h>
int main(void)
{
    fftw_complex *x = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 4);
    if (!x) return 2;
    fftw_plan p = fftw_plan_dft_1d(4, x, x, FFTW_FORWARD, FFTW_ESTIMATE);
    if (!p)
    {
        fftw_free(x);
        return 3;
    }
    fftw_destroy_plan(p);
    fftw_free(x);
    return 0;
}
EOF

cat > "$WORKDIR/double_threads.c" <<'EOF'
#include <fftw3.h>
int main(void)
{
    if (!fftw_init_threads()) return 2;
    fftw_plan_with_nthreads(2);
    fftw_complex *x = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 4);
    if (!x) return 3;
    fftw_plan p = fftw_plan_dft_1d(4, x, x, FFTW_FORWARD, FFTW_ESTIMATE);
    if (!p)
    {
        fftw_free(x);
        fftw_cleanup_threads();
        return 4;
    }
    fftw_destroy_plan(p);
    fftw_free(x);
    fftw_cleanup_threads();
    return 0;
}
EOF

cat > "$WORKDIR/float_core.c" <<'EOF'
#include <fftw3.h>
int main(void)
{
    fftwf_complex *x = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * 4);
    if (!x) return 2;
    fftwf_plan p = fftwf_plan_dft_1d(4, x, x, FFTW_FORWARD, FFTW_ESTIMATE);
    if (!p)
    {
        fftwf_free(x);
        return 3;
    }
    fftwf_destroy_plan(p);
    fftwf_free(x);
    return 0;
}
EOF

cat > "$WORKDIR/float_threads.c" <<'EOF'
#include <fftw3.h>
int main(void)
{
    if (!fftwf_init_threads()) return 2;
    fftwf_plan_with_nthreads(2);
    fftwf_complex *x = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * 4);
    if (!x) return 3;
    fftwf_plan p = fftwf_plan_dft_1d(4, x, x, FFTW_FORWARD, FFTW_ESTIMATE);
    if (!p)
    {
        fftwf_free(x);
        fftwf_cleanup_threads();
        return 4;
    }
    fftwf_destroy_plan(p);
    fftwf_free(x);
    fftwf_cleanup_threads();
    return 0;
}
EOF

section "Compiler tests"

if [ -z "$CCBIN" ]; then
    line "No compiler found, skipping compile/link tests."
    exit 0
fi

run_cmd_test "header_only_no_flags" \
    "$CCBIN" -c "$WORKDIR/header_only.c" -o "$WORKDIR/header_only.o"

run_cmd_test "direct_double_core" \
    "$CCBIN" "$WORKDIR/double_core.c" -o "$WORKDIR/direct_double_core.bin" -lfftw3 -lm

run_cmd_test "direct_double_threads" \
    "$CCBIN" "$WORKDIR/double_threads.c" -o "$WORKDIR/direct_double_threads.bin" -lfftw3_threads -lfftw3 -lm

run_cmd_test "direct_float_core" \
    "$CCBIN" "$WORKDIR/float_core.c" -o "$WORKDIR/direct_float_core.bin" -lfftw3f -lm

run_cmd_test "direct_float_threads" \
    "$CCBIN" "$WORKDIR/float_threads.c" -o "$WORKDIR/direct_float_threads.bin" -lfftw3f_threads -lfftw3f -lm

if command -v pkg-config >/dev/null 2>&1; then
    for m in "${PC_MODULES[@]}"; do
        if pkg-config --exists "$m"; then
            flags=()
            # shellcheck disable=SC2206
            flags=( $(pkg-config --cflags --libs "$m") )

            case "$m" in
                fftw3)
                    run_cmd_test "pkgconfig_${m}" \
                        "$CCBIN" "$WORKDIR/double_core.c" -o "$WORKDIR/pkg_${m}.bin" "${flags[@]}"
                    ;;
                fftw3_threads)
                    run_cmd_test "pkgconfig_${m}" \
                        "$CCBIN" "$WORKDIR/double_threads.c" -o "$WORKDIR/pkg_${m}.bin" "${flags[@]}"
                    ;;
                fftw3f)
                    run_cmd_test "pkgconfig_${m}" \
                        "$CCBIN" "$WORKDIR/float_core.c" -o "$WORKDIR/pkg_${m}.bin" "${flags[@]}"
                    ;;
                fftw3f_threads)
                    run_cmd_test "pkgconfig_${m}" \
                        "$CCBIN" "$WORKDIR/float_threads.c" -o "$WORKDIR/pkg_${m}.bin" "${flags[@]}"
                    ;;
            esac
        fi
    done
fi

section "Done"
line "Please paste the full output of this script."
