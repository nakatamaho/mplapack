#!/bin/sh
set -u

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
backends="mpfr gmp binary128 binary80 double dd qd"
passed=0
failed=0
skipped=0
seen=0

run_example() {
    exe=$1
    rel=$2
    case $rel in
        *readfromfile*)
            printf 'SKIP %s
' "$rel"
            skipped=$((skipped + 1))
            return
            ;;
    esac
    seen=$((seen + 1))
    if "$exe" >/dev/null 2>&1; then
        printf 'PASS %s
' "$rel"
        passed=$((passed + 1))
    else
        printf 'FAIL %s
' "$rel"
        failed=$((failed + 1))
    fi
}

run_category_tree() {
    base=$1
    for category in "$base"/[0-9][0-9]_* "$base"/90_PrecisionComparison; do
        [ -d "$category" ] || continue
        for backend in $backends; do
            for exe in "$category"/*_"$backend"; do
                [ -f "$exe" ] || continue
                [ -x "$exe" ] || continue
                run_example "$exe" "${exe#"$base"/}"
            done
        done
    done
}

run_cmake_tree() {
    base=$1
    for backend in $backends; do
        [ -d "$base/$backend" ] || continue
        for exe in "$base/$backend"/example_mplapack_*_"$backend"; do
            [ -f "$exe" ] || continue
            [ -x "$exe" ] || continue
            run_example "$exe" "${exe#"$base"/}"
        done
    done
}

if [ -d "examples/mplapack/double" ] || [ -d "examples/mplapack/mpfr" ]; then
    run_cmake_tree "$(pwd)/examples/mplapack"
elif [ -d "double" ] || [ -d "mpfr" ]; then
    run_cmake_tree "$(pwd)"
elif [ -d "00_GeneralLinearEquations" ]; then
    run_category_tree "$(pwd)"
elif [ -d "examples/mplapack/00_GeneralLinearEquations" ]; then
    run_category_tree "$(pwd)/examples/mplapack"
else
    run_category_tree "$script_dir"
fi

printf 'summary: passed=%d failed=%d skipped=%d total=%d
' "$passed" "$failed" "$skipped" "$seen"
[ "$failed" -eq 0 ]
