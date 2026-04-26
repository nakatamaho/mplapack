fable_clang_format_style='{
    BasedOnStyle: llvm,
    IndentWidth: 4,
    ColumnLimit: 10000,
    SortIncludes: false,
    AlignEscapedNewlines: LeftWithLastLine,
    SpaceBeforeRangeBasedForLoopColon: false,
    PointerAlignment: Right,
    NamespaceIndentation: Inner,
    AlwaysBreakTemplateDeclarations: No,
    BreakBeforeConceptDeclarations: Never,
}'

fable_clang_format_header_style='{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }'

fable_clang_format_version_major() {
    "$1" --version 2>/dev/null | sed -n 's/.*version \([0-9][0-9]*\).*/\1/p'
}

fable_check_clang_format_version() {
    fable_clang_format_version="$(fable_clang_format_version_major "$1")"
    if [ -n "$fable_clang_format_version" ] && [ "$fable_clang_format_version" -ge 19 ]; then
        return 0
    fi

    if [ "${2:-}" != "quiet" ]; then
        if [ -n "$fable_clang_format_version" ]; then
            echo "Error: $1 version $fable_clang_format_version is too old. clang-format 19 or newer is required." >&2
        else
            echo "Error: could not determine $1 version. clang-format 19 or newer is required." >&2
        fi
    fi
    return 1
}

fable_find_clang_format() {
    if [ -n "${CLANG_FORMAT:-}" ]; then
        if command -v "$CLANG_FORMAT" >/dev/null 2>&1; then
            fable_check_clang_format_version "$CLANG_FORMAT" && return 0
            exit 1
        fi
        echo "Error: CLANG_FORMAT='$CLANG_FORMAT' is not executable." >&2
        exit 1
    fi

    if command -v clang-format >/dev/null 2>&1; then
        CLANG_FORMAT=clang-format
        fable_check_clang_format_version "$CLANG_FORMAT" quiet && return 0
    fi

    if command -v clang-format-19 >/dev/null 2>&1; then
        CLANG_FORMAT=clang-format-19
        fable_check_clang_format_version "$CLANG_FORMAT" quiet && return 0
    fi

    echo "Error: clang-format 19 or newer is not installed." >&2
    exit 1
}

fable_clang_format_inplace() {
    fable_find_clang_format
    "$CLANG_FORMAT" -i -style "$fable_clang_format_style" "$@"
}

fable_clang_format_stdout() {
    fable_find_clang_format
    "$CLANG_FORMAT" -style "$fable_clang_format_header_style" "$@"
}
