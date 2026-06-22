#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}"
CONF_FILE="${CONF_FILE:-$SCRIPT_DIR/build-matrix.conf}"
LOGDIR="${LOGDIR:-$SCRIPT_DIR/logs/$(date +%Y%m%d_%H%M%S)}"
SUCCESS_DIR="${SUCCESS_DIR:-$SCRIPT_DIR/success}"
MACOS_REMOTE_SSH="${MACOS_REMOTE_SSH:-ssh}"
MPLAPACK_REF="${MPLAPACK_REF:-$(git -C "$PROJECT_ROOT" rev-parse HEAD)}"
MPLAPACK_SOURCE_MODE="${MPLAPACK_SOURCE_MODE:-dist}"

name="${1:?Usage: run-remote-macos.sh <matrix-name>}"

mkdir -p "$LOGDIR" "$SUCCESS_DIR"

row="$(awk -F'|' -v target="$name" '$1 == target { print; exit }' "$CONF_FILE")"
if [[ -z "$row" ]]; then
    echo "ERROR: $name not found in $CONF_FILE" >&2
    exit 1
fi

IFS='|' read -r matrix_name host target_dir script_rel remote_cmd source_type <<< "$row"
remote_cmd="${remote_cmd:-bash}"

if [[ "$source_type" != "remote-macos" ]]; then
    echo "ERROR: $name is source_type=$source_type, expected remote-macos" >&2
    exit 1
fi
if [[ -z "$host" || -z "$target_dir" || -z "$script_rel" ]]; then
    echo "ERROR: invalid remote-macos row: $row" >&2
    exit 1
fi

script_path="$PROJECT_ROOT/$script_rel"
if [[ ! -f "$script_path" ]]; then
    echo "ERROR: build script not found: $script_path" >&2
    exit 1
fi

tier="${matrix_name%%-*}"
arch="${matrix_name##*-}"
macos_version="$($MACOS_REMOTE_SSH "$host" "sw_vers -productVersion 2>/dev/null | awk -F. '{print \$1}'" 2>/dev/null || true)"
if [[ -n "$macos_version" ]]; then
    os_label="macos${macos_version}"
else
    os_label="macos_unknown"
fi
stamp_prefix="$SUCCESS_DIR/${tier}-${os_label}-${arch}-"
stamp_suffix="-macos.ok"
logfile="$LOGDIR/${tier}-${os_label}-${arch}-macos.log"
resultfile="$LOGDIR/results_${tier}-${os_label}-${arch}.csv"
stamp=""
COLLECTED_RESULTS_STAGE=""

cleanup_stale_success_links() {
    find "$SUCCESS_DIR" -maxdepth 1 -xtype l -name '*.ok' -exec rm -f {} +
}

detect_compiler_label_from_results() {
    local stage="$1"
    local outdir base compiler

    if [[ -d "$stage" ]]; then
        outdir="$(find "$stage" -path '*/results/*' -type d -mindepth 4 -maxdepth 5 | head -n 1 || true)"
        if [[ -n "$outdir" ]]; then
            base="$(basename "$outdir")"
            compiler="$(printf '%s' "$base" | sed -n 's/^.*_\([A-Za-z][A-Za-z0-9+.-]*-[0-9][A-Za-z0-9_.-]*\)$/\1/p')"
            if [[ -n "$compiler" ]]; then
                printf '%s' "$compiler" | sed 's/-//; s/[^A-Za-z0-9_+.-]/_/g'
                return 0
            fi
        fi
    fi
    printf 'compiler_unknown'
}

copy_tree_contents() {
    local src="$1"
    local dst="$2"
    mkdir -p "$dst"
    if command -v rsync >/dev/null 2>&1; then
        rsync -a "$src/" "$dst/"
    else
        cp -a "$src/." "$dst/"
    fi
}

collect_remote_test_results() {
    local remote_root="${target_dir}.distcheck-results"
    local local_stage="$LOGDIR/${matrix_name}_distcheck-results"
    local version_dir version suite src dst copied
    COLLECTED_RESULTS_STAGE="$local_stage"

    if ! "$MACOS_REMOTE_SSH" "$host" "test -d '$remote_root'" >/dev/null 2>&1; then
        echo "No remote distcheck results found: $host:$remote_root" >&2
        return 0
    fi

    rm -rf "$local_stage"
    mkdir -p "$local_stage"
    if command -v rsync >/dev/null 2>&1 && rsync -a "$host:$remote_root/" "$local_stage/"; then
        :
    else
        "$MACOS_REMOTE_SSH" "$host" "cd '$remote_root' && tar -cf - ." | tar -C "$local_stage" -xf -
    fi

    copied=0
    for version_dir in "$local_stage"/*; do
        [[ -d "$version_dir" ]] || continue
        version="$(basename "$version_dir")"
        for suite in lin eig; do
            src="$version_dir/$suite/results"
            [[ -d "$src" ]] || continue
            dst="$PROJECT_ROOT/mplapack/test/$suite/results/$version"
            copy_tree_contents "$src" "$dst"
            copied=$((copied + 1))
            echo "Collected remote $suite results: $host:$remote_root/$version/$suite/results -> $dst" >&2
        done
    done

    if [[ "$copied" -eq 0 ]]; then
        echo "No staged lin/eig results found under $host:$remote_root" >&2
    fi
}

lang_c_date() {
    LANG=C LC_ALL=C date
}

format_elapsed() {
    local seconds="$1"
    printf '%02d:%02d:%02d' "$((seconds / 3600))" "$(((seconds % 3600) / 60))" "$((seconds % 60))"
}

write_log_start() {
    local epoch="$1"
    {
        echo "=== LOG START: ${matrix_name} ==="
        echo "LOG_START_EPOCH=${epoch}"
        printf 'LOG_START_DATE='
        lang_c_date
        echo
    } >> "$logfile"
}

write_log_end() {
    local rc="$1"
    local end_epoch="$2"
    local elapsed="$3"
    {
        echo
        echo "=== LOG END: ${matrix_name} ==="
        echo "LOG_END_EPOCH=${end_epoch}"
        printf 'LOG_END_DATE='
        lang_c_date
        echo "LOG_ELAPSED_SECONDS=${elapsed}"
        echo "LOG_ELAPSED_HMS=$(format_elapsed "$elapsed")"
        echo "=== LOG ELAPSED: ${elapsed}s ($(format_elapsed "$elapsed")) | rc: ${rc} ==="
    } >> "$logfile"
}

echo "name,arch,base,stage,result,elapsed,source_type" > "$resultfile"
if [[ "$MPLAPACK_SOURCE_MODE" == "dist" && -z "${MPLAPACK_SOURCE_TARBALL:-}" ]]; then
    source_info_file="$LOGDIR/${tier}-${arch}_source.env"
    PROJECT_ROOT="$PROJECT_ROOT" LOGDIR="$LOGDIR" MPLAPACK_SOURCE_INFO_FILE="$source_info_file" \
        "$SCRIPT_DIR/make-source-snapshot.sh"
    # shellcheck disable=SC1090
    source "$source_info_file"
elif [[ "$MPLAPACK_SOURCE_MODE" != "dist" && "$MPLAPACK_SOURCE_MODE" != "ref" ]]; then
    echo "ERROR: MPLAPACK_SOURCE_MODE must be 'dist' or 'ref' (got '$MPLAPACK_SOURCE_MODE')" >&2
    exit 1
fi

if [[ "$MPLAPACK_SOURCE_MODE" == "dist" ]]; then
    stamp_suffix="-${MPLAPACK_SOURCE_LABEL:-dist}-macos.ok"
fi

cleanup_stale_success_links

existing_stamp="$(find "$SUCCESS_DIR" -maxdepth 1 -type l -name "$(basename "$stamp_prefix")*${stamp_suffix}" | head -n 1 || true)"
if [[ -n "$existing_stamp" && ! -e "$existing_stamp" ]]; then
    rm -f "$existing_stamp"
    existing_stamp=""
fi
if [[ -n "$existing_stamp" && -e "$existing_stamp" ]]; then
    echo "Already passed; skipping $matrix_name" >&2
    echo "Success log: $(readlink "$existing_stamp")" >&2
    echo "$matrix_name,$arch,$host,test,SKIPPED,0,$source_type" | tee -a "$resultfile"
    exit 0
fi


echo "Running $script_rel on $host:$target_dir" >&2
echo "MPLAPACK_REF: $MPLAPACK_REF" >&2
echo "MPLAPACK_SOURCE_MODE: $MPLAPACK_SOURCE_MODE" >&2
echo "Log: $logfile" >&2

context_name="${tier}-${os_label}-${arch}"
remote_context_tar="$(dirname "$target_dir")/${context_name}.context.tar.gz"
context_tar="$LOGDIR/${context_name}_context.tar.gz"
if [[ "$MPLAPACK_SOURCE_MODE" == "dist" ]]; then
    if [[ -z "${MPLAPACK_SOURCE_TARBALL:-}" || ! -f "$MPLAPACK_SOURCE_TARBALL" ]]; then
        echo "ERROR: MPLAPACK_SOURCE_TARBALL is not set or missing: ${MPLAPACK_SOURCE_TARBALL:-<unset>}" >&2
        exit 1
    fi
    tar_args=(-C "$(dirname "$MPLAPACK_SOURCE_TARBALL")" "$(basename "$MPLAPACK_SOURCE_TARBALL")")
    for source_extra in "${MPLAPACK_SOURCE_METADATA:-}" "${MPLAPACK_SOURCE_PATCH:-}" "${MPLAPACK_SOURCE_STATUS:-}"; do
        if [[ -n "$source_extra" && -f "$source_extra" ]]; then
            tar_args+=(-C "$(dirname "$source_extra")" "$(basename "$source_extra")")
        fi
    done
    tar -czf "$context_tar" "${tar_args[@]}"
fi

start="$(date +%s)"
set +e
write_log_start "$start"
if [[ "$MPLAPACK_SOURCE_MODE" == "dist" ]]; then
    "$MACOS_REMOTE_SSH" "$host" "mkdir -p '$(dirname "$target_dir")' && cat > '$remote_context_tar'" < "$context_tar" && \
        "$MACOS_REMOTE_SSH" "$host" "MPLAPACK_REMOTE_WORKDIR='$target_dir' MPLAPACK_CONTEXT_TARBALL='$remote_context_tar' MPLAPACK_REF='$MPLAPACK_REF' MPLAPACK_SOURCE_MODE='$MPLAPACK_SOURCE_MODE' $remote_cmd -s" < "$script_path" >> "$logfile" 2>&1
    rc=$?
else
    "$MACOS_REMOTE_SSH" "$host" "MPLAPACK_REMOTE_WORKDIR='$target_dir' MPLAPACK_REF='$MPLAPACK_REF' MPLAPACK_SOURCE_MODE='$MPLAPACK_SOURCE_MODE' $remote_cmd -s" < "$script_path" >> "$logfile" 2>&1
    rc=$?
fi
set -e

end="$(date +%s)"
elapsed=$((end - start))
write_log_end "$rc" "$end" "$elapsed"
if [[ "$rc" -eq 0 ]]; then
    collect_remote_test_results
    compiler_label="$(detect_compiler_label_from_results "$COLLECTED_RESULTS_STAGE")"
    if [[ "$compiler_label" == "compiler_unknown" ]]; then
        echo "ERROR: compiler label could not be determined from collected results; refusing to create success stamp." >&2
        echo "$matrix_name,$arch,${remote_docker_base:-${host:-}},test,FAILED,$elapsed,$source_type" | tee -a "$resultfile"
        exit 1
    fi
    final_logfile="$LOGDIR/${tier}-${os_label}-${arch}-${compiler_label}-macos.log"
    if [[ "$final_logfile" != "$logfile" ]]; then
        mv "$logfile" "$final_logfile"
        logfile="$final_logfile"
    fi
    stamp="${stamp_prefix}${compiler_label}${stamp_suffix}"
    ln -sfn "$logfile" "$stamp"
    echo "$matrix_name,$arch,$host,test,OK,$elapsed,$source_type" | tee -a "$resultfile"
else
    echo "$matrix_name,$arch,$host,test,FAILED,$elapsed,$source_type" | tee -a "$resultfile"
    exit "$rc"
fi
