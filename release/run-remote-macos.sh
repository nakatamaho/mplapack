#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}"
CONF_FILE="${CONF_FILE:-$SCRIPT_DIR/build-matrix.conf}"
LOGDIR="${LOGDIR:-$SCRIPT_DIR/logs/$(date +%Y%m%d_%H%M%S)}"
SUCCESS_DIR="${SUCCESS_DIR:-$SCRIPT_DIR/success}"
MACOS_REMOTE_SSH="${MACOS_REMOTE_SSH:-ssh}"

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

echo "name,arch,base,stage,result,elapsed,source_type" > "$resultfile"
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

start="$(date +%s)"
echo "Running $script_rel on $host:$target_dir" >&2
echo "Log: $logfile" >&2

set +e
"$MACOS_REMOTE_SSH" "$host" "MPLAPACK_REMOTE_WORKDIR=$target_dir $remote_cmd -s" < "$script_path" > "$logfile" 2>&1
rc=$?
set -e

elapsed=$(($(date +%s) - start))
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
