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

arch="${matrix_name%%-*}"
logfile="$LOGDIR/${matrix_name}_${arch}_remote-macos.log"
resultfile="$LOGDIR/results_${matrix_name}.csv"
stamp="$SUCCESS_DIR/${matrix_name}_${arch}_remote-macos.ok"

echo "name,arch,base,stage,result,elapsed,source_type" > "$resultfile"

if [[ -L "$stamp" && ! -e "$stamp" ]]; then
    rm -f "$stamp"
fi
if [[ -L "$stamp" && -e "$stamp" ]]; then
    echo "Already passed; skipping $matrix_name" >&2
    echo "Success log: $(readlink "$stamp")" >&2
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
    ln -sfn "$logfile" "$stamp"
    echo "$matrix_name,$arch,$host,test,OK,$elapsed,$source_type" | tee -a "$resultfile"
else
    echo "$matrix_name,$arch,$host,test,FAILED,$elapsed,$source_type" | tee -a "$resultfile"
    exit "$rc"
fi
