#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}"
CONF_FILE="${CONF_FILE:-$SCRIPT_DIR/build-matrix.conf}"
LOGDIR="${LOGDIR:-$SCRIPT_DIR/logs/$(date +%Y%m%d_%H%M%S)}"
SUCCESS_DIR="${SUCCESS_DIR:-$SCRIPT_DIR/success}"
REMOTE_SSH="${REMOTE_LINUX_SSH:-ssh}"
MPLAPACK_REF="${MPLAPACK_REF:-$(git -C "$PROJECT_ROOT" rev-parse HEAD)}"

name="${1:?Usage: run-remote-linux-docker.sh <matrix-name>}"

mkdir -p "$LOGDIR" "$SUCCESS_DIR"

row="$(awk -F'|' -v target="$name" '$1 == target { print; exit }' "$CONF_FILE")"
if [[ -z "$row" ]]; then
    echo "ERROR: $name not found in $CONF_FILE" >&2
    exit 1
fi

IFS='|' read -r matrix_name host target_dir script_rel remote_cmd source_type matrix_ubuntu_version matrix_docker_base <<< "$row"
remote_cmd="${remote_cmd:-bash}"
matrix_ubuntu_version="${matrix_ubuntu_version:-}"
matrix_docker_base="${matrix_docker_base:-}"
remote_ubuntu_version="${MPLAPACK_UBUNTU_VERSION:-$matrix_ubuntu_version}"
remote_docker_base="${MPLAPACK_DOCKER_BASE:-$matrix_docker_base}"
if [[ -z "$remote_docker_base" && -n "$remote_ubuntu_version" ]]; then
    remote_docker_base="ubuntu:$remote_ubuntu_version"
fi

if [[ "$source_type" != "remote-linux-docker" ]]; then
    echo "ERROR: $name is source_type=$source_type, expected remote-linux-docker" >&2
    exit 1
fi
if [[ -z "$host" || -z "$target_dir" || -z "$script_rel" ]]; then
    echo "ERROR: invalid remote-linux-docker row: $row" >&2
    exit 1
fi

script_path="$PROJECT_ROOT/$script_rel"
if [[ ! -f "$script_path" ]]; then
    echo "ERROR: build script not found: $script_path" >&2
    exit 1
fi

arch="${matrix_name%%-*}"
docker_label="${remote_docker_base:-ubuntu:${remote_ubuntu_version:-default}}"
docker_label="$(printf '%s' "$docker_label" | sed 's/[^A-Za-z0-9._+-]/_/g')"
logfile="$LOGDIR/${matrix_name}_${docker_label}_${arch}_remote-linux-docker.log"
resultfile="$LOGDIR/results_${matrix_name}.csv"
stamp="$SUCCESS_DIR/${matrix_name}_${docker_label}_${arch}_remote-linux-docker.ok"

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

    if ! "$REMOTE_SSH" "$host" "test -d '$remote_root'" >/dev/null 2>&1; then
        echo "No remote distcheck results found: $host:$remote_root" >&2
        return 0
    fi

    rm -rf "$local_stage"
    mkdir -p "$local_stage"
    if command -v rsync >/dev/null 2>&1 && rsync -a "$host:$remote_root/" "$local_stage/"; then
        :
    else
        "$REMOTE_SSH" "$host" "cd '$remote_root' && tar -cf - ." | tar -C "$local_stage" -xf -
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

if [[ -L "$stamp" && ! -e "$stamp" ]]; then
    rm -f "$stamp"
fi
if [[ -L "$stamp" && -e "$stamp" ]]; then
    echo "Already passed; skipping $matrix_name" >&2
    echo "Success log: $(readlink "$stamp")" >&2
    echo "$matrix_name,$arch,${remote_docker_base:-$host},test,SKIPPED,0,$source_type" | tee -a "$resultfile"
    exit 0
fi

start="$(date +%s)"
echo "Running $script_rel on $host:$target_dir" >&2
echo "MPLAPACK_REF: $MPLAPACK_REF" >&2
echo "MPLAPACK_UBUNTU_VERSION: ${remote_ubuntu_version:-<default>}" >&2
echo "MPLAPACK_DOCKER_BASE: ${remote_docker_base:-<default>}" >&2
echo "Log: $logfile" >&2

context_tar="$LOGDIR/${matrix_name}_context.tar.gz"
remote_context_tar="${target_dir}.context.tar.gz"
tar -C "$PROJECT_ROOT" -czf "$context_tar" release/docker

set +e
{
    "$REMOTE_SSH" "$host" "mkdir -p '$(dirname "$target_dir")' && cat > '$remote_context_tar'" < "$context_tar" && \
    "$REMOTE_SSH" "$host" "MPLAPACK_REMOTE_WORKDIR='$target_dir' MPLAPACK_CONTEXT_TARBALL='$remote_context_tar' MPLAPACK_REF='$MPLAPACK_REF' MPLAPACK_UBUNTU_VERSION='$remote_ubuntu_version' MPLAPACK_DOCKER_BASE='$remote_docker_base' $remote_cmd -s" < "$script_path"
} > "$logfile" 2>&1
rc=$?
set -e

elapsed=$(($(date +%s) - start))
if [[ "$rc" -eq 0 ]]; then
    collect_remote_test_results
    ln -sfn "$logfile" "$stamp"
    echo "$matrix_name,$arch,${remote_docker_base:-$host},test,OK,$elapsed,$source_type" | tee -a "$resultfile"
else
    echo "$matrix_name,$arch,${remote_docker_base:-$host},test,FAILED,$elapsed,$source_type" | tee -a "$resultfile"
    exit "$rc"
fi
