#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}"
CONF_FILE="${CONF_FILE:-$SCRIPT_DIR/build-matrix.conf}"
LOGDIR="${LOGDIR:-$SCRIPT_DIR/logs/$(date +%Y%m%d_%H%M%S)}"
SUCCESS_DIR="${SUCCESS_DIR:-$SCRIPT_DIR/success}"
MACOS_REMOTE_SSH="${MACOS_REMOTE_SSH:-ssh}"
REMOTE_LINUX_SSH="${REMOTE_LINUX_SSH:-ssh}"

if [[ "$#" -eq 0 ]]; then
    echo "Usage: run-remote-targets-by-host.sh <target>..." >&2
    exit 1
fi
if [[ ! -f "$CONF_FILE" ]]; then
    echo "ERROR: build matrix not found: $CONF_FILE" >&2
    exit 1
fi

declare -A host_seen=()
declare -A host_target_files=()
host_order=()
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

matrix_row() {
    local target="$1"
    awk -F'|' -v target="$target" '$1 == target { print; exit }' "$CONF_FILE"
}

for target in "$@"; do
    row="$(matrix_row "$target")"
    if [[ -z "$row" ]]; then
        echo "ERROR: $target not found in $CONF_FILE" >&2
        exit 1
    fi
    IFS='|' read -r matrix_name host target_dir script_rel remote_cmd source_type rest <<< "$row"
    case "$source_type" in
        remote-macos|remote-linux-docker) ;;
        *)
            echo "ERROR: $target is source_type=$source_type, expected remote target" >&2
            exit 1
            ;;
    esac
    if [[ -z "$host" ]]; then
        echo "ERROR: $target has empty host in $CONF_FILE" >&2
        exit 1
    fi
    if [[ -z "${host_seen[$host]:-}" ]]; then
        host_seen[$host]=1
        host_order+=("$host")
        host_target_files[$host]="$tmpdir/host_${#host_order[@]}.targets"
        : > "${host_target_files[$host]}"
    fi
    printf '%s\n' "$target" >> "${host_target_files[$host]}"
done

run_target() {
    local target="$1"
    local row matrix_name host target_dir script_rel remote_cmd source_type rest

    row="$(matrix_row "$target")"
    IFS='|' read -r matrix_name host target_dir script_rel remote_cmd source_type rest <<< "$row"
    case "$source_type" in
        remote-macos)
            CONF_FILE="$CONF_FILE" LOGDIR="$LOGDIR" PROJECT_ROOT="$PROJECT_ROOT" SUCCESS_DIR="$SUCCESS_DIR" \
                MACOS_REMOTE_SSH="$MACOS_REMOTE_SSH" "$SCRIPT_DIR/run-remote-macos.sh" "$target"
            ;;
        remote-linux-docker)
            CONF_FILE="$CONF_FILE" LOGDIR="$LOGDIR" PROJECT_ROOT="$PROJECT_ROOT" SUCCESS_DIR="$SUCCESS_DIR" \
                REMOTE_LINUX_SSH="$REMOTE_LINUX_SSH" "$SCRIPT_DIR/run-remote-linux-docker.sh" "$target"
            ;;
    esac
}

run_host_targets() {
    local host="$1"
    local host_rc target target_file target_rc

    host_rc=0
    target_file="${host_target_files[$host]}"
    echo "=== Running remote targets on $host ===" >&2
    while IFS= read -r target; do
        [[ -n "$target" ]] || continue
        echo "=== Running remote target on $host: $target ===" >&2
        # Keep ssh in child scripts from consuming the remaining target list.
        if run_target "$target" </dev/null; then
            echo "=== Finished remote target on $host: $target ===" >&2
        else
            target_rc="$?"
            echo "ERROR: remote target failed on $host: $target (rc=$target_rc)" >&2
            echo "=== Finished remote target on $host: $target (FAILED) ===" >&2
            host_rc=1
        fi
    done < "$target_file"
    return "$host_rc"
}

pids=()
pid_hosts=()
for host in "${host_order[@]}"; do
    run_host_targets "$host" &
    pids+=("$!")
    pid_hosts+=("$host")
done

rc=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
        echo "ERROR: remote target group failed on ${pid_hosts[$i]}" >&2
        rc=1
    fi
done
exit "$rc"
