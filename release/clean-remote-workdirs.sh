#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONF_FILE="${CONF_FILE:-$SCRIPT_DIR/build-matrix.conf}"
MACOS_REMOTE_SSH="${MACOS_REMOTE_SSH:-ssh}"
REMOTE_LINUX_SSH="${REMOTE_LINUX_SSH:-ssh}"

usage() {
    echo "Usage: clean-remote-workdirs.sh" >&2
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [[ ! -f "$CONF_FILE" ]]; then
    echo "ERROR: build matrix not found: $CONF_FILE" >&2
    exit 1
fi

sanitize_label() {
    printf '%s' "$1" | sed 's#^.*/##; s/[^A-Za-z0-9]/_/g'
}

linux_context_name() {
    local matrix_name="$1"
    local docker_base="$2"
    local tier arch docker_name docker_version os_label name_middle flavor

    tier="${matrix_name%%-*}"
    arch="${matrix_name##*-}"
    docker_name="${docker_base%%:*}"
    docker_version="${docker_base#*:}"
    if [[ "$docker_name" == "$docker_base" ]]; then
        docker_version=""
    fi
    os_label="$(sanitize_label "$docker_name")$(sanitize_label "$docker_version")"
    if [[ -z "$os_label" ]]; then
        os_label="linux_unknown"
    fi
    if [[ "$matrix_name" == *mingw* ]]; then
        os_label="windows"
    fi

    name_middle="${matrix_name#${tier}-}"
    name_middle="${name_middle%-${arch}}"
    flavor=""
    if [[ "$os_label" == "windows" ]]; then
        flavor="$name_middle"
    elif [[ "$name_middle" == *-* ]]; then
        flavor="${name_middle#*-}"
    fi
    printf '%s-%s%s-%s' "$tier" "$os_label" "${flavor:+-$flavor}" "$arch"
}

is_safe_remote_workdir() {
    case "$1" in
        /Users/*/tmp/*|/home/*/tmp/*) return 0 ;;
        *) return 1 ;;
    esac
}

clean_one() {
    local name="$1"
    local host="$2"
    local target_dir="$3"
    local source_type="$4"
    local docker_base="${5:-}"
    local ssh_cmd context_tar context_name

    if [[ -z "$host" || -z "$target_dir" ]]; then
        echo "SKIP $name: missing host or target_dir" >&2
        return 0
    fi
    if ! is_safe_remote_workdir "$target_dir"; then
        echo "SKIP $name: refusing unsafe remote workdir: $host:$target_dir" >&2
        return 0
    fi

    ssh_cmd="$MACOS_REMOTE_SSH"
    context_tar="${target_dir}.context.tar.gz"
    if [[ "$source_type" == "remote-linux-docker" ]]; then
        ssh_cmd="$REMOTE_LINUX_SSH"
        context_name="$(linux_context_name "$name" "$docker_base")"
        context_tar="$(dirname "$target_dir")/${context_name}.context.tar.gz"
    elif [[ "$source_type" == "remote-tarball-docker" ]]; then
        ssh_cmd="$REMOTE_LINUX_SSH"
    fi

    echo "Cleaning $name on $host:$target_dir" >&2
    remote_workdir_q="$(printf '%q' "$target_dir")"
    remote_context_tar_q="$(printf '%q' "$context_tar")"
    remote_source_type_q="$(printf '%q' "$source_type")"
    remote_docker_base_q="$(printf '%q' "$docker_base")"
    "$ssh_cmd" "$host" "MPLAPACK_CLEAN_WORKDIR=$remote_workdir_q MPLAPACK_CLEAN_CONTEXT_TAR=$remote_context_tar_q MPLAPACK_CLEAN_SOURCE_TYPE=$remote_source_type_q MPLAPACK_CLEAN_DOCKER_BASE=$remote_docker_base_q bash -s" <<'REMOTE_CLEAN'
set -eu
workdir="${MPLAPACK_CLEAN_WORKDIR:?}"
context_tar="${MPLAPACK_CLEAN_CONTEXT_TAR:?}"
source_type="${MPLAPACK_CLEAN_SOURCE_TYPE:?}"
docker_base="${MPLAPACK_CLEAN_DOCKER_BASE:-}"
lockdir="${workdir}.lock"
results_dir="${workdir}.distcheck-results"
legacy_context_tar="${workdir}.context.tar.gz"

case "$workdir" in
    /Users/*/tmp/*|/home/*/tmp/*) ;;
    *)
        echo "ERROR: refusing unsafe remote workdir: $workdir" >&2
        exit 1
        ;;
esac

if [ -f "${lockdir}/pid" ]; then
    old_pid="$(cat "${lockdir}/pid" 2>/dev/null || true)"
    if [ -n "$old_pid" ] && kill -0 "$old_pid" 2>/dev/null; then
        echo "Stopping remote build pid $old_pid for $workdir" >&2
        kill "$old_pid" 2>/dev/null || true
        i=0
        while [ "$i" -lt 60 ] && kill -0 "$old_pid" 2>/dev/null; do
            sleep 1
            i=$((i + 1))
        done
        if kill -0 "$old_pid" 2>/dev/null; then
            echo "Killing remote build pid $old_pid for $workdir" >&2
            kill -KILL "$old_pid" 2>/dev/null || true
            sleep 1
        fi
    fi
fi

stop_remote_docker_jobs() {
    if ! command -v docker >/dev/null 2>&1; then
        return 0
    fi

    ids="$(docker ps -q --filter label=org.mplapack.project=mplapack --filter label=org.mplapack.purpose=mplapack-qa 2>/dev/null || true)"
    if [ -n "$ids" ]; then
        echo "Stopping labeled MPLAPACK Docker containers: $ids" >&2
        docker stop $ids >/dev/null 2>&1 || true
    fi

    for image in \
        mplapack-tier1-ubuntu-amd64:latest \
        mplapack-tier1-ubuntu-arm64:latest \
        mplapack-tier1-debian-i386:latest \
        mplapack-tier1-mingw64-amd64:latest \
        mplapack-tier1-ubuntu-inteloneapi-amd64:latest; do
        ids="$(docker ps -q --filter ancestor=$image 2>/dev/null || true)"
        if [ -n "$ids" ]; then
            echo "Stopping legacy MPLAPACK Docker containers for $image: $ids" >&2
            docker stop $ids >/dev/null 2>&1 || true
        fi
    done
}

cleanup_with_docker() {
    parent_dir="$(dirname "$workdir")"
    work_base="$(basename "$workdir")"
    context_base="$(basename "$context_tar")"
    context_parent="$(dirname "$context_tar")"
    image="${docker_base:-ubuntu:24.04}"
    docker_platform=""

    if ! command -v docker >/dev/null 2>&1; then
        return 1
    fi
    if [ "$context_parent" != "$parent_dir" ]; then
        return 1
    fi
    case "$(uname -m 2>/dev/null || true)" in
        x86_64|amd64) docker_platform="linux/amd64" ;;
        aarch64|arm64) docker_platform="linux/arm64" ;;
        i386|i686) docker_platform="linux/386" ;;
    esac

    echo "Retrying cleanup via Docker root helper: $image" >&2
    if [ -n "$docker_platform" ]; then
        docker run --rm \
            --platform "$docker_platform" \
            -v "${parent_dir}:/cleanup:rw" \
            "$image" \
            sh -c 'chmod -R u+rwX "/cleanup/$1" "/cleanup/$2" 2>/dev/null || true; rm -rf -- "/cleanup/$1" "/cleanup/$1.lock" "/cleanup/$2" "/cleanup/$3" "/cleanup/$4"' \
            sh "$work_base" "${work_base}.distcheck-results" "${work_base}.context.tar.gz" "$context_base"
    else
        docker run --rm \
            -v "${parent_dir}:/cleanup:rw" \
            "$image" \
            sh -c 'chmod -R u+rwX "/cleanup/$1" "/cleanup/$2" 2>/dev/null || true; rm -rf -- "/cleanup/$1" "/cleanup/$1.lock" "/cleanup/$2" "/cleanup/$3" "/cleanup/$4"' \
            sh "$work_base" "${work_base}.distcheck-results" "${work_base}.context.tar.gz" "$context_base"
    fi
}

stop_remote_docker_jobs

chmod -R u+rwX "$workdir" "$results_dir" 2>/dev/null || true
if ! rm -rf -- "$workdir" "$lockdir" "$results_dir" "$legacy_context_tar" "$context_tar" 2>/dev/null; then
    if [ "$source_type" = "remote-linux-docker" ] || [ "$source_type" = "remote-tarball-docker" ]; then
        echo "Regular cleanup could not remove every file; using Docker root helper." >&2
        cleanup_with_docker
    else
        echo "ERROR: failed to remove remote workdir: $workdir" >&2
        exit 1
    fi
fi
REMOTE_CLEAN
}

while IFS='|' read -r name host target_dir script_rel remote_cmd source_type docker_base rest; do
    case "$name" in
        ''|'#'*) continue ;;
    esac
    case "$source_type" in
        remote-macos|remote-linux-docker|remote-tarball-docker) clean_one "$name" "$host" "$target_dir" "$source_type" "${docker_base:-}" ;;
    esac
done < "$CONF_FILE"
