#!/bin/bash
set -euo pipefail
export PATH="/opt/local/bin:/opt/local/sbin:${PATH}"

log() {
    echo "$*"
}

docker_run() {
    if [ -n "${MPLAPACK_DOCKER_PLATFORM:-}" ] && [ -n "${MPLAPACK_DOCKER_RUN_RUNTIME:-}" ]; then
        docker run --platform "${MPLAPACK_DOCKER_PLATFORM}" --runtime "${MPLAPACK_DOCKER_RUN_RUNTIME}" "$@"
    elif [ -n "${MPLAPACK_DOCKER_PLATFORM:-}" ]; then
        docker run --platform "${MPLAPACK_DOCKER_PLATFORM}" "$@"
    elif [ -n "${MPLAPACK_DOCKER_RUN_RUNTIME:-}" ]; then
        docker run --runtime "${MPLAPACK_DOCKER_RUN_RUNTIME}" "$@"
    else
        docker run "$@"
    fi
}

safe_workdir() {
    local target="$1"
    if [ -z "${HOME:-}" ] || [ "${HOME}" = "/" ]; then
        log "ERROR: HOME is '${HOME:-<unset>}', refusing to use '${target}'."
        exit 1
    fi
    case "${target}" in
        "${HOME}/"*) ;;
        *)
            log "ERROR: '${target}' is not under HOME '${HOME}', refusing to continue."
            exit 1
            ;;
    esac
}

log_ccache_stats() {
    local label="$1"
    log "=== CCACHE STATS (${label}) ==="
    docker_run --rm \
        -v "${MPLAPACK_CCACHE_DIR}:/ccache:rw" \
        "${MPLAPACK_IMAGE_TAG}" \
        ccache -s || true
    log "=== END CCACHE STATS (${label}) ==="
}

: "${MPLAPACK_REMOTE_WORKDIR:?}"
: "${MPLAPACK_DOCKER_BASE:?}"
: "${MPLAPACK_DOCKERFILE:?}"
: "${MPLAPACK_IMAGE_TAG:?}"
: "${MPLAPACK_TARBALL_NAME:?}"
: "${MPLAPACK_CCACHE_DIR:=${HOME}/.ccache}"
: "${MPLAPACK_CCACHE_MAXSIZE:=80G}"
: "${MPLAPACK_COLIMA_CPUS:=$(sysctl -n hw.ncpu 2>/dev/null || echo 10)}"
: "${MPLAPACK_COLIMA_MEMORY_GB:=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 34359738368) / 1024 / 1024 / 1024 / 2 ))}"
: "${MPLAPACK_COLIMA_DISK_GB:=100}"

WORKDIR="${MPLAPACK_REMOTE_WORKDIR}"
CONTEXT_DIR="${WORKDIR}/context"
DOCKER_ROOT="${CONTEXT_DIR}/release/docker"
INPUT_DIR="${WORKDIR}/input"
BUILD_WORK_DIR="${WORKDIR}/work"
LOCKDIR="${WORKDIR}.lock"

safe_workdir "${WORKDIR}"
mkdir -p "$(dirname "${WORKDIR}")"

if mkdir "${LOCKDIR}" 2>/dev/null; then
    printf "%s\n" "$$" > "${LOCKDIR}/pid"
else
    old_pid=""
    [ -f "${LOCKDIR}/pid" ] && old_pid="$(cat "${LOCKDIR}/pid" 2>/dev/null || true)"
    if [ -n "${old_pid}" ] && [ "${old_pid}" != "$$" ] && kill -0 "${old_pid}" 2>/dev/null; then
        log "Another remote tarball build is running (pid: ${old_pid}); stopping it."
        kill "${old_pid}" 2>/dev/null || true
        for _wait_i in $(seq 1 60); do
            kill -0 "${old_pid}" 2>/dev/null || break
            sleep 1
        done
        if kill -0 "${old_pid}" 2>/dev/null; then
            log "Previous remote tarball build pid ${old_pid} did not stop after 60s; killing it."
            kill -KILL "${old_pid}" 2>/dev/null || true
            sleep 1
        fi
    else
        log "WARNING: Stale lock detected; removing it: ${LOCKDIR}"
    fi
    rm -rf "${LOCKDIR}"
    mkdir "${LOCKDIR}"
    printf "%s\n" "$$" > "${LOCKDIR}/pid"
fi

cleanup() {
    local ids
    if command -v docker >/dev/null 2>&1; then
        ids="$(docker ps -q --filter "ancestor=${MPLAPACK_IMAGE_TAG}" 2>/dev/null || true)"
        if [ -n "${ids}" ]; then
            docker stop ${ids} >/dev/null 2>&1 || true
        fi
    fi
    rm -rf "${LOCKDIR}"
}
trap cleanup EXIT INT TERM HUP

if command -v colima >/dev/null 2>&1; then
    colima_status=""
    colima_cpus=""
    colima_memory=""
    colima_disk=""
    if colima list 2>/dev/null | awk '$1 == "default" { found=1 } END { exit found ? 0 : 1 }'; then
        colima_status="$(colima list | awk '$1 == "default" { print $2; exit }')"
        colima_cpus="$(colima list | awk '$1 == "default" { print $4; exit }')"
        colima_memory="$(colima list | awk '$1 == "default" { print $5; exit }')"
        colima_disk="$(colima list | awk '$1 == "default" { print $6; exit }')"
    fi

    desired_memory="${MPLAPACK_COLIMA_MEMORY_GB}GiB"
    desired_disk="${MPLAPACK_COLIMA_DISK_GB}GiB"
    if [ "${colima_status}" != "Running" ] || \
        [ "${colima_cpus}" != "${MPLAPACK_COLIMA_CPUS}" ] || \
        [ "${colima_memory}" != "${desired_memory}" ] || \
        [ "${colima_disk}" != "${desired_disk}" ]; then
        log "Configuring Colima: cpu=${MPLAPACK_COLIMA_CPUS}, memory=${MPLAPACK_COLIMA_MEMORY_GB}GiB, disk=${MPLAPACK_COLIMA_DISK_GB}GiB"
        if [ "${colima_status}" = "Running" ]; then
            colima stop
        fi
        colima start --cpu "${MPLAPACK_COLIMA_CPUS}" --memory "${MPLAPACK_COLIMA_MEMORY_GB}" --disk "${MPLAPACK_COLIMA_DISK_GB}"
    fi
fi

if ! command -v docker >/dev/null 2>&1; then
    log "ERROR: docker command not found. PATH=${PATH}"
    exit 1
fi

docker info >/dev/null
mkdir -p "${MPLAPACK_CCACHE_DIR}" "${BUILD_WORK_DIR}"

log "WORKDIR: ${WORKDIR}"
log "MPLAPACK_DOCKER_BASE: ${MPLAPACK_DOCKER_BASE}"
log "MPLAPACK_DOCKERFILE: ${MPLAPACK_DOCKERFILE}"
log "MPLAPACK_IMAGE_TAG: ${MPLAPACK_IMAGE_TAG}"
log "MPLAPACK_CCACHE_DIR: ${MPLAPACK_CCACHE_DIR}"
log "MPLAPACK_CCACHE_MAXSIZE: ${MPLAPACK_CCACHE_MAXSIZE}"
log "MPLAPACK_TARBALL: ${INPUT_DIR}/${MPLAPACK_TARBALL_NAME}"

if [ ! -f "${DOCKER_ROOT}/${MPLAPACK_DOCKERFILE}" ]; then
    log "ERROR: Dockerfile not found in context: ${DOCKER_ROOT}/${MPLAPACK_DOCKERFILE}"
    exit 1
fi
if [ ! -f "${INPUT_DIR}/${MPLAPACK_TARBALL_NAME}" ]; then
    log "ERROR: tarball not found in input: ${INPUT_DIR}/${MPLAPACK_TARBALL_NAME}"
    exit 1
fi

docker build \
    --build-arg BASE="${MPLAPACK_DOCKER_BASE}" \
    --label org.mplapack.project=mplapack \
    --label org.mplapack.purpose=mplapack-qa \
    --label org.mplapack.source_type=remote-tarball-docker \
    -t "${MPLAPACK_IMAGE_TAG}" \
    -f "${DOCKER_ROOT}/${MPLAPACK_DOCKERFILE}" \
    "${DOCKER_ROOT}"

docker_run --rm \
    -v "${MPLAPACK_CCACHE_DIR}:/ccache:rw" \
    "${MPLAPACK_IMAGE_TAG}" \
    ccache -M "${MPLAPACK_CCACHE_MAXSIZE}"

log_ccache_stats START

docker_run --rm \
    -e CCACHE_DIR=/ccache \
    -e CCACHE_MAXSIZE="${MPLAPACK_CCACHE_MAXSIZE}" \
    -v "${MPLAPACK_CCACHE_DIR}:/ccache:rw" \
    -v "${INPUT_DIR}:/input:ro" \
    -v "${BUILD_WORK_DIR}:/work:rw" \
    "${MPLAPACK_IMAGE_TAG}"

log_ccache_stats END
log "=== REMOTE TARBALL BUILD COMPLETED SUCCESSFULLY ==="