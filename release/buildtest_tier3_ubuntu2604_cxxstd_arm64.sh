#!/bin/bash
set -euo pipefail
export PATH="/opt/local/bin:/opt/local/sbin:${PATH}"

log() {
    echo "$*"
}

docker_run() {
    local labels=(
        --label org.mplapack.project=mplapack
        --label org.mplapack.purpose=mplapack-qa
    )
    if [ -n "${MPLAPACK_DOCKER_PLATFORM:-}" ] && [ -n "${MPLAPACK_DOCKER_RUN_RUNTIME:-}" ]; then
        docker run --platform "${MPLAPACK_DOCKER_PLATFORM}" --runtime "${MPLAPACK_DOCKER_RUN_RUNTIME}" "${labels[@]}" "$@"
    elif [ -n "${MPLAPACK_DOCKER_PLATFORM:-}" ]; then
        docker run --platform "${MPLAPACK_DOCKER_PLATFORM}" "${labels[@]}" "$@"
    elif [ -n "${MPLAPACK_DOCKER_RUN_RUNTIME:-}" ]; then
        docker run --runtime "${MPLAPACK_DOCKER_RUN_RUNTIME}" "${labels[@]}" "$@"
    else
        docker run "${labels[@]}" "$@"
    fi
}

log_ccache_stats() {
    local label="$1"

    log "=== CCACHE STATS (${label}) ==="
    docker_run --rm \
        -v "${CCACHE_DIR_HOST}:/ccache:rw" \
        "${MPLAPACK_IMAGE_TAG}" \
        ccache -s || true
    log "=== END CCACHE STATS (${label}) ==="
}

log_ccache_start() {
    log_ccache_stats "START"
    CCACHE_STATS_STARTED=1
}

log_ccache_end_once() {
    if [ "${CCACHE_STATS_STARTED:-0}" -eq 1 ] && [ "${CCACHE_STATS_ENDED:-0}" -eq 0 ]; then
        log_ccache_stats "END"
        CCACHE_STATS_ENDED=1
    fi
}

safe_rmdir() {
    local target="$1"
    if [ -z "${HOME:-}" ] || [ "${HOME}" = "/" ]; then
        log "ERROR: HOME is '${HOME:-<unset>}', refusing to rm -rf '${target}'."
        exit 1
    fi
    case "${target}" in
        "${HOME}/"*) rm -rf "${target}" ;;
        *)
            log "ERROR: '${target}' is not under HOME '${HOME}', refusing to rm -rf."
            exit 1
            ;;
    esac
}

: "${MPLAPACK_REMOTE_WORKDIR:=${HOME}/tmp/mplapack-tier3-ubuntu2604-cxxstd-arm64}"
: "${MPLAPACK_REF:=master}"
: "${MPLAPACK_DISTRO_VERSION:=26.04}"
: "${MPLAPACK_DOCKER_BASE:=ubuntu:${MPLAPACK_DISTRO_VERSION}}"
: "${MPLAPACK_DOCKERFILE:=release/docker/tier3/Dockerfile.ubuntu-cxxstd}"
: "${MPLAPACK_DOCKER_CONTEXT:=release/docker}"
: "${MPLAPACK_IMAGE_TAG:=mplapack-tier3-ubuntu2604-cxxstd-arm64:latest}"
: "${MPLAPACK_CCACHE_DIR:=/Users/maho/.ccache}"
: "${MPLAPACK_CCACHE_MAXSIZE:=80G}"
: "${MPLAPACK_COLIMA_CPUS:=$(sysctl -n hw.ncpu 2>/dev/null || echo 10)}"
: "${MPLAPACK_COLIMA_MEMORY_GB:=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 34359738368) / 1024 / 1024 / 1024 / 2 ))}"
: "${MPLAPACK_COLIMA_DISK_GB:=100}"
: "${MPLAPACK_CONTEXT_TARBALL:=${MPLAPACK_REMOTE_WORKDIR}.context.tar.gz}"

WORKDIR="${MPLAPACK_REMOTE_WORKDIR}"
LOCKDIR="${WORKDIR}.lock"
CONTEXT_DIR="${WORKDIR}/context"
CCACHE_DIR_HOST="${MPLAPACK_CCACHE_DIR}"
CCACHE_STATS_STARTED=0
CCACHE_STATS_ENDED=0

case "${WORKDIR}" in
    "${HOME}/"*) ;;
    *)
        log "ERROR: Refusing to use WORKDIR outside HOME: ${WORKDIR}"
        exit 1
        ;;
esac

mkdir -p "$(dirname "${WORKDIR}")"

if mkdir "${LOCKDIR}" 2>/dev/null; then
    printf "%s\n" "$$" > "${LOCKDIR}/pid"
else
    old_pid=""
    [ -f "${LOCKDIR}/pid" ] && old_pid="$(cat "${LOCKDIR}/pid" 2>/dev/null || true)"
    if [ -n "${old_pid}" ] && [ "${old_pid}" != "$$" ] && kill -0 "${old_pid}" 2>/dev/null; then
        log "Another tier3-ubuntu2604-cxxstd-arm64 build is running (pid: ${old_pid}); stopping it."
        kill "${old_pid}" 2>/dev/null || true
        for _wait_i in $(seq 1 60); do
            kill -0 "${old_pid}" 2>/dev/null || break
            sleep 1
        done
        if kill -0 "${old_pid}" 2>/dev/null; then
            log "Previous build pid ${old_pid} did not stop after 60s; killing it."
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
    log_ccache_end_once
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
safe_rmdir "${WORKDIR}"
mkdir -p "${CONTEXT_DIR}" "${CCACHE_DIR_HOST}"

log "WORKDIR: ${WORKDIR}"
log "MPLAPACK_REF: ${MPLAPACK_REF}"
log "MPLAPACK_DISTRO_VERSION: ${MPLAPACK_DISTRO_VERSION}"
log "MPLAPACK_DOCKER_BASE: ${MPLAPACK_DOCKER_BASE}"
log "MPLAPACK_DOCKERFILE: ${MPLAPACK_DOCKERFILE}"
log "MPLAPACK_DOCKER_CONTEXT: ${MPLAPACK_DOCKER_CONTEXT}"
log "MPLAPACK_IMAGE_TAG: ${MPLAPACK_IMAGE_TAG}"
log "MPLAPACK_CCACHE_DIR: ${CCACHE_DIR_HOST}"
log "MPLAPACK_CCACHE_MAXSIZE: ${MPLAPACK_CCACHE_MAXSIZE}"

if [ ! -f "${MPLAPACK_CONTEXT_TARBALL}" ]; then
    log "ERROR: context tarball not found: ${MPLAPACK_CONTEXT_TARBALL}"
    exit 1
fi

tar -xzf "${MPLAPACK_CONTEXT_TARBALL}" -C "${CONTEXT_DIR}"

if [ ! -f "${CONTEXT_DIR}/${MPLAPACK_DOCKERFILE}" ]; then
    log "ERROR: Dockerfile not found in context: ${CONTEXT_DIR}/${MPLAPACK_DOCKERFILE}"
    exit 1
fi

SOURCE_BUNDLE="$(find "${CONTEXT_DIR}" -maxdepth 1 -type f -name '*_source.bundle' | head -n 1 || true)"
if [ -n "${SOURCE_BUNDLE}" ]; then
    log "MPLAPACK_SOURCE_BUNDLE: ${SOURCE_BUNDLE}"
fi
SOURCE_TARBALL="$(find "${CONTEXT_DIR}" -maxdepth 1 -type f -name 'mplapack-*.tar.*' | head -n 1 || true)"
SOURCE_METADATA="$(find "${CONTEXT_DIR}" -maxdepth 1 -type f -name 'source-metadata.txt' | head -n 1 || true)"
if [ -n "${SOURCE_TARBALL}" ]; then
    log "MPLAPACK_SOURCE_TARBALL: ${SOURCE_TARBALL}"
fi
if [ -n "${SOURCE_METADATA}" ]; then
    log "=== MPLAPACK SOURCE METADATA ==="
    cat "${SOURCE_METADATA}"
    log "=== END MPLAPACK SOURCE METADATA ==="
fi

docker build \
    --build-arg BASE="${MPLAPACK_DOCKER_BASE}" \
    --label org.mplapack.project=mplapack \
    --label org.mplapack.purpose=mplapack-qa \
    -t "${MPLAPACK_IMAGE_TAG}" \
    -f "${CONTEXT_DIR}/${MPLAPACK_DOCKERFILE}" \
    "${CONTEXT_DIR}/${MPLAPACK_DOCKER_CONTEXT}"

docker_run --rm \
    -v "${CCACHE_DIR_HOST}:/ccache:rw" \
    "${MPLAPACK_IMAGE_TAG}" \
    ccache -M "${MPLAPACK_CCACHE_MAXSIZE}"

log_ccache_start

if [ -n "${SOURCE_TARBALL}" ]; then
    docker_run --rm \
        -e MPLAPACK_REF="${MPLAPACK_REF}" \
        -e MPLAPACK_DISTRO_VERSION="${MPLAPACK_DISTRO_VERSION:-}" \
        -e CCACHE_MAXSIZE="${MPLAPACK_CCACHE_MAXSIZE}" \
        -e MPLAPACK_SOURCE_TARBALL=/source.tarball \
        -v "${SOURCE_TARBALL}:/source.tarball:ro" \
        -v "${CCACHE_DIR_HOST}:/ccache:rw" \
        "${MPLAPACK_IMAGE_TAG}"
elif [ -n "${SOURCE_BUNDLE}" ]; then
    docker_run --rm \
        -e MPLAPACK_REF="${MPLAPACK_REF}" \
        -e CCACHE_MAXSIZE="${MPLAPACK_CCACHE_MAXSIZE}" \
        -e MPLAPACK_REPO=/source.bundle \
        -v "${SOURCE_BUNDLE}:/source.bundle:ro" \
        -v "${CCACHE_DIR_HOST}:/ccache:rw" \
        "${MPLAPACK_IMAGE_TAG}"
else
    docker_run --rm \
        -e MPLAPACK_REF="${MPLAPACK_REF}" \
        -e CCACHE_MAXSIZE="${MPLAPACK_CCACHE_MAXSIZE}" \
        -v "${CCACHE_DIR_HOST}:/ccache:rw" \
        "${MPLAPACK_IMAGE_TAG}"
fi

log_ccache_end_once

log "=== ALL TIER3 DOCKER BUILD STEPS COMPLETED SUCCESSFULLY ==="
