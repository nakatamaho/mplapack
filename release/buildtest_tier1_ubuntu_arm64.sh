#!/bin/bash
set -euo pipefail
export PATH="/opt/local/bin:/opt/local/sbin:${PATH}"

log() {
    echo "$*"
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

: "${MPLAPACK_REMOTE_WORKDIR:=${HOME}/tmp/mplapack-ubuntu-arm64}"
: "${MPLAPACK_REF:=master}"
: "${MPLAPACK_DISTRO_VERSION:=24.04}"
: "${MPLAPACK_DOCKER_BASE:=ubuntu:${MPLAPACK_DISTRO_VERSION}}"
: "${MPLAPACK_DOCKERFILE:=release/docker/Dockerfile.ubuntu}"
: "${MPLAPACK_DOCKER_CONTEXT:=release/docker}"
: "${MPLAPACK_IMAGE_TAG:=mplapack-tier1-ubuntu-arm64:latest}"
: "${MPLAPACK_CCACHE_DIR:=/Users/maho/.ccache}"
: "${MPLAPACK_CCACHE_MAXSIZE:=80G}"
: "${MPLAPACK_COLIMA_CPUS:=$(sysctl -n hw.ncpu 2>/dev/null || echo 10)}"
: "${MPLAPACK_COLIMA_MEMORY_GB:=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 34359738368) / 1024 / 1024 / 1024 / 2 ))}"
: "${MPLAPACK_COLIMA_DISK_GB:=100}"
: "${MPLAPACK_CPU_MODEL_OVERRIDE:=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || true)}"
: "${MPLAPACK_RESULTS_DIR:=${MPLAPACK_REMOTE_WORKDIR}.distcheck-results}"
: "${MPLAPACK_CONTEXT_TARBALL:=${MPLAPACK_REMOTE_WORKDIR}.context.tar.gz}"

WORKDIR="${MPLAPACK_REMOTE_WORKDIR}"
LOCKDIR="${WORKDIR}.lock"
CONTEXT_DIR="${WORKDIR}/context"
RESULTS_DIR="${MPLAPACK_RESULTS_DIR}"
CCACHE_DIR_HOST="${MPLAPACK_CCACHE_DIR}"

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
        log "Another tier1-ubuntu-arm64 buildtest is running (pid: ${old_pid}); stopping it."
        kill "${old_pid}" 2>/dev/null || true
        for _wait_i in $(seq 1 60); do
            kill -0 "${old_pid}" 2>/dev/null || break
            sleep 1
        done
        if kill -0 "${old_pid}" 2>/dev/null; then
            log "Previous buildtest pid ${old_pid} did not stop after 60s; killing it."
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
mkdir -p "${CONTEXT_DIR}" "${RESULTS_DIR}" "${CCACHE_DIR_HOST}"

log "WORKDIR: ${WORKDIR}"
log "MPLAPACK_REF: ${MPLAPACK_REF}"
log "MPLAPACK_DISTRO_VERSION: ${MPLAPACK_DISTRO_VERSION}"
log "MPLAPACK_DOCKER_BASE: ${MPLAPACK_DOCKER_BASE}"
log "MPLAPACK_DOCKERFILE: ${MPLAPACK_DOCKERFILE}"
log "MPLAPACK_DOCKER_CONTEXT: ${MPLAPACK_DOCKER_CONTEXT}"
log "MPLAPACK_IMAGE_TAG: ${MPLAPACK_IMAGE_TAG}"
log "MPLAPACK_RESULTS_DIR: ${RESULTS_DIR}"
log "MPLAPACK_CCACHE_DIR: ${CCACHE_DIR_HOST}"
log "MPLAPACK_CCACHE_MAXSIZE: ${MPLAPACK_CCACHE_MAXSIZE}"
log "MPLAPACK_CPU_MODEL_OVERRIDE: ${MPLAPACK_CPU_MODEL_OVERRIDE}"

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

docker build \
    --build-arg BASE="${MPLAPACK_DOCKER_BASE}" \
    -t "${MPLAPACK_IMAGE_TAG}" \
    -f "${CONTEXT_DIR}/${MPLAPACK_DOCKERFILE}" \
    "${CONTEXT_DIR}/${MPLAPACK_DOCKER_CONTEXT}"

docker run --rm \
    -v "${CCACHE_DIR_HOST}:/ccache:rw" \
    "${MPLAPACK_IMAGE_TAG}" \
    ccache -M "${MPLAPACK_CCACHE_MAXSIZE}"

docker_source_args=()
if [ -n "${SOURCE_BUNDLE}" ]; then
    docker_source_args=(-e MPLAPACK_REPO=/source.bundle -v "${SOURCE_BUNDLE}:/source.bundle:ro")
fi

docker run --rm \
    -e MPLAPACK_REF="${MPLAPACK_REF}" \
    -e MPLAPACK_TEST_RESULTS_BASE=/results \
    -e CCACHE_MAXSIZE="${MPLAPACK_CCACHE_MAXSIZE}" \
    -e CPU_MODEL_OVERRIDE="${MPLAPACK_CPU_MODEL_OVERRIDE}" \
    "${docker_source_args[@]}" \
    -v "${CCACHE_DIR_HOST}:/ccache:rw" \
    -v "${RESULTS_DIR}:/results:rw" \
    "${MPLAPACK_IMAGE_TAG}"

log "=== ALL STEPS COMPLETED SUCCESSFULLY ==="
