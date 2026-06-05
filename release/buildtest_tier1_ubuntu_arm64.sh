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
: "${MPLAPACK_DOCKER_BASE:=ubuntu:24.04}"
: "${MPLAPACK_DOCKERFILE:=release/docker/Dockerfile.ubuntu}"
: "${MPLAPACK_IMAGE_TAG:=mplapack-tier1-ubuntu-arm64:latest}"
: "${MPLAPACK_CCACHE_DIR:=${HOME}/tmp/mplapack-ccache}"
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

if ! command -v docker >/dev/null 2>&1; then
    log "ERROR: docker command not found. PATH=${PATH}"
    exit 1
fi

docker info >/dev/null
safe_rmdir "${WORKDIR}"
mkdir -p "${CONTEXT_DIR}" "${RESULTS_DIR}" "${CCACHE_DIR_HOST}"

log "WORKDIR: ${WORKDIR}"
log "MPLAPACK_REF: ${MPLAPACK_REF}"
log "MPLAPACK_DOCKER_BASE: ${MPLAPACK_DOCKER_BASE}"
log "MPLAPACK_DOCKERFILE: ${MPLAPACK_DOCKERFILE}"
log "MPLAPACK_IMAGE_TAG: ${MPLAPACK_IMAGE_TAG}"
log "MPLAPACK_RESULTS_DIR: ${RESULTS_DIR}"

if [ ! -f "${MPLAPACK_CONTEXT_TARBALL}" ]; then
    log "ERROR: context tarball not found: ${MPLAPACK_CONTEXT_TARBALL}"
    exit 1
fi

tar -xzf "${MPLAPACK_CONTEXT_TARBALL}" -C "${CONTEXT_DIR}"

if [ ! -f "${CONTEXT_DIR}/${MPLAPACK_DOCKERFILE}" ]; then
    log "ERROR: Dockerfile not found in context: ${CONTEXT_DIR}/${MPLAPACK_DOCKERFILE}"
    exit 1
fi

docker build \
    --build-arg BASE="${MPLAPACK_DOCKER_BASE}" \
    -t "${MPLAPACK_IMAGE_TAG}" \
    -f "${CONTEXT_DIR}/${MPLAPACK_DOCKERFILE}" \
    "${CONTEXT_DIR}/release/docker"

docker run --rm \
    -e MPLAPACK_REF="${MPLAPACK_REF}" \
    -e MPLAPACK_TEST_RESULTS_BASE=/results \
    -v "${CCACHE_DIR_HOST}:/ccache:rw" \
    -v "${RESULTS_DIR}:/results:rw" \
    "${MPLAPACK_IMAGE_TAG}"

log "=== ALL STEPS COMPLETED SUCCESSFULLY ==="
