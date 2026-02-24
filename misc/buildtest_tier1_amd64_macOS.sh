#!/bin/bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
LOG_DIR="${HOME}/mplapack_build_logs/$(LANG=C date +%Y%m%d_%H%M%S)_$$_tier1_amd64_macOS"
mkdir -p "${LOG_DIR}"

log() {
    echo "$*" | tee -a "${LOG_DIR}/summary.log"
}

# ERR trap: fires on any uncaught error outside run_step, prints line/command
trap 'echo "[ERR trap] line ${LINENO}: ${BASH_COMMAND}" \
      | tee -a "${LOG_DIR}/summary.log"; exit 1' ERR

# ---------------------------------------------------------------------------
# run_step: runs a command, logs start/end/elapsed/rc; always writes END log
# ---------------------------------------------------------------------------
run_step() {
    local name="$1"; shift
    local logfile="${LOG_DIR}/${name}.log"
    local rc t_start t_end

    log ""
    log "=== START: ${name} ==="
    LANG=C date | tee -a "${logfile}" | tee -a "${LOG_DIR}/summary.log"
    t_start=$(date +%s)

    # Disable errexit inside run_step so rc is always captured before exit
    set +e
    "$@" 2>&1 | tee -a "${logfile}"
    rc=${PIPESTATUS[0]}
    set -e

    t_end=$(date +%s)
    LANG=C date | tee -a "${logfile}" | tee -a "${LOG_DIR}/summary.log"
    log "=== END: ${name} | elapsed: $((t_end - t_start))s | rc: ${rc} ==="

    if [ "${rc}" -ne 0 ]; then
        log "FAILED at step: ${name}. Aborting."
        exit "${rc}"
    fi
}

# ---------------------------------------------------------------------------
# Safe directory / prefix removal
# ---------------------------------------------------------------------------
safe_rmdir() {
    local target="$1"
    if [ -z "${HOME:-}" ] || [ "${HOME}" = "/" ]; then
        log "ERROR: HOME is '${HOME:-<unset>}', refusing to rm -rf '${target}'."
        exit 1
    fi
    case "${target}" in
        "${HOME}/"*)
            rm -rf "${target}"
            ;;
        *)
            log "ERROR: '${target}' is not under HOME '${HOME}', refusing to rm -rf."
            exit 1
            ;;
    esac
}

# ---------------------------------------------------------------------------
# Environment snapshot
# ---------------------------------------------------------------------------
log_env() {
    {
        echo '=== Environment snapshot ==='
        LANG=C uname -a
        sw_vers 2>/dev/null || true
        echo '--- compiler versions ---'
        clang   --version 2>/dev/null || true
        clang++ --version 2>/dev/null || true
        gfortran --version 2>/dev/null || true
        brew --version 2>/dev/null || true
        echo '--- PATH ---'
        echo "${PATH}"
        echo '=== End of environment snapshot ===='
    } 2>&1 | tee "${LOG_DIR}/env.log" | tee -a "${LOG_DIR}/summary.log"
}

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PREFIX_DIR="${HOME}/MPLAPACK"
# NOTE: On case-insensitive filesystems (default macOS APFS), paths that differ
# only by letter case collide (e.g., "$HOME/mplapack" vs "$HOME/MPLAPACK").
# WORKDIR is created per-run via mktemp in the Main section to avoid collisions.

# ---------------------------------------------------------------------------
# DISTCHECK_CONFIGURE_FLAGS: feature flags only, no --prefix
# Note: Apple Silicon (arm64) does not have x87, so binary80 is x86_64 only
# ---------------------------------------------------------------------------
COMMON_FLAGS="--enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-test=yes --enable-benchmark=yes"
ARCH=$(uname -m)
case "${ARCH}" in
    x86_64|i686|i586|i386)
        DISTCHECK_CONFIGURE_FLAGS="${COMMON_FLAGS} --enable-binary80=yes"
        ;;
    *)
        DISTCHECK_CONFIGURE_FLAGS="${COMMON_FLAGS}"
        ;;
esac
export DISTCHECK_CONFIGURE_FLAGS
log "ARCH: ${ARCH}"
log "DISTCHECK_CONFIGURE_FLAGS: ${DISTCHECK_CONFIGURE_FLAGS}"

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
log_env
safe_rmdir "${PREFIX_DIR}"

# Create a unique per-run WORKDIR under $HOME/tmp.
# Include a timestamp so the directory is identifiable when kept for debugging.
: "${HOME:?HOME is not set}"
mkdir -p "${HOME}/tmp"
ts="$(LANG=C date +%Y%m%d_%H%M%S)"
WORKDIR="$(mktemp -d "${HOME}/tmp/mplapack.${ts}.XXXXXX")"
log "WORKDIR: ${WORKDIR}"

cleanup() {
    rc=$?
    cd "${HOME}" || true
    if [ "${rc}" -eq 0 ]; then
        rm -rf "${WORKDIR}" || true
    else
        echo "Keeping WORKDIR for debugging: ${WORKDIR}" | tee -a "${LOG_DIR}/summary.log"
    fi
}
trap cleanup EXIT

git clone --depth 1 --branch release/2.1 git@github.com:nakatamaho/mplapack.git "${WORKDIR}"
cd "${WORKDIR}"
git --no-pager log -1 | tee "${LOG_DIR}/git_log.log" | tee -a "${LOG_DIR}/summary.log"

run_step "reconfig"       bash misc/reconfig.macOS.sh
run_step "make"           make -j4
run_step "make_install"   make install

cp config.log "${LOG_DIR}/config.log" 2>/dev/null || true
grep "^  \$ \./configure" "${LOG_DIR}/config.log" 2>/dev/null \
    | tee -a "${LOG_DIR}/summary.log" || true

safe_rmdir "${PREFIX_DIR}"
run_step "make_distcheck" make distcheck DISTCHECK_CONFIGURE_FLAGS="${DISTCHECK_CONFIGURE_FLAGS}"
run_step "make_check"     make check -j4

log ""
log "=== ALL STEPS COMPLETED SUCCESSFULLY ==="
log "Logs: ${LOG_DIR}"
