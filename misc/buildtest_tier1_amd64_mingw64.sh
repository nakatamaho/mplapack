#!/bin/bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
LOG_DIR="${HOME}/mplapack_build_logs/$(LANG=C date +%Y%m%d_%H%M%S)_$$_tier1_amd64_mingw"
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
    # Confirm target starts with $HOME to prevent accidental wide deletion
    case "${target}" in
        "${HOME}/"*)
            # make distcheck leaves read-only files in extracted tarball trees;
            # restore write permission before removal so rm -rf succeeds on rerun.
            chmod -R u+rwX "${target}" 2>/dev/null || true
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
        if [ -f /etc/os-release ]; then cat /etc/os-release; fi
        echo '--- compiler versions ---'
        x86_64-w64-mingw32-gcc   --version 2>/dev/null || true
        x86_64-w64-mingw32-g++   --version 2>/dev/null || true
        x86_64-w64-mingw32-gfortran --version 2>/dev/null || true
        echo '--- PATH ---'
        echo "${PATH}"
        echo '=== End of environment snapshot ===='
    } 2>&1 | tee "${LOG_DIR}/env.log" | tee -a "${LOG_DIR}/summary.log"
}

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PREFIX_DIR="${HOME}/MPLAPACK_MINGW"
WORKDIR="${HOME}/mplapack"

# ---------------------------------------------------------------------------
# DISTCHECK_CONFIGURE_FLAGS: feature flags only, no --prefix
# ---------------------------------------------------------------------------
COMMON_FLAGS="--enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-test=yes"
ARCH=$(uname -m)
case "${ARCH}" in
    x86_64|i686|i586|i486|i386)
        DISTCHECK_CONFIGURE_FLAGS="${COMMON_FLAGS} --enable-binary80=yes --enable-benchmark=yes"
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
safe_rmdir "${WORKDIR}"

git clone --depth 1 git@github.com:nakatamaho/mplapack.git "${WORKDIR}"
cd "${WORKDIR}"
git --no-pager log -1 | tee "${LOG_DIR}/git_log.log" | tee -a "${LOG_DIR}/summary.log"

run_step "reconfig"       bash misc/reconfig.ubuntu24.04.mingw64.sh
run_step "make"           make -j32
run_step "make_install"   make install

# Copy config.log (records actual configure invocation and detected settings)
cp config.log "${LOG_DIR}/config.log" 2>/dev/null || true
grep "^  \$ \./configure" "${LOG_DIR}/config.log" 2>/dev/null \
    | tee -a "${LOG_DIR}/summary.log" || true

# make distcheck internally does: dist -> extract -> configure -> make -> check -> install -> uninstall
# --prefix must NOT be in DISTCHECK_CONFIGURE_FLAGS (distcheck uses its own isolated prefix)
safe_rmdir "${PREFIX_DIR}"
run_step "autoreconf"     autoreconf -fi
run_step "make_distcheck" env CC="ccache x86_64-w64-mingw32-gcc" CXX="ccache x86_64-w64-mingw32-g++" FC="ccache x86_64-w64-mingw32-gfortran" \
                          make distcheck LOG_COMPILER=wine MAKEFLAGS="-j32" DISTCHECK_CONFIGURE_FLAGS="--host=x86_64-w64-mingw32 ${DISTCHECK_CONFIGURE_FLAGS}"

log ""
log "=== ALL STEPS COMPLETED SUCCESSFULLY ==="
log "Logs: ${LOG_DIR}"
