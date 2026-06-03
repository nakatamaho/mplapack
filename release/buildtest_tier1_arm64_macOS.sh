#!/bin/bash
set -euo pipefail
export PATH="/opt/local/bin:/opt/local/sbin:${PATH}"

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
LOG_DIR="${HOME}/mplapack_build_logs/$(LANG=C date +%Y%m%d_%H%M%S)_$$_tier1_arm64_macOS"
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
            # restore write permission before removal so rm -rf succeeds.
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
        sw_vers 2>/dev/null || true
        echo '--- compiler versions ---'
        gcc-mp-14   --version 2>/dev/null || true
        g++-mp-14 --version 2>/dev/null || true
        gfortran-mp-14 --version 2>/dev/null || true
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

# Create a unique per-run WORKDIR under $HOME/tmp.
# Include a timestamp so the directory is identifiable when kept for debugging.
: "${HOME:?HOME is not set}"

realpath_safe() {
    # Resolve to a canonical absolute path (no symlinks).
    # Prefer realpath(1), then python3, then perl.
    local p="$1"
    if command -v realpath >/dev/null 2>&1; then
        realpath "$p"
        return
    fi
    if command -v python3 >/dev/null 2>&1; then
        python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "$p"
        return
    fi
    if command -v perl >/dev/null 2>&1; then
        perl -MCwd=realpath -e 'print realpath($ARGV[0])' "$p"
        return
    fi
    # Last resort: best-effort physical path for an existing directory.
    (cd "$p" 2>/dev/null && pwd -P) || return 1
}

# Use a fixed WORKDIR to improve ccache hit rate (avoid per-run path variability).
WORKDIR="${MPLAPACK_REMOTE_WORKDIR:-${HOME}/tmp/mplapack}"
case "${WORKDIR}" in
    "~/"*) WORKDIR="${HOME}/${WORKDIR#~/}" ;;
esac
BASE_TMP_DIR="$(dirname "${WORKDIR}")"

# Prevent concurrent runs from clobbering each other.
# Use an atomic mkdir lock to avoid relying on flock(1) availability.
LOCKDIR="${WORKDIR}.lock"
# Safety guard: refuse to operate outside HOME.
case "${WORKDIR}" in
    "${HOME}/"*) ;;
    *)
        echo "ERROR: Refusing to use WORKDIR outside HOME: ${WORKDIR}" >&2
        exit 1
        ;;
esac

# Symlink guard: refuse if BASE_TMP_DIR or WORKDIR is a symlink (prevents deleting through links).
if [ -L "${BASE_TMP_DIR}" ]; then
    echo "ERROR: Refusing symlink BASE_TMP_DIR: ${BASE_TMP_DIR}" >&2
    exit 1
fi
if [ -L "${WORKDIR}" ]; then
    echo "ERROR: Refusing symlink WORKDIR: ${WORKDIR}" >&2
    exit 1
fi

# Ensure WORKDIR exists before realpath resolution (first-run case: directory not yet created).
mkdir -p "${WORKDIR}"

# Realpath guard: ensure resolved WORKDIR matches resolved expected path.
EXPECTED_WORKDIR="${WORKDIR}"
WORKDIR_REAL="$(realpath_safe "${WORKDIR}")" || { echo "ERROR: Failed to resolve realpath: ${WORKDIR}" >&2; exit 1; }
EXPECTED_REAL="$(realpath_safe "${EXPECTED_WORKDIR}")" || { echo "ERROR: Failed to resolve realpath: ${EXPECTED_WORKDIR}" >&2; exit 1; }
if [ "${WORKDIR_REAL}" != "${EXPECTED_REAL}" ]; then
    echo "ERROR: WORKDIR realpath mismatch: '${WORKDIR_REAL}' != '${EXPECTED_REAL}'" >&2
    exit 1
fi

# Acquire lock (with stale-lock detection).
if mkdir "${LOCKDIR}" 2>/dev/null; then
    printf '%s\n' "$$" > "${LOCKDIR}/pid"
else
    if [ -f "${LOCKDIR}/pid" ]; then
        old_pid="$(cat "${LOCKDIR}/pid" 2>/dev/null || true)"
    else
        old_pid=""
    fi

    if [ -n "${old_pid}" ] && kill -0 "${old_pid}" 2>/dev/null; then
        echo "ERROR: Another buildtest instance is running (lock: ${LOCKDIR}, pid: ${old_pid})" >&2
        exit 1
    fi

    echo "WARNING: Stale lock detected; removing it: ${LOCKDIR}" >&2
    rm -rf "${LOCKDIR}"
    if ! mkdir "${LOCKDIR}" 2>/dev/null; then
        echo "ERROR: Failed to acquire lock after removing stale lock: ${LOCKDIR}" >&2
        exit 1
    fi
    printf '%s\n' "$$" > "${LOCKDIR}/pid"
fi

# Always release the lock on exit.
cleanup_lock() {
    # make distcheck intentionally makes extracted tree read-only;
    # restore write bits so the directory can be removed on any exit path.
    chmod -R u+rwX "${WORKDIR}" 2>/dev/null || true
    rm -rf "${LOCKDIR}"
}
trap cleanup_lock EXIT INT TERM HUP

mkdir -p "${WORKDIR}"
# make distcheck leaves read-only files in the extracted tarball tree;
# restore write permission before removal so rm -rf succeeds on rerun.
chmod -R u+rwX "${WORKDIR}" 2>/dev/null || true
find "${WORKDIR}" -mindepth 1 -maxdepth 1 -exec rm -rf -- {} +

log "WORKDIR: ${WORKDIR}"

git clone --depth 1 git@github.com:nakatamaho/mplapack.git "${WORKDIR}"
cd "${WORKDIR}"
git --no-pager log -1 | tee "${LOG_DIR}/git_log.log" | tee -a "${LOG_DIR}/summary.log"

run_step "reconfig"       bash misc/reconfig.macOS.sh
run_step "make"           make -j4
run_step "make_install"   make install

# Copy config.log (records actual configure invocation and detected settings)
cp config.log "${LOG_DIR}/config.log" 2>/dev/null || true
grep "^  \$ \./configure" "${LOG_DIR}/config.log" 2>/dev/null \
    | tee -a "${LOG_DIR}/summary.log" || true

# make distcheck internally does: dist -> extract -> configure -> make -> check -> install -> uninstall
# --prefix must NOT be in DISTCHECK_CONFIGURE_FLAGS (distcheck uses its own isolated prefix)
safe_rmdir "${PREFIX_DIR}"
run_step "autoreconf"     autoreconf -fi
run_step "make_distcheck" env CC="ccache gcc" CXX="ccache g++" FC="ccache gfortran" \
                          make distcheck MAKEFLAGS="-j4" DISTCHECK_CONFIGURE_FLAGS="${DISTCHECK_CONFIGURE_FLAGS}"

log ""
log "=== ALL STEPS COMPLETED SUCCESSFULLY ==="
log "Logs: ${LOG_DIR}"
