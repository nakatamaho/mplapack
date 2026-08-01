#!/bin/bash
set -euo pipefail
export PATH="/opt/local/bin:/opt/local/sbin:${PATH}"

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
LOG_DIR="${HOME}/mplapack_build_logs/$(LANG=C LC_ALL=C date +%Y%m%d_%H%M%S)_$$_tier1_macos_arm64"
mkdir -p "${LOG_DIR}"
CCACHE_STATS_STARTED=0
CCACHE_STATS_ENDED=0

log() {
    echo "$*" | tee -a "${LOG_DIR}/summary.log"
}

log_ccache_stats() {
    local label="$1"

    log "=== CCACHE STATS (${label}) ==="
    if command -v ccache >/dev/null 2>&1; then
        ccache -s 2>&1 | tee -a "${LOG_DIR}/summary.log"
    else
        log "ccache command not found"
    fi
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
    LANG=C LC_ALL=C date | tee -a "${logfile}" | tee -a "${LOG_DIR}/summary.log"
    t_start=$(date +%s)

    # Disable errexit inside run_step so rc is always captured before exit
    set +e
    "$@" 2>&1 | tee -a "${logfile}"
    rc=${PIPESTATUS[0]}
    set -e

    t_end=$(date +%s)
    LANG=C LC_ALL=C date | tee -a "${logfile}" | tee -a "${LOG_DIR}/summary.log"
    log "=== END: ${name} | elapsed: $((t_end - t_start))s | rc: ${rc} ==="

    if [ "${rc}" -ne 0 ]; then
        log "FAILED at step: ${name}. Aborting."
        exit "${rc}"
    fi
}

get_make_jobs() {
    if [ -n "${MPLAPACK_MAKE_JOBS:-}" ]; then
        printf '%s\n' "${MPLAPACK_MAKE_JOBS}"
        return
    fi
    if command -v sysctl >/dev/null 2>&1; then
        jobs="$(sysctl -n hw.physicalcpu 2>/dev/null || true)"
        if [ -n "${jobs}" ] && [ "${jobs}" -gt 0 ] 2>/dev/null; then
            printf '%s\n' "${jobs}"
            return
        fi
    fi
    if command -v getconf >/dev/null 2>&1; then
        jobs="$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)"
        if [ -n "${jobs}" ] && [ "${jobs}" -gt 0 ] 2>/dev/null; then
            printf '%s\n' "${jobs}"
            return
        fi
    fi
    printf '4\n'
}

select_make_cmd() {
    if [ -n "${MPLAPACK_MAKE:-}" ]; then
        printf '%s\n' "${MPLAPACK_MAKE}"
        return
    fi
    if command -v gmake >/dev/null 2>&1; then
        printf 'gmake\n'
        return
    fi
    echo "ERROR: gmake not found; macOS tier1 parallel distcheck requires GNU make. Install gmake or set MPLAPACK_MAKE." >&2
    exit 1
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
        gcc-mp-15   --version 2>/dev/null || true
        g++-mp-15 --version 2>/dev/null || true
        gfortran-mp-15 --version 2>/dev/null || true
        brew --version 2>/dev/null || true
        echo '--- PATH ---'
        echo "${PATH}"
        echo '=== End of environment snapshot ===='
    } 2>&1 | tee "${LOG_DIR}/env.log" | tee -a "${LOG_DIR}/summary.log"
}

# ---------------------------------------------------------------------------
# Test result collection
# ---------------------------------------------------------------------------
get_mplapack_version() {
    local major minor patch
    major="$(sed -n 's/^m4_define(\[MPLAPACK_VER_MAJOR\], \[\([^]]*\)\])/\1/p' configure.ac | tr -d '[:space:]')"
    minor="$(sed -n 's/^m4_define(\[MPLAPACK_VER_MINOR\], \[\([^]]*\)\])/\1/p' configure.ac | tr -d '[:space:]')"
    patch="$(sed -n 's/^m4_define(\[MPLAPACK_VER_PATCH\], \[\([^]]*\)\])/\1/p' configure.ac | tr -d '[:space:]')"
    if [ -n "${major}" ] && [ -n "${minor}" ] && [ -n "${patch}" ]; then
        printf '%s.%s.%s\n' "${major}" "${minor}" "${patch}"
    fi
}

collect_test_results() {
    local version="${MPLAPACK_RESULTS_VERSION:-}"
    local staging_root="${MPLAPACK_TEST_RESULTS_STAGING:-}"
    local local_root="${MPLAPACK_RESULTS_ROOT:-${HOME}/mplapack}"
    local suite src_base src_dir dst_dir outdir copied

    if [ -z "${version}" ]; then
        version="$(get_mplapack_version)"
    fi
    if [ -z "${version}" ]; then
        log "ERROR: Failed to determine MPLAPACK results version."
        exit 1
    fi

    if [ -z "${staging_root}" ]; then
        staging_root="${WORKDIR}.distcheck-results/${version}"
    fi

    case "${local_root}" in
        "${HOME}"|"${HOME}/"*) ;;
        *)
            log "ERROR: Refusing MPLAPACK_RESULTS_ROOT outside HOME: ${local_root}"
            exit 1
            ;;
    esac

    for suite in lin eig; do
        src_base="${staging_root}/${suite}/results"
        dst_dir="${local_root}/mplapack/test/${suite}/results/${version}"
        mkdir -p "${dst_dir}"
        copied=0

        if [ -d "${src_base}" ]; then
            while IFS= read -r src_dir; do
                [ -d "${src_dir}" ] || continue
                outdir="$(basename "${src_dir}")"
                if command -v rsync >/dev/null 2>&1; then
                    rsync -a "${src_dir}/" "${dst_dir}/${outdir}/"
                else
                    mkdir -p "${dst_dir}/${outdir}"
                    cp -a "${src_dir}/." "${dst_dir}/${outdir}/"
                fi
                copied=$((copied + 1))
                log "Collected ${suite} results: ${src_dir} -> ${dst_dir}/${outdir}"
            done < <(find "${src_base}" -mindepth 1 -maxdepth 1 -type d -print)
        fi

        if [ "${copied}" -eq 0 ]; then
            log "No ${suite} results found under ${src_base}"
        fi
    done
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
    x86_64|i686|i586|i486|i386)
        DISTCHECK_CONFIGURE_FLAGS="${COMMON_FLAGS} --enable-binary80=yes"
        ;;
    *)
        DISTCHECK_CONFIGURE_FLAGS="${COMMON_FLAGS}"
        ;;
esac
export DISTCHECK_CONFIGURE_FLAGS
log "ARCH: ${ARCH}"
log "DISTCHECK_CONFIGURE_FLAGS: ${DISTCHECK_CONFIGURE_FLAGS}"
MAKE_JOBS="$(get_make_jobs)"
if [ -n "${MPLAPACK_MAKE_JOBS:-}" ]; then
    log "MAKE_JOBS (MPLAPACK_MAKE_JOBS override): ${MAKE_JOBS}"
else
    log "MAKE_JOBS (physical cores): ${MAKE_JOBS}"
fi
MAKE_CMD="$(select_make_cmd)"
export MAKE="${MAKE_CMD}"
log "MAKE_CMD: ${MAKE_CMD}"

: "${MPLAPACK_CCACHE_DIR:=/Users/maho/.ccache}"
: "${MPLAPACK_CCACHE_MAXSIZE:=80G}"
export CCACHE_DIR="${MPLAPACK_CCACHE_DIR}"
mkdir -p "${CCACHE_DIR}"
if command -v ccache >/dev/null 2>&1; then
    ccache -M "${MPLAPACK_CCACHE_MAXSIZE}"
fi
log "CCACHE_DIR: ${CCACHE_DIR}"
log "MPLAPACK_CCACHE_MAXSIZE: ${MPLAPACK_CCACHE_MAXSIZE}"
log_ccache_start
trap log_ccache_end_once EXIT

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
log_env

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

# Acquire lock. If another buildtest is running, stop it and take over.
if mkdir "${LOCKDIR}" 2>/dev/null; then
    printf "%s\n" "$$" > "${LOCKDIR}/pid"
else
    if [ -f "${LOCKDIR}/pid" ]; then
        old_pid="$(cat "${LOCKDIR}/pid" 2>/dev/null || true)"
    else
        old_pid=""
    fi

    if [ -n "${old_pid}" ] && [ "${old_pid}" != "$$" ] && kill -0 "${old_pid}" 2>/dev/null; then
        log "Another buildtest instance is running (lock: ${LOCKDIR}, pid: ${old_pid}); stopping it."
        kill "${old_pid}" 2>/dev/null || true
        for _wait_i in $(seq 1 60); do
            if ! kill -0 "${old_pid}" 2>/dev/null; then
                break
            fi
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
    if ! mkdir "${LOCKDIR}" 2>/dev/null; then
        echo "ERROR: Failed to acquire lock after removing stale lock: ${LOCKDIR}" >&2
        exit 1
    fi
    printf "%s\n" "$$" > "${LOCKDIR}/pid"
fi

# Always release the lock on exit.
cleanup_lock() {
    log_ccache_end_once
    # make distcheck intentionally makes extracted tree read-only;
    # restore write bits so the directory can be removed on any exit path.
    chmod -R u+rwX "${WORKDIR}" 2>/dev/null || true
    rm -rf "${LOCKDIR}"
}
trap cleanup_lock EXIT INT TERM HUP

safe_rmdir "${PREFIX_DIR}"
mkdir -p "${WORKDIR}"
# make distcheck leaves read-only files in the extracted tarball tree;
# restore write permission before removal so rm -rf succeeds on rerun.
chmod -R u+rwX "${WORKDIR}" 2>/dev/null || true
find "${WORKDIR}" -mindepth 1 -maxdepth 1 -exec rm -rf -- {} +

log "WORKDIR: ${WORKDIR}"
log "MPLAPACK_REF: ${MPLAPACK_REF:-master}"
log "MPLAPACK_SOURCE_MODE: ${MPLAPACK_SOURCE_MODE:-dist}"

SOURCE_KIND=git
if [ "${MPLAPACK_SOURCE_MODE:-dist}" = "dist" ]; then
    if [ -z "${MPLAPACK_CONTEXT_TARBALL:-}" ] || [ ! -f "${MPLAPACK_CONTEXT_TARBALL}" ]; then
        log "ERROR: context tarball not found: ${MPLAPACK_CONTEXT_TARBALL:-<unset>}"
        exit 1
    fi
    CONTEXT_DIR="${WORKDIR}.context"
    safe_rmdir "${CONTEXT_DIR}"
    mkdir -p "${CONTEXT_DIR}"
    tar -xzf "${MPLAPACK_CONTEXT_TARBALL}" -C "${CONTEXT_DIR}"
    SOURCE_TARBALL="$(find "${CONTEXT_DIR}" -maxdepth 1 -type f -name 'mplapack-*.tar.*' | head -n 1 || true)"
    SOURCE_METADATA="$(find "${CONTEXT_DIR}" -maxdepth 1 -type f -name 'source-metadata.txt' | head -n 1 || true)"
    if [ -z "${SOURCE_TARBALL}" ]; then
        log "ERROR: no mplapack source tarball found in context: ${CONTEXT_DIR}"
        exit 1
    fi
    log "MPLAPACK_SOURCE_TARBALL: ${SOURCE_TARBALL}"
    if [ -n "${SOURCE_METADATA}" ]; then
        log "=== MPLAPACK SOURCE METADATA ==="
        cat "${SOURCE_METADATA}" | tee -a "${LOG_DIR}/summary.log"
        log "=== END MPLAPACK SOURCE METADATA ==="
    fi
    tar xf "${SOURCE_TARBALL}" -C "${WORKDIR}" --strip-components 1
    cd "${WORKDIR}"
    SOURCE_KIND=tarball
else
    git clone git@github.com:nakatamaho/mplapack.git "${WORKDIR}"
    cd "${WORKDIR}"
    git checkout "${MPLAPACK_REF:-master}"
    git --no-pager log -1 | tee "${LOG_DIR}/git_log.log" | tee -a "${LOG_DIR}/summary.log"
fi

RESULTS_VERSION="$(get_mplapack_version)"
if [ -z "${RESULTS_VERSION}" ]; then
    log "ERROR: Failed to determine MPLAPACK results version."
    exit 1
fi
DISTCHECK_RESULTS_STAGING="${WORKDIR}.distcheck-results/${RESULTS_VERSION}"
safe_rmdir "${WORKDIR}.distcheck-results"
mkdir -p "${DISTCHECK_RESULTS_STAGING}"
export MPLAPACK_RESULTS_VERSION="${RESULTS_VERSION}"
export MPLAPACK_TEST_RESULTS_STAGING="${DISTCHECK_RESULTS_STAGING}"
log "MPLAPACK_RESULTS_VERSION: ${MPLAPACK_RESULTS_VERSION}"
log "MPLAPACK_TEST_RESULTS_STAGING: ${MPLAPACK_TEST_RESULTS_STAGING}"

if [ "$SOURCE_KIND" = "git" ]; then
    run_step "reconfig"       bash misc/reconfig.macOS.sh
else
    run_step "configure"      env CC="ccache gcc-mp-15" CXX="ccache g++-mp-15" FC="ccache gfortran-mp-15" \
                              ./configure --prefix="${PREFIX_DIR}" ${DISTCHECK_CONFIGURE_FLAGS}
fi
run_step "make"           "${MAKE_CMD}" -j"${MAKE_JOBS}"
run_step "make_install"   "${MAKE_CMD}" install
run_step "check_installed_examples" bash release/check-installed-examples.sh "${PREFIX_DIR}" Makefile.macos "${MAKE_JOBS}"
run_step "check_installed_benchmarks" bash release/check-installed-benchmarks.sh "${PREFIX_DIR}"

# Copy config.log (records actual configure invocation and detected settings)
cp config.log "${LOG_DIR}/config.log" 2>/dev/null || true
grep "^  \$ \./configure" "${LOG_DIR}/config.log" 2>/dev/null \
    | tee -a "${LOG_DIR}/summary.log" || true

# make distcheck internally does: dist -> extract -> configure -> make -> check -> install -> uninstall
# --prefix must NOT be in DISTCHECK_CONFIGURE_FLAGS (distcheck uses its own isolated prefix)
safe_rmdir "${PREFIX_DIR}"
if [ "$SOURCE_KIND" = "git" ]; then
    run_step "autoreconf"     autoreconf -fi
else
    log "Using distributed configure files from source tarball; skipping autoreconf."
fi
run_step "make_distcheck" env CC="ccache gcc" CXX="ccache g++" FC="ccache gfortran-mp-15" \
                          "${MAKE_CMD}" -j"${MAKE_JOBS}" distcheck DISTCHECK_CONFIGURE_FLAGS="${DISTCHECK_CONFIGURE_FLAGS}"
run_step "collect_test_results" collect_test_results

log_ccache_end_once

log ""
log "=== ALL STEPS COMPLETED SUCCESSFULLY ==="
log "Logs: ${LOG_DIR}"
