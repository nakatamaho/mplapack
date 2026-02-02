#!/bin/bash
# build-all.sh
# Multi-environment build and test script for MPLAPACK release validation
# Uses host-mounted ccache directory

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

CONF_FILE="${CONF_FILE:-$SCRIPT_DIR/build-matrix.conf}"
LOGDIR="${LOGDIR:-$SCRIPT_DIR/logs/$(date +%Y%m%d_%H%M%S)}"
DOCKER_DIR="$SCRIPT_DIR/docker"
TARBALL="${TARBALL:-}"
PHASE="${PHASE:-all}"
FILTER_NAME="${FILTER_NAME:-}"
FILTER_ARCH="${FILTER_ARCH:-}"
USE_GPU="${USE_GPU:-auto}"
CCACHE_MAXSIZE="${CCACHE_MAXSIZE:-40G}"

# Host directories
HOST_WORK_DIR="${HOST_WORK_DIR:-$PROJECT_ROOT}"
HOST_CCACHE_DIR="${HOST_CCACHE_DIR:-/work/ccache}"

mkdir -p "$LOGDIR"
mkdir -p "$HOST_CCACHE_DIR"

# Initialize buildx builder if not exists
docker buildx inspect mplapack-builder >/dev/null 2>&1 || \
    docker buildx create --name mplapack-builder --use >/dev/null

log() {
    echo "[$(date +%H:%M:%S)] $*"
}

setup_ccache() {
    local conf="$HOST_CCACHE_DIR/ccache.conf"

    log "Configuring ccache: dir=$HOST_CCACHE_DIR max_size=$CCACHE_MAXSIZE"

    # Prefer using ccache itself if available on the host.
    if command -v ccache >/dev/null 2>&1; then
        CCACHE_DIR="$HOST_CCACHE_DIR" ccache -M "$CCACHE_MAXSIZE" >/dev/null 2>&1 || true
    fi

    # Ensure the cache-specific config file exists and contains max_size.
    if [[ -f "$conf" ]]; then
        if grep -qE '^[[:space:]]*max_size[[:space:]]*=' "$conf"; then
            sed -i.bak -E "s|^[[:space:]]*max_size[[:space:]]*=.*|max_size = ${CCACHE_MAXSIZE}|" "$conf"
        else
            printf "max_size = %s\n" "$CCACHE_MAXSIZE" >> "$conf"
        fi
    else
        printf "max_size = %s\n" "$CCACHE_MAXSIZE" > "$conf"
    fi
}

setup_ccache


# Detect host Docker platform (best-effort)
detect_host_platform() {
    local m
    m="$(uname -m)"
    case "$m" in
        x86_64|amd64) echo "linux/amd64" ;;
        aarch64|arm64) echo "linux/arm64" ;;
        armv7l|armv7*) echo "linux/arm/v7" ;;
        armv6l|armv6*) echo "linux/arm/v6" ;;
        ppc64le) echo "linux/ppc64le" ;;
        s390x) echo "linux/s390x" ;;
        *) echo "" ;;
    esac
}

HOST_PLATFORM="$(detect_host_platform)"
EMULATE="${EMULATE:-auto}"   # auto|yes|no

# Ensure binfmt/qemu is installed for cross-arch runs (Docker Desktop often has this already)
ensure_binfmt() {
    local target="$1"
    if [[ "$EMULATE" == "no" ]]; then
        return 1
    fi

    # If we can't detect host platform, play it safe and attempt to install binfmt in auto mode.
    if [[ -n "${HOST_PLATFORM:-}" && "$target" == "$HOST_PLATFORM" ]]; then
        return 0
    fi

    # Quick sanity check: can we run a trivial container for target arch?
    if docker run --rm --platform "$target" alpine:3.19 true >/dev/null 2>&1; then
        return 0
    fi

    if [[ "$EMULATE" == "auto" || "$EMULATE" == "yes" ]]; then
        log "Enabling binfmt/qemu for cross-arch ($HOST_PLATFORM -> $target)"
        docker run --privileged --rm tonistiigi/binfmt:latest --install all >/dev/null 2>&1 || true
        # Re-test after install
        docker run --rm --platform "$target" alpine:3.19 true >/dev/null 2>&1
    else
        return 1
    fi
}

# Check if NVIDIA GPU is available
check_gpu() {
    if [[ "$USE_GPU" == "no" ]]; then
        return 1
    elif [[ "$USE_GPU" == "yes" ]]; then
        return 0
    else
        docker run --rm --gpus all nvidia/cuda:12.4.0-base-ubuntu22.04 nvidia-smi >/dev/null 2>&1
    fi
}

# Generate release tarball using autotools
make_dist() {
    log "=== Phase 0: make dist ==="
    local distdir="$LOGDIR/dist"
    mkdir -p "$distdir"

    cd "$PROJECT_ROOT"

    if [[ ! -f configure.ac ]]; then
        log "ERROR: configure.ac not found in $PROJECT_ROOT"
        exit 1
    fi

    log "Running autoreconf..."
    autoreconf -fi > "$distdir/autoreconf.log" 2>&1

    log "Running configure..."
    ./configure > "$distdir/configure.log" 2>&1

    log "Running make dist..."
    make dist > "$distdir/make_dist.log" 2>&1

    cp mplapack-*.tar.gz "$distdir/"
    TARBALL=$(ls "$distdir"/mplapack-*.tar.gz | head -1)
    log "Created: $TARBALL"

    cd "$SCRIPT_DIR"
}

# Build and test single configuration
build_one() {
    local name=$1
    local base=$2
    local arch=$3
    local dockerfile=$4
    local test_cmd=$5
    local source_type=$6

    local arch_short=${arch##*/}
    local tag="mplapack-${name}-${arch_short}"
    local logprefix="$LOGDIR/${name}_${arch_short}"
    local start=$(date +%s)
    local gpu_flag=""


# Cross-architecture support (qemu/binfmt)
if [[ -n "${HOST_PLATFORM:-}" && "$arch" != "$HOST_PLATFORM" ]]; then
    if ! ensure_binfmt "$arch"; then
        log "  (Cross-arch disabled or binfmt unavailable; skipping $name / ${arch_short})"
        echo "$name,$arch_short,$base,build,SKIPPED,0,$source_type"
        return 0
    fi
    log "  (Cross-arch enabled: $HOST_PLATFORM -> $arch)"
fi

# Detect CUDA build
    if [[ "$name" == *cuda* ]]; then
        if check_gpu; then
            gpu_flag="--gpus all"
            log "  (GPU detected)"
        else
            log "  (No GPU detected, tests may be limited)"
        fi
    fi

    # Prepare work directory
    local work_mount="$HOST_WORK_DIR"

    # For tarball tests, create temp directory with tarball
    if [[ "$source_type" == "tarball" ]]; then
        if [[ -z "$TARBALL" || ! -f "$TARBALL" ]]; then
            echo "$name,$arch_short,$base,build,SKIPPED,0,$source_type"
            return 0
        fi
        local tmpdir="$LOGDIR/workdir_${name}_${arch_short}"
        mkdir -p "$tmpdir"
        cp "$TARBALL" "$tmpdir/"
        work_mount="$tmpdir"
    fi

    # Build image (environment setup only, no compilation)
    if ! docker buildx build \
        --platform "$arch" \
        --build-arg BASE="$base" \
        --load \
        -t "$tag" \
        -f "$DOCKER_DIR/$dockerfile" \
        "$DOCKER_DIR" \
        > "${logprefix}_image.log" 2>&1; then

        local elapsed=$(($(date +%s) - start))
        echo "$name,$arch_short,$base,image,FAILED,$elapsed,$source_type"
        return 1
    fi

    # Run build with host directories mounted
    if ! docker run --rm \
        --platform "$arch" \
        -v "$work_mount:/work:rw" \
        -v "$HOST_CCACHE_DIR:/ccache:rw" \
        -e CCACHE_DIR=/ccache \
        -e CCACHE_MAXSIZE="$CCACHE_MAXSIZE" \
        $gpu_flag \
        "$tag" \
        > "${logprefix}_build.log" 2>&1; then

        local elapsed=$(($(date +%s) - start))
        echo "$name,$arch_short,$base,build,FAILED,$elapsed,$source_type"
        return 1
    fi

    # Run tests
    if ! docker run --rm \
        --platform "$arch" \
        -v "$work_mount:/work:rw" \
        -v "$HOST_CCACHE_DIR:/ccache:rw" \
        -e CCACHE_DIR=/ccache \
        -e CCACHE_MAXSIZE="$CCACHE_MAXSIZE" \
        $gpu_flag \
        "$tag" \
        sh -c "$test_cmd" \
        > "${logprefix}_test.log" 2>&1; then

        local elapsed=$(($(date +%s) - start))
        echo "$name,$arch_short,$base,test,FAILED,$elapsed,$source_type"
        return 1
    fi

    local elapsed=$(($(date +%s) - start))
    echo "$name,$arch_short,$base,test,OK,$elapsed,$source_type"
}

# Process build matrix
run_matrix() {
    echo "name,arch,base,stage,result,elapsed,source_type"

    while IFS='|' read -r name base archs dockerfile test_cmd source_type; do
        # Skip comments and empty lines
        [[ -z "$name" || "$name" =~ ^# ]] && continue

        # Apply name filter
        [[ -n "$FILTER_NAME" && "$name" != *"$FILTER_NAME"* ]] && continue

        # Set defaults
        source_type="${source_type:-branch}"
        test_cmd="${test_cmd:-make check}"

        # Apply phase filter
        [[ "$PHASE" != "all" && "$source_type" != "$PHASE" ]] && continue

        # Process each architecture
        IFS=',' read -ra arch_array <<< "$archs"
        for arch in "${arch_array[@]}"; do
            local arch_short=${arch##*/}

            # Apply arch filter
            [[ -n "$FILTER_ARCH" && "$arch_short" != "$FILTER_ARCH" ]] && continue

            log "START $name / $arch_short ($source_type)"
            build_one "$name" "$base" "$arch" "$dockerfile" "$test_cmd" "$source_type" || true
        done
    done < "$CONF_FILE"
}

# Display summary of results
show_summary() {
    [[ -f "$LOGDIR/results.csv" ]] || return

    echo ""
    log "=== SUMMARY ==="
    column -t -s, "$LOGDIR/results.csv"

    local total=$(( $(wc -l < "$LOGDIR/results.csv") - 1 ))
    local passed=$(grep -c ',OK,' "$LOGDIR/results.csv" || echo 0)
    local failed=$(grep -c ',FAILED,' "$LOGDIR/results.csv" || echo 0)
    local skipped=$(grep -c ',SKIPPED,' "$LOGDIR/results.csv" || echo 0)

    echo ""
    echo "Total: $total | Passed: $passed | Failed: $failed | Skipped: $skipped"
    echo "Logs:   $LOGDIR/"
    echo "ccache: $HOST_CCACHE_DIR"

    [[ $failed -eq 0 ]] || exit 1
}

# Main entry point
case "${1:-}" in
    dist)
        make_dist
        ;;
    branch)
        PHASE=branch run_matrix | tee "$LOGDIR/results.csv"
        show_summary
        ;;
    tarball)
        [[ -z "$TARBALL" ]] && make_dist
        PHASE=tarball run_matrix | tee "$LOGDIR/results.csv"
        show_summary
        ;;
    *)
        # Full cycle: branch -> dist -> tarball
        log "=== Phase 1: Branch tests ==="
        PHASE=branch run_matrix | tee "$LOGDIR/results_branch.csv"

        echo ""
        make_dist

        echo ""
        log "=== Phase 2: Tarball tests ==="
        PHASE=tarball run_matrix | tee "$LOGDIR/results_tarball.csv"

        # Merge results
        cat "$LOGDIR/results_branch.csv" > "$LOGDIR/results.csv"
        tail -n +2 "$LOGDIR/results_tarball.csv" >> "$LOGDIR/results.csv"

        show_summary
        ;;
esac
