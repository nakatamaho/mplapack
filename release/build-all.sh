#!/bin/bash
# build-all.sh
# Multi-environment build and test script for MPLAPACK release validation

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

mkdir -p "$LOGDIR"

# Initialize buildx builder if not exists
docker buildx inspect mplapack-builder >/dev/null 2>&1 || \
    docker buildx create --name mplapack-builder --use >/dev/null

log() {
    echo "[$(date +%H:%M:%S)] $*"
}

# Check if NVIDIA GPU is available
check_gpu() {
    if [[ "$USE_GPU" == "no" ]]; then
        return 1
    elif [[ "$USE_GPU" == "yes" ]]; then
        return 0
    else
        # auto-detect
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
    local is_cuda=false
    local gpu_flag=""

    # Detect CUDA build
    if [[ "$name" == *cuda* ]]; then
        is_cuda=true
        if check_gpu; then
            gpu_flag="--gpus all"
            log "  (GPU detected, will use --gpus all for tests)"
        else
            log "  (No GPU detected, tests may be limited)"
        fi
    fi

    # Copy tarball to docker context if tarball test
    if [[ "$source_type" == "tarball" ]]; then
        if [[ -z "$TARBALL" || ! -f "$TARBALL" ]]; then
            echo "$name,$arch_short,$base,build,SKIPPED,0,$source_type"
            return 0
        fi
        cp "$TARBALL" "$DOCKER_DIR/"
    fi

    # Build image
    if ! docker buildx build \
        --platform "$arch" \
        --build-arg BASE="$base" \
        --load \
        -t "$tag" \
        -f "$DOCKER_DIR/$dockerfile" \
        "$PROJECT_ROOT" \
        > "${logprefix}_build.log" 2>&1; then

        local elapsed=$(($(date +%s) - start))
        echo "$name,$arch_short,$base,build,FAILED,$elapsed,$source_type"
        return 1
    fi

    # Run tests (with GPU if available for CUDA builds)
    if ! docker run --rm --platform "$arch" $gpu_flag "$tag" sh -c "$test_cmd" \
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
    echo "Logs:  $LOGDIR/"

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
