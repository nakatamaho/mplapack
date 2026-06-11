#!/bin/bash
# build-all.sh
# Multi-environment build and test script for MPLAPACK release validation
# Uses host-mounted ccache directory; supports cross-arch via binfmt/qemu

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

CONF_FILE="${CONF_FILE:-$SCRIPT_DIR/build-matrix.conf}"
LOGDIR="${LOGDIR:-$SCRIPT_DIR/logs/$(date +%Y%m%d_%H%M%S)}"
SUCCESS_DIR="${SUCCESS_DIR:-$SCRIPT_DIR/success}"
DOCKER_DIR="$SCRIPT_DIR/docker"
TARBALL="${TARBALL:-}"
PHASE="${PHASE:-all}"
FILTER_NAME="${FILTER_NAME:-}"
FILTER_ARCH="${FILTER_ARCH:-}"
USE_GPU="${USE_GPU:-auto}"
CCACHE_MAXSIZE="${CCACHE_MAXSIZE:-200G}"
WORK_MOUNT_MODE="${WORK_MOUNT_MODE:-bind}"  # bind|tmpfs

# Host directories
HOST_WORK_DIR="${HOST_WORK_DIR:-/work/mplapack-work}"
HOST_CCACHE_DIR="${HOST_CCACHE_DIR:-/work/ccache}"

mkdir -p "$LOGDIR"
mkdir -p "$SUCCESS_DIR"
mkdir -p "$HOST_CCACHE_DIR"
if [[ "$WORK_MOUNT_MODE" == "bind" ]]; then
    mkdir -p "$HOST_WORK_DIR"
fi

if [[ "$WORK_MOUNT_MODE" != "bind" && "$WORK_MOUNT_MODE" != "tmpfs" ]]; then
    echo "ERROR: WORK_MOUNT_MODE must be 'bind' or 'tmpfs'." >&2
    exit 1
fi

# Initialize buildx builder if not exists
docker buildx inspect mplapack-builder >/dev/null 2>&1 || \
    docker buildx create --name mplapack-builder --use >/dev/null

log() {
    echo "[$(date +%H:%M:%S)] $*" >&2
}


# Image naming (safer than mplapack-*)
IMAGE_PREFIX="${IMAGE_PREFIX:-mplapack-qa}"   # e.g., mplapack-qa, mplapack-ci
IMAGE_TTL_HOURS="${IMAGE_TTL_HOURS:-}"        # optional: used only for label metadata

sanitize_token() {
    # Keep only [a-z0-9._-], convert others to "-". Lowercase output.
    local s="$1"
    s="$(echo "$s" | tr "[:upper:]" "[:lower:]" )"
    s="$(echo "$s" | sed -E 's/[^a-z0-9._-]+/-/g; s/^-+//; s/-+$//; s/-+/-/g')"
    echo "$s"
}

git_meta() {
    # Print: <branch> <sha7>. If not a git repo, print: nogit nogit
    local dir="$1"
    local br="nogit" sha="nogit"
    if command -v git >/dev/null 2>&1 && git -C "$dir" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        br="$(git -C "$dir" rev-parse --abbrev-ref HEAD 2>/dev/null || echo nogit)"
        sha="$(git -C "$dir" rev-parse --short=7 HEAD 2>/dev/null || echo nogit)"
    fi
    echo "$br $sha"
}

make_image_tag() {
    # Usage: make_image_tag <name> <arch_short> <source_type>
    local name_raw="$1" arch_short="$2" source_type="$3"
    local name arch src br sha stamp extra
    name="$(sanitize_token "$name_raw")"
    arch="$(sanitize_token "$arch_short")"
    src="$(sanitize_token "$source_type")"
    read -r br sha <<< "$(git_meta "$PROJECT_ROOT")"
    br="$(sanitize_token "$br")"
    sha="$(sanitize_token "$sha")"
    stamp="$(date +%Y%m%d_%H%M%S)"
    extra="${EXTRA_TAG:-}"
    if [[ -n "$extra" ]]; then extra="$(sanitize_token "$extra")"; fi

    # Format: <prefix>:<name>-<arch>-<src>-<branch>-<sha>-<timestamp>[-<extra>]
    local tag="${IMAGE_PREFIX}:${name}-${arch}-${src}-${br}-${sha}-${stamp}"
    if [[ -n "$extra" ]]; then tag="${tag}-${extra}"; fi
    echo "$tag"
}

make_success_stamp() {
    # Usage: make_success_stamp <name> <arch> <source_type>
    local name_raw="$1" arch_raw="$2" source_type_raw="$3"
    local name arch src
    name="$(sanitize_token "$name_raw")"
    arch="$(sanitize_token "$arch_raw")"
    src="$(sanitize_token "$source_type_raw")"
    echo "$SUCCESS_DIR/${name}_${arch}_${src}.ok"
}

abspath() {
    # Print an absolute path without relying on GNU-specific realpath/readlink options.
    local path="$1"
    if [[ -d "$path" ]]; then
        (
            cd "$path"
            pwd
        )
    else
        local dir base
        dir="$(dirname "$path")"
        base="$(basename "$path")"
        (cd "$dir" && printf '%s/%s\n' "$(pwd)" "$base")
    fi
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

# Return 0 if IMAGE has a manifest entry for PLATFORM (e.g., linux/mips64le)
has_manifest() {
    local image="$1" platform="$2" pat
    case "$platform" in
        linux/arm64)
            # buildx prints "linux/arm64" or "linux/arm64/v8"; both are arm64
            pat='linux/arm64(/v8)?'
            ;;
        linux/arm/v7)
            pat='linux/arm/v7'
            ;;
        linux/arm/v6)
            pat='linux/arm(/v6)?'
            ;;
        *)
            pat="$(echo "$platform" | sed 's|/|\\/|g')"
            ;;
    esac
    # "docker buildx imagetools inspect" outputs "Platform: linux/arm64/v8" lines.
    # Match the Platform field directly so trailing version suffixes don't cause misses.
    docker buildx imagetools inspect "$image" 2>/dev/null \
        | grep -qE "Platform:[[:space:]]+${pat}([[:space:]]|/|$)"
}

# Ordered list of candidate probe images for multi-arch manifest checks.
# These are well-known images that publish manifests for most platforms.
PROBE_CANDIDATES=(
    "debian:12"
    "ubuntu:24.04"
    "alpine:3.19"
    "busybox:1.36"
)

# Pick a probe image that actually supports the target platform (best-effort).
pick_probe_image() {
    local target="$1"
    local img
    for img in "${PROBE_CANDIDATES[@]}"; do
        if has_manifest "$img" "$target"; then
            echo "$img"
            return 0
        fi
    done
    # Fallback: nothing found
    return 1
}

# Ensure binfmt/qemu is installed for cross-arch runs
ensure_binfmt() {
    local target="$1"

    # Respect explicit disable
    if [[ "${EMULATE:-auto}" == "no" ]]; then
        return 1
    fi

    # Native run: no emulation needed
    if [[ -n "${HOST_PLATFORM:-}" && "$target" == "$HOST_PLATFORM" ]]; then
        return 0
    fi

    # Choose an appropriate probe image for this target
    local probe_img=""
    if probe_img="$(pick_probe_image "$target")"; then
        :
    else
        # We cannot even find a probe image that has a manifest for target.
        # This is NOT a binfmt problem.
        log "No probe image found with manifest for $target (probe image coverage / manifest-detection issue)."
        log "  (This is probe selection/manifest-detection failure, not a binfmt problem.)"
        log "  (Probe candidates: ${PROBE_CANDIDATES[*]})"
        # Emit a compact per-candidate verdict for debugging without flooding logs.
        for img in "${PROBE_CANDIDATES[@]}"; do
            if has_manifest "$img" "$target"; then
                log "  (probe check: OK  image=$img target=$target)"
            else
                log "  (probe check: NG  image=$img target=$target)"
            fi
        done
        if command -v docker >/dev/null 2>&1; then
            log "  (docker buildx version: $(docker buildx version 2>/dev/null | tr '\n' ' '))"
        fi
        log "  (Tip: arm64 may appear as linux/arm64/v8 in buildx output.)"
    fi

    # Quick sanity check: can we run a trivial container for target arch?
    if docker run --rm --platform "$target" "$probe_img" true >/dev/null 2>&1; then
        return 0
    fi

    # If not runnable, try enabling binfmt (if allowed)
    if [[ "${EMULATE:-auto}" == "auto" || "${EMULATE:-auto}" == "yes" ]]; then
        log "Enabling binfmt/qemu for cross-arch (${HOST_PLATFORM:-unknown} -> $target)"
        docker run --privileged --rm tonistiigi/binfmt:latest --install all >/dev/null 2>&1 || true

        # Re-test after install with the SAME probe image
        if docker run --rm --platform "$target" "$probe_img" true >/dev/null 2>&1; then
            return 0
        fi

        # Still failing: now it really is an emulation/runtime limitation.
        log "Cross-arch run still failing for $target even after binfmt install (probe=$probe_img)."
        return 1
    fi

    return 1
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
    local source_type=$5

    # Refuse editor backup Dockerfiles (very common footgun: Dockerfile.debian~)
    if [[ "$dockerfile" == *"~" ]]; then
        log "ERROR: Refusing backup Dockerfile: $DOCKER_DIR/$dockerfile"
        log "       Fix build-matrix.conf or remove the backup file."
        echo "$name,${arch##*/},$base,image,FAILED,0,$source_type"
        return 1
    fi
    if [[ ! -f "$DOCKER_DIR/$dockerfile" ]]; then
        log "ERROR: Dockerfile not found: $DOCKER_DIR/$dockerfile"
        echo "$name,${arch##*/},$base,image,FAILED,0,$source_type"
        return 1
    fi

    local arch_short=${arch##*/}
    local tag="$(make_image_tag "$name" "$arch_short" "$source_type")"
    local logprefix="$LOGDIR/${name}_${arch_short}_${source_type}"
    local logfile="${logprefix}_build.log"
    local stamp
    stamp="$(make_success_stamp "$name" "$arch" "$source_type")"
    local start=$(date +%s)
    local gpu_flag=""

    if [[ -L "$stamp" && ! -e "$stamp" ]]; then
        log "Removing broken success link: $stamp"
        rm -f "$stamp"
    fi

    if [[ -L "$stamp" && -e "$stamp" ]]; then
        log "  (Already passed; skipping $name / ${arch_short} / $source_type)"
        log "  (Success log: $(readlink "$stamp"))"
        echo "$name,$arch_short,$base,build,SKIPPED,0,$source_type"
        return 0
    fi

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

    # Prepare input directory.
    # If a tarball is used, mount it to /input (read-only) instead of copying
    # it into /work before the container starts.
    local input_mount=""

    # For tarball tests, create temp directory with tarball
    if [[ "$source_type" == "tarball" ]]; then
        if [[ -z "$TARBALL" || ! -f "$TARBALL" ]]; then
            echo "$name,$arch_short,$base,build,SKIPPED,0,$source_type"
            return 0
        fi
        local tmpdir="$LOGDIR/workdir_${name}_${arch_short}"
        mkdir -p "$tmpdir"
        cp "$TARBALL" "$tmpdir/"
        input_mount="$tmpdir"
    fi

    # Common docker run arguments.
    # - Default: bind-mount /work from the host so large release builds do not
    #   exhaust tmpfs during heavy MinGW linking.
    # - Optional: set WORK_MOUNT_MODE=tmpfs for small, fully memory-backed builds.
    # - For tarball mode, the tarball directory is mounted to /input:ro.
    local -a docker_run_args
    docker_run_args=(
        --rm
        --platform "$arch"
        -v "$HOST_CCACHE_DIR:/ccache:rw"
        -e CCACHE_DIR=/ccache
        -e CCACHE_MAXSIZE="$CCACHE_MAXSIZE"
    )
    if [[ "$WORK_MOUNT_MODE" == "bind" ]]; then
        docker_run_args+=( -v "$HOST_WORK_DIR:/work:rw" )
    else
        docker_run_args+=( --tmpfs /work:rw,exec )
    fi

    # Personality note: linux/386 on x86_64 host is handled by
    # entrypoint.sh inside the Dockerfile (auto-detects and wraps with linux32).
    # Do NOT use --entrypoint here; it would override the Dockerfile ENTRYPOINT.
    if [[ "$arch" == "linux/386" && "$HOST_PLATFORM" == "linux/amd64" ]]; then
        log "  (i386 on amd64: entrypoint.sh will apply linux32 personality)"
    fi

    if [[ -n "$gpu_flag" ]]; then
        # gpu_flag may contain multiple tokens (e.g., "--gpus all")
        docker_run_args+=( $gpu_flag )
    fi
    if [[ -n "$input_mount" ]]; then
        docker_run_args+=( -v "$input_mount:/input:ro" )
    fi

    # Build image (environment setup only, no compilation)
    if ! docker buildx build \
        --platform "$arch" \
        --build-arg BASE="$base" \
        --label org.mplapack.project=mplapack \
        --label org.mplapack.purpose="${IMAGE_PREFIX}" \
        --label org.mplapack.matrix_name="$name" \
        --label org.mplapack.arch="$arch" \
        --label org.mplapack.source_type="$source_type" \
        --load \
        -t "$tag" \
        -f "$DOCKER_DIR/$dockerfile" \
        "$DOCKER_DIR" \
        > "${logprefix}_image.log" 2>&1; then

        local elapsed=$(($(date +%s) - start))
        echo "$name,$arch_short,$base,image,FAILED,$elapsed,$source_type"
        return 1
    fi

    # Run the branch/tarball build in a single container.
    # /work is tmpfs so a second docker run would lose all build artifacts.
    if ! docker run "${docker_run_args[@]}" \
        "$tag" > "$logfile" 2>&1; then

        local elapsed=$(($(date +%s) - start))
        echo "$name,$arch_short,$base,build,FAILED,$elapsed,$source_type"
        return 1
    fi

    local elapsed=$(($(date +%s) - start))
    ln -sfn "$(abspath "$logfile")" "$stamp"
    echo "$name,$arch_short,$base,build,OK,$elapsed,$source_type"
}

# Process build matrix
run_matrix() {
    echo "name,arch,base,stage,result,elapsed,source_type"

    while IFS= read -r line; do
        # Skip comments and empty lines
        [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue

        IFS='|' read -r -a fields <<< "$line"
        local name="${fields[0]:-}"
        local base="${fields[1]:-}"
        local archs="${fields[2]:-}"
        local dockerfile="${fields[3]:-}"
        local source_type=""

        if [[ "${fields[5]:-}" == remote-* ]]; then
            source_type="${fields[5]}"
        else
            source_type="${fields[4]:-branch}"
        fi

        # Apply name filter
        [[ -n "$FILTER_NAME" && "$name" != *"$FILTER_NAME"* ]] && continue

        # Set defaults
        source_type="${source_type:-branch}"
        # Remote targets are launched by explicit Makefile targets, not Docker matrix runs.
        [[ "$source_type" == remote-* ]] && continue

        # Apply phase filter
        [[ "$PHASE" != "all" && "$source_type" != "$PHASE" ]] && continue

        # Process each architecture
        IFS=',' read -ra arch_array <<< "$archs"
        for arch in "${arch_array[@]}"; do
            local arch_short=${arch##*/}

            # Apply arch filter
            [[ -n "$FILTER_ARCH" && "$arch_short" != "$FILTER_ARCH" ]] && continue

            log "START $name / $arch_short ($source_type)"
            build_one "$name" "$base" "$arch" "$dockerfile" "$source_type" || true
        done
    done < "$CONF_FILE"
}

# Display summary of results
show_summary() {
    [[ -f "$LOGDIR/results.csv" ]] || return

    echo "" >&2
    log "=== SUMMARY ==="
    column -t -s, "$LOGDIR/results.csv" >&2

    local total=$(( $(wc -l < "$LOGDIR/results.csv") - 1 ))
    local passed; passed=$(grep -c ',OK,' "$LOGDIR/results.csv") || passed=0
    local failed; failed=$(grep -c ',FAILED,' "$LOGDIR/results.csv") || failed=0
    local skipped; skipped=$(grep -c ',SKIPPED,' "$LOGDIR/results.csv") || skipped=0

    echo "" >&2
    echo "Total: $total | Passed: $passed | Failed: $failed | Skipped: $skipped" >&2
    echo "Logs:   $LOGDIR/" >&2
    echo "Passes: $SUCCESS_DIR/" >&2
    echo "ccache: $HOST_CCACHE_DIR" >&2

    [[ $failed -eq 0 ]] || exit 1
}

# Main entry point
case "${1:-}" in
    clean-success)
        rm -rf "$SUCCESS_DIR"
        mkdir -p "$SUCCESS_DIR"
        log "Cleaned success links: $SUCCESS_DIR"
        ;;
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
        log "=== Phase 1: Branch matrix builds ==="
        PHASE=branch run_matrix | tee "$LOGDIR/results_branch.csv"

        echo "" >&2
        make_dist

        echo "" >&2
        log "=== Phase 2: Tarball tests ==="
        PHASE=tarball run_matrix | tee "$LOGDIR/results_tarball.csv"

        # Merge results
        cat "$LOGDIR/results_branch.csv" > "$LOGDIR/results.csv"
        tail -n +2 "$LOGDIR/results_tarball.csv" >> "$LOGDIR/results.csv"

        show_summary
        ;;
esac
