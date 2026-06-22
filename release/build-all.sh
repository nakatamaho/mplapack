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
REMOTE_TARBALL_SSH="${REMOTE_TARBALL_SSH:-${REMOTE_LINUX_SSH:-ssh}}"
REMOTE_TARBALL_CCACHE_DIR="${REMOTE_TARBALL_CCACHE_DIR:-}"
REMOTE_TARBALL_CCACHE_MAXSIZE="${REMOTE_TARBALL_CCACHE_MAXSIZE:-}"
MPLAPACK_REF="${MPLAPACK_REF:-$(git -C "$PROJECT_ROOT" rev-parse HEAD 2>/dev/null || echo master)}"

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

lang_c_date() {
    LANG=C LC_ALL=C date
}

format_elapsed() {
    local seconds="$1"
    printf '%02d:%02d:%02d' "$((seconds / 3600))" "$(((seconds % 3600) / 60))" "$((seconds % 60))"
}

write_log_start() {
    local logfile="$1"
    local label="$2"
    local epoch="$3"
    {
        echo "=== LOG START: ${label} ==="
        echo "LOG_START_EPOCH=${epoch}"
        printf 'LOG_START_DATE='
        lang_c_date
        echo
    } >> "$logfile"
}

write_log_end() {
    local logfile="$1"
    local label="$2"
    local rc="$3"
    local end_epoch="$4"
    local elapsed="$5"
    {
        echo
        echo "=== LOG END: ${label} ==="
        echo "LOG_END_EPOCH=${end_epoch}"
        printf 'LOG_END_DATE='
        lang_c_date
        echo "LOG_ELAPSED_SECONDS=${elapsed}"
        echo "LOG_ELAPSED_HMS=$(format_elapsed "$elapsed")"
        echo "=== LOG ELAPSED: ${elapsed}s ($(format_elapsed "$elapsed")) | rc: ${rc} ==="
    } >> "$logfile"
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
    # Print: <ref> <sha7>. If not a git repo, print: nogit nogit
    local dir="$1"
    local ref="${MPLAPACK_REF:-nogit}" sha="nogit"
    if command -v git >/dev/null 2>&1 && git -C "$dir" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        if git -C "$dir" rev-parse --verify -q "$ref^{commit}" >/dev/null 2>&1; then
            sha="$(git -C "$dir" rev-parse --short=7 "$ref^{commit}" 2>/dev/null || echo nogit)"
        fi
    fi
    echo "$ref $sha"
}

make_image_tag() {
    # Usage: make_image_tag <name> <arch_short> <source_type>
    local name_raw="$1" arch_short="$2" source_type="$3"
    local name arch src ref sha stamp extra
    name="$(sanitize_token "$name_raw")"
    arch="$(sanitize_token "$arch_short")"
    src="$(sanitize_token "$source_type")"
    read -r ref sha <<< "$(git_meta "$PROJECT_ROOT")"
    ref="$(sanitize_token "$ref")"
    sha="$(sanitize_token "$sha")"
    stamp="$(date +%Y%m%d_%H%M%S)"
    extra="${EXTRA_TAG:-}"
    if [[ -n "$extra" ]]; then extra="$(sanitize_token "$extra")"; fi

    # Format: <prefix>:<name>-<arch>-<src>-<ref>-<sha>-<timestamp>[-<extra>]
    local tag="${IMAGE_PREFIX}:${name}-${arch}-${src}-${ref}-${sha}-${stamp}"
    if [[ -n "$extra" ]]; then tag="${tag}-${extra}"; fi
    echo "$tag"
}

make_success_stamp() {
    # Usage: make_success_stamp <name> <arch> <source_type>
    local name_raw="$1" arch_raw="$2" source_type_raw="$3"
    local name arch src ref sha ref_part
    name="$(sanitize_token "$name_raw")"
    arch="$(sanitize_token "$arch_raw")"
    src="$(sanitize_token "$source_type_raw")"
    if [[ "$source_type_raw" == "ref" ]]; then
        read -r ref sha <<< "$(git_meta "$PROJECT_ROOT")"
        ref="$(sanitize_token "$ref")"
        sha="$(sanitize_token "$sha")"
        ref_part="_${ref}_${sha}"
    else
        ref_part=""
    fi
    echo "$SUCCESS_DIR/${name}_${arch}_${src}${ref_part}.ok"
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


write_checksum_file() {
    # Usage: write_checksum_file <sha256|md5> <tarball> <output-file>
    local algorithm="$1" tarball="$2" output="$3"
    local dir base
    dir="$(dirname "$tarball")"
    base="$(basename "$tarball")"

    (
        cd "$dir"
        case "$algorithm" in
            sha256)
                if command -v sha256sum >/dev/null 2>&1; then
                    sha256sum "$base"
                elif command -v shasum >/dev/null 2>&1; then
                    shasum -a 256 "$base"
                elif command -v openssl >/dev/null 2>&1; then
                    openssl dgst -sha256 -r "$base"
                else
                    echo "ERROR: no sha256 checksum command found" >&2
                    exit 1
                fi
                ;;
            md5)
                if command -v md5sum >/dev/null 2>&1; then
                    md5sum "$base"
                elif command -v md5 >/dev/null 2>&1; then
                    md5 -r "$base"
                elif command -v openssl >/dev/null 2>&1; then
                    openssl dgst -md5 -r "$base"
                else
                    echo "ERROR: no md5 checksum command found" >&2
                    exit 1
                fi
                ;;
            *)
                echo "ERROR: unknown checksum algorithm: $algorithm" >&2
                exit 1
                ;;
        esac
    ) > "$output"
}

report_tarball() {
    # Create checksum sidecars next to the tarball and print the release artifact summary.
    local tarball="$1"
    local tarball_abs sha_file md5_file sha md5

    if [[ -z "$tarball" || ! -f "$tarball" ]]; then
        log "ERROR: tarball not found: ${tarball:-<empty>}"
        exit 1
    fi

    tarball_abs="$(abspath "$tarball")"
    sha_file="${tarball_abs}.sha256sum"
    md5_file="${tarball_abs}.md5sum"

    write_checksum_file sha256 "$tarball_abs" "$sha_file"
    write_checksum_file md5 "$tarball_abs" "$md5_file"

    read -r sha _ < "$sha_file"
    read -r md5 _ < "$md5_file"

    log "Tarball: $tarball_abs"
    log "SHA256: $sha"
    log "SHA256 file: $sha_file"
    log "MD5: $md5"
    log "MD5 file: $md5_file"
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

    log "Generating compare Makefile.am files..."
    (
        cd mplapack/test/compare
        bash gen.Makefile.am.sh
    ) > "$distdir/gen_compare_makefiles.log" 2>&1

    log "Running autoreconf..."
    autoreconf -fi > "$distdir/autoreconf.log" 2>&1

    log "Running configure..."
    ./configure > "$distdir/configure.log" 2>&1

    log "Running make dist..."
    make dist > "$distdir/make_dist.log" 2>&1

    TARBALL=$(ls -t mplapack-*.tar.xz | head -1)
    TARBALL=$(abspath "$TARBALL")
    log "Created: $TARBALL"
    report_tarball "$TARBALL"

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
    local image_logfile="${logprefix}_image.log"
    local stamp
    stamp="$(make_success_stamp "$name" "$arch" "$source_type")"
    local start
    local gpu_flag=""
    start=$(date +%s)

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
        -e MPLAPACK_REF="$MPLAPACK_REF"
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
    local image_start image_end image_elapsed
    image_start="$(date +%s)"
    : > "$image_logfile"
    write_log_start "$image_logfile" "docker image ${name} / ${arch_short} (${source_type})" "$image_start"
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
        >> "$image_logfile" 2>&1; then

        image_end="$(date +%s)"
        image_elapsed=$((image_end - image_start))
        write_log_end "$image_logfile" "docker image ${name} / ${arch_short} (${source_type})" 1 "$image_end" "$image_elapsed"
        local elapsed=$(($(date +%s) - start))
        echo "$name,$arch_short,$base,image,FAILED,$elapsed,$source_type"
        return 1
    fi
    image_end="$(date +%s)"
    image_elapsed=$((image_end - image_start))
    write_log_end "$image_logfile" "docker image ${name} / ${arch_short} (${source_type})" 0 "$image_end" "$image_elapsed"

    # Run the ref/tarball build in a single container.
    # /work is tmpfs so a second docker run would lose all build artifacts.
    local run_start run_end run_elapsed
    run_start="$(date +%s)"
    : > "$logfile"
    write_log_start "$logfile" "docker run ${name} / ${arch_short} (${source_type})" "$run_start"
    if ! docker run "${docker_run_args[@]}" \
        "$tag" >> "$logfile" 2>&1; then

        run_end="$(date +%s)"
        run_elapsed=$((run_end - run_start))
        write_log_end "$logfile" "docker run ${name} / ${arch_short} (${source_type})" 1 "$run_end" "$run_elapsed"
        local elapsed=$(($(date +%s) - start))
        echo "$name,$arch_short,$base,build,FAILED,$elapsed,$source_type"
        return 1
    fi

    run_end="$(date +%s)"
    run_elapsed=$((run_end - run_start))
    write_log_end "$logfile" "docker run ${name} / ${arch_short} (${source_type})" 0 "$run_end" "$run_elapsed"
    local elapsed=$(($(date +%s) - start))
    ln -sfn "$(abspath "$logfile")" "$stamp"
    echo "$name,$arch_short,$base,build,OK,$elapsed,$source_type"
}


shell_quote() {
    printf '%q' "$1"
}

remote_tarball_default_ccache_dir() {
    local target_dir="$1" rest user
    case "$target_dir" in
        /Users/*/tmp/*)
            rest="${target_dir#/Users/}"
            user="${rest%%/*}"
            echo "/Users/${user}/.ccache"
            ;;
        /home/*/tmp/*)
            rest="${target_dir#/home/}"
            user="${rest%%/*}"
            echo "/home/${user}/.ccache"
            ;;
        *)
            echo "${HOME:-/tmp}/.ccache"
            ;;
    esac
}

remote_tarball_default_ccache_maxsize() {
    local target_dir="$1"
    case "$target_dir" in
        /home/*/tmp/*) echo '200G' ;;
        *) echo '80G' ;;
    esac
}

# Build and test a tarball on a native remote Docker host.
build_remote_tarball_one() {
    local name=$1
    local host=$2
    local target_dir=$3
    local dockerfile=$4
    local remote_cmd=$5
    local source_type=$6
    local base=$7
    local arch=$8
    local matrix_ccache_dir=${9:-}
    local matrix_ccache_maxsize=${10:-}

    local remote_ccache_dir="${REMOTE_TARBALL_CCACHE_DIR:-${matrix_ccache_dir:-$(remote_tarball_default_ccache_dir "$target_dir")}}"
    local remote_ccache_maxsize="${REMOTE_TARBALL_CCACHE_MAXSIZE:-${matrix_ccache_maxsize:-$(remote_tarball_default_ccache_maxsize "$target_dir")}}"

    local arch_short=${arch##*/}
    local logprefix="$LOGDIR/${name}_${arch_short}_${source_type}"
    local logfile="${logprefix}_build.log"
    local stamp
    stamp="$(make_success_stamp "$name" "$arch" "$source_type")"
    local start elapsed rc
    local image_tag="$(make_image_tag "$name" "$arch_short" "$source_type")"
    image_tag="${image_tag/:/-remote:}"

    if [[ -z "$TARBALL" || ! -f "$TARBALL" ]]; then
        echo "$name,$arch_short,$base,build,SKIPPED,0,$source_type"
        return 0
    fi
    if [[ ! -f "$DOCKER_DIR/$dockerfile" ]]; then
        log "ERROR: Dockerfile not found: $DOCKER_DIR/$dockerfile"
        echo "$name,$arch_short,$base,image,FAILED,0,$source_type"
        return 1
    fi

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

    local tarball_name
    tarball_name="$(basename "$TARBALL")"
    local remote_tarball="$target_dir/input/$tarball_name"
    local remote_script="$SCRIPT_DIR/run-remote-tarball-docker-worker.sh"
    if [[ ! -f "$remote_script" ]]; then
        log "ERROR: remote tarball worker not found: $remote_script"
        echo "$name,$arch_short,$base,build,FAILED,0,$source_type"
        return 1
    fi

    start=$(date +%s)
    log "START remote $name / $arch_short ($source_type) on $host"
    log "Remote tarball log: $logfile"
    : > "$logfile"
    write_log_start "$logfile" "remote tarball ${name} / ${arch_short} (${source_type})" "$start"

    set +e
    {
        "$REMOTE_TARBALL_SSH" "$host" \
            "MPLAPACK_REMOTE_WORKDIR=$(shell_quote "$target_dir") bash -s" <<'REMOTE_PREP'
set -euo pipefail
export PATH="/opt/local/bin:/opt/local/sbin:${PATH}"
: "${MPLAPACK_REMOTE_WORKDIR:?}"
if [ -z "${HOME:-}" ] || [ "${HOME}" = "/" ]; then
    echo "ERROR: invalid HOME: ${HOME:-<unset>}"
    exit 1
fi
case "${MPLAPACK_REMOTE_WORKDIR}" in
    "${HOME}/"*) ;;
    *)
        echo "ERROR: refusing remote workdir outside HOME: ${MPLAPACK_REMOTE_WORKDIR}"
        exit 1
        ;;
esac
rm -rf "${MPLAPACK_REMOTE_WORKDIR}/context" "${MPLAPACK_REMOTE_WORKDIR}/input" "${MPLAPACK_REMOTE_WORKDIR}/work"
mkdir -p "${MPLAPACK_REMOTE_WORKDIR}/context" "${MPLAPACK_REMOTE_WORKDIR}/input" "${MPLAPACK_REMOTE_WORKDIR}/work"
REMOTE_PREP
        tar -C "$PROJECT_ROOT" -czf - release/docker | \
            "$REMOTE_TARBALL_SSH" "$host" "tar -C $(shell_quote "$target_dir/context") -xzf -"
        "$REMOTE_TARBALL_SSH" "$host" "cat > $(shell_quote "$remote_tarball")" < "$TARBALL"
        "$REMOTE_TARBALL_SSH" "$host" \
            "MPLAPACK_REMOTE_WORKDIR=$(shell_quote "$target_dir") MPLAPACK_DOCKER_BASE=$(shell_quote "$base") MPLAPACK_DOCKERFILE=$(shell_quote "$dockerfile") MPLAPACK_IMAGE_TAG=$(shell_quote "$image_tag") MPLAPACK_TARBALL_NAME=$(shell_quote "$tarball_name") MPLAPACK_CCACHE_DIR=$(shell_quote "$remote_ccache_dir") MPLAPACK_CCACHE_MAXSIZE=$(shell_quote "$remote_ccache_maxsize") MPLAPACK_DOCKER_PLATFORM=$(shell_quote "$arch") $remote_cmd -s" < "$remote_script"
    } >> "$logfile" 2>&1
    rc=$?
    set -e

    end="$(date +%s)"
    elapsed=$((end - start))
    write_log_end "$logfile" "remote tarball ${name} / ${arch_short} (${source_type})" "$rc" "$end" "$elapsed"
    if [[ "$rc" -eq 0 ]]; then
        ln -sfn "$(abspath "$logfile")" "$stamp"
        echo "$name,$arch_short,$base,build,OK,$elapsed,$source_type"
        return 0
    fi

    echo "$name,$arch_short,$base,build,FAILED,$elapsed,$source_type"
    return "$rc"
}

# Process build matrix
run_matrix() {
    echo "name,arch,base,stage,result,elapsed,source_type"

    local remote_tmpdir
    local -a remote_pids=()
    local -a remote_result_files=()
    local -a remote_fallback_rows=()
    remote_tmpdir="$(mktemp -d "$LOGDIR/remote_tarball_results.XXXXXX")"

    cleanup_remote_tarball_jobs() {
        local rc="${1:-130}" pid
        trap - INT TERM HUP
        for pid in "${remote_pids[@]:-}"; do
            if kill -0 "$pid" 2>/dev/null; then
                kill "$pid" 2>/dev/null || true
            fi
        done
        for pid in "${remote_pids[@]:-}"; do
            wait "$pid" 2>/dev/null || true
        done
        rm -rf "$remote_tmpdir"
        exit "$rc"
    }

    trap 'cleanup_remote_tarball_jobs 130' INT
    trap 'cleanup_remote_tarball_jobs 143' TERM HUP

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
            source_type="${fields[4]:-ref}"
        fi

        # Apply name filter
        [[ -n "$FILTER_NAME" && "$name" != *"$FILTER_NAME"* ]] && continue

        # Set defaults
        source_type="${source_type:-ref}"

        if [[ "$source_type" == "remote-tarball-docker" ]]; then
            local host="$base"
            local target_dir="$archs"
            local remote_cmd="${fields[4]:-bash}"
            local docker_base="${fields[6]:-}"
            local remote_arch="${fields[7]:-linux/arm64}"
            local remote_ccache_dir="${fields[8]:-}"
            local remote_ccache_maxsize="${fields[9]:-}"
            local arch_short="${remote_arch##*/}"
            local result_file

            [[ "$PHASE" != "all" && "$PHASE" != "tarball" ]] && continue
            [[ -n "$FILTER_ARCH" && "$arch_short" != "$FILTER_ARCH" ]] && continue

            log "START $name / $arch_short ($source_type) on $host"
            result_file="$remote_tmpdir/${#remote_pids[@]}_${name}_${arch_short}.csv"
            (
                build_remote_tarball_one "$name" "$host" "$target_dir" "$dockerfile" "$remote_cmd" "$source_type" "$docker_base" "$remote_arch" "$remote_ccache_dir" "$remote_ccache_maxsize"
            ) > "$result_file" &
            remote_pids+=("$!")
            remote_result_files+=("$result_file")
            remote_fallback_rows+=("$name,$arch_short,$docker_base,build,FAILED,0,$source_type")
            continue
        fi

        # Remote distcheck targets are launched by explicit Makefile targets, not Docker matrix runs.
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

    local i
    for i in "${!remote_pids[@]}"; do
        wait "${remote_pids[$i]}" || true
    done

    for i in "${!remote_result_files[@]}"; do
        if [[ -s "${remote_result_files[$i]}" ]]; then
            cat "${remote_result_files[$i]}"
        else
            echo "${remote_fallback_rows[$i]}"
        fi
    done

    rm -rf "$remote_tmpdir"
    trap - INT TERM HUP
    return 0
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
    ref)
        PHASE=ref run_matrix | tee "$LOGDIR/results.csv"
        show_summary
        ;;
    tarball)
        if [[ -z "$TARBALL" ]]; then
            make_dist
        else
            report_tarball "$TARBALL"
        fi
        PHASE=tarball run_matrix | tee "$LOGDIR/results.csv"
        show_summary
        ;;
    *)
        # Full cycle: ref -> dist -> tarball
        log "=== Phase 1: Git ref matrix builds ==="
        PHASE=ref run_matrix | tee "$LOGDIR/results_ref.csv"

        echo "" >&2
        make_dist

        echo "" >&2
        log "=== Phase 2: Tarball tests ==="
        PHASE=tarball run_matrix | tee "$LOGDIR/results_tarball.csv"

        # Merge results
        cat "$LOGDIR/results_ref.csv" > "$LOGDIR/results.csv"
        tail -n +2 "$LOGDIR/results_tarball.csv" >> "$LOGDIR/results.csv"

        show_summary
        ;;
esac
