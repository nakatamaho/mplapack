#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}"
CONF_FILE="${CONF_FILE:-$SCRIPT_DIR/build-matrix.conf}"
LOGDIR="${LOGDIR:-$SCRIPT_DIR/logs/$(date +%Y%m%d_%H%M%S)}"
SUCCESS_DIR="${SUCCESS_DIR:-$SCRIPT_DIR/success}"
FAILED_DIR="${FAILED_DIR:-$SCRIPT_DIR/failed}"
REMOTE_SSH="${REMOTE_LINUX_SSH:-ssh}"
MPLAPACK_REF="${MPLAPACK_REF:-$(git -C "$PROJECT_ROOT" rev-parse HEAD)}"
MPLAPACK_SOURCE_MODE="${MPLAPACK_SOURCE_MODE:-dist}"

name="${1:?Usage: run-remote-linux-docker.sh <matrix-name>}"

mkdir -p "$LOGDIR" "$SUCCESS_DIR" "$FAILED_DIR"

row="$(awk -F'|' -v target="$name" '$1 == target { print; exit }' "$CONF_FILE")"
if [[ -z "$row" ]]; then
    echo "ERROR: $name not found in $CONF_FILE" >&2
    exit 1
fi

IFS='|' read -r matrix_name host target_dir script_rel remote_cmd source_type matrix_docker_base <<< "$row"
remote_cmd="${remote_cmd:-bash}"
matrix_docker_base="${matrix_docker_base:-}"
remote_docker_base="${MPLAPACK_DOCKER_BASE:-$matrix_docker_base}"
remote_distro_version="${MPLAPACK_DISTRO_VERSION:-}"
if [[ -z "$remote_distro_version" && "$remote_docker_base" == *:* ]]; then
    remote_distro_version="${remote_docker_base##*:}"
fi

if [[ "$source_type" != "remote-linux-docker" ]]; then
    echo "ERROR: $name is source_type=$source_type, expected remote-linux-docker" >&2
    exit 1
fi
if [[ -z "$host" || -z "$target_dir" || -z "$script_rel" ]]; then
    echo "ERROR: invalid remote-linux-docker row: $row" >&2
    exit 1
fi

script_path="$PROJECT_ROOT/$script_rel"
if [[ ! -f "$script_path" ]]; then
    echo "ERROR: build script not found: $script_path" >&2
    exit 1
fi

tier="${matrix_name%%-*}"
arch="${matrix_name##*-}"
docker_label="${remote_docker_base:-ubuntu:${remote_distro_version:-default}}"
os_label="$(printf '%s' "$docker_label" | awk -F: '{name=$1; version=$2; sub(/^.*\//, "", name); gsub(/[^A-Za-z0-9]+/, "_", name); gsub(/[^A-Za-z0-9]+/, "_", version); print name version}')"
if [[ -z "$os_label" ]]; then
    os_label="linux_unknown"
fi
if [[ "$matrix_name" == *mingw* ]]; then
    os_label="windows"
fi
name_middle="${matrix_name#${tier}-}"
name_middle="${name_middle%-${arch}}"
flavor=""
if [[ "$os_label" == "windows" ]]; then
    flavor="$name_middle"
elif [[ "$name_middle" == *-* ]]; then
    flavor="${name_middle#*-}"
fi
name_prefix="${tier}-${os_label}${flavor:+-${flavor}}-${arch}"
stamp_prefix="$SUCCESS_DIR/${name_prefix}-"
stamp_suffix="-linux-docker.ok"
logfile="$LOGDIR/${name_prefix}-linux-docker.log"
source_log_part=""
resultfile="$LOGDIR/results_${name_prefix}.csv"
stamp=""
COLLECTED_RESULTS_STAGE=""

cleanup_stale_success_links() {
    find "$SUCCESS_DIR" -maxdepth 1 -xtype l -name '*.ok' -exec rm -f {} +
    find "$FAILED_DIR" -maxdepth 1 -xtype l -name '*.failed' -exec rm -f {} +
}

failed_stamp_path() {
    local label="${1:-failed}"
    printf '%s/%s-%s%s-linux-docker.failed\n' "$FAILED_DIR" "$name_prefix" "$label" "$source_log_part"
}

link_failed_log() {
    local label="${1:-failed}"
    mkdir -p "$FAILED_DIR"
    ln -sfn "$logfile" "$(failed_stamp_path "$label")"
}

clear_failed_logs() {
    find "$FAILED_DIR" -maxdepth 1 -type l -name "${name_prefix}-*${source_log_part}-linux-docker.failed" -exec rm -f {} +
}

detect_compiler_label_from_results() {
    local stage="$1"
    local outdir base compiler

    if [[ -d "$stage" ]]; then
        outdir="$(find "$stage" -path '*/results/*' -type d -mindepth 4 -maxdepth 5 | head -n 1 || true)"
        if [[ -n "$outdir" ]]; then
            base="$(basename "$outdir")"
            compiler="$(printf '%s' "$base" | sed -n 's/^.*_\([A-Za-z][A-Za-z0-9+.-]*-[0-9][A-Za-z0-9_.-]*\)$/\1/p')"
            if [[ -n "$compiler" ]]; then
                printf '%s' "$compiler" | sed 's/-//; s/[^A-Za-z0-9_+.-]/_/g'
                return 0
            fi
        fi
    fi
    printf 'compiler_unknown'
}

copy_tree_contents() {
    local src="$1"
    local dst="$2"
    mkdir -p "$dst"
    if command -v rsync >/dev/null 2>&1; then
        rsync -a "$src/" "$dst/"
    else
        cp -a "$src/." "$dst/"
    fi
}

collect_remote_test_results() {
    local remote_root="${target_dir}.distcheck-results"
    local local_stage="$LOGDIR/${matrix_name}_distcheck-results"
    local version_dir version suite src dst copied
    COLLECTED_RESULTS_STAGE="$local_stage"

    if ! "$REMOTE_SSH" "$host" "test -d '$remote_root'" >/dev/null 2>&1; then
        echo "No remote distcheck results found: $host:$remote_root" >&2
        return 0
    fi

    rm -rf "$local_stage"
    mkdir -p "$local_stage"
    if command -v rsync >/dev/null 2>&1 && rsync -a "$host:$remote_root/" "$local_stage/"; then
        :
    else
        "$REMOTE_SSH" "$host" "cd '$remote_root' && tar -cf - ." | tar -C "$local_stage" -xf -
    fi

    copied=0
    for version_dir in "$local_stage"/*; do
        [[ -d "$version_dir" ]] || continue
        version="$(basename "$version_dir")"
        for suite in lin eig; do
            src="$version_dir/$suite/results"
            [[ -d "$src" ]] || continue
            dst="$PROJECT_ROOT/mplapack/test/$suite/results/$version"
            copy_tree_contents "$src" "$dst"
            copied=$((copied + 1))
            echo "Collected remote $suite results: $host:$remote_root/$version/$suite/results -> $dst" >&2
        done
    done

    if [[ "$copied" -eq 0 ]]; then
        echo "No staged lin/eig results found under $host:$remote_root" >&2
    fi
}

lang_c_date() {
    LANG=C LC_ALL=C date
}

format_elapsed() {
    local seconds="$1"
    printf '%02d:%02d:%02d' "$((seconds / 3600))" "$(((seconds % 3600) / 60))" "$((seconds % 60))"
}


append_final_ccache_stats() {
    local src_logfile="$1"
    local tmp

    [[ -f "$src_logfile" ]] || return 0
    tmp="$(mktemp "${src_logfile}.ccache-summary.XXXXXX")" || return 0
    awk '
        /=== CCACHE STATS \(START\) ===/ { in_start=1; start=$0 ORS; next }
        in_start {
            start=start $0 ORS
            if ($0 ~ /=== END CCACHE STATS \(START\) ===/) in_start=0
            next
        }
        /=== CCACHE STATS \(END\) ===/ { in_end=1; end=$0 ORS; next }
        in_end {
            end=end $0 ORS
            if ($0 ~ /=== END CCACHE STATS \(END\) ===/) in_end=0
            next
        }
        END {
            if (start != "" || end != "") {
                print ""
                print "=== FINAL CCACHE STATS SUMMARY ==="
                if (start != "") printf "%s", start
                if (end != "") printf "%s", end
                print "=== END FINAL CCACHE STATS SUMMARY ==="
            }
        }
    ' "$src_logfile" > "$tmp"
    if [[ -s "$tmp" ]]; then
        cat "$tmp" >> "$src_logfile"
    fi
    rm -f "$tmp"
}

write_log_start() {
    local epoch="$1"
    {
        echo "=== LOG START: ${matrix_name} ==="
        echo "LOG_START_EPOCH=${epoch}"
        printf 'LOG_START_DATE='
        lang_c_date
        echo
    } >> "$logfile"
}

write_log_end() {
    local rc="$1"
    local end_epoch="$2"
    local elapsed="$3"
    {
        echo
        echo "=== LOG END: ${matrix_name} ==="
        echo "LOG_END_EPOCH=${end_epoch}"
        printf 'LOG_END_DATE='
        lang_c_date
        echo "LOG_ELAPSED_SECONDS=${elapsed}"
        echo "LOG_ELAPSED_HMS=$(format_elapsed "$elapsed")"
        echo "=== LOG ELAPSED: ${elapsed}s ($(format_elapsed "$elapsed")) | rc: ${rc} ==="
    } >> "$logfile"
    append_final_ccache_stats "$logfile"
}

mkdir -p "$(dirname "$resultfile")"
echo "name,arch,base,stage,result,elapsed,source_type" > "$resultfile"
if [[ "$MPLAPACK_SOURCE_MODE" == "dist" && -z "${MPLAPACK_SOURCE_TARBALL:-}" ]]; then
    source_info_file="$LOGDIR/source/source.env"
    PROJECT_ROOT="$PROJECT_ROOT" LOGDIR="$LOGDIR" MPLAPACK_SOURCE_INFO_FILE="$source_info_file" \
        "$SCRIPT_DIR/make-source-snapshot.sh"
    # shellcheck disable=SC1090
    source "$source_info_file"
elif [[ "$MPLAPACK_SOURCE_MODE" != "dist" && "$MPLAPACK_SOURCE_MODE" != "ref" ]]; then
    echo "ERROR: MPLAPACK_SOURCE_MODE must be 'dist' or 'ref' (got '$MPLAPACK_SOURCE_MODE')" >&2
    exit 1
fi

if [[ "$MPLAPACK_SOURCE_MODE" == "dist" ]]; then
    source_log_part="-${MPLAPACK_SOURCE_LABEL:-dist}"
    stamp_suffix="${source_log_part}-linux-docker.ok"
    logfile="$LOGDIR/${name_prefix}${source_log_part}-linux-docker.log"
fi

cleanup_stale_success_links

existing_stamp="$(find "$SUCCESS_DIR" -maxdepth 1 -type l -name "$(basename "$stamp_prefix")*${stamp_suffix}" | head -n 1 || true)"
if [[ -n "$existing_stamp" && ! -e "$existing_stamp" ]]; then
    rm -f "$existing_stamp"
    existing_stamp=""
fi
if [[ -n "$existing_stamp" && -e "$existing_stamp" ]]; then
    echo "Already passed; skipping $matrix_name" >&2
    echo "Success log: $(readlink "$existing_stamp")" >&2
    echo "$matrix_name,$arch,${remote_docker_base:-$host},test,SKIPPED,0,$source_type" | tee -a "$resultfile"
    exit 0
fi


echo "Running $script_rel on $host:$target_dir" >&2
echo "MPLAPACK_REF: $MPLAPACK_REF" >&2
echo "MPLAPACK_SOURCE_MODE: $MPLAPACK_SOURCE_MODE" >&2
if [[ "$MPLAPACK_SOURCE_MODE" == "dist" ]]; then
    echo "MPLAPACK_SOURCE_TARBALL: ${MPLAPACK_SOURCE_TARBALL:-<unset>}" >&2
    echo "MPLAPACK_SOURCE_LABEL: ${MPLAPACK_SOURCE_LABEL:-<unset>}" >&2
fi
echo "MPLAPACK_DISTRO_VERSION: ${remote_distro_version:-<default>}" >&2
echo "MPLAPACK_DOCKER_BASE: ${remote_docker_base:-<default>}" >&2
echo "Log: $logfile" >&2

context_name="${name_prefix}${source_log_part}"
context_tar="$LOGDIR/${context_name}_context.tar.gz"
source_bundle="$LOGDIR/${context_name}_source.bundle"
remote_context_tar="$(dirname "$target_dir")/${context_name}.context.tar.gz"
if [[ "$MPLAPACK_SOURCE_MODE" == "dist" ]]; then
    if [[ -z "${MPLAPACK_SOURCE_TARBALL:-}" || ! -f "$MPLAPACK_SOURCE_TARBALL" ]]; then
        echo "ERROR: MPLAPACK_SOURCE_TARBALL is not set or missing: ${MPLAPACK_SOURCE_TARBALL:-<unset>}" >&2
        exit 1
    fi
    tar_args=(-C "$PROJECT_ROOT" release/docker -C "$(dirname "$MPLAPACK_SOURCE_TARBALL")" "$(basename "$MPLAPACK_SOURCE_TARBALL")")
    for source_extra in "${MPLAPACK_SOURCE_METADATA:-}" "${MPLAPACK_SOURCE_PATCH:-}" "${MPLAPACK_SOURCE_STATUS:-}"; do
        if [[ -n "$source_extra" && -f "$source_extra" ]]; then
            tar_args+=(-C "$(dirname "$source_extra")" "$(basename "$source_extra")")
        fi
    done
    tar -czf "$context_tar" "${tar_args[@]}"
else
    bundle_ref="refs/heads/mplapack-buildtest-${context_name}"
    git -C "$PROJECT_ROOT" update-ref "$bundle_ref" "$MPLAPACK_REF"
    git -C "$PROJECT_ROOT" bundle create "$source_bundle" "$bundle_ref"
    git -C "$PROJECT_ROOT" update-ref -d "$bundle_ref"
    tar -czf "$context_tar" \
        -C "$PROJECT_ROOT" release/docker \
        -C "$LOGDIR" "$(basename "$source_bundle")"
fi

remote_pid=""
cleanup_remote_command() {
    local rc=130
    local end_epoch start_epoch
    trap - INT TERM HUP
    if [[ -n "${remote_pid:-}" ]] && kill -0 "$remote_pid" 2>/dev/null; then
        kill "$remote_pid" 2>/dev/null || true
        wait "$remote_pid" 2>/dev/null || true
    fi
    end_epoch="$(date +%s)"
    start_epoch="${start:-$end_epoch}"
    write_log_end "$rc" "$end_epoch" "$((end_epoch - start_epoch))"
    link_failed_log "interrupted"
    exit "$rc"
}
trap cleanup_remote_command INT TERM HUP

start="$(date +%s)"
set +e
write_log_start "$start"
(
    "$REMOTE_SSH" "$host" "mkdir -p '$(dirname "$target_dir")' && cat > '$remote_context_tar'" < "$context_tar" &&     "$REMOTE_SSH" "$host" "MPLAPACK_REMOTE_WORKDIR='$target_dir' MPLAPACK_CONTEXT_TARBALL='$remote_context_tar' MPLAPACK_REF='$MPLAPACK_REF' MPLAPACK_SOURCE_MODE='$MPLAPACK_SOURCE_MODE' MPLAPACK_DISTRO_VERSION='$remote_distro_version' MPLAPACK_DOCKER_BASE='$remote_docker_base' $remote_cmd -s" < "$script_path"
) >> "$logfile" 2>&1 &
remote_pid="$!"
wait "$remote_pid"
rc=$?
remote_pid=""
trap - INT TERM HUP
set -e

end="$(date +%s)"
elapsed=$((end - start))
write_log_end "$rc" "$end" "$elapsed"
if [[ "$rc" -eq 0 ]]; then
    collect_remote_test_results
    compiler_label="$(detect_compiler_label_from_results "$COLLECTED_RESULTS_STAGE")"
    if [[ "$compiler_label" == "compiler_unknown" ]]; then
        echo "ERROR: compiler label could not be determined from collected results; refusing to create success stamp." >&2
        link_failed_log "compiler_unknown"
        mkdir -p "$(dirname "$resultfile")"
        echo "$matrix_name,$arch,${remote_docker_base:-${host:-}},test,FAILED,$elapsed,$source_type" | tee -a "$resultfile"
        exit 1
    fi
    final_logfile="$LOGDIR/${name_prefix}-${compiler_label}${source_log_part}-linux-docker.log"
    if [[ "$final_logfile" != "$logfile" ]]; then
        mv "$logfile" "$final_logfile"
        logfile="$final_logfile"
    fi
    stamp="${stamp_prefix}${compiler_label}${stamp_suffix}"
    clear_failed_logs
    ln -sfn "$logfile" "$stamp"
    mkdir -p "$(dirname "$resultfile")"
    echo "$matrix_name,$arch,${remote_docker_base:-$host},test,OK,$elapsed,$source_type" | tee -a "$resultfile"
else
    link_failed_log "failed"
    mkdir -p "$(dirname "$resultfile")"
    echo "$matrix_name,$arch,${remote_docker_base:-$host},test,FAILED,$elapsed,$source_type" | tee -a "$resultfile"
    exit "$rc"
fi
