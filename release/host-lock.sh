#!/bin/bash

mplapack_host_lock_key() {
    printf '%s' "$1" | sed 's/[^A-Za-z0-9_.-]/_/g'
}

mplapack_host_lock_acquire() {
    local host="$1"
    local label="${2:-$1}"
    local lock_root="${MPLAPACK_HOST_LOCK_DIR:-${TMPDIR:-/tmp}/mplapack-release-host-locks-${USER:-user}}"
    local poll_seconds="${MPLAPACK_HOST_LOCK_POLL_SECONDS:-30}"
    local lock_key lockdir old_pid old_label waited

    mkdir -p "$lock_root"
    lock_key="$(mplapack_host_lock_key "$host")"
    lockdir="$lock_root/${lock_key}.lock"
    waited=0

    while ! mkdir "$lockdir" 2>/dev/null; do
        old_pid=""
        old_label=""
        [ -f "$lockdir/pid" ] && old_pid="$(sed -n '1p' "$lockdir/pid" 2>/dev/null || true)"
        [ -f "$lockdir/label" ] && old_label="$(sed -n '1p' "$lockdir/label" 2>/dev/null || true)"

        if [ -z "$old_pid" ] || ! kill -0 "$old_pid" 2>/dev/null; then
            echo "WARNING: removing stale host lock for $host: $lockdir" >&2
            rm -rf "$lockdir"
            continue
        fi

        if [ "$waited" -eq 0 ]; then
            echo "Host $host is busy with ${old_label:-pid $old_pid}; waiting for host lock." >&2
        fi
        waited=$((waited + poll_seconds))
        sleep "$poll_seconds"
    done

    printf '%s\n' "$$" > "$lockdir/pid"
    printf '%s\n' "$label" > "$lockdir/label"
    date '+%Y-%m-%d %H:%M:%S %Z' > "$lockdir/started_at" 2>/dev/null || true

    MPLAPACK_HELD_HOST_LOCKS="${MPLAPACK_HELD_HOST_LOCKS:-}${MPLAPACK_HELD_HOST_LOCKS:+
}$lockdir"
    echo "Acquired host lock for $host: $label" >&2
}

mplapack_host_lock_release_all() {
    local lockdir pid

    [ -n "${MPLAPACK_HELD_HOST_LOCKS:-}" ] || return 0
    while IFS= read -r lockdir; do
        [ -n "$lockdir" ] || continue
        pid=""
        [ -f "$lockdir/pid" ] && pid="$(sed -n '1p' "$lockdir/pid" 2>/dev/null || true)"
        if [ "$pid" = "$$" ]; then
            rm -rf "$lockdir"
        fi
    done <<EOF
$MPLAPACK_HELD_HOST_LOCKS
EOF
    MPLAPACK_HELD_HOST_LOCKS=""
}
