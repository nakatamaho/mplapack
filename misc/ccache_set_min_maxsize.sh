#!/bin/bash
set -euo pipefail

min_size="${1:-80G}"

if ! command -v ccache >/dev/null 2>&1; then
    exit 0
fi

ccache_size_mib() {
    awk -v value="$1" '
    BEGIN {
        s = value
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
        lower = tolower(s)
        gsub(/[[:space:]]+/, "", lower)
        if (lower == "" || lower == "0" || lower == "unlimited") {
            print 1099511627776
            exit
        }
        n = lower
        sub(/[^0-9.].*$/, "", n)
        u = lower
        sub(/^[0-9.]+/, "", u)
        if (n == "" || n !~ /^[0-9]+(\.[0-9]+)?$/) {
            print -1
            exit
        }
        n += 0
        if (u == "" || u == "b") print n / 1024 / 1024
        else if (u == "k" || u == "kb" || u == "kib") print n / 1024
        else if (u == "m" || u == "mb" || u == "mib") print n
        else if (u == "g" || u == "gb" || u == "gib") print n * 1024
        else if (u == "t" || u == "tb" || u == "tib") print n * 1024 * 1024
        else print -1
    }'
}

current_size="$(ccache --get-config max_size 2>/dev/null || true)"
if [ -z "$current_size" ] && ccache -p >/dev/null 2>&1; then
    current_size="$(ccache -p 2>/dev/null | sed -n 's/^[[:space:]]*max_size[[:space:]]*=[[:space:]]*//p' | tail -n 1)"
fi

min_mib="$(ccache_size_mib "$min_size")"
current_mib="$(ccache_size_mib "$current_size")"

if awk -v c="$current_mib" -v m="$min_mib" 'BEGIN { exit !(c < 0 || c < m) }'; then
    ccache -M "$min_size"
else
    echo "ccache max_size already >= ${min_size}: ${current_size}"
fi
