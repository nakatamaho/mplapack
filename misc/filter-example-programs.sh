#!/bin/bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
    echo "usage: $0 <makefile> [disabled-suffix ...]" >&2
    exit 2
fi

makefile="$1"
shift

if [ "$#" -eq 0 ]; then
    exit 0
fi

tmp="${makefile}.tmp.$$"
awk -v disabled_suffixes="$*" '
BEGIN {
    ndisabled = split(disabled_suffixes, suffixes, /[[:space:]]+/)
    for (i = 1; i <= ndisabled; i++) {
        disabled[suffixes[i]] = 1
    }
}

function program_suffix(program, base, parts, n) {
    base = program
    sub(/\.exe$/, "", base)
    sub(/_opt$/, "", base)
    n = split(base, parts, "_")
    return parts[n]
}

function emit_programs(text, tokens, ntokens, i, output, suffix) {
    ntokens = split(text, tokens, /[[:space:]]+/)
    output = ""
    for (i = 1; i <= ntokens; i++) {
        if (tokens[i] == "") {
            continue
        }
        suffix = program_suffix(tokens[i])
        if (suffix in disabled) {
            continue
        }
        output = output (output == "" ? "" : " ") tokens[i]
    }
    print "programs=" output
}

/^programs[[:space:]]*=/ {
    text = $0
    sub(/^[^=]*=/, "", text)
    while ($0 ~ /\\[[:space:]]*$/) {
        sub(/\\[[:space:]]*$/, "", text)
        if ((getline) <= 0) {
            break
        }
        text = text " " $0
    }
    emit_programs(text)
    next
}

{
    print
}
' "${makefile}" > "${tmp}"
mv "${tmp}" "${makefile}"
