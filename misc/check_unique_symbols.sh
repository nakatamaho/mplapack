#!/bin/sh
# usage: check_unique_symbols.sh <lib...>
# Fails if any STRONG global defined symbol appears more than once in one
# library. Weak (W/w/V/v) and GNU-unique (u) definitions may legitimately
# repeat across TUs (C++ COMDAT) and are excluded.
status=0
nm_tmp=${TMPDIR:-/tmp}/mplapack-check-symbols.$$
trap 'rm -f "$nm_tmp"' EXIT HUP INT TERM

for lib in "$@"; do
  if nm -m -g "$lib" >"$nm_tmp" 2>/dev/null && grep -q ' external ' "$nm_tmp"; then
    # Darwin/Mach-O reports C++ COMDAT/header-inline definitions as
    # "weak external" only in -m output; plain nm can show them as T.
    dups=$(awk '/ external / && $0 !~ /\(undefined\)/ && $0 !~ / weak / {print $NF}' "$nm_tmp" \
          | grep -v '^$' | sort | uniq -d)
  else
    dups=$(nm -g --defined-only "$lib" 2>/dev/null \
          | awk '$2 !~ /^[WwVvu]$/ {print $3}' | grep -v '^$' | sort | uniq -d)
  fi
  if [ -n "$dups" ]; then
    echo "DUPLICATE SYMBOLS in $lib:"; echo "$dups"; status=1
  fi
done
exit $status
