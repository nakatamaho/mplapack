#!/bin/sh
# usage: check_unique_symbols.sh <lib...>
# Fails if any STRONG global defined symbol appears more than once in one
# library. Weak (W/w/V/v) and GNU-unique (u) definitions may legitimately
# repeat across TUs (C++ COMDAT) and are excluded.
status=0
for lib in "$@"; do
  dups=$(nm -g --defined-only "$lib" 2>/dev/null \
        | awk '$2 !~ /^[WwVvu]$/ {print $3}' | grep -v '^$' | sort | uniq -d)
  if [ -n "$dups" ]; then
    echo "DUPLICATE SYMBOLS in $lib:"; echo "$dups"; status=1
  fi
done
exit $status
