#!/bin/sh
# usage: check_unique_symbols.sh <lib...>
# Fails if any strong global defined symbol appears more than once in one library.
status=0
for lib in "$@"; do
  dups=$(nm -g --defined-only "$lib" 2>/dev/null \
        | awk '$2 !~ /^[WwVv]$/ {print $3}' | grep -v '^$' | sort | uniq -d)
  if [ -n "$dups" ]; then
    echo "DUPLICATE SYMBOLS in $lib:"; echo "$dups"; status=1
  fi
done
exit $status
