#!/bin/sh
# usage: check_unique_symbols.sh <lib...>
# Fails if any STRONG global defined symbol appears more than once in one
# library. Weak (W/w/V/v) and GNU-unique (u) definitions may legitimately
# repeat across TUs (C++ COMDAT) and are excluded.  PE/COFF nm output does not
# mark COMDAT definitions as weak, so MinGW archives are checked by asking the
# linker to do a relocatable whole-archive link.
status=0
NM_CMD=${NM:-nm}
OBJDUMP_CMD=${OBJDUMP:-objdump}
CXX_CMD=${CXX:-c++}
tmp_base=${TMPDIR:-/tmp}/mplapack-check-symbols.$$
nm_tmp=${tmp_base}.nm
objdump_tmp=${tmp_base}.objdump
link_tmp=${tmp_base}.o
link_err=${tmp_base}.err
trap 'rm -f "$nm_tmp" "$objdump_tmp" "$link_tmp" "$link_err"' EXIT HUP INT TERM

is_pe_coff_archive() {
  "$OBJDUMP_CMD" -f "$1" >"$objdump_tmp" 2>/dev/null &&
    grep -q '^In archive ' "$objdump_tmp" &&
    grep -Eqi 'file format (pe|pei)' "$objdump_tmp"
}

check_pe_coff_archive() {
  lib=$1
  rm -f "$link_tmp" "$link_err"
  # Intentionally leave CXX_CMD unquoted: autotools may set CXX to a command
  # with arguments, e.g. "ccache x86_64-w64-mingw32-g++-posix".
  if $CXX_CMD -no-pie -nostdlib -Wl,-r -Wl,--whole-archive "$lib" \
      -Wl,--no-whole-archive -o "$link_tmp" >"$link_err" 2>&1; then
    return 0
  fi

  if grep -qi 'multiple definition' "$link_err"; then
    echo "DUPLICATE SYMBOLS in $lib:"
    grep -i 'multiple definition' "$link_err"
    return 1
  fi

  echo "ERROR: unable to check duplicate symbols in $lib using CXX='$CXX_CMD'" >&2
  sed -n '1,40p' "$link_err" >&2
  return 1
}

for lib in "$@"; do
  if is_pe_coff_archive "$lib"; then
    if ! check_pe_coff_archive "$lib"; then
      status=1
    fi
    continue
  fi

  if "$NM_CMD" -m -g "$lib" >"$nm_tmp" 2>/dev/null && grep -q ' external ' "$nm_tmp"; then
    # Darwin/Mach-O reports C++ COMDAT/header-inline definitions as
    # "weak external" only in -m output; plain nm can show them as T.
    dups=$(awk '/ external / && $0 !~ /\(undefined\)/ && $0 !~ / weak / {print $NF}' "$nm_tmp" \
          | grep -v '^$' | sort | uniq -d)
  else
    dups=$("$NM_CMD" -g --defined-only "$lib" 2>/dev/null \
          | awk '$2 !~ /^[WwVvu]$/ {print $3}' | grep -v '^$' | sort | uniq -d)
  fi
  if [ -n "$dups" ]; then
    echo "DUPLICATE SYMBOLS in $lib:"; echo "$dups"; status=1
  fi
done
exit $status
