#!/bin/sh
# usage: check_source_manifests.sh <top_srcdir>
# The 972-entry LAPACK reference source list exists twice in autotools
# (mplapack/reference/sources.am and accelerator-sources.am differ only in
# relative path prefix) while CMake uses file(GLOB). This check fails if the
# two manifests and the actual directory contents disagree by basename, so an
# upstream LAPACK sync cannot silently update one copy but not the other.
top="${1:?usage: check_source_manifests.sh <top_srcdir>}"
ref="$top/mplapack/reference"
tmp="${TMPDIR:-/tmp}/mplapack_manifest_$$"
mkdir -p "$tmp" || exit 1
trap 'rm -rf "$tmp"' EXIT

extract() { # extract() <fragment.am> -> sorted basenames of listed .cpp files
  sed -e 's/\\$//' "$1" | tr ' ' '\n' | grep '\.cpp$' | xargs -n1 basename | sort
}

extract "$ref/sources.am"             > "$tmp/opt"
extract "$ref/accelerator-sources.am" > "$tmp/accel"
ls "$ref"/*.cpp | xargs -n1 basename | sort > "$tmp/dir"

status=0
if ! cmp -s "$tmp/opt" "$tmp/accel"; then
  echo "MANIFEST DRIFT: sources.am vs accelerator-sources.am"
  diff "$tmp/opt" "$tmp/accel"; status=1
fi
if ! cmp -s "$tmp/opt" "$tmp/dir"; then
  echo "MANIFEST DRIFT: sources.am vs mplapack/reference/*.cpp"
  diff "$tmp/opt" "$tmp/dir"; status=1
fi
exit $status
