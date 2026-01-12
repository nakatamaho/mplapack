#!/usr/bin/env bash
set -euo pipefail

USE_CCACHE=${USE_CCACHE:-yes}
PREFIX=${PREFIX:-"$HOME/MPLAPACK"}

# Move to repository root (directory containing configure.ac)
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir/.."

# Optional: ccache
if [[ "$USE_CCACHE" == "yes" ]]; then
  command -v ccache >/dev/null 2>&1 || { echo "ccache not found"; exit 1; }
  export CC="ccache gcc"
  export CXX="ccache g++"
  export FC="ccache gfortran"
  ccache -M 80G
else
  export CC="gcc"
  export CXX="g++"
  export FC="gfortran"
fi

# Generate Makefile.am (if your project needs it)
pushd mplapack/test/compare >/dev/null
bash gen.Makefile.am.sh
popd >/dev/null

# Regenerate autotools files (libtoolize must happen before automake)
rm -rf autom4te.cache
mkdir -p m4

libtoolize --force --copy
autoreconf -fi

# Configure options
common_opts=(
  "--prefix=$PREFIX"
  "--enable-gmp=yes"
  "--enable-mpfr=yes"
  "--enable-_Float128=yes"
  "--enable-qd=yes"
  "--enable-dd=yes"
  "--enable-double=yes"
  "--enable-test=yes"
  "--enable-benchmark=yes"
)

# _Float64x is typically x86_64-only in your logic
if [[ "$(uname -m)" == "x86_64" ]]; then
  common_opts+=("--enable-_Float64x=yes")
fi

./configure "${common_opts[@]}"
