# GOAL SERIES: Unify MPLAPACK library layout (branch: topic/unified-libs)

This file is the shared context for goals 01-05. Read it before executing any
numbered goal. Do NOT implement anything from this file alone; it only defines
the target architecture and the invariants every goal must preserve.

## Motivation

Today each precision backend ships two libraries (libmpblas_<b> and
libmplapack_<b>), and the optimized layer ships yet more (libmpblas_<b>_opt,
libmpblas_dd_cuda, libmpblas_binary128_opt_opencl). The split has no benefit
and causes duplicate-symbol hazards: libmpblas_<b>_opt already re-bundles all
reference BLAS objects, so linking it together with libmpblas_<b> multiply
defines Rgemm, Rdot, etc. The CUDA library defines the public symbol Rgemm
again, conflicting with the _opt dispatcher.

## Target library taxonomy

Exactly ONE library per "flavor". A consumer links exactly one of these:

| Library                              | Contents                                                        |
|--------------------------------------|-----------------------------------------------------------------|
| libmplapack_<b>                      | mplapack/reference + mpblas/reference                           |
| libmplapack_<b>_opt                  | reference set, with basename-shadowed files replaced by mpblas/optimized/<b> (+ openmp helpers) and mplapack/optimized/<b> (skeleton, empty for now) |
| libmplapack_dd_opt_cuda              | the dd _opt set, with cuda-shadowed dispatchers replaced by mpblas/optimized/dd/cuda/*.cu |
| libmplapack_binary128_opt_opencl     | the binary128 _opt set, with opencl-shadowed dispatchers replaced by mpblas/optimized/binary128/opencl sources |

Backends <b>: gmp, mpfr, qd, dd, double, binary80, binary128 (same option
gating as today: MPLAPACK_ENABLE_* in CMake, --enable-* in configure.ac).

libmpblas_<b>, libmpblas_<b>_opt, libmpblas_dd_cuda are REMOVED. No stub
libraries, no compatibility symlinks. Migration is documented (goal 05).

## Symbol / shadowing policy (the fix for duplicate symbols)

Layer precedence, highest first, resolved per source-file stem (basename
without extension):

    accelerator (cuda/ or opencl/)  >  optimized  >  reference

Rules:

1. A library's source list is computed by starting from the reference list
   and DELETING every file whose stem also exists in a higher layer that is
   included in that library. Example: optimized/dd/Rgemm.cpp shadows
   reference/Rgemm.cpp in libmplapack_dd_opt; cuda/Rgemm.cu shadows
   optimized/dd/Rgemm.cpp in libmplapack_dd_opt_cuda.
2. Helper files never shadow and are never shadowed: anything matching
   `*_ref.cpp`, `*_omp.cpp`, `*_cuda.cu` (internal entry), and split-kernel
   files (`*_fermi_*`, `*_tesla_*`, `Rsyrk_[NT][LU]_*`). These define private
   symbols only. KNOWN EXCEPTION: `binary128/opencl/Rgemm_opencl.cpp` defines
   the PUBLIC symbol `Rgemm` despite its helper-looking name; goal 03 must
   either rename it to `Rgemm.cpp` (preferred, restores stem-shadowing) or
   register it as an explicit shadow of the reference `Rgemm`. Decide by `nm`
   on compiled objects, never by filename alone.
3. Every public API symbol (Rgemm, Rdot, ...) must be defined exactly once
   inside each installed library. Enforced by a check (see below), not by
   convention.
4. If an accelerator layer today defines the public symbol directly (cuda
   Rgemm.cu defines Rgemm), that is acceptable as the top layer, but goal 03
   must verify only one definition survives per library.

## Mandatory verification gate (used by every goal)

Add `misc/check_unique_symbols.sh` (goal 01 creates it) and run it on every
produced static archive and shared object:

```sh
#!/bin/sh
# usage: check_unique_symbols.sh <lib...>
# Fails if any STRONG global defined symbol appears more than once in one
# library. Weak (W/w/V/v) and GNU-unique (u) definitions may legitimately
# repeat across TUs (C++ COMDAT: inline functions, templates, vtables,
# function-local statics) and are excluded.
status=0
for lib in "$@"; do
  dups=$(nm -g --defined-only "$lib" 2>/dev/null \
        | awk '$2 !~ /^[WwVvu]$/ {print $3}' | grep -v '^$' | sort | uniq -d)
  if [ -n "$dups" ]; then
    echo "DUPLICATE SYMBOLS in $lib:"; echo "$dups"; status=1
  fi
done
exit $status
```

(For shared objects, duplicate definitions are already a link error; the
script mainly guards static archives where the linker silently picks one.)

## Invariants for ALL goals

- Branch: work only on `topic/unified-libs`.
- Both build systems (autotools AND CMake) must stay green after each goal;
  CMake parity is finalized in goal 04 but do not knowingly break it earlier.
- Do not modify any generated LAPACK routine source in mplapack/reference or
  mpblas/reference (they are FABLE-generated); only build files, new
  dispatchers, new directories, docs, and tests may change.
- Keep -ffp-contract=off handling for dd, OpenMP flags, macOS
  flat_namespace flags, and MinGW no-undefined flags exactly as they are
  attached to the current libraries, transplanted to the merged targets.
- ABI: this is a deliberate ABI/link-interface break. Bump the package
  version minor (2.2.x -> 2.3.0-dev) and libtool -version-info current
  (goal 05 finalizes numbers). Never pretend compatibility.
- Every goal ends with: `autoreconf` path builds (`./gen_configure.sh` if
  that is the repo convention), CMake configure+build for at least backends
  dd and double with tests ON, the unique-symbol check passing, and a
  focused commit with a message referencing the goal file.
- pkg-config (`mplapack.pc.in`, `cmake/mplapack.pc.cmake.in`), examples,
  benchmark, docker, and packaging link lines are updated in goal 05, but if
  a goal breaks their link lines, fix them in that goal.

## Execution order

01 -> 02 -> 03 -> 04 -> 05. Do not reorder; each assumes the previous state.
