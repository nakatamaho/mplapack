# GOAL 03: Dedicated CPU LU (Rgetrf / Cgetrf, Rgetrf2 / Cgetrf2)

Prerequisites: goals 01, 02. Read 00-overview.md first, including "Known
defects found by goal 01": in the accelerator libraries the LU task
schedule must not issue concurrent GPU BLAS calls until the bug-fix track
is merged (same rule as goal 02, task 4).

## Scope

mplapack/optimized/common/{R,C}getrf.cpp, new {R,C}getrf2.cpp (and getf2
if used) in the same directory, Rlookahead_common.h.

## Tasks

1. Default mode (bit-identical): parallel recursive panel. Rgetrf2's
   recursion stays; its column searches (iRamax / iCamax, first index on
   ties), swaps and scalings are split across threads by row blocks with a
   reduction that returns exactly the reference pivot; its Rtrsm / Rgemm
   calls run on the optimized BLAS. Apply each tile's row interchanges,
   Rtrsm and Rgemm as one task (fused per tile, still in reference order).
2. Fast mode (MPLAPACK_LAPACK_FAST=1): OpenBLAS-style LU
   (lapack/getrf/getrf_parallel.c as the model): recursive panel of its own
   width, dynamic look-ahead, updates fused with the BLAS engine's packing
   where the optimized Rgemm exposes it.
3. Rgetrs keeps the reference code (it calls the optimized Rlaswp / Rtrsm);
   check that Rgesv picks up the new Rgetrf.

## Acceptance criteria

- Rfactor_test: bit-identical default mode (A, ipiv, info), backward error
  ≤ 30 and identical pivots on non-tie matrices in fast mode; singular and
  rectangular cases included.
- mplapack/test/lin Rchkge / Rdrvge (and complex) pass for dd, qd, mpfr,
  double.
- Timings (commit body) as in goal 02; default mode never slower than
  before beyond noise.
- check_unique_symbols.sh passes.
- One commit: "Rgetrf, Cgetrf: dedicated CPU LU with parallel panel (factorization goal 03)".
