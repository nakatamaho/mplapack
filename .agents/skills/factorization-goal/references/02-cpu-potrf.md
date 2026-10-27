# GOAL 02: Dedicated CPU Cholesky (Rpotrf / Cpotrf)

Prerequisites: goal 01 (commit 0cb1fc78). Read 00-overview.md first,
including "Known defects found by goal 01".

## Scope

mplapack/optimized/common/{R,C}potrf.cpp and new helpers next to them
(Rpotrf2 / Cpotrf2 replacements may be added as files of the same name in
mplapack/optimized/common, shadowing the reference by basename).

## Tasks

1. Default mode (bit-identical): keep the reference operation order and
   improve only the schedule — parallel recursive diagonal block (the
   Rpotrf2 recursion's Rtrsm / Rsyrk / Rherk calls as tasks on tiles), no
   barrier between the diagonal block and the panel's Rtrsm tiles, tile
   size tuned per backend (dd/td/qd vs gmp/mpfr), and the existing
   look-ahead.
2. Fast mode (MPLAPACK_LAPACK_FAST=1): an OpenBLAS-style recursive Cholesky
   on the whole matrix (lapack/potrf/potrf_[UL]_parallel.c as the model):
   recursive split, Rsyrk / Rherk and Rtrsm on the halves with the
   optimized BLAS, its own block size.
3. Document the modes and variables in mplapack/optimized/common/README.md.
4. Accelerator libraries: mplapack/optimized/common/{R,C}potrf.cpp are also
   compiled into libmplapack_{dd,td,qd,mpfr}_opt_cuda, whose GPU BLAS calls
   are not safe to issue concurrently yet (Known defect 1). In those
   libraries keep Cholesky free of concurrent BLAS calls: no task schedule
   (the current production behaviour, la_lookahead(false, ...) without
   MPLAPACK_LAPACK_LOOKAHEAD=2) unless the bug-fix track has been merged
   (check git log for it and say which case applies in the report). Find
   or add a compile-time way to tell an accelerator build (e.g. a define
   passed only by the cuda/ Makefile.am and the CMake accelerator targets)
   rather than relying on the environment. Do not edit CUDA sources.

## Acceptance criteria

- Rfactor_test: bit-identical in default mode, backward error ≤ 30 in fast
  mode, for every backend, both UPLO, including info = j cases; the CUDA
  variants (dd, td, qd, mpfr; default and fast) still pass unchanged.
- mplapack/test/lin Rchkpo / Rdrvpo (and complex) pass for dd, qd, mpfr,
  double.
- Timings (commit body) for dd, qd, mpfr 512 bits at n = 512, 1024, 2048 on
  this machine's 32 cores: before (goal 01 numbers) and after, default and
  fast mode. The default mode must not be slower than before beyond noise
  (0.9-1.1) at any measured size.
- check_unique_symbols.sh passes.
- One commit: "Rpotrf, Cpotrf: dedicated CPU Cholesky (factorization goal 02)".
