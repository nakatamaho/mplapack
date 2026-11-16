# GOAL 02: Dedicated CPU Cholesky (Rpotrf / Cpotrf)

Prerequisites: goal 01 (commit 0cb1fc78). Read 00-overview.md first,
including "Known defects found by goals 01 and 02".

Status: done in f9216152. CPU correctness gates passed (goal 02 commit).
Performance gate: passed on interleaved re-measurement (2026-10-27,
A100 host shared with other workloads, load average about 100 on 64
hardware threads, so absolute times are not calibrated): Rpotrf lower,
default mode, OMP_NUM_THREADS=32, before = 1a9e8a7f, after = f9216152,
5 rounds alternating the order; median after/before dd 1.09 / 0.99 / 0.92,
qd 0.99 / 1.04 / 0.99, mpfr 1.05 / 1.07 / 1.00 at n = 512 / 1024 / 2048
(all <= 1.1; per-round ratios 0.72-1.24 without a consistent sign).
FAST=1 medians 0.95-1.06 except dd n = 512 (1.17, 0.07 s runs, ratios
0.81-1.68). Raw data: bench_local/results/goal02_potrf_interleaved.csv
in the working tree of that measurement (not committed).

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
   compiled into libmplapack_{dd,td,qd,mpfr}_opt_cuda. (Done in f9216152 with
   MPLAPACK_CHOLESKY_CPU_TASKS granted to the CPU libraries only; since the
   fix b626121f the accelerator libraries define MPLAPACK_ACCELERATOR_LIBRARY,
   and goal 03 replaces the former by the latter.)

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
