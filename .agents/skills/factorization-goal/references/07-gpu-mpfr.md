# GOAL 07: MPFR GPU Cholesky and LU (Rgetrf / Rpotrf / C*, mpfr)

Prerequisites: goals 05, 06, 07a. Read 00-overview.md (including "Known
defects found by goal 01"; the bug-fix track must be merged first) and
mpblas/optimized/mpfr/cuda/README.md first.

## Scope

Rgetrf_mpfr_cuda, Rpotrf_mpfr_cuda, Cgetrf_mpfr_cuda, Cpotrf_mpfr_cuda in
mpblas/optimized/mpfr/cuda (and mplapack-side dispatchers there).

## Tasks

1. Eligibility exactly like the MPFR BLAS entry points (same precision for
   every element and the default, MPFR_RNDN, finite values, exponent bounds;
   fixed precisions 256..2048 on cu_freal, others on cu_mpfr).
2. Bit identity: every operation correctly rounded as the CPU. Divisions
   and square roots of the panels on the fixed precisions use cu_fdiv_rn /
   cu_sqrt_rn (goal 07a); compare by measurement with (a) cu_mpfr on the
   device and (b) the panel on the host with the trailing update on the
   device (hybrid, MAGMA style), and document the choice. Updates (Rgemm, Rtrsm with diag = U, Rsyrk) on
   cu_freal.
3. Device memory check (an n = 10000, 512-bit matrix is 8 GB on the device)
   and `MPLAPACK_MPFR_LAPACK_CUDA*` variables as in goal 05.

## Acceptance criteria

- Rfactor_test CUDA variant for mpfr at 256, 512, 1024 bits and one runtime
  precision (e.g. 200 bits): bit-identical A, ipiv, info (zeros may differ
  only in sign on the cu_freal path, as documented for Rgemm — count them).
- Host emulation variant passes.
- Timings, mpfr 512 bits, n = 512, 1024, 2048, 4096: CPU default, current
  _opt_cuda, new entry point.
- nm -D, check_unique_symbols.sh.
- One commit: "MPFR CUDA: device-resident getrf/potrf (factorization goal 07)".
