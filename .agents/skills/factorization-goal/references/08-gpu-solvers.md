# GOAL 08: GPU solvers (getrs, potrs, gesv, posv)

Prerequisites: goals 05, 06 (and 07 for mpfr). Read 00-overview.md first.

## Scope

Rgetrs/Rpotrs/Rgesv/Rposv and the complex ones as `<R>_<b>_cuda` for dd,
td, qd (mpfr if goal 07 is done), dispatchers, headers.

## Tasks

1. getrs / potrs entry points: upload factors and B, device Rlaswp and
   Rtrsm in reference order, download B.
2. gesv / posv entry points: one residency for factorization and solve (no
   download of the factors before the solve; download A, ipiv and B once).
3. Benchmarks for gesv / posv (nrhs = 1 and n) in the goal 01 programs.

## Acceptance criteria

- Bit-identical X, A, ipiv, info against the reference drivers; Rdrvge /
  Rdrvpo of mplapack/test/lin pass for dd, qd (mpfr if built) with the
  accelerator library.
- Timings for gesv / posv: CPU, per-call GPU BLAS path, entry point.
- nm -D, check_unique_symbols.sh.
- One commit: "CUDA getrs/potrs/gesv/posv entry points (factorization goal 08)".
