# GOAL 06: GPU LU entry points (Rgetrf / Cgetrf, dd td qd)

Prerequisites: goals 03, 04, 05. Read 00-overview.md (including "Known
defects found by goal 01"; the bug-fix track must be merged first) first.

## Scope

Rgetrf_<b>_cuda, Cgetrf_<b>_cuda (b = dd, td, qd), dispatchers, headers,
both build systems.

## Tasks

1. Entry point as in goal 05, LU in reference order on the device: Rgetrf2
   recursion on the panel (device iRamax, swaps, scaling), device Rlaswp on
   the left and right blocks, Rtrsm, Rgemm; look-ahead with two streams;
   ipiv downloaded with A; info as the reference (singular pivots do not
   stop the factorization).
2. Fast mode `<P>_LAPACK_FAST=1` as in goal 05.

## Acceptance criteria

- Rfactor_test CUDA variant: bit-identical A, ipiv, info for the whole
  matrix set (rectangular, singular, ties, lda > m); false for invalid
  arguments and quick returns.
- Host emulation variant passes.
- Timings as in goal 05 for getrf.
- nm -D, check_unique_symbols.sh.
- One commit: "Rgetrf/Cgetrf_<b>_cuda: device-resident LU (factorization goal 06)".
