# GOAL 11: GPU/CPU kernels — MPFR GPU inverses, BLAS1, BLAS2 triangular solves, and their square-size benchmarks

Prerequisites: goals 05-09 (GPU factors, solvers, dd/td/qd GPU inverses),
07a (cu_fdiv_rn), the kernel benchmark 5b5728451. Read 00-overview.md
("Result modes", "Naming and API", "GPU design") and 09-inverses.md first.
(Goal number 10 is reserved for the Cholesky speedup.)

The maintainer's plan has three stages: (1) implement the kernels, (2)
benchmark them alone at square sizes, (3) connect them (gates, LAPACK paths,
dispatcher defaults). **This goal is stages (1) and (2) only.** Stage (3) is a
later goal, after the maintainer has seen the benchmark.

## Scope

1. MPFR GPU inverses: `Rtrtri_mpfr_cuda`, `Rlauum_mpfr_cuda`,
   `Rpotri_mpfr_cuda`, `Rgetri_mpfr_cuda` and the complex `C*` ones
   (09-inverses.md covered dd/td/qd only; this goal adds MPFR).
2. BLAS1 GPU entry points for dd/td/qd/mpfr that are still missing:
   element-wise real `Raxpy`, `Rcopy` (and any other missing element-wise
   routine), `iRamax`, `iCamax`; and the reductions `Rdot`, `Cdotc`,
   `Cdotu`, `Rasum`, `RCasum`, `Rnrm2`, `RCnrm2` with an opt-in fast mode on
   the GPU and on the CPU (below). Not the scalar routines (rotg, rotm, rotmg).
3. BLAS2 GPU entry points `Rtrsv`, `Rtbsv`, `Rtpsv`, `Ctrsv`, `Ctbsv`,
   `Ctpsv` for dd/td/qd/mpfr. (The CPU `*_omp` versions exist.)
4. Square-size kernel benchmarks and graphs of 1-3.

Also in scope: `misc/gen_accelerator_headers.py` and the "Accelerator entry
points" rule of AGENTS.md (the result argument below), the "Result modes"
section of 00-overview.md (the fast switch), one Scope line in
09-inverses.md, the benchmark programs/scripts in `benchmark/`.

Out of scope: dispatcher defaults and gates (stage 3), optimizing existing
kernels, existing arithmetic, the known zero-sign allowance of the legacy
fixed-precision MPFR BLAS, the CPU small-n LU on a loaded host.

## Rules specific to this goal

- **Entry points of functions with a result** (dot, asum, nrm2, iamax, ...):
  the arguments of `<Routine>` in the same order, then **one reference
  argument `result`** of the routine's return type; the return value is
  `bool`. Example: `bool Rdot_dd_cuda(mplapackint const n, dd_real *dx,
  mplapackint const incx, dd_real *dy, mplapackint const incy, dd_real
  &result)`; `bool iRamax_dd_cuda(..., mplapackint &result)`. `result` is
  written only when true is returned. This is the form of the existing CPU
  `Rdot_par(..., REAL &result)`. Extend gen_accelerator_headers.py to
  generate these from non-void declarations without changing any existing
  generated line, and add the extension to AGENTS.md.
- **Fast mode of the reductions** (the only results allowed to differ from
  the reference): the same routine, switched by environment variables read
  on every call — CPU `MPLAPACK_MPBLAS_FAST=1`, GPU `<P>_MPBLAS_FAST=1`
  (`<P>` = MPLAPACK_DD, _TD, _QD, _MPFR). Off (default): CPU as today (bit
  identical), GPU entry returns false. On: parallel reductions, with a fixed
  summation order so that the same input and settings give the same result;
  nrm2 keeps the overflow/underflow-safe scaling. The default
  (bit-identical) mode must stay: it is the QA reference. No separately
  named routines.
- Everything else is bit-identical to the reference (zero signs, info,
  ipiv), as in 00-overview.md. Triangular solves: solve the columns in the
  reference order and update the independent elements of a column in
  parallel; where the reference forms a dot product (trans = T, C), keep its
  order and parallelize only the products.
- MPFR: fixed precisions (cu_freal, `cu_fdiv_rn`) and runtime precision
  (cu_mpfr), the eligibility, lock, stack (mpfr_cuda_stack.h) and memory
  policy of goals 07/07b.
- New GPU entry points are opt-in in the dispatchers (default: not used).

## Acceptance criteria

- Both build systems; host-emulation tests and GPU quick tests: bit-identical
  to the reference for everything except the fast reductions; fast
  reductions within a documented error bound (about c·n·u·Σ|x_i y_i|/|result|
  for dot/asum, a few u for nrm2) and reproducible run to run.
- nm -D: exactly the declared entry points added; check_unique_symbols.sh.
- Benchmarks: n = multiples of 7 and of 31 and powers of two (vectors up to
  about 2^24), GPU kernel only (data on the device, CUDA events, no
  transfers) and the CPU library (32 threads), graphs in the style of
  `benchmark/plot_sweep.py` (model: results/2026/<host>/rgemm/mpfr512_sweep10k.png),
  results in `benchmark/results/2026/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/kernel-bench-2/`.
- One commit per item (four): "MPFR CUDA: device-resident trtri/lauum/potri/getri (factorization goal 11)",
  "BLAS1: GPU entry points and CPU kernels (element-wise, iamax; dot/asum/nrm2 with an opt-in fast mode) (factorization goal 11)",
  "BLAS2: GPU entry points for trsv/tbsv/tpsv (factorization goal 11)",
  "Kernel benchmark: square-size sweeps of the new kernels (factorization goal 11)".
  A session instruction from the maintainer may give more detail; follow it where it does not contradict this file.
