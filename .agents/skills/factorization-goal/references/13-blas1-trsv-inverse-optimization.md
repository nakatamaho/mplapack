# GOAL 13: optimize the goal 11 kernels — iamax, the fast BLAS1 reductions, the triangular solves (trsv/tbsv/tpsv) and the inverses (trtri/lauum/potri/getri)

Prerequisites: goal 11 (6caad21c7 BLAS1, fc85906fc BLAS2, c1e4c603b MPFR
inverses, e1bba005b kernel benchmark, 1577e80c5 roofline), goal 09 (dd/td/qd
inverses), goal 10 (its GEMM/SYRK/HERK/TRSM kernels), goal 12 finished or
recorded as stalled (it changes cu_freal, which the MPFR measurements here
depend on). Read 00-overview.md ("Result modes", "GPU design"),
09-inverses.md, 11-gpu-kernels.md and mpblas/optimized/common/README.md
("Opt-in BLAS1 reductions (goal 11)") first.

This is stage 2.5 of the maintainer's kernel plan: the kernels exist and are
benchmarked; this goal makes them fast. Connecting them (gates, dispatcher
defaults) stays a later goal.

## Why (measured, A100 + Threadripper 3970X, goal 11 sweep and roofline)

Roofline (results/2026/<host>/roofline/goal11/README.md, MPFR on the level 2
peaks): element-wise axpy/copy reach 84-100% of the bandwidth for dd/td, but

| kernel | dd | qd | MPFR 512 | why (from the source) |
|---|---:|---:|---:|---|
| Rdot fast (2^21) | 6% | 6% | 3% | `Rblas1_jobs_common.h`: an O(n) array of `state` (re, im, scale, index), one global-memory pass per tree level, O(log n) launches |
| Rasum / Rnrm2 fast (2^21) | 3% | 3-6% | 1-2% | same |
| iRamax (2^21) | 2.7% | 2.7% | 1.7% | same tree for a reduction that needs no rounding |
| Rtrsv / Rtpsv (2048) | 0.2% | 0.4% | 0.2% | `Rtrisolve_jobs.h`: two kernel launches per column (pivot job + column job), 4096 launches at n = 2048; trans forms sum the column with one thread |
| Rtbsv (2048, k = 64) | 0.0% | 0.0% | 0.0% | same, with only k elements per launch |
| trtri / lauum (dd 1280/2048, qd 928/1280, MPFR 1022) | 0.1 / 1.5% | 0.2 / 1.1% | 0.1 / 0.3% | reference-blocked schedule, but TRMM/TRSM distribute only columns/rows and the leaves (trti2, lauu2) run per column |
| getri / potri (MPFR 1022 / 511) | not measured | not measured | 0.3 / 0.1% | same |

Wall times at the representative sizes (kernel only on the GPU, CPU library
with 32 threads; benchmark/results/2026/<host>/kernel-bench-2/kernels.csv):

| | GPU | CPU fast | CPU default |
|---|---:|---:|---:|
| dd Rdot 2^21 | 0.79 ms | **80.6 ms** | 4.4 ms |
| qd Rdot 2^21 | 1.63 ms | **145 ms** | 12.9 ms |
| MPFR512 Rdot 2^21 | 10.4 ms | **612 ms** | 47.8 ms |
| dd iRamax 2^21 | 0.88 ms | 28.1 ms (one thread path) | |
| dd Rtrsv 2048 | **15.0 ms** | 5.1 ms (CPU *_omp) | |
| dd Rtbsv 2048, k = 64 | **13.7 ms** | 0.41 ms | |
| MPFR512 Rtrtri 1022 | **33.5 s** | 12.4 s | |
| MPFR512 Rgetri 1022 | **35.1 s** | 15.6 s | |

So: the CPU fast mode is 10-20x SLOWER than the CPU default it is supposed
to beat; the GPU triangular solves and the MPFR GPU inverses are slower than
the CPU; the GPU reductions and iamax are an order of magnitude from the
bandwidth.

## Scope (in this order; one commit each)

1. **iamax (iRamax, iCamax), GPU and CPU, bit-identical by default.** The
   result is exact (comparisons of exact scores; iCamax's score
   |re| + |im| is computed per element exactly as the reference), so a
   parallel reduction gives the reference index when ties keep the first
   index. GPU: one pass over x with per-thread/per-warp/per-block (score,
   index) in registers and shared memory, then one small pass over the block
   results (no O(n) workspace, no per-level launches). CPU: OpenMP chunks
   with the same rule, combined in chunk order. Keep the NaN/positive
   increment contracts of goal 11. Target: >= 50% of the bandwidth on the
   GPU for dd at n >= 2^21 (MPFR: of the bandwidth or the level 2 bound,
   whichever binds); the CPU at least 8x faster than today at 2^21.
2. **Fast reductions** (Rdot, Cdotc, Cdotu, Rasum, RCasum, Rnrm2, RCnrm2;
   only with `MPLAPACK_MPBLAS_FAST=1` / `<P>_MPBLAS_FAST=1`; the default
   stays the reference order on the CPU and false on the GPU). Replace the
   O(n)-state pairwise tree by a **fixed two-level order that depends on n
   only**: contiguous chunks of a fixed length C (e.g. 1024 or 4096; one
   constant, documented) are summed in a fixed order (sequential or a fixed
   in-chunk tree), and the chunk results are combined by the existing
   adjacent-pair tree. The same order on the CPU and on the GPU, so the CPU
   and GPU fast results are identical and do not depend on thread counts or
   launch shapes (keep this property, and test it). nrm2 keeps the
   (scale, sumsq) safe scaling. The fast results may differ from goal 11's
   fast results (the order changes): update the documented order and error
   bound in mpblas/optimized/common/README.md. Targets: GPU >= 50% of the
   bandwidth for dd/td/qd at n >= 2^21 (MPFR: the binding bound); CPU fast
   faster than CPU default by at least 8x for dd/qd/MPFR at 2^21 with 32
   threads.
3. **Triangular solves** (Rtrsv, Rtbsv, Rtpsv, Ctrsv, Ctbsv, Ctpsv; all
   uplo/trans/diag; dd/td/qd/mpfr; bit-identical). In the reference every
   element receives its terms in a fixed order: NoTrans (axpy form) — x(i)
   gets x(i) - x(j)*A(i,j) for j in column order (and the zero skip when
   x(j) == 0), then its division when its own column comes; Trans/ConjTrans
   (dot form) — temp = x(j), minus A(i,j)*x(i) for i in the reference order,
   then the division. A term becomes available as soon as the x it uses is
   final. Hence a **wavefront**: one persistent kernel (one block, or a
   grid with device-side progress flags; no launch per column) in which the
   owner of each element applies its terms in the reference order as they
   become available; the column solves proceed one after another but all
   pending updates of all elements run concurrently. Every element keeps
   exactly the reference sequence, so the result is bit-identical. Band and
   packed storage use the same scheme with their index maps. Target: dd
   Rtrsv n = 2048 at least 10x faster than today and faster than the CPU
   *_omp; Rtbsv proportionally (O(n·k) work); MPFR likewise. Also apply the
   wavefront to the CPU *_omp versions if they are slower than the GPU
   result suggests is possible (measure first; optional).
4. **Inverses** (Rtrtri, Rlauum, Rpotri, Rgetri and C*; dd/td/qd/mpfr;
   bit-identical; the reference iMlaenv blocking stays). The reference
   blocked algorithms call Rtrmm/Rtrsm/Rtrti2 (trtri), Rtrmm/Rgemm/Rsyrk or
   Rherk/Rlauu2 (lauum), Rgemv/Rgemm/Rtrsm/Rswap (getri). Profile first
   (nsys) which part dominates; then give every call the level of its
   goal 10 counterpart: TRSM with the goal 10 blocked team kernels, GEMM and
   SYRK/HERK tiled, TRMM tiled the same way (each output element keeps its
   reference sequence of terms), the leaves (trti2, lauu2) with one block per
   diagonal block instead of one launch per column, GEMV as one kernel over
   the output elements. Target: the GPU faster than the CPU library (32
   threads) at n >= 1024 for every type, and the % of roofline before and
   after at n = 512, 1024, 2048 (MPFR 1024: 512, 1024).
5. **Benchmark and roofline** of 1-4 before and after: the kernel benchmark
   (benchmark/Rkernel_bench.cpp / kernel_bench.py; vectors n = multiples of
   7 and 31 and powers of two up to 2^24, matrices up to 2048-4096), CPU and
   GPU, graphs like the goal 11 ones, roofline with benchmark/roofline.py
   (peaks.csv: MPFR level 2 peaks; also the level 1 peaks), results in
   `benchmark/results/2026/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/kernel-bench-3/`.

Out of scope: gates and dispatcher defaults; cu_freal arithmetic (goal 12);
the default-mode GPU reductions (they keep returning false: a reference-order
sum is sequential); new fast modes for other routines; the known zero-sign
allowance of the legacy fixed-precision MPFR BLAS; the mpc_cuda complex
division defect (MPFR complex memcheck stays UNVERIFIED as decided in goal 10).

## Rules

- Everything except the fast reductions is bit-identical to the reference
  (values, zero signs, info, the iamax index), host emulation and GPU.
- The fast reductions: deterministic, the same result on the CPU and the
  GPU, independent of thread counts, block counts and launch shapes; error
  bound documented (about c·n·u·Σ|x_i y_i| / |result| for dot/asum with the
  new c; a few u for nrm2); tested.
- Entry point signatures, eligibility variables and the goal 11 contracts
  (false leaves all outputs untouched, no CPU fallback inside an entry)
  stay. No new public symbols unless needed (then generated by
  misc/gen_accelerator_headers.py and declared).
- MPFR: fixed precisions (cu_freal, cu_fdiv_rn) and runtime precision
  (cu_mpfr); the eligibility, lock, stack (mpfr_cuda_stack.h), arena and
  memory policies of goals 07/07b. **Runtime precision has bounded scratch
  per lane: any change that runs more lanes concurrently must stay within
  it** (goal 10 item 2 broke this; see goal 12 addendum 3).
- Both build systems.

## Acceptance criteria (per item)

- Bit identity as above on all forms touched: n = 0, 1, block/chunk
  boundaries (C-1, C, C+1, multiples), lda > n, negative and non-unit
  increments where the contract allows them, zeros (zero skips, signed
  zeros), complex. For iamax: ties at block/chunk boundaries, all-equal
  vectors, the maximum at the first and last index.
- Fast reductions: CPU result == GPU result bit for bit for the same n and
  input, for 1, 2, 7, 32 threads and several launch shapes; error bound
  tests.
- All existing tests that use the kernels pass, **run to completion**,
  including the MPFR runtime-precision solver/inverse tests at 768/1024/2048
  bits (one to two hours each). compute-sanitizer memcheck on one size per
  item (MPFR complex: UNVERIFIED as decided).
- Performance before/after with the kernel benchmark and the roofline at the
  sizes of item 5; report the targets met and missed (with a profile for the
  missed ones).
- nm -D unchanged (or exactly the declared additions); check_unique_symbols.sh.
- Commits: "BLAS1: one-pass parallel iamax on the GPU and the CPU, bit-identical (factorization goal 13)",
  "BLAS1 fast reductions: fixed two-level order, identical on the CPU and the GPU (factorization goal 13)",
  "BLAS2: wavefront trsv/tbsv/tpsv on the GPU, bit-identical (factorization goal 13)",
  "Inverses: trtri/lauum/potri/getri on the tiled GPU kernels, bit-identical (factorization goal 13)",
  "Kernel benchmark: goal 13 before/after and roofline (factorization goal 13)".
  Each body: before/after times and % of roofline, tests, UNVERIFIED items.
  A session instruction from the maintainer may give more detail; follow it
  where it does not contradict this file.
