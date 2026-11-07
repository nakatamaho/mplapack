# GOAL 12: fast cu_freal addition, subtraction and multiplication (mpc_cuda patch, bit-identical)

Prerequisites: goal 07a (the mpc_cuda patch mechanism,
external/mpc_cuda/patches/), goal 10 finished (its MPFR GEMM item 5 and the
upper Cholesky item 4b are committed or recorded as stalled), the two-level
MPFR roofline (commits 3e2ae9b67 and f8cc4b441:
benchmark/mpfr_path_bench.cu, benchmark/mpfr_gemm_paths.py,
benchmark/mpfr_level_peaks.py,
benchmark/results/2026/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/roofline/mpfr/README.md).
Read 00-overview.md and 07a-cu-freal-correct-rounding.md first.

## Why

Every MPLAPACK fixed-precision MPFR GPU kernel (PB = 256, 512, 768, 1024,
2048; MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS) computes with
cu_fp::cu_freal<PB> of the bundled mpc_cuda (include/mpc_cuda/cu_freal.cuh).
The two-level roofline measured on the A100 (roofline/mpfr/paths.csv):

| PB | cu_fmul / hardware full-product bound | cu_fadd (GEMM path mix) / hardware add bound | GEMM: level 2 / level 1 |
|---:|---:|---:|---:|
| 256 | 50% | 10% | 30% |
| 512 | 61% | 12% | 44% |
| 1024 | 62% | 13% | 52% |

Level 1 = the minimum 64-bit integer work (N^2 limb products with carries for
a product, N = PB/64; N funnel shifts + N adds with carry for a sum) at the
measured A100 instruction rates; level 2 = cu_freal as implemented. An
addition costs as much as a multiplication at PB = 256 although it needs
about 5 times fewer instructions. The causes visible in the source:

- `cu_freal<PB>::finalize` (normalize + round) and `cu_submag` index their
  limb arrays (D, A, T of length N+1 .. N+2) at run-time positions (the MSB
  position, the shift by the exponent gap), so nvcc puts the arrays in local
  memory: stack frames of 112 / 208 / 688 bytes at PB = 256 / 512 / 1024 in
  benchmark/mpfr_path_bench.cu (`-Xptxas -v`).
- The top-bit search and the sticky computation test bits / limbs one at a
  time in loops with run-time bounds.
- `cu_mul_full` forms the full 2N-limb product although rounding to PB bits
  needs only the upper part (plus a rare exact check).

No MPLAPACK-side change can recover this (goal 10 brought the MPFR GEMM to
about 54% of level 2 at 256/512 bits). Hence this goal: rewrite the
internals of cu_freal's +, -, x so that the results stay bit for bit the same.

## Why bit identity is preserved

cu_freal's +, - and x are correctly rounded (round to nearest, ties to even,
bit-exact with MPFR RNDN; 07a background). For a fixed PB the correctly
rounded result is unique, so any correct implementation returns the same
bits. The data layout, the zero convention (a flag, exp ==
CU_FREAL_ZERO_EXP, zero is +0, sign = 1), the absence of Inf/NaN, the
exponent convention and the unused low SB bits (SB = N*64 - PB, 0 or 32)
must stay exactly as they are. Every MPLAPACK result therefore stays
unchanged, which the existing bit-identity tests confirm.

## Maintainer decision (November 7, 2026): zero results

The goal's first run found that the old cu_fsub negates the sign of a zero b,
so +0 - +0 returns a zero with sign = -1 (MPFR RNDN: +0). This is an
upstream bug. It is fixed FIRST, by a separate small patch
`patch-cu_freal-zero-sign` (applied after the 07a patch, before the fast
patch; its own commit "mpc_cuda patch: canonical +0 for zero results of
cu_freal (factorization goal 12)"; a version against plain abd9407 for an
upstream PR independent of PR #1). Every zero result is then the canonical
+0 (set_zero()); non-zero results are unchanged. MPLAPACK kernels that keep
MPFR signed-zero semantics in the unused sign field of zero values must keep
their results (verify; stop if any MPLAPACK result changes). The fast
patch's oracle compares against cu_freal WITH the zero-sign fix (all fields)
and against MPFR including the sign of zero. Details:
bench_local/pi_instructions/goal12_addendum_zero.md.

## Policy: a second mpc_cuda patch, kept in MPLAPACK, offered upstream

As in 07a:

- New patch file `external/mpc_cuda/patches/patch-cu_freal-fast-add-mul`, a
  unified diff with an explanatory header, applied AFTER
  patch-cu_freal-correctly-rounded-div-sqrt (it may touch cu_freal.cuh only,
  plus a test/benchmark tool if useful upstream). Both build systems apply it
  at extract time (external/mpc_cuda/Makefile.am with the stamp depending on
  every patch file; the CMake ExternalProject PATCH_COMMAND); a failing patch
  stops the build. Do not change the tarball or its SHA256.
- The patch changes the IMPLEMENTATION of `operator+`, `operator-`,
  `operator*` of cu_freal (cu_addmag, cu_submag, cu_mul_full / the product,
  finalize), keeping every signature, every public name and every result.
  Functions that call finalize (cu_fdiv_rn, cu_sqrt_rn, conversions,
  cu_fround, precision conversion) must keep their results too; if finalize
  keeps a general version for them, that is fine.
- Update external/mpc_cuda/patches/README.md (what the patch does, the base
  commit abd9407 + the 07a patch, status "not yet submitted", the PR text).
  Do not open the upstream PR (needs the maintainer's account).
- No change to MPLAPACK kernel code is required; MPLAPACK kernels benefit
  through the operators. MPLAPACK-side changes are allowed only where a
  kernel works around the old costs in a way that is now harmful (measure
  first, record it).

## Implementation guidance (not prescriptive; results are what counts)

- Compile-time indexing: write the shifts and the normalization as
  unrolled loops over compile-time limb indices with run-time
  shift amounts selected by predication / funnel shifts (`__funnelshift_r`,
  `shf.r.wrap`), not as array accesses at a run-time index. Keep the
  operands and the working limbs in registers for PB <= 1024 (check
  `-Xptxas -v`: 0 bytes stack frame for the add and mul paths at
  256/512/1024 in the benchmark; 2048 may spill, record it).
- Addition paths by exponent gap d = ea - eb (|a| >= |b|): d = 0;
  1 <= d < 64; 64 <= d < PB (limb shift + bit shift); d >= PB + 2 (b only
  affects rounding as a sticky bit: the result is a or its neighbour,
  decided without the full shift). Same sign: one carry chain, at most one
  bit of renormalization. Opposite sign: borrow chain; for d >= 2 at most one
  bit of renormalization; for d in {0, 1} a leading-zero count over limbs
  (`__clzll`, limb-wise skip) and one funnel-shift pass.
- Rounding: round bit and sticky from the guard limbs with ties to even, as
  cu_freal does now (SB = 32: round at bit 32 of m[0]).
- Multiplication: the product only needs the bits that decide rounding. A
  truncated product (limb products with i + j >= N - 2, say, plus the
  carries into them) has a bounded error; when the truncated value is
  farther than that bound from a rounding boundary (the midpoint between two
  PB-bit numbers for RN), the rounding is decided; otherwise compute the
  full product (rare). Write down the error bound and the decision rule in
  the patch header and test the boundary cases explicitly. Use the PTX
  mad.lo.cc / madc.hi.cc chains as cu_mul_full does.
- Keep everything `__host__ __device__` (the host emulation tests use the
  same code).

## Tasks

1. Baseline: rebuild benchmark/mpfr_path_bench.cu against the current
   mpc_cuda and keep the numbers (roofline/mpfr/paths.csv is the reference;
   re-measure on a quiet GPU), plus the `-Xptxas -v` register/stack report of
   the add and mul paths.
2. Oracle test (new, kept in the patch as an mpc_cuda tool if it fits, and in
   MPLAPACK under mpblas/optimized/mpfr/cuda/test/): for PB in {256, 512,
   768, 1024, 2048} and two SB = 32 precisions (e.g. 224 and 992), host and
   device, compare the new +, -, x bit for bit against the UNPATCHED
   cu_freal (a copy of the old functions under another name inside the test)
   and against MPFR (mpfr_add / mpfr_sub / mpfr_mul, RNDN, precision PB):
   - random full-width significands, random signs, exponent gaps covering
     0, 1, 2, 62, 63, 64, 65, 127, 128, PB-2 .. PB+2, large;
   - exact ties (halfway cases), all-ones significands (carry out of every
     limb), powers of two, results that cancel to every length (1 bit .. PB
     bits, exact zero), operands equal in magnitude with opposite signs,
     products near a rounding boundary (constructed: x*y whose low part is
     0 or just above / below the midpoint), the smallest and largest
     exponents MPLAPACK's tests use;
   - at least 10^7 random pairs per PB and operation on the device,
     10^6 on the host; 0 mismatches required (sign, exponent, every limb,
     zero flag).
3. Implement the patch (task 2's test must pass before any performance
   claim).
4. Measure after: mpfr_path_bench (every path), the GEMM sweep
   (benchmark/run_gemm_sweep.sh, MPFR 256/512/1024, n = multiples of 7 and
   31 and powers of two), the kernel benchmark (benchmark/Rkernel_bench.cpp /
   kernel_bench.py) for gemm, syrk, cgemm, cherk, trsm and the panels, and
   Rgetrf/Rpotrf/Cgetrf/Cpotrf MPFR 512 through the GPU entry points
   (benchmark/Rfactor.cpp --entry, n = 1024, 2048, 4096; Upper and Lower),
   before and after, warmup 1, best of 3, alternating before/after runs.
5. Roofline: regenerate the level 2 peaks from the new paths.csv
   (benchmark/mpfr_level_peaks.py) into a new file next to the old one
   (keep both; the old level 2 is the "before" reference), and report
   kernel / level 1 before and after, and level 2 (new) / level 1.

## Acceptance criteria

- Task 2's oracle: 0 mismatches, host and device, all PB listed.
- All existing MPFR CUDA tests (fixed and runtime precision; BLAS, factor,
  solver, inverse, Rgemm/Winograd, cu_fround tests) pass in both build
  systems, bit-identical results unchanged; compute-sanitizer memcheck on one
  size of the MPFR real factor test (MPFR complex stays UNVERIFIED as decided
  in goal 10, bench_local/pi_instructions/goal10_addendum_mpc.md).
- Performance targets (guides; if missed, show with a profile what limits):
  cu_fadd (GEMM path mix) >= 40% of the level 1 add bound at 256/512/1024;
  cu_fmul >= 75% of the level 1 full-product bound (or report the truncated
  product's own bound); MPFR GEMM at n = 2048 (1024 bits: n = 1024) at least
  1.5x faster at 256 and 512 bits.
- Both patches apply in order on a fresh extract in both build systems;
  editing the patch file rebuilds mpc_cuda; nm -D of the MPLAPACK libraries
  unchanged; misc/check_unique_symbols.sh.
- Commit: "mpc_cuda patch: fast bit-identical cu_freal add/sub/mul
  (factorization goal 12)", body with the oracle counts, the path rates and
  the GEMM / factorization timings before and after, the roofline table, and
  the UNVERIFIED items.

Out of scope: the runtime-precision cu_mpfr (device MPFR port), cu_fdiv_rn /
cu_sqrt_rn algorithms (only keep their results), complex cu_fc arithmetic
beyond what follows from the real operators, fused multiply-add (the
reference rounds the product before the sum; a single-rounding fma would
change results), gates and dispatcher defaults.
