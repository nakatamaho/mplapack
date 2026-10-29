# GOAL SERIES: dedicated LU and Cholesky (CPU and GPU)

Shared context for goals 01-09 and 07a. Read it before any numbered goal. Do not
implement anything from this file alone.

## Motivation

OpenBLAS does not run the reference LAPACK getrf/potrf: it has its own
lapack/ directory (getrf, getf2, getrs, laswp, potrf, potf2, lauum, lauu2,
trtri, trti2, trtrs) with recursive, parallel panels and update steps
fused with its packed GEMM, and the gesv/posv drivers on top of them.

MPLAPACK today (branch claude/mpfr-cuda-rgemm):

- mplapack/optimized/common/{R,C}getrf.cpp and {R,C}potrf.cpp run the
  reference blocked algorithm (same nb from iMlaenv, Rgetrf2 / Rpotrf2
  panels, Rlaswp, Rtrsm, Rgemm, Rsyrk) with OpenMP task look-ahead
  (Rlookahead_common.h); results are the reference results bit for bit.
  MPLAPACK_LAPACK_FAST=1 runs the recursive Rgetrf2 / Rpotrf2 on the whole
  matrix instead (other rounding). {R,C}laswp.cpp are optimized too.
- The accelerator libraries (libmplapack_{dd,td,qd,mpfr}_opt_cuda) compile
  those same files: the factorization runs on the CPU, and each BLAS call it
  makes (Rgemm, Rtrsm, Rsyrk, ...) copies its operands to the GPU and back.
  There is no GPU panel, no GPU Rlaswp, and the matrix never stays on the
  device across calls.
- Every other routine of the table below runs the reference code.

In multiple precision the factorizations are dominated by the arithmetic of
their Level 3 updates, so the CPU gains of dedicated code come from the
panel (parallel pivot search and scaling, recursion), fewer barriers and
better use of the cores for small and medium n; the GPU gains come from
keeping the matrix on the device for the whole factorization (one upload,
one download) and running the panel there too.

## Known defects found by goals 01 and 02 (fixed in b626121f)

Goals 01 and 02 found defects of the accelerator libraries; they were fixed
outside this series in commit b626121f ("CUDA/OpenCL: make GPU BLAS calls
safe under the LAPACK look-ahead tasks"):

1. dd/td/qd CUDA: each GPU BLAS call wrote back the whole memory span of its
   operands, overwriting rows that concurrent look-ahead tasks were updating
   (lower Rpotrf/Cpotrf with MPLAPACK_LAPACK_LOOKAHEAD=2). Now only the
   referenced elements are transferred (gpu_operand rows/cols/region,
   Rblas_gpu_operand_common.h).
2. MPFR CUDA: concurrent calls raced on mpc_cuda's arena globals, per-call
   device limits and the complex bridge's shared store (illegal memory
   access in Rgetrf n = 1024, 2048). The MPFR GPU entry points are
   serialized by a mutex.
3. MPFR CUDA complex BLAS: the fixed-precision Ctrsm overflowed the default
   1 KB device stack even in serial calls (found by goal 02). The bridge
   raises the stack limit before its launches.
4. MPFR CUDA Rswap/Rrot: x was not copied back from the device.

Consequences for this series:
- The accelerator libraries (compiled with MPLAPACK_ACCELERATOR_LIBRARY)
  run the LAPACK look-ahead tasks only with MPLAPACK_LAPACK_LOOKAHEAD=2; by
  default their factorizations are serial (each task's BLAS call would copy
  its operands to the GPU and wait for the others: dd LU n = 4096 26.5 s
  with the old default, 5.4 s now). Concurrent GPU BLAS calls are correct
  but serialized; do not design GPU goals around them.
- MPLAPACK_ACCELERATOR_LIBRARY is the one compile-time switch for "this is
  an accelerator library". Goal 02 introduced MPLAPACK_CHOLESKY_CPU_TASKS
  (granted to the CPU libraries only) before the fix existed; goal 03
  replaces it by !MPLAPACK_ACCELERATOR_LIBRARY (one switch, same effect).
- The CUDA files of the fix (mpblas/optimized/common/*cuda*, blas_cuda/,
  {dd,td,qd,mpfr}/cuda/, Rlookahead_common.h) are no longer off limits.

Provisional `<P>_LAPACK_CUDA_MIN` planning defaults from goal 01 (GPU BLAS
per call; retune after device residency, goals 05-07): MPLAPACK_DD 4096,
MPLAPACK_QD 512, MPLAPACK_TD 512 (QD proxy, not measured), MPLAPACK_MPFR 2048
(MPFR LU timings above 512 were unverified because of defect 2; remeasure).

## Routine table (what to implement)

Backends: CPU paths are generic (mplapack/optimized/common, every backend);
GPU paths are dd, td, qd (CUDA, BF arithmetic) and mpfr (CUDA, cu_freal /
cu_mpfr). P0 = first, P1 = next, P2 = later.

| Routine (R and C) | Role | Today | CPU: implement | GPU: implement | Prio | Goal |
|---|---|---|---|---|---|---|
| getrf | LU with partial pivoting | look-ahead, bit-identical | dedicated recursive LU (parallel panel, fused swap/trsm/gemm per tile), opt-in fast mode | `Rgetrf_<b>_cuda`: device-resident, bit-identical default | P0 | 03, 06, 07 |
| getrf2 / getf2 | LU panel | reference | parallel recursive panel (pivot search, swaps, scaling across threads) | panel kernels on the device (iamax, swap, scal/div, rank-1 / small gemm) | P0 | 03, 04, 06 |
| laswp | row interchanges | optimized | (done) | device kernel, in place on the device matrix | P0 | 04 |
| potrf | Cholesky | look-ahead, bit-identical | dedicated recursive Cholesky, opt-in fast mode | `Rpotrf_<b>_cuda`: device-resident, bit-identical default | P0 | 02, 05, 07 |
| potrf2 / potf2 | diagonal block | reference | parallel recursive diagonal block | device kernels (sqrt, scaling, syrk/herk on the block) | P0 | 02, 04, 05 |
| getrs / potrs | solve with the factors | reference (calls opt BLAS) | keep (BLAS-bound) | `Rgetrs_<b>_cuda`, `Rpotrs_<b>_cuda` | P1 | 08 |
| gesv / posv | factor + solve | reference | keep (calls the optimized getrf/potrf) | `Rgesv_<b>_cuda`, `Rposv_<b>_cuda`: factor and solve without leaving the device | P1 | 08 |
| trtri, lauum (trti2, lauu2) | triangular inverse, U·U^T | reference | dedicated recursive versions | `Rtrtri_<b>_cuda`, `Rlauum_<b>_cuda` | P2 | 09 |
| potri, getri | inverse from the factors | reference | keep (calls trtri / lauum) | `Rpotri_<b>_cuda`, `Rgetri_<b>_cuda` | P2 | 09 |
| trtrs | triangular solve with checks | reference | keep | (Rtrsm is already on the GPU) | — | — |
| (support) cu_freal div / sqrt | correctly rounded division and square root for the MPFR fixed-precision device type | additive cu_fround.cuh patch provides RNDN APIs; cu_fmath stays approximate | — | `cu_fdiv_rn`, `cu_sqrt_rn` available as a patch to the bundled mpc_cuda (upstream proposal prepared, not submitted, carried on 3.1); MPFR Rtrsm on cu_freal | P0 for mpfr | 07a |

Out of scope: QR, eigenvalue and SVD routines, pivoted Cholesky (pstrf),
band and packed storage, mixed precision refinement.

## Result modes (normative)

1. **Default = reference results bit for bit.** Every public routine
   (`Rgetrf`, `Rpotrf`, ...) and every GPU entry point in its default mode
   performs, element by element, the operations of the reference routine
   in their order: same nb (iMlaenv), same recursion of Rgetrf2 / Rpotrf2,
   same Rlaswp, Rtrsm, Rgemm, Rsyrk calls. Reordering the *scheduling*
   (tasks, look-ahead, streams, which device computes a tile) is allowed;
   reordering *arithmetic* is not. The dd/td/qd GPU kernels are already
   bit-identical to the CPU (BF arithmetic, `-fmad=false`, accurate_div);
   this is what makes a bit-identical GPU factorization possible.
2. **Fast mode = opt-in, other rounding.** `MPLAPACK_LAPACK_FAST=1` (CPU,
   existing switch) and `<P>_LAPACK_FAST=1` for a GPU entry point (`<P>` =
   MPLAPACK_DD, _TD, _QD, _MPFR) may use any algorithm (recursion on the
   whole matrix, other block sizes, OpenBLAS-style fused updates, Winograd
   or sloppy GEMM). Gate: backward error (below), and the same pivots as
   the default mode on the test matrices whenever the pivot choice is not
   a tie.
3. No silent CPU fallback inside a GPU entry point (AGENTS.md): it returns
   false and leaves every output untouched; the dispatcher then computes on
   the CPU.

## Naming and API (AGENTS.md, "Accelerator entry points")

- GPU entry points: `<Routine>_<backend>_cuda` with the arguments of the
  routine and a bool result, e.g.
  `bool Rgetrf_dd_cuda(mplapackint const m, mplapackint const n, dd_real *a,
  mplapackint const lda, mplapackint *ipiv, mplapackint &info);`
  (the declaration of Rgetrf in include/mplapack_dd.h with `bool` and the
  suffix).
  Declared in include/mpblas_<backend>_cuda.h (rename to the LAPACK side if
  a goal decides so; keep one header per library), generated by
  misc/gen_accelerator_headers.py — extend its routine lists and teach it to
  read include/mplapack_<backend>.h for LAPACK routines.
- info semantics: identical to the reference (singular U(i,i) = 0 gives
  info = i, the factorization completes; potrf stops at the first non
  positive pivot with info = j). Argument errors: return false (the
  dispatcher's reference call reports them through Mxerbla).
- The dispatcher (`Rgetrf` in the accelerator library, a file under the
  library's cuda/ directory shadowing mplapack/optimized/common/Rgetrf.cpp)
  calls the entry point, then the CPU path.
- Internal helpers live in namespaces (mplapack_lapack_gpu, ...), never
  with names of the entry point form.
- Environment variables (read at every call, documented next to the code
  and in the README of the directory): `<P>_LAPACK_CUDA=0` (never),
  `<P>_LAPACK_CUDA_MIN=N` (smallest n, default chosen from the goal 01
  measurements), `<P>_LAPACK_FAST=1`, `<P>_CUDA_VERBOSE=1`.

## GPU design (normative for goals 04-08)

- One call = one device residency: allocate, upload A (and B), factor and
  solve on the device, download, free. Device memory per call: the matrix
  plus O(n·nb) workspace; check it against cudaMemGetInfo and return false
  when it does not fit.
- Element storage on the device: the planar / packed layouts the existing
  kernels use (dd/td/qd: common/Rgemm_cuda_common.cuh, Rblas_cuda_common.cuh;
  mpfr: cu_freal<PB> for 256..2048 bits, cu_mpfr otherwise). Reuse the
  device BLAS bodies (blas_dev/, Rblas_cuda_common.cuh, gemm_cuda tiled
  kernels) through device-pointer variants; do not duplicate arithmetic.
- Panels on the device. Pivot search: |a| comparisons exactly as the
  reference (iRamax / iCamax semantics, first index on ties). Row swaps:
  device Rlaswp. Scaling: the reference test `|pivot| >= sfmin` then
  multiply by the reciprocal, else divide — dd/td/qd divide with
  accurate_div (bit-identical); for mpfr, the bundled additive mpc_cuda
  patch provides correctly rounded cu_fdiv_rn / cu_sqrt_rn (goal 07a;
  see mpblas/optimized/mpfr/cuda/README.md and external/mpc_cuda/patches).
  Real Rtrsm uses cu_fdiv_rn at fixed precision, and the panels can stay
  on cu_freal (goal 07 still measures the alternatives: cu_mpfr on the
  device, the panel on the host).
- Streams: a panel stream and an update stream with events, so that the
  look-ahead of Rlookahead_common.h maps onto the GPU (panel j+1 while the
  trailing update of step j runs). Scheduling only; arithmetic order as in
  mode 1.
- The future device-resident public API (matrices that stay on the GPU
  across calls, handles, streams) is NOT part of this series; keep the
  device-pointer internals clean enough to expose later.

## Verification gates (every goal that touches a routine)

- Bit identity (default mode): compare with the reference CPU routine
  (`<R>_reference_cpu` or a build of libmplapack_<b>) with memcmp of A,
  ipiv and info, on: square and rectangular (m > n, m < n), n in {1, 2, 3,
  nb-1, nb, nb+1, 2nb+3, odd sizes up to 600}, lda > m, singular matrices
  (zero column, exactly singular), matrices that need pivoting at every
  step, tie pivots, for potrf: SPD, not SPD at step j (info = j), both UPLO.
- Backward error (fast mode, and as a sanity check of the default):
  ‖P·A − L·U‖_1 / (n·‖A‖_1·eps) ≤ 30 and ‖A − U^H·U‖_1 / (n·‖A‖_1·eps)
  ≤ 30, eps of the backend (Rlamch("E")).
- LAPACK test suite: the drivers of mplapack/test/lin (common/Rchkge.cpp,
  Rchkpo.cpp, Rdrvge.cpp, Rdrvpo.cpp and the complex ones) must still pass
  for the backend.
- GPU device tests are skipped with exit code 77 without a device
  (CMake SKIP_RETURN_CODE 77); a host emulation variant
  (MPLAPACK_GEMM_CUDA_HOST_EMULATION / MPLAPACK_MPFR_CUDA_HOST_EMULATION)
  runs the same kernels on the host in CI without a GPU.
- Both build systems; misc/check_unique_symbols.sh on every changed
  archive; `nm -D` shows each declared entry point.

## Benchmarks (goal 01 creates them; later goals report against them)

benchmark/ programs for Rgetrf, Rpotrf (and Rgesv, Rposv from goal 08):
n = m in powers of two and multiples of 13 up to a size limit, best of 3
for n ≤ 2000, one run above; output n, seconds, GFLOPS with the LAPACK
counts (getrf: 2n³/3, potrf: n³/3, real; ×4 complex). Report for dd, qd and
mpfr 512 bits at least: reference, current optimized (look-ahead), the
goal's new path, and for GPU goals the time split (upload, factor,
download). Reference points from this machine (A100 80GB PCIe,
Threadripper 3970X 32 cores): MPFR 512-bit GPU Rgemm 3.6 GFLOPS at n=1024,
4.6 at 2048, 6.0 at 8192 (bench_local, not in the repository); dd GPU
Rgemm 200 GFLOPS at n=2048 (theoretical dd peak of the BF arithmetic on
A100: 336 GFLOPS).

## Environment notes (this repository, Ubuntu 26.04, CUDA 13.1)

- autotools: misc/reconfig.ubuntu26.04.nvidia.sh enables dd/td/qd/mpfr CUDA
  (`--with-mpfr-cuda-arch=80,90`). The CUDA device tests are defined only in
  CMake (`ctest -R cuda`); autotools `make check` runs the host emulation
  tests.
- nvcc with glibc 2.42 needs `-U_GNU_SOURCE -D_XOPEN_SOURCE=700` (already in
  the Makefile.am / CMake flags; keep it on new .cu files).
- dd needs -ffp-contract=off on every TU and `-fmad=false` on every .cu.
- The bundled mpc_cuda builds with LC_ALL=C (sort | comm in its build).
- Test binaries linked against the uninstalled libraries need the QD / MPFR
  library directories on the RPATH (`-Wl,--disable-new-dtags`), the
  libraries' own RUNPATH points at the install prefix.
