# GOAL 10: optimize the GPU kernels under getrf / potrf (stage 3 of the kernel plan)

Prerequisites: goal 11 (all four commits, including its kernel benchmark), the
roofline analysis 89019cc3c (benchmark/peak_bench.cu, benchmark/roofline.py,
results/2026/<host>/roofline/). Read 00-overview.md ("Result modes", "GPU
design") first.

## Why

The roofline (results/2026/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/roofline/README.md)
shows that the dd/td/qd GEMM, its k = 64 update form and SYRK Upper/Trans are
at 93-101% of the arithmetic peak, while the kernels that getrf and potrf also
need are far from it (n = 2048): TRSM 1-7%, the LU panel (getrf2, n x 64)
0.7-3%, the Cholesky diagonal block (potrf2) 0.1-0.3%, SYRK Lower/NoTrans about
50%, complex GEMM about 60%, complex HERK 2-5%, MPFR GEMM 14-46%. The GPU
Cholesky is slower than the GPU LU for that reason, and the complex Cholesky
entry is opt-in.

## Scope (in this order; one commit each)

1. **TRSM** (real and complex, dd/td/qd/mpfr; the forms used by getrf, potrf,
   getrs, potrs, trtri and their k = 64 panel forms). Today one GPU thread
   solves one column (Left) or one row (Right) of B, so only m or n threads run
   and each runs a long sequential loop. Parallelize over the elements of B as
   well: in the reference loops (mpblas/reference/Rtrsm.cpp, Ctrsm.cpp), for a
   fixed step k the updates of the remaining elements of a column (Left,
   axpy forms) or of a row (Right) are independent, so a block of threads can
   own a set of columns (rows), loop over k, update the independent elements in
   parallel and synchronize between steps. Every element keeps the reference
   sequence of operations (and its zero-skips and reciprocal multiplications),
   so the result stays bit-identical. Where the reference forms a dot product
   per element (Left Upper Trans / Left Lower Trans), keep its order and
   parallelize across the independent columns and the products only.
2. **SYRK Lower/NoTrans and HERK** (real and complex). First find with a
   profiler why Lower/NoTrans runs at half the speed of Upper/Trans with the
   same tiled kernel (k_syrk_tiled) (memory access pattern, the reference's
   zero-skip branch, the triangle mask); then bring both to the GEMM level.
   HERK: the same, plus the real diagonal of the reference.
3. **Panels**: the LU panel (getrf2: iamax, swap, scaling, the inner updates)
   and the Cholesky diagonal block (potrf2), real and complex. Use the
   parallel TRSM of item 1 and the GEMM inside the recursion, avoid one kernel
   launch per column where a single block (or a cooperative launch) can run
   the column loop on the device, and keep the reference order of every
   element and the pivot choice (first index on ties).
4. **Cholesky Lower right-looking** (bench_local/pi_instructions/goal10.md has
   the full argument; summary): the reference Lower potrf uses Rsyrk("L","N")
   / Rgemm("N","T") (Cherk("L","N") / Cgemm("N","C")) that add the terms of
   l = 1..k one by one into C, so applying the contributions of each block
   column to the whole trailing matrix (right-looking) gives every element the
   same sequence of operations. Per element: the zero-skip rule of the routine
   that updates it in the reference (SYRK/HERK for the diagonal block, GEMM
   below it) and Cherk's real diagonal. On info > 0 the output must equal the
   reference (the reference leaves the columns to the right of the failing
   block, and the part of the failing block column below its diagonal block,
   untouched). Upper stays left-looking (its updates are dot products).
5. **MPFR GEMM** (cu_freal 256/512/1024): find why it reaches only 14-46% of
   the madd peak (shared-memory tiling of large elements, register pressure,
   occupancy) and improve it; bit identity unchanged.

Out of scope: gates and dispatcher defaults (measure only; the reconnection
is a later goal), the fast modes, CPU kernels, the known zero-sign allowance
of the legacy fixed-precision MPFR BLAS.

## Acceptance criteria (per item)

- Bit-identical to the reference in default mode (zero signs, info, ipiv),
  host emulation and GPU, all forms of the routine touched, n = 0, 1, block
  boundaries, lda > n; the existing GPU BLAS / factorization / solver / inverse
  tests that use the kernel still pass. compute-sanitizer memcheck on one size.
  The MPFR runtime-precision device and host solver tests at 768, 1024 and
  2048 bits (mpfr_cuda_solver_{device,host}_runtime_*; one to two hours each)
  are part of "existing tests" and must run to completion: ba73551e8 (item 2)
  broke mpfr_cuda_solver_device_runtime_{768,1024,2048} (illegal memory
  access in the runtime complex HERK lanes) and the item's validation did not
  run them; fixed under goal 12, addendum 3. (Policy since November 8: host-only tests in
  parallel; the runtime 2048-bit tests once at the end of a goal.)
- Performance, measured before and after with the kernel benchmark
  (benchmark/Rkernel_bench.cpp / kernel_bench.py, square sizes n = multiples of
  7 and 31 and powers of two, GPU kernel only) and the roofline
  (benchmark/roofline.py with results/.../roofline/peaks.csv): report the % of
  roofline before and after at n = 512, 1024, 2048 (MPFR 1024: 512, 1024).
- After items 1-4: getrf and potrf (real and complex, dd/td/qd/mpfr 512,
  Lower and Upper, n = 1024, 2048, 4096) through the GPU entry points, before
  and after (benchmark/Rfactor.cpp --entry), and the CPU library for reference.
- nm -D unchanged (no new public symbols), check_unique_symbols.sh, both build
  systems.
- Commits: "GPU TRSM: parallel over the elements of B, bit-identical (factorization goal 10)",
  "GPU SYRK/HERK: Lower/NoTrans at the GEMM level (factorization goal 10)",
  "GPU getrf2/potrf2 panels on the parallel TRSM/GEMM (factorization goal 10)",
  "GPU potrf: right-looking lower Cholesky, bit-identical (factorization goal 10)",
  "MPFR CUDA GEMM: closer to the cu_freal madd peak (factorization goal 10)".
  If an item cannot reach bit identity or stalls, leave it uncommitted, record
  why, and continue with the next item. A session instruction from the
  maintainer may give more detail; follow it where it does not contradict this file.
