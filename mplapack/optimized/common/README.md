# Optimized LAPACK routines shared by every backend

`R*.cpp` / `C*.cpp` here replace the files of the same name in
`mplapack/reference` in every optimized library. They run the reference
code (included unchanged under another name) on independent parts of the
problem in parallel (`mpblas/optimized/common/Rblas_split_common.h`), so the
results are the reference results bit for bit.

| Routine | Parallel over |
|---|---|
| `Rlaswp`, `Claswp` | ranges of columns |
| `Rgetrf`, `Cgetrf`, `Rpotrf`, `Cpotrf` | tasks on column (row) tiles, below |

`Rgetrf`, `Cgetrf`, `Rpotrf`, `Cpotrf` (`Rlookahead_common.h`): the
reference blocked loop (same block size, `Rgetrf2`/`Rpotrf2` on the panels and
diagonal blocks, then `Rlaswp`, `Rtrsm`, `Rgemm`, `Rsyrk`) split into OpenMP
tasks on tiles of nb columns (rows), with dependencies on the tiles each task
reads and writes, so that the next panel is factored while the other tiles
are being updated (look-ahead). Every tile receives the operations of the
reference in their order, so the results, info and pivots included, are the
reference ones bit for bit (`mpblas/optimized/common/test/Rlapack_omp_test.cpp`).

| Variable | Meaning |
|---|---|
| `MPLAPACK_LAPACK_LOOKAHEAD=1` | default: tasks for LU with min(m,n) >= 12 nb (768) |
| `MPLAPACK_LAPACK_LOOKAHEAD=2` | tasks also for smaller LU and for Cholesky |
| `MPLAPACK_LAPACK_LOOKAHEAD=0` | the reference loop |
| `MPLAPACK_LAPACK_FAST=1` | LU: reference whole-matrix recursion; Cholesky: dedicated whole-matrix recursion described below; different rounding |

Measured on 4 cores (median of 5), against the reference loop on the same
parallel BLAS: `Rgetrf` double 1.6x (n = 1000), 1.23x (2000), 1.10x (3000),
dd 1.0x (600), 1.08x (1000); slower below n = 768 (0.65-0.74x), hence the
default. Cholesky with tasks: 0.66-1.17x, because the left-looking order of
the reference puts the diagonal blocks and their updates on a serial path,
which the reference loop runs on the parallel BLAS; it is therefore off by
default. The old whole-matrix reference recursion measured 0.85-1.2x of the blocked routine.
Most of the speed of these routines comes from the BLAS below.

## Dedicated CPU Cholesky

`Rpotrf.cpp` and `Cpotrf.cpp` reuse the unchanged generated blocked loop
under a private name, separately from the `*_reference_cpu` oracle. Default
mode retains its `iMlaenv` block size, binary `potrf2` recursion, full-length
BLAS reductions, failure stopping point, and untouched triangle/padding.
There is no split of the reduction dimension and no replacement of HERK by
GEMM. `Rpotrf_common.h` schedules independent triangular-solve tiles and
SYRK/HERK output-column tiles in the recursive diagonal block. Its local
taskgroups wait only for that recursion node, not for the whole team; the
existing look-ahead releases panel GEMM/TRSM tiles immediately when the
diagonal is ready. The outer team installs the caller's MPFR/GMP precision
and rounding context on every worker.

As before, Cholesky task look-ahead requires `MPLAPACK_LAPACK_LOOKAHEAD=2`;
mode 0/1, a single thread, or a call from an existing parallel region uses
the blocked schedule with the parallel BLAS. Small blocks do not create
tasks. Keeping task scheduling opt-in avoids imposing its overhead on
small matrices and preserves the established production schedule.

`MPLAPACK_LAPACK_FAST=1` instead recursively factors the whole matrix,
solves the off-diagonal half, updates the trailing half with the optimized
SYRK/HERK, then recurses. The split is aligned to a separate leaf size
(not `iMlaenv`); leaves use the unchanged reference `potrf2`. These BLAS
calls use their own CPU parallelism, not concurrent factorization tasks.
Fast mode may round differently; its correctness gate is backward error
<= 30, not bit identity. Both modes support both triangles and propagate
`info` at the first nonpositive pivot.

| Variable | Meaning (read on every call) |
|---|---|
| `MPLAPACK_POTRF_TILE` | task tile width; default 16 for fixed-precision backends, 8 for GMP/MPFR; positive integer, clamped to at least 1 |
| `MPLAPACK_POTRF_LEAF` | fast recursion leaf size/alignment; default 64 for fixed-precision backends, 32 for GMP/MPFR; clamped to at least 1 |

**Accelerator safety:** only CPU optimized targets define
`MPLAPACK_CHOLESKY_CPU_TASKS`. Autotools grants this in the generated CPU
Makefiles; CMake grants it only to the two Cholesky sources of CPU targets.
CUDA (and OpenCL) libraries do not define it: even `LOOKAHEAD=2` cannot
activate Cholesky tasks there, in either result mode. This is a compile-time
permission, not an environment-based guess at the linked BLAS. It is
necessary while the independent CUDA-concurrency fix is unmerged. It does
not fix LU look-ahead or MPFR CUDA's illegal-memory-access defect, nor make
simultaneous external calls safe.

The other LAPACK routines (`Rgetrs`, `Rtrtri`, `Rlauum`,
`Rsygst`, ...) are the reference ones; they spend their time in the Level 3
BLAS (`Rgemm`, `Rtrsm`, `Rsyrk`, `Rtrmm`, `Cherk`, ...), which the optimized
libraries run in parallel (and with AVX2 for dd, td, qd) with the reference
results bit for bit, so these routines are faster with unchanged results.
