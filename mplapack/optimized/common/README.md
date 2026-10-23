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
| `MPLAPACK_LAPACK_FAST=1` | the recursive `Rgetrf2`/`Rpotrf2` (`Cgetrf2`/`Cpotrf2`) of the reference on the whole matrix; different rounding |

Measured on 4 cores (median of 5), against the reference loop on the same
parallel BLAS: `Rgetrf` double 1.6x (n = 1000), 1.23x (2000), 1.10x (3000),
dd 1.0x (600), 1.08x (1000); slower below n = 768 (0.65-0.74x), hence the
default. Cholesky with tasks: 0.66-1.17x, because the left-looking order of
the reference puts the diagonal blocks and their updates on a serial path,
which the reference loop runs on the parallel BLAS; it is therefore off by
default. The recursive variant measured 0.85-1.2x of the blocked routine.
Most of the speed of these routines comes from the BLAS below.

The other LAPACK routines (`Rgetrs`, `Rtrtri`, `Rlauum`,
`Rsygst`, ...) are the reference ones; they spend their time in the Level 3
BLAS (`Rgemm`, `Rtrsm`, `Rsyrk`, `Rtrmm`, `Cherk`, ...), which the optimized
libraries run in parallel (and with AVX2 for dd, td, qd) with the reference
results bit for bit, so these routines are faster with unchanged results.
