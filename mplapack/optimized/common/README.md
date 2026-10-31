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
| `MPLAPACK_LAPACK_LOOKAHEAD=1` | CPU default: LU tasks at min(m,n) >= 12 nb (768), or real LU >= 2 nb with more than eight requested threads; also the CPU-only fallback of dd/td/qd/MPFR CUDA factors |
| `MPLAPACK_LAPACK_LOOKAHEAD=2` | tasks also for smaller LU and for Cholesky |
| `MPLAPACK_LAPACK_LOOKAHEAD=0` | the reference loop |
| `MPLAPACK_LAPACK_FAST=1` | dedicated LU / Cholesky algorithms described below; different rounding |

Measured on 4 cores (median of 5), against the reference loop on the same
parallel BLAS: `Rgetrf` double 1.6x (n = 1000), 1.23x (2000), 1.10x (3000),
dd 1.0x (600), 1.08x (1000); slower below n = 768 (0.65-0.74x), hence the
original default. On a loaded 32-core host, however, repeatedly launching
32-thread BLAS teams below that boundary costs far more than the arithmetic:
dd LU n=512 takes about 1.4 s, versus 0.15 s at n=1024. With more than eight
requested threads, small real LU now shares a single task team, bounded by the
number of column tiles (or 1024-row panel chunks for tall matrices). The
n=512 time becomes about 0.034 s with unchanged arithmetic. No process-wide
environment or OpenMP ICV is changed; explicit LOOKAHEAD=0 remains serial
outer scheduling. `Rfactor_test --scheduling` covers n=256/384/512/768/1024
in production scheduling mode, including real/complex LU and both Cholesky
triangles. Complex LU did not exhibit the same extreme slowdown and retains
its existing default gate; Cholesky scheduling itself is unchanged.

Cholesky with tasks: 0.66-1.17x, because the left-looking order of
the reference puts the diagonal blocks and their updates on a serial path,
which the reference loop runs on the parallel BLAS; it is therefore off by
default. The old whole-matrix reference recursion measured 0.85-1.2x of the blocked routine.
Most of the speed of these routines comes from the BLAS below.

## Dedicated CPU LU

`Rgetrf2.cpp` / `Cgetrf2.cpp` retain the reference binary recursion and
safe-minimum scaling test. Pivot search uses contiguous row blocks with a
strict-greater, increasing-index reduction (complex magnitude is abs(real)
+ abs(imag), not Euclidean magnitude). Scaling uses exactly the reference
reciprocal, or elementwise division below the safe minimum. Row blocks
of 1024 become tasks only at 4096 rows; smaller columns stay serial.
Ordinary panels (at most 64 columns and fewer than 4096 rows) reuse the
unchanged reference recursion as a serial leaf, avoiding scheduler overhead.
Update coefficients retain the reference `-one` representation, including
signed-zero complex components and expansion limbs; constructing a fresh
`COMPLEX(-1)` is not interchangeable for bit identity. Regression tests cover
signed zeros in wide panels and in the serial blocked dispatcher.
The pivot interchange precedes scaling; swaps across independent columns
use the existing optimized LASWP. NaN block candidates trigger the serial
pivot search to retain its comparison behavior.

Within an existing LU team, recursive updates fuse LASWP/TRSM/GEMM on
independent column tiles with local taskgroups. No reduction dimension is
split. The default outer loop keeps `iMlaenv` and the look-ahead policy above;
standalone tall panels can create their own team. Workers inherit the
caller's MPFR/GMP precision, rounding, and exponent range. LOOKAHEAD=0
also disables panel tasks. The independent `*_reference_cpu` routines call
unchanged, separately named reference panels, not the new implementation.

Fast LU uses its own panel width (`MPLAPACK_GETRF_TILE`, default 32,
clamped to at least 1), recursive panels and dependency-driven look-ahead
when the normal task gate permits it. Each trailing tile fuses swaps,
TRSM and GEMM. The optimized GEMM exposes no public packing handle, so
packing remains inside the BLAS engine rather than being duplicated here.
The different blocking may change rounding; backward error <= 30 and
non-tie pivot identity are required. `Rgetrs` / `Cgetrs` and `Rgesv` /
`Cgesv` remain unchanged and resolve the new public LU through normal linking.

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

**Accelerator scheduling:** both build systems define only
`MPLAPACK_ACCELERATOR_LIBRARY` for accelerator targets. The dd/td/qd/MPFR
CUDA factor dispatchers now use CPU BLAS throughout their rejected-entry
path, including recursive panels (`Rfactor_cpu_blas.h`, goal 08b). Their
CPU scheduling therefore matches the CPU library, without GPU transfers or
process-wide environment overrides. Other accelerator flavors still require
explicit `MPLAPACK_LAPACK_LOOKAHEAD=2`: transfer/serialization overhead makes
per-call GPU BLAS tasks unsuitable as their default schedule. CPU libraries
and resident GPU arithmetic are unchanged.

The CUDA libraries additionally provide [resident real/complex solvers](Rsolver_cuda.md)
for dd/td/qd/mpfr, and [resident real/complex inverses](Rinverse_cuda.md) for
dd/td/qd (opt-in until crossovers are measured); the CPU solvers are unchanged.

## Dedicated CPU triangular inverses and products

`Rtrtri` / `Ctrtri` and `Rlauum` / `Clauum` use `Rinverse_cpu_common.h`.
Default mode keeps the exact reference block width, diagonal TRTI2/LAUU2
leaves, singularity scan, flags, scalar representations and full reductions.
Recursive tasks split only independent output lanes: left TRMM columns,
right TRSM rows, and the independent rows/columns of the LAUUM updates.
Workers inherit the caller's arithmetic context. No triangular dependency
or GEMM/SYRK/HERK reduction axis is split. The normal LAPACK look-ahead
controls select the shared team; LOOKAHEAD=0 uses CPU BLAS directly.

`MPLAPACK_LAPACK_FAST=1` selects a whole-matrix binary recursion with 32-wide
reference leaves. This may round differently and is checked by inverse or
product residuals, not bit identity. Unit diagonal storage, the opposite
triangle and padding remain untouched. A non-unit singular TRTRI returns the
first zero diagonal index before changing A in either mode.

POTRI and GETRI retain their generated reference outer algorithms and call
the public optimized TRTRI/LAUUM routines. CUDA-library CPU branches use
`Rinverse_cpu_blas.h`, extending the factor CPU-BLAS binding to TRMM and the
leaf BLAS too; they do not route fallback arithmetic through GPU BLAS.
`Rfactor_test --inverse-quick` exercises all four inverse families in real
and complex arithmetic; `--inverses` extends the sizes through 599. Default
tests compare the complete storage and GETRI workspace, with zero-sign-sensitive
component checks; fast tests check products/inverse residuals.

The other CPU LAPACK routines (`Rgetrs`, `Rsygst`, ...) are the reference ones; they spend their time in the Level 3
BLAS (`Rgemm`, `Rtrsm`, `Rsyrk`, `Rtrmm`, `Cherk`, ...), which the optimized
libraries run in parallel (and with AVX2 for dd, td, qd) with the reference
results bit for bit, so these routines are faster with unchanged results.
