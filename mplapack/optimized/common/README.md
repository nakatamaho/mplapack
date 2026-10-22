# Optimized LAPACK routines shared by every backend

`R*.cpp` / `C*.cpp` here replace the files of the same name in
`mplapack/reference` in every optimized library. They run the reference
code (included unchanged under another name) on independent parts of the
problem in parallel (`mpblas/optimized/common/Rblas_split_common.h`), so the
results are the reference results bit for bit.

| Routine | Parallel over |
|---|---|
| `Rlaswp`, `Claswp` | ranges of columns |

The other LAPACK routines (`Rgetrf`, `Rpotrf`, `Rgetrs`, `Rtrtri`, `Rlauum`,
`Rsygst`, ...) are the reference ones; they spend their time in the Level 3
BLAS (`Rgemm`, `Rtrsm`, `Rsyrk`, `Rtrmm`, `Cherk`, ...), which the optimized
libraries run in parallel (and with AVX2 for dd, td, qd) with the reference
results bit for bit, so these routines are faster with unchanged results.
