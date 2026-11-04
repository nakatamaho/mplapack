# Goal 09 CPU: implementation-stage checks

Dedicated R/C TRTRI and LAUUM use a shared recursive output-lane scheduler
and unchanged reference leaves. Default keeps iMlaenv blocking and each
full reduction; fast mode uses a whole-matrix binary recursion (32-wide
leaves). POTRI/GETRI retain the reference outer algorithms. CUDA flavors
compile the same CPU implementation with optimized CPU-BLAS bindings,
including the leaf routines, rather than private serial GEMM fallbacks.
Autotools source/shadow manifests and generated backend makefiles were updated;
CMake discovers/shadows these four common sources automatically.

## Checks actually completed before the CPU commit

| Check | Result |
|---|---|
| CMake libraries | dd/td/qd CPU and CUDA builds PASS |
| Autotools libraries | dd CPU/CUDA builds PASS; `gen_configure.sh`, `config.status` PASS |
| CPU inverse quick, CMake | dd/td/qd default: **360 checks each, 0 failures**; fast: **232 each, 0 failures** |
| CUDA-library CPU fallback inverse quick, CMake | dd default **361**, fast **233**, **0 failures** |
| CPU inverse quick, autotools | dd default **360**, fast **232**, **0 failures** |
| Source manifests | PASS |
| `check_unique_symbols.sh` | **8 archives PASS**: CMake dd/td/qd CPU/CUDA and autotools dd CPU/CUDA |
| `nm -D` | all four ordinary inverse symbols defined once in each autotools dd CPU/CUDA shared library |

Quick sizes are n=0/3/65, lda=n+3, real/complex, both triangles, unit/nonunit
TRTRI, zero diagonals, signed-zero off-diagonals, pivoted GETRI, workspace
queries. Default comparisons include complete A and GETRI work storage;
fast checks use product/inverse residuals. No comparison tolerance was
relaxed. CPU quick tests use OMP_NUM_THREADS=8; the standard CPU test policy
forces task look-ahead to exercise the recursive output tasks at small n.

**Not claimed:** GPU inverse entry points (the next item), MPFR/all-other-
backend build completion, full inverse size sweeps, `test/lin` inverse cases,
memcheck, or performance crossovers. Those are UNVERIFIED here. No inverse
benchmark matrix was run, in accordance with the maintainer's stage-1-only
amendment. QA は段階 2 でまとめて行う (a later session); no push.
