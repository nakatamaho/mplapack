# Goal 09 GPU — implementation-stage receipts

Added 8 bool entries per backend (dd/td/qd): real/complex TRTRI, LAUUM,
POTRI, GETRI, with ordinary basename-shadowing dispatchers. The common device
engine retains reference blocking and reductions and reuses/generated-adapts
unchanged reference arithmetic. GETRI includes reduced-workspace paths and
preserves untouched work slots. Unit-diagonal storage is neither read nor
published. Rejection is transactional. CPU fallbacks use the same optimized
CPU sources and CPU BLAS, including GETRI/POTRI's outer composition; they do
not re-enter per-call GPU BLAS or a serial reference GEMM.

All new inverse entries are **opt-in until measured** (LONG_MAX default;
common or inverse-specific minimum overrides). Existing factor/solver gates
were not retuned. MPFR GPU inverses are out of scope. API/controls and the
MPFR follow-up direction: [Rinverse_cuda.md](../../../mplapack/optimized/common/Rinverse_cuda.md).

## Final precommit checks actually completed

| Check | Result |
|---|---|
| CMake libraries | dd/td/qd CUDA + their CPU dependencies built PASS |
| Autotools libraries | dd/td/qd CPU + CUDA built PASS |
| CMake direct device quick, default | dd/td/qd: **694 checks + 62 rejection cases each**, **0 failures**, 168 accepted numerical entries each |
| CMake direct host quick, default | dd: same **694 + 62**, **0 failures** |
| CMake direct device/host fast quick | dd: **462 checks + 62 rejection cases each**, **0 failures**, 168 accepted entries each |
| Autotools direct host quick | dd default **694 + 62**, fast **462 + 62**, **0 failures** |
| CMake ordinary public CUDA-library inverse quick | dd disabled resident entry / forced GPU x default / fast: **695 / 463 checks per configuration**, **0 failures** |
| Autotools actual shared-library inverse quick | dd disabled-entry / forced GPU, default: **695 checks each**, **0 failures**; all 8 entry families traced in forced mode |
| `check_unique_symbols.sh` | **12 archives PASS**: CMake + autotools x dd/td/qd x CPU/CUDA |
| `nm -D` | **8/8 new entries defined** in each dd/td/qd autotools CUDA shared library |
| `ldd -r` | all three CUDA shared libraries: **0 missing libraries / unresolved symbols**, with QD/CUDA runtime library paths |
| Source/generated manifests; whitespace | PASS; new CMake quick test registrations inspected (not an aggregate CTest run) |

Quick sizes: 0/3/65/131, real/complex, U/L, unit/nonunit, singular
TRTRI/POTRI/GETRI, signed-zero off-diagonals, lda=n+3, pivoted GETRI,
workspace queries and 1/2/optimal workspace columns. Default comparisons are
complete-storage, component-byte comparisons (including zero signs and
untouched A/work), not numerical tolerances. Fast checks use inverse/product
residuals. No threshold/tolerance was relaxed to obtain these results.

The direct harness requires every eligible n>0 numeric entry to return true;
all runs above recorded **0 unexpected rejections**. n=0/query numerical
wrapper cases use CPU no-computation semantics; separate direct contracts
require those entries to reject without publication. Per direct suite, the
62 rejection cases cover disabled/threshold/invalid arguments, null pointers,
invalid pivots, quick returns/queries, and size overflow. The public-library
runs additionally exercise actual dispatch/CPU-fallback composition.

Precommit build issues were corrected before these final runs: autotools
CUDA translation units need explicit backend defines (not merely host
compiler flags), and new host oracles need rebuilt parent CPU archives.
These failed intermediate builds are not counted as passes. Device generator
integer MOD calls were adapted to device-safe C++ remainder; no FABLE body
was edited. All verification above was rerun on the final relevant paths.

## Deferred / UNVERIFIED

No aggregate stage-2 QA or stage-3 benchmark/gate tuning was run, following
the maintainer's amendment. Full-size numerical sweeps, test/lin inverse
and GE/PO drivers, memcheck, hidden-GPU runs, MPFR/other-backend build and
numerical coverage, cold/warm timings and CPU/GPU crossovers remain deferred.
Only A100/sm_80 executed; sm_90 was compiled by autotools, not executed.
QA は段階 2 でまとめて行う, in a later session. No push.
