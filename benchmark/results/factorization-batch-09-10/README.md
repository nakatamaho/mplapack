# Batch 09/10 — implementation-stage receipts

Base: `573645db3`, branch `claude/mpfr-cuda-rgemm`. The maintainer subsequently
requested **stage 1 only**: no aggregate stage-2 QA, stage-3 benchmark matrix,
gate retuning, or push. The measurements below diagnose item 0; they are not
crossover measurements or a substitute for deferred QA.

Stage-1 status: CPU LU follow-up and goal-09 CPU/GPU implementations are
committed separately. Goal 10 was **paused after its mixed initial profile**,
not reported as an implementation success. Receipts:
[Cholesky initial profile](chol-profile.md),
[goal-09 CPU](cpu-inverses.md), [goal-09 GPU](gpu-inverses.md).
Stage 2 and stage 3 were not run; existing gates were not retuned.

## Item 0: small CPU LU

Reproduced without CUDA: with 32 threads, dd Rgetrf n=256/384/512 took
0.956/1.634/1.388 s, whereas n=768/1024 took 0.091/0.145 s (best of two).
At n=512, eight threads took 0.0286 s, one thread 0.0884 s. Passive OpenMP
waiting reduced the 32-thread case to 0.161 s; forced task look-ahead reduced
it to 0.0446 s. CUDA initialization cannot explain this CPU-only discontinuity.

The default outer schedule changes at 12*nb=768: below that boundary it
repeatedly starts small parallel BLAS calls, including recursive-panel calls,
with the entire requested team. On this loaded host, active waits and team
barriers overwhelm their arithmetic. There is no nested active parallel team
in this default path. The fix shares one bounded task team for **real** LU
at >=2*nb when more than eight threads were requested. Small task teams are
bounded by column tiles / tall-panel row chunks. No arithmetic, nb, environment
variable, OpenMP ICV, or Cholesky default policy changes.

The initial candidate applied the small gate to complex LU as well. The
360 diagnostic points in `cpu-lu-diagnosis.csv` retain that **candidate**
(`path=after`), not a falsely labelled final complex implementation. All
three backends, R/C getrf/potrf, n=256/384/512/768/1024, and 1/8/32 threads
are retained. Complex LU did not show the extreme real-LU discontinuity;
its default gate was restored. Cholesky was never changed. Single-process
best-of-two diagnostics are not statistically calibrated medians; host load
was approximately 65. An existing GPU quick test briefly overlapped the
initial diagnosis; these are diagnostic, not isolated benchmark receipts.

`cpu-lu-final-smoke.csv` is a separate final-code check (32 threads):

| Backend | Rgetrf 512, old (s) | Rgetrf 512, final (s) | Rgetrf 1024, final (s) |
|---|---:|---:|---:|
| dd | 1.27455 | 0.03513 | 0.15157 |
| td | 1.60328 | 0.12069 | 0.45075 |
| qd | 2.23714 | 0.23261 | 0.86854 |

CUDA-library CPU fallback, dd, disabled resident entries, best of two:
512 **0.03370 s**, 1024 **0.14360 s**. This is the same shared CPU fix,
not a serial reference substitute. Final Cgetrf smoke points are retained
in the CSV too; no complex-LU speedup is claimed.

### Item-0 precommit checks actually run

* CMake: dd CPU/CUDA libraries, td/qd CPU libraries, their factor benchmarks,
  and dd CPU/CUDA factor tests built successfully.
* Autotools: dd CPU and CUDA libraries built successfully; CPU factor test built.
* CMake dd CPU production-default factor test: **2300 checks, 0 failures**.
* Added `Rfactor_test --scheduling`: dd CPU **100 checks**, CUDA-library CPU
  fallback **101 checks**, **0 failures**; n=256/384/512/768/1024, R/C LU,
  both R/C Cholesky triangles, exact default-mode comparison and residuals.
* dd solver host/device quick tests: **1600 checks / 1,538,784 component
  comparisons each**, **0 signed-zero errors**. These use unchanged GPU engines.
* Six rebuilt archives: `check_unique_symbols.sh` PASS (CMake dd CPU/CUDA,
  td/qd CPU, autotools dd CPU/CUDA).
* `nm -D`: autotools dd CUDA's complete defined dynamic-symbol set is unchanged;
  CPU Rgetrf/Cgetrf definitions also inspected.
* One full, unchanged dd GPU-getrf sweep was attempted but cancelled at the
  tool timeout. It is **not** counted as a pass. The full host-getrf sweep
  completed, but does not replace deferred device QA.

**Deferred / UNVERIFIED:** aggregate QA, remaining backend builds and full
numerical suites, sanitizer checks, calibrated performance comparisons,
threshold changes, other architectures. QA は段階 2 でまとめて行う (in a later
session under the maintainer's stage-1-only amendment).
