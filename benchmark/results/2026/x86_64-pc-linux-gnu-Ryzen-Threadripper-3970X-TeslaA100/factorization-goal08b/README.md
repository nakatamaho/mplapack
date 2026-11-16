# Factorization goal 08b: CPU-library-speed fallback

Goal-08 base: `a86ee4fc5`. During the session HEAD independently advanced to
`f5bd68e28` (documentation only); that commit is preserved, not rewritten.
GPU arithmetic, transfers, precision/stack policies and public signatures are
unchanged. No generated reference routine was edited.

## Changes

* Rejected gesv/posv call the **public factor dispatcher + CPU solve**.
  Remove the private reference LU/Cholesky and serial reference GEMM from
  the driver fallback. getrs/potrs CPU fallbacks are unchanged.
* The CPU branches of dd/td/qd/MPFR factors and their recursive panels use
  CPU BLAS, including the optimized backend GEMM dispatcher. They use the
  CPU library's look-ahead policy, with no temporary environment overrides.
* Autotools' MPFR CUDA library, unlike CMake's, did not contain
  `Rgemm_blocked_omp`. Reuse that implementation privately (including CPU
  Winograd leaves); inlining/hidden constructor instantiation preserves the
  existing dynamic symbol set rather than adding a public GEMM dependency.
* The benchmark's `--legacy` comparator keeps private generated factors so
  it still measures per-call GPU BLAS after public CPU factors stop doing so.

## Eligibility: old -> new

Factor defaults (real/complex), common `<P>_LAPACK_CUDA_MIN` override unchanged:

| Backend | getrf real / complex | potrf real / complex |
|---|---|---|
| dd | 2048 / 2048 (unchanged) | 4096 / 4096 (unchanged) |
| td | 1024 / 1024 (unchanged) | 2048 / 2048 (unchanged) |
| qd | 1024 / 1024 (unchanged) | 2048 / 4096 (unchanged) |
| mpfr | 512 / 1024 (unchanged) | 2048 / **2048 -> 4096** |

MPFR Cpotrf L at 2048: CPU **69.6981932 s**, forced GPU **76.8520319 s**.
At 4096: CPU **490.306155 s**, forced GPU **315.619880 s**. Upper at 4096
is independently **810.859391 / 332.048643 s** (CPU/GPU), not inferred from
lower. After retuning, the 2048 default is **69.993494 s** versus CPU
**70.158712 s**; forced GPU remains **77.231652 s**.
Other factor gates remain conservative: resident calls at the gates beat
CPU; below-gate calls now use the CPU library's BLAS/schedule. Small warm
wins alone do not justify extrapolating a new universal cold crossover.

Fused-driver defaults (real/complex):

| Backend | gesv old -> new | posv old -> new |
|---|---|---|
| dd | 2048/2048 -> opt-in | 4096/4096 -> opt-in |
| td | 1024/1024 -> opt-in | 2048/2048 -> opt-in |
| qd | 1024/1024 -> opt-in | 2048/4096 -> opt-in |
| mpfr | 512/1024 unchanged | 2048/2048 -> opt-in |

Opt-in means an explicit common or driver-specific minimum enables the
**fused entry**; its default minimum is `LONG_MAX`. Ordinary drivers still
offload their **public factor** where eligible, then solve on the CPU.
Standalone solve controls/defaults are unchanged. No CPU fallback is added
inside a bool entry, and rejection still leaves every output untouched.

The goal-08 interleaved CPU/entry results already show why factor gates are
not driver gates: dd Rgesv 2048/rhs1 **0.754/1.730 s**, dd Rposv 4096/rhs1
**4.013/5.122 s**, td Rgesv 1024/rhs1 **0.504/0.820 s**, qd Rgesv 2048/rhs1
**3.705/4.818 s**, MPFR Cposv 2048/rhs1 **75.953/217.613 s**. See the retained
[goal-08 table](../factorization-goal08/solvers.csv). These fused engines and
CPU-library algorithms have not changed. Marginal lower-only wins do not
establish a safe n gate for all orientations/RHS; no unmeasured larger gate
is invented. MPFR LU retains its measured CPU-library crossovers.

## Measurements

A100 80 GB PCIe, Threadripper 3970X (32 cores / 64 threads), CUDA 13.1.115,
GCC 15.2.0, `OMP_NUM_THREADS=32`, `OMP_DYNAMIC=FALSE`, MPFR **512 bits**.
Shared-host load was approximately 72-77; absolute CPU times are not calibrated.
No correctness jobs, builds or other benchmarks ran concurrently with samples.

Each point uses three fresh-process rounds, reversing mode order on alternate
rounds, and reports their median. The benchmark's normal policy is retained:
best of three calls at n<=2000, one call above; initialization/input copies
are outside the timer. There is no silent substitution of historical timings
for the requested four-backend/four-size factor comparison.

* `factors-before-gates.csv`: **64 points / 192 medians / 576 retained results**,
  CPU library / CUDA default / forced CUDA, after CPU-binding fixes and before
  eligibility retuning. `factors-before-gates-samples.csv` retains every result.
* `factors.csv`: final factor comparison; only the changed MPFR Cpotrf 2048
  path is replaced by a fresh post-retuning three-mode comparison. The other
  63 points execute unchanged factor paths. Raw final results are in
  `factors-samples.csv`; the changed point is also in `retuned-factors*.csv`.
* `drivers.csv`, `drivers-samples.csv`: **16 points / 48 medians / 144 results**,
  CPU / parent executable / new fallback, rhs=1, n=1024 BF / n=256 MPFR.
  A high driver minimum forces fallback; public factor controls remain normal.
  Thus td/qd LU can use their resident public factors in the new composition.
* `upper-crossover*.csv`: extra MPFR Cpotrf U 4096 CPU/forced comparison.
* `driver-gates*.csv`: fresh CPU/default/forced comparisons for dd Rgesv 2048
  and Rposv 4096 after fused-entry defaults become opt-in:
  Rgesv **0.750172 / 0.669455 / 1.709880 s** and Rposv
  **3.815835 / 1.553601 / 5.100371 s** (CPU/default/forced).
* `gpu-phase-samples.csv`: **552** observed per-call upload/factor(or solve)/download
  intervals. For best-of-three processes these are **all** observed calls,
  not falsely paired with the selected minimum total. Intervals include host
  staging/publication where the entry does; they need not sum to total time.

Selected fallback driver medians, seconds (all 16 points retained above):

| Backend / call | n | CPU | parent fallback | new fallback |
|---|---:|---:|---:|---:|
| dd Rgesv | 1024 | 0.170147 | 3.150783 | 0.162293 |
| dd Rposv | 1024 | 0.683023 | 1.433871 | 0.642856 |
| td Rgesv | 1024 | 0.486015 | 10.256936 | 0.161708 |
| qd Rgesv | 1024 | 0.964976 | 21.030337 | 0.294031 |

The hypothesized equally large complex-driver regression was **not** observed:
dd Cgesv was 0.626 -> 0.517 s (CPU 0.531), Cposv 0.680 -> 0.688 s
(CPU 0.671). Complex tables are not omitted when there is no large win.

### CPU-only stability / benchmark synchronization artifact

A first timed `cudaDeviceSynchronize` in the CUDA benchmark can initialize
CUDA even after a CPU-only library call. This is **not** GPU BLAS activity
inside the factor. It contributes roughly 0.16 s to some single-call
below-gate totals (notably dd Rpotrf 2048). The benchmark was not silently
changed to hide it.

Seven interleaved rounds with explicit `LAPACK_CUDA=0` and repeat=3 isolate
the CPU path. `cpu-fallback-stability*.csv` retains all **56 results**:

| Call | n | CPU | CPU fallback | fallback / CPU |
|---|---:|---:|---:|---:|
| dd Rpotrf | 1024 | 0.628801 | 0.612936 | 0.9748 |
| dd Rpotrf | 2048 | 1.267107 | 1.288574 | 1.0169 |
| td Cpotrf | 1024 | 1.321018 | 1.356329 | 1.0267 |
| mpfr Rgesv | 256 | 0.614402 | 0.588817 | 0.9584 |

Initial three-round MPFR Rgesv 256 medians were CPU 0.478 / new 0.571 s;
those results are retained, not discarded. CPU samples ranged 0.403-1.002 s
across the initial/follow-up rounds; the seven-round comparison does not
support a persistent regression. The CPU binding/scheduling is the same
algorithm as the CPU library, not a slower serial substitute.

## Verification

Completed numerical checks use the original criteria: BF component bytes;
MPFR precision, exact value and zero sign; full A/B including padding and
unused triangles, pivots and info. The driver oracle is **public factor +
CPU solve**, not a private factor selected to mask a public-factor convention.
Dense pivoting, singular/nonpositive matrices, signed zeros, n=0/3/65/257,
rhs=0/1/3, N/T/C and U/L extend the original standalone zero probe.

| Gate | Result |
|---|---|
| CMake strict public fallbacks | 9/9 configurations; 352 cases / 14,149,260 component comparisons each; **0 bit / zero-sign errors** |
| Autotools-library strict public fallbacks | 9/9 configurations, same cases/comparisons, **0 errors** |
| Post visibility-only relink rechecks | MPFR six dynamic configurations + CMake 512, 7/7, **0 errors** |
| BF direct solvers host/device | 6/6 full suites, unchanged strict criteria |
| BF public factors default/fast | 6/6; bit identity/default and original backward-error/fast criteria |
| BF factor/BLAS device regressions | 18/18, all dd/td/qd families |
| MPFR public factors default/fast | 2/2; six-precision fallback configurations separately covered above |
| MPFR factor/solver device regressions | 512 bits, real/complex factors and eight solvers, automatic/runtime, 6/6 quick suites; **0 signed-zero errors** |
| MPFR BLAS device regressions | Rgemm, fixed complex BLAS, fixed-only Level 3: 3/3 |
| Autotools host checks | BF solver 3/3 full; MPFR factor 2/2 quick + solver quick pass at 512 bits |
| CUDA-library-linked GE/PO | 8/8 CPU-fallback runs (R/C x four backends), plus 6/6 forced BF runs; each run passes all four GE/PO routine/driver summaries |
| memcheck | BF quick + MPFR 12 required frame sweeps (six precisions x automatic/runtime), **0 sanitizer errors** |
| Final post-retuning public fallbacks | BF 3/3 + MPFR 512 1/1, **0 bit / zero-sign errors** |
| No visible GPU | CMake/autotools host emulation 8/8; actual CUDA-library CPU fallbacks dd/MPFR 512 2/2, **0 errors** |
| Final memcheck replay | MPFR 12/12 six-precision frame sweeps again, **0 sanitizer errors** |
| Builds / packaging manifests | Autotools and CMake PASS; both CMake source-manifest checks PASS; runner CLI smoke PASS |
| Archive / shared ABI | 16 CPU/accelerator archives PASS; all four full dynamic defined-symbol sets unchanged; 65 entries/library (8 solvers), 0 unresolved imports |

Legacy fixed BLAS regressions retain their existing zero-only convention:
Rgemm recorded **2**, complex BLAS **1096**, and fixed Level 3 **580**
signed-zero differences, with **0 failures**. These are not driver-composition
mismatches and no tolerance was widened. Repairing that legacy BLAS convention
(as done for real TRSM in goal 07a) remains separate work, not a reason to
restore a serial driver fallback.

Long unmodified GPU arithmetic sweeps were deliberately narrowed: redundant
BF host-factor tails (covered by full host/device fused solvers), full MPFR
host/runtime BLAS sweeps, and forced MPFR full GE/PO runs were cancelled and
**not counted as passes**. MPFR fallback coverage spans all six required
precisions; unchanged GPU factor/solve arithmetic is rechecked at 512 bits,
and all six precisions retain frame/stack memcheck coverage. No transfer,
stream-dependency, memory-budget, precision-restoration or device arithmetic
implementation changed, so delayed-copy negative controls and every old
failure-injection suite were not redundantly rerun.

**UNVERIFIED:** other GPU architectures, other precisions' performance
crossovers, unmeasured shapes/thread counts/sizes, and the full unchanged GPU
numerical sweeps beyond the selected 512-bit regressions. Eligibility is not
claimed to prove speedup on unsampled workloads.

## Reproduction

```sh
python3 benchmark/run_factor_dispatch_bench.py \
  --bin-dir /tmp/factor04-cmake/benchmark \
  --mpfr-bin-dir /tmp/factor01-cmake/benchmark \
  --sizes 512 1024 2048 4096 --rounds 3 --threads 32 --output factors.jsonl
python3 benchmark/run_factor_dispatch_bench.py \
  --bin-dir /tmp/factor01-cmake/benchmark --backends mpfr \
  --routines potrf --kinds complex --upper --sizes 4096 \
  --modes cpu forced --rounds 3 --threads 32 --output upper.jsonl
```

Use a fresh output file after changing a build/gate. No push is performed.
