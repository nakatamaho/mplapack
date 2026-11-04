# Maintainer-requested pause, 2026-11-02

**Follow-up closed 2026-11-04:** resumed measurements were subsequently stopped
at the maintainer's request; [completed data and final gate decisions](../2026/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/gate-retune/README.md)
are archived. The statements below about unstarted benchmarks/unchanged gates
are historical, not the final status.

**Ended by explicit instruction on 2026-11-03.** HEAD is now `93d1a70fb`,
which reverts the small-n team cap. All benchmark executables were rebuilt
in both CMake trees and autotools before new GPU measurements. The CPU
comparison is not repeated; its timing concern is historical and was
superseded by the revert. Current work uses `/tmp/batch0910-retune-93d1/`,
not the old `phase3.py`. The remainder of this file is the preserved pause
record, not a current instruction to remain idle.

GPU benchmarks have **not started**. No CUDA LAPACK gate has been changed.
The retuning commit requested for stage 3 has not been made.

At 18:24 JST the active CPU comparison and its automatic continuation were
terminated. There are no remaining comparison, CUDA benchmark or sanitizer
processes from this session. Resume GPU computation **only after an explicit
maintainer instruction**; the maintainer is using the GPU for OCR.

Only `nvidia-smi`/NVML process-list recording remains active (no CUDA/OpenCL
computation). Observations during this hold are expected maintainer OCR
usage, **not measurement anomalies or contaminated benchmark samples**.
There are no GPU benchmark samples to invalidate.

## CPU comparison scope amendment

The maintainer accepts n=256/512/1024, threads=8/32 for the CPU-LU follow-up;
in particular all six required MPFR points are complete. Across all eight
CPU backends this is **48 points / 288 paired R/C comparisons**, with zero
A/pivot/info identity failures.

The original request for n>=2048 and one thread has been withdrawn:
**修正の対象外の経路なので省略**. The in-flight MPFR 4096/1 point was
cancelled; MPFR 4096/8 and /32 were never started. See `cpu-lu-omitted.csv`.
GMP 4096 and MPFR 2048 had already completed before the amendment arrived;
their existing data are retained as extra references, not rerun or required
acceptance points. No completed samples were deleted or relabelled as omitted.

`cpu-lu-comparison-samples.csv` retains **522 completed paired comparisons**:
288 required, 234 extra. `cpu-lu-comparison.csv` contains their medians and
scope labels. The 69 extra points are not prerequisites for GPU benchmarking.

### Timing interpretation / pending CPU performance concern

Exact-output checks passed; this does **not** automatically mean a timing
nonregression passed. At MPFR 512 bits, Rgetrf n=512/32 threads has medians
**3.17298448 -> 3.45718416 s** (1.08957x). This is an apparent ~9% slowdown
candidate, not a resolved performance PASS. It needs CPU-only follow-up before
claiming all-backend nonregression; no further CPU measurements or numerical
implementation changes were made after the hold instruction. The 8-thread control is
1.95505395 -> 1.96420628 s. All original samples are retained.

For comparison, dd Rgetrf 512/32 improves **1.51390501 -> 0.036979789 s**;
its final 1024/32 is **0.171058439 s**. See the CSV for exact retained values.

Method: reconstructed pre-cdc38b323 Rgetrf scheduling, with the other CPU
callees held at the final numerical implementation, compared in one process
against the current CPU library. Cgetrf is unchanged and is measured twice
as a noise/control comparison using the same function. Inputs are seeded
16-bit dyadics; comparisons include the full padded A, pivots and info using
the existing backend representation checks. Three alternating-order pairs
were used through n=1024; already completed larger extras used one pair.
These are shared-host timings, not calibrated CPU-core isolation.

## Resume state

Local state:

* `/tmp/batch0910-qa/GPU_HOLD`: explicit no-GPU guard.
* `/tmp/batch0910-qa/phase3.py`: stopped; updated to exclude withdrawn CPU
  sizes/thread counts and to exit when the hold guard exists.
* `/tmp/batch0910-measurements/`: original CPU samples and completion markers;
  no GPU benchmark JSONL exists.
* `/tmp/batch0910-measurements/gpu-hold-observations.jsonl`: OCR-period
  process-list records; their policy is record-only, not anomaly detection.
* `/tmp/batch0910-qa/gpu-hold-observer.py`: NVML recorder, exits when the
  guard is removed. It does not remove the guard or start any benchmark.

Do not remove the guard, launch CUDA tests/benchmarks, or infer permission
from an idle GPU. Await explicit permission. After permission, review the
CPU timing concern, stop/coordinate the NVML recorder as appropriate, and
continue the interleaved CPU/default/forced cold/warm GPU matrix. The
measurement runner and benchmark inverse/sample extensions remain uncommitted
preparation; no retuning or speedup claim follows from building them.
