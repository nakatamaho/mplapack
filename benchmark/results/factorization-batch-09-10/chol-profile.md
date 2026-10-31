# Goal-10 initial profile: mixed result; implementation paused

Nsight Systems 2025.5.2, A100 80 GB PCIe, CUDA 13.1, resident Lower
Cholesky forced with `<P>_LAPACK_CUDA_MIN=0`. dd/qd n=2048, MPFR 512-bit
n=1024, one call per process. CUDA tracing only (CPU sampling unavailable).
These are diagnostic kernel sums, not isolated wall-time speedup benchmarks.
The profile was collected after item 0; resident arithmetic is still the base
`573645db3` implementation. Builds/CPU checks overlapped part of collection;
no other GPU computation ran concurrently with the six recorded profiles.

| Type | Real kernel time (%) | Complex kernel time (%) |
|---|---|---|
| dd | TRSM lanes 52.5; SYRK 28.0; GEMM 15.8; scalar base 3.6 | HERK 83.2; TRSM 9.6; GEMM 7.0; scalar base 0.3 |
| qd | TRSM lanes 47.2; SYRK 24.6; GEMM 22.2; scalar base 5.9 | HERK 76.4; TRSM 8.2; GEMM 15.1; scalar base 0.3 |
| mpfr 512 | panel/TRSM job 74.6; update job 25.4 | update job 84.4; panel/TRSM job 15.6 |

MPFR's single update kernel combines SYRK/HERK and GEMM, and its panel job
combines TRSM and scalar panel operations. No finer MPFR split is inferred
from those names. BF separates those kernel families. Totals are sums of
kernel durations, not critical-path times; streams can overlap.

The hypothesized update bottleneck is supported for **complex** Cholesky,
but is not the common dominant bottleneck claimed for real Cholesky. Real
TRSM/panel work dominates, especially MPFR, and right-looking updates alone
do not remove it. Under the instruction to stop goal 10 if the initial
hypothesis is not borne out, no right-looking code or eligibility changes
were made. Goal 09 proceeds independently. A future Cholesky session should
separate real TRSM/panel optimization from complex HERK parallelism, and
retain the stated per-element arithmetic and nonpositive-publication gates.

Raw kernel-summary CSVs are in `profiles/`. Reproduction:

```sh
OMP_NUM_THREADS=32 MPLAPACK_DD_LAPACK_CUDA_MIN=0 \
  /usr/local/cuda/bin/nsys profile --trace=cuda --sample=none \
  --cpuctxsw=none -o lower-dd \
  /tmp/factor04-cmake/benchmark/Rfactor.dd_cuda_total \
  --potrf --n 2048 --repeat 1
/usr/local/cuda/bin/nsys stats --report cuda_gpu_kern_sum --format csv lower-dd.nsys-rep
```

Add `--complex` for complex, select qd instead of dd, or use
`/tmp/factor01-cmake/benchmark/Rfactor.mpfr_cuda_total` at n=1024 with the
MPFR gate override. Trace binaries remain outside git. No goal-10 success
commit, bit-identity claim, or before/after speedup is asserted.
