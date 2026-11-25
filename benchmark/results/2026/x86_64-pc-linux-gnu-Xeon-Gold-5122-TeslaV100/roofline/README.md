# Roofline of the dd/td/qd GEMM on the V100 (November 22, 2026)

Tesla V100-PCIE-32GB (sm_70). The roofline compares the kernel-only GEMM
sweep with the arithmetic bound of two FLOPs per measured madd. No
kernel-bench run was included.

- `peaks.csv`: `benchmark/peak_bench.cu` built for sm_70 (`-fmad=false`),
  GPU alone. Only the dd/td/qd rows are used. The MPFR rows overestimate
  the rate because the operands come from double, so they are excluded.
- `roofline.csv`, `roofline_summary.md`, `roofline_{dd,td,qd}.png`:
  `benchmark/roofline.py peaks.csv <empty kernels.csv> ../gemm-sweep .`

| type | madd peak (/s) | GEMM bound 2·madd (GFLOPS) | GEMM NN n=2048 | % of roofline |
|---|---:|---:|---:|---:|
| dd | 1.138e11 | 227.7 | 231.6 | 101.7 |
| td | 3.566e10 | 71.3 | 71.1 | 99.7 |
| qd | 1.674e10 | 33.5 | 33.0 | 98.7 |

GEMM throughput is close to the arithmetic bound from the madd
microbenchmark. Device memory (read + write): 0.47 TB/s.
