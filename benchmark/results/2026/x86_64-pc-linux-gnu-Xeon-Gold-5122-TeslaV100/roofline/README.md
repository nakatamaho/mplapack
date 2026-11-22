# Roofline of the dd/td/qd GEMM on the V100 (November 22, 2026)

Tesla V100-PCIE-32GB (sm_70). Method as the A100 roofline
(`../../x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/roofline`),
for the GEMM sweep only (no kernel-bench run on the V100 yet).

- `peaks.csv`: `benchmark/peak_bench.cu` built for sm_70 (`-fmad=false`),
  GPU alone. Only the dd/td/qd rows are used; the MPFR rows are the
  peak_bench values with operands from double, which overestimate the rate
  (see the A100 README), and are not used here.
- `roofline.csv`, `roofline_summary.md`, `roofline_{dd,td,qd}.png`:
  `benchmark/roofline.py peaks.csv <empty kernels.csv> ../gemm-sweep .`

| type | madd peak (/s) | GEMM bound 2·madd (GFLOPS) | GEMM NN n=2048 | % of roofline |
|---|---:|---:|---:|---:|
| dd | 1.138e11 | 227.7 | 231.6 | 101.7 |
| td | 3.566e10 | 71.3 | 71.1 | 99.7 |
| qd | 1.674e10 | 33.5 | 33.0 | 98.7 |

As on the A100 (99-101%), the GEMM runs at the madd rate of the
microbenchmark. Device memory (read + write): 0.47 TB/s.
