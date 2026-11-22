# dd/td/qd Rgemm on the V100 (November 22, 2026)

Tesla V100-PCIE-32GB (sm_70, 1380 MHz, CUDA 12.6 of HPC SDK 24.11), commit
9d8716b9f, GPU 0 alone (`gpu-processes-at-*.txt`). Kernel only:
`benchmark/run_gemm_sweep.sh` (data on the GPU, CUDA events, default tiled
kernel), the sizes of the MPFR 512 sweep (`../mpfr512-sweep`) and of the A100
`../../x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/gemm-sweep`:
multiples of 7 and 31 and powers of two up to 4096, best of 3 after a warm-up
for n <= 2000, one run above; NN, NT, TN, TT. Run with `./run.sh`
(about 2 hours), graphs with `python3 plot_gemm.py`.

GFLOPS (2 n^3 / t), V100 against the A100:

| type | form | n=1024 | n=2048 | n=4096 | A100 n=4096 | V100/A100 |
|---|---|---:|---:|---:|---:|---:|
| dd | NN | 227.4 | 231.6 | 233.3 | 321.7 | 0.725 |
| dd | TT | | | 233.8 | 321.5 | 0.727 |
| td | NN | 70.6 | 71.1 | 71.9 | 99.8 | 0.721 |
| td | TT | | | 72.4 | 99.7 | 0.727 |
| qd | NN | 33.0 | 33.0 | 33.5 | 45.9 | 0.731 |
| qd | TT | | | 33.5 | 45.9 | 0.731 |

The four forms are within 1% of each other above n = 1000. The ratio to the
A100 is the ratio of the FP64 peaks (7.0 / 9.7 TFLOPS = 0.72): the kernels
are at the arithmetic bound on both GPUs (`../roofline`: 99-102% of the
measured madd peak).

Graphs: `gemm_{dd,td,qd}_v100_a100.png` (solid: V100, dashed: A100).
