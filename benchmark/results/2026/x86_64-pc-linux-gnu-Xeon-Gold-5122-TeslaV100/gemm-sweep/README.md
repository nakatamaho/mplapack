# dd/td/qd Rgemm on the V100 (November 22, 2026)

Tesla V100-PCIE-32GB (sm_70, 1380 MHz, CUDA 12.6 of HPC SDK 24.11), commit
9d8716b9f, GPU 0 alone (`gpu-processes-at-*.txt`). Kernel only:
`benchmark/run_gemm_sweep.sh` (data on the GPU, CUDA events, default tiled
kernel). Matrix sizes are multiples of 7 and 31 and powers of two up to 4096, best of 3 after a warm-up
for n <= 2000, one run above; NN, NT, TN, TT. Run with `./run.sh`
(about 2 hours), graphs with `python3 plot_gemm.py`.

GFLOPS (2 n^3 / t):

| type | form | n=1024 | n=2048 | n=4096 |
|---|---|---:|---:|---:|
| dd | NN | 227.4 | 231.6 | 233.3 |
| dd | TT | | | 233.8 |
| td | NN | 70.6 | 71.1 | 71.9 |
| td | TT | | | 72.4 |
| qd | NN | 33.0 | 33.0 | 33.5 |
| qd | TT | | | 33.5 |

The four forms are within 1% of each other above n = 1000. The kernels
run near the arithmetic bound (99-102% of the measured madd peak; see
[`../roofline`](../roofline/README.md)).

Graphs: [dd](gemm_dd.png), [td](gemm_td.png), [qd](gemm_qd.png).
