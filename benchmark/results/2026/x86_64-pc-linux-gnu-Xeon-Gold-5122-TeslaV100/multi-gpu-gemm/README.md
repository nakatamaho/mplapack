# dd/td/qd Rgemm on one and two V100s (November 25, 2026)

Two Tesla V100-PCIE-32GB (sm_70) on one Xeon Gold 5122, both on NUMA node 0,
connected through PCIe only (`gpu-topo.txt`: NODE, no NVLink). CUDA 12.6 of
HPC SDK 24.11, commit in `commit.txt`. GPUs otherwise idle
(`gpu-processes-at-*.txt`, `run.log`).

The public `Rgemm` of libmplapack_{dd,td,qd}_opt_cuda (`benchmark/Rgemm.cpp`,
`Rgemm.{dd,td,qd}_cuda_total`): the whole call, host to device transfers,
kernels and transfers back. `MPLAPACK_{DD,TD,QD}_CUDA_NGPUS=1` against `=2`
(`mpblas/optimized/common/cuda_multi_gpu.h`: the columns of C are cut into one
block per device, each computed as on one device, so the result does not
depend on the number of devices). NN, m = n = k, one warm-up, then 3 runs;
median below, all statistics in `{dd,td,qd}_ngpus{1,2}.csv` (MFLOPS). No
reference check in the timed runs (`-NOCHECK`); a separate check at n = 600
against the reference library gave error 0 for both settings (dd and qd).
Run with `./run.sh` (about 14 min).

GFLOPS (2 n^3 / t):

| type | n | 1 GPU | 2 GPUs | speedup |
|---|---:|---:|---:|---:|
| dd | 1024 | 104.3 | 152.7 | 1.46 |
| dd | 2048 | 156.8 | 264.4 | 1.69 |
| dd | 4096 | 189.0 | 355.8 | 1.88 |
| dd | 8192 | 208.9 | 404.9 | 1.94 |
| td | 1024 | 48.1 | 73.8 | 1.53 |
| td | 2048 | 58.4 | 109.2 | 1.87 |
| td | 4096 | 64.7 | 125.8 | 1.94 |
| td | 8192 | 68.1 | 133.9 | 1.97 |
| qd | 1024 | 25.7 | 41.4 | 1.61 |
| qd | 2048 | 29.2 | 55.9 | 1.92 |
| qd | 4096 | 31.3 | 61.4 | 1.96 |
| qd | 8192 | 32.4 | 64.0 | 1.98 |

Two GPUs reach 1.9-2.0x from n = 4096 on. The scaling is lowest for dd,
whose kernel is the fastest relative to the PCIe transfers: on one GPU the
whole call at n = 4096 runs at 81% of the kernel-only rate of
`../gemm-sweep` (dd 233.3, td 71.9, qd 33.5 GFLOPS), against 90% for td and
93% for qd. At n = 1024 dd with two GPUs varies between runs (cv 14%).
