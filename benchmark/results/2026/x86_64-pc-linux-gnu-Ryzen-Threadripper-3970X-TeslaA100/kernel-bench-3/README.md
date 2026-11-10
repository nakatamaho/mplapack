# Kernel benchmark: goal 13 before/after (iamax, fast reductions, wavefront triangular solves)

Measured on the A100 80GB PCIe of this host, CUDA events around the
device-pointer schedules of `Rblas1_jobs_common.h` / `Rtrisolve_jobs.h`
(kernel only: staging, transfers and allocation outside the interval), CPU
library with OpenMP. "Before" values are the goal 11 sweep
(`../kernel-bench-2/kernels.csv`). Per the session instruction, MPFR numbers
use the current mpc_cuda v0.0.2-dev (goal 12) state.

**CPU caveat:** every CPU timing window of this session carried an external
host load average of 100-118 (64 logical cores, other tenants); CPU wall
times are therefore not trustworthy and are marked `cpu-load` below. Serial
(1 thread) parity measurements are included where they are informative
(dd fast ≈ default per-element cost; the fixed order adds no work).
The >= 8x CPU target of the goal could not be verified honestly.

## Results (type,kernel,form,n,method,seconds[,notes])

```
# goal 13 after (GPU kernel only)
dd,iRamax,inc1,2097152,gpu,0.0000471        # before 0.000881  (18.7x; 712 GB/s ~ 53% of ~1.35 TB/s)
dd,iRamax,inc1,4194304,gpu,0.0000748        # 1.04 TB/s
dd,iRamax,inc1,16777216,gpu,0.000204        # 1.31 TB/s (~97%)
dd,iCamax,inc1,2097152,gpu,0.0000645        # complex 32 B/elem: 1.04 TB/s
td,iRamax,inc1,2097152,gpu,0.0000768
td,iCamax,inc1,2097152,gpu,0.000108
qd,iRamax,inc1,2097152,gpu,0.0000932
qd,iCamax,inc1,2097152,gpu,0.000144
mpfr512,iRamax,inc1,2097152,gpu,0.001137    # before 0.006884
mpfr512,iCamax,inc1,2097152,gpu,0.002070    # before 0.007333

dd,Rdot,inc1,2097152,gpu,0.000206           # before 0.000791; 326 GB/s ~ 24% (target 50% missed)
dd,Rdot,inc1,4194304,gpu,0.000397
dd,Rdot,inc1,8388608,gpu,0.000643
dd,Cdotc,inc1,2097152,gpu,0.000377
dd,Rasum,inc1,2097152,gpu,0.000196
dd,Rnrm2,inc1,2097152,gpu,0.000309
qd,Rdot,inc1,2097152,gpu,0.000498
mpfr512,Rdot,inc1,2097152,gpu,0.004294      # before 0.0104
mpfr512,RCnrm2,inc1,2097152,gpu,0.060220

dd,Rtrsv,UNN_k64,2048,gpu,0.008179          # before 0.0150 (1.8x; target 10x missed, see README)
dd,Rtrsv,UNN_k64,4096,gpu,0.023307          # before 0.035
dd,Rtbsv,UNN_k64,2048,gpu,0.004805          # before 0.0137 (2.9x)
dd,Rtpsv,UNN_k64,2048,gpu,0.007768
mpfr512,Rtrsv,UNN_k64,2048,gpu,0.183217     # before 0.173895 (no gain; documented)
mpfr512,Rtbsv,UNN_k64,2048,gpu,0.082561     # before 0.174091

# goal 13 after (CPU; load-polluted, indicative only)
dd,iRamax,inc1,2097152,cpu-load,0.00134     # before (clean, goal 11) 0.0281
qd,iRamax,inc1,2097152,cpu-load,0.00473     # before 0.0868
mpfr512,iRamax,inc1,2097152,cpu-load,0.0108 # before 0.1468
dd,Rdot,inc1,2097152,cpu-load,0.0061-0.032
mpfr512,Rdot,inc1,2097152,cpu-load,0.24-0.26
dd,Rtrsv,UNN_k64,2048,cpu-load,~0.005       # CPU *_omp, before ~0.0051
```

## Roofline percentages (item 1 / item 2, GPU)

Bandwidth reference ~1.35 TB/s effective (A100 80GB PCIe HBM2e). Bytes =
n x element size (x + y where applicable).

| kernel | n | bytes | after | GB/s | % | target |
|---|---:|---:|---:|---:|---:|---|
| dd iRamax | 2^21 | 33.6 MB | 47.1 us | 712 | ~53% | >= 50% met |
| dd iCamax | 2^21 | 67.1 MB | 64.5 us | 1040 | ~77% | met |
| td iRamax | 2^21 | 50.3 MB | 76.8 us | 655 | ~49% | (dd is the goal's binding size) |
| qd iRamax | 2^21 | 67.1 MB | 93.2 us | 720 | ~53% | met |
| dd Rdot | 2^21 | 67.1 MB | 205.8 us | 326 | ~24% | >= 50% MISSED |
| dd Rdot | 2^24 | 537 MB | 0.75 ms | 716 | ~53% | |
| MPFR 512 | 2^21 | 134 MB (x) | 4.29 ms | — | — | bound is arithmetic (level 2 peaks), not bandwidth |

## Item 4 (inverses) status

Not implemented in this session: after items 1-3 (including the MPFR device
debugging and the repeated stale-binary detours) the session budget was
exhausted. The inverses still run the goal 11 schedule. No numbers.

## Graphs / roofline.py

Not produced for this session (see item 4 note); the CSV above carries the
measured points. The goal 11 plotting pipeline
(`benchmark/kernel_bench.py new --phase plot`) can consume this directory
once a full sweep is run.
