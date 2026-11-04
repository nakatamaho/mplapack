# Representative measured times and graph index

CPU/GPU is CPU seconds divided by GPU seconds: >1 means faster GPU.
CPU is FAST for the seven reductions, default otherwise; OMP_NUM_THREADS=32.
Prefer n=2048, then1024, otherwise the largest common measured power of two.
† marks a sampled load>=16 or other-GPU-app flag. See README for scope and omissions.
The CSV also contains the CPU-default reduction times and vector n=2^21 comparisons.

| Type | Routine | n | GPU ms | CPU ms | CPU/GPU | Flags | Graph |
|---|---|---:|---:|---:|---:|---|---|
| dd | Cdotc | 2048 | 0.04915 | 0.06002 | 1.22× |  | [PNG](dd/Cdotc.png) |
| dd | Cdotu | 2048 | 0.04915 | 0.06619 | 1.35× |  | [PNG](dd/Cdotu.png) |
| dd | Ctbsv | 2048 | 17.06 | 1.199 | 0.0703× |  | [PNG](dd/Ctbsv.png) |
| dd | Ctpsv | 2048 | 18.67 | 16.19 | 0.867× |  | [PNG](dd/Ctpsv.png) |
| dd | Ctrsv | 2048 | 18.78 | 16.41 | 0.874× |  | [PNG](dd/Ctrsv.png) |
| dd | RCasum | 2048 | 0.04915 | 0.06083 | 1.24× |  | [PNG](dd/RCasum.png) |
| dd | RCnrm2 | 2048 | 0.07987 | 0.1541 | 1.93× |  | [PNG](dd/RCnrm2.png) |
| dd | Rasum | 2048 | 0.04915 | 0.06205 | 1.26× |  | [PNG](dd/Rasum.png) |
| dd | Raxpy | 2048 | 0.00512 | 0.01041 | 2.03× |  | [PNG](dd/Raxpy.png) |
| dd | Rcopy | 2048 | 0.00512 | 0.008647 | 1.69× |  | [PNG](dd/Rcopy.png) |
| dd | Rdot | 2048 | 0.04915 | 0.06128 | 1.25× |  | [PNG](dd/Rdot.png) |
| dd | Rnrm2 | 2048 | 0.0727 | 0.09017 | 1.24× |  | [PNG](dd/Rnrm2.png) |
| dd | Rtbsv | 2048 | 13.67 | 0.4101 | 0.03× |  | [PNG](dd/Rtbsv.png) |
| dd | Rtpsv | 2048 | 14.74 | 4.98 | 0.338× |  | [PNG](dd/Rtpsv.png) |
| dd | Rtrsv | 2048 | 14.99 | 5.149 | 0.343× |  | [PNG](dd/Rtrsv.png) |
| dd | iCamax | 2048 | 0.0512 | 0.05013 | 0.979× |  | [PNG](dd/iCamax.png) |
| dd | iRamax | 2048 | 0.05222 | 0.01735 | 0.332× |  | [PNG](dd/iRamax.png) |
| mpfr512 | Cdotc | 2048 | 0.2703 | 1.148 | 4.24× |  | [PNG](mpfr512/Cdotc.png) |
| mpfr512 | Cdotu | 2048 | 0.2703 | 0.8602 | 3.18× |  | [PNG](mpfr512/Cdotu.png) |
| mpfr512 | Cgetri | 256 | 8753 | 391.5 | 0.0447× |  | [PNG](mpfr512/Cgetri.png) |
| mpfr512 | Clauum | 256 | 8352 | 219.1 | 0.0262× |  | [PNG](mpfr512/Clauum.png) |
| mpfr512 | Cpotri | 128 | 7204 | 135 | 0.0187× |  | [PNG](mpfr512/Cpotri.png) |
| mpfr512 | Ctbsv | 2048 | 308.4 | 74.78 | 0.242× |  | [PNG](mpfr512/Ctbsv.png) |
| mpfr512 | Ctpsv | 2048 | 310.1 | 1236 | 3.99× |  | [PNG](mpfr512/Ctpsv.png) |
| mpfr512 | Ctrsv | 2048 | 310.5 | 1253 | 4.03× |  | [PNG](mpfr512/Ctrsv.png) |
| mpfr512 | Ctrtri | 256 | 8259 | 181.5 | 0.022× |  | [PNG](mpfr512/Ctrtri.png) |
| mpfr512 | RCasum | 2048 | 0.1475 | 0.7915 | 5.37× |  | [PNG](mpfr512/RCasum.png) |
| mpfr512 | RCnrm2 | 2048 | 1.986 | 1.929 | 0.971× |  | [PNG](mpfr512/RCnrm2.png) |
| mpfr512 | Rasum | 2048 | 0.1454 | 0.9941 | 6.84× |  | [PNG](mpfr512/Rasum.png) |
| mpfr512 | Raxpy | 2048 | 0.0215 | 0.03054 | 1.42× |  | [PNG](mpfr512/Raxpy.png) |
| mpfr512 | Rcopy | 2048 | 0.007168 | 0.0168 | 2.34× |  | [PNG](mpfr512/Rcopy.png) |
| mpfr512 | Rdot | 2048 | 0.1884 | 0.9766 | 5.18× |  | [PNG](mpfr512/Rdot.png) |
| mpfr512 | Rgetri | 512 | 6456 | 542.1 | 0.084× |  | [PNG](mpfr512/Rgetri.png) |
| mpfr512 | Rlauum | 512 | 4281 | 229.1 | 0.0535× |  | [PNG](mpfr512/Rlauum.png) |
| mpfr512 | Rnrm2 | 2048 | 1.825 | 0.9693 | 0.531× |  | [PNG](mpfr512/Rnrm2.png) |
| mpfr512 | Rpotri | 256 | 3554 | 103.2 | 0.029× |  | [PNG](mpfr512/Rpotri.png) |
| mpfr512 | Rtbsv | 2048 | 174.1 | 17.84 | 0.102× |  | [PNG](mpfr512/Rtbsv.png) |
| mpfr512 | Rtpsv | 2048 | 174 | 282.6 | 1.62× |  | [PNG](mpfr512/Rtpsv.png) |
| mpfr512 | Rtrsv | 2048 | 173.9 | 270 | 1.55× |  | [PNG](mpfr512/Rtrsv.png) |
| mpfr512 | Rtrtri | 512 | 6125 | 186 | 0.0304× |  | [PNG](mpfr512/Rtrtri.png) |
| mpfr512 | iCamax | 2048 | 0.09626 | 1.424 | 14.8× |  | [PNG](mpfr512/iCamax.png) |
| mpfr512 | iRamax | 2048 | 0.09114 | 0.3578 | 3.93× |  | [PNG](mpfr512/iRamax.png) |
| qd | Cdotc | 2048 | 0.06656 | 0.1033 | 1.55× |  | [PNG](qd/Cdotc.png) |
| qd | Cdotu | 2048 | 0.06656 | 0.1036 | 1.56× |  | [PNG](qd/Cdotu.png) |
| qd | Ctbsv | 2048 | 38.35 | 7.756 | 0.202× |  | [PNG](qd/Ctbsv.png) |
| qd | Ctpsv | 2048 | 40.51 | 102.1 | 2.52× |  | [PNG](qd/Ctpsv.png) |
| qd | Ctrsv | 2048 | 40.52 | 101.3 | 2.5× |  | [PNG](qd/Ctrsv.png) |
| qd | RCasum | 2048 | 0.06451 | 0.1015 | 1.57× |  | [PNG](qd/RCasum.png) |
| qd | RCnrm2 | 2048 | 0.1894 | 0.3364 | 1.78× |  | [PNG](qd/RCnrm2.png) |
| qd | Rasum | 2048 | 0.06349 | 0.09744 | 1.53× |  | [PNG](qd/Rasum.png) |
| qd | Raxpy | 2048 | 0.006144 | 0.01382 | 2.25× |  | [PNG](qd/Raxpy.png) |
| qd | Rcopy | 2048 | 0.00512 | 0.009368 | 1.83× |  | [PNG](qd/Rcopy.png) |
| qd | Rdot | 2048 | 0.06451 | 0.1131 | 1.75× |  | [PNG](qd/Rdot.png) |
| qd | Rnrm2 | 2048 | 0.172 | 0.2465 | 1.43× |  | [PNG](qd/Rnrm2.png) |
| qd | Rtbsv | 2048 | 20.71 | 2.189 | 0.106× |  | [PNG](qd/Rtbsv.png) |
| qd | Rtpsv | 2048 | 22.63 | 26.92 | 1.19× |  | [PNG](qd/Rtpsv.png) |
| qd | Rtrsv | 2048 | 22.63 | 27.26 | 1.2× |  | [PNG](qd/Rtrsv.png) |
| qd | iCamax | 2048 | 0.06554 | 0.09509 | 1.45× |  | [PNG](qd/iCamax.png) |
| qd | iRamax | 2048 | 0.06554 | 0.05035 | 0.768× |  | [PNG](qd/iRamax.png) |
| td | Cdotc | 2048 | 0.05939 | 0.1098 | 1.85× |  | [PNG](td/Cdotc.png) |
| td | Cdotu | 2048 | 0.05837 | 0.1048 | 1.79× |  | [PNG](td/Cdotu.png) |
| td | Ctbsv | 2048 | 25.02 | 5.508 | 0.22× |  | [PNG](td/Ctbsv.png) |
| td | Ctpsv | 2048 | 27.63 | 78.33 | 2.83× |  | [PNG](td/Ctpsv.png) |
| td | Ctrsv | 2048 | 27.42 | 78.9 | 2.88× |  | [PNG](td/Ctrsv.png) |
| td | RCasum | 2048 | 0.05734 | 0.1031 | 1.8× |  | [PNG](td/RCasum.png) |
| td | RCnrm2 | 2048 | 0.128 | 0.2046 | 1.6× |  | [PNG](td/RCnrm2.png) |
| td | Rasum | 2048 | 0.05734 | 0.08823 | 1.54× |  | [PNG](td/Rasum.png) |
| td | Raxpy | 2048 | 0.006144 | 0.01106 | 1.8× |  | [PNG](td/Raxpy.png) |
| td | Rcopy | 2048 | 0.00512 | 0.009929 | 1.94× |  | [PNG](td/Rcopy.png) |
| td | Rdot | 2048 | 0.05734 | 0.09478 | 1.65× |  | [PNG](td/Rdot.png) |
| td | Rnrm2 | 2048 | 0.1147 | 0.2367 | 2.06× |  | [PNG](td/Rnrm2.png) |
| td | Rtbsv | 2048 | 16.57 | 1.313 | 0.0792× |  | [PNG](td/Rtbsv.png) |
| td | Rtpsv | 2048 | 17.94 | 17.18 | 0.958× |  | [PNG](td/Rtpsv.png) |
| td | Rtrsv | 2048 | 18.05 | 17.33 | 0.96× |  | [PNG](td/Rtrsv.png) |
| td | iCamax | 2048 | 0.06042 | 0.05586 | 0.925× |  | [PNG](td/iCamax.png) |
| td | iRamax | 2048 | 0.06042 | 0.01909 | 0.316× |  | [PNG](td/iRamax.png) |
