| type | kernel | form | n | bound | % of roofline | % of GEMM plateau |
|---|---|---|---:|---|---:|---:|
| dd | Cdotc | inc1 | 2097152 | memory | 12.5 | 6.8 |
| dd | Cdotu | inc1 | 2097152 | memory | 12.5 | 6.8 |
| dd | Ctbsv | UNN_k64 | 2048 | arith | 0.0 | 0.0 |
| dd | Ctpsv | UNN_k64 | 2048 | arith | 0.3 | 0.3 |
| dd | Ctrsv | UNN_k64 | 2048 | arith | 0.3 | 0.3 |
| dd | RCasum | inc1 | 2097152 | memory | 6.5 | 1.8 |
| dd | RCnrm2 | inc1 | 2097152 | memory | 3.0 | 2.5 |
| dd | Rasum | inc1 | 2097152 | memory | 3.0 | 0.8 |
| dd | Raxpy | inc1 | 2097152 | memory | 93.9 | 17.0 |
| dd | Rcopy | inc1 | 2097152 | memory | 83.7 |  |
| dd | Rdot | inc1 | 2097152 | memory | 6.1 | 1.6 |
| dd | Rnrm2 | inc1 | 2097152 | memory | 2.7 | 2.2 |
| dd | Rtbsv | UNN_k64 | 2048 | memory | 0.0 | 0.0 |
| dd | Rtpsv | UNN_k64 | 2048 | memory | 0.2 | 0.1 |
| dd | Rtrsv | UNN_k64 | 2048 | memory | 0.2 | 0.1 |
| dd | gemm_sweep | NN | 2048 | arith | 99.8 |  |
| dd | iCamax | inc1 | 2097152 | memory | 5.8 | 0.8 |
| dd | iRamax | inc1 | 2097152 | memory | 2.7 |  |
| td | Cdotc | inc1 | 2097152 | arith | 14.7 | 14.5 |
| td | Cdotu | inc1 | 2097152 | arith | 14.7 | 14.5 |
| td | Ctbsv | UNN_k64 | 2048 | arith | 0.0 | 0.0 |
| td | Ctpsv | UNN_k64 | 2048 | arith | 0.6 | 0.6 |
| td | Ctrsv | UNN_k64 | 2048 | arith | 0.6 | 0.6 |
| td | RCasum | inc1 | 2097152 | memory | 6.5 | 3.8 |
| td | RCnrm2 | inc1 | 2097152 | arith | 4.5 | 4.7 |
| td | Rasum | inc1 | 2097152 | memory | 3.1 | 1.8 |
| td | Raxpy | inc1 | 2097152 | memory | 100.1 | 39.0 |
| td | Rcopy | inc1 | 2097152 | memory | 95.1 |  |
| td | Rdot | inc1 | 2097152 | memory | 6.0 | 3.5 |
| td | Rnrm2 | inc1 | 2097152 | arith | 4.1 | 4.2 |
| td | Rtbsv | UNN_k64 | 2048 | arith | 0.0 | 0.0 |
| td | Rtpsv | UNN_k64 | 2048 | arith | 0.2 | 0.2 |
| td | Rtrsv | UNN_k64 | 2048 | arith | 0.2 | 0.2 |
| td | gemm_sweep | NN | 2048 | arith | 101.1 |  |
| td | iCamax | inc1 | 2097152 | memory | 5.8 | 1.7 |
| td | iRamax | inc1 | 2097152 | memory | 2.8 |  |
| qd | Cdotc | inc1 | 2097152 | arith | 28.5 | 28.6 |
| qd | Cdotu | inc1 | 2097152 | arith | 28.5 | 28.6 |
| qd | Ctbsv | UNN_k64 | 2048 | arith | 0.1 | 0.1 |
| qd | Ctpsv | UNN_k64 | 2048 | arith | 0.9 | 0.9 |
| qd | Ctrsv | UNN_k64 | 2048 | arith | 0.9 | 0.9 |
| qd | RCasum | inc1 | 2097152 | memory | 6.5 | 6.2 |
| qd | RCnrm2 | inc1 | 2097152 | arith | 6.7 | 6.7 |
| qd | Rasum | inc1 | 2097152 | memory | 3.0 | 2.9 |
| qd | Raxpy | inc1 | 2097152 | memory | 51.6 | 32.8 |
| qd | Rcopy | inc1 | 2097152 | memory | 37.3 |  |
| qd | Rdot | inc1 | 2097152 | memory | 5.9 | 5.6 |
| qd | Rnrm2 | inc1 | 2097152 | arith | 6.1 | 6.1 |
| qd | Rtbsv | UNN_k64 | 2048 | arith | 0.0 | 0.0 |
| qd | Rtpsv | UNN_k64 | 2048 | arith | 0.4 | 0.4 |
| qd | Rtrsv | UNN_k64 | 2048 | arith | 0.4 | 0.4 |
| qd | gemm_sweep | NN | 2048 | arith | 99.4 |  |
| qd | iCamax | inc1 | 2097152 | memory | 5.8 | 2.8 |
| qd | iRamax | inc1 | 2097152 | memory | 2.7 |  |
| mpfr256 | gemm_sweep | NN | 2048 | arith | 54.2 |  |
| mpfr512 | Cdotc | inc1 | 2097152 | arith | 7.7 | 14.3 |
| mpfr512 | Cdotu | inc1 | 2097152 | arith | 7.5 | 14.0 |
| mpfr512 | Cgetri | UN | 496 | arith | 0.2 | 0.4 |
| mpfr512 | Clauum | UN | 496 | arith | 0.1 | 0.1 |
| mpfr512 | Cpotri | UN | 248 | arith | 0.0 | 0.0 |
| mpfr512 | Ctbsv | UNN_k64 | 2048 | arith | 0.0 | 0.1 |
| mpfr512 | Ctpsv | UNN_k64 | 2048 | arith | 0.5 | 0.9 |
| mpfr512 | Ctrsv | UNN_k64 | 2048 | arith | 0.5 | 0.9 |
| mpfr512 | Ctrtri | UN | 496 | arith | 0.1 | 0.1 |
| mpfr512 | RCasum | inc1 | 2097152 | arith | 2.8 | 7.5 |
| mpfr512 | RCnrm2 | inc1 | 2097152 | arith | 1.0 | 1.7 |
| mpfr512 | Rasum | inc1 | 2097152 | arith | 1.5 | 3.9 |
| mpfr512 | Raxpy | inc1 | 2097152 | memory | 7.8 | 14.5 |
| mpfr512 | Rcopy | inc1 | 2097152 | memory | 56.8 |  |
| mpfr512 | Rdot | inc1 | 2097152 | arith | 3.4 | 6.4 |
| mpfr512 | Rgetri | UN | 1022 | arith | 0.3 | 0.6 |
| mpfr512 | Rlauum | UN | 1022 | arith | 0.3 | 0.5 |
| mpfr512 | Rnrm2 | inc1 | 2097152 | arith | 1.0 | 1.7 |
| mpfr512 | Rpotri | UN | 511 | arith | 0.1 | 0.1 |
| mpfr512 | Rtbsv | UNN_k64 | 2048 | arith | 0.0 | 0.0 |
| mpfr512 | Rtpsv | UNN_k64 | 2048 | arith | 0.2 | 0.4 |
| mpfr512 | Rtrsv | UNN_k64 | 2048 | arith | 0.2 | 0.4 |
| mpfr512 | Rtrtri | UN | 1022 | arith | 0.1 | 0.2 |
| mpfr512 | gemm_sweep | NN | 2048 | arith | 52.8 |  |
| mpfr512 | iCamax | inc1 | 2097152 | memory | 3.3 | 4.5 |
| mpfr512 | iRamax | inc1 | 2097152 | memory | 1.7 |  |
| mpfr1024 | gemm_sweep | NN | 1024 | arith | 26.5 |  |
