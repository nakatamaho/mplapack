| type | kernel | form | n | bound | % of roofline | % of GEMM plateau |
|---|---|---|---:|---|---:|---:|
| dd | Cdotc | inc1 | 2097152 | memory | 66.9 | 36.4 |
| dd | Cdotu | inc1 | 2097152 | memory | 66.4 | 36.1 |
| dd | Ctbsv | UNN_k64 | 2048 | arith | 0.0 | 0.0 |
| dd | Ctpsv | UNN_k64 | 2048 | arith | 0.7 | 0.7 |
| dd | Ctrsv | UNN_k64 | 2048 | arith | 0.7 | 0.7 |
| dd | RCasum | inc1 | 2097152 | memory | 38.1 | 10.3 |
| dd | RCnrm2 | inc1 | 2097152 | memory | 5.9 | 4.8 |
| dd | Rasum | inc1 | 2097152 | memory | 33.0 | 9.0 |
| dd | Rdot | inc1 | 2097152 | memory | 53.2 | 14.5 |
| dd | Rnrm2 | inc1 | 2097152 | memory | 8.0 | 6.6 |
| dd | Rtbsv | UNN_k64 | 2048 | memory | 0.0 | 0.0 |
| dd | Rtpsv | UNN_k64 | 2048 | memory | 0.6 | 0.3 |
| dd | Rtrsv | UNN_k64 | 2048 | memory | 0.6 | 0.3 |
| dd | gemm_sweep | NN | 2048 | arith | 99.8 |  |
| dd | iCamax | inc1 | 2097152 | memory | 74.3 | 10.1 |
| dd | iRamax | inc1 | 2097152 | memory | 52.0 |  |
| td | Rdot | inc1 | 2097152 | memory | 57.6 | 33.7 |
| td | Rtrsv | UNN_k64 | 2048 | arith | 0.6 | 0.6 |
| td | gemm_sweep | NN | 2048 | arith | 101.1 |  |
| td | iRamax | inc1 | 2097152 | memory | 46.8 |  |
| qd | Rdot | inc1 | 2097152 | memory | 45.0 | 42.9 |
| qd | Rtrsv | UNN_k64 | 2048 | arith | 0.8 | 0.8 |
| qd | gemm_sweep | NN | 2048 | arith | 99.4 |  |
| qd | iRamax | inc1 | 2097152 | memory | 50.9 |  |
| mpfr256 | gemm_sweep | NN | 2048 | arith | 24.3 |  |
| mpfr512 | Cdotc | inc1 | 2097152 | arith | 6.0 | 26.7 |
| mpfr512 | RCasum | inc1 | 2097152 | memory | 7.2 | 20.0 |
| mpfr512 | RCnrm2 | inc1 | 2097152 | arith | 0.9 | 3.3 |
| mpfr512 | Rasum | inc1 | 2097152 | memory | 4.8 | 13.4 |
| mpfr512 | Rdot | inc1 | 2097152 | memory | 5.6 | 15.5 |
| mpfr512 | Rgetri | UN | 1024 | arith | 0.2 | 1.1 |
| mpfr512 | Rlauum | UN | 1024 | arith | 0.2 | 0.7 |
| mpfr512 | Rnrm2 | inc1 | 2097152 | arith | 0.9 | 3.3 |
| mpfr512 | Rpotri | UN | 1024 | arith | 0.1 | 0.4 |
| mpfr512 | Rtbsv | UNN_k64 | 2048 | arith | 0.0 | 0.1 |
| mpfr512 | Rtrsv | UNN_k64 | 2048 | arith | 0.2 | 0.8 |
| mpfr512 | Rtrtri | UN | 1024 | arith | 0.1 | 0.3 |
| mpfr512 | gemm_sweep | NN | 2048 | arith | 21.9 |  |
| mpfr512 | iCamax | inc1 | 2097152 | memory | 11.6 | 16.0 |
| mpfr512 | iRamax | inc1 | 2097152 | memory | 10.5 |  |
| mpfr1024 | gemm_sweep | NN | 1024 | arith | 11.2 |  |
