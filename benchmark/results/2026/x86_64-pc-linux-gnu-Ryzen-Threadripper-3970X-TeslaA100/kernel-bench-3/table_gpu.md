| point | before (goal 11 / goal-13 part 1) | after (this session) | speedup |
|---|---|---|---|
| dd iRamax inc1 2097152 | 0.882 ms | 0.046 ms | 19.13x |
| dd iCamax inc1 2097152 | 0.825 ms | 0.065 ms | 12.79x |
| dd iRamax inc1 16777216 | 6.316 ms | 0.204 ms | 30.99x |
| dd Rdot inc1 2097152 | 0.791 ms | 0.090 ms | 8.77x |
| dd Rdot inc1 16777216 | 5.637 ms | 0.434 ms | 12.98x |
| dd Cdotc inc1 2097152 | 0.767 ms | 0.143 ms | 5.35x |
| dd Rasum inc1 2097152 | 0.794 ms | 0.073 ms | 10.92x |
| dd Rnrm2 inc1 2097152 | 0.882 ms | 0.298 ms | 2.96x |
| dd Rtrsv UNN_k64 1024 | 6.898 ms | 1.835 ms | 3.76x |
| dd Rtrsv UNN_k64 2048 | 14.993 ms | 4.008 ms | 3.74x |
| dd Rtrsv UNN_k64 4096 | 30.230 ms | 8.273 ms | 3.65x |
| dd Rtbsv UNN_k64 2048 | 13.667 ms | 3.615 ms | 3.78x |
| dd Rtpsv UNN_k64 2048 | 14.743 ms | 4.090 ms | 3.60x |
| dd Ctrsv UNN_k64 2048 | 18.779 ms | 7.481 ms | 2.51x |
| td iRamax inc1 2097152 | 1.300 ms | 0.077 ms | 16.93x |
| td Rdot inc1 2097152 | 1.191 ms | 0.125 ms | 9.53x |
| td Rtrsv UNN_k64 2048 | 18.049 ms | 6.918 ms | 2.61x |
| qd iRamax inc1 2097152 | 1.754 ms | 0.094 ms | 18.62x |
| qd Rdot inc1 2097152 | 1.634 ms | 0.213 ms | 7.67x |
| qd Rtrsv UNN_k64 2048 | 22.634 ms | 11.111 ms | 2.04x |
| mpfr512 iRamax inc1 2097152 | 6.884 ms | 1.141 ms | 6.04x |
| mpfr512 Rdot inc1 2097152 | 10.407 ms | 4.294 ms | 2.42x |
| mpfr512 Rtrsv UNN_k64 2048 | 173.895 ms | 86.714 ms | 2.01x |
