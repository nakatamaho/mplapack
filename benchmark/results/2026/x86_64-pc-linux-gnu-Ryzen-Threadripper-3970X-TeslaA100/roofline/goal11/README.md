# Roofline of the goal 11 kernels (BLAS1, BLAS2 triangular solves, MPFR inverses)

Same method as ../README.md (peaks ../peaks.csv, LAWN 41 counts, % = max(t_arith, t_mem) / t_measured), applied to
../../kernel-bench-2/kernels.csv (phase = sweep only). Representative sizes: vectors n = 2^21, matrices n = 2048
(MPFR inverses: the largest n measured, 248-1022). `roofline_summary.md` has every kernel; `roofline_<type>.png` the curves.

Counts used for the new kernels: axpy n mul + n add; dot 2n (complex: 4n mul + 4n add); asum n add (complex 2n);
nrm2 2n mul + n add (the scaled square and sum; complex 4n + 2n); copy, iamax no arithmetic (memory only; iCamax
n additions); trsv / tpsv n(n-1)/2 mul + add + n div; tbsv (band k = 64) n·k mul + add + n div; getri, trtri, lauum,
potri (= trtri + lauum) from LAWN 41.

| kernel (n) | dd | td | qd | MPFR 512 |
|---|---:|---:|---:|---:|
| Raxpy (2^21, memory bound) | 94 | 100 | 52 | 7.8 |
| Rcopy (2^21, memory bound) | 84 | 95 | 37 | 57 |
| Rdot fast (2^21) | 6.1 | 6.0 | 5.9 | 2.3 |
| Cdotc fast (2^21) | 12.5 | 14.7 | 28.5 | 4.9 |
| Rasum fast (2^21) | 3.0 | 3.1 | 3.0 | 1.4 |
| Rnrm2 fast (2^21) | 2.7 | 4.1 | 6.1 | 0.6 |
| iRamax (2^21, memory bound) | 2.7 | 2.8 | 2.7 | 1.7 |
| Rtrsv / Rtpsv (2048) | 0.2 | 0.2 | 0.4 | 0.1 |
| Ctrsv / Ctpsv (2048) | 0.3 | 0.6 | 0.9 | 0.3 |
| Rtbsv / Ctbsv (2048, k = 64) | 0.0 | 0.0 | 0.0-0.1 | 0.0 |
| MPFR Rgetri / Rtrtri / Rlauum (1022) | | | | 0.2 / 0.1 / 0.2 |
| MPFR Rpotri (511), Cgetri / Ctrtri / Clauum (496) | | | | 0.0-0.1 |

Reading: the element-wise kernels reach the memory bandwidth for dd/td; the fast reductions and iamax are far from
it (a few %, i.e. the reduction tree / launches dominate); the triangular solves and the MPFR inverses are
latency bound (one dependent step per column), which a throughput roofline does not capture: their % is the
distance from what the throughput would allow, not a measure of a fixable constant factor alone.
