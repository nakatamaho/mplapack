# Stage-2 fix: MPFR LIN generator scaling

Base numerical code: `b236594c4`. MPFR 512 bits, default MPFR exponent
range, OMP_NUM_THREADS=4 (8 in the initial optimized run). No exponent
clamping, precision reduction, test-type exclusion or tolerance change.

The newly required TR tests exposed a pre-existing **test-input** defect:

| Original full GE/PO/TR run | DTR | ZTR | GE/PO |
|---|---:|---:|---|
| Optimized CPU | 62/8008 failures, 24 error messages | 18/8008 failures | all four summaries PASS |
| Reference CPU | identical | identical | all four summaries PASS |

The complete reference/optimized output is identical after excluding elapsed
time. Failures are n=2, type 6 (large-norm triangular input), not the allowed
legacy CUDA BLAS zero-sign discrepancies.

A standalone `Rlattr(6,"U","N",...,n=2,seed={1988,1989,1990,1991})`
probe against the **reference** archives returned info=0 but generated an
all-zero A before TRTRI was called. Its requested norm was
`1.25215321e+323228343`; its two requested singular values were nonzero.
Reference Rtrtri correctly returned info=1 for that invalid test input.

`Rlatb4` called the existing wide-exponent `Rlabad` policy only for GMP;
`Clatb4` only for binary80/binary128. Include MPFR in these guards and in
both active FABLE 3.12.1 patches. The older 3.9.1 patches already call
Rlabad unconditionally. No production BLAS/LAPACK routine was changed.
This is the existing generator scaling policy, not a change to MPFR's
actual exponent range or a weakened residual threshold.

After the fix, the same real probe requests norm `8.64360355e+161614094`,
generates nonzero finite diagonal entries, and reference TRTRI returns 0.
The input remains far outside hardware floating-point ranges.

| Affected verification | Result |
|---|---|
| Autotools reference/optimized x real/complex, full GE/PO/TR input | 4/4 PASS |
| CMake reference/optimized x real/complex, same input and full MPFR exponent range | 4/4 PASS |
| Per run | GE 3653 + driver 5748; PO 1628 + driver 1910; TR 8008 = **20947** tests, zero failures |
| Aggregate | **167576 tests**, zero failures |
| Minimal reference DTR type-6 reproducer | 22/36 failures before; 36/36 PASS after |
| Changed LIN archives, both build systems | 4/4 `check_unique_symbols.sh` PASS |
| Active FABLE patches | reverse/apply round-trip reproduces both helpers byte for byte |

Reproduce the minimal profile (both before and after the fix):

```sh
OMP_NUM_THREADS=4 MPFRXX_DEFAULT_PRECISION_BITS=512 \
  mplapack/test/lin/mpfr/xlintstR_mpfr < \
  benchmark/results/factorization-batch-09-10/mpfr-tr-type6.in
```

Full input is the standard Rtest/Ctest parameter header with just GE 11,
PO 9 and TR 18 selected (all test types), n/m=0/1/2/3/5/10/50,
NB=1/3/3/3/20, NRHS=1/2/15, threshold 30. CUDA-library LIN runs are checked
later in the aggregate QA; this receipt does not claim those pending runs.
