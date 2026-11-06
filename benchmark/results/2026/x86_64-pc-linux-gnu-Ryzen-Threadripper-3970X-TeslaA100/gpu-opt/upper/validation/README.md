# Upper Cholesky validation

These are supplementary checks, not replacements for the repository's factor,
solver, inverse and public-BLAS suites. See the parent report and
`../provenance.json` for acceptance receipts and build/source identities.

## Reproduce the supplementary checks

Use a checkout of baseline `7e0b1f55c`, a checkout of the candidate, a configured
MPFR CMake build, and separately built baseline/candidate CUDA libraries. The
candidate checkout supplies the installed dependencies under `external/i/`.
For example, from this directory:

```sh
python3 run.py \
  --source-before /path/to/baseline \
  --source-after /path/to/candidate \
  --build /path/to/mpfr-cmake-build \
  --libs-before /path/to/baseline-libraries \
  --libs-after /path/to/candidate-libraries \
  --output /path/to/new-validation-output \
  --memcheck-updates
```

CUDA, GMP, MPFR, MPC and GMPFRXX dependencies are required. Run on an otherwise
idle GPU. Optional `--timings` also runs the isolated resident-kernel timings;
these are not whole-factorization measurements. The runner saves commands,
outputs, exit receipts, binary hashes and boundary-dump comparisons. It asserts
actual component precision, not just an environment label.

- `resident-dots.cu`: fixed-complex512 resident C/N GEMM and Upper/C HERK;
  compares baseline/candidate dumps at7,63,64,65, including padding. These are
  direct resident jobs, not the unchanged public MPFR Cgemm. Timing sizes are
 512,1024,2048, one warmup and best of three CUDA-event samples.
- `herk-boundaries.cpp`: independently renamed reference Cherk; U/C and L/C,
  sizes1,7,63,64,65 with k17, precisions256,512,768,1024,2048. Full-width mixed-sign
  operands, alpha/beta cases, diagonals and padding;50 direct calls and263100
  components per library stage. No ordinary or signed-zero allowance.
- `trsm-canary.cpp`:15 real transpose-panel cases across the same five fixed
  precisions, including short/odd right-hand-side counts.
- `ctrsm-canary.cpp`:15 complex transpose-panel cases with an explicitly non-real
  diagonal1.25+0.125i, checking actual division rather than a reciprocal rewrite.

The source/binary pairs must match. This runner does not rebuild the entire
project, validate public symbol sets, or run the full consumer corpus.

## Full acceptance corpus

The campaign also checks both CMake and autotools builds,60 BF and40 MPFR
host/device consumers,26 autotools-built host consumers, six MPFR public-BLAS
suites,22 extended triangular/HERK suites and10 n1000 factor suites. Existing
near-singular, sparse-zero and failure-output cases remain part of that corpus;
they are not presented as newly added tests here. Factor comparisons require
exact results, including `info`, pivots, signed zeros and padding.

After the final **host-launch-order-only** correction, all40 MPFR consumers,
all26 autotools hosts, all four large MPFR factors and all three MPFR-real
memchecks were rerun. The unchanged BF/public-BLAS/edge records were retained
explicitly; BF shared libraries were verified byte-identical. This is not a
claim that every retained record came from a new execution.

The pre-existing public fixed-complex BLAS zero-sign allowance was not enlarged.
Paired full suites each reported45996 checks, zero failures and1096 zero-sign
exceptions; paired quick TRSM suites each reported504 checks, zero failures and
22194 exceptions. These are matching **counts**, not an allowance for factor or
resident-boundary differences.

## Sanitizer limits

BF dd/td/qd potrf and MPFR-real512/200/2048 memchecks pass. The optional runner
memchecks only the new value-only resident complex GEMM/HERK helpers at65.
**Full MPFR complex memcheck remains UNVERIFIED under the previously approved
exception.** Passing these two helpers does not lift that exception. No stack
policy or external MPC-CUDA implementation was changed to obtain a pass.

Initial helper type/accessor compilation failures, a denied debugger attach,
rejected real-SYRK tiling and an unbuildable register-cap probe are retained in
the local campaign evidence; none supplies accepted performance data.
