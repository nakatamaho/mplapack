# MPLAPACK 2.4.0: MinGW/GMP numerical QA investigation

Investigation date: 2026-10-18. No production sources or build rules were
changed during this investigation; the tested release artifact is preserved.

## Findings

The collected 22 MinGW/GMP EIG failures are genuine failure summaries, not
an aggregation-only error. They occur in Csep, Cse2, Rsep, and Rse2. Several
calls to Cheevx/Chpevx/Cheevr, Rspevx, and Cstein/Rstein return INFO=1 for
clustered-eigenvalue matrices of order 25. Ratios of approximately
9.2233720368547758e18 in these records are failure sentinels, not merely
small threshold overruns.

More importantly, the GMP tests do not use the intended precision:

| Results | Reported machine precision | Effective GMP precision |
| --- | --- | ---: |
| 2.3.0 MinGW GMP Csep | 7.4583407312002067e-155 | 512 bits |
| 2.4.0 MinGW GMP Csep (both host versions) | 5.4210108624275222e-20 | 64 bits |
| 2.4.0 other platforms GMP Csep | 7.4583407312002067e-155 | 512 bits |

The same 64-bit precision is reported throughout the collected MinGW GMP
LIN outputs, even though those tests report no numerical failures.
Consequently, the LIN all-pass count is not evidence of a 512-bit MinGW/GMP
validation run.

Both original host result directories were inspected, not just the combined
repository result tree:

- `/home/maho/tmp/mplapack-ubuntu2404-mingw64-amd64.distcheck-results/`
- `/home/maho/tmp/mplapack-ubuntu2604-mingw64-amd64.distcheck-results/`

Both contain the same Csep precision and failing-summary pattern. This rules
out the shared triplet destination being the sole explanation for these
failures, although the destination-name collision remains an audit limitation.

## Linkage mechanism

`include/mpblas_gmp.h` specifies a 512-bit default.
`mpblas/reference/mplapackinit.cpp` contains a constructor,
`mplapack_initialize_gmp()`, that sets this precision and optionally honors
`MPLAPACK_GMP_PRECISION`.

The 2.3.0 MinGW test link rules used `--whole-archive` for the backend archive.
The current rules in `mplapack/test/eig/gmp/Makefile.am` and the actual QA
link commands instead pass `libmplapack_gmp.a` as an ordinary static archive.
The initialization function has no ordinary caller and can therefore be
omitted from archive extraction. A constructor attribute does not by itself
cause its archive member to be selected.

An isolated Linux linker experiment used the actual, unchanged
`mpblas/reference/mplapackinit.cpp`, compiled with the GMP backend macro,
then placed that object in an archive. A small main printed
`mpf_get_default_prec()` with `MPLAPACK_GMP_PRECISION` unset:

```text
Ordinary archive link:                          64 bits
Link with -Wl,-u,_Z23mplapack_initialize_gmpv:   512 bits
```

`nm` confirms that the initialization symbol is absent from the ordinary
executable and present in the forced-symbol executable. This reproduces the
archive-extraction mechanism on Linux; it is not a completed Windows numerical
retest. Temporary experiment files are under
`/tmp/mplapack-240-init-check.hN9L0Q/`.

## Disposition

Treat missing backend initialization as a release-blocking test-link issue,
not as an accepted numerical tolerance difference. It is a strong explanation
for the precision loss and a plausible explanation for the new GMP failures,
but **the disappearance of all 22 failures has not yet been verified**.

Before publication:

1. Ensure the GMP initializer is retained in MinGW static test links (ordinary,
   optimized, LIN, EIG, and DMD as appropriate), without restoring duplicate
   override/common-array ownership. Also audit other backend initializers and
   CMake static consumers for the same archive-extraction issue.
2. Add a runtime precision check so job-level OK cannot hide a missing
   initializer. Setting only an environment variable cannot help if the
   initializer that reads it was never linked.
3. Rebuild/relink and rerun the affected MinGW tests, checking 512-bit machine
   precision and the numerical summaries. Do not simply raise thresholds.
4. If build rules change, produce a new source snapshot and tarball with new
   provenance. The current QA artifact must not be represented as containing
   such a fix.

## Subsequent fix and focused validation

After approval, all seven MinGW backend rules (GMP, MPFR, QD, DD, double,
binary80, and binary128) in LIN and EIG now retain their matching initializer
using `-Wl,-u,<C++ ABI symbol>`. The shared `test_ldflags` covers reference,
optimized, and DMD targets. The production routine bodies are unchanged.
CMake static test executables use equivalent retention; shared-library test
links are unchanged. The link-rule checker now checks initializer retention,
and the GMP override probe verifies the effective runtime precision.

Validation performed:

- Regenerated Automake outputs; all 14 backend/family rules and their generated
  files contain the matching initializer symbol.
- Built GMP 6.3.0 for MinGW in an isolated Docker container. Compiled the actual
  MPLAPACK GMP initializer and linked the probe as a Windows executable.
  Wine confirmed 64 bits for the unretained negative control, 512 bits with
  retention, and 128 bits with `MPLAPACK_GMP_PRECISION=128`.
- MinGW/Wine LIN and EIG override-order regression probes passed, including
  their archive-first negative controls.
- CMake static double LIN, EIG, and DMD override probes built and passed (3/3).
  Generated link lines contain the double initializer retention option.
- Duplicate-symbol checks passed for the produced archives; `git diff --check`
  passed. Cross-test logs: `/tmp/mplapack-240-cross-init-output/`.

The full MinGW numerical retest, including confirmation that all 22 GMP
failures disappear, remains pending. The old QA tarball has not been changed
and does not include this fix. A new source snapshot is required for release QA.

No tag, public release, or announcement was published during this check.

## Release QA retest (resolved)

A new source snapshot was made from `00e0b659d01d` (label
`00e0b659d01d-bfdd8c64498a`, run `release/logs/20261020_143347/`), which
contains the initializer-retention fix. All 12 Tier 1 and 4 Tier 2 jobs
report OK.

- Both MinGW hosts (Ubuntu 24.04 and 26.04) report machine precision
  `7.4583407312002067e-155` (512 bits) in every GMP LIN and EIG output.
  The same value is reported by all 11 result triplets.
- The 22 MinGW/GMP EIG failures (Csep, Cse2, Rsep, Rse2) are gone; GMP has
  no failures on any platform. Both original per-host result directories
  were checked, not only the combined repository tree.
- The remaining EIG failures (8 tests in 6 records: Debian i386 double
  Rse2/Rsvd, macOS 15 and MinGW binary80 Rsep) match those published for 2.3.0.
- LIN: 0 failures in 126,515,712 recognized tests.

The missing initializer was the cause of the 22 failures. The issue is closed.
