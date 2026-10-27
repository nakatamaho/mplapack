# MPLAPACK 2.4.0 - Changes since 2.3.0

MPLAPACK 2.4.0 is a **maintenance release of the 2.x series**, based on
LAPACK 3.12.1. It focuses on dependency linkage, installation, test-driver
correctness, and reproducible release QA. It retains the unified per-flavor
library layout introduced in 2.3.0; it does not introduce another library-layout
or ABI transition.

This summary covers the 27 commits in
`v2.3.0..00e0b659d01d073a36a58fa5851fae067fdd3a10`.
The latest published 2.3.x release is 2.3.0.

## Backend dependencies and installed consumers

- Corrected GMP, MPFR/MPC, DD, and QD dependency linkage in the Autotools
  libraries and public per-flavor pkg-config metadata. GMP consumers receive
  gmpxx and gmp; MPFR consumers receive mpc, mpfr, gmpxx, and gmp.
- Propagated static QD math-library dependencies and optimized-library OpenMP
  runtime dependencies. `OPENMP_LIBS` can override runtime discovery;
  `--disable-openmp` leaves that interface empty.
- Improved CMake GMP/QD discovery and exported dependency metadata, including
  custom prefixes without matching pkg-config files. A host pkg-config file
  for a different library must not supply the selected library's dependencies.
- Corrected installation order for MATGEN, LIN, and EIG support libraries,
  and added installed-consumer, shared-loader, dependency, and link-order probes.
- Avoided duplicate installation of generated binary80/binary128 utility headers.

## Test-driver ownership and MinGW fixes

- Made LIN/EIG test overrides (`Mxerbla`, `iMlaenv`, and `Mxlaenv`) belong to
  each test executable instead of separately linked override libraries. This
  avoids duplicate ownership and initialization/finalization of test common
  arrays. Both Autotools and CMake use the executable-owned model.
- DMD drivers own only `Mxerbla` and `Mxlaenv`, retaining the production
  `iMlaenv` tuning behavior.
- Corrected MinGW LIN/EIG override link order and added focused error-exit,
  override ownership, and lifetime regression checks.
- Retained each backend's precision initializer in static LIN/EIG test links
  (`-Wl,-u,<initializer>` for all seven MinGW backends; equivalent retention
  for CMake static test executables). An ordinary static archive link could
  omit the constructor-only initializer, so the MinGW GMP tests ran at 64 bits
  instead of 512 bits. The GMP override probe now checks the effective runtime
  precision. Production routine bodies are unchanged.
- Preserved the private dlfcn-win32 headers and libraries during `uninstall`.
  Previously, distcheck's subsequent DESTDIR installation could try to rebuild
  benchmarks after `dlfcn.h` had been deleted while its install stamp remained.
  Private dependency cleanup now occurs at `distclean` instead.

## MPFR behavior and header compatibility

- Removed automatic precision-dependent narrowing of the process-wide MPFR
  exponent range during MPLAPACK initialization. Without explicit overrides,
  the existing MPFR exponent range is preserved.
- `MPLAPACK_MPFR_EMIN` and `MPLAPACK_MPFR_EMAX` remain available for applications
  that need explicit bounds. Applications relying on the old automatic limits
  should set their intended bounds explicitly.
- Pinned the MPFR compare-test exponent range to `emin=-32765`, `emax=32768`
  for reproducible reference-output comparisons.
- Included `mplapack_config.h` before `<complex>` in the binary80 BLAS header.

## Build and release engineering

- Enabled disabled-by-default Automake maintainer mode for distributed sources.
- Avoided duplicate Automake `check-recursive` recipes and corrected installed
  example directory checks.
- Made the macOS dylib fixup skip directories rather than treating them as
  candidate libraries.
- Added Tier 1 Fable source-reproduction comparison and regression checks for
  generated sources, headers, source lists, and inputs.
- Improved shared source-snapshot preparation for parallel tier runs.
- CMake: compile with `-fp-model=strict` under Intel compilers, matching the
  Autotools build; Intel's default value-unsafe model breaks the DD/QD
  expansion arithmetic. The option is also exported with the no-FMA backend
  options.
- CMake tests: use the bundled, patched reference LAPACK 3.12.1 as the
  comparison oracle, as Autotools does, and link its Fortran runtime. An
  optimized system BLAS is not a valid oracle.
- CMake: no longer build the mpblas test support file as a test.
- Benchmarks: fixed plot log names, stopped linking binary80 with
  `-lquadmath`, added `misc/gen_benchmark_header.py`, and regenerated the
  benchmark headers from the per-flavor libraries.
- Fable: removed dead duplicate hunks from `patch-Claror.cpp`.
- Organized release validation into Tier 1 numerical/distcheck jobs, Tier 2
  tarball and C++ standard buildability jobs, and non-gating Tier 3 checks.
  C++17, C++20, C++23, and C++26 buildability is checked on Ubuntu 26.04
  amd64 and arm64.

## Version and compatibility

- Package version: **2.4.0** in both Autotools and CMake.
- Libtool version-info remains **2:0:0**; CMake library SONAME remains **2**.
- Continue linking one `libmplapack_<flavor>` library, as with 2.3.0.
- This release retains the 2.x MPFR/MPC backend interface; it is not a migration
  to the 3.x series.
- The MPFR exponent-range behavior change above is relevant even though the
  public library layout and ABI version remain unchanged.

## Release artifact and QA provenance

The QA source is a clean checkout (no tracked source patch):

- Git commit: `00e0b659d01d073a36a58fa5851fae067fdd3a10`
- Source label: `00e0b659d01d-bfdd8c64498a`
- QA run: `release/logs/20261020_143347/`
- Artifact: `mplapack-2.4.0.tar.xz`
- SHA-256: `bfdd8c64498a1b9472b98d3f40bcf94a51fb3e54f2921996cf1eb2c85c6093df`

All **12 Tier 1 jobs** and **4 Tier 2 jobs** report `OK` in this run's CSV
records. Tier 1 covers macOS 15 amd64, macOS 26 arm64, Ubuntu 24.04/26.04
amd64 and arm64, MinGW-w64 via both Ubuntu versions, Intel oneAPI via both
Ubuntu versions, and Debian 12/13 i386. Tier 2 covers Ubuntu 24.04 tarball
builds and Ubuntu 26.04 C++ standard builds, each on amd64 and arm64.

Job completion is not equivalent to zero numerical threshold failures.
The collected numerical result trees summarize as follows:

| Category | Triplets | Precision groups | Recognized tests | Failed tests | Failed records |
| --- | ---: | ---: | ---: | ---: | ---: |
| EIG | 11 | 96 | 116,310,045 | 8 | 6 |
| LIN | 11 | 96 | 126,515,712 | 0 | 0 |

These are counts of recognized summary records, not a claim of unique test
cases. Two MinGW host configurations share the same result triplet name, so
the repository result tree does not independently preserve both; the original
per-host result directories were checked separately and agree.

EIG failure records:

| Platform/compiler | Backend | File | Suite | Failed/Total |
| --- | --- | --- | --- | ---: |
| Debian 12 GCC 12.2.0 | double | Rse2.out | DST | 1/4,440 |
| Debian 13 GCC 14.2.0 | double | Rse2.out | DST | 1/4,440 |
| Debian 12 GCC 12.2.0 | double | Rsvd.out | DBD | 1/10,260 |
| Debian 13 GCC 14.2.0 | double | Rsvd.out | DBD | 1/10,260 |
| macOS 15 GCC 15.2.0 | binary80 | Rsep.out | DST drivers | 2/13,464 |
| MinGW GCC 13.x | binary80 | Rsep.out | DST drivers | 2/13,464 |

These match the platform/backend/file/suite failure counts reported in the
2.3.0 release notes. All GMP results, on every platform including both MinGW
hosts, report the intended 512-bit machine precision
(`7.4583407312002067e-155`) and have no failures.

An earlier QA run of `6de413e4a75f` found 22 additional MinGW/GMP EIG failures.
They were caused by the MinGW GMP tests running at 64 bits because the static
test link omitted the backend initializer; they disappear with the
initializer-retention fix above. See
[the investigation report](doc/QA.2.4.0-mingw-gmp.md).

Reproduce the summaries from the repository root:

```sh
python3 misc/summarize_mplapack_tests.py 2.4.0 mplapack/test/eig --only-fail
python3 misc/summarize_mplapack_tests.py 2.4.0 mplapack/test/lin --only-fail
```
