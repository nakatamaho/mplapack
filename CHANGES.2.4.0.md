# MPLAPACK 2.4.0 - Changes since 2.3.0

MPLAPACK 2.4.0 is a **small maintenance release of the 2.x series**, based on
LAPACK 3.12.1. The only user-visible behavior change is that MPLAPACK no longer
sets the MPFR exponent range: it is left at the MPFR default unless explicitly
specified. All other changes are QA and release-engineering matters (dependency
linkage, installation, test drivers, and release validation). The unified
per-flavor library layout of 2.3.0, the ABI, and SONAME 2 are unchanged.

This summary covers the 27 commits in
`v2.3.0..00e0b659d01d073a36a58fa5851fae067fdd3a10`.

## User-visible change: MPFR exponent range

- MPLAPACK initialization no longer narrows the process-wide MPFR exponent
  range (`emin`/`emax`) according to the precision. Unless specified, the MPFR
  default exponent range (or whatever range the application has already set)
  is used as is.
- `MPLAPACK_MPFR_EMIN` and `MPLAPACK_MPFR_EMAX` remain available to specify the
  range explicitly. Applications relying on the old automatic limits should
  set their intended bounds explicitly.

## QA and release-engineering changes

These do not change the numerical routines or the library interface.

Dependencies and installation:

- Corrected GMP, MPFR/MPC, DD, and QD dependency linkage in the Autotools
  libraries and the per-flavor pkg-config metadata, propagated static QD
  math-library and optimized-library OpenMP dependencies, and improved CMake
  GMP/QD discovery and exported dependency metadata.
- Corrected the installation order of the MATGEN, LIN, and EIG support
  libraries, avoided duplicate installation of generated binary80/binary128
  utility headers, and added installed-consumer and link-order probes.
- Included `mplapack_config.h` before `<complex>` in the binary80 BLAS header.

Test drivers:

- LIN/EIG test overrides (`Mxerbla`, `iMlaenv`, `Mxlaenv`) are owned by each
  test executable instead of separately linked override libraries; DMD drivers
  own only `Mxerbla` and `Mxlaenv`. Corrected the MinGW override link order.
- Retained each backend's precision initializer in static LIN/EIG test links.
  Without it, the MinGW GMP tests ran at 64 bits instead of 512 bits.
- Pinned the MPFR compare-test exponent range to `emin=-32765`, `emax=32768`
  for reproducible reference-output comparisons.
- CMake tests use the bundled, patched reference LAPACK 3.12.1 as the
  comparison oracle, as Autotools does.

Build and release engineering:

- CMake compiles with `-fp-model=strict` under Intel compilers, matching the
  Autotools build.
- Preserved the private dlfcn-win32 files during `uninstall` so MinGW
  distcheck reinstalls work; Automake maintainer mode, `check-recursive`, and
  macOS dylib fixup cleanups.
- Benchmark log-name and header-generation fixes; Fable source-reproduction
  checks.
- Release validation organized into Tier 1 numerical/distcheck jobs, Tier 2
  tarball and C++17/20/23/26 buildability jobs, and non-gating Tier 3 checks,
  with shared source snapshots for parallel tier runs.

## Version and compatibility

- Package version: **2.4.0** in both Autotools and CMake.
- Libtool version-info remains **2:0:0**; CMake library SONAME remains **2**.
- Continue linking one `libmplapack_<flavor>` library, as with 2.3.0.
- This release retains the 2.x MPFR/MPC backend interface; it is not a migration
  to the 3.x series.

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
[the investigation report](https://github.com/nakatamaho/mplapack/blob/2.4/doc/QA.2.4.0-mingw-gmp.md).

Reproduce the summaries from the repository root:

```sh
python3 misc/summarize_mplapack_tests.py 2.4.0 mplapack/test/eig --only-fail
python3 misc/summarize_mplapack_tests.py 2.4.0 mplapack/test/lin --only-fail
```
