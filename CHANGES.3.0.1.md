# MPLAPACK 3.0.1 - Release notes

Release date: 2026-09-15 (release baseline: `953d7a4916554546937a753a30b0619691072841`)

MPLAPACK 3.0.1 is a build, packaging, and release-engineering patch release
following 3.0.0.  It does not change the numerical API or the libtool ABI.

## Release artifact

- Archive: `mplapack-3.0.1.tar.xz`
- Source commit: `953d7a4916554546937a753a30b0619691072841`
- SHA256: `47ebb653b21f0c62e8144c1e515e94d76965d9d0d4b7ba216b034d778570cbaa`
- MD5: `8963437912a012d9f00e972c4d2dfde4`
- Archive source metadata: `release/logs/20260910_175329/source/source-metadata.txt`

The archive was generated from a clean tracked source tree at the release
baseline.  The final release tag also contains the result snapshot and the
release documentation added after that source archive was generated.

## Build and dependency fixes

- Corrected shared and static dependency metadata for the MPFR, GMP, QD, and
  DD backends across Autotools/libtool, CMake, and pkg-config.
- Preserved the complete dependency interface discovered by the CMake GMP and
  QD find modules, including QD's required math-library dependency.
- Added clean-process shared-library load, relocation, and external consumer
  checks for backend dependency regressions.
- Added the MPFR precision-scope public header to installation and release
  manifests, and exposed the external `gmpfrxx_mkII` include path in the MPFR
  pkg-config metadata.
- Fixed OpenMP runtime dependency propagation for optimized libraries and the
  corresponding generated Autoconf case syntax.
- Hardened macOS QD/DD shared-library load checks and portable pkg-config
  generation.
- Added bundled CMake bootstrapping for legacy environments and disabled
  OpenSSL for that bootstrap because MPLAPACK does not require it.

The package version is 3.0.1 in both Autotools and CMake.  The libtool ABI
version remains `3:0:0`; no SONAME or public numerical API change is intended.

## Test result summary

Results were collected under:

- `mplapack/test/eig/results/3.0.1/`
- `mplapack/test/lin/results/3.0.1/`

The snapshot covers 11 result triplets:

- `Apple_M4_macos26_gcc-15_2_0`
- `Apple_M4_ubuntu24_04_gcc-13_3_0`
- `Apple_M4_ubuntu26_04_gcc-15_2_0`
- `Core_i7-6920HQ_macos15_gcc-15_2_0`
- `Ryzen_Threadripper_3970X_debian12_gcc-12_2_0`
- `Ryzen_Threadripper_3970X_debian13_gcc-14_2_0`
- `Ryzen_Threadripper_3970X_ubuntu24_04_gcc-13_3_0`
- `Ryzen_Threadripper_3970X_ubuntu24_04_icx-2026_1_1`
- `Ryzen_Threadripper_3970X_ubuntu26_04_gcc-15_2_0`
- `Ryzen_Threadripper_3970X_ubuntu26_04_icx-2026_1_1`
- `Ryzen_Threadripper_3970X_windows_gcc-13_x_x`

| Category | Triplets | `.out` files | Recognized tests | Failed tests | Failed records |
| --- | ---: | ---: | ---: | ---: | ---: |
| eig | 11 | 4,032 | 116,811,306 | 30 | 28 |
| lin | 11 | 384 | 126,515,712 | 0 | 0 |

The recurring `eig` threshold-edge results are unchanged from 3.0.0:

| Backend/result | Environments | Failed/total per environment | Aggregate failed/total |
| --- | ---: | ---: | ---: |
| GMP `Cgd.out` / ZGV drivers | 11 | 1/1,092 | 11/12,012 |
| MPFR default `Rgg.default.out` / DGG | 11 | 1/3,120 | 11/34,320 |
| double `Rse2.out` / DST | 2 | 1/4,440 | 2/8,880 |
| double `Rsvd.out` / DBD | 2 | 1/10,260 | 2/20,520 |
| binary80 `Rsep.out` / DST drivers | 2 | 2/13,464 | 4/26,928 |

The GMP and MPFR cases are identical one-test threshold exceedances on all
11 triplets.  The double and binary80 cases are the previously documented
platform-specific differences.  These 30 threshold-edge results are retained
for traceability and are not release blockers.

The result summary can be reproduced with:

```sh
cd mplapack/test
python3 ../../misc/summarize_mplapack_tests.py eig eig/results/3.0.1/* --only-fail
python3 ../../misc/summarize_mplapack_tests.py lin lin/results/3.0.1/* --only-fail
```

## Release QA

The release QA snapshot is in `release/logs/20260910_175329/`:

| QA area | Result |
| --- | --- |
| Tier1 build and test | 12/12 `OK` |
| Tier3 C++ standard build | 2/2 `OK` |
| Ubuntu 24.04 tarball smoke build | amd64 and arm64 completed |
| Tier2 branch matrix | 15 entries skipped; Ubuntu 22 MinGW64 and Debian 13 MinGW64 reference entries failed |

Tier2 is not a release-gating tier in the current project policy.  The
legacy skipped entries and the two MinGW reference-build failures remain
visible in `results.csv`; they are not represented as successful tests.

## Known limitations

- The known non-harmful eig threshold-edge results listed above remain.
- MinGW/Wine long-double trigonometric limitations remain documented for the
  affected wide-precision paths.
- Legacy Tier2 environments may lack a sufficiently new CMake or otherwise
  fall outside the supported release build matrix.
