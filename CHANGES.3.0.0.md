# MPLAPACK 3.0.0 - Release notes

Release date: 2026-08-29 (release baseline: `eec51cd962d452563f6bc8703a441d0c3d410c88`)

This summary covers the 3.0.0 release baseline through
`eec51cd962d452563f6bc8703a441d0c3d410c88`.

MPLAPACK 3.0.0 completes the GMP/MPFR C++ wrapper migration to
`gmpfrxx_mkII`, updates the bundled high-precision dependencies, and adds
portability and release-test hardening for current C++, CMake, MinGW, macOS,
and wide-precision environments.

## Release artifact

The release archive is generated from the release baseline above.  The
published checksum files are the authoritative identifiers for the archive.

- Archive: `mplapack-3.0.0.tar.xz`
- Source commit: `eec51cd962d452563f6bc8703a441d0c3d410c88`
- SHA256: see `mplapack-3.0.0.tar.xz.sha256sum`
- Source state: clean (`MPLAPACK_SOURCE_DIRTY=no`)
- Dist-cache key: `git-eec51cd962d452563f6bc8703a441d0c3d410c88`

## Version and ABI

- Package version is now 3.0.0 in Autotools, CMake, and `iMlaver`.
- The shared-library ABI advances independently to libtool `3:0:0`.
- CMake shared libraries use version `3.0.0` and SONAME `3`.
- The ABI change is intentional: GMP and MPFR wrappers now use the
  `gmpfrxx_mkII` type and API family instead of the legacy in-tree `mpfrc++`
  implementation.
- Existing GMP/MPFR consumers using the old wrapper headers or old ABI must
  update their includes, type names, and link configuration.

## Dependency projects and licensing

- The GMP and MPFR C++ backend migration is based on the separately maintained
  [`gmpfrxx_mkII`](https://github.com/nakatamaho/gmpfrxx_mkII) project. Its
  wrapper code is distributed under the **BSD-2-Clause** license. This is the
  wrapper's license; the upstream GMP, MPFR, and MPC licenses remain
  applicable to those libraries.
- The QD backend is based on the separately maintained
  [`libQD3`](https://github.com/nakatamaho/libQD3) project, which retains the
  original DD/QD functionality and adds `td_real` and `edd_real` types.
- Although libQD3 provides `td_real` and `edd_real`, MPLAPACK 3.0.0 does not
  support them as MPLAPACK backends. MPLAPACK exposes the `dd_real` and
  `qd_real` backends only; `td_real` and `edd_real` are dependency-level
  capabilities and are not part of the MPLAPACK 3.0.0 API or test matrix.

## GMP and MPFR backend migration

- GMP and MPFR integration was migrated from the legacy `mpfrc++` code to the
  bundled `gmpfrxx_mkII` integration contract.
- The legacy `mpreal` and `mpcomplex` wrapper types were retired in favor of
  `mpfr_class` and `mpc_class`.
- On the GMP side, the former `mpc_class` complex wrapper was renamed to
  `mpfc_class`.
- The new wrapper uses expression templates for arithmetic expressions,
  reducing unnecessary intermediate objects and providing an expected
  performance improvement for high-precision compound expressions.
- The migration covers Autotools, CMake, benchmarks, examples, MPBLAS,
  MPLAPACK, compare tests, lin/eig tests, and installed headers.
- Backend-specific numeric, random-state, exponent, and binary-float adapter
  compatibility layers were added or centralized.
- GMP now honors the external default precision context instead of imposing a
  separate wrapper default.
- MPFR now preserves the default MPFR exponent range in normal builds.
  Test configurations explicitly apply the required emulation profiles when
  testing binary64, binary80, and binary128 behavior.
- MPFR result comparison and `Rlaruv` reference handling were regenerated for
  the new wrapper and random-state behavior.
- A GMP/MPFR result-diff helper was added for diagnosing backend drift.
- Binary128 adapter tests are guarded when the installed MPFR does not expose
  the optional float128 conversion API.
- GCC binary128 limit macros and non-reserved MPLAPACK configuration macro
  names are handled across supported compiler modes.

## gmpfrxx_mkII validation and packaging

- The repository now contains explicit P0, P1, and P1R migration gates,
  consumer validation, ABI inventories, baseline comparisons, and spike tests.
- Local compatibility patches cover typed integer accessors, MPFR exponent
  accessors, MPC/MPFC component accessors, semantic binary adapters, random
  state, pointer conversion, and integer power operations.
- The bundled `gmpfrxx_mkII` dependency was updated through the 1.2.x, 1.3.x,
  and 1.4.0 releases; MPLAPACK 3.0.0 uses 1.4.0.
- External CMake checks receive the same compiler, linker, and launcher flags
  as the parent build, including ccache where configured.
- MinGW dependency discovery, GCC-version parsing, checksum verification, and
  DLL execution through libtool/Wine wrappers were hardened.
- Bundled gmpfrxx_mkII headers are removed correctly during uninstall and
  distcheck.

## QD and libQD3 backend

- The bundled QD implementation migrated to
  [`libQD3`](https://github.com/nakatamaho/libQD3) and was updated through
  1.3.0, 1.3.1, and finally 1.4.0.
- QD patches add or harden integer arithmetic and comparison overloads for
  `dd_real`, `td_real`, `qd_real`, and `edd_real`.
- The `td_real` and `edd_real` capabilities described above belong to libQD3;
  they are not MPLAPACK-supported backends in this release.
- QD CMake gained an uninstall target, shared-library SONAME handling, and
  cross-configure support for MinGW.
- Runtime search paths for QD shared libraries are propagated to DD/QD test
  binaries and external QD builds; bundled QD now respects the MPLAPACK
  install prefix.
- FMA selection includes a residual correctness probe rather than relying
  only on compiler feature detection.
- x87 FPU-mode handling separates DD/TD/QD arithmetic from EDD's verified
  binary80 mode.
- MinGW/Wine EDD long-double trigonometric argument-reduction limitations are
  documented. DD/TD/QD tests remain supported; EDD trigonometric failures
  under Wine remain a known platform limitation.

## Public headers and C++ compatibility

- Public BLAS declarations include `<complex>` directly instead of relying on
  transitive includes.
- Compare-debug support includes the standard headers it uses.
- Strict C++ modes avoid exposing unsupported C complex binary128 fallbacks.
- FABLE I/O helpers now support guarded `long long` read/write overloads.
- `mplapackint` width selection is configurable in both Autotools and CMake.
- Binary80 and binary128 integer-exponent `pow` overloads were added for
  `mplapackint` exponents.
- The generated examples were adjusted for QD integer-conversion ambiguity,
  expression-template conditional-expression issues, and high-precision
  numeric literal handling.

## Examples

- The examples tree was reorganized into named problem categories, including
  general linear equations and precision-comparison examples.
- New generated examples cover generalized symmetric, generalized
  nonsymmetric, and generalized singular-value problems.
- The `90_PrecisionComparison` suite includes Hilbert, Pascal, Vandermonde,
  Wilkinson, Frank, Kahan, GEJSV/GESVD, and related examples.
- Generic example templates are generated for all supported backend flavors.
- Generated examples are self-contained and no longer depend on an implicit
  working-directory stream or missing wrapper include path.
- Example directories, generated Makefiles, install checks, and source
  manifests are included correctly in release tarballs.

## Build, CMake, and release tooling

- CMake is installed in tier1, tier2, and tarball smoke-test images so the
  vendored gmpfrxx_mkII CMake checks can run during distcheck.
- CMake checks propagate the parent compiler and ccache launcher settings.
- macOS distcheck uses GNU make where required by the generated dependency
  fragments.
- MinGW/GMP test execution uses libtool execute wrappers and the correct
  runtime DLL lookup path.
- Remote release scripts now handle MinGW compiler suffixes, compiler-version
  log namespaces, host locks, tarball reuse, and source/dist labels more
  consistently.
- Tier2 and tarball logs use explicit namespaces, and release cleanup removes
  local and remote build logs, temporary build directories, QA-tagged Docker
  images, and generated test-result directories.
- Installed example-directory checks and dist manifests were hardened after
  the examples reorganization.
- CMake and Autotools remain supported as parallel build systems.

## Validation and known limitations

- The preceding 3.0.0 QA snapshot was exercised through the Tier1, Tier2,
  Tier3, and tarball workflows, including Linux, macOS, MinGW/Wine, Intel
  oneAPI, CMake, and distcheck paths.  Tier1 covered 11 platform/compiler
  triplets; the tarball smoke tests passed on Ubuntu 24.04 amd64 and arm64.
- The release commit adds the bundled-QD internal `lib` install-directory
  fix.  Its targeted internal QD build and all 10 QD tests pass; the full
  Tier1 matrix was not rerun for this release-only build-system fix.
- The 3.0.0 `lin` result set passed all 126,515,712 recognized tests.
- The `eig` result set contains 30 threshold exceedances in 28 output records
  out of 116,811,306 recognized tests.  The recurring cases are GMP `Cgd`
  (1/1092 on all 11 triplets) and MPFR default `Rgg` (1/3120 on all 11
  triplets); the remaining cases are the documented double/binary80
  platform-specific differences.
- These eig cases are known non-harmful threshold-edge results and are not
  treated as release blockers.  They are retained in the result set for
  traceability and follow-up.
- Upstream PRs #1365 and #1366 are still pending independently of the
  MPLAPACK-side patches.
- The MinGW/Wine EDD long-double trigonometric limitation remains documented
  and is not treated as a general MPLAPACK algorithm failure.
- Cross-platform last-bit differences and DD `Claror` optimization sensitivity
  remain documented engineering limitations rather than ABI changes.
- Ubuntu 18 and Ubuntu 20 reference builds use CMake versions below the
  bundled libQD3 minimum (and Ubuntu 18 also lacks the modern `-S`/`-B`
  invocation).  These legacy Tier2 environments are outside the release
  support set and are not release blockers.

## 3.0.0 test result summary

The 3.0.0 result trees were summarized from:

- `mplapack/test/eig/results/3.0.0/`
- `mplapack/test/lin/results/3.0.0/`

The result set covers the following 11 platform/compiler triplets:

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

| Category | Triplets | Groups | `.out` files | Recognized tests | Failed tests | Failed records |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| eig | 11 | 96 | 4,032 | 116,811,306 | 30 | 28 |
| lin | 11 | 96 | 384 | 126,515,712 | 0 | 0 |

The `lin` tests passed for all 126,515,712 recognized tests.  The `eig`
results contain 30 threshold exceedances in 28 output records, corresponding
to approximately 0.000026% of the 116,811,306 recognized tests.

The `eig` failure records are:

| Triplet | Precision | File | Suite | Failed/Total | Fail rate |
| --- | --- | --- | --- | ---: | ---: |
| Apple_M4_macos26_gcc-15_2_0 | gmp | gmp/Cgd.out | ZGV drivers | 1/1,092 | 0.09% |
| Apple_M4_ubuntu24_04_gcc-13_3_0 | gmp | gmp/Cgd.out | ZGV drivers | 1/1,092 | 0.09% |
| Apple_M4_ubuntu26_04_gcc-15_2_0 | gmp | gmp/Cgd.out | ZGV drivers | 1/1,092 | 0.09% |
| Core_i7-6920HQ_macos15_gcc-15_2_0 | gmp | gmp/Cgd.out | ZGV drivers | 1/1,092 | 0.09% |
| Ryzen_Threadripper_3970X_debian12_gcc-12_2_0 | gmp | gmp/Cgd.out | ZGV drivers | 1/1,092 | 0.09% |
| Ryzen_Threadripper_3970X_debian13_gcc-14_2_0 | gmp | gmp/Cgd.out | ZGV drivers | 1/1,092 | 0.09% |
| Ryzen_Threadripper_3970X_ubuntu24_04_gcc-13_3_0 | gmp | gmp/Cgd.out | ZGV drivers | 1/1,092 | 0.09% |
| Ryzen_Threadripper_3970X_ubuntu24_04_icx-2026_1_1 | gmp | gmp/Cgd.out | ZGV drivers | 1/1,092 | 0.09% |
| Ryzen_Threadripper_3970X_ubuntu26_04_gcc-15_2_0 | gmp | gmp/Cgd.out | ZGV drivers | 1/1,092 | 0.09% |
| Ryzen_Threadripper_3970X_ubuntu26_04_icx-2026_1_1 | gmp | gmp/Cgd.out | ZGV drivers | 1/1,092 | 0.09% |
| Ryzen_Threadripper_3970X_windows_gcc-13_x_x | gmp | gmp/Cgd.out | ZGV drivers | 1/1,092 | 0.09% |
| Apple_M4_macos26_gcc-15_2_0 | mpfr (default) | mpfr/Rgg.default.out | DGG | 1/3,120 | 0.03% |
| Apple_M4_ubuntu24_04_gcc-13_3_0 | mpfr (default) | mpfr/Rgg.default.out | DGG | 1/3,120 | 0.03% |
| Apple_M4_ubuntu26_04_gcc-15_2_0 | mpfr (default) | mpfr/Rgg.default.out | DGG | 1/3,120 | 0.03% |
| Core_i7-6920HQ_macos15_gcc-15_2_0 | mpfr (default) | mpfr/Rgg.default.out | DGG | 1/3,120 | 0.03% |
| Ryzen_Threadripper_3970X_debian12_gcc-12_2_0 | mpfr (default) | mpfr/Rgg.default.out | DGG | 1/3,120 | 0.03% |
| Ryzen_Threadripper_3970X_debian13_gcc-14_2_0 | mpfr (default) | mpfr/Rgg.default.out | DGG | 1/3,120 | 0.03% |
| Ryzen_Threadripper_3970X_ubuntu24_04_gcc-13_3_0 | mpfr (default) | mpfr/Rgg.default.out | DGG | 1/3,120 | 0.03% |
| Ryzen_Threadripper_3970X_ubuntu24_04_icx-2026_1_1 | mpfr (default) | mpfr/Rgg.default.out | DGG | 1/3,120 | 0.03% |
| Ryzen_Threadripper_3970X_ubuntu26_04_gcc-15_2_0 | mpfr (default) | mpfr/Rgg.default.out | DGG | 1/3,120 | 0.03% |
| Ryzen_Threadripper_3970X_ubuntu26_04_icx-2026_1_1 | mpfr (default) | mpfr/Rgg.default.out | DGG | 1/3,120 | 0.03% |
| Ryzen_Threadripper_3970X_windows_gcc-13_x_x | mpfr (default) | mpfr/Rgg.default.out | DGG | 1/3,120 | 0.03% |
| Ryzen_Threadripper_3970X_debian12_gcc-12_2_0 | double | double/Rse2.out | DST | 1/4,440 | 0.02% |
| Ryzen_Threadripper_3970X_debian13_gcc-14_2_0 | double | double/Rse2.out | DST | 1/4,440 | 0.02% |
| Core_i7-6920HQ_macos15_gcc-15_2_0 | binary80 | binary80/Rsep.out | DST drivers | 2/13,464 | 0.01% |
| Ryzen_Threadripper_3970X_windows_gcc-13_x_x | binary80 | binary80/Rsep.out | DST drivers | 2/13,464 | 0.01% |
| Ryzen_Threadripper_3970X_debian12_gcc-12_2_0 | double | double/Rsvd.out | DBD | 1/10,260 | 0.01% |
| Ryzen_Threadripper_3970X_debian13_gcc-14_2_0 | double | double/Rsvd.out | DBD | 1/10,260 | 0.01% |

The 11 GMP `Cgd.out` records and the 11 MPFR default
`Rgg.default.out` records are identical recurring one-test threshold
exceedances across every triplet.  The remaining records are the documented
double/binary80 platform-specific differences and are retained for
traceability.

The result summary can be reproduced with:

```sh
cd mplapack/test
python3 ../../misc/summarize_mplapack_tests.py eig eig/results/3.0.0/* --only-fail
python3 ../../misc/summarize_mplapack_tests.py lin lin/results/3.0.0/* --only-fail
```
