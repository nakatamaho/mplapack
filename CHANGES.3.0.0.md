# MPLAPACK 3.0.0 - Changes since 2.3.0

Release date: TBD (review draft; base commit: `a67000765864c8bcddb3a2017c20ca73f1869e1a`)

This summary is based on 96 commits in `v2.3.0..a67000765864c8bcddb3a2017c20ca73f1869e1a`.

MPLAPACK 3.0.0 completes the GMP/MPFR C++ wrapper migration to
`gmpfrxx_mkII`, updates the bundled high-precision dependencies, and adds
portability and release-test hardening for current C++, CMake, MinGW, macOS,
and wide-precision environments.

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

- The 3.0.0 release candidate was exercised through the tier1, tier2, tier3,
  and tarball workflows, including Linux, macOS, MinGW/Wine, Intel oneAPI,
  CMake, and distcheck paths.
- Upstream PRs #1365 and #1366 are still pending independently of the
  MPLAPACK-side patches.
- The MinGW/Wine EDD long-double trigonometric limitation remains documented
  and is not treated as a general MPLAPACK algorithm failure.
- Cross-platform last-bit differences and DD `Claror` optimization sensitivity
  remain documented engineering limitations rather than ABI changes.
