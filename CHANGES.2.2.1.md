# MPLAPACK 2.2.1 - Changes since 2.2.0

Release date: 2026-07-21 (base commit: `8fb8c83427c3561e003fab19202cd8c0f7833bbb`)

MPLAPACK 2.2.1 is a patch release on top of 2.2.0. It keeps the
LAPACK 3.12.1-based reference code from 2.2.0 and focuses on GMP
transcendental correctness, complex helper coverage, binary128/binary80
C++ standard compatibility, OpenCL binary128 work, an experimental CMake
build, and release engineering for the expanded build matrix.

This summary is based on 168 commits in `v2.2.0..8fb8c834`.

No public BLAS/LAPACK routine set is intentionally changed in this
release. The notable user-visible changes are in backend helper behavior,
headers, configure/CMake feature detection, examples, benchmarks, and
release/test infrastructure.

---

## 1. GMP backend transcendental functions

- Added native GMP-backed helpers for `log2`, `log1p`, `log`, `exp`,
  `expm1`, `sin`/`cos`, `atan`, `atan2`, and `pow`, with MPFR-free
  primary tests and an MPFR oracle test where appropriate
  (`1764beada` through `3265bda91`, `566680bd4`, `b1c515a88`).
- Added AGM-based constants and logarithms, plus a `quo_rem` argument
  reduction primitive used by the trigonometric path.
- Hardened the transcendental implementation after merge:
  denominator-overflow handling in `sincos`, tighter `quo_rem` precision
  handling, `expm1` cancellation fixes, and additional edge-case tests
  (`4e0325007`, `f346408eb`, `52e3f3aa7`, `d15693a9a`,
  `636d8f959`).
- Re-enabled the GMP transcendental path in `LATMS`/`LATMT`, enabled GMP
  CSD tests, and removed the temporary GMP-specific `BBCSD` tolerance
  override (`ef3cd6993`, `56285b344`, `4e9653233`).
- Fixed GMP singleton interval checks in `Rstebz` (`567344252`,
  `6c8c842f3`).
- Relaxed one GMP cosine-vs-MPFR comparison tolerance on i386
  (`41129b5bf`).

## 2. Multiprecision complex and predicate hardening

- Completed multiprecision complex helper coverage and preserved
  `mpreal` precision when copying (`4c58126d3`, `2d963aa24`).
- Added MPC pointer conversions for `mpcomplex` (`a823c3a19`).
- Applied safer `mpcomplex` constructor behavior (`f038f019e`).
- Fixed NaN and infinity predicates across generated headers
  (`65ac082fb`).
- Added `mpcomplex` overloads for binary complex addition and arithmetic
  with `std::complex<_Float128>` / `std::complex<__float128>` /
  `std::complex<_Float64x>` forms, avoiding ambiguous overloads in
  C++23-era builds (`4e7992b7c`, `b54ef2506`).

## 3. C++ standard and binary128/binary80 feature detection

- Added `--with-cxx-standard` support for C++17, C++20, C++23, and
  C++26, including GNU-extension and strict-standard modes where
  supported (`4957620ca`, `3d25873f7`).
- Hardened binary128 `_Float128` math selection:
  `std::abs`, scalar math (`sin`, `cos`, `atan2`, `exp`, `log`, etc.),
  `std::complex`, and C complex support are now controlled by configure
  probes rather than only by language-mode assumptions
  (`8f8cdf62e`, `f5de81ace`, `b0351c2d5`).
- Added analogous binary80 `_Float64x` math and `abs` probing
  (`b0351c2d5`).
- Extended configure summaries to report binary128/binary80 `std::abs`,
  `std::math`, `std::complex`, and C-complex availability, including the
  probed function groups.
- Added Ubuntu 26.04 C++ standard build scripts and Tier 3 C++ standard
  coverage for both amd64 and arm64 (`50a181857`, `e8d16eb99`,
  `dadfbb45a`).
- Documented the binary128/binary80 type-support and math-support
  situation in `binary128_binary80_type_support.md` (`b0351c2d5`,
  `2353d0352`).

## 4. Header include-order fixes

- `mplapack_binary128.h` now includes `<complex>` itself so users can
  include the installed LAPACK declarations directly (`0786e001b`).
- `mpblas_binary128.h` was hardened so `mplapack_config.h` appears before
  standard headers that may depend on feature-test macros
  (`2b92d5a21`).
- `mpblas_binary80.h` received the same include-order hardening
  (`308abcccd`).

## 5. Experimental CMake build

- Added an experimental modern CMake build system with FetchContent-based
  dependency support (`073560fa8`).
- Added install/package support, pkg-config/CMake package files, and CPack
  source/binary packaging as the CMake-side equivalent of `make dist`
  (`6d88fee92`).
- Ported the MPBLAS correctness tests and `mplapack/test/compare` suite to
  CMake, including link-group/whole-archive fixes for the test override
  libraries (`19a041488`, `7a69e0279`, `5f7ecb460`).
- Added CMake benchmark coverage for reference BLAS vs. OpenBLAS
  comparison builds (`38a302ddf`, `d0ded909d`).
- Added CMake-side binary128/binary80 math probes so the CMake config
  header mirrors the autotools probe results for the affected math and
  `abs` helpers (`2353d0352`).
- Documented the experimental CMake build in `README.md`; the autotools
  build remains the primary release path (`2353d0352`).

## 6. OpenCL binary128 work

- Merged and ported OpenCL support for binary128 RGEMM, including older
  OpenCL platform/device fixes and moving scalar parameters to the GPU
  path (`353b090b2`, `b88af945e`, `f1a7462be`, `dc079c806`).
- Removed obsolete `_Float128` OpenCL sources after the port
  (`2770aec9e`).
- Added binary128 OpenCL benchmark targets and an `Rgemm` OpenCL benchmark
  runner (`c46c82e18`, `caed88118`).
- Made OpenCL enablement depend on a usable OpenCL configuration and
  gated OpenCL binary128 in remote builds (`58e369c27`, `ba7fae1cf`).
- Used a portable standard-library include in the binary128 OpenCL path
  (`70ac89c78`).

## 7. External libraries and toolchain fixes

- Updated MPC to 1.4.1 and added checksum verification for external
  distfiles (`c3daa1f8d`, `bd4f4e25f`).
- Updated OpenBLAS to 0.3.33, then pruned obsolete external archives and
  retained the required OpenBLAS version in the external tree
  (`e64ecc16e`, `540fc0f9f`, `8ab77298f`, `8201aaabc`).
- Fixed Intel oneAPI OpenBLAS configure handling and switched oneAPI
  benchmark BLAS selection to MKL where appropriate (`ce438a2a6`,
  `0ce19284c`).
- Propagated quadmath flags to oneAPI examples and distributed Intel
  oneAPI example makefiles (`5cc3a32d3`, `1d1616c47`, `dd31b81ff`).
- Fixed macOS SDK Mach-O include handling with GNU C++ (`9d65b086d`).
- Avoided a glibc `rsqrt` conflict in CUDA builds and generated CUDA
  libtool objects without depending on system libtool
  (`59b1637c4`, `57593f553`).
- Enforced a single `-std=gnu++17` in the legacy configure path so stale
  environment flags do not accumulate conflicting `-std=` options
  (`ef772cc19`).

## 8. Release engineering

- Version bumped to 2.2.1 in `configure.ac` (`0c4c3c609`,
  `ee73cd2c6`).
- Release buildtests were moved to remote tier targets, with remote macOS
  and Ubuntu arm64 Docker coverage added (`1bf24c411`, `ee09860db`,
  `67f90f0c4`, `dd7f2ba7a`).
- Release target names now include distro versions and use `ref`
  terminology consistently (`cba67fb1d`, `2c5f98430`).
- Release jobs run by host and use stable host target files so SSH does
  not consume target lists (`d0f833661`, `65c8b01c3`, `4f4a2ed53`).
- Remote jobs are serialized per host, so `make -j tier1 tier2 tier3
  tarball` can prepare independent targets while only one job runs on a
  given host (`fe8f283be`).
- Release source snapshots are now preserved across tier runs, dist
  tarballs can be cached by source key, and fresh source snapshot
  generation was fixed (`af13a20b8`, `e8fc3c424`, `05ed567a1`).
- `make realclean` and remote cleanup were tightened to remove local and
  remote release logs, temporary workdirs, success/failure stamps, and
  QA-tagged Docker containers/images without deleting the local checkout
  (`d4fb6eaa7`, `8fb8c8342`).
- Release logging now records timings, failures, and ccache statistics
  (`a9fd658a5`, `0026d55d9`, `a4a94f266`).
- Release Docker matrix layout was refactored; tarball smoke tests now
  run on remote native Docker hosts and use xz artifacts
  (`3ce0deb13`, `42350f31f`, `0649360b6`).
- Installed examples and benchmark scripts are checked in Tier 1
  release testing, with example targets filtered by enabled backends
  (`6c2de4952`, `0915833f0`).
- `ccache` use was made more consistent in release scripts, including
  Fortran reconfig scripts (`afc842962`, `cc3ed8f9c`).

## 9. Build matrix and developer containers

- Added an Ubuntu 26.04 development container and subsequent tool updates,
  including CUDA/WineHQ setup, clang, and portable clang-format wrappers
  (`a5c804d00`, `cae44a64c`, `015612c82`, `9eb78d462`,
  `8820b5441`).
- Refined Tier 1 / Tier 2 / Tier 3 release coverage: Ubuntu amd64 NVIDIA
  distcheck, Ubuntu arm64 Docker, Debian i386 split image, MinGW updates,
  Intel oneAPI entries, and selected native benchmarks on non-x86
  (`bca934bf5`, `5f6a670d9`, `839dae437`, `3ab6e8e26`).
- Added branch-matrix and Dockerfile cleanups, including removal of old
  tarball Dockerfiles and stale release distcheck logs (`1b78b86e0`,
  `9084dede2`, `5e38066cb`).

## 10. Test results and documentation

- Added and reorganized 2.2.0 test result logs for
  `Core_i5-9600_ubuntu22_04_gcc-11_4_0` and updated README release news
  after the 2.2.0 release (`e2cc45c13`, `a7bf9090d`).
- Added Git checkout bootstrap requirements documentation (`f8decdeea`).
- Updated MPLAPACK 3.0 TODO notes (`c9957c07d`).
- Ignored new autotools, CMake, and release distcheck artifacts in
  `.gitignore` (`58a65ce5c`, `b12ae7ecc`).

## 11. Compatibility notes

- Re-run `configure` or CMake when changing the C++ language standard.
  The binary128/binary80 math helper decisions are probe results and can
  differ between `c++17`, `gnu++17`, `c++20`, `c++23`, and `c++26`.
- Code that directly includes installed binary128/binary80 BLAS/LAPACK
  headers no longer needs to include `<complex>` first.
- The CMake build is still experimental in 2.2.1. Use autotools for the
  primary release build and use the CMake path for targeted validation and
  developer testing.
