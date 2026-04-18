# Changes in MPLAPACK 2.1.1

Released: April 10, 2026

MPLAPACK 2.1.1 is a patch release on top of 2.1.0. It is ABI- and
API-compatible with 2.1.0. The release focuses on GCC 15 compatibility,
portability fixes across musl / MinGW / macOS arm64 / riscv64, a latent
DD-backend miscompilation risk, and an expanded build matrix. No public
BLAS / LAPACK routines were changed.

## Compiler and runtime fixes

* **GCC 15 miscompilation workaround in `{C,R}drgvx`.**
  `g++-15` at `-O2` (the test-driver default) miscompiles
  `mplapack/test/eig/common/{C,R}drgvx.cpp` via the new IPA modref
  pass, causing a segfault on the second explicit-example iteration.
  `-fno-ipa-modref` alone fixes it, but a per-function
  `optimize("no-ipa-modref")` attribute does not, because IPA passes
  are driven from the call site. The minimum effective workaround is
  `__attribute__((optimize("O1")))` on the `{C,R}drgvx` functions
  themselves. The attribute is unconditional; clang ignores unknown
  optimize strings and older GCC accepts `"O1"`. Only the test driver
  is affected; `g++-14 -O2` and `clang -O1` with ASan/UBSan both run
  cleanly.

* **`binary128`: unqualified `abs()` in `sign()`.**
  `std::abs(mplapack_binary128_t)` was ambiguous on platforms where
  `__float128` has no `std::abs` overload (e.g. MacPorts GCC 15 on
  `aarch64-apple-darwin`). Each `MPLAPACK_BINARY128_MATH` branch
  already defines `abs()` in the global namespace (`fabsq` for
  `MATH_QUADMATH`, `fabsf128` for `MATH_F128`, `using std::abs` for
  `MATH_LDBL`); using unqualified `abs()` lets normal lookup resolve
  the correct overload in all three cases.

* **MinGW build failure in `lapack/laed3`.**
  `common_interface.h` declares `slamc3` with return type `FLOATRET`
  (when `NEED_F2CCONV` is active), but `laed3_single.c` and
  `laed3_parallel.c` redeclared `LAMC3` as returning `FLOAT`, causing
  a `conflicting types for 'slamc3_'` error on MinGW. The local
  prototypes are now `FLOATRET`, matching the shared declaration.
  `laed3_parallel.c` also undefines the Windows `max` macro before its
  local `max` definition to avoid redefinition warnings.

* **Alpine / musl ARM64: fallback for `M_PIl`.**
  `M_PIl` is a glibc extension and is not provided by musl libc.
  A fallback definition is now emitted when system headers do not
  supply it, fixing the build on Alpine ARM64.

* **External GMP: GCC 15 (C23 default) build fix.**
  GCC 15 switches the default language standard to C23, which no
  longer accepts the untyped `void g(){}` in GMP's configure test.
  The upstream patch replacing it with a properly typed prototype is
  now applied to both `internal/` and `install/` GMP trees, but only
  when the C compiler reports major version >= 15 (detected at
  configure time as `GMP_NEEDS_GCC15_PATCH`). Older GCC releases skip
  the patch. The patch file (`gcc-15.patch`) is included in
  `EXTRA_DIST`.

* **External OpenBLAS bumped to 0.3.32.**
  Also fixes the external OpenBLAS build on macOS arm64.

* **`DD` backend: `-ffp-contract=off` (and `-fp-model strict` for
  icpx) is now propagated to all DD build targets.**
  The QD library itself was already built with `-ffp-contract=off`
  via `QD_EXTRA_CXXFLAGS`, but the MPLAPACK DD backend (reference
  libraries, optimized BLAS, tests, benchmarks) was not receiving the
  flag. This left a latent miscompilation risk on any target where
  the compiler enables FMA by default (x86-64 with AVX2/FMA3,
  AArch64, POWER). `DD_CXXFLAGS` is now defined from the already-
  probed `qd_extra_cxxflags` and appended to:
  `mpblas/optimized/dd`, `mpblas/reference`, `mplapack/reference`,
  `mpblas/test/dd`, `mplapack/test/compare/dd`,
  `mplapack/test/{eig,lin}/dd`, and `benchmark/Makefile.dd.am`.
  CUDA DD targets are intentionally excluded: `nvcc` does not accept
  `-ffp-contract=off`.

## Configure / build system

* **`configure.ac`: removed stale `CXXFLAGS="$SAVE_CXXFLAGS"`.**
  `SAVE_CXXFLAGS` was never assigned anywhere in `configure.ac`, so
  the restore line unconditionally expanded to `CXXFLAGS=""` and
  silently discarded any flags passed by the user via the environment
  (e.g. `-fsanitize=address`). The binary128 detection block already
  uses its own local `save_CXXFLAGS` pattern and cleans up before
  reaching this point, so the line was both incorrect and unnecessary.

* **`configure`: exact-overload-match probe for `std::abs(__float128)`.**
  The probe no longer relies on `__SIZEOF_FLOAT128__`. It now runs
  only when `__float128` support has already been detected, and
  requires an exact match for `&std::abs` with signature
  `__float128(__float128)`. This avoids reporting success when
  `std::abs(x)` merely happens to compile through an implicit
  conversion, and makes the probe reflect what MPLAPACK actually
  needs before suppressing the fallback helper.

* **`configure`: human-readable build summary.**
  At the end of configure, MPLAPACK now reports the C/C++ compiler
  commands and detected versions, target CPU, integer model, OpenMP
  status and inferred runtime, enabled backends, binary128 / binary80
  details (type, I/O, math functions, literal suffix, interop), and
  `std::abs(__float128)` availability.

* **`configure`: benchmark flag typo fix.**
  The internal LAPACK condition incorrectly used `enable_bench`
  instead of `enable_benchmark`.

* **`autoconf`: `external/openblas/Makefile` moved to the
  benchmark-conditional `AC_CONFIG_FILES` block.**
  It was unconditionally registered in `AC_CONFIG_FILES`, but
  `external/openblas` is only added to `SUBDIRS` when the benchmark
  is enabled. `make distcheck` (which runs `./configure` without
  options) generated the Makefile but never removed it in
  `make distclean`, failing the `distcleancheck` step. The Makefile
  is now generated conditionally, matching the actual `SUBDIRS`
  conditionality.

## Test harness

* **dd / qd comparison tests now use numerical tolerance.**
  The old scripts compared outputs with `diff` (plus brittle `sed`
  rewrites for qd), which failed on 1-bit discrepancies from
  hardware-dependent rounding such as FMA. The dd and qd test
  runners now invoke a new `misc/num_diff.py` helper with relative
  tolerance `1e-30` (dd) and `1e-60` (qd). Legacy `sed`-based output
  mangling was removed from the qd path. Tests fail only when the
  numerical difference exceeds the tolerance.

* **`mplapack/test/compare`: VPATH build fixes.**
  `run_Rlaruv_test.sh` / `run_Rlamch_test.sh` used `script_dir`
  (which resolves to the build directory) to locate
  `_reference.txt` files; under `make distcheck`'s out-of-tree
  build these files only exist in the source directory. The scripts
  now use `${srcdir}`, with the shebang switched to `/bin/sh` and
  `TOLERANCE` inlined. `NUM_DIFF_PY` in the dd/qd
  `Makefile.am.part` files now uses `$(top_srcdir)/misc/num_diff.py`
  instead of a fragile relative path, so the generated shell scripts
  carry the correct absolute path under distcheck.

* **GMP/MPFR/QD test drivers no longer link external LAPACK(E)/BLAS.**
  The external LAPACKE / LAPACK / BLAS / gfortran flags have been
  dropped from both the `eig` and `lin` GMP / MPFR / QD test
  drivers. These libraries were not needed by those test binaries
  and were masking the intended dependency structure.

## Distribution (what ships in the tarball)

* `misc/num_diff.py` is now shipped in release tarballs via
  `EXTRA_DIST` in the top-level `Makefile.am`.
* `gcc-15.patch` is now shipped in `EXTRA_DIST`.
* The following `misc/` helper scripts are now shipped in
  `EXTRA_DIST`:
  `buildtest_tier1_arm64_macOS.sh`,
  `buildtest_tier2_amd64_ubuntu_intel.sh`,
  `buildtest_tier2_i386_debian.sh`,
  `reconfig.bullseye.sh`,
  `reconfig.rocky10.sh`,
  `reconfig.ubuntu22.04.intel.sh`,
  `reconfig.ubuntu22.04.mingw64.sh`,
  `reconfig.ubuntu22.04.nvidia.sh`,
  `reconfig.ubuntu22.04.sh`,
  `reconfig.ubuntu24.04.intel.sh`.
* `misc/summarize_mplapack_tests.py` and related helpers restored
  Python 3.6 compatibility (PEP 585 `list[...]` annotations replaced
  with `typing.List`, `typing.Dict`, `typing.Tuple`), so they run on
  legacy distributions.

## Platform support

* **arm64 Ubuntu promoted from tier 2 to tier 1.**
  `buildtest_tier2_arm64_ubuntu.sh` renamed to
  `buildtest_tier1_arm64_ubuntu.sh`.

* **Release build matrix overhaul.**
  The distros and architectures covered by the release build matrix
  have been significantly expanded relative to 2.1.0:
  * Alpine: 3.19 → 3.22 / 3.23; riscv64 added.
  * Fedora: 42 and 43 added.
  * Ubuntu: 18.04 and 20.04 added; Intel oneAPI entries
    (`ubuntu22-intel`, `ubuntu24-intel`) promoted.
  * Debian: 11 added; riscv64 added to Debian 13; `debian13-mingw64`
    added.
  * openSUSE: Leap 15.6, 16.0, and Tumbleweed promoted from
    commented-out to active.
  * MinGW release image: base upgraded from Ubuntu 22.04 to 24.04;
    switched to the `x86_64-w64-mingw32 *-posix` toolchain; the
    native build triple is now detected dynamically; `make check`
    runs the generated `.exe` test binaries through Wine
    (`LOG_COMPILER=wine`).

* **`cpu_os_detection.sh`: Apple Silicon support for Linux
  containers.**
  The virtualization layer (Apple VZ) masks CPU part codes to
  `0x000`, making chip generation undetectable from inside a
  container. A new `CPU_MODEL_OVERRIDE` environment variable
  bypasses detection and passes the value through `normalize_cpu()`
  (e.g. `"Apple M4 Pro"` → `"Apple_M4_Pro"`). When `lscpu` returns
  Vendor ID `Apple` with an empty Model name, ISA flags are used as
  a best-effort fallback: `bf16=0, afp=0` → M1; `bf16=1, afp=0`
  → M2; `bf16=1, afp=1` → M3 or M4 (indistinguishable in a
  container, so `CPU_MODEL_OVERRIDE` is recommended when exact
  generation is required).

## Documentation

* **Manual bumped from 1.0.0 to 2.0.1.**
  `manual.tex` and `manual.bib` refreshed. Added EPS benchmark
  figures for Raxpy, Rdot, Rgemm, Rgemv, Rgetrf, Rpotrf, and Rsyrk.

* **New binary128 / binary80 type support matrix.**
  A new document was produced against a 77-configuration build
  matrix covering Alpine 3.19-3.23, Debian 11-13, Ubuntu 18.04-24.04,
  Rocky 8-10, Fedora 42/43, openSUSE Leap 15/16 and Tumbleweed, on
  x86_64, i386/i686, aarch64, ppc64le, s390x, riscv64, and mips64le,
  plus Intel oneAPI (icpx 2025.3.2) and MinGW-w64 cross targets.
  Key points for downstream users:
  * The rule for which binary128 path is selected is now
    `(GCC >= 13) AND (libc provides strfromf128)` — i.e. Alpine with
    GCC 14.2 / 15.2 stays on `__float128 + libquadmath` because musl
    does not ship `strfromf128` / `strfromf64x`.
  * RISC-V is the only non-x86 architecture that swings between
    combinations D and E on libc alone (musl vs. glibc), making it a
    useful regression target.
  * `mips64le` (Debian 12) is combination E. The configure summary
    reports CPU as `mips64` instead of `mips64le` / `mips64el`.
  * `ppc64le` long double is IBM double-double or IEEE-128 depending
    on ABI flags, so `ppc64le` never takes combination E.
  * Rocky 8/9, openSUSE Leap 15, and Ubuntu 18 / 20 / 22 remain on
    the libquadmath path for their full support window, so
    libquadmath support cannot be dropped.
  * A guard-expression appendix gives the correct
    architecture-first, libc-second, GCC-version-third ordering for
    conditional compilation.

* `README.md` updated for the 2.1.1 release (GCC 15 support,
  supported platforms).

## Version

* `MPLAPACK_VER_PATCH` bumped from 0 to 1 in `configure.ac` and
  `iMlaver.cpp`. The reported version is now 2.1.1.
