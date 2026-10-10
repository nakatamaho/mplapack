# Adding an arithmetic backend

This is the recipe for adding a backend (a REAL/COMPLEX type pair such as
`td_real`/`td_complex`) to MPLAPACK, written as a sequence of
transformations from a model backend.  It was derived from the td backend
(3.1 branch, model dd, partly qd) and is meant to be replayed for ts, qs,
ds, edd and float.

Read `AGENTS.md` first.  Sources under `mpblas/reference/` and
`mplapack/reference/` come from the fable pipeline; both build systems
(autotools and CMake) must keep working.

## 0. Notation

| Symbol | Meaning | td |
|---|---|---|
| `B` | new backend name, `[a-z][a-z0-9]*` | `td` |
| `M` | model backend: the existing backend closest to `B` | `dd` (`qd` for lin/eig `Makefile.am`) |
| `B↑`, `M↑` | upper case | `TD`, `DD` |
| `REAL`, `COMPLEX` | C++ types | `td_real`, `td_complex` |
| `L` | number of limbs (expansion types) | 3 |
| `T` | `REAL` | |
| `R(M→B)` | `misc/backend_rename.py M B <file>` | |
| `SITES(M,B)` | `misc/backend_sites.py M B` | |

Model choice: dd, td or qd for a double expansion; td for ds/ts/qs
(`single_real<N>` follows the same libQD3 path); binary80 for edd where it
concerns platform conditions and dd for the expansion code; double for
float.

## 1. Operations

Every step below is one of these operations.

| Op | Input | Effect | Done when |
|---|---|---|---|
| `TABLE` | table file | add a row or section for `B` | row present |
| `GEN` | command | regenerate; never edit the output | `git diff` shows only `B` |
| `COPY` | `src` of `M` | `dst` := `R(M→B)(src)`, then the listed `VALUE` edits | file builds |
| `BLOCK` | `#if/#elif MPLAPACK_BUILD_WITH_M↑` alone | add a sibling `#elif defined MPLAPACK_BUILD_WITH_B↑` block := `R(M→B)(block)`, then `VALUE` edits | `SITES` no longer lists it |
| `JOIN` | `#if` testing `M↑` with other backends | `DECIDE`; if yes, add `\|\| defined MPLAPACK_BUILD_WITH_B↑` | decision recorded in section 5 |
| `LIST` | hand-written list of backend names | insert `B` next to `M` (or to the group `B` joins) | `SITES` no longer lists it |
| `PAIR` | fable patch + generated source | apply the same `BLOCK`/`JOIN` to `fable/3.12.1/lapack/patch-X.cpp` and to the generated `X.cpp` | `patch -R` check in section 5 |
| `VALUE` | per-type constant | value from section 3 | test in section 8 passes |
| `DECIDE` | choice | make it, record it in this file | recorded |

`R(M→B)` replaces `M`/`M↑` wherever it is a token (bounded by anything but a
letter or digit: `dd_real`, `Rlamch_dd`, `_dd_opt`, `Rgemm.dd`,
`optimized/dd/`, `ENABLE_DD`).  It never touches the libQD3 names that every
libQD3 backend shares (KEEP): `QD_LIBDIR`, `QD_INCLUDEDIR`,
`QD_RPATH_LDFLAGS`, `QD_EXTRA_CXXFLAGS`, `@qd_extra_cxxflags@`, `-lqd`,
`libqd`, `<qd/...>`, `include/qd`, and `DD_CXXFLAGS` (the no-FMA flags of
every `nofma` backend).  A plain `sed s/dd/td/` is wrong: it produced
`LDATD` from `LDADD` and `top_builtdir` from `top_builddir` during td.

## 2. Procedure

Run the steps in order.  `f` is a file, `<b>` stands for `B` in paths.

### S1. Table

| Op | File | Content (td) |
|---|---|---|
| `TABLE` | `backends.txt` | `td td_real td_complex OFF qd,nofma`; libQD3 expansions take `qd,nofma` |
| `TABLE` | `misc/backend_makefiles/backends.conf` | `[td]` `compare_tolerance = 1e-45` (section 3) |

### S2. Per-backend files

| Op | `dst` (from `M`'s file) | Residual edits after `R(M→B)` (td) |
|---|---|---|
| `COPY` | `mpblas/optimized/<b>/{Raxpy,Rcopy,Rdot,Rgemm}{,_ref}.cpp`, `openmp/*.cpp` | none |
| `COPY` | `mplapack/optimized/<b>/Makefile.am`, `.gitkeep` | none |
| `COPY` | `mpblas/reference/mpblas_<b>.h.in`, `mplapack/reference/mplapack_<b>.h.in`, `mplapack/test/{lin,eig}/common/mplapack_{lin,eig}_<b>.h.in`, `mplapack/test/matgen/mplapack_matgen_<b>.h.in` | none |
| `COPY` | `examples/{mpblas,mplapack}/generic/header_<b>{,_complex}` | none |
| `COPY` | `mplapack/test/{lin,eig}/<b>/Makefile.am`, model qd | add `$(DD_CXXFLAGS)` to every `*_CXXFLAGS` line (trait `nofma`); copy the current file (libraries in `LDADD`) |
| `COPY` | `include/mplapack_utils_<b>.h` | `B↑_PRECISION`; limb-wise hex printing (`L` limbs); `castREAL` sets limbs `2..L-1` to 0 |
| `COPY` | `include/mplapack_arithmetic_params_<b>.h` | `eps`, emin source, Blue exponents (section 3) |

### S3. Generated files

| Op | Command | Output |
|---|---|---|
| `GEN` | `python3 misc/gen_backend_makefiles.py` | `mpblas/optimized/<b>/Makefile.am`, `mpblas/test/<b>/Makefile.am`, `mplapack/test/compare/<b>/Makefile.am`, `benchmark/Makefile.<b>.am` |
| `GEN` | `bash fable/gen_include_{mpblas,mplapack,mplapack_lin,mplapack_eig,mplapack_matgen}.sh` (universal-ctags, `CLANG_FORMAT=clang-format-19`) | `include/{mpblas,mplapack,mplapack_lin,mplapack_eig,mplapack_matgen}_<b>.h`; review: other backends' headers must not change |
| `GEN` | `misc/gen_benchmark_header.py <b> libmplapack_<b>.a` (after S6 builds the library) | `include/mplapack_benchmark_<b>.h`; renaming another header is wrong (libQD3 complex types are `qd3_complex<T>` aliases) |
| `GEN` | run `Rlamch.test_<b>`, `Rlaruv.test_<b>` in `mplapack/test/compare/<b>` | `Rlamch_reference.txt`, `Rlaruv_reference.txt` (from the backend's own `Rlamch.txt`, `Rlaruv.txt`; review before commit) |
| `GEN` | `cd examples && bash gen_all.sh` (after S7 edits the templates) | `examples/**/*_<b>.cpp`, the example `Makefile.am`/`Makefile.*.in`; run once before editing to confirm it reproduces the tree |
| `GEN` | `./gen_configure.sh` | autotools files; fails if S3's makefiles are stale |

### S4. Shared sources: `SITES(M,B)` loop

Run `misc/backend_sites.py M B` and handle every line until only the
intended exclusions (section 6) remain.  Kind → operation:

| Kind | Operation |
|---|---|
| `PP-OWN` | `BLOCK` (or `DECIDE` not to follow `M`) |
| `PP-GROUP` | `JOIN` + `DECIDE` |
| `ENABLE` | `LIST` (configure, Makefile.am, packaging, release); configurations that disable `M` need nothing for a default-off `B` |
| `LIST` | `LIST` |
| `NAME` | an unguarded per-backend definition → `BLOCK` without `#if` (as `Mxerbla_<b>`); a file name list → `LIST` |
| `PLOT` | `DECIDE` which plot (S7) |

The files td needed, by operation (in `fable/3.12.1/lapack/patch-*.cpp` and
their generated sources, every change is `PAIR`):

| Op | Files |
|---|---|
| `BLOCK` | `include/{mpblas,mplapack,mplapack_arithmetic_params,mplapack_utils,mplapack_benchmark,mplapack_compare_debug,mplapack_lin,mplapack_eig,mplapack_matgen}.h`; `mpblas/reference/mplapackinit.cpp`; `mplapack/reference/{Rlamch,Rlaruv,Risnan,Risinf,Mexponent}.cpp`; `mpblas/test/common/{arithmetic,complex,mplapack}.test.cpp`; `mplapack/test/compare/common/Rlamch.test.cpp`; `fable/fem/{read,write}.hpp` |
| `BLOCK` (no `#if`) | `mpblas/test/common/Mxerbla.override.cpp` (`Mxerbla_<b>`) |
| `JOIN` + `PAIR` | `Cgejsv`, `Cgesvj`, `Rgejsv`, `Rgesvj`, `iMieeeck` (reference); `Rlatb4`, `Clatb4` (lin); `Rchkbd`, `Cchkbd` (eig); `Claror` (matgen) |
| `JOIN` | the `||` groups inside the `BLOCK` files above (for example the QD/DD/TD groups in `mplapack_compare_debug.h`, `mplapackinit.cpp`, `fem/*.hpp`) |
| `LIST` | `fable/convert_blas_all.sh`, `fable/convert_lapack_all.sh`, `fable/go_testing.sh` (`KEEP_HAND_WRITTEN_FILES`: the five `*_<b>.h.in`); `misc/check_backend_support_libraries.sh` (`<b> <b>_opt`, dependency case); `misc/check_test_override_link.sh` (second loop: backends with an external library) |

Scripts that read `backends.txt` and need nothing: CMake
(`cmake/MplapackBackendTable.cmake`), the fable header generators,
`misc/check_source_manifests.sh`, `misc/check_test_link_order.sh`, the first
loop of `misc/check_test_override_link.sh`,
`release/docker/common/compare-fable-outputs.py`;
`misc/check_backend_support_install.sh` takes the built
`mpblas/optimized/*/`.

### S5. Build wiring (`LIST`, by hand)

`SITES` lists these through `ENABLE`; the content to add:

| File | Content |
|---|---|
| `configure.ac` | `AC_ARG_ENABLE(<b>)`, `AM_CONDITIONAL(ENABLE_B↑)`; for trait `qd`: the `--with-system-qd` and `need_internal_qdlib` conditions; enable/summary messages; type detection if platform-dependent (as binary80); `AC_CONFIG_FILES` for `mpblas/optimized/<b>`, `mplapack/optimized/<b>`, `mpblas/test/<b>`, `mplapack/test/compare/<b>`, `mplapack/test/lin/<b>`, `mplapack/test/eig/<b>` |
| `Makefile.am` | unique-symbol archive list, `include/*_<b>.h`, `mplapack_<b>.pc`, `mplapack_<b>_opt.pc`, the libQD3 shared-load loop |
| `mplapack/reference/Makefile.am` | `libmplapack_<b>.la` (copy `M`'s block) |
| `mpblas/optimized/Makefile.am`, `mpblas/test/Makefile.am`, `mplapack/test/compare/Makefile.am`, `mplapack/optimized/Makefile.am` | `SUBDIRS += <b>` under `if ENABLE_B↑` |
| `mplapack/test/{lin,eig,matgen}/Makefile.am` | `lin_<b>`, `eig_<b>`, `matgen_<b>` libraries, `SUBDIRS`, `CHECK_BACKENDS` |
| `benchmark/Makefile.am` | `include Makefile.<b>.am` under `if ENABLE_B↑` |

### S6. Build, then the remaining `GEN` of S3

Build with autotools and CMake; generate the benchmark header and the
compare reference data.

### S7. Examples, benchmark scripts, packaging, documentation

| Op | File | Content |
|---|---|---|
| `LIST` | `examples/{mpblas,mplapack}/generic/generate.sh` | `MPLIBS` lists (mplapack: three, the `Cgeev_NPR` list without gmp, the `for _mplib` loop and its `case` → `B↑LIBS`); `REAL`/`COMPLEX` substitution branch; `if ENABLE_B↑` block for `Makefile.am`; `DISABLED_EXAMPLE_SUFFIXES` |
| `LIST` | `examples/{mpblas,mplapack}/generic/Makefile.{freebsd,linux,linux.inteloneAPI,macos,mingw}.in` | `B↑LIBS` (mpblas also `B↑OPTLIBS`), `programs=`, one rule per program; `Makefile.linux_cuda.in` is dd only |
| `LIST` | `examples/mplapack/run_smoke.sh` | `backends=` |
| `GEN` | `examples/gen_all.sh` | S3 |
| `LIST` | `benchmark/go.<Routine>.sh.in` | one of the two `for _mplib` loops |
| `DECIDE` + `LIST` | `benchmark/<Routine>{1,2}.plt.in` | plot 1 = fast types (binary80, binary128, dd), plot 2 = high precision (MPFR, GMP, qd); plain and `_opt` curve |
| `LIST` | `packaging/PKGBUILD.in`, `packaging/debian/rules`, `packaging/mplapack.spec.in` (+ `Provides:`), `misc/reconfig.*.sh`, `release/buildtest_tier1_macos_*.sh`, `release/docker/common/tarball-smoke.sh`, `release/docker/distcheck/*.sh`, `release/docker/matrix/*` | `--enable-<b>` (`=yes`) next to `--enable-M` |
| `LIST` | `README.md` | "Supported Precision Backends"; CMake option if default off |
| — | `MIGRATION.md` | nothing (breaking changes only) |
| — | `doc/manual/manual.tex` | not maintained (still 2.0.1, committed with its PDF) |

### S8. Release tests

The release QA (`release/README.md`: Tier 1 gating with `make distcheck`
and the full lin/eig suites, Tier 2 tarball and C++ standard builds, Tier 3
matrix) covers `B` only where `B` is enabled and only through checks that
know `B`.

| Op | File | Content |
|---|---|---|
| `LIST` | release configure options (S7) | Tier 1 `run-*-distcheck.sh` and macOS scripts, Tier 2 via `misc/reconfig.ubuntu26.04.*.sh`, Tier 3 `docker/matrix/*`, `tarball-smoke.sh` |
| `COPY` | `misc/backend_<b>_pkgconfig_consumer.cpp` from `M` | none |
| `LIST` | `cmake/tests/CMakeLists.txt` | consume source branch (`flavor MATCHES "^<b>"`), `mplapack_add_consume_tests(<b>)` and `(<b>_opt)`, `_pkgconfig_flavors`, the `(gmp\|qd\|dd\|...)` regex and its source branch |
| `LIST` | `Makefile.am` `EXTRA_DIST` | `misc/backend_<b>_pkgconfig_consumer.cpp` |
| `LIST` | `misc/check_backend_support_libraries.sh`, `misc/check_test_override_link.sh` | S4 |

Nothing per backend: `release/check-installed-examples.sh`,
`release/check-installed-benchmarks.sh`, `misc/summarize_mplapack_tests.py`
(backend from the directory name), the Tier 1 Fable reproduction (uses
`compare-fable-outputs.py`, which reads `backends.txt`).

## 3. Values

How each per-type constant is obtained.  `DERIVED`: computed from the type;
`MEASURED`: produced by running code; `CHOSEN`: picked, then validated by the
named test.

| Value | Where | Kind | Rule | dd | td | qd |
|---|---|---|---|---|---|---|
| `t` (digits) | arithmetic params | DERIVED | `std::numeric_limits<T>::digits` | 106 | 157 | 212 |
| `eps` | arithmetic params | DERIVED | `2^-t` as an exact literal (dd/qd keep historic QD literals) | — | `0x1p-157` | — |
| `prec` | arithmetic params | DERIVED | `eps * 2` | | | |
| `emin` | arithmetic params | DERIVED | exponent of `frexp(T::_min_normalized)`: all limbs stay normalized, so it is above the limb type's emin | -968 | -915 | -862 |
| `emax` | arithmetic params | DERIVED | `numeric_limits<limb>::max_exponent` | 1024 | 1024 | 1024 |
| Blue `exp_tsml` | arithmetic params | DERIVED | `ceil((emin-1)/2)` | | -458 | |
| Blue `exp_tbig` | arithmetic params | DERIVED | `floor((emax-t+1)/2)` | | 434 | |
| Blue `exp_ssml` | arithmetic params | DERIVED | `-floor((emin-t)/2)` | | 536 | |
| Blue `exp_sbig` | arithmetic params | DERIVED | `-ceil((emax+t-1)/2)` | | -590 | |
| `B↑_PRECISION` | `mplapack_utils_<b>.h` | DERIVED | `T::_ndigits - 1` | 30 | 46 | 61 |
| digits compared | `arithmetic.test.cpp` | DERIVED | `B↑_PRECISION` | 30 | 46 | 61 |
| write digits | `fable/fem/write.hpp` | DERIVED | `16·L` (limb decimal digits × L) | 32 | 48 | 64 |
| Rlaruv limb scales | `Rlaruv.cpp` | DERIVED | limb `k` = uniform × `2^(-53k)`, `k = 0..L-1` | | `1, 2^-53, 2^-106` | |
| cast to MPFR | `mplapack_compare_debug.h` (`cast2ref`) | DERIVED | `mpfr_set_d` limb 0, `mpfr_add_d` limbs `1..L-1` | | 3 limbs | |
| `compare_tolerance` | `backends.conf` | DERIVED | `1e-(15·L)` for double limbs | 1e-30 | 1e-45 | 1e-60 |
| `EPSILON` | `mplapack_compare_debug.h` | CHOSEN | about `compare_tolerance`; validated by the compare routine tests | 1e-31 | 1e-45 | 1e-61 |
| `EPSILON2`..`EPSILON12` | `mplapack_compare_debug.h` | CHOSEN | between dd and qd by digits; validated by the compare routine tests | | 1e-42 … 1e-36 | |
| random numbers | `mplapack.test.cpp` | — | the type's own generator (`tdrand()`) | `ddrand` | `tdrand` | `qdrand` |
| reference data | `compare/<b>/*_reference.txt` | MEASURED | S3 | | | |
| benchmark symbols | `mplapack_benchmark_<b>.h` | MEASURED | S3 | | | |

For float limbs (ds/ts/qs) the "16·L" and "15·L" rules become their float
analogues (limb decimal digits 7–8); measure them instead of assuming, as
emin rises with the number of limbs.

## 4. Fable `PAIR` check

When only `+` lines of a patch change (adding `B↑` to an `#if`), editing the
patch and the generated source identically is equivalent to regenerating.
Check: `patch -R -o out <source> < <patch>` on the edited pair must give the
same result as on the old pair.  `patch-Claror.cpp` contains its hunks
twice and already reverse-applies with offsets and a leftover line; that
predates td.  `fable/3.9.1/` is unused.

## 5. Decisions taken for td

| Site | Decision |
|---|---|
| `JOIN` sites that guard against non-IEEE arithmetic or the double exponent range (`Cgejsv`, `Cgesvj`, `Rgejsv`, `Rgesvj`, `iMieeeck`, `Rlatb4`, `Clatb4`, `Rchkbd`, `Claror`) | join, with dd/qd |
| `Cchkbd` singular value range | qd's narrower range `-(half*half)*log(ulp)`, not dd's |
| `Rget32`, `Rget34` (qd only) | do not follow |
| `Rchkee` threshold increase (dd only) | do not follow |
| benchmark loop / plot | high-precision loop; plot 2 before qd |
| CMake default | `OFF`; `--enable-td` default `no` |
| packaging, release, reconfig | enabled wherever dd is enabled |

Subset reference at the time of writing (`grep -lr MPLAPACK_BUILD_WITH_
mplapack/reference mplapack/test/*/common mplapack/test/matgen
fable/3.12.1/lapack` lists them):

| Subset | Files |
|---|---|
| binary80, binary128 | `Cgeev`, `Cgges`, `Cgges3`, `Clatrs`; tests `Cget23`, `Cqrt13`, `Cqrt15`, `Clatb4` |
| dd, td, qd | `Cgejsv`, `Cgesvj`, `Rgejsv`, `Rgesvj` (with gmp, mpfr), `iMieeeck` (with gmp); tests `Cchkbd`, `Rchkbd`, `Clatb4`, `Rlatb4` (dd, td); `Rget32`, `Rget34` (qd only); matgen `Claror` (dd, td) |
| gmp, mpfr | `Rlarrb`, `Rlarrd`, `Rlarrk`, `Rstebz`, `Cgges3`; tests `Rdrgev3`, `Rget38` |
| gmp only | `Cbdsqr`, `Rbdsqr`, `Rbdsvdx`; test `Cchkee` |
| binary80, double, gmp | tests `Cchktsqr`, `Rchktsqr` |
| binary80, dd, double, gmp, mpfr | test `Rchkee` |

`MPLAPACK_BUILD_WITH_MPFR` in `mpblas/test/common/*.test.cpp` and
`mplapack/test/compare/common/*.test.cpp` is the test oracle; no change.

## 6. Expected residual of `SITES(dd,td)`

After td, `misc/backend_sites.py dd td` prints only these; for a new backend
the residual should be of the same kind.

| Kind | Site | Why it stays |
|---|---|---|
| `PP-OWN` | `Cchkbd` (patch and source), `Rchkee` (patch and source) | section 5 |
| `PP-GROUP` | `Rchkee` gmp/dd group | section 5 |
| `ENABLE` | `misc/reconfig.macOS.clang.sanitize.sh`, `misc/reconfig.ubuntu24.04.clang.asan.sh`, `release/docker/common/check-fable-reproduction.sh` | they disable dd; td is off by default |
| `LIST` | `benchmark/go.*.sh.in` fast loop | td is in the other loop |
| `LIST` | fallback lists in `check_test_link_order.sh`, `check_test_override_link.sh`, `compare-fable-outputs.py` | used only without `backends.txt` |
| `LIST`, `NAME` | `examples/README` | historical text |
| `NAME` | `cmake/tests/consume/main.cpp` | single-backend consumer sample |
| `NAME` | `examples/mpblas/generic/Makefile.linux_cuda.in` | dd-only CUDA |
| `NAME` | comments in `include/mplapack_symbol_resolver.h`, `mplapack/reference/Rlamch.cpp` | comments |
| `PLOT` | `benchmark/*1.plt.in` | td is in plot 2 |

## 7. Failures seen during td

Each was found by a failing test or build, not by review; `SITES` now
reports all of them.

| Site | Symptom | Kind that reports it |
|---|---|---|
| `mpblas/test/common/Mxerbla.override.cpp` has no `Mxerbla_<b>` | invalid-argument tests reach the library `Mxerbla` and exit; 98/150 mpblas tests fail | `NAME` |
| `fable/fem/{read,write}.hpp` lack a branch | lin/eig/matgen test build stops at `#error`; library and mpblas tests build | `PP-OWN` |
| `KEEP_HAND_WRITTEN_FILES` in the fable scripts lack `*_<b>.h.in` | nothing fails until a fable re-conversion deletes them | `NAME` |
| lin/eig `Makefile.am` copied before master's link-order fix | old link order kept after the merge (new files, no conflict) | `misc/check_test_link_order.sh` |
| release QA checks with their own backend lists (`check_test_override_link.sh`, `check_backend_support_install.sh`, `compare-fable-outputs.py`, `check_backend_support_libraries.sh`) | td silently not checked | `LIST` |
| `cmake/tests/CMakeLists.txt` had no consume or pkg-config test for td | the installed td CMake package and `mplapack_td.pc` were never consumed | `LIST`, `NAME` (since `CMakeLists.txt` is scanned) |
| `Rlamch.test.cpp` td section left with dd names | compile error | `R(M→B)` applied by hand instead of the script |
| `sed s/dd/td/` on Makefiles | `LDATD`, `top_builtdir` | use `R(M→B)` |

## 8. Verification

| # | Command | Pass condition |
|---|---|---|
| V1 | `python3 misc/gen_backend_makefiles.py --check`; fable header generators | no diff for other backends |
| V2 | `misc/backend_sites.py M B` | only section 6-type residuals |
| V3 | CMake configure + build with `B` on and off | builds |
| V4 | `./gen_configure.sh`; configure with and without `--enable-<b>`; build | builds |
| V5 | `misc/check_unique_symbols.sh libmplapack_<b>.a libmplapack_<b>_opt.a` | no duplicates |
| V6 | `make check-all` in `mpblas/test/<b>` | all pass (td: 150/150) |
| V7 | `make check` in `mplapack/test/compare/<b>` | Rlamch, Rlaruv pass |
| V8 | `make check-all` there builds the routine tests but runs none (no log rule for `*.test_<b>`, every backend); run `for t in *.test_<b>; do ./$t > $t.out 2>&1 \|\| echo FAIL $t; done` | failures only where dd/qd fail too (Rpotri, Classq); these are ignored |
| V9 | `sh misc/check_test_link_order.sh .`; `misc/check_test_override_link.sh <build> .`; `misc/check_backend_support_libraries.sh <build> . <qd libdir>` | PASS |
| V10 | lin suite; compare pass lines with `mplapack/test/lin/results/` of the model | same pass lines |
| V11 | eig suite (hours; full machine, keep logs) | as V10 |
| V12 | examples: build with both systems; run `*_<b>` (bundled libQD3 needs `LD_LIBRARY_PATH`); `03_SymmetricEigenproblems/Rsyev*_Frank`, `Rsyevd_DingDong` sweep n=5..1000, use `-STARTN 5 -ENDN 200 -STEPN 15` | exit 0 (td: 67/67) |
| V13 | `python3 release/test-fable-comparison.py`; `python3 release/test-source-preparation.py` | OK / PASS |
| V14 | `make dist`; per-backend files of `M` and `B` in the tarball, compared after renaming | differ only by `M`-only files (dd: CUDA) |
| V15 | tarball smoke: unpack, `./configure --enable-<b>=yes` (others off), build | builds |
| V16 | CMake with tests: `ctest -R 'consume\|pkgconfig'` | `consume_*_<b>_*`, `pkgconfig_mplapack_<b>_*` pass |
| V17 | Fable reproduction as in Tier 1 (`release/docker/common/check-fable-reproduction.sh` steps: fresh clone, `fable/go.sh`, `fable/go_testing.sh`, `compare-fable-outputs.py <tree> <clone>`; needs universal-ctags, clang-format, GNU parallel, the LAPACK tarball) | no content difference |

## 9. Plan for the 3.1 backends

libQD3 1.6.0 or later provides all of them (the complex types for ds/ts/qs
first appear in 1.6.0).

1. **td**: done (model dd; lin/eig `Makefile.am` from qd).  eig suite (V11)
   still to be run on a full machine.
2. **ds, ts, qs** together, model td: all three are `single_real<N>` from
   one header.  Numerically the hardest: `float` exponent range, emin rising
   with `L`, so every `DERIVED` rule in section 3 that assumes double limbs
   must be re-derived or measured.  Mainly the CPU reference for GPU
   kernels.
3. **edd** (`_Float64x`/long double pairs), model dd for the expansion code
   and binary80 for the platform conditions: x86 only.  The bundled libQD3
   build disables edd on MinGW (`external/qd/Makefile.am`).
4. **float**, model double; can go in at any point.
