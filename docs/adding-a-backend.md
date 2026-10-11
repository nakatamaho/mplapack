# Adding an arithmetic backend

How to add a backend (a REAL/COMPLEX type pair) to MPLAPACK, written from
the td backend (`td_real`/`td_complex`, libQD3 1.6.0, branch 3.1; commits
745aaa91a … d2d0e43b2).  It is meant to be replayed for ds, ts, qs, edd and
float.

Read `AGENTS.md` first: sources under `mpblas/reference/` and
`mplapack/reference/` come from the fable pipeline, and autotools and CMake
must both keep working.

Contents

1. Workflow and tools
2. Files, by kind of change
3. Magic numbers: where each value comes from
4. Fable: where it gets stuck
5. Tests
6. Other precisions (ds, ts, qs, edd, float)
7. Record of td: decisions, failures, residual

---

## 1. Workflow and tools

Notation: `B` new backend (`td`), `M` model backend (`dd`), `B↑`/`M↑` upper
case, `L` number of limbs, `<b>` = `B` in paths.

| Step | What | Section |
|---|---|---|
| 1 | Add the table row; let generators produce what they can | 2.1 |
| 2 | Copy the per-backend files from `M` with the rename tool | 2.2 |
| 3 | Set the per-type numbers | 3 |
| 4 | Handle every site the site lister reports (shared sources, fable special cases, build wiring, scripts) | 2.3–2.6, 4 |
| 5 | Build with autotools and CMake; generate benchmark header and compare reference data | 2.1 |
| 6 | Tests | 5 |
| 7 | Examples, benchmark scripts, packaging, release, README | 2.7 |

Tools (git tree only, not shipped):

| Tool | Use |
|---|---|
| `misc/backend_rename.py M B <file>` | Copy-rename.  Replaces `M`/`M↑` where it is a token (`dd_real`, `Rlamch_dd`, `_dd_opt`, `Rgemm.dd`, `optimized/dd/`, `ENABLE_DD`), never inside words (`LDADD`, `top_builddir` stay intact), and never the names all libQD3 backends share: `QD_LIBDIR`, `QD_INCLUDEDIR`, `QD_RPATH_LDFLAGS`, `QD_EXTRA_CXXFLAGS`, `@qd_extra_cxxflags@`, `-lqd`, `libqd`, `<qd/…>`, `include/qd`, `DD_CXXFLAGS` (the no-FMA flags of every `nofma` backend).  `sed s/dd/td/` is wrong: it produced `LDATD` and `top_builtdir` during td. |
| `misc/backend_sites.py M B` | To-do list: every place that names `M` but not `B`.  Kinds: `PP-OWN` (an `#if` for `M` alone → add a block for `B`), `PP-GROUP` (`M` in an `\|\|` group → decide whether `B` joins), `ENABLE` (`ENABLE_M`, `--enable-m` lists), `LIST` (`M` in a list of backend names), `NAME` (an unguarded `…_m` identifier such as `Mxerbla_dd`), `PLOT` (benchmark plot).  Run it on the tree before td, it reports every shared file td had to change, including the ones found only by failing tests.  Repeat until only intended differences remain (section 7.3). |
| `misc/gen_backend_makefiles.py` | Per-backend `Makefile.am` files from `backends.txt` + `misc/backend_makefiles/` |
| `fable/gen_include_*.sh` | Prototype headers from `backends.txt` + `.h.in` |
| `misc/gen_benchmark_header.py B libmplapack_<b>.a` | Benchmark symbol header from the built library |
| `examples/gen_all.sh` | Example sources and makefiles |

Model choice: dd/td/qd for a double expansion; td for ds/ts/qs; dd for the
edd expansion code and binary80 for its platform conditions; double for
float.  td used dd, except `mplapack/test/{lin,eig}/<b>/Makefile.am` (from
qd, see 2.2).

## 2. Files, by kind of change

### 2.1 Table and generated files

| Kind | File / command | td |
|---|---|---|
| table | `backends.txt` | `td td_real td_complex OFF qd,nofma` |
| table | `misc/backend_makefiles/backends.conf` | `[td] compare_tolerance = 1e-45` |
| generated | `python3 misc/gen_backend_makefiles.py` → `mpblas/optimized/<b>/Makefile.am`, `mpblas/test/<b>/Makefile.am`, `mplapack/test/compare/<b>/Makefile.am`, `benchmark/Makefile.<b>.am` | |
| generated | `fable/gen_include_{mpblas,mplapack,mplapack_lin,mplapack_eig,mplapack_matgen}.sh` → `include/*_<b>.h` | equal to dd's headers with names replaced (section 4) |
| generated | `misc/gen_benchmark_header.py` → `include/mplapack_benchmark_<b>.h` | after the library builds; renaming dd's header gives wrong names (`td_complex` is `qd3_complex<td_real>`) |
| generated | `Rlamch.test_<b>`, `Rlaruv.test_<b>` → `mplapack/test/compare/<b>/{Rlamch,Rlaruv}_reference.txt` | review before committing (section 5.2) |
| generated | `examples/gen_all.sh` | section 2.7 |
| generated | `./gen_configure.sh` | fails when the generated makefiles are stale |
| automatic | CMake (`cmake/MplapackBackendTable.cmake`): option, macro, libraries, QD link, `-ffp-contract=off` from the traits | nothing to edit |

### 2.2 Per-backend files (copy from `M` with `backend_rename.py`)

| Files | Edits after the rename (td) |
|---|---|
| `mpblas/optimized/<b>/{Raxpy,Rcopy,Rdot,Rgemm}{,_ref}.cpp`, `openmp/*.cpp` | none |
| `mplapack/optimized/<b>/Makefile.am`, `.gitkeep` | none |
| `mpblas/reference/mpblas_<b>.h.in`, `mplapack/reference/mplapack_<b>.h.in`, `mplapack/test/{lin,eig}/common/mplapack_{lin,eig}_<b>.h.in`, `mplapack/test/matgen/mplapack_matgen_<b>.h.in` | none |
| `examples/{mpblas,mplapack}/generic/header_<b>{,_complex}` | none |
| `misc/backend_<b>_pkgconfig_consumer.cpp` | none |
| `mplapack/test/{lin,eig}/<b>/Makefile.am`, from qd | add `$(DD_CXXFLAGS)` to every `*_CXXFLAGS` line; copy the current file (libraries in `LDADD`, after the objects) |
| `include/mplapack_utils_<b>.h` | `B↑_PRECISION`; hex printing of `L` limbs; `castREAL` zeroes limbs `1..L-1` |
| `include/mplapack_arithmetic_params_<b>.h` | section 3 |

### 2.3 Shared sources: one `#if` block per backend (`PP-OWN`)

Add a block for `B` next to `M`'s, made with the rename tool, then set its
numbers (section 3).

| Area | Files |
|---|---|
| public headers | `include/{mpblas,mplapack,mplapack_arithmetic_params,mplapack_utils,mplapack_benchmark,mplapack_lin,mplapack_eig,mplapack_matgen}.h` |
| library, MPLAPACK's own code | `mpblas/reference/mplapackinit.cpp` (x87 precision), `mplapack/reference/{Rlamch,Rlaruv,Risnan,Risinf,Mexponent}.cpp` |
| test support | `include/mplapack_compare_debug.h` (EPSILON, cast to MPFR), `mpblas/test/common/{arithmetic,complex,mplapack}.test.cpp`, `mplapack/test/compare/common/Rlamch.test.cpp` |
| Fortran I/O emulation | `fable/fem/read.hpp`, `fable/fem/write.hpp` |
| unguarded (`NAME`) | `mpblas/test/common/Mxerbla.override.cpp`: add `Mxerbla_<b>` (no `#if`) |

Inside the same files, `||` groups (`QD/DD/TD …`) take `B` where `B` behaves
like `M` (`PP-GROUP`).

### 2.4 Fable special cases (`PP-GROUP` in generated code)

Edit the patch `fable/3.12.1/lapack/patch-X.cpp` and the generated `X.cpp`
the same way; section 4.

| td joined | Files |
|---|---|
| dd/qd | `Cgejsv`, `Cgesvj`, `Rgejsv`, `Rgesvj`, `iMieeeck` (reference); `Rchkbd` (eig) |
| dd | `Rlatb4`, `Clatb4` (lin); `Claror` (matgen) |
| qd's branch | `Cchkbd` (eig) |

### 2.5 Build wiring (`ENABLE`, by hand)

| File | Content |
|---|---|
| `configure.ac` | `AC_ARG_ENABLE(<b>)` (td: default no), `AM_CONDITIONAL(ENABLE_B↑)`; trait `qd`: the `--with-system-qd` and `need_internal_qdlib` conditions; messages; platform detection if needed (binary80 style); `AC_CONFIG_FILES` for `mpblas/optimized/<b>`, `mplapack/optimized/<b>`, `mpblas/test/<b>`, `mplapack/test/compare/<b>`, `mplapack/test/lin/<b>`, `mplapack/test/eig/<b>` |
| `Makefile.am` | unique-symbol archives, `include/*_<b>.h`, `mplapack_<b>{,_opt}.pc`, libQD3 shared-load loop, `EXTRA_DIST` (consumer source) |
| `mplapack/reference/Makefile.am` | `libmplapack_<b>.la` with `$(DD_CXXFLAGS)` |
| `mpblas/optimized/Makefile.am`, `mpblas/test/Makefile.am`, `mplapack/test/compare/Makefile.am`, `mplapack/optimized/Makefile.am` | `SUBDIRS` |
| `mplapack/test/{lin,eig,matgen}/Makefile.am` | `lin_<b>`, `eig_<b>`, `matgen_<b>`, `SUBDIRS`, `CHECK_BACKENDS` |
| `benchmark/Makefile.am` | `include Makefile.<b>.am` |

### 2.6 Scripts with their own backend lists (`LIST`)

| File | Content |
|---|---|
| `fable/convert_blas_all.sh`, `fable/convert_lapack_all.sh`, `fable/go_testing.sh` | `KEEP_HAND_WRITTEN_FILES`: the five `*_<b>.h.in` |
| `cmake/tests/CMakeLists.txt` | consume source branch, `mplapack_add_consume_tests(<b>)`/`(<b>_opt)`, `_pkgconfig_flavors`, the `(gmp\|qd\|dd\|td)` regex and its branch |
| `misc/check_backend_support_libraries.sh` | `<b> <b>_opt`, dependency case |
| `misc/check_test_override_link.sh` | second loop (backends whose `.pc` pulls an external library) |

These read `backends.txt` and need nothing: CMake, the fable header
generators, `misc/check_source_manifests.sh`, `misc/check_test_link_order.sh`,
the first loop of `check_test_override_link.sh`,
`release/docker/common/compare-fable-outputs.py`;
`misc/check_backend_support_install.sh` takes the built `mpblas/optimized/*/`.

### 2.7 Examples, benchmark scripts, packaging, release, documentation

| Area | File | Content |
|---|---|---|
| examples | `examples/{mpblas,mplapack}/generic/generate.sh` | `MPLIBS` lists (mplapack: three, the `Cgeev_NPR` list without gmp, the `for _mplib` loop and `case` → `B↑LIBS`), `REAL`/`COMPLEX` substitution branch, `if ENABLE_B↑` block, `DISABLED_EXAMPLE_SUFFIXES` |
| examples | `generic/Makefile.{freebsd,linux,linux.inteloneAPI,macos,mingw}.in` | `B↑LIBS` (mpblas: `B↑OPTLIBS`), `programs=`, rules; not `Makefile.linux_cuda.in` (dd only) |
| examples | `examples/mplapack/run_smoke.sh` | `backends=` |
| examples | `cd examples && bash gen_all.sh` | run it once before editing: it reproduces the tree, so afterwards `git diff` shows only `B` |
| benchmark | `benchmark/go.<Routine>.sh.in` | one of the two `for _mplib` loops |
| benchmark | `benchmark/<Routine>{1,2}.plt.in` | plot 1 fast types (binary80, binary128, dd), plot 2 high precision (MPFR, GMP, qd); td: plot 2 before qd |
| packaging | `packaging/PKGBUILD.in`, `packaging/debian/rules`, `packaging/mplapack.spec.in` (+ `Provides:`) | `--enable-<b>` |
| release | `misc/reconfig.*.sh` (Tier 2 uses the `ubuntu26.04.*` ones), `release/buildtest_tier1_macos_*.sh`, `release/docker/distcheck/run-*-distcheck.sh` (Tier 1), `release/docker/common/tarball-smoke.sh`, `release/docker/matrix/*` (Tier 3) | `--enable-<b>=yes` next to `--enable-dd=yes`; configurations with `--enable-dd=no` need nothing for a default-off backend |
| docs | `README.md` | "Supported Precision Backends", CMake option if default off |
| docs | `MIGRATION.md` | nothing (breaking changes only) |
| docs | `doc/manual/manual.tex` | not maintained (2.0.1, committed with its PDF) |

## 3. Magic numbers: where each value comes from

Kinds: **library** (read from the type at compile or run time; nothing to
choose), **derived** (formula from the type's constants), **chosen**
(picked, then validated by the named test), **measured** (output of running
code).

libQD3 constants (printed from the bundled libQD3 1.6.0):

| | dd | td | qd | ds | ts | qs | edd |
|---|---|---|---|---|---|---|---|
| limb type, its digits / max exponent | double 53/1024 | double 53/1024 | double 53/1024 | float 24/128 | float 24/128 | float 24/128 | long double 64/16384 |
| `L` | 2 | 3 | 4 | 2 | 3 | 4 | 2 |
| `numeric_limits<T>::digits` | 104 | 157 | 209 | 46 | 70 | 94 | 126 |
| `T::_eps` | 2^-104 | 2^-157 | 2^-209 | 2^-46 | 2^-70 | 2^-94 | 2^-126 |
| `T::_min_normalized` | 2^-969 | 2^-916 | 2^-863 | 2^-102 | 2^-78 | 2^-54 | 2^-16318 |
| `T::_ndigits` | 31 | 47 | 62 | 14 | 21 | 28 | 38 |

`_min_normalized = 2^(limb emin + limb digits × (L-1))` (limb emin: -1022
double, -126 float, -16382 long double): every limb must stay normalized, so
each extra limb raises the minimum by one limb's digits.

| Value | File | Kind | Rule | dd | td | qd |
|---|---|---|---|---|---|---|
| `t` | `arithmetic_params_<b>.h` | library | `numeric_limits<T>::digits` | 104 | 157 | 209 |
| `eps` (`E`) | `arithmetic_params_<b>.h`, `Rlamch.test.cpp` | derived | exact literal `2^-t`; dd/qd keep the historic QD literals (`0x1.fff…p-105`, `…p-209`) | | `0x1p-157` | |
| `prec` (`P`) | same | derived | `2·E` | | `0x1p-156` | |
| `emin` | `arithmetic_params_<b>.h` | library | exponent of `frexp(T::_min_normalized)` | -968 | -915 | -862 |
| `emax` | same | library | `numeric_limits<limb>::max_exponent` | 1024 | 1024 | 1024 |
| Blue `exp_tsml` | same | derived | `ceil((emin-1)/2)` | | -458 | |
| Blue `exp_tbig` | same | derived | `floor((emax-t+1)/2)` | | 434 | |
| Blue `exp_ssml` | same | derived | `-floor((emin-t)/2)` | | 536 | |
| Blue `exp_sbig` | same | derived | `-ceil((emax+t-1)/2)` | | -590 | |
| `B↑_PRECISION` | `mplapack_utils_<b>.h` | derived | `T::_ndigits - 1` | 30 | 46 | 61 |
| compared digits | `arithmetic.test.cpp` | derived | `B↑_PRECISION` | 30 | 46 | 61 |
| `B↑_PRECISION_SHORT`, `MPLAPACK_BUFLEN` | `mplapack_utils_<b>.h` | copied | 16, 1024 (all libQD3 backends) | | | |
| hex print buffers | `mplapack_utils_<b>.h` | chosen | 128 bytes per limb | | | |
| write digits | `fable/fem/write.hpp` | derived | `16·L` for double limbs (`≥ _ndigits+1`) | 32 | 48 | 64 |
| Rlaruv limb scales | `Rlaruv.cpp` | derived | limb `k` = uniform × `2^(-53k)` (53 = limb digits) | | 1, 2^-53, 2^-106 | |
| cast to MPFR | `mplapack_compare_debug.h` | derived | `mpfr_set_d` limb 0, `mpfr_add_d` limbs 1..L-1 | | 3 limbs | |
| random numbers | `mplapack.test.cpp` | library | the type's generator | `ddrand` | `tdrand` | `qdrand` |
| `compare_tolerance` | `backends.conf` | chosen | `1e-(15·L)` for double limbs, about `10^-(_ndigits-1)`; validated by compare Rlamch/Rlaruv | 1e-30 | 1e-45 | 1e-60 |
| `EPSILON` | `mplapack_compare_debug.h` | chosen | about `compare_tolerance` | 1e-31 | 1e-45 | 1e-61 |
| `EPSILON2..12` | `mplapack_compare_debug.h` | chosen | dd's ladder moved to start at 1e-45 with a smaller step (td: 1e-42, 1e-42, 1e-41, 1e-40, 1e-39, 1e-39, 1e-38, 1e-37, 1e-36; dd: 1e-27 … 1e-20); validated by the compare routine tests (103/105, the two failures also fail for dd/qd) | | | |
| `Cchkbd` singular value range | `Cchkbd.cpp` (patch) | chosen | qd's `-(half*half)*log(ulp)` | | | |
| reference data | `compare/<b>/*_reference.txt` | measured | section 5.2 | | | |
| benchmark symbols | `mplapack_benchmark_<b>.h` | measured | section 2.1 | | | |

The **chosen** values are the ones to revisit for a new type.  Derive the
**library**/**derived** ones from the table above; never copy them from the
model.

## 4. Fable: where it gets stuck

| # | Point | What happens | What to do |
|---|---|---|---|
| F1 | Prototype headers | `fable/gen_include_*.sh` needs universal-ctags and clang-format 19+ (`CLANG_FORMAT=clang-format-19`), and rewrites the headers of every backend | review the whole diff: a different ctags can add or drop prototypes of existing backends.  td's headers came out equal to dd's with names replaced, in fable's sort order. |
| F2 | `.h.in` templates | `convert_blas_all.sh`, `convert_lapack_all.sh`, `go_testing.sh` delete generated files except `KEEP_HAND_WRITTEN_FILES` | add the five `*_<b>.h.in`, or the next conversion deletes them (nothing fails before that) |
| F3 | Special cases in generated code | `#if` in `mplapack/reference/*.cpp` and the test sources come from `fable/3.12.1/lapack/patch-*.cpp` | never edit the generated file alone.  When only `+` lines change (adding `B↑` to a condition), edit patch and source identically; check that `patch -R -o out <source> < <patch>` gives the same result as on the old pair |
| F4 | dead text in a patch | `patch-Claror.cpp` carried a stale copy of its last three hunks after a blank line; `patch` stops at the blank line, so the copy was never applied but made reverse checks noisy | removed (master, 2.4, 3.1).  Any `@@` line right after a blank line in a patch is dead text: `awk 'p=="" && /^@@/' patch-*.cpp` |
| F5 | Fortran I/O emulation | `fable/fem/{read,write}.hpp` end in `#error` for an unknown backend | add the branch; the library and mpblas tests build without it, lin/eig/matgen stop |
| F6 | Release reproduction (Tier 1, Ubuntu 26.04 amd64) | `check-fable-reproduction.sh` regenerates everything from a fresh checkout and compares with the tree | run it before release: fresh clone, `autoreconf`, configure (double only), `bash fable/go.sh`, `ROOT=$PWD bash fable/go_testing.sh`, `python3 release/docker/common/compare-fable-outputs.py <tree> <clone>`; needs GNU parallel and `external/lapack/download/lapack-3.12.1.tar.gz`.  td: 1763 files, 0 differences (ctags 5.9.0, clang-format 19). |
| F7 | `fable/3.9.1/` | old LAPACK patches | unused; ignore |

## 5. Tests

### 5.1 What each suite needs and how to run it

| Suite | Needs for `B` | Run | td result |
|---|---|---|---|
| mpblas | generated `mpblas/test/<b>/Makefile.am`; `Mxerbla_<b>` in `Mxerbla.override.cpp`; branches in `arithmetic/complex/mplapack.test.cpp`; `EPSILON`, cast in `compare_debug.h` | `make check-all` in `mpblas/test/<b>` (`make check` runs only the smoke tests) | 150/150 |
| compare | generated `Makefile.am`; `compare_tolerance`; `Rlamch.test.cpp` branch; reference data | `make check` runs Rlamch/Rlaruv only; `make check-all` builds but runs no routine test (no log rule for `*.test_<b>`, every backend), so `for t in *.test_<b>; do ./$t > $t.out 2>&1 \|\| echo FAIL $t; done` | Rlamch, Rlaruv pass; routines 103/105 (Rpotri, Classq; they also fail for dd/qd and are ignored) |
| lin | `lin/<b>/Makefile.am` (from qd), fem branches, fable special cases | `make check` in `mplapack/test/lin/<b>` | Rtest, Ctest, Rtest_rfp, Ctest_rfp: 0 failures, same pass lines as dd's committed results |
| eig | `eig/<b>/Makefile.am`, fable special cases | `make check` in `mplapack/test/eig/<b>`; hours | **not run** |
| matgen | `matgen/Makefile.am`, `Claror` | built by lin/eig | builds |
| link checks | `{lin,eig}/<b>/Makefile.am` in the current form | `sh misc/check_test_link_order.sh .`; `misc/check_test_override_link.sh <build> .`; `misc/check_backend_support_libraries.sh <build> . <qd libdir>` | pass |
| symbols | | `misc/check_unique_symbols.sh libmplapack_<b>.a libmplapack_<b>_opt.a` | pass |
| CMake | nothing (table) | `ctest` | library; mpblas 76/76; `compare_td_Rlamch`, `compare_td_Rgesv` pass |
| CMake consume / pkg-config | section 2.6 | `ctest -R 'consume\|pkgconfig'` | pass, when `CMAKE_INSTALL_PREFIX` is the consume prefix (below) |
| examples | section 2.7 | build with both systems; run `*_<b>` | 67/67 exit 0 |
| release | sections 2.7, 4 F6 | `python3 release/test-fable-comparison.py`, `python3 release/test-source-preparation.py`; `make dist` and compare `M`'s and `B`'s files in the tarball; td-only build from the tarball | pass; tarball differs from dd only by dd's CUDA files |

### 5.2 Generating the compare reference data

`Rlamch_reference.txt` and `Rlaruv_reference.txt` are the backend's own
`Rlamch.txt`/`Rlaruv.txt`.  Accept them only if `Rlamch.test_<b>` reports
"Testing Rlamch successful": it computes the expected values with MPFR from
`t`, `emin`, `emax` (td: E=2^-157, P=2^-156, N=157, M=-915, L=1024,
U=S=2^-916, Blue -458/434/536/-590).  For Rlaruv td checked: 108128 values
in [0,1), all distinct, mean 0.4991, variance 0.08315 (1/12 = 0.08333),
chi-square 92.8 over 100 bins.

### 5.3 Test pitfalls

| Pitfall | Effect | Fix |
|---|---|---|
| no `Mxerbla_<b>` | invalid-argument tests reach the library `Mxerbla` and exit (td: 98/150 failed) | add it (2.3) |
| `Rlamch.test.cpp` branch copied by hand | leftover `Rlamch_dd_test`, two-limb probes | use the rename tool, then build the lowest-limb probes with `L` limbs |
| lin/eig `Makefile.am` from an old copy | old link order after a merge, no conflict (new files) | copy the current file; `check_test_link_order.sh` |
| bundled libQD3 not installed | test and example binaries fail to load `libqd.so.2` | `LD_LIBRARY_PATH=<build>/external/i/QD/lib` |
| long examples | `03_SymmetricEigenproblems/Rsyev*_Frank`, `Rsyevd_DingDong` sweep n=5..1000 | `-STARTN 5 -ENDN 200 -STEPN 15` |
| CMake pkg-config tests | `.pc` files carry `CMAKE_INSTALL_PREFIX` from configure time; with the default `/usr/local` they fail for every backend | configure with `-DCMAKE_INSTALL_PREFIX=<build>/cmake/tests/consume/install` |
| CMake mpblas tests | `mplapack.test.cpp` was globbed as a test without `main()`, so no CMake test build ever completed (fixed on master, 2.4, 3.1) | — |
| CMake mpfr mpblas tests | after that fix, `Rnrm2`, `RCnrm2`, `Crotg` fail for mpfr under CMake but pass under autotools | open, not td-specific |

## 6. Other precisions

What changes relative to td, from the constants in section 3.

| | ds / ts / qs | edd | float |
|---|---|---|---|
| model | td (same libQD3 path; `single_real<N>`, one header) | dd (expansion code), binary80 (platform) | double |
| traits | `qd,nofma` | `qd,nofma` | none |
| complex | `ds/ts/qs_complex` (new in libQD3 1.6.0) | `edd_complex` | `std::complex<float>` |
| `t` | 46 / 70 / 94 | 126 | 24 |
| `emin` | -101 / -77 / -53 | -16317 | -125 |
| `emax` | 128 | 16384 | 128 |
| `_PRECISION` | 13 / 20 / 27 | 37 | from double's header, float values |
| Rlaruv limb scale | `2^(-24k)` | `2^(-64k)` | — |
| write digits | not `16·L`; at least `_ndigits+1` (15 / 22 / 29), check by reading back lin/eig input | at least 39 | as double |
| `compare_tolerance`, `EPSILON*` | not `1e-(15·L)`; measure | measure | from double's, scaled |
| platform | all | x86 only; the bundled libQD3 disables edd on MinGW (`external/qd/Makefile.am`); condition like binary80 | all |
| Blue exponents tsml/tbig/ssml/sbig | ds -51/41/74/-87; ts -39/29/74/-99; qs -27/17/74/-111 | -8159/8129/8222/-8255 | as double |
| risk | `emin` rises by 24 per limb.  For ts and qs, `sbig` (2^-99, 2^-111) is below the type's smallest normalized value (2^-78, 2^-54), so the Blue scaling in `Rnrm2`/`Rlassq` uses a constant the type cannot hold with all limbs normalized; decide how `arithmetic_params` handles it (clamp, or accept a leading-limb-only constant) before the tests, and expect lin/eig threshold failures that need the special cases of 2.4 revisited | long double formatting and `ldexp` on `_Float64x`; MinGW | small: no expansion; `Rlamch` and tolerances |
| order (section 7.4) | after td, together | after ds/ts/qs | any time |

For every new type: compute `t`, `emin`, `emax`, `_ndigits` from libQD3
(as in section 3, by printing them), derive the **derived** rows, choose the
**chosen** rows, and confirm with compare Rlamch (MPFR-derived) before
running lin/eig.

## 7. Record of td

### 7.1 Decisions

| Site | Decision |
|---|---|
| fable special cases guarding non-IEEE arithmetic or the double range | join dd/qd (2.4) |
| `Cchkbd` | qd's range |
| `Rget32`, `Rget34` (qd only), `Rchkee` threshold (dd only) | do not follow |
| defaults | CMake `OFF`, `--enable-td` default `no` |
| packaging, release, reconfig | enabled wherever dd is |
| benchmark | high-precision loop, plot 2 |

### 7.2 Failures during td, and what reports them now

| Failure | Found by | Now reported by |
|---|---|---|
| no `Mxerbla_td` | 98/150 mpblas failures | `backend_sites.py` (`NAME`) |
| no fem branch | `#error` in the matgen build | `PP-OWN` |
| `.h.in` missing from KEEP lists | review after the fact | `NAME` |
| lin/eig `Makefile.am` in the pre-merge form | merging master | `check_test_link_order.sh` |
| release checks with own lists (`check_test_override_link.sh`, `check_backend_support_install.sh`, `compare-fable-outputs.py`, `check_backend_support_libraries.sh`) | `backend_sites.py` | `LIST` |
| no CMake consume / pkg-config test | release review | `LIST`, `NAME` |
| `Rlamch.test.cpp` dd leftovers | compile error | use the rename tool |
| `LDATD`, `top_builtdir` | build | use the rename tool |

### 7.3 Expected residual of `backend_sites.py dd td`

| Kind | Site | Why |
|---|---|---|
| `PP-OWN`, `PP-GROUP` | `Cchkbd`, `Rchkee` (patch and source) | 7.1 |
| `ENABLE` | `misc/reconfig.macOS.clang.sanitize.sh`, `misc/reconfig.ubuntu24.04.clang.asan.sh`, `release/docker/common/check-fable-reproduction.sh` | they disable dd |
| `LIST` | `benchmark/go.*.sh.in` fast loop | td is in the other loop |
| `LIST` | fallback lists in `check_test_link_order.sh`, `check_test_override_link.sh`, `compare-fable-outputs.py` | used only without `backends.txt` |
| `LIST`, `NAME` | `examples/README` | historical |
| `NAME` | `cmake/tests/consume/main.cpp` | single-backend sample |
| `NAME` | `examples/mpblas/generic/Makefile.linux_cuda.in` | dd-only CUDA |
| `NAME` | comments in `include/mplapack_symbol_resolver.h`, `mplapack/reference/Rlamch.cpp` | comments |
| `PLOT` | `benchmark/*1.plt.in` | td in plot 2 |

### 7.4 Status and plan

| Backend | Status |
|---|---|
| td | done; eig suite and platforms other than Linux x86_64 GCC not run |
| ds, ts, qs | next, together, model td |
| edd | after ds/ts/qs |
| float | any time, model double |
