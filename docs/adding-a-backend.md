# Adding an arithmetic backend

This is the checklist for adding a backend (a REAL/COMPLEX type pair such as
`td_real`/`td_complex`) to MPLAPACK.  It reflects the tree as of the 3.1
branch.  Parts of the work are generated from `backends.txt`; the rest is
still by hand and is listed here so nothing is missed.

Read `AGENTS.md` first.  In particular: sources under `mpblas/reference/`
and `mplapack/reference/` come from the fable pipeline and are never edited
by hand, and both build systems (autotools and CMake) must keep working.

## 1. What is generated and what is not

| Area | Source of truth | Status |
|---|---|---|
| CMake options, library targets, dependency discovery, tests, examples, benchmarks | `backends.txt` | generated |
| Prototype headers `include/{mpblas,mplapack,mplapack_lin,mplapack_eig,mplapack_matgen}_<b>.h` | `backends.txt` + `.h.in` templates, via `fable/gen_include_*.sh` | generated |
| `mpblas/optimized/<b>/Makefile.am`, `mpblas/test/<b>/Makefile.am`, `mplapack/test/compare/<b>/Makefile.am`, `benchmark/Makefile.<b>.am` | `backends.txt` + `misc/backend_makefiles/`, via `misc/gen_backend_makefiles.py` | generated |
| Type-specific headers, optimized BLAS sources, reference data | — | by hand (section 3) |
| Shared sources with one branch per backend | — | by hand (section 4) |
| `configure.ac`, the other `Makefile.am` files | — | by hand (section 6) |
| Example sources and their autotools makefiles | `examples/*/generic/` templates, via `examples/gen_all.sh` | generated; the templates by hand (section 7) |
| Benchmark driver and plot scripts, packaging, release scripts, README | — | by hand (section 7) |

## 2. Table entries

1. Add a row to `backends.txt`: name, REAL type, COMPLEX type, default of
   the CMake option `MPLAPACK_ENABLE_<NAME>`, and build traits.  The file
   documents the columns and the traits (`gmp`, `mpfr`, `qd`, `nofma`,
   `binary128libs`).  A libQD3 expansion type (dd, qd, td, ds, ts, qs, edd)
   takes `qd,nofma`: it links libQD3 and must be compiled with
   `-ffp-contract=off`.
2. If the backend needs autotools settings that the traits do not cover,
   add a section to `misc/backend_makefiles/backends.conf`.  Every backend
   whose compare tests are not bit-exact needs `compare_tolerance`.
3. Regenerate:

   ```sh
   python3 misc/gen_backend_makefiles.py
   bash fable/gen_include_mpblas.sh
   bash fable/gen_include_mplapack.sh
   bash fable/gen_include_mplapack_lin.sh
   bash fable/gen_include_mplapack_eig.sh
   bash fable/gen_include_mplapack_matgen.sh
   ```

   The fable scripts need universal-ctags and clang-format 19 or newer
   (`CLANG_FORMAT=clang-format-19`).  They rewrite the headers of every
   backend; review the diff, since a different ctags version can add or drop
   prototypes of existing backends too.

## 3. Per-backend files written by hand

Start from the closest existing backend (dd or qd for a libQD3 type, double
for `float`) and adapt.  `<b>` is the backend name.

Headers in `include/`:

- `mplapack_arithmetic_params_<b>.h`: digits, emin, emax and the Blue
  scaling exponents used by `Rlamch`, `Rnrm2`, `Rlassq`.  Derive them from
  the type's own constants (for libQD3: `_eps`, `_min_normalized`, `_max`),
  not from the underlying `float`/`double`: the lower limbs of an expansion
  must stay normalized, which raises the effective emin (dd has emin -968,
  not -1021).
- `mplapack_utils_<b>.h`: conversions, printing and helper functions.
- `mplapack_benchmark_<b>.h`: mangled routine names for the benchmarks.
  Generate it after the library builds with
  `misc/gen_benchmark_header.py <b> path/to/libmplapack_<b>.a`; renaming
  another backend's header gives wrong names, because libQD3 complex types
  are aliases of `qd3_complex<T>`.

Templates for the generated prototype headers (section 2, step 3):

- `mpblas/reference/mpblas_<b>.h.in`
- `mplapack/reference/mplapack_<b>.h.in`
- `mplapack/test/lin/common/mplapack_lin_<b>.h.in`
- `mplapack/test/eig/common/mplapack_eig_<b>.h.in`
- `mplapack/test/matgen/mplapack_matgen_<b>.h.in`

Optimized BLAS in `mpblas/optimized/<b>/`: `Raxpy.cpp`, `Rcopy.cpp`,
`Rdot.cpp`, `Rgemm.cpp`, their `*_ref.cpp`, and `openmp/` (`Raxpy_omp.cpp`,
`Rcopy_omp.cpp`, `Rdot_omp.cpp`, `Rgemm_omp.cpp`, `Rgemm_{NN,NT,TN,TT}_omp.cpp`).
The generated `Makefile.am` lists exactly these files.

`mplapack/optimized/<b>/Makefile.am` (two lines) and `.gitkeep`; copy from
another backend.

Test directories with their own `Makefile.am` (not generated yet):
`mplapack/test/lin/<b>/` and `mplapack/test/eig/<b>/`.  These differ between
backends by more than the name; compare the existing ones before copying.

Reference data for the compare tests,
`mplapack/test/compare/<b>/Rlamch_reference.txt` and `Rlaruv_reference.txt`:
these are the `Rlamch.txt` and `Rlaruv.txt` that the backend's own
`Rlamch.test_<b>` and `Rlaruv.test_<b>` write.  Generate them once the
backend builds and commit them after review.  `Rlamch.test.cpp` needs a
branch for the new type (section 4): it derives the expected Rlamch values
with MPFR from the type's digits, emin and emax and reports
"Testing Rlamch successful" only when the computed values match.  The comparison uses `compare_tolerance` from
`backends.conf`, or an exact diff when none is set.

## 4. Shared sources with one branch per backend

Each of these has an `#if defined MPLAPACK_BUILD_WITH_<NAME>` block for
every backend.  A new backend needs its own block in all of them.

Public and test headers in `include/`:

- `mpblas.h`, `mplapack.h`: REAL/COMPLEX typedefs and the `_<b>` name
  mappings (`Rlamch_<b>`, `iMlaenv_<b>`, ...).
- `mplapack_arithmetic_params.h`, `mplapack_utils.h`,
  `mplapack_benchmark.h`, `mplapack_compare_debug.h`,
  `mplapack_lin.h`, `mplapack_eig.h`, `mplapack_matgen.h`.

Library sources that are MPLAPACK's own (not fable output):

- `mpblas/reference/mplapackinit.cpp`
- `mplapack/reference/Rlamch.cpp`, `Rlaruv.cpp`, `Risnan.cpp`,
  `Risinf.cpp`, `Mexponent.cpp`

Test sources:

- `mpblas/test/common/arithmetic.test.cpp`, `complex.test.cpp`,
  `mplapack.test.cpp`
- `mplapack/test/compare/common/Rlamch.test.cpp`

Fortran I/O emulation used by the lin/eig/matgen test code:

- `fable/fem/read.hpp`, `fable/fem/write.hpp`: backend utility header,
  type headers, and `read_loop`/`write_loop` overloads for REAL and COMPLEX.
  These end in `#error` for an unknown backend, so a library-only build does
  not reveal a missing branch; only the test code includes them.

To find what still lacks a branch, list the files that mention a similar
backend but not the new one (here for td, modelled on dd):

```sh
git ls-files | grep -v '^external/\|/results/\|^examples/' \
  | xargs grep -l MPLAPACK_BUILD_WITH_DD | xargs grep -L MPLAPACK_BUILD_WITH_TD
```

`fable/3.9.1/` holds patches for an older LAPACK and is not used.

### Places that were missed for td

All three were found only by failing tests, not by review.  The first and
third are not found by the grep above.

- `mpblas/test/common/Mxerbla.override.cpp` defines `Mxerbla_<b>` for every
  backend without any `MPLAPACK_BUILD_WITH_*` guard.  Without a definition
  for the new backend, the invalid-argument tests reach the library's own
  `Mxerbla`, which exits; for td 98 of 150 mpblas tests failed.
- `fable/fem/read.hpp`, `fable/fem/write.hpp` (listed above): without a
  branch the lin/eig/matgen test build stops at `#error`, while the library
  and the mpblas tests build normally.
- `fable/convert_blas_all.sh`, `fable/convert_lapack_all.sh`,
  `fable/go_testing.sh`: the `KEEP_HAND_WRITTEN_FILES` lists must include the
  new `*_<b>.h.in` files.  Otherwise a fable re-conversion deletes them; the
  current build is unaffected, so nothing fails until the next conversion.

Found while doing section 7 for td:

- `mplapack/test/{lin,eig}/td/Makefile.am` were copied before master's
  link-order fix and kept the old form after the merge; the merge itself
  had no conflict, because the td files are new on 3.1.
  `misc/check_test_link_order.sh` did not look at td, since its backend list
  was written by hand.
- The release QA checks `misc/check_test_override_link.sh`,
  `misc/check_backend_support_install.sh` and
  `release/docker/common/compare-fable-outputs.py` also had their own
  backend lists without td.

## 5. Special cases to decide per backend

Some fable-generated sources treat a subset of backends differently.  For
a new backend, decide for each one whether it joins the subset.  The
subsets at the time of writing:

| Subset | Files (reference / test) |
|---|---|
| binary80, binary128 | `Cgeev`, `Cgges`, `Cgges3`, `Clatrs`; tests `Cget23`, `Cqrt13`, `Cqrt15`, `Clatb4` |
| dd, td, qd | `Cgejsv`, `Cgesvj`, `Rgejsv`, `Rgesvj` (with gmp, mpfr), `iMieeeck` (with gmp); tests `Cchkbd`, `Rchkbd`, `Clatb4`, `Rlatb4` (dd, td); `Rget32`, `Rget34` (qd only); matgen `Claror` (dd, td) |
| gmp, mpfr | `Rlarrb`, `Rlarrd`, `Rlarrk`, `Rstebz`, `Cgges3`; tests `Rdrgev3`, `Rget38` |
| gmp only | `Cbdsqr`, `Rbdsqr`, `Rbdsvdx`; test `Cchkee` |
| binary80, double, gmp | tests `Cchktsqr`, `Rchktsqr` |
| binary80, dd, double, gmp, mpfr | test `Rchkee` |

Change these through the patches in `fable/3.12.1/lapack/patch-*.cpp` and
regenerate; never edit `mpblas/reference/`, `mplapack/reference/` or the
generated test sources directly.  When only `+` lines of a patch change
(adding a backend to an `#if` condition), editing the patch and the
generated source the same way is equivalent to regenerating; check it by
reverse-applying the patch to the edited source and comparing with the
reverse of the old pair (`patch -R -o out source < patch`).
`patch-Claror.cpp` contains its hunks twice and already reverse-applies with
offsets and a leftover line; that predates the td backend.

Decisions taken for td, as a reference for the next expansion types: td
joins every dd/qd case above that guards against non-IEEE arithmetic or
the double exponent range; in `Cchkbd` it takes qd's narrower singular
value range (`-(half * half) * log(ulp)`); it does not take the qd-only
`Rget32`/`Rget34` workarounds or the dd-only threshold increase in
`Rchkee`.  To list the current subsets:

```sh
grep -lr MPLAPACK_BUILD_WITH_ mplapack/reference mplapack/test/*/common \
    mplapack/test/matgen fable/3.12.1/lapack
```

`mpblas/test/common/*.test.cpp` and `mplapack/test/compare/common/*.test.cpp`
mention `MPLAPACK_BUILD_WITH_MPFR` because MPFR is the test oracle; those
need no change.

## 6. Wiring that is still by hand

`configure.ac`:

- `AC_ARG_ENABLE(<b>, ...)` and `AM_CONDITIONAL(ENABLE_<B>, ...)`, next to
  the existing backends.
- For a libQD3 backend: the conditions that decide whether QD is needed
  (`--with-system-qd` block and `need_internal_qdlib`) test
  `enable_qd`/`enable_dd`; add the new backend there.
- The "Enable/Disable ... version" and summary messages.
- Type detection, if the type depends on the compiler or platform (as
  binary80 and binary128 do).
- `AC_CONFIG_FILES`: `mpblas/optimized/<b>/Makefile`,
  `mplapack/optimized/<b>/Makefile`, `mpblas/test/<b>/Makefile`,
  `mplapack/test/compare/<b>/Makefile`, `mplapack/test/lin/<b>/Makefile`,
  `mplapack/test/eig/<b>/Makefile`.

`Makefile.am` files:

- `Makefile.am`: unique-symbol archive list, installed headers
  (`include/*_<b>.h`), pkg-config files `mplapack_<b>.pc`,
  `mplapack_<b>_opt.pc`.
- `mplapack/reference/Makefile.am`: the `libmplapack_<b>.la` library.
- `mpblas/optimized/Makefile.am`, `mpblas/test/Makefile.am`,
  `mplapack/test/compare/Makefile.am`: `SUBDIRS += <b>` under
  `if ENABLE_<B>`.
- `mplapack/optimized/Makefile.am`: `SUBDIRS`.
- `mplapack/test/lin/Makefile.am`, `mplapack/test/eig/Makefile.am`,
  `mplapack/test/matgen/Makefile.am`: the `lin_<b>`, `eig_<b>`,
  `matgen_<b>` libraries, `SUBDIRS`, `CHECK_BACKENDS`.
- `benchmark/Makefile.am`: `include Makefile.<b>.am` under `if ENABLE_<B>`.

When creating `mplapack/test/{lin,eig}/<b>/Makefile.am` by copying another
backend's, copy the current file: master changed their link order (libraries
in `LDADD`, after the executable's objects; no `--whole-archive` on MinGW).
`misc/check_test_link_order.sh .` fails on a file in the old form.

Scripts that list backends by hand: `cmake/tests/CMakeLists.txt`
(pkg-config consumer tests), and the second loop of
`misc/check_test_override_link.sh` (backends whose pkg-config file pulls in
an external library).  `misc/check_source_manifests.sh`,
`misc/check_test_link_order.sh`, the first loop of
`misc/check_test_override_link.sh` and
`release/docker/common/compare-fable-outputs.py` read `backends.txt`;
`misc/check_backend_support_install.sh` takes whatever
`mpblas/optimized/*/` was built.

## 7. Examples, benchmark scripts, packaging and README

Do these after the library and tests pass.  None of them reads
`backends.txt`.

Examples.  The sources under `examples/` are generated from
`examples/{mpblas,mplapack}/generic/`; CMake builds every
`examples/**/*_<b>.cpp` it finds, autotools needs the generated
`Makefile.am`.

- `generic/generate.sh` (mpblas and mplapack): the `MPLIBS` lists (the
  mplapack one has three, plus the `Cgeev_NPR` list without gmp and the
  `for _mplib in ...` loop with its `case` that maps a backend to
  `<B>LIBS`), the `REAL`/`COMPLEX` substitution branch, the
  `if ENABLE_<B>` block written to `Makefile.am`, and
  `DISABLED_EXAMPLE_SUFFIXES`.
- `generic/header_<b>` and `generic/header_<b>_complex`: copy from a
  similar backend.
- `generic/Makefile.{freebsd,linux,linux.inteloneAPI,macos,mingw}.in`:
  the `<B>LIBS` (and for mpblas `<B>OPTLIBS`) line, the `programs=` lists
  and one rule per program.  `Makefile.linux_cuda.in` is dd only.
- `examples/mplapack/run_smoke.sh`: `backends=`.
- Regenerate with `cd examples && bash gen_all.sh`.  Run it once before the
  change: it reproduces the committed files exactly, so afterwards
  `git diff examples` must show only the new backend.

Benchmark scripts (`benchmark/Makefile.<b>.am` is generated, section 1):

- `go.<Routine>.sh.in`: add the backend to one of the two `for _mplib`
  loops.
- `<Routine>1.plt.in` holds the fast types (binary80, binary128, dd),
  `<Routine>2.plt.in` the high-precision ones (MPFR, GMP, qd); add a plain
  and an `_opt` curve to one of them.  td went into the second, before qd.

Packaging and release:

- `--enable-<b>` in `packaging/PKGBUILD.in`, `packaging/debian/rules`,
  `packaging/mplapack.spec.in` (and its `Provides:` lines),
  `misc/reconfig.*.sh` and the configure options in `release/`
  (`buildtest_tier1_macos_*.sh`, `docker/common/tarball-smoke.sh`,
  `docker/distcheck/*.sh`, `docker/matrix/*`).  Configurations that turn dd
  off (sanitizer builds, `check-fable-reproduction.sh`) need nothing for a
  backend that is off by default.

Documentation:

- `README.md`: "Supported Precision Backends", and the CMake option if the
  backend is off by default.  The install snippets that download an older
  release tarball stay as they are.
- `MIGRATION.md` lists breaking changes only; a new backend adds nothing
  there.
- `doc/manual/manual.tex` still describes 2.0.1 and is committed with its
  PDF; it was not updated for td.

## 8. Verification

1. `python3 misc/gen_backend_makefiles.py --check` passes, and the fable
   header generators leave the other backends' headers unchanged.
2. CMake configures and builds with the new backend enabled and with it
   disabled.
3. autotools: `./gen_configure.sh`, configure with and without
   `--enable-<b>`, build.
4. `misc/check_unique_symbols.sh` on `libmplapack_<b>.a` and
   `libmplapack_<b>_opt.a`.
5. `make check-all` in `mpblas/test/<b>` (every BLAS test; plain
   `make check` runs only the smoke tests).
6. `make check` in `mplapack/test/compare/<b>` runs only the Rlamch and
   Rlaruv reference comparisons.  `make check-all` there currently runs
   none of the routine tests (no log rule matches `*.test_<b>`, for every
   backend), so run the programs directly:

   ```sh
   cd mplapack/test/compare/<b> && make check-all   # builds the programs
   for t in *.test_<b>; do ./$t > $t.out 2>&1 || echo "FAIL $t"; done
   ```

   Rpotri fails for dd and Classq fails for dd and qd as well; compare a new
   backend's failures with those before treating them as its own.
7. The lin and eig suites (`mplapack/test/lin`, `mplapack/test/eig`), and
   compare the pass lines with the committed results of a similar backend
   under `mplapack/test/*/results/`.  eig takes hours; run it on a full
   machine and keep the logs.

## 9. Plan for the 3.1 backends

The order chosen for float, ds, ts, qs, td and edd (libQD3 1.6.0 or later
provides all of them; the complex types for ds/ts/qs first appear in 1.6.0):

1. **td** first, end to end.  It goes through the same libQD3 path as
   ds/ts/qs/edd (`qd3_complex`, `ldexp`, `numeric_limits`) and is
   numerically between dd and qd, so it is the lowest-risk way to find what
   is still hard-coded.  The lin/eig `Makefile.am` files and the remaining
   `configure.ac`/`Makefile.am` lists are generalized as td needs them.
2. **ds, ts, qs** together: all three are `single_real<N>` from the same
   header.  This is the hardest step numerically: the exponent range is that
   of `float` and the effective emin rises with the number of limbs, so
   `arithmetic_params`, `Rlamch` and the test tolerances must come from
   measured values.  Their use is mainly as the CPU reference for GPU
   kernels.
3. **edd** (`_Float64x`/long double pairs): x86 only, conditioned like
   binary80.  The bundled libQD3 build disables edd on MinGW
   (`external/qd/Makefile.am`).
4. **float**: close to double; can go in at any point.
