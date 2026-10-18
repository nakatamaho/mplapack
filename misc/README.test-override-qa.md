# Test override dependency QA

Ordinary LIN and EIG drivers compile Mxerbla, iMlaenv and Mxlaenv into the
executable so production definitions cannot replace the test overrides under
shared-library lookup or --as-needed. DMD drivers remain separate: they own
Mxerbla and Mxlaenv, but use the production iMlaenv (ISPEC=9 defaults to 25).
Both build systems build no separate LIN/EIG override libraries: the drivers are the
sole owners of those routines and common arrays. Linking another copy could
initialize and finalize the MPFR arrays twice.

CMake also registers test_override_lin/eig/dmd_<backend>[_opt] smoke tests.
Run them with ctest --test-dir <builddir> -R '^test_override_' after building.
They check error handling, tuning ownership and normal process finalization.

MinGW static LIN/EIG links retain each constructor-only backend initializer
with `-Wl,-u,<symbol>`, including reference, optimized and DMD drivers. CMake
static test executables do the same. GMP is different on 3.x: its normally
referenced context provider replaces the old GMP initializer, so do not add
an undefined `mplapack_initialize_gmp` reference. MPFR's initializer is a no-op;
DD/TD/QD initialization remains necessary for the x87 precision control word.

Check both link order and initializer retention with:

    sh misc/check_test_link_order.sh .

An isolated runtime regression compiles the actual backend initializer into
an archive. It presets x87 extended precision before the backend constructor
and checks that retaining only the initializer selects double precision;
its unretained negative control must fail:

    CPPFLAGS='-I<qd-prefix>/include' LIBS='-L<qd-prefix>/lib -lqd' \
      sh misc/check_backend_initializer.sh . dd

Repeat for td and qd. For MinGW set CXX/AR/NM to the cross tools,
LOG_COMPILER=wine, and WINEPATH to the dependency DLL directories. Supply an
rpath or temporary LD_LIBRARY_PATH for a native non-system QD installation.
The probe skips non-x86 hosts; it does not replace numerical QA.

After building the configured GMP, MPFR, QD or DD reference and optimized
support libraries, run:

    sh misc/check_test_override_link.sh <top-builddir> <top-srcdir>

The Linux probe checks error handling, test parameter exchange, and the DMD
production default without running release numerical suites. If actual MPFR
LIN drivers exist, it also runs their DGE error-exit checks with zero-size
matrices. MPFR uses the master MPFRXX_DEFAULT_EMIN/EMAX environment interface.

Run the installation regression with an isolated configured prefix:

    sh misc/check_backend_support_install.sh <top-builddir>

It stages already built primary, MATGEN, LIN and EIG libraries into a fresh
temporary DESTDIR using the real Automake install rules. No production
installation is changed. Support-library harness globals are supplied by
the final executable; the existing support relocation check distinguishes
these from unresolved external backend symbols.

Run lightweight release regressions independently of active QA:

    python3 release/test-fable-comparison.py
    python3 release/test-source-preparation.py

Full Fable conversion, native macOS/MinGW execution and numerical validation
remain part of release QA. These smoke checks do not replace that QA.
