# Test override dependency QA

Ordinary LIN and EIG drivers compile Mxerbla, iMlaenv and Mxlaenv into the
executable so production definitions cannot replace the test overrides under
shared-library lookup or --as-needed. DMD drivers remain separate: they own
Mxerbla and Mxlaenv, but use the production iMlaenv (ISPEC=9 defaults to 25).
Neither LIN nor EIG support libraries depend on the override libraries.

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
