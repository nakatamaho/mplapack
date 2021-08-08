if ENABLE_MPFR
mplapackexamples_PROGRAMS = Rgges_test_mpfr
mpfrcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends =  -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Rgges_test_mpfr_SOURCES= Rgges_test_mpfr.cpp
Rgges_test_mpfr_CXXFLAGS= $(mpfrcxxflags)
Rgges_test_mpfr_LDFLAGS= $(mpfr_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/mplapack/examples/04_NonsymmetricEigenproblems
