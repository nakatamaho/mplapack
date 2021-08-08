if ENABLE_MPFR
mplapackexamples_PROGRAMS = Rgesvd_test_mpfr Rgesvd_readfromfile_mpfr

mpfrcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends =  -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Rgesvd_test_mpfr_SOURCES= Rgesvd_test_mpfr.cpp
Rgesvd_test_mpfr_CXXFLAGS= $(mpfrcxxflags)
Rgesvd_test_mpfr_LDFLAGS= $(mpfr_libdepends)

Rgesvd_readfromfile_mpfr_SOURCES= Rgesvd_readfromfile_mpfr.cpp
Rgesvd_readfromfile_mpfr_CXXFLAGS= $(mpfrcxxflags)
Rgesvd_readfromfile_mpfr_LDFLAGS= $(mpfr_libdepends)

endif


mplapackexamplesdir=$(prefix)/share/mplapack/examples/05_SingularValueDecomposition

