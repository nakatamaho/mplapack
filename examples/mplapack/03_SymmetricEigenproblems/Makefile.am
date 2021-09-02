mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Rsyev_test_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Rsyev_test_mpfr_SOURCES = Rsyev_test_mpfr.cpp
Rsyev_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rsyev_test_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Rsyev_test_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp -L$(GMP_LIBDIR) -lgmp

Rsyev_test_gmp_SOURCES = Rsyev_test_gmp.cpp
Rsyev_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Rsyev_test_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE__FLOAT128
mplapackexamples_PROGRAMS += Rsyev_test__Float128

_Float128_cxxflags = $(OPENMP_CXXFLAGS)
if WANT_QUADMATH
_Float128_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128 -lquadmath
else
_Float128_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128
endif

Rsyev_test__Float128_SOURCES = Rsyev_test__Float128.cpp
Rsyev_test__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rsyev_test__Float128_LDFLAGS = $(_Float128_libdepends)

endif

if ENABLE__FLOAT64X
mplapackexamples_PROGRAMS += Rsyev_test__Float64x

_Float64x_cxxflags = $(OPENMP_CXXFLAGS)
_Float64x_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack__Float64x -L$(top_builddir)/mpblas/reference -lmpblas__Float64x

Rsyev_test__Float64x_SOURCES = Rsyev_test__Float64x.cpp
Rsyev_test__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rsyev_test__Float64x_LDFLAGS = $(_Float64x_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Rsyev_test_double

double_cxxflags = $(OPENMP_CXXFLAGS)
double_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_double -L$(top_builddir)/mpblas/reference -lmpblas_double

Rsyev_test_double_SOURCES = Rsyev_test_double.cpp
Rsyev_test_double_CXXFLAGS = $(double_cxxflags)
Rsyev_test_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Rsyev_test_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(top_builddir)/mpblas/reference -lmpblas_dd -L$(QD_LIBDIR) -lqd

Rsyev_test_dd_SOURCES = Rsyev_test_dd.cpp
Rsyev_test_dd_CXXFLAGS = $(dd_cxxflags)
Rsyev_test_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Rsyev_test_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(top_builddir)/mpblas/reference -lmpblas_qd -L$(QD_LIBDIR) -lqd

Rsyev_test_qd_SOURCES = Rsyev_test_qd.cpp
Rsyev_test_qd_CXXFLAGS = $(qd_cxxflags)
Rsyev_test_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/mplapack/examples/03_SymmetricEigenproblems
