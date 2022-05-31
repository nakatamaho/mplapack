mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Cgetri_test_mpfr Rgesv_test_mpfr Rgetri_Hilbert_mpfr Rgetri_test2_mpfr Rgetri_test3_mpfr Rgetri_test_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Cgetri_test_mpfr_SOURCES = Cgetri_test_mpfr.cpp
Cgetri_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgetri_test_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgesv_test_mpfr_SOURCES = Rgesv_test_mpfr.cpp
Rgesv_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgesv_test_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgetri_Hilbert_mpfr_SOURCES = Rgetri_Hilbert_mpfr.cpp
Rgetri_Hilbert_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgetri_Hilbert_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgetri_test2_mpfr_SOURCES = Rgetri_test2_mpfr.cpp
Rgetri_test2_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgetri_test2_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgetri_test3_mpfr_SOURCES = Rgetri_test3_mpfr.cpp
Rgetri_test3_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgetri_test3_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgetri_test_mpfr_SOURCES = Rgetri_test_mpfr.cpp
Rgetri_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgetri_test_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Cgetri_test_gmp Rgesv_test_gmp Rgetri_Hilbert_gmp Rgetri_test2_gmp Rgetri_test3_gmp Rgetri_test_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp -L$(GMP_LIBDIR) -lgmp

Cgetri_test_gmp_SOURCES = Cgetri_test_gmp.cpp
Cgetri_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Cgetri_test_gmp_LDFLAGS = $(gmp_libdepends)

Rgesv_test_gmp_SOURCES = Rgesv_test_gmp.cpp
Rgesv_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgesv_test_gmp_LDFLAGS = $(gmp_libdepends)

Rgetri_Hilbert_gmp_SOURCES = Rgetri_Hilbert_gmp.cpp
Rgetri_Hilbert_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgetri_Hilbert_gmp_LDFLAGS = $(gmp_libdepends)

Rgetri_test2_gmp_SOURCES = Rgetri_test2_gmp.cpp
Rgetri_test2_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgetri_test2_gmp_LDFLAGS = $(gmp_libdepends)

Rgetri_test3_gmp_SOURCES = Rgetri_test3_gmp.cpp
Rgetri_test3_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgetri_test3_gmp_LDFLAGS = $(gmp_libdepends)

Rgetri_test_gmp_SOURCES = Rgetri_test_gmp.cpp
Rgetri_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgetri_test_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE__FLOAT128
mplapackexamples_PROGRAMS += Cgetri_test__Float128 Rgesv_test__Float128 Rgetri_Hilbert__Float128 Rgetri_test2__Float128 Rgetri_test3__Float128 Rgetri_test__Float128

_Float128_cxxflags = $(OPENMP_CXXFLAGS)
if WANT_QUADMATH
_Float128_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128 -lquadmath
else
_Float128_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128
endif

Cgetri_test__Float128_SOURCES = Cgetri_test__Float128.cpp
Cgetri_test__Float128_CXXFLAGS = $(_Float128_cxxflags)
Cgetri_test__Float128_LDFLAGS = $(_Float128_libdepends)

Rgesv_test__Float128_SOURCES = Rgesv_test__Float128.cpp
Rgesv_test__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgesv_test__Float128_LDFLAGS = $(_Float128_libdepends)

Rgetri_Hilbert__Float128_SOURCES = Rgetri_Hilbert__Float128.cpp
Rgetri_Hilbert__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgetri_Hilbert__Float128_LDFLAGS = $(_Float128_libdepends)

Rgetri_test2__Float128_SOURCES = Rgetri_test2__Float128.cpp
Rgetri_test2__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgetri_test2__Float128_LDFLAGS = $(_Float128_libdepends)

Rgetri_test3__Float128_SOURCES = Rgetri_test3__Float128.cpp
Rgetri_test3__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgetri_test3__Float128_LDFLAGS = $(_Float128_libdepends)

Rgetri_test__Float128_SOURCES = Rgetri_test__Float128.cpp
Rgetri_test__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgetri_test__Float128_LDFLAGS = $(_Float128_libdepends)

endif

if ENABLE__FLOAT64X
mplapackexamples_PROGRAMS += Cgetri_test__Float64x Rgesv_test__Float64x Rgetri_Hilbert__Float64x Rgetri_test2__Float64x Rgetri_test3__Float64x Rgetri_test__Float64x

_Float64x_cxxflags = $(OPENMP_CXXFLAGS)
_Float64x_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack__Float64x -L$(top_builddir)/mpblas/reference -lmpblas__Float64x

Cgetri_test__Float64x_SOURCES = Cgetri_test__Float64x.cpp
Cgetri_test__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Cgetri_test__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgesv_test__Float64x_SOURCES = Rgesv_test__Float64x.cpp
Rgesv_test__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgesv_test__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgetri_Hilbert__Float64x_SOURCES = Rgetri_Hilbert__Float64x.cpp
Rgetri_Hilbert__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgetri_Hilbert__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgetri_test2__Float64x_SOURCES = Rgetri_test2__Float64x.cpp
Rgetri_test2__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgetri_test2__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgetri_test3__Float64x_SOURCES = Rgetri_test3__Float64x.cpp
Rgetri_test3__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgetri_test3__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgetri_test__Float64x_SOURCES = Rgetri_test__Float64x.cpp
Rgetri_test__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgetri_test__Float64x_LDFLAGS = $(_Float64x_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Cgetri_test_double Rgesv_test_double Rgetri_Hilbert_double Rgetri_test2_double Rgetri_test3_double Rgetri_test_double

double_cxxflags = $(OPENMP_CXXFLAGS)
double_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_double -L$(top_builddir)/mpblas/reference -lmpblas_double

Cgetri_test_double_SOURCES = Cgetri_test_double.cpp
Cgetri_test_double_CXXFLAGS = $(double_cxxflags)
Cgetri_test_double_LDFLAGS = $(double_libdepends)

Rgesv_test_double_SOURCES = Rgesv_test_double.cpp
Rgesv_test_double_CXXFLAGS = $(double_cxxflags)
Rgesv_test_double_LDFLAGS = $(double_libdepends)

Rgetri_Hilbert_double_SOURCES = Rgetri_Hilbert_double.cpp
Rgetri_Hilbert_double_CXXFLAGS = $(double_cxxflags)
Rgetri_Hilbert_double_LDFLAGS = $(double_libdepends)

Rgetri_test2_double_SOURCES = Rgetri_test2_double.cpp
Rgetri_test2_double_CXXFLAGS = $(double_cxxflags)
Rgetri_test2_double_LDFLAGS = $(double_libdepends)

Rgetri_test3_double_SOURCES = Rgetri_test3_double.cpp
Rgetri_test3_double_CXXFLAGS = $(double_cxxflags)
Rgetri_test3_double_LDFLAGS = $(double_libdepends)

Rgetri_test_double_SOURCES = Rgetri_test_double.cpp
Rgetri_test_double_CXXFLAGS = $(double_cxxflags)
Rgetri_test_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Cgetri_test_dd Rgesv_test_dd Rgetri_Hilbert_dd Rgetri_test2_dd Rgetri_test3_dd Rgetri_test_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(top_builddir)/mpblas/reference -lmpblas_dd -L$(QD_LIBDIR) -lqd

Cgetri_test_dd_SOURCES = Cgetri_test_dd.cpp
Cgetri_test_dd_CXXFLAGS = $(dd_cxxflags)
Cgetri_test_dd_LDFLAGS = $(dd_libdepends)

Rgesv_test_dd_SOURCES = Rgesv_test_dd.cpp
Rgesv_test_dd_CXXFLAGS = $(dd_cxxflags)
Rgesv_test_dd_LDFLAGS = $(dd_libdepends)

Rgetri_Hilbert_dd_SOURCES = Rgetri_Hilbert_dd.cpp
Rgetri_Hilbert_dd_CXXFLAGS = $(dd_cxxflags)
Rgetri_Hilbert_dd_LDFLAGS = $(dd_libdepends)

Rgetri_test2_dd_SOURCES = Rgetri_test2_dd.cpp
Rgetri_test2_dd_CXXFLAGS = $(dd_cxxflags)
Rgetri_test2_dd_LDFLAGS = $(dd_libdepends)

Rgetri_test3_dd_SOURCES = Rgetri_test3_dd.cpp
Rgetri_test3_dd_CXXFLAGS = $(dd_cxxflags)
Rgetri_test3_dd_LDFLAGS = $(dd_libdepends)

Rgetri_test_dd_SOURCES = Rgetri_test_dd.cpp
Rgetri_test_dd_CXXFLAGS = $(dd_cxxflags)
Rgetri_test_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Cgetri_test_qd Rgesv_test_qd Rgetri_Hilbert_qd Rgetri_test2_qd Rgetri_test3_qd Rgetri_test_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(top_builddir)/mpblas/reference -lmpblas_qd -L$(QD_LIBDIR) -lqd

Cgetri_test_qd_SOURCES = Cgetri_test_qd.cpp
Cgetri_test_qd_CXXFLAGS = $(qd_cxxflags)
Cgetri_test_qd_LDFLAGS = $(qd_libdepends)

Rgesv_test_qd_SOURCES = Rgesv_test_qd.cpp
Rgesv_test_qd_CXXFLAGS = $(qd_cxxflags)
Rgesv_test_qd_LDFLAGS = $(qd_libdepends)

Rgetri_Hilbert_qd_SOURCES = Rgetri_Hilbert_qd.cpp
Rgetri_Hilbert_qd_CXXFLAGS = $(qd_cxxflags)
Rgetri_Hilbert_qd_LDFLAGS = $(qd_libdepends)

Rgetri_test2_qd_SOURCES = Rgetri_test2_qd.cpp
Rgetri_test2_qd_CXXFLAGS = $(qd_cxxflags)
Rgetri_test2_qd_LDFLAGS = $(qd_libdepends)

Rgetri_test3_qd_SOURCES = Rgetri_test3_qd.cpp
Rgetri_test3_qd_CXXFLAGS = $(qd_cxxflags)
Rgetri_test3_qd_LDFLAGS = $(qd_libdepends)

Rgetri_test_qd_SOURCES = Rgetri_test_qd.cpp
Rgetri_test_qd_CXXFLAGS = $(qd_cxxflags)
Rgetri_test_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/mplapack/examples/mplapack/00_LinearEquations
