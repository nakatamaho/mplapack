mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Cheev_test_mpfr Rsyev_Frank_mpfr Rsyev_test_mpfr Rsyevd_DingDong_mpfr Rsyevd_Frank_mpfr Rsyevr_Frank_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Cheev_test_mpfr_SOURCES = Cheev_test_mpfr.cpp
Cheev_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cheev_test_mpfr_LDFLAGS = $(mpfr_libdepends)

Rsyev_Frank_mpfr_SOURCES = Rsyev_Frank_mpfr.cpp
Rsyev_Frank_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rsyev_Frank_mpfr_LDFLAGS = $(mpfr_libdepends)

Rsyev_test_mpfr_SOURCES = Rsyev_test_mpfr.cpp
Rsyev_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rsyev_test_mpfr_LDFLAGS = $(mpfr_libdepends)

Rsyevd_DingDong_mpfr_SOURCES = Rsyevd_DingDong_mpfr.cpp
Rsyevd_DingDong_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rsyevd_DingDong_mpfr_LDFLAGS = $(mpfr_libdepends)

Rsyevd_Frank_mpfr_SOURCES = Rsyevd_Frank_mpfr.cpp
Rsyevd_Frank_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rsyevd_Frank_mpfr_LDFLAGS = $(mpfr_libdepends)

Rsyevr_Frank_mpfr_SOURCES = Rsyevr_Frank_mpfr.cpp
Rsyevr_Frank_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rsyevr_Frank_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Cheev_test_gmp Rsyev_Frank_gmp Rsyev_test_gmp Rsyevd_DingDong_gmp Rsyevd_Frank_gmp Rsyevr_Frank_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp -L$(GMP_LIBDIR) -lgmp

Cheev_test_gmp_SOURCES = Cheev_test_gmp.cpp
Cheev_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Cheev_test_gmp_LDFLAGS = $(gmp_libdepends)

Rsyev_Frank_gmp_SOURCES = Rsyev_Frank_gmp.cpp
Rsyev_Frank_gmp_CXXFLAGS = $(gmp_cxxflags)
Rsyev_Frank_gmp_LDFLAGS = $(gmp_libdepends)

Rsyev_test_gmp_SOURCES = Rsyev_test_gmp.cpp
Rsyev_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Rsyev_test_gmp_LDFLAGS = $(gmp_libdepends)

Rsyevd_DingDong_gmp_SOURCES = Rsyevd_DingDong_gmp.cpp
Rsyevd_DingDong_gmp_CXXFLAGS = $(gmp_cxxflags)
Rsyevd_DingDong_gmp_LDFLAGS = $(gmp_libdepends)

Rsyevd_Frank_gmp_SOURCES = Rsyevd_Frank_gmp.cpp
Rsyevd_Frank_gmp_CXXFLAGS = $(gmp_cxxflags)
Rsyevd_Frank_gmp_LDFLAGS = $(gmp_libdepends)

Rsyevr_Frank_gmp_SOURCES = Rsyevr_Frank_gmp.cpp
Rsyevr_Frank_gmp_CXXFLAGS = $(gmp_cxxflags)
Rsyevr_Frank_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE__FLOAT128
mplapackexamples_PROGRAMS += Cheev_test__Float128 Rsyev_Frank__Float128 Rsyev_test__Float128 Rsyevd_DingDong__Float128 Rsyevd_Frank__Float128 Rsyevr_Frank__Float128

_Float128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if WANT_QUADMATH
_Float128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128 -lquadmath
else
_Float128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128
endif

Cheev_test__Float128_SOURCES = Cheev_test__Float128.cpp
Cheev_test__Float128_CXXFLAGS = $(_Float128_cxxflags)
Cheev_test__Float128_LDFLAGS = $(_Float128_libdepends)

Rsyev_Frank__Float128_SOURCES = Rsyev_Frank__Float128.cpp
Rsyev_Frank__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rsyev_Frank__Float128_LDFLAGS = $(_Float128_libdepends)

Rsyev_test__Float128_SOURCES = Rsyev_test__Float128.cpp
Rsyev_test__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rsyev_test__Float128_LDFLAGS = $(_Float128_libdepends)

Rsyevd_DingDong__Float128_SOURCES = Rsyevd_DingDong__Float128.cpp
Rsyevd_DingDong__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rsyevd_DingDong__Float128_LDFLAGS = $(_Float128_libdepends)

Rsyevd_Frank__Float128_SOURCES = Rsyevd_Frank__Float128.cpp
Rsyevd_Frank__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rsyevd_Frank__Float128_LDFLAGS = $(_Float128_libdepends)

Rsyevr_Frank__Float128_SOURCES = Rsyevr_Frank__Float128.cpp
Rsyevr_Frank__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rsyevr_Frank__Float128_LDFLAGS = $(_Float128_libdepends)

endif

if ENABLE__FLOAT64X
mplapackexamples_PROGRAMS += Cheev_test__Float64x Rsyev_Frank__Float64x Rsyev_test__Float64x Rsyevd_DingDong__Float64x Rsyevd_Frank__Float64x Rsyevr_Frank__Float64x

_Float64x_cxxflags = $(OPENMP_CXXFLAGS)
_Float64x_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack__Float64x -L$(top_builddir)/mpblas/reference -lmpblas__Float64x

Cheev_test__Float64x_SOURCES = Cheev_test__Float64x.cpp
Cheev_test__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Cheev_test__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rsyev_Frank__Float64x_SOURCES = Rsyev_Frank__Float64x.cpp
Rsyev_Frank__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rsyev_Frank__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rsyev_test__Float64x_SOURCES = Rsyev_test__Float64x.cpp
Rsyev_test__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rsyev_test__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rsyevd_DingDong__Float64x_SOURCES = Rsyevd_DingDong__Float64x.cpp
Rsyevd_DingDong__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rsyevd_DingDong__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rsyevd_Frank__Float64x_SOURCES = Rsyevd_Frank__Float64x.cpp
Rsyevd_Frank__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rsyevd_Frank__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rsyevr_Frank__Float64x_SOURCES = Rsyevr_Frank__Float64x.cpp
Rsyevr_Frank__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rsyevr_Frank__Float64x_LDFLAGS = $(_Float64x_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Cheev_test_double Rsyev_Frank_double Rsyev_test_double Rsyevd_DingDong_double Rsyevd_Frank_double Rsyevr_Frank_double

double_cxxflags = $(OPENMP_CXXFLAGS)
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double -L$(top_builddir)/mpblas/reference -lmpblas_double

Cheev_test_double_SOURCES = Cheev_test_double.cpp
Cheev_test_double_CXXFLAGS = $(double_cxxflags)
Cheev_test_double_LDFLAGS = $(double_libdepends)

Rsyev_Frank_double_SOURCES = Rsyev_Frank_double.cpp
Rsyev_Frank_double_CXXFLAGS = $(double_cxxflags)
Rsyev_Frank_double_LDFLAGS = $(double_libdepends)

Rsyev_test_double_SOURCES = Rsyev_test_double.cpp
Rsyev_test_double_CXXFLAGS = $(double_cxxflags)
Rsyev_test_double_LDFLAGS = $(double_libdepends)

Rsyevd_DingDong_double_SOURCES = Rsyevd_DingDong_double.cpp
Rsyevd_DingDong_double_CXXFLAGS = $(double_cxxflags)
Rsyevd_DingDong_double_LDFLAGS = $(double_libdepends)

Rsyevd_Frank_double_SOURCES = Rsyevd_Frank_double.cpp
Rsyevd_Frank_double_CXXFLAGS = $(double_cxxflags)
Rsyevd_Frank_double_LDFLAGS = $(double_libdepends)

Rsyevr_Frank_double_SOURCES = Rsyevr_Frank_double.cpp
Rsyevr_Frank_double_CXXFLAGS = $(double_cxxflags)
Rsyevr_Frank_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Cheev_test_dd Rsyev_Frank_dd Rsyev_test_dd Rsyevd_DingDong_dd Rsyevd_Frank_dd Rsyevr_Frank_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(top_builddir)/mpblas/reference -lmpblas_dd -L$(QD_LIBDIR) -lqd

Cheev_test_dd_SOURCES = Cheev_test_dd.cpp
Cheev_test_dd_CXXFLAGS = $(dd_cxxflags)
Cheev_test_dd_LDFLAGS = $(dd_libdepends)

Rsyev_Frank_dd_SOURCES = Rsyev_Frank_dd.cpp
Rsyev_Frank_dd_CXXFLAGS = $(dd_cxxflags)
Rsyev_Frank_dd_LDFLAGS = $(dd_libdepends)

Rsyev_test_dd_SOURCES = Rsyev_test_dd.cpp
Rsyev_test_dd_CXXFLAGS = $(dd_cxxflags)
Rsyev_test_dd_LDFLAGS = $(dd_libdepends)

Rsyevd_DingDong_dd_SOURCES = Rsyevd_DingDong_dd.cpp
Rsyevd_DingDong_dd_CXXFLAGS = $(dd_cxxflags)
Rsyevd_DingDong_dd_LDFLAGS = $(dd_libdepends)

Rsyevd_Frank_dd_SOURCES = Rsyevd_Frank_dd.cpp
Rsyevd_Frank_dd_CXXFLAGS = $(dd_cxxflags)
Rsyevd_Frank_dd_LDFLAGS = $(dd_libdepends)

Rsyevr_Frank_dd_SOURCES = Rsyevr_Frank_dd.cpp
Rsyevr_Frank_dd_CXXFLAGS = $(dd_cxxflags)
Rsyevr_Frank_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Cheev_test_qd Rsyev_Frank_qd Rsyev_test_qd Rsyevd_DingDong_qd Rsyevd_Frank_qd Rsyevr_Frank_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(top_builddir)/mpblas/reference -lmpblas_qd -L$(QD_LIBDIR) -lqd

Cheev_test_qd_SOURCES = Cheev_test_qd.cpp
Cheev_test_qd_CXXFLAGS = $(qd_cxxflags)
Cheev_test_qd_LDFLAGS = $(qd_libdepends)

Rsyev_Frank_qd_SOURCES = Rsyev_Frank_qd.cpp
Rsyev_Frank_qd_CXXFLAGS = $(qd_cxxflags)
Rsyev_Frank_qd_LDFLAGS = $(qd_libdepends)

Rsyev_test_qd_SOURCES = Rsyev_test_qd.cpp
Rsyev_test_qd_CXXFLAGS = $(qd_cxxflags)
Rsyev_test_qd_LDFLAGS = $(qd_libdepends)

Rsyevd_DingDong_qd_SOURCES = Rsyevd_DingDong_qd.cpp
Rsyevd_DingDong_qd_CXXFLAGS = $(qd_cxxflags)
Rsyevd_DingDong_qd_LDFLAGS = $(qd_libdepends)

Rsyevd_Frank_qd_SOURCES = Rsyevd_Frank_qd.cpp
Rsyevd_Frank_qd_CXXFLAGS = $(qd_cxxflags)
Rsyevd_Frank_qd_LDFLAGS = $(qd_libdepends)

Rsyevr_Frank_qd_SOURCES = Rsyevr_Frank_qd.cpp
Rsyevr_Frank_qd_CXXFLAGS = $(qd_cxxflags)
Rsyevr_Frank_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/examples/mplapack/03_SymmetricEigenproblems

mplapackexamples_DATA = Cheev_test_mpfr.cpp Cheev_test_gmp.cpp Cheev_test__Float128.cpp Cheev_test__Float64x.cpp Cheev_test_double.cpp Cheev_test_dd.cpp Cheev_test_qd.cpp Rsyev_Frank_mpfr.cpp Rsyev_Frank_gmp.cpp Rsyev_Frank__Float128.cpp Rsyev_Frank__Float64x.cpp Rsyev_Frank_double.cpp Rsyev_Frank_dd.cpp Rsyev_Frank_qd.cpp Rsyev_test_mpfr.cpp Rsyev_test_gmp.cpp Rsyev_test__Float128.cpp Rsyev_test__Float64x.cpp Rsyev_test_double.cpp Rsyev_test_dd.cpp Rsyev_test_qd.cpp Rsyevd_DingDong_mpfr.cpp Rsyevd_DingDong_gmp.cpp Rsyevd_DingDong__Float128.cpp Rsyevd_DingDong__Float64x.cpp Rsyevd_DingDong_double.cpp Rsyevd_DingDong_dd.cpp Rsyevd_DingDong_qd.cpp Rsyevd_Frank_mpfr.cpp Rsyevd_Frank_gmp.cpp Rsyevd_Frank__Float128.cpp Rsyevd_Frank__Float64x.cpp Rsyevd_Frank_double.cpp Rsyevd_Frank_dd.cpp Rsyevd_Frank_qd.cpp Rsyevr_Frank_mpfr.cpp Rsyevr_Frank_gmp.cpp Rsyevr_Frank__Float128.cpp Rsyevr_Frank__Float64x.cpp Rsyevr_Frank_double.cpp Rsyevr_Frank_dd.cpp Rsyevr_Frank_qd.cpp

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapackexamplesdir) $(prefix)
endif
