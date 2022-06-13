mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Rgesvd_random_highcond_mpfr Rgesvd_readfromfile_mpfr Rgesvd_test2_mpfr Rgesvd_test_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Rgesvd_random_highcond_mpfr_SOURCES = Rgesvd_random_highcond_mpfr.cpp
Rgesvd_random_highcond_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgesvd_random_highcond_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgesvd_readfromfile_mpfr_SOURCES = Rgesvd_readfromfile_mpfr.cpp
Rgesvd_readfromfile_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgesvd_readfromfile_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgesvd_test2_mpfr_SOURCES = Rgesvd_test2_mpfr.cpp
Rgesvd_test2_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgesvd_test2_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgesvd_test_mpfr_SOURCES = Rgesvd_test_mpfr.cpp
Rgesvd_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgesvd_test_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Rgesvd_random_highcond_gmp Rgesvd_readfromfile_gmp Rgesvd_test2_gmp Rgesvd_test_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp -L$(GMP_LIBDIR) -lgmp

Rgesvd_random_highcond_gmp_SOURCES = Rgesvd_random_highcond_gmp.cpp
Rgesvd_random_highcond_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgesvd_random_highcond_gmp_LDFLAGS = $(gmp_libdepends)

Rgesvd_readfromfile_gmp_SOURCES = Rgesvd_readfromfile_gmp.cpp
Rgesvd_readfromfile_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgesvd_readfromfile_gmp_LDFLAGS = $(gmp_libdepends)

Rgesvd_test2_gmp_SOURCES = Rgesvd_test2_gmp.cpp
Rgesvd_test2_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgesvd_test2_gmp_LDFLAGS = $(gmp_libdepends)

Rgesvd_test_gmp_SOURCES = Rgesvd_test_gmp.cpp
Rgesvd_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgesvd_test_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE__FLOAT128
mplapackexamples_PROGRAMS += Rgesvd_random_highcond__Float128 Rgesvd_readfromfile__Float128 Rgesvd_test2__Float128 Rgesvd_test__Float128

_Float128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if WANT_QUADMATH
_Float128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128 -lquadmath
else
_Float128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128
endif

Rgesvd_random_highcond__Float128_SOURCES = Rgesvd_random_highcond__Float128.cpp
Rgesvd_random_highcond__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgesvd_random_highcond__Float128_LDFLAGS = $(_Float128_libdepends)

Rgesvd_readfromfile__Float128_SOURCES = Rgesvd_readfromfile__Float128.cpp
Rgesvd_readfromfile__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgesvd_readfromfile__Float128_LDFLAGS = $(_Float128_libdepends)

Rgesvd_test2__Float128_SOURCES = Rgesvd_test2__Float128.cpp
Rgesvd_test2__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgesvd_test2__Float128_LDFLAGS = $(_Float128_libdepends)

Rgesvd_test__Float128_SOURCES = Rgesvd_test__Float128.cpp
Rgesvd_test__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgesvd_test__Float128_LDFLAGS = $(_Float128_libdepends)

endif

if ENABLE__FLOAT64X
mplapackexamples_PROGRAMS += Rgesvd_random_highcond__Float64x Rgesvd_readfromfile__Float64x Rgesvd_test2__Float64x Rgesvd_test__Float64x

_Float64x_cxxflags = $(OPENMP_CXXFLAGS)
_Float64x_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack__Float64x -L$(top_builddir)/mpblas/reference -lmpblas__Float64x

Rgesvd_random_highcond__Float64x_SOURCES = Rgesvd_random_highcond__Float64x.cpp
Rgesvd_random_highcond__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgesvd_random_highcond__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgesvd_readfromfile__Float64x_SOURCES = Rgesvd_readfromfile__Float64x.cpp
Rgesvd_readfromfile__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgesvd_readfromfile__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgesvd_test2__Float64x_SOURCES = Rgesvd_test2__Float64x.cpp
Rgesvd_test2__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgesvd_test2__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgesvd_test__Float64x_SOURCES = Rgesvd_test__Float64x.cpp
Rgesvd_test__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgesvd_test__Float64x_LDFLAGS = $(_Float64x_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Rgesvd_random_highcond_double Rgesvd_readfromfile_double Rgesvd_test2_double Rgesvd_test_double

double_cxxflags = $(OPENMP_CXXFLAGS)
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double -L$(top_builddir)/mpblas/reference -lmpblas_double

Rgesvd_random_highcond_double_SOURCES = Rgesvd_random_highcond_double.cpp
Rgesvd_random_highcond_double_CXXFLAGS = $(double_cxxflags)
Rgesvd_random_highcond_double_LDFLAGS = $(double_libdepends)

Rgesvd_readfromfile_double_SOURCES = Rgesvd_readfromfile_double.cpp
Rgesvd_readfromfile_double_CXXFLAGS = $(double_cxxflags)
Rgesvd_readfromfile_double_LDFLAGS = $(double_libdepends)

Rgesvd_test2_double_SOURCES = Rgesvd_test2_double.cpp
Rgesvd_test2_double_CXXFLAGS = $(double_cxxflags)
Rgesvd_test2_double_LDFLAGS = $(double_libdepends)

Rgesvd_test_double_SOURCES = Rgesvd_test_double.cpp
Rgesvd_test_double_CXXFLAGS = $(double_cxxflags)
Rgesvd_test_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Rgesvd_random_highcond_dd Rgesvd_readfromfile_dd Rgesvd_test2_dd Rgesvd_test_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(top_builddir)/mpblas/reference -lmpblas_dd -L$(QD_LIBDIR) -lqd

Rgesvd_random_highcond_dd_SOURCES = Rgesvd_random_highcond_dd.cpp
Rgesvd_random_highcond_dd_CXXFLAGS = $(dd_cxxflags)
Rgesvd_random_highcond_dd_LDFLAGS = $(dd_libdepends)

Rgesvd_readfromfile_dd_SOURCES = Rgesvd_readfromfile_dd.cpp
Rgesvd_readfromfile_dd_CXXFLAGS = $(dd_cxxflags)
Rgesvd_readfromfile_dd_LDFLAGS = $(dd_libdepends)

Rgesvd_test2_dd_SOURCES = Rgesvd_test2_dd.cpp
Rgesvd_test2_dd_CXXFLAGS = $(dd_cxxflags)
Rgesvd_test2_dd_LDFLAGS = $(dd_libdepends)

Rgesvd_test_dd_SOURCES = Rgesvd_test_dd.cpp
Rgesvd_test_dd_CXXFLAGS = $(dd_cxxflags)
Rgesvd_test_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Rgesvd_random_highcond_qd Rgesvd_readfromfile_qd Rgesvd_test2_qd Rgesvd_test_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(top_builddir)/mpblas/reference -lmpblas_qd -L$(QD_LIBDIR) -lqd

Rgesvd_random_highcond_qd_SOURCES = Rgesvd_random_highcond_qd.cpp
Rgesvd_random_highcond_qd_CXXFLAGS = $(qd_cxxflags)
Rgesvd_random_highcond_qd_LDFLAGS = $(qd_libdepends)

Rgesvd_readfromfile_qd_SOURCES = Rgesvd_readfromfile_qd.cpp
Rgesvd_readfromfile_qd_CXXFLAGS = $(qd_cxxflags)
Rgesvd_readfromfile_qd_LDFLAGS = $(qd_libdepends)

Rgesvd_test2_qd_SOURCES = Rgesvd_test2_qd.cpp
Rgesvd_test2_qd_CXXFLAGS = $(qd_cxxflags)
Rgesvd_test2_qd_LDFLAGS = $(qd_libdepends)

Rgesvd_test_qd_SOURCES = Rgesvd_test_qd.cpp
Rgesvd_test_qd_CXXFLAGS = $(qd_cxxflags)
Rgesvd_test_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/mplapack/examples/05_SingularValueDecomposition

mplapackexamples_DATA = Rgesvd_random_highcond_mpfr.cpp Rgesvd_random_highcond_gmp.cpp Rgesvd_random_highcond__Float128.cpp Rgesvd_random_highcond__Float64x.cpp Rgesvd_random_highcond_double.cpp Rgesvd_random_highcond_dd.cpp Rgesvd_random_highcond_qd.cpp Rgesvd_readfromfile_mpfr.cpp Rgesvd_readfromfile_gmp.cpp Rgesvd_readfromfile__Float128.cpp Rgesvd_readfromfile__Float64x.cpp Rgesvd_readfromfile_double.cpp Rgesvd_readfromfile_dd.cpp Rgesvd_readfromfile_qd.cpp Rgesvd_test2_mpfr.cpp Rgesvd_test2_gmp.cpp Rgesvd_test2__Float128.cpp Rgesvd_test2__Float64x.cpp Rgesvd_test2_double.cpp Rgesvd_test2_dd.cpp Rgesvd_test2_qd.cpp Rgesvd_test_mpfr.cpp Rgesvd_test_gmp.cpp Rgesvd_test__Float128.cpp Rgesvd_test__Float64x.cpp Rgesvd_test_double.cpp Rgesvd_test_dd.cpp Rgesvd_test_qd.cpp

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapackexamplesdir) $(prefix)
endif
