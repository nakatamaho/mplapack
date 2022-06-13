mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Cgeev_test2_mpfr Cgeev_test_mpfr Rgees_Grcar_mpfr Rgees_readfromfile_mpfr Rgees_test_mpfr Rgeev_Frank_mpfr Rgeev_random_highcond_mpfr Rgeev_readfromfile_mpfr Rgeev_test_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Cgeev_test2_mpfr_SOURCES = Cgeev_test2_mpfr.cpp
Cgeev_test2_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgeev_test2_mpfr_LDFLAGS = $(mpfr_libdepends)

Cgeev_test_mpfr_SOURCES = Cgeev_test_mpfr.cpp
Cgeev_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgeev_test_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgees_Grcar_mpfr_SOURCES = Rgees_Grcar_mpfr.cpp
Rgees_Grcar_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgees_Grcar_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgees_readfromfile_mpfr_SOURCES = Rgees_readfromfile_mpfr.cpp
Rgees_readfromfile_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgees_readfromfile_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgees_test_mpfr_SOURCES = Rgees_test_mpfr.cpp
Rgees_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgees_test_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgeev_Frank_mpfr_SOURCES = Rgeev_Frank_mpfr.cpp
Rgeev_Frank_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgeev_Frank_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgeev_random_highcond_mpfr_SOURCES = Rgeev_random_highcond_mpfr.cpp
Rgeev_random_highcond_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgeev_random_highcond_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgeev_readfromfile_mpfr_SOURCES = Rgeev_readfromfile_mpfr.cpp
Rgeev_readfromfile_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgeev_readfromfile_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgeev_test_mpfr_SOURCES = Rgeev_test_mpfr.cpp
Rgeev_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgeev_test_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Cgeev_test2_gmp Cgeev_test_gmp Rgees_Grcar_gmp Rgees_readfromfile_gmp Rgees_test_gmp Rgeev_Frank_gmp Rgeev_random_highcond_gmp Rgeev_readfromfile_gmp Rgeev_test_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp -L$(GMP_LIBDIR) -lgmp

Cgeev_test2_gmp_SOURCES = Cgeev_test2_gmp.cpp
Cgeev_test2_gmp_CXXFLAGS = $(gmp_cxxflags)
Cgeev_test2_gmp_LDFLAGS = $(gmp_libdepends)

Cgeev_test_gmp_SOURCES = Cgeev_test_gmp.cpp
Cgeev_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Cgeev_test_gmp_LDFLAGS = $(gmp_libdepends)

Rgees_Grcar_gmp_SOURCES = Rgees_Grcar_gmp.cpp
Rgees_Grcar_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgees_Grcar_gmp_LDFLAGS = $(gmp_libdepends)

Rgees_readfromfile_gmp_SOURCES = Rgees_readfromfile_gmp.cpp
Rgees_readfromfile_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgees_readfromfile_gmp_LDFLAGS = $(gmp_libdepends)

Rgees_test_gmp_SOURCES = Rgees_test_gmp.cpp
Rgees_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgees_test_gmp_LDFLAGS = $(gmp_libdepends)

Rgeev_Frank_gmp_SOURCES = Rgeev_Frank_gmp.cpp
Rgeev_Frank_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgeev_Frank_gmp_LDFLAGS = $(gmp_libdepends)

Rgeev_random_highcond_gmp_SOURCES = Rgeev_random_highcond_gmp.cpp
Rgeev_random_highcond_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgeev_random_highcond_gmp_LDFLAGS = $(gmp_libdepends)

Rgeev_readfromfile_gmp_SOURCES = Rgeev_readfromfile_gmp.cpp
Rgeev_readfromfile_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgeev_readfromfile_gmp_LDFLAGS = $(gmp_libdepends)

Rgeev_test_gmp_SOURCES = Rgeev_test_gmp.cpp
Rgeev_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgeev_test_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE__FLOAT128
mplapackexamples_PROGRAMS += Cgeev_test2__Float128 Cgeev_test__Float128 Rgees_Grcar__Float128 Rgees_readfromfile__Float128 Rgees_test__Float128 Rgeev_Frank__Float128 Rgeev_random_highcond__Float128 Rgeev_readfromfile__Float128 Rgeev_test__Float128

_Float128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if WANT_QUADMATH
_Float128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128 -lquadmath
else
_Float128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128
endif

Cgeev_test2__Float128_SOURCES = Cgeev_test2__Float128.cpp
Cgeev_test2__Float128_CXXFLAGS = $(_Float128_cxxflags)
Cgeev_test2__Float128_LDFLAGS = $(_Float128_libdepends)

Cgeev_test__Float128_SOURCES = Cgeev_test__Float128.cpp
Cgeev_test__Float128_CXXFLAGS = $(_Float128_cxxflags)
Cgeev_test__Float128_LDFLAGS = $(_Float128_libdepends)

Rgees_Grcar__Float128_SOURCES = Rgees_Grcar__Float128.cpp
Rgees_Grcar__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgees_Grcar__Float128_LDFLAGS = $(_Float128_libdepends)

Rgees_readfromfile__Float128_SOURCES = Rgees_readfromfile__Float128.cpp
Rgees_readfromfile__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgees_readfromfile__Float128_LDFLAGS = $(_Float128_libdepends)

Rgees_test__Float128_SOURCES = Rgees_test__Float128.cpp
Rgees_test__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgees_test__Float128_LDFLAGS = $(_Float128_libdepends)

Rgeev_Frank__Float128_SOURCES = Rgeev_Frank__Float128.cpp
Rgeev_Frank__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgeev_Frank__Float128_LDFLAGS = $(_Float128_libdepends)

Rgeev_random_highcond__Float128_SOURCES = Rgeev_random_highcond__Float128.cpp
Rgeev_random_highcond__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgeev_random_highcond__Float128_LDFLAGS = $(_Float128_libdepends)

Rgeev_readfromfile__Float128_SOURCES = Rgeev_readfromfile__Float128.cpp
Rgeev_readfromfile__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgeev_readfromfile__Float128_LDFLAGS = $(_Float128_libdepends)

Rgeev_test__Float128_SOURCES = Rgeev_test__Float128.cpp
Rgeev_test__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgeev_test__Float128_LDFLAGS = $(_Float128_libdepends)

endif

if ENABLE__FLOAT64X
mplapackexamples_PROGRAMS += Cgeev_test2__Float64x Cgeev_test__Float64x Rgees_Grcar__Float64x Rgees_readfromfile__Float64x Rgees_test__Float64x Rgeev_Frank__Float64x Rgeev_random_highcond__Float64x Rgeev_readfromfile__Float64x Rgeev_test__Float64x

_Float64x_cxxflags = $(OPENMP_CXXFLAGS)
_Float64x_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack__Float64x -L$(top_builddir)/mpblas/reference -lmpblas__Float64x

Cgeev_test2__Float64x_SOURCES = Cgeev_test2__Float64x.cpp
Cgeev_test2__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Cgeev_test2__Float64x_LDFLAGS = $(_Float64x_libdepends)

Cgeev_test__Float64x_SOURCES = Cgeev_test__Float64x.cpp
Cgeev_test__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Cgeev_test__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgees_Grcar__Float64x_SOURCES = Rgees_Grcar__Float64x.cpp
Rgees_Grcar__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgees_Grcar__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgees_readfromfile__Float64x_SOURCES = Rgees_readfromfile__Float64x.cpp
Rgees_readfromfile__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgees_readfromfile__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgees_test__Float64x_SOURCES = Rgees_test__Float64x.cpp
Rgees_test__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgees_test__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgeev_Frank__Float64x_SOURCES = Rgeev_Frank__Float64x.cpp
Rgeev_Frank__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgeev_Frank__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgeev_random_highcond__Float64x_SOURCES = Rgeev_random_highcond__Float64x.cpp
Rgeev_random_highcond__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgeev_random_highcond__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgeev_readfromfile__Float64x_SOURCES = Rgeev_readfromfile__Float64x.cpp
Rgeev_readfromfile__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgeev_readfromfile__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgeev_test__Float64x_SOURCES = Rgeev_test__Float64x.cpp
Rgeev_test__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgeev_test__Float64x_LDFLAGS = $(_Float64x_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Cgeev_test2_double Cgeev_test_double Rgees_Grcar_double Rgees_readfromfile_double Rgees_test_double Rgeev_Frank_double Rgeev_random_highcond_double Rgeev_readfromfile_double Rgeev_test_double

double_cxxflags = $(OPENMP_CXXFLAGS)
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double -L$(top_builddir)/mpblas/reference -lmpblas_double

Cgeev_test2_double_SOURCES = Cgeev_test2_double.cpp
Cgeev_test2_double_CXXFLAGS = $(double_cxxflags)
Cgeev_test2_double_LDFLAGS = $(double_libdepends)

Cgeev_test_double_SOURCES = Cgeev_test_double.cpp
Cgeev_test_double_CXXFLAGS = $(double_cxxflags)
Cgeev_test_double_LDFLAGS = $(double_libdepends)

Rgees_Grcar_double_SOURCES = Rgees_Grcar_double.cpp
Rgees_Grcar_double_CXXFLAGS = $(double_cxxflags)
Rgees_Grcar_double_LDFLAGS = $(double_libdepends)

Rgees_readfromfile_double_SOURCES = Rgees_readfromfile_double.cpp
Rgees_readfromfile_double_CXXFLAGS = $(double_cxxflags)
Rgees_readfromfile_double_LDFLAGS = $(double_libdepends)

Rgees_test_double_SOURCES = Rgees_test_double.cpp
Rgees_test_double_CXXFLAGS = $(double_cxxflags)
Rgees_test_double_LDFLAGS = $(double_libdepends)

Rgeev_Frank_double_SOURCES = Rgeev_Frank_double.cpp
Rgeev_Frank_double_CXXFLAGS = $(double_cxxflags)
Rgeev_Frank_double_LDFLAGS = $(double_libdepends)

Rgeev_random_highcond_double_SOURCES = Rgeev_random_highcond_double.cpp
Rgeev_random_highcond_double_CXXFLAGS = $(double_cxxflags)
Rgeev_random_highcond_double_LDFLAGS = $(double_libdepends)

Rgeev_readfromfile_double_SOURCES = Rgeev_readfromfile_double.cpp
Rgeev_readfromfile_double_CXXFLAGS = $(double_cxxflags)
Rgeev_readfromfile_double_LDFLAGS = $(double_libdepends)

Rgeev_test_double_SOURCES = Rgeev_test_double.cpp
Rgeev_test_double_CXXFLAGS = $(double_cxxflags)
Rgeev_test_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Cgeev_test2_dd Cgeev_test_dd Rgees_Grcar_dd Rgees_readfromfile_dd Rgees_test_dd Rgeev_Frank_dd Rgeev_random_highcond_dd Rgeev_readfromfile_dd Rgeev_test_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(top_builddir)/mpblas/reference -lmpblas_dd -L$(QD_LIBDIR) -lqd

Cgeev_test2_dd_SOURCES = Cgeev_test2_dd.cpp
Cgeev_test2_dd_CXXFLAGS = $(dd_cxxflags)
Cgeev_test2_dd_LDFLAGS = $(dd_libdepends)

Cgeev_test_dd_SOURCES = Cgeev_test_dd.cpp
Cgeev_test_dd_CXXFLAGS = $(dd_cxxflags)
Cgeev_test_dd_LDFLAGS = $(dd_libdepends)

Rgees_Grcar_dd_SOURCES = Rgees_Grcar_dd.cpp
Rgees_Grcar_dd_CXXFLAGS = $(dd_cxxflags)
Rgees_Grcar_dd_LDFLAGS = $(dd_libdepends)

Rgees_readfromfile_dd_SOURCES = Rgees_readfromfile_dd.cpp
Rgees_readfromfile_dd_CXXFLAGS = $(dd_cxxflags)
Rgees_readfromfile_dd_LDFLAGS = $(dd_libdepends)

Rgees_test_dd_SOURCES = Rgees_test_dd.cpp
Rgees_test_dd_CXXFLAGS = $(dd_cxxflags)
Rgees_test_dd_LDFLAGS = $(dd_libdepends)

Rgeev_Frank_dd_SOURCES = Rgeev_Frank_dd.cpp
Rgeev_Frank_dd_CXXFLAGS = $(dd_cxxflags)
Rgeev_Frank_dd_LDFLAGS = $(dd_libdepends)

Rgeev_random_highcond_dd_SOURCES = Rgeev_random_highcond_dd.cpp
Rgeev_random_highcond_dd_CXXFLAGS = $(dd_cxxflags)
Rgeev_random_highcond_dd_LDFLAGS = $(dd_libdepends)

Rgeev_readfromfile_dd_SOURCES = Rgeev_readfromfile_dd.cpp
Rgeev_readfromfile_dd_CXXFLAGS = $(dd_cxxflags)
Rgeev_readfromfile_dd_LDFLAGS = $(dd_libdepends)

Rgeev_test_dd_SOURCES = Rgeev_test_dd.cpp
Rgeev_test_dd_CXXFLAGS = $(dd_cxxflags)
Rgeev_test_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Cgeev_test2_qd Cgeev_test_qd Rgees_Grcar_qd Rgees_readfromfile_qd Rgees_test_qd Rgeev_Frank_qd Rgeev_random_highcond_qd Rgeev_readfromfile_qd Rgeev_test_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(top_builddir)/mpblas/reference -lmpblas_qd -L$(QD_LIBDIR) -lqd

Cgeev_test2_qd_SOURCES = Cgeev_test2_qd.cpp
Cgeev_test2_qd_CXXFLAGS = $(qd_cxxflags)
Cgeev_test2_qd_LDFLAGS = $(qd_libdepends)

Cgeev_test_qd_SOURCES = Cgeev_test_qd.cpp
Cgeev_test_qd_CXXFLAGS = $(qd_cxxflags)
Cgeev_test_qd_LDFLAGS = $(qd_libdepends)

Rgees_Grcar_qd_SOURCES = Rgees_Grcar_qd.cpp
Rgees_Grcar_qd_CXXFLAGS = $(qd_cxxflags)
Rgees_Grcar_qd_LDFLAGS = $(qd_libdepends)

Rgees_readfromfile_qd_SOURCES = Rgees_readfromfile_qd.cpp
Rgees_readfromfile_qd_CXXFLAGS = $(qd_cxxflags)
Rgees_readfromfile_qd_LDFLAGS = $(qd_libdepends)

Rgees_test_qd_SOURCES = Rgees_test_qd.cpp
Rgees_test_qd_CXXFLAGS = $(qd_cxxflags)
Rgees_test_qd_LDFLAGS = $(qd_libdepends)

Rgeev_Frank_qd_SOURCES = Rgeev_Frank_qd.cpp
Rgeev_Frank_qd_CXXFLAGS = $(qd_cxxflags)
Rgeev_Frank_qd_LDFLAGS = $(qd_libdepends)

Rgeev_random_highcond_qd_SOURCES = Rgeev_random_highcond_qd.cpp
Rgeev_random_highcond_qd_CXXFLAGS = $(qd_cxxflags)
Rgeev_random_highcond_qd_LDFLAGS = $(qd_libdepends)

Rgeev_readfromfile_qd_SOURCES = Rgeev_readfromfile_qd.cpp
Rgeev_readfromfile_qd_CXXFLAGS = $(qd_cxxflags)
Rgeev_readfromfile_qd_LDFLAGS = $(qd_libdepends)

Rgeev_test_qd_SOURCES = Rgeev_test_qd.cpp
Rgeev_test_qd_CXXFLAGS = $(qd_cxxflags)
Rgeev_test_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/examples/mplapack/04_NonsymmetricEigenproblems

mplapackexamples_DATA = Cgeev_test2_mpfr.cpp Cgeev_test2_gmp.cpp Cgeev_test2__Float128.cpp Cgeev_test2__Float64x.cpp Cgeev_test2_double.cpp Cgeev_test2_dd.cpp Cgeev_test2_qd.cpp Cgeev_test_mpfr.cpp Cgeev_test_gmp.cpp Cgeev_test__Float128.cpp Cgeev_test__Float64x.cpp Cgeev_test_double.cpp Cgeev_test_dd.cpp Cgeev_test_qd.cpp Rgees_Grcar_mpfr.cpp Rgees_Grcar_gmp.cpp Rgees_Grcar__Float128.cpp Rgees_Grcar__Float64x.cpp Rgees_Grcar_double.cpp Rgees_Grcar_dd.cpp Rgees_Grcar_qd.cpp Rgees_readfromfile_mpfr.cpp Rgees_readfromfile_gmp.cpp Rgees_readfromfile__Float128.cpp Rgees_readfromfile__Float64x.cpp Rgees_readfromfile_double.cpp Rgees_readfromfile_dd.cpp Rgees_readfromfile_qd.cpp Rgees_test_mpfr.cpp Rgees_test_gmp.cpp Rgees_test__Float128.cpp Rgees_test__Float64x.cpp Rgees_test_double.cpp Rgees_test_dd.cpp Rgees_test_qd.cpp Rgeev_Frank_mpfr.cpp Rgeev_Frank_gmp.cpp Rgeev_Frank__Float128.cpp Rgeev_Frank__Float64x.cpp Rgeev_Frank_double.cpp Rgeev_Frank_dd.cpp Rgeev_Frank_qd.cpp Rgeev_random_highcond_mpfr.cpp Rgeev_random_highcond_gmp.cpp Rgeev_random_highcond__Float128.cpp Rgeev_random_highcond__Float64x.cpp Rgeev_random_highcond_double.cpp Rgeev_random_highcond_dd.cpp Rgeev_random_highcond_qd.cpp Rgeev_readfromfile_mpfr.cpp Rgeev_readfromfile_gmp.cpp Rgeev_readfromfile__Float128.cpp Rgeev_readfromfile__Float64x.cpp Rgeev_readfromfile_double.cpp Rgeev_readfromfile_dd.cpp Rgeev_readfromfile_qd.cpp Rgeev_test_mpfr.cpp Rgeev_test_gmp.cpp Rgeev_test__Float128.cpp Rgeev_test__Float64x.cpp Rgeev_test_double.cpp Rgeev_test_dd.cpp Rgeev_test_qd.cpp

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapackexamplesdir) $(prefix)
endif
