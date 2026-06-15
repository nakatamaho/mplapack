mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Cgees_test_mpfr Cgeev_NPR_mpfr Cgeev_test2_mpfr Cgeev_test_mpfr Rgees_Grcar_mpfr Rgees_readfromfile_mpfr Rgees_test_mpfr Rgeev_Frank_mpfr Rgeev_random_highcond_mpfr Rgeev_readfromfile_mpfr Rgeev_test_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Cgees_test_mpfr_SOURCES = Cgees_test_mpfr.cpp
Cgees_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgees_test_mpfr_LDFLAGS = $(mpfr_libdepends)

Cgeev_NPR_mpfr_SOURCES = Cgeev_NPR_mpfr.cpp
Cgeev_NPR_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgeev_NPR_mpfr_LDFLAGS = $(mpfr_libdepends)

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
mplapackexamples_PROGRAMS += Cgees_test_gmp Cgeev_test2_gmp Cgeev_test_gmp Rgees_Grcar_gmp Rgees_readfromfile_gmp Rgees_test_gmp Rgeev_Frank_gmp Rgeev_random_highcond_gmp Rgeev_readfromfile_gmp Rgeev_test_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp -L$(GMP_LIBDIR) -lgmp

Cgees_test_gmp_SOURCES = Cgees_test_gmp.cpp
Cgees_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Cgees_test_gmp_LDFLAGS = $(gmp_libdepends)

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

if ENABLE_BINARY128
mplapackexamples_PROGRAMS += Cgees_test_binary128 Cgeev_NPR_binary128 Cgeev_test2_binary128 Cgeev_test_binary128 Rgees_Grcar_binary128 Rgees_readfromfile_binary128 Rgees_test_binary128 Rgeev_Frank_binary128 Rgeev_random_highcond_binary128 Rgeev_readfromfile_binary128 Rgeev_test_binary128

binary128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if MPLAPACK_BINARY128_MODE_QUADMATH
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -L$(top_builddir)/mpblas/reference -lmpblas_binary128 -lquadmath
else
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -L$(top_builddir)/mpblas/reference -lmpblas_binary128
endif

Cgees_test_binary128_SOURCES = Cgees_test_binary128.cpp
Cgees_test_binary128_CXXFLAGS = $(binary128_cxxflags)
Cgees_test_binary128_LDFLAGS = $(binary128_libdepends)

Cgeev_NPR_binary128_SOURCES = Cgeev_NPR_binary128.cpp
Cgeev_NPR_binary128_CXXFLAGS = $(binary128_cxxflags)
Cgeev_NPR_binary128_LDFLAGS = $(binary128_libdepends)

Cgeev_test2_binary128_SOURCES = Cgeev_test2_binary128.cpp
Cgeev_test2_binary128_CXXFLAGS = $(binary128_cxxflags)
Cgeev_test2_binary128_LDFLAGS = $(binary128_libdepends)

Cgeev_test_binary128_SOURCES = Cgeev_test_binary128.cpp
Cgeev_test_binary128_CXXFLAGS = $(binary128_cxxflags)
Cgeev_test_binary128_LDFLAGS = $(binary128_libdepends)

Rgees_Grcar_binary128_SOURCES = Rgees_Grcar_binary128.cpp
Rgees_Grcar_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgees_Grcar_binary128_LDFLAGS = $(binary128_libdepends)

Rgees_readfromfile_binary128_SOURCES = Rgees_readfromfile_binary128.cpp
Rgees_readfromfile_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgees_readfromfile_binary128_LDFLAGS = $(binary128_libdepends)

Rgees_test_binary128_SOURCES = Rgees_test_binary128.cpp
Rgees_test_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgees_test_binary128_LDFLAGS = $(binary128_libdepends)

Rgeev_Frank_binary128_SOURCES = Rgeev_Frank_binary128.cpp
Rgeev_Frank_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgeev_Frank_binary128_LDFLAGS = $(binary128_libdepends)

Rgeev_random_highcond_binary128_SOURCES = Rgeev_random_highcond_binary128.cpp
Rgeev_random_highcond_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgeev_random_highcond_binary128_LDFLAGS = $(binary128_libdepends)

Rgeev_readfromfile_binary128_SOURCES = Rgeev_readfromfile_binary128.cpp
Rgeev_readfromfile_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgeev_readfromfile_binary128_LDFLAGS = $(binary128_libdepends)

Rgeev_test_binary128_SOURCES = Rgeev_test_binary128.cpp
Rgeev_test_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgeev_test_binary128_LDFLAGS = $(binary128_libdepends)

endif

if ENABLE_BINARY80
mplapackexamples_PROGRAMS += Cgees_test_binary80 Cgeev_NPR_binary80 Cgeev_test2_binary80 Cgeev_test_binary80 Rgees_Grcar_binary80 Rgees_readfromfile_binary80 Rgees_test_binary80 Rgeev_Frank_binary80 Rgeev_random_highcond_binary80 Rgeev_readfromfile_binary80 Rgeev_test_binary80

binary80_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
binary80_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary80 -L$(top_builddir)/mpblas/reference -lmpblas_binary80

Cgees_test_binary80_SOURCES = Cgees_test_binary80.cpp
Cgees_test_binary80_CXXFLAGS = $(binary80_cxxflags)
Cgees_test_binary80_LDFLAGS = $(binary80_libdepends)

Cgeev_NPR_binary80_SOURCES = Cgeev_NPR_binary80.cpp
Cgeev_NPR_binary80_CXXFLAGS = $(binary80_cxxflags)
Cgeev_NPR_binary80_LDFLAGS = $(binary80_libdepends)

Cgeev_test2_binary80_SOURCES = Cgeev_test2_binary80.cpp
Cgeev_test2_binary80_CXXFLAGS = $(binary80_cxxflags)
Cgeev_test2_binary80_LDFLAGS = $(binary80_libdepends)

Cgeev_test_binary80_SOURCES = Cgeev_test_binary80.cpp
Cgeev_test_binary80_CXXFLAGS = $(binary80_cxxflags)
Cgeev_test_binary80_LDFLAGS = $(binary80_libdepends)

Rgees_Grcar_binary80_SOURCES = Rgees_Grcar_binary80.cpp
Rgees_Grcar_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgees_Grcar_binary80_LDFLAGS = $(binary80_libdepends)

Rgees_readfromfile_binary80_SOURCES = Rgees_readfromfile_binary80.cpp
Rgees_readfromfile_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgees_readfromfile_binary80_LDFLAGS = $(binary80_libdepends)

Rgees_test_binary80_SOURCES = Rgees_test_binary80.cpp
Rgees_test_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgees_test_binary80_LDFLAGS = $(binary80_libdepends)

Rgeev_Frank_binary80_SOURCES = Rgeev_Frank_binary80.cpp
Rgeev_Frank_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgeev_Frank_binary80_LDFLAGS = $(binary80_libdepends)

Rgeev_random_highcond_binary80_SOURCES = Rgeev_random_highcond_binary80.cpp
Rgeev_random_highcond_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgeev_random_highcond_binary80_LDFLAGS = $(binary80_libdepends)

Rgeev_readfromfile_binary80_SOURCES = Rgeev_readfromfile_binary80.cpp
Rgeev_readfromfile_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgeev_readfromfile_binary80_LDFLAGS = $(binary80_libdepends)

Rgeev_test_binary80_SOURCES = Rgeev_test_binary80.cpp
Rgeev_test_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgeev_test_binary80_LDFLAGS = $(binary80_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Cgees_test_double Cgeev_NPR_double Cgeev_test2_double Cgeev_test_double Rgees_Grcar_double Rgees_readfromfile_double Rgees_test_double Rgeev_Frank_double Rgeev_random_highcond_double Rgeev_readfromfile_double Rgeev_test_double

double_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double -L$(top_builddir)/mpblas/reference -lmpblas_double

Cgees_test_double_SOURCES = Cgees_test_double.cpp
Cgees_test_double_CXXFLAGS = $(double_cxxflags)
Cgees_test_double_LDFLAGS = $(double_libdepends)

Cgeev_NPR_double_SOURCES = Cgeev_NPR_double.cpp
Cgeev_NPR_double_CXXFLAGS = $(double_cxxflags)
Cgeev_NPR_double_LDFLAGS = $(double_libdepends)

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
mplapackexamples_PROGRAMS += Cgees_test_dd Cgeev_NPR_dd Cgeev_test2_dd Cgeev_test_dd Rgees_Grcar_dd Rgees_readfromfile_dd Rgees_test_dd Rgeev_Frank_dd Rgeev_random_highcond_dd Rgeev_readfromfile_dd Rgeev_test_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(top_builddir)/mpblas/reference -lmpblas_dd -L$(QD_LIBDIR) -lqd

Cgees_test_dd_SOURCES = Cgees_test_dd.cpp
Cgees_test_dd_CXXFLAGS = $(dd_cxxflags)
Cgees_test_dd_LDFLAGS = $(dd_libdepends)

Cgeev_NPR_dd_SOURCES = Cgeev_NPR_dd.cpp
Cgeev_NPR_dd_CXXFLAGS = $(dd_cxxflags)
Cgeev_NPR_dd_LDFLAGS = $(dd_libdepends)

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
mplapackexamples_PROGRAMS += Cgees_test_qd Cgeev_NPR_qd Cgeev_test2_qd Cgeev_test_qd Rgees_Grcar_qd Rgees_readfromfile_qd Rgees_test_qd Rgeev_Frank_qd Rgeev_random_highcond_qd Rgeev_readfromfile_qd Rgeev_test_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(top_builddir)/mpblas/reference -lmpblas_qd -L$(QD_LIBDIR) -lqd

Cgees_test_qd_SOURCES = Cgees_test_qd.cpp
Cgees_test_qd_CXXFLAGS = $(qd_cxxflags)
Cgees_test_qd_LDFLAGS = $(qd_libdepends)

Cgeev_NPR_qd_SOURCES = Cgeev_NPR_qd.cpp
Cgeev_NPR_qd_CXXFLAGS = $(qd_cxxflags)
Cgeev_NPR_qd_LDFLAGS = $(qd_libdepends)

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

mplapackexamples_DATA = Cgees_test_mpfr.cpp Cgees_test_gmp.cpp Cgees_test_binary128.cpp Cgees_test_binary80.cpp Cgees_test_double.cpp Cgees_test_dd.cpp Cgees_test_qd.cpp Cgeev_NPR_mpfr.cpp Cgeev_NPR_binary128.cpp Cgeev_NPR_binary80.cpp Cgeev_NPR_double.cpp Cgeev_NPR_dd.cpp Cgeev_NPR_qd.cpp Cgeev_test2_mpfr.cpp Cgeev_test2_gmp.cpp Cgeev_test2_binary128.cpp Cgeev_test2_binary80.cpp Cgeev_test2_double.cpp Cgeev_test2_dd.cpp Cgeev_test2_qd.cpp Cgeev_test_mpfr.cpp Cgeev_test_gmp.cpp Cgeev_test_binary128.cpp Cgeev_test_binary80.cpp Cgeev_test_double.cpp Cgeev_test_dd.cpp Cgeev_test_qd.cpp Rgees_Grcar_mpfr.cpp Rgees_Grcar_gmp.cpp Rgees_Grcar_binary128.cpp Rgees_Grcar_binary80.cpp Rgees_Grcar_double.cpp Rgees_Grcar_dd.cpp Rgees_Grcar_qd.cpp Rgees_readfromfile_mpfr.cpp Rgees_readfromfile_gmp.cpp Rgees_readfromfile_binary128.cpp Rgees_readfromfile_binary80.cpp Rgees_readfromfile_double.cpp Rgees_readfromfile_dd.cpp Rgees_readfromfile_qd.cpp Rgees_test_mpfr.cpp Rgees_test_gmp.cpp Rgees_test_binary128.cpp Rgees_test_binary80.cpp Rgees_test_double.cpp Rgees_test_dd.cpp Rgees_test_qd.cpp Rgeev_Frank_mpfr.cpp Rgeev_Frank_gmp.cpp Rgeev_Frank_binary128.cpp Rgeev_Frank_binary80.cpp Rgeev_Frank_double.cpp Rgeev_Frank_dd.cpp Rgeev_Frank_qd.cpp Rgeev_random_highcond_mpfr.cpp Rgeev_random_highcond_gmp.cpp Rgeev_random_highcond_binary128.cpp Rgeev_random_highcond_binary80.cpp Rgeev_random_highcond_double.cpp Rgeev_random_highcond_dd.cpp Rgeev_random_highcond_qd.cpp Rgeev_readfromfile_mpfr.cpp Rgeev_readfromfile_gmp.cpp Rgeev_readfromfile_binary128.cpp Rgeev_readfromfile_binary80.cpp Rgeev_readfromfile_double.cpp Rgeev_readfromfile_dd.cpp Rgeev_readfromfile_qd.cpp Rgeev_test_mpfr.cpp Rgeev_test_gmp.cpp Rgeev_test_binary128.cpp Rgeev_test_binary80.cpp Rgeev_test_double.cpp Rgeev_test_dd.cpp Rgeev_test_qd.cpp \
Matrix_01.txt Matrix_02.txt Matrix_03.txt Matrix_04.txt Matrix_05.txt Matrix_06.txt Matrix_07.txt Makefile.freebsd Makefile.linux Makefile.macos Makefile.mingw

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapackexamplesdir) $(prefix)
endif

EXTRA_DIST = \
	Makefile.freebsd.in \
	Makefile.linux.in \
	Makefile.macos.in \
	Makefile.mingw.in

GENERATED_MAKEFILES = \
	Makefile.freebsd \
	Makefile.linux \
	Makefile.macos \
	Makefile.mingw

# Resolve %%LIBQUADMATH%% based on whether quadmath backend is selected.
if MPLAPACK_BINARY128_MODE_QUADMATH
LIBQUADMATH_SUBST = -lquadmath
else
LIBQUADMATH_SUBST =
endif

if MPLAPACK_BINARY128_MODE_QUADMATH
MPLAPACK_SED_SUBST = \
        -e 's|%%MPLAPACKDIR%%|$(prefix)|g' \
        -e 's|%%LIBQUADMATH%%| -lquadmath|g'
else
MPLAPACK_SED_SUBST = \
        -e 's|%%MPLAPACKDIR%%|$(prefix)|g' \
        -e 's|%%LIBQUADMATH%%||g'
endif

DISABLED_EXAMPLE_SUFFIXES =
if ENABLE_MPFR
else
DISABLED_EXAMPLE_SUFFIXES += mpfr
endif
if ENABLE_GMP
else
DISABLED_EXAMPLE_SUFFIXES += gmp
endif
if ENABLE_QD
else
DISABLED_EXAMPLE_SUFFIXES += qd
endif
if ENABLE_DD
else
DISABLED_EXAMPLE_SUFFIXES += dd
endif
if ENABLE_DOUBLE
else
DISABLED_EXAMPLE_SUFFIXES += double
endif
if ENABLE_BINARY128
else
DISABLED_EXAMPLE_SUFFIXES += binary128
endif
if ENABLE_BINARY80
else
DISABLED_EXAMPLE_SUFFIXES += binary80
endif

Makefile.freebsd: Makefile.freebsd.in
	sed $(MPLAPACK_SED_SUBST) $< > $@
	bash $(top_srcdir)/misc/filter-example-programs.sh $@ $(DISABLED_EXAMPLE_SUFFIXES)

Makefile.linux: Makefile.linux.in
	sed $(MPLAPACK_SED_SUBST) $< > $@
	bash $(top_srcdir)/misc/filter-example-programs.sh $@ $(DISABLED_EXAMPLE_SUFFIXES)

Makefile.macos: Makefile.macos.in
	sed $(MPLAPACK_SED_SUBST) $< > $@
	bash $(top_srcdir)/misc/filter-example-programs.sh $@ $(DISABLED_EXAMPLE_SUFFIXES)

Makefile.mingw: Makefile.mingw.in
	sed $(MPLAPACK_SED_SUBST) $< > $@
	bash $(top_srcdir)/misc/filter-example-programs.sh $@ $(DISABLED_EXAMPLE_SUFFIXES)

all-local: $(GENERATED_MAKEFILES)

clean-local:
	rm -f $(GENERATED_MAKEFILES)
