mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Cgesv_example_mpfr Cgetrf_getrs_example_mpfr Cgetri_example_mpfr Rgbsv_example_mpfr Rgecon_example_mpfr Rgesv_example_mpfr Rgesvx_example_mpfr Rgetrf_getrs_example_mpfr Rgetri_example2_mpfr Rgetri_example3_mpfr Rgetri_example_mpfr Rgtsv_example_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Cgesv_example_mpfr_SOURCES = Cgesv_example_mpfr.cpp
Cgesv_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgesv_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Cgetrf_getrs_example_mpfr_SOURCES = Cgetrf_getrs_example_mpfr.cpp
Cgetrf_getrs_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgetrf_getrs_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Cgetri_example_mpfr_SOURCES = Cgetri_example_mpfr.cpp
Cgetri_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgetri_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgbsv_example_mpfr_SOURCES = Rgbsv_example_mpfr.cpp
Rgbsv_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgbsv_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgecon_example_mpfr_SOURCES = Rgecon_example_mpfr.cpp
Rgecon_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgecon_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgesv_example_mpfr_SOURCES = Rgesv_example_mpfr.cpp
Rgesv_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgesv_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgesvx_example_mpfr_SOURCES = Rgesvx_example_mpfr.cpp
Rgesvx_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgesvx_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgetrf_getrs_example_mpfr_SOURCES = Rgetrf_getrs_example_mpfr.cpp
Rgetrf_getrs_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgetrf_getrs_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgetri_example2_mpfr_SOURCES = Rgetri_example2_mpfr.cpp
Rgetri_example2_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgetri_example2_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgetri_example3_mpfr_SOURCES = Rgetri_example3_mpfr.cpp
Rgetri_example3_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgetri_example3_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgetri_example_mpfr_SOURCES = Rgetri_example_mpfr.cpp
Rgetri_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgetri_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgtsv_example_mpfr_SOURCES = Rgtsv_example_mpfr.cpp
Rgtsv_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgtsv_example_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Cgesv_example_gmp Cgetrf_getrs_example_gmp Cgetri_example_gmp Rgbsv_example_gmp Rgecon_example_gmp Rgesv_example_gmp Rgesvx_example_gmp Rgetrf_getrs_example_gmp Rgetri_example2_gmp Rgetri_example3_gmp Rgetri_example_gmp Rgtsv_example_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(GMP_LIBDIR) -lgmp

Cgesv_example_gmp_SOURCES = Cgesv_example_gmp.cpp
Cgesv_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Cgesv_example_gmp_LDFLAGS = $(gmp_libdepends)

Cgetrf_getrs_example_gmp_SOURCES = Cgetrf_getrs_example_gmp.cpp
Cgetrf_getrs_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Cgetrf_getrs_example_gmp_LDFLAGS = $(gmp_libdepends)

Cgetri_example_gmp_SOURCES = Cgetri_example_gmp.cpp
Cgetri_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Cgetri_example_gmp_LDFLAGS = $(gmp_libdepends)

Rgbsv_example_gmp_SOURCES = Rgbsv_example_gmp.cpp
Rgbsv_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgbsv_example_gmp_LDFLAGS = $(gmp_libdepends)

Rgecon_example_gmp_SOURCES = Rgecon_example_gmp.cpp
Rgecon_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgecon_example_gmp_LDFLAGS = $(gmp_libdepends)

Rgesv_example_gmp_SOURCES = Rgesv_example_gmp.cpp
Rgesv_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgesv_example_gmp_LDFLAGS = $(gmp_libdepends)

Rgesvx_example_gmp_SOURCES = Rgesvx_example_gmp.cpp
Rgesvx_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgesvx_example_gmp_LDFLAGS = $(gmp_libdepends)

Rgetrf_getrs_example_gmp_SOURCES = Rgetrf_getrs_example_gmp.cpp
Rgetrf_getrs_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgetrf_getrs_example_gmp_LDFLAGS = $(gmp_libdepends)

Rgetri_example2_gmp_SOURCES = Rgetri_example2_gmp.cpp
Rgetri_example2_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgetri_example2_gmp_LDFLAGS = $(gmp_libdepends)

Rgetri_example3_gmp_SOURCES = Rgetri_example3_gmp.cpp
Rgetri_example3_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgetri_example3_gmp_LDFLAGS = $(gmp_libdepends)

Rgetri_example_gmp_SOURCES = Rgetri_example_gmp.cpp
Rgetri_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgetri_example_gmp_LDFLAGS = $(gmp_libdepends)

Rgtsv_example_gmp_SOURCES = Rgtsv_example_gmp.cpp
Rgtsv_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgtsv_example_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE_BINARY128
mplapackexamples_PROGRAMS += Cgesv_example_binary128 Cgetrf_getrs_example_binary128 Cgetri_example_binary128 Rgbsv_example_binary128 Rgecon_example_binary128 Rgesv_example_binary128 Rgesvx_example_binary128 Rgetrf_getrs_example_binary128 Rgetri_example2_binary128 Rgetri_example3_binary128 Rgetri_example_binary128 Rgtsv_example_binary128

binary128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if MPLAPACK_BINARY128_MODE_QUADMATH
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -lquadmath
else
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128
endif

Cgesv_example_binary128_SOURCES = Cgesv_example_binary128.cpp
Cgesv_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Cgesv_example_binary128_LDFLAGS = $(binary128_libdepends)

Cgetrf_getrs_example_binary128_SOURCES = Cgetrf_getrs_example_binary128.cpp
Cgetrf_getrs_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Cgetrf_getrs_example_binary128_LDFLAGS = $(binary128_libdepends)

Cgetri_example_binary128_SOURCES = Cgetri_example_binary128.cpp
Cgetri_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Cgetri_example_binary128_LDFLAGS = $(binary128_libdepends)

Rgbsv_example_binary128_SOURCES = Rgbsv_example_binary128.cpp
Rgbsv_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgbsv_example_binary128_LDFLAGS = $(binary128_libdepends)

Rgecon_example_binary128_SOURCES = Rgecon_example_binary128.cpp
Rgecon_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgecon_example_binary128_LDFLAGS = $(binary128_libdepends)

Rgesv_example_binary128_SOURCES = Rgesv_example_binary128.cpp
Rgesv_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgesv_example_binary128_LDFLAGS = $(binary128_libdepends)

Rgesvx_example_binary128_SOURCES = Rgesvx_example_binary128.cpp
Rgesvx_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgesvx_example_binary128_LDFLAGS = $(binary128_libdepends)

Rgetrf_getrs_example_binary128_SOURCES = Rgetrf_getrs_example_binary128.cpp
Rgetrf_getrs_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgetrf_getrs_example_binary128_LDFLAGS = $(binary128_libdepends)

Rgetri_example2_binary128_SOURCES = Rgetri_example2_binary128.cpp
Rgetri_example2_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgetri_example2_binary128_LDFLAGS = $(binary128_libdepends)

Rgetri_example3_binary128_SOURCES = Rgetri_example3_binary128.cpp
Rgetri_example3_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgetri_example3_binary128_LDFLAGS = $(binary128_libdepends)

Rgetri_example_binary128_SOURCES = Rgetri_example_binary128.cpp
Rgetri_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgetri_example_binary128_LDFLAGS = $(binary128_libdepends)

Rgtsv_example_binary128_SOURCES = Rgtsv_example_binary128.cpp
Rgtsv_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgtsv_example_binary128_LDFLAGS = $(binary128_libdepends)

endif

if ENABLE_BINARY80
mplapackexamples_PROGRAMS += Cgesv_example_binary80 Cgetrf_getrs_example_binary80 Cgetri_example_binary80 Rgbsv_example_binary80 Rgecon_example_binary80 Rgesv_example_binary80 Rgesvx_example_binary80 Rgetrf_getrs_example_binary80 Rgetri_example2_binary80 Rgetri_example3_binary80 Rgetri_example_binary80 Rgtsv_example_binary80

binary80_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
binary80_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary80

Cgesv_example_binary80_SOURCES = Cgesv_example_binary80.cpp
Cgesv_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Cgesv_example_binary80_LDFLAGS = $(binary80_libdepends)

Cgetrf_getrs_example_binary80_SOURCES = Cgetrf_getrs_example_binary80.cpp
Cgetrf_getrs_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Cgetrf_getrs_example_binary80_LDFLAGS = $(binary80_libdepends)

Cgetri_example_binary80_SOURCES = Cgetri_example_binary80.cpp
Cgetri_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Cgetri_example_binary80_LDFLAGS = $(binary80_libdepends)

Rgbsv_example_binary80_SOURCES = Rgbsv_example_binary80.cpp
Rgbsv_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgbsv_example_binary80_LDFLAGS = $(binary80_libdepends)

Rgecon_example_binary80_SOURCES = Rgecon_example_binary80.cpp
Rgecon_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgecon_example_binary80_LDFLAGS = $(binary80_libdepends)

Rgesv_example_binary80_SOURCES = Rgesv_example_binary80.cpp
Rgesv_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgesv_example_binary80_LDFLAGS = $(binary80_libdepends)

Rgesvx_example_binary80_SOURCES = Rgesvx_example_binary80.cpp
Rgesvx_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgesvx_example_binary80_LDFLAGS = $(binary80_libdepends)

Rgetrf_getrs_example_binary80_SOURCES = Rgetrf_getrs_example_binary80.cpp
Rgetrf_getrs_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgetrf_getrs_example_binary80_LDFLAGS = $(binary80_libdepends)

Rgetri_example2_binary80_SOURCES = Rgetri_example2_binary80.cpp
Rgetri_example2_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgetri_example2_binary80_LDFLAGS = $(binary80_libdepends)

Rgetri_example3_binary80_SOURCES = Rgetri_example3_binary80.cpp
Rgetri_example3_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgetri_example3_binary80_LDFLAGS = $(binary80_libdepends)

Rgetri_example_binary80_SOURCES = Rgetri_example_binary80.cpp
Rgetri_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgetri_example_binary80_LDFLAGS = $(binary80_libdepends)

Rgtsv_example_binary80_SOURCES = Rgtsv_example_binary80.cpp
Rgtsv_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgtsv_example_binary80_LDFLAGS = $(binary80_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Cgesv_example_double Cgetrf_getrs_example_double Cgetri_example_double Rgbsv_example_double Rgecon_example_double Rgesv_example_double Rgesvx_example_double Rgetrf_getrs_example_double Rgetri_example2_double Rgetri_example3_double Rgetri_example_double Rgtsv_example_double

double_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double

Cgesv_example_double_SOURCES = Cgesv_example_double.cpp
Cgesv_example_double_CXXFLAGS = $(double_cxxflags)
Cgesv_example_double_LDFLAGS = $(double_libdepends)

Cgetrf_getrs_example_double_SOURCES = Cgetrf_getrs_example_double.cpp
Cgetrf_getrs_example_double_CXXFLAGS = $(double_cxxflags)
Cgetrf_getrs_example_double_LDFLAGS = $(double_libdepends)

Cgetri_example_double_SOURCES = Cgetri_example_double.cpp
Cgetri_example_double_CXXFLAGS = $(double_cxxflags)
Cgetri_example_double_LDFLAGS = $(double_libdepends)

Rgbsv_example_double_SOURCES = Rgbsv_example_double.cpp
Rgbsv_example_double_CXXFLAGS = $(double_cxxflags)
Rgbsv_example_double_LDFLAGS = $(double_libdepends)

Rgecon_example_double_SOURCES = Rgecon_example_double.cpp
Rgecon_example_double_CXXFLAGS = $(double_cxxflags)
Rgecon_example_double_LDFLAGS = $(double_libdepends)

Rgesv_example_double_SOURCES = Rgesv_example_double.cpp
Rgesv_example_double_CXXFLAGS = $(double_cxxflags)
Rgesv_example_double_LDFLAGS = $(double_libdepends)

Rgesvx_example_double_SOURCES = Rgesvx_example_double.cpp
Rgesvx_example_double_CXXFLAGS = $(double_cxxflags)
Rgesvx_example_double_LDFLAGS = $(double_libdepends)

Rgetrf_getrs_example_double_SOURCES = Rgetrf_getrs_example_double.cpp
Rgetrf_getrs_example_double_CXXFLAGS = $(double_cxxflags)
Rgetrf_getrs_example_double_LDFLAGS = $(double_libdepends)

Rgetri_example2_double_SOURCES = Rgetri_example2_double.cpp
Rgetri_example2_double_CXXFLAGS = $(double_cxxflags)
Rgetri_example2_double_LDFLAGS = $(double_libdepends)

Rgetri_example3_double_SOURCES = Rgetri_example3_double.cpp
Rgetri_example3_double_CXXFLAGS = $(double_cxxflags)
Rgetri_example3_double_LDFLAGS = $(double_libdepends)

Rgetri_example_double_SOURCES = Rgetri_example_double.cpp
Rgetri_example_double_CXXFLAGS = $(double_cxxflags)
Rgetri_example_double_LDFLAGS = $(double_libdepends)

Rgtsv_example_double_SOURCES = Rgtsv_example_double.cpp
Rgtsv_example_double_CXXFLAGS = $(double_cxxflags)
Rgtsv_example_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Cgesv_example_dd Cgetrf_getrs_example_dd Cgetri_example_dd Rgbsv_example_dd Rgecon_example_dd Rgesv_example_dd Rgesvx_example_dd Rgetrf_getrs_example_dd Rgetri_example2_dd Rgetri_example3_dd Rgetri_example_dd Rgtsv_example_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(QD_LIBDIR) -lqd

Cgesv_example_dd_SOURCES = Cgesv_example_dd.cpp
Cgesv_example_dd_CXXFLAGS = $(dd_cxxflags)
Cgesv_example_dd_LDFLAGS = $(dd_libdepends)

Cgetrf_getrs_example_dd_SOURCES = Cgetrf_getrs_example_dd.cpp
Cgetrf_getrs_example_dd_CXXFLAGS = $(dd_cxxflags)
Cgetrf_getrs_example_dd_LDFLAGS = $(dd_libdepends)

Cgetri_example_dd_SOURCES = Cgetri_example_dd.cpp
Cgetri_example_dd_CXXFLAGS = $(dd_cxxflags)
Cgetri_example_dd_LDFLAGS = $(dd_libdepends)

Rgbsv_example_dd_SOURCES = Rgbsv_example_dd.cpp
Rgbsv_example_dd_CXXFLAGS = $(dd_cxxflags)
Rgbsv_example_dd_LDFLAGS = $(dd_libdepends)

Rgecon_example_dd_SOURCES = Rgecon_example_dd.cpp
Rgecon_example_dd_CXXFLAGS = $(dd_cxxflags)
Rgecon_example_dd_LDFLAGS = $(dd_libdepends)

Rgesv_example_dd_SOURCES = Rgesv_example_dd.cpp
Rgesv_example_dd_CXXFLAGS = $(dd_cxxflags)
Rgesv_example_dd_LDFLAGS = $(dd_libdepends)

Rgesvx_example_dd_SOURCES = Rgesvx_example_dd.cpp
Rgesvx_example_dd_CXXFLAGS = $(dd_cxxflags)
Rgesvx_example_dd_LDFLAGS = $(dd_libdepends)

Rgetrf_getrs_example_dd_SOURCES = Rgetrf_getrs_example_dd.cpp
Rgetrf_getrs_example_dd_CXXFLAGS = $(dd_cxxflags)
Rgetrf_getrs_example_dd_LDFLAGS = $(dd_libdepends)

Rgetri_example2_dd_SOURCES = Rgetri_example2_dd.cpp
Rgetri_example2_dd_CXXFLAGS = $(dd_cxxflags)
Rgetri_example2_dd_LDFLAGS = $(dd_libdepends)

Rgetri_example3_dd_SOURCES = Rgetri_example3_dd.cpp
Rgetri_example3_dd_CXXFLAGS = $(dd_cxxflags)
Rgetri_example3_dd_LDFLAGS = $(dd_libdepends)

Rgetri_example_dd_SOURCES = Rgetri_example_dd.cpp
Rgetri_example_dd_CXXFLAGS = $(dd_cxxflags)
Rgetri_example_dd_LDFLAGS = $(dd_libdepends)

Rgtsv_example_dd_SOURCES = Rgtsv_example_dd.cpp
Rgtsv_example_dd_CXXFLAGS = $(dd_cxxflags)
Rgtsv_example_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Cgesv_example_qd Cgetrf_getrs_example_qd Cgetri_example_qd Rgbsv_example_qd Rgecon_example_qd Rgesv_example_qd Rgesvx_example_qd Rgetrf_getrs_example_qd Rgetri_example2_qd Rgetri_example3_qd Rgetri_example_qd Rgtsv_example_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(QD_LIBDIR) -lqd

Cgesv_example_qd_SOURCES = Cgesv_example_qd.cpp
Cgesv_example_qd_CXXFLAGS = $(qd_cxxflags)
Cgesv_example_qd_LDFLAGS = $(qd_libdepends)

Cgetrf_getrs_example_qd_SOURCES = Cgetrf_getrs_example_qd.cpp
Cgetrf_getrs_example_qd_CXXFLAGS = $(qd_cxxflags)
Cgetrf_getrs_example_qd_LDFLAGS = $(qd_libdepends)

Cgetri_example_qd_SOURCES = Cgetri_example_qd.cpp
Cgetri_example_qd_CXXFLAGS = $(qd_cxxflags)
Cgetri_example_qd_LDFLAGS = $(qd_libdepends)

Rgbsv_example_qd_SOURCES = Rgbsv_example_qd.cpp
Rgbsv_example_qd_CXXFLAGS = $(qd_cxxflags)
Rgbsv_example_qd_LDFLAGS = $(qd_libdepends)

Rgecon_example_qd_SOURCES = Rgecon_example_qd.cpp
Rgecon_example_qd_CXXFLAGS = $(qd_cxxflags)
Rgecon_example_qd_LDFLAGS = $(qd_libdepends)

Rgesv_example_qd_SOURCES = Rgesv_example_qd.cpp
Rgesv_example_qd_CXXFLAGS = $(qd_cxxflags)
Rgesv_example_qd_LDFLAGS = $(qd_libdepends)

Rgesvx_example_qd_SOURCES = Rgesvx_example_qd.cpp
Rgesvx_example_qd_CXXFLAGS = $(qd_cxxflags)
Rgesvx_example_qd_LDFLAGS = $(qd_libdepends)

Rgetrf_getrs_example_qd_SOURCES = Rgetrf_getrs_example_qd.cpp
Rgetrf_getrs_example_qd_CXXFLAGS = $(qd_cxxflags)
Rgetrf_getrs_example_qd_LDFLAGS = $(qd_libdepends)

Rgetri_example2_qd_SOURCES = Rgetri_example2_qd.cpp
Rgetri_example2_qd_CXXFLAGS = $(qd_cxxflags)
Rgetri_example2_qd_LDFLAGS = $(qd_libdepends)

Rgetri_example3_qd_SOURCES = Rgetri_example3_qd.cpp
Rgetri_example3_qd_CXXFLAGS = $(qd_cxxflags)
Rgetri_example3_qd_LDFLAGS = $(qd_libdepends)

Rgetri_example_qd_SOURCES = Rgetri_example_qd.cpp
Rgetri_example_qd_CXXFLAGS = $(qd_cxxflags)
Rgetri_example_qd_LDFLAGS = $(qd_libdepends)

Rgtsv_example_qd_SOURCES = Rgtsv_example_qd.cpp
Rgtsv_example_qd_CXXFLAGS = $(qd_cxxflags)
Rgtsv_example_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/examples/mplapack/00_GeneralLinearEquations

mplapackexamples_DATA = Cgesv_example_mpfr.cpp Cgesv_example_gmp.cpp Cgesv_example_binary128.cpp Cgesv_example_binary80.cpp Cgesv_example_double.cpp Cgesv_example_dd.cpp Cgesv_example_qd.cpp Cgetrf_getrs_example_mpfr.cpp Cgetrf_getrs_example_gmp.cpp Cgetrf_getrs_example_binary128.cpp Cgetrf_getrs_example_binary80.cpp Cgetrf_getrs_example_double.cpp Cgetrf_getrs_example_dd.cpp Cgetrf_getrs_example_qd.cpp Cgetri_example_mpfr.cpp Cgetri_example_gmp.cpp Cgetri_example_binary128.cpp Cgetri_example_binary80.cpp Cgetri_example_double.cpp Cgetri_example_dd.cpp Cgetri_example_qd.cpp Rgbsv_example_mpfr.cpp Rgbsv_example_gmp.cpp Rgbsv_example_binary128.cpp Rgbsv_example_binary80.cpp Rgbsv_example_double.cpp Rgbsv_example_dd.cpp Rgbsv_example_qd.cpp Rgecon_example_mpfr.cpp Rgecon_example_gmp.cpp Rgecon_example_binary128.cpp Rgecon_example_binary80.cpp Rgecon_example_double.cpp Rgecon_example_dd.cpp Rgecon_example_qd.cpp Rgesv_example_mpfr.cpp Rgesv_example_gmp.cpp Rgesv_example_binary128.cpp Rgesv_example_binary80.cpp Rgesv_example_double.cpp Rgesv_example_dd.cpp Rgesv_example_qd.cpp Rgesvx_example_mpfr.cpp Rgesvx_example_gmp.cpp Rgesvx_example_binary128.cpp Rgesvx_example_binary80.cpp Rgesvx_example_double.cpp Rgesvx_example_dd.cpp Rgesvx_example_qd.cpp Rgetrf_getrs_example_mpfr.cpp Rgetrf_getrs_example_gmp.cpp Rgetrf_getrs_example_binary128.cpp Rgetrf_getrs_example_binary80.cpp Rgetrf_getrs_example_double.cpp Rgetrf_getrs_example_dd.cpp Rgetrf_getrs_example_qd.cpp Rgetri_example2_mpfr.cpp Rgetri_example2_gmp.cpp Rgetri_example2_binary128.cpp Rgetri_example2_binary80.cpp Rgetri_example2_double.cpp Rgetri_example2_dd.cpp Rgetri_example2_qd.cpp Rgetri_example3_mpfr.cpp Rgetri_example3_gmp.cpp Rgetri_example3_binary128.cpp Rgetri_example3_binary80.cpp Rgetri_example3_double.cpp Rgetri_example3_dd.cpp Rgetri_example3_qd.cpp Rgetri_example_mpfr.cpp Rgetri_example_gmp.cpp Rgetri_example_binary128.cpp Rgetri_example_binary80.cpp Rgetri_example_double.cpp Rgetri_example_dd.cpp Rgetri_example_qd.cpp Rgtsv_example_mpfr.cpp Rgtsv_example_gmp.cpp Rgtsv_example_binary128.cpp Rgtsv_example_binary80.cpp Rgtsv_example_double.cpp Rgtsv_example_dd.cpp Rgtsv_example_qd.cpp \
 Makefile.freebsd Makefile.linux Makefile.linux.inteloneAPI Makefile.macos Makefile.mingw

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapackexamplesdir) $(prefix)
endif

EXTRA_DIST = \
	Makefile.freebsd.in \
	Makefile.linux.in \
	Makefile.linux.inteloneAPI.in \
	Makefile.macos.in \
	Makefile.mingw.in

GENERATED_MAKEFILES = \
	Makefile.freebsd \
	Makefile.linux \
	Makefile.linux.inteloneAPI \
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
        -e 's|%%LIBQUADMATH%%| -lquadmath|g' \
        -e 's|%%QUADMATH_CPPFLAGS%%|$(QUADMATH_CPPFLAGS)|g' \
        -e 's|%%QUADMATH_LDFLAGS%%|$(QUADMATH_LDFLAGS)|g'
else
MPLAPACK_SED_SUBST = \
        -e 's|%%MPLAPACKDIR%%|$(prefix)|g' \
        -e 's|%%LIBQUADMATH%%||g' \
        -e 's|%%QUADMATH_CPPFLAGS%%||g' \
        -e 's|%%QUADMATH_LDFLAGS%%||g'
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

Makefile.linux.inteloneAPI: Makefile.linux.inteloneAPI.in
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
