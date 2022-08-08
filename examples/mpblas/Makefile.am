mpblasexamples_PROGRAMS =

if ENABLE_MPFR
mpblasexamples_PROGRAMS += Cgemm_mpfr Rgemm_mpfr Cgemm_mpfr_opt Rgemm_mpfr_opt

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp
mpfr_opt_libdepends = -L$(top_builddir)/mpblas/optimized/mpfr -lmpblas_mpfr_opt -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Cgemm_mpfr_SOURCES = Cgemm_mpfr.cpp
Cgemm_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgemm_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgemm_mpfr_SOURCES = Rgemm_mpfr.cpp
Rgemm_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgemm_mpfr_LDFLAGS = $(mpfr_libdepends)

Cgemm_mpfr_opt_SOURCES = Cgemm_mpfr.cpp
Cgemm_mpfr_opt_CXXFLAGS = $(mpfr_cxxflags)
Cgemm_mpfr_opt_LDFLAGS = $(mpfr_opt_libdepends)

Rgemm_mpfr_opt_SOURCES = Rgemm_mpfr.cpp
Rgemm_mpfr_opt_CXXFLAGS = $(mpfr_cxxflags)
Rgemm_mpfr_opt_LDFLAGS = $(mpfr_opt_libdepends)

endif

if ENABLE_GMP
mpblasexamples_PROGRAMS += Cgemm_gmp Rgemm_gmp Cgemm_gmp_opt Rgemm_gmp_opt

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_gmp -L$(GMP_LIBDIR) -lgmp
gmp_opt_libdepends = -L$(top_builddir)/mpblas/optimized/gmp -lmpblas_gmp_opt -L$(GMP_LIBDIR) -lgmp

Cgemm_gmp_SOURCES = Cgemm_gmp.cpp
Cgemm_gmp_CXXFLAGS = $(gmp_cxxflags)
Cgemm_gmp_LDFLAGS = $(gmp_libdepends)

Rgemm_gmp_SOURCES = Rgemm_gmp.cpp
Rgemm_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgemm_gmp_LDFLAGS = $(gmp_libdepends)

Cgemm_gmp_opt_SOURCES = Cgemm_gmp.cpp
Cgemm_gmp_opt_CXXFLAGS = $(gmp_cxxflags)
Cgemm_gmp_opt_LDFLAGS = $(gmp_opt_libdepends)

Rgemm_gmp_opt_SOURCES = Rgemm_gmp.cpp
Rgemm_gmp_opt_CXXFLAGS = $(gmp_cxxflags)
Rgemm_gmp_opt_LDFLAGS = $(gmp_opt_libdepends)

endif

if ENABLE__FLOAT128
mpblasexamples_PROGRAMS += Cgemm__Float128 Rgemm__Float128 Cgemm__Float128_opt Rgemm__Float128_opt

_Float128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if WANT_QUADMATH
_Float128_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas__Float128 -lquadmath
_Float128_opt_libdepends = -L$(top_builddir)/mpblas/optimized/_Float128 -lmpblas__Float128_opt -lquadmath
else
_Float128_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas__Float128
_Float128_opt_libdepends = -L$(top_builddir)/mpblas/optimized/_Float128 -lmpblas__Float128_opt
endif

Cgemm__Float128_SOURCES = Cgemm__Float128.cpp
Cgemm__Float128_CXXFLAGS = $(_Float128_cxxflags)
Cgemm__Float128_LDFLAGS = $(_Float128_libdepends)

Rgemm__Float128_SOURCES = Rgemm__Float128.cpp
Rgemm__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgemm__Float128_LDFLAGS = $(_Float128_libdepends)

Cgemm__Float128_opt_SOURCES = Cgemm__Float128.cpp
Cgemm__Float128_opt_CXXFLAGS = $(_Float128_cxxflags)
Cgemm__Float128_opt_LDFLAGS = $(_Float128_opt_libdepends)

Rgemm__Float128_opt_SOURCES = Rgemm__Float128.cpp
Rgemm__Float128_opt_CXXFLAGS = $(_Float128_cxxflags)
Rgemm__Float128_opt_LDFLAGS = $(_Float128_opt_libdepends)

endif

if ENABLE__FLOAT64X
mpblasexamples_PROGRAMS += Cgemm__Float64x Rgemm__Float64x Cgemm__Float64x_opt Rgemm__Float64x_opt

_Float64x_cxxflags = $(OPENMP_CXXFLAGS)
_Float64x_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas__Float64x
_Float64x_opt_libdepends = -L$(top_builddir)/mpblas/optimized/_Float64x -lmpblas__Float64x_opt

Cgemm__Float64x_SOURCES = Cgemm__Float64x.cpp
Cgemm__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Cgemm__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgemm__Float64x_SOURCES = Rgemm__Float64x.cpp
Rgemm__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgemm__Float64x_LDFLAGS = $(_Float64x_libdepends)

Cgemm__Float64x_opt_SOURCES = Cgemm__Float64x.cpp
Cgemm__Float64x_opt_CXXFLAGS = $(_Float64x_cxxflags)
Cgemm__Float64x_opt_LDFLAGS = $(_Float64x_opt_libdepends)

Rgemm__Float64x_opt_SOURCES = Rgemm__Float64x.cpp
Rgemm__Float64x_opt_CXXFLAGS = $(_Float64x_cxxflags)
Rgemm__Float64x_opt_LDFLAGS = $(_Float64x_opt_libdepends)

endif

if ENABLE_DOUBLE
mpblasexamples_PROGRAMS += Cgemm_double Rgemm_double Cgemm_double_opt Rgemm_double_opt

double_cxxflags = $(OPENMP_CXXFLAGS)
double_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_double
double_opt_libdepends = -L$(top_builddir)/mpblas/optimized/double -lmpblas_double_opt

Cgemm_double_SOURCES = Cgemm_double.cpp
Cgemm_double_CXXFLAGS = $(double_cxxflags)
Cgemm_double_LDFLAGS = $(double_libdepends)

Rgemm_double_SOURCES = Rgemm_double.cpp
Rgemm_double_CXXFLAGS = $(double_cxxflags)
Rgemm_double_LDFLAGS = $(double_libdepends)

Cgemm_double_opt_SOURCES = Cgemm_double.cpp
Cgemm_double_opt_CXXFLAGS = $(double_cxxflags)
Cgemm_double_opt_LDFLAGS = $(double_opt_libdepends)

Rgemm_double_opt_SOURCES = Rgemm_double.cpp
Rgemm_double_opt_CXXFLAGS = $(double_cxxflags)
Rgemm_double_opt_LDFLAGS = $(double_opt_libdepends)

endif

if ENABLE_DD
mpblasexamples_PROGRAMS += Cgemm_dd Rgemm_dd Cgemm_dd_opt Rgemm_dd_opt

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_dd -L$(QD_LIBDIR) -lqd
dd_opt_libdepends = -L$(top_builddir)/mpblas/optimized/dd -lmpblas_dd_opt -L$(QD_LIBDIR) -lqd

Cgemm_dd_SOURCES = Cgemm_dd.cpp
Cgemm_dd_CXXFLAGS = $(dd_cxxflags)
Cgemm_dd_LDFLAGS = $(dd_libdepends)

Rgemm_dd_SOURCES = Rgemm_dd.cpp
Rgemm_dd_CXXFLAGS = $(dd_cxxflags)
Rgemm_dd_LDFLAGS = $(dd_libdepends)

Cgemm_dd_opt_SOURCES = Cgemm_dd.cpp
Cgemm_dd_opt_CXXFLAGS = $(dd_cxxflags)
Cgemm_dd_opt_LDFLAGS = $(dd_opt_libdepends)

Rgemm_dd_opt_SOURCES = Rgemm_dd.cpp
Rgemm_dd_opt_CXXFLAGS = $(dd_cxxflags)
Rgemm_dd_opt_LDFLAGS = $(dd_opt_libdepends)

endif

if ENABLE_QD
mpblasexamples_PROGRAMS += Cgemm_qd Rgemm_qd Cgemm_qd_opt Rgemm_qd_opt

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_qd -L$(QD_LIBDIR) -lqd
qd_opt_libdepends = -L$(top_builddir)/mpblas/optimized/qd -lmpblas_qd_opt -L$(QD_LIBDIR) -lqd

Cgemm_qd_SOURCES = Cgemm_qd.cpp
Cgemm_qd_CXXFLAGS = $(qd_cxxflags)
Cgemm_qd_LDFLAGS = $(qd_libdepends)

Rgemm_qd_SOURCES = Rgemm_qd.cpp
Rgemm_qd_CXXFLAGS = $(qd_cxxflags)
Rgemm_qd_LDFLAGS = $(qd_libdepends)

Cgemm_qd_opt_SOURCES = Cgemm_qd.cpp
Cgemm_qd_opt_CXXFLAGS = $(qd_cxxflags)
Cgemm_qd_opt_LDFLAGS = $(qd_opt_libdepends)

Rgemm_qd_opt_SOURCES = Rgemm_qd.cpp
Rgemm_qd_opt_CXXFLAGS = $(qd_cxxflags)
Rgemm_qd_opt_LDFLAGS = $(qd_opt_libdepends)

endif

mpblasexamplesdir=$(prefix)/share/examples/mpblas/

mpblasexamples_DATA = Cgemm_mpfr.cpp Cgemm_gmp.cpp Cgemm__Float128.cpp Cgemm__Float64x.cpp Cgemm_double.cpp Cgemm_dd.cpp Cgemm_qd.cpp Rgemm_mpfr.cpp Rgemm_gmp.cpp Rgemm__Float128.cpp Rgemm__Float64x.cpp Rgemm_double.cpp Rgemm_dd.cpp Rgemm_qd.cpp \
Makefile.freebsd Makefile.linux Makefile.linux_cuda Makefile.macos Makefile.mingw

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mpblasexamplesdir) $(prefix)
endif
