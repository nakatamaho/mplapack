mpblasexamples_PROGRAMS =

if ENABLE_MPFR
mpblasexamples_PROGRAMS += Cgemm_mpfr Rgemm_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Cgemm_mpfr_SOURCES = Cgemm_mpfr.cpp
Cgemm_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgemm_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgemm_mpfr_SOURCES = Rgemm_mpfr.cpp
Rgemm_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgemm_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mpblasexamples_PROGRAMS += Cgemm_gmp Rgemm_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_gmp -L$(GMP_LIBDIR) -lgmp

Cgemm_gmp_SOURCES = Cgemm_gmp.cpp
Cgemm_gmp_CXXFLAGS = $(gmp_cxxflags)
Cgemm_gmp_LDFLAGS = $(gmp_libdepends)

Rgemm_gmp_SOURCES = Rgemm_gmp.cpp
Rgemm_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgemm_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE__FLOAT128
mpblasexamples_PROGRAMS += Cgemm__Float128 Rgemm__Float128

_Float128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if WANT_QUADMATH
_Float128_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas__Float128 -lquadmath
else
_Float128_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas__Float128
endif

Cgemm__Float128_SOURCES = Cgemm__Float128.cpp
Cgemm__Float128_CXXFLAGS = $(_Float128_cxxflags)
Cgemm__Float128_LDFLAGS = $(_Float128_libdepends)

Rgemm__Float128_SOURCES = Rgemm__Float128.cpp
Rgemm__Float128_CXXFLAGS = $(_Float128_cxxflags)
Rgemm__Float128_LDFLAGS = $(_Float128_libdepends)

endif

if ENABLE__FLOAT64X
mpblasexamples_PROGRAMS += Cgemm__Float64x Rgemm__Float64x

_Float64x_cxxflags = $(OPENMP_CXXFLAGS)
_Float64x_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas__Float64x

Cgemm__Float64x_SOURCES = Cgemm__Float64x.cpp
Cgemm__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Cgemm__Float64x_LDFLAGS = $(_Float64x_libdepends)

Rgemm__Float64x_SOURCES = Rgemm__Float64x.cpp
Rgemm__Float64x_CXXFLAGS = $(_Float64x_cxxflags)
Rgemm__Float64x_LDFLAGS = $(_Float64x_libdepends)

endif

if ENABLE_DOUBLE
mpblasexamples_PROGRAMS += Cgemm_double Rgemm_double

double_cxxflags = $(OPENMP_CXXFLAGS)
double_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_double

Cgemm_double_SOURCES = Cgemm_double.cpp
Cgemm_double_CXXFLAGS = $(double_cxxflags)
Cgemm_double_LDFLAGS = $(double_libdepends)

Rgemm_double_SOURCES = Rgemm_double.cpp
Rgemm_double_CXXFLAGS = $(double_cxxflags)
Rgemm_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mpblasexamples_PROGRAMS += Cgemm_dd Rgemm_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_dd -L$(QD_LIBDIR) -lqd

Cgemm_dd_SOURCES = Cgemm_dd.cpp
Cgemm_dd_CXXFLAGS = $(dd_cxxflags)
Cgemm_dd_LDFLAGS = $(dd_libdepends)

Rgemm_dd_SOURCES = Rgemm_dd.cpp
Rgemm_dd_CXXFLAGS = $(dd_cxxflags)
Rgemm_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mpblasexamples_PROGRAMS += Cgemm_qd Rgemm_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_qd -L$(QD_LIBDIR) -lqd

Cgemm_qd_SOURCES = Cgemm_qd.cpp
Cgemm_qd_CXXFLAGS = $(qd_cxxflags)
Cgemm_qd_LDFLAGS = $(qd_libdepends)

Rgemm_qd_SOURCES = Rgemm_qd.cpp
Rgemm_qd_CXXFLAGS = $(qd_cxxflags)
Rgemm_qd_LDFLAGS = $(qd_libdepends)

endif

mpblasexamplesdir=$(prefix)/share/mpblas/examples/mpblas

mpblasexamples_DATA = Cgemm_mpfr.cpp Cgemm_gmp.cpp Cgemm__Float128.cpp Cgemm__Float64x.cpp Cgemm_double.cpp Cgemm_dd.cpp Cgemm_qd.cpp Rgemm_mpfr.cpp Rgemm_gmp.cpp Rgemm__Float128.cpp Rgemm__Float64x.cpp Rgemm_double.cpp Rgemm_dd.cpp Rgemm_qd.cpp

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mpblasexamplesdir) $(prefix)
endif
