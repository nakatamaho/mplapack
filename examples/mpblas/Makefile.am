blasexamples_PROGRAMS =

blasexamplesdir=$(prefix)/share/mplapack/examples/mpblas/

blasexamples_DATA = Rgemm__Float128.cpp Rgemm__Float64x.cpp Rgemm_dd.cpp Rgemm_double.cpp Rgemm_gmp.cpp Rgemm_mpfr.cpp Rgemm_qd.cpp Makefile.freebsd Makefile.linux Makefile.macos Makefile.mingw

if ENABLE_GMP
blasexamples_PROGRAMS += Rgemm_gmp
gmpcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) 
gmp_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_gmp -L$(GMP_LIBDIR) -lgmpxx -lgmp
Rgemm_gmp_SOURCES= Rgemm_gmp.cpp
Rgemm_gmp_CXXFLAGS= $(gmpcxxflags)
Rgemm_gmp_LDFLAGS= -Wl,-rpath,$(libdir) $(gmp_libdepends)
blasexamples_DATA +=
endif

if ENABLE_MPFR
blasexamples_PROGRAMS += Rgemm_mpfr
mpfrcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp
Rgemm_mpfr_SOURCES= Rgemm_mpfr.cpp
Rgemm_mpfr_CXXFLAGS= $(mpfrcxxflags)
Rgemm_mpfr_LDFLAGS= -Wl,-rpath,$(libdir) $(mpfr_libdepends)
endif

if ENABLE_QD
blasexamples_PROGRAMS += Rgemm_qd
qdcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) 
qd_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_qd -L$(QD_LIBDIR) -lqd
Rgemm_qd_SOURCES= Rgemm_qd.cpp
Rgemm_qd_CXXFLAGS= $(qdcxxflags)
Rgemm_qd_LDFLAGS= -Wl,-rpath,$(libdir) $(qd_libdepends)
endif

if ENABLE_DD
blasexamples_PROGRAMS += Rgemm_dd
ddcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) 
dd_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_dd -L$(QD_LIBDIR) -lqd
Rgemm_dd_SOURCES= Rgemm_dd.cpp
Rgemm_dd_CXXFLAGS= $(ddcxxflags)
Rgemm_dd_LDFLAGS= -Wl,-rpath,$(libdir) $(dd_libdepends)
endif

if ENABLE_DOUBLE
blasexamples_PROGRAMS +=  Rgemm_double
doublecxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_double
Rgemm_double_SOURCES= Rgemm_double.cpp
Rgemm_double_CXXFLAGS= $(doublecxxflags)
Rgemm_double_LDFLAGS= -Wl,-rpath,$(libdir) $(double_libdepends)
endif

if ENABLE__FLOAT128
blasexamples_PROGRAMS += Rgemm__Float128
_Float128cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if WANT_QUADMATH
_Float128_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas__Float128 -lquadmath
else
_Float128_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas__Float128
endif
Rgemm__Float128_SOURCES= Rgemm__Float128.cpp
Rgemm__Float128_CXXFLAGS= $(_Float128cxxflags)
Rgemm__Float128_LDFLAGS= -Wl,-rpath,$(libdir) $(_Float128_libdepends)
endif

if ENABLE__FLOAT64X
blasexamples_PROGRAMS += Rgemm__Float64x
_Float64xcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
_Float64x_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas__Float64x
Rgemm__Float64x_SOURCES= Rgemm__Float64x.cpp
Rgemm__Float64x_CXXFLAGS= $(_Float64xcxxflags)
Rgemm__Float64x_LDFLAGS= -Wl,-rpath,$(libdir) $(_Float64x_libdepends)
endif

