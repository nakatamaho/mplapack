blasexamples_PROGRAMS =

if ENABLE_GMP
blasexamples_PROGRAMS += Rgemm.gmp
gmpcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) 
gmp_libdepends = -L$(top_builddir)/mblas/optimized/gmp -lmblas_gmp -L$(GMP_LIBDIR) -lgmpxx -lgmp
Rgemm_gmp_SOURCES= Rgemm_gmp.cpp
Rgemm_gmp_CXXFLAGS= $(gmpcxxflags)
Rgemm_gmp_LDFLAGS= $(gmp_libdepends)
endif

if ENABLE_MPFR
blasexamples_PROGRAMS += Rgemm.mpfr
mpfrcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -L$(top_builddir)/mblas/optimized/mpfr -lmblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp
Rgemm_mpfr_SOURCES= Rgemm_mpfr.cpp
Rgemm_mpfr_CXXFLAGS= $(mpfrcxxflags)
Rgemm_mpfr_LDFLAGS= $(mpfr_libdepends)
endif

if ENABLE_QD
blasexamples_PROGRAMS += Rgemm.qd
qdcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) 
qd_libdepends = -L$(top_builddir)/mblas/optimized/qd -lmblas_qd -L$(QD_LIBDIR) -lqd
Rgemm_qd_SOURCES= Rgemm_qd.cpp
Rgemm_qd_CXXFLAGS= $(qdcxxflags)
Rgemm_qd_LDFLAGS= $(qd_libdepends)
endif

if ENABLE_DD
blasexamples_PROGRAMS += Rgemm.dd
ddcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) 
dd_libdepends = -L$(top_builddir)/mblas/optimized/dd -lmblas_dd -L$(QD_LIBDIR) -lqd
Rgemm_dd_SOURCES= Rgemm_dd.cpp
Rgemm_dd_CXXFLAGS= $(ddcxxflags)
Rgemm_dd_LDFLAGS= $(dd_libdepends)
endif

if ENABLE_DOUBLE
blasexamples_PROGRAMS +=  Rgemm.double
doublecxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -L$(top_builddir)/mblas/optimized/double -lmblas_double
Rgemm_double_SOURCES= Rgemm_double.cpp
Rgemm_double_CXXFLAGS= $(doublecxxflags)
Rgemm_double_LDFLAGS= $(double_libdepends)
endif

if ENABLE___FLOAT128
blasexamples_PROGRAMS += Rgemm.__float128
__float128cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
__float128_libdepends = -L$(top_builddir)/mblas/optimized/__float128 -lmblas___float128 -lquadmath
Rgemm___float128_SOURCES= Rgemm___float128.cpp
Rgemm___float128_CXXFLAGS= $(__float128cxxflags)
Rgemm___float128_LDFLAGS= $(__float128_libdepends)
endif

blasexamplesdir=$(prefix)/share/mpack/examples
