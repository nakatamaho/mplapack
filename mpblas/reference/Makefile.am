lib_LTLIBRARIES =

if ENABLE_GMP
lib_LTLIBRARIES += libmpblas_gmp_ref.la
endif

if ENABLE_MPFR
lib_LTLIBRARIES += libmpblas_mpfr_ref.la
endif

if ENABLE_QD
lib_LTLIBRARIES += libmpblas_qd_ref.la
endif

if ENABLE_DD
lib_LTLIBRARIES += libmpblas_dd_ref.la
endif

if ENABLE_DOUBLE
lib_LTLIBRARIES += libmpblas_double_ref.la
endif

if ENABLE__FLOAT64X
lib_LTLIBRARIES += libmpblas__Float64x_ref.la
endif

if ENABLE__FLOAT128
lib_LTLIBRARIES += libmpblas__Float128_ref.la
endif

MPBLAS_SOURCES = Caxpy.cpp Ccopy.cpp Cdotc.cpp Cdotu.cpp Cgbmv.cpp \
Cgemm.cpp Cgemv.cpp Cgerc.cpp Cgeru.cpp Chbmv.cpp Chemm.cpp Chemv.cpp \
Cher.cpp Cher2.cpp Cher2k.cpp Cherk.cpp Chpmv.cpp Chpr.cpp Chpr2.cpp \
Crotg.cpp Cscal.cpp CRrot.cpp CRscal.cpp Cswap.cpp Csymm.cpp Csyr2k.cpp \
Csyrk.cpp Ctbmv.cpp Ctbsv.cpp Ctpmv.cpp Ctpsv.cpp Ctrmm.cpp Ctrmv.cpp \
Ctrsm.cpp Ctrsv.cpp RCabs1.cpp RCasum.cpp RCnrm2.cpp Rasum.cpp Raxpy.cpp\
Rcopy.cpp Rdot.cpp Rgbmv.cpp Rgemm.cpp Rgemv.cpp Rger.cpp Rnrm2.cpp Rrot.cpp\
Rrotg.cpp Rrotm.cpp Rrotmg.cpp Rsbmv.cpp Rscal.cpp Rspmv.cpp \
Rspr.cpp Rspr2.cpp Rswap.cpp Rsymm.cpp Rsymv.cpp Rsyr.cpp Rsyr2.cpp \
Rsyr2k.cpp Rsyrk.cpp Rtbmv.cpp Rtbsv.cpp Rtpmv.cpp Rtpsv.cpp Rtrmm.cpp \
Rtrmv.cpp Rtrsm.cpp Rtrsv.cpp iCamax.cpp iRamax.cpp \
Mlsame.cpp Mxerbla.cpp mplapackinit.cpp

if ENABLE_GMP
libmpblas_gmp_ref_la_CPPFLAGS = -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___
libmpblas_gmp_ref_la_SOURCES = $(MPBLAS_SOURCES)
if IS_MACOSX
libmpblas_gmp_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_MPFR
libmpblas_mpfr_ref_la_CPPFLAGS = -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___
libmpblas_mpfr_ref_la_SOURCES = $(MPBLAS_SOURCES)
if IS_MACOSX
libmpblas_mpfr_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_QD
libmpblas_qd_ref_la_CPPFLAGS = -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___
libmpblas_qd_ref_la_SOURCES = $(MPBLAS_SOURCES)
if IS_MACOSX
libmpblas_qd_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DD
libmpblas_dd_ref_la_CPPFLAGS = -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___
libmpblas_dd_ref_la_SOURCES = $(MPBLAS_SOURCES)
if IS_MACOSX
libmpblas_dd_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DOUBLE
libmpblas_double_ref_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_DOUBLE___
libmpblas_double_ref_la_SOURCES = $(MPBLAS_SOURCES)
if IS_MACOSX
libmpblas_double_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE__FLOAT64X
libmpblas__Float64x_ref_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT64X___
libmpblas__Float64x_ref_la_SOURCES = $(MPBLAS_SOURCES)
if IS_MACOSX
libmpblas__Float64x_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE__FLOAT128
libmpblas__Float128_ref_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT128___
libmpblas__Float128_ref_la_SOURCES = $(MPBLAS_SOURCES)
if IS_MACOSX
libmpblas__Float128_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif
