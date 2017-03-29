lib_LTLIBRARIES =

if ENABLE_GMP
lib_LTLIBRARIES += libmblas_gmp_ref.la
endif

if ENABLE_MPFR
lib_LTLIBRARIES += libmblas_mpfr_ref.la
endif

if ENABLE_QD
lib_LTLIBRARIES += libmblas_qd_ref.la
endif

if ENABLE_DD
lib_LTLIBRARIES += libmblas_dd_ref.la
endif

if ENABLE_DOUBLE
lib_LTLIBRARIES += libmblas_double_ref.la
endif

if ENABLE___FLOAT128
lib_LTLIBRARIES += libmblas___float128_ref.la
endif

MBLAS_SOURCES = Caxpy.cpp Ccopy.cpp Cdotc.cpp Cdotu.cpp Cgbmv.cpp \
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
Rset.cpp Mlsame.cpp Mxerbla.cpp mpackinit.cpp

if ENABLE_GMP
libmblas_gmp_ref_la_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/$(GMP_INCLUDEDIR) -D___MPACK_BUILD_WITH_GMP___
libmblas_gmp_ref_la_SOURCES = $(MBLAS_SOURCES)
if IS_MACOSX
libmblas_gmp_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_MPFR
libmblas_mpfr_ref_la_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/$(GMP_INCLUDEDIR) -I$(top_srcdir)/$(MPFR_INCLUDEDIR) -I$(top_srcdir)/$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPACK_BUILD_WITH_MPFR___
libmblas_mpfr_ref_la_SOURCES = $(MBLAS_SOURCES)
if IS_MACOSX
libmblas_mpfr_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_QD
libmblas_qd_ref_la_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/$(QD_INCLUDEDIR) -D___MPACK_BUILD_WITH_QD___
libmblas_qd_ref_la_SOURCES = $(MBLAS_SOURCES)
if IS_MACOSX
libmblas_qd_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DD
libmblas_dd_ref_la_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/$(QD_INCLUDEDIR) -D___MPACK_BUILD_WITH_DD___
libmblas_dd_ref_la_SOURCES = $(MBLAS_SOURCES)
if IS_MACOSX
libmblas_dd_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DOUBLE
libmblas_double_ref_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPACK_BUILD_WITH_DOUBLE___
libmblas_double_ref_la_SOURCES = $(MBLAS_SOURCES)
if IS_MACOSX
libmblas_double_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE___FLOAT128
libmblas___float128_ref_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPACK_BUILD_WITH___FLOAT128___
libmblas___float128_ref_la_SOURCES = $(MBLAS_SOURCES)
if IS_MACOSX
libmblas___float128_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif
