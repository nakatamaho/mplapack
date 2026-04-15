lib_LTLIBRARIES =

if ENABLE_GMP
lib_LTLIBRARIES += libmpblas_gmp.la
endif

if ENABLE_MPFR
lib_LTLIBRARIES += libmpblas_mpfr.la
endif

if ENABLE_QD
lib_LTLIBRARIES += libmpblas_qd.la
endif

if ENABLE_DD
lib_LTLIBRARIES += libmpblas_dd.la
endif

if ENABLE_DOUBLE
lib_LTLIBRARIES += libmpblas_double.la
endif

if ENABLE_BINARY80
lib_LTLIBRARIES += libmpblas_binary80.la
endif

if ENABLE_BINARY128
lib_LTLIBRARIES += libmpblas_binary128.la
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

MPLAPACK_LT_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
MPLAPACK_RPATH_LDFLAGS = -Wl,-rpath,$(libdir)

if ENABLE_GMP
libmpblas_gmp_la_CPPFLAGS = -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___
libmpblas_gmp_la_SOURCES = $(MPBLAS_SOURCES)
libmpblas_gmp_la_LDFLAGS = $(MPLAPACK_LT_LDFLAGS) $(MPLAPACK_RPATH_LDFLAGS)
if IS_MACOS
libmpblas_gmp_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmpblas_gmp_la_LDFLAGS += -no-undefined -L$(GMP_LIBDIR) -lgmp
endif
endif

if ENABLE_MPFR
libmpblas_mpfr_la_CPPFLAGS = -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___
libmpblas_mpfr_la_SOURCES = $(MPBLAS_SOURCES)
libmpblas_mpfr_la_LDFLAGS = $(MPLAPACK_LT_LDFLAGS) $(MPLAPACK_RPATH_LDFLAGS)
if IS_MACOS
libmpblas_mpfr_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmpblas_mpfr_la_LDFLAGS += -no-undefined -L$(MPC_LIBDIR) -lmpc -L$(MPFR_LIBDIR) -lmpfr -L$(GMP_LIBDIR) -lgmp
endif
endif

if ENABLE_QD
libmpblas_qd_la_CPPFLAGS = -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___
libmpblas_qd_la_SOURCES = $(MPBLAS_SOURCES)
libmpblas_qd_la_LDFLAGS = $(MPLAPACK_LT_LDFLAGS) $(MPLAPACK_RPATH_LDFLAGS)
if IS_MACOS
libmpblas_qd_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmpblas_qd_la_LDFLAGS += -no-undefined -L$(QD_LIBDIR) -lqd
endif
endif

if ENABLE_DD
libmpblas_dd_la_CPPFLAGS = -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___
libmpblas_dd_la_SOURCES = $(MPBLAS_SOURCES)
libmpblas_dd_la_LDFLAGS = $(MPLAPACK_LT_LDFLAGS) $(MPLAPACK_RPATH_LDFLAGS)
libmpblas_dd_la_CXXFLAGS  = $(DD_CXXFLAGS)
if IS_MACOS
libmpblas_dd_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmpblas_dd_la_LDFLAGS += -no-undefined -L$(QD_LIBDIR) -lqd
endif
endif

if ENABLE_DOUBLE
libmpblas_double_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_DOUBLE___
libmpblas_double_la_SOURCES = $(MPBLAS_SOURCES)
libmpblas_double_la_LDFLAGS = $(MPLAPACK_LT_LDFLAGS)
if IS_MACOS
libmpblas_double_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmpblas_double_la_LDFLAGS += -no-undefined
endif
endif

if ENABLE_BINARY80
libmpblas_binary80_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_BINARY80___
libmpblas_binary80_la_SOURCES = $(MPBLAS_SOURCES)
libmpblas_binary80_la_LDFLAGS = $(MPLAPACK_LT_LDFLAGS)
if IS_MACOS
libmpblas_binary80_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmpblas_binary80_la_LDFLAGS += -no-undefined
endif
endif

if ENABLE_BINARY128
libmpblas_binary128_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_BINARY128___
libmpblas_binary128_la_SOURCES = $(MPBLAS_SOURCES)
libmpblas_binary128_la_LDFLAGS = $(MPLAPACK_LT_LDFLAGS)
if IS_MACOS
libmpblas_binary128_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmpblas_binary128_la_LDFLAGS += -no-undefined -lquadmath
endif
endif

