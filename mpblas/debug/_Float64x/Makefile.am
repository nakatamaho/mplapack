#CTBSV, CTPSV might have bugs or intrinsically inaccurate
#Rnrm2 should be checked more
check_PROGRAMS = $(mpblas__Float64x_test_PROGRAMS) 
TESTS = $(check_PROGRAMS)
mpblas__Float64x_test_PROGRAMS = \
arithmetic.debug__Float64x \
complex.debug__Float64x \
iCamax.debug._Float64x_ref iCamax.debug._Float64x   \
iRamax.debug._Float64x_ref iRamax.debug._Float64x   \
Mlsame.debug._Float64x_ref Mlsame.debug._Float64x   \
Mxerbla.debug._Float64x_ref Mxerbla.debug._Float64x \
CRrot.debug._Float64x_ref  CRrot.debug._Float64x    \
CRscal.debug._Float64x_ref CRscal.debug._Float64x   \
Caxpy.debug._Float64x_ref  Caxpy.debug._Float64x    \
Ccopy.debug._Float64x_ref  Ccopy.debug._Float64x    \
Cdotc.debug._Float64x_ref  Cdotc.debug._Float64x	\
Cdotu.debug._Float64x_ref  Cdotu.debug._Float64x	\
Cgbmv.debug._Float64x_ref  Cgbmv.debug._Float64x	\
Cgemm.debug._Float64x_ref  Cgemm.debug._Float64x	\
Cgemv.debug._Float64x_ref  Cgemv.debug._Float64x	\
Cgerc.debug._Float64x_ref  Cgerc.debug._Float64x	\
Cgeru.debug._Float64x_ref  Cgeru.debug._Float64x    \
Chbmv.debug._Float64x_ref  Chbmv.debug._Float64x	\
Chemm.debug._Float64x_ref  Chemm.debug._Float64x	\
Chemv.debug._Float64x_ref  Chemv.debug._Float64x	\
Cher.debug._Float64x_ref   Cher.debug._Float64x 	\
Cher2.debug._Float64x_ref  Cher2.debug._Float64x	\
Cher2k.debug._Float64x_ref Cher2k.debug._Float64x	\
Cherk.debug._Float64x_ref  Cherk.debug._Float64x	\
Chpmv.debug._Float64x_ref  Chpmv.debug._Float64x	\
Chpr.debug._Float64x_ref   Chpr.debug._Float64x	\
Chpr2.debug._Float64x_ref  Chpr2.debug._Float64x	\
Crotg.debug._Float64x_ref  Crotg.debug._Float64x	\
Cscal.debug._Float64x_ref  Cscal.debug._Float64x	\
Cswap.debug._Float64x_ref  Cswap.debug._Float64x	\
Csymm.debug._Float64x_ref  Csymm.debug._Float64x    \
Csyr2k.debug._Float64x_ref Csyr2k.debug._Float64x	\
Csyrk.debug._Float64x_ref  Csyrk.debug._Float64x	\
Ctbmv.debug._Float64x_ref  Ctbmv.debug._Float64x	\
Ctbsv.debug._Float64x_ref  Ctbsv.debug._Float64x	\
Ctpmv.debug._Float64x_ref  Ctpmv.debug._Float64x	\
Ctpsv.debug._Float64x_ref  Ctpsv.debug._Float64x	\
Ctrmm.debug._Float64x_ref  Ctrmm.debug._Float64x	\
Ctrmv.debug._Float64x_ref  Ctrmv.debug._Float64x	\
Ctrsm.debug._Float64x_ref  Ctrsm.debug._Float64x	\
Ctrsv.debug._Float64x_ref  Ctrsv.debug._Float64x	\
RCabs1.debug._Float64x_ref RCabs1.debug._Float64x	\
RCasum.debug._Float64x_ref RCasum.debug._Float64x	\
RCnrm2.debug._Float64x_ref RCnrm2.debug._Float64x	\
Rasum.debug._Float64x_ref  Rasum.debug._Float64x    \
Raxpy.debug._Float64x_ref  Raxpy.debug._Float64x	\
Rcopy.debug._Float64x_ref  Rcopy.debug._Float64x	\
Rdot.debug._Float64x_ref   Rdot.debug._Float64x	\
Rgbmv.debug._Float64x_ref  Rgbmv.debug._Float64x	\
Rgemm.debug._Float64x_ref  Rgemm.debug._Float64x	\
Rgemv.debug._Float64x_ref  Rgemv.debug._Float64x	\
Rger.debug._Float64x_ref   Rger.debug._Float64x	\
Rnrm2.debug._Float64x_ref  Rnrm2.debug._Float64x	\
Rrot.debug._Float64x_ref   Rrot.debug._Float64x	\
Rrotg.debug._Float64x_ref  Rrotg.debug._Float64x	\
Rrotm.debug._Float64x_ref  Rrotm.debug._Float64x	\
Rsbmv.debug._Float64x_ref  Rsbmv.debug._Float64x	\
Rscal.debug._Float64x_ref  Rscal.debug._Float64x	\
Rspmv.debug._Float64x_ref  Rspmv.debug._Float64x    \
Rspr.debug._Float64x_ref   Rspr.debug._Float64x	\
Rspr2.debug._Float64x_ref  Rspr2.debug._Float64x	\
Rswap.debug._Float64x_ref  Rswap.debug._Float64x	\
Rsymm.debug._Float64x_ref  Rsymm.debug._Float64x	\
Rsymv.debug._Float64x_ref  Rsymv.debug._Float64x	\
Rsyr.debug._Float64x_ref   Rsyr.debug._Float64x	\
Rsyr2.debug._Float64x_ref  Rsyr2.debug._Float64x	\
Rsyr2k.debug._Float64x_ref Rsyr2k.debug._Float64x	\
Rsyrk.debug._Float64x_ref  Rsyrk.debug._Float64x	\
Rtbmv.debug._Float64x_ref  Rtbmv.debug._Float64x	\
Rtbsv.debug._Float64x_ref  Rtbsv.debug._Float64x	\
Rtpmv.debug._Float64x_ref  Rtpmv.debug._Float64x	\
Rtpsv.debug._Float64x_ref  Rtpsv.debug._Float64x	\
Rtrmm.debug._Float64x_ref  Rtrmm.debug._Float64x    \
Rtrmv.debug._Float64x_ref  Rtrmv.debug._Float64x	\
Rtrsm.debug._Float64x_ref  Rtrsm.debug._Float64x	\
Rtrsv.debug._Float64x_ref  Rtrsv.debug._Float64x

mpblas__Float64x_testdir = $(prefix)/share/mplapack/testing/_Float64x
ref_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas__Float64x_ref -lmpblas_mpfr_ref
opt_libdepends = -L$(top_builddir)/mpblas/optimized/_Float64x -lmpblas__Float64x -L$(top_builddir)/mpblas/reference -lmpblas_mpfr_ref
mplibs=-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp -lquadmath
refblas=-L$(top_builddir)/external/i/LAPACK/lib -lblas $(FCLIBS)

debug_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH__FLOAT64X___

if !IS_MACOSX
debug_srcdepends = ../common/xerbla.cpp ../common/mplapack.debug.cpp ../common/Mxerbla.override.cpp
endif

if IS_MACOSX
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT64X___
libmxerbla_override_la_SOURCES =  ../common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS= -Wl,-flat_namespace,-undefined,dynamic_lookup
debug_srcdepends = ../common/xerbla.cpp ../common/mplapack.debug.cpp
mplibs+=-Wl,-flat_namespace,-undefined,dynamic_lookup
TESTS_ENVIRONMENT=DYLD_INSERT_LIBRARIES=.libs/libmxerbla_override.0.dylib
endif

##
arithmetic_debug__Float64x_SOURCES=../common/arithmetic.debug.cpp $(debug_srcdepends)
arithmetic_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
arithmetic_debug__Float64x_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

complex_debug__Float64x_SOURCES=../common/complex.debug.cpp $(debug_srcdepends)
complex_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
complex_debug__Float64x_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_debug__Float64x_ref_SOURCES=../common/iCamax.debug.cpp $(debug_srcdepends)
iCamax_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
iCamax_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_debug__Float64x_SOURCES=../common/iCamax.debug.cpp $(debug_srcdepends)
iCamax_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
iCamax_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iRamax_debug__Float64x_ref_SOURCES=../common/iRamax.debug.cpp $(debug_srcdepends)
iRamax_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
iRamax_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iRamax_debug__Float64x_SOURCES=../common/iRamax.debug.cpp $(debug_srcdepends)
iRamax_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
iRamax_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mlsame_debug__Float64x_ref_SOURCES=../common/Mlsame.debug.cpp $(debug_srcdepends)
Mlsame_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Mlsame_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mlsame_debug__Float64x_SOURCES=../common/Mlsame.debug.cpp $(debug_srcdepends)
Mlsame_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Mlsame_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mxerbla_debug__Float64x_ref_SOURCES=../common/Mxerbla.debug.cpp $(debug_srcdepends)
Mxerbla_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Mxerbla_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mxerbla_debug__Float64x_SOURCES=../common/Mxerbla.debug.cpp $(debug_srcdepends)
Mxerbla_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Mxerbla_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRrot_debug__Float64x_ref_SOURCES=../common/CRrot.debug.cpp $(debug_srcdepends)
CRrot_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
CRrot_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRrot_debug__Float64x_SOURCES=../common/CRrot.debug.cpp $(debug_srcdepends)
CRrot_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
CRrot_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRscal_debug__Float64x_ref_SOURCES=../common/CRscal.debug.cpp $(debug_srcdepends)
CRscal_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
CRscal_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRscal_debug__Float64x_SOURCES=../common/CRscal.debug.cpp $(debug_srcdepends)
CRscal_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
CRscal_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_debug__Float64x_SOURCES=../common/Caxpy.debug.cpp $(debug_srcdepends)
Caxpy_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Caxpy_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_debug__Float64x_ref_SOURCES=../common/Caxpy.debug.cpp $(debug_srcdepends)
Caxpy_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Caxpy_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ccopy_debug__Float64x_SOURCES=../common/Ccopy.debug.cpp $(debug_srcdepends)
Ccopy_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Ccopy_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ccopy_debug__Float64x_ref_SOURCES=../common/Ccopy.debug.cpp $(debug_srcdepends)
Ccopy_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Ccopy_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotc_debug__Float64x_SOURCES=../common/Cdotc.debug.cpp $(debug_srcdepends)
Cdotc_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Cdotc_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotc_debug__Float64x_ref_SOURCES=../common/Cdotc.debug.cpp $(debug_srcdepends)
Cdotc_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Cdotc_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotu_debug__Float64x_SOURCES=../common/Cdotu.debug.cpp $(debug_srcdepends)
Cdotu_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Cdotu_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotu_debug__Float64x_ref_SOURCES=../common/Cdotu.debug.cpp $(debug_srcdepends)
Cdotu_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Cdotu_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgbmv_debug__Float64x_SOURCES=../common/Cgbmv.debug.cpp $(debug_srcdepends)
Cgbmv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Cgbmv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgbmv_debug__Float64x_ref_SOURCES=../common/Cgbmv.debug.cpp $(debug_srcdepends)
Cgbmv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Cgbmv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemm_debug__Float64x_SOURCES=../common/Cgemm.debug.cpp $(debug_srcdepends)
Cgemm_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Cgemm_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemm_debug__Float64x_ref_SOURCES=../common/Cgemm.debug.cpp $(debug_srcdepends)
Cgemm_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Cgemm_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemv_debug__Float64x_SOURCES=../common/Cgemv.debug.cpp $(debug_srcdepends)
Cgemv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Cgemv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemv_debug__Float64x_ref_SOURCES=../common/Cgemv.debug.cpp $(debug_srcdepends)
Cgemv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Cgemv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgerc_debug__Float64x_SOURCES=../common/Cgerc.debug.cpp $(debug_srcdepends)
Cgerc_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Cgerc_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgerc_debug__Float64x_ref_SOURCES=../common/Cgerc.debug.cpp $(debug_srcdepends)
Cgerc_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Cgerc_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgeru_debug__Float64x_SOURCES=../common/Cgeru.debug.cpp $(debug_srcdepends)
Cgeru_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Cgeru_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgeru_debug__Float64x_ref_SOURCES=../common/Cgeru.debug.cpp $(debug_srcdepends)
Cgeru_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Cgeru_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chbmv_debug__Float64x_SOURCES=../common/Chbmv.debug.cpp $(debug_srcdepends)
Chbmv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Chbmv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chbmv_debug__Float64x_ref_SOURCES=../common/Chbmv.debug.cpp $(debug_srcdepends)
Chbmv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Chbmv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemm_debug__Float64x_SOURCES=../common/Chemm.debug.cpp $(debug_srcdepends)
Chemm_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Chemm_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemm_debug__Float64x_ref_SOURCES=../common/Chemm.debug.cpp $(debug_srcdepends)
Chemm_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Chemm_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemv_debug__Float64x_SOURCES=../common/Chemv.debug.cpp $(debug_srcdepends)
Chemv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Chemv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemv_debug__Float64x_ref_SOURCES=../common/Chemv.debug.cpp $(debug_srcdepends)
Chemv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Chemv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher_debug__Float64x_SOURCES=../common/Cher.debug.cpp $(debug_srcdepends)
Cher_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Cher_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher_debug__Float64x_ref_SOURCES=../common/Cher.debug.cpp $(debug_srcdepends)
Cher_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Cher_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2_debug__Float64x_SOURCES=../common/Cher2.debug.cpp $(debug_srcdepends)
Cher2_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Cher2_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2_debug__Float64x_ref_SOURCES=../common/Cher2.debug.cpp $(debug_srcdepends)
Cher2_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Cher2_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2k_debug__Float64x_SOURCES=../common/Cher2k.debug.cpp $(debug_srcdepends)
Cher2k_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Cher2k_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2k_debug__Float64x_ref_SOURCES=../common/Cher2k.debug.cpp $(debug_srcdepends)
Cher2k_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Cher2k_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cherk_debug__Float64x_SOURCES=../common/Cherk.debug.cpp $(debug_srcdepends)
Cherk_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Cherk_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cherk_debug__Float64x_ref_SOURCES=../common/Cherk.debug.cpp $(debug_srcdepends)
Cherk_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Cherk_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpmv_debug__Float64x_SOURCES=../common/Chpmv.debug.cpp $(debug_srcdepends)
Chpmv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Chpmv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpmv_debug__Float64x_ref_SOURCES=../common/Chpmv.debug.cpp $(debug_srcdepends)
Chpmv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Chpmv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr_debug__Float64x_SOURCES=../common/Chpr.debug.cpp $(debug_srcdepends)
Chpr_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Chpr_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr_debug__Float64x_ref_SOURCES=../common/Chpr.debug.cpp $(debug_srcdepends)
Chpr_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Chpr_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr2_debug__Float64x_SOURCES=../common/Chpr2.debug.cpp $(debug_srcdepends)
Chpr2_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Chpr2_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr2_debug__Float64x_ref_SOURCES=../common/Chpr2.debug.cpp $(debug_srcdepends)
Chpr2_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Chpr2_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Crotg_debug__Float64x_SOURCES=../common/Crotg.debug.cpp $(debug_srcdepends)
Crotg_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Crotg_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Crotg_debug__Float64x_ref_SOURCES=../common/Crotg.debug.cpp $(debug_srcdepends)
Crotg_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Crotg_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cscal_debug__Float64x_SOURCES=../common/Cscal.debug.cpp $(debug_srcdepends)
Cscal_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Cscal_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cscal_debug__Float64x_ref_SOURCES=../common/Cscal.debug.cpp $(debug_srcdepends)
Cscal_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Cscal_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cswap_debug__Float64x_SOURCES=../common/Cswap.debug.cpp $(debug_srcdepends)
Cswap_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Cswap_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cswap_debug__Float64x_ref_SOURCES=../common/Cswap.debug.cpp $(debug_srcdepends)
Cswap_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Cswap_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csymm_debug__Float64x_SOURCES=../common/Csymm.debug.cpp $(debug_srcdepends)
Csymm_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Csymm_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csymm_debug__Float64x_ref_SOURCES=../common/Csymm.debug.cpp $(debug_srcdepends)
Csymm_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Csymm_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyr2k_debug__Float64x_SOURCES=../common/Csyr2k.debug.cpp $(debug_srcdepends)
Csyr2k_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Csyr2k_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyr2k_debug__Float64x_ref_SOURCES=../common/Csyr2k.debug.cpp $(debug_srcdepends)
Csyr2k_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Csyr2k_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyrk_debug__Float64x_SOURCES=../common/Csyrk.debug.cpp $(debug_srcdepends)
Csyrk_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Csyrk_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyrk_debug__Float64x_ref_SOURCES=../common/Csyrk.debug.cpp $(debug_srcdepends)
Csyrk_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Csyrk_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbmv_debug__Float64x_SOURCES=../common/Ctbmv.debug.cpp $(debug_srcdepends)
Ctbmv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Ctbmv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbmv_debug__Float64x_ref_SOURCES=../common/Ctbmv.debug.cpp $(debug_srcdepends)
Ctbmv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Ctbmv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbsv_debug__Float64x_SOURCES=../common/Ctbsv.debug.cpp $(debug_srcdepends)
Ctbsv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Ctbsv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbsv_debug__Float64x_ref_SOURCES=../common/Ctbsv.debug.cpp $(debug_srcdepends)
Ctbsv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Ctbsv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpmv_debug__Float64x_SOURCES=../common/Ctpmv.debug.cpp $(debug_srcdepends)
Ctpmv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Ctpmv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpmv_debug__Float64x_ref_SOURCES=../common/Ctpmv.debug.cpp $(debug_srcdepends)
Ctpmv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Ctpmv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpsv_debug__Float64x_SOURCES=../common/Ctpsv.debug.cpp $(debug_srcdepends)
Ctpsv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Ctpsv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpsv_debug__Float64x_ref_SOURCES=../common/Ctpsv.debug.cpp $(debug_srcdepends)
Ctpsv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Ctpsv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmm_debug__Float64x_SOURCES=../common/Ctrmm.debug.cpp $(debug_srcdepends)
Ctrmm_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Ctrmm_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmm_debug__Float64x_ref_SOURCES=../common/Ctrmm.debug.cpp $(debug_srcdepends)
Ctrmm_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Ctrmm_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmv_debug__Float64x_SOURCES=../common/Ctrmv.debug.cpp $(debug_srcdepends)
Ctrmv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Ctrmv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmv_debug__Float64x_ref_SOURCES=../common/Ctrmv.debug.cpp $(debug_srcdepends)
Ctrmv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Ctrmv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsm_debug__Float64x_SOURCES=../common/Ctrsm.debug.cpp $(debug_srcdepends)
Ctrsm_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Ctrsm_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsm_debug__Float64x_ref_SOURCES=../common/Ctrsm.debug.cpp $(debug_srcdepends)
Ctrsm_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Ctrsm_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsv_debug__Float64x_SOURCES=../common/Ctrsv.debug.cpp $(debug_srcdepends)
Ctrsv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Ctrsv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsv_debug__Float64x_ref_SOURCES=../common/Ctrsv.debug.cpp $(debug_srcdepends)
Ctrsv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Ctrsv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_debug__Float64x_ref_SOURCES=../common/RCabs1.debug.cpp $(debug_srcdepends)
RCabs1_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
RCabs1_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_debug__Float64x_SOURCES=../common/RCabs1.debug.cpp $(debug_srcdepends)
RCabs1_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
RCabs1_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCasum_debug__Float64x_ref_SOURCES=../common/RCasum.debug.cpp $(debug_srcdepends)
RCasum_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
RCasum_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCasum_debug__Float64x_SOURCES=../common/RCasum.debug.cpp $(debug_srcdepends)
RCasum_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
RCasum_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCnrm2_debug__Float64x_ref_SOURCES=../common/RCnrm2.debug.cpp $(debug_srcdepends)
RCnrm2_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
RCnrm2_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCnrm2_debug__Float64x_SOURCES=../common/RCnrm2.debug.cpp $(debug_srcdepends)
RCnrm2_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
RCnrm2_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rasum_debug__Float64x_ref_SOURCES=../common/Rasum.debug.cpp $(debug_srcdepends)
Rasum_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rasum_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rasum_debug__Float64x_SOURCES=../common/Rasum.debug.cpp $(debug_srcdepends)
Rasum_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rasum_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Raxpy_debug__Float64x_ref_SOURCES=../common/Raxpy.debug.cpp $(debug_srcdepends)
Raxpy_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Raxpy_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Raxpy_debug__Float64x_SOURCES=../common/Raxpy.debug.cpp $(debug_srcdepends)
Raxpy_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Raxpy_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rcopy_debug__Float64x_ref_SOURCES=../common/Rcopy.debug.cpp $(debug_srcdepends)
Rcopy_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rcopy_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rcopy_debug__Float64x_SOURCES=../common/Rcopy.debug.cpp $(debug_srcdepends)
Rcopy_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rcopy_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rdot_debug__Float64x_ref_SOURCES=../common/Rdot.debug.cpp $(debug_srcdepends)
Rdot_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rdot_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rdot_debug__Float64x_SOURCES=../common/Rdot.debug.cpp $(debug_srcdepends)
Rdot_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rdot_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgbmv_debug__Float64x_ref_SOURCES=../common/Rgbmv.debug.cpp $(debug_srcdepends)
Rgbmv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rgbmv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgbmv_debug__Float64x_SOURCES=../common/Rgbmv.debug.cpp $(debug_srcdepends)
Rgbmv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rgbmv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemm_debug__Float64x_ref_SOURCES=../common/Rgemm.debug.cpp $(debug_srcdepends)
Rgemm_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rgemm_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemm_debug__Float64x_SOURCES=../common/Rgemm.debug.cpp $(debug_srcdepends)
Rgemm_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rgemm_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemv_debug__Float64x_ref_SOURCES=../common/Rgemv.debug.cpp $(debug_srcdepends)
Rgemv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rgemv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemv_debug__Float64x_SOURCES=../common/Rgemv.debug.cpp $(debug_srcdepends)
Rgemv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rgemv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rger_debug__Float64x_ref_SOURCES=../common/Rger.debug.cpp $(debug_srcdepends)
Rger_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rger_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rger_debug__Float64x_SOURCES=../common/Rger.debug.cpp $(debug_srcdepends)
Rger_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rger_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rnrm2_debug__Float64x_ref_SOURCES=../common/Rnrm2.debug.cpp $(debug_srcdepends)
Rnrm2_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rnrm2_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rnrm2_debug__Float64x_SOURCES=../common/Rnrm2.debug.cpp $(debug_srcdepends)
Rnrm2_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rnrm2_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrot_debug__Float64x_ref_SOURCES=../common/Rrot.debug.cpp $(debug_srcdepends)
Rrot_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rrot_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrot_debug__Float64x_SOURCES=../common/Rrot.debug.cpp $(debug_srcdepends)
Rrot_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rrot_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotg_debug__Float64x_ref_SOURCES=../common/Rrotg.debug.cpp $(debug_srcdepends)
Rrotg_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rrotg_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotg_debug__Float64x_SOURCES=../common/Rrotg.debug.cpp $(debug_srcdepends)
Rrotg_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rrotg_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotm_debug__Float64x_ref_SOURCES=../common/Rrotm.debug.cpp $(debug_srcdepends)
Rrotm_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rrotm_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotm_debug__Float64x_SOURCES=../common/Rrotm.debug.cpp $(debug_srcdepends)
Rrotm_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rrotm_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsbmv_debug__Float64x_ref_SOURCES=../common/Rsbmv.debug.cpp $(debug_srcdepends)
Rsbmv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rsbmv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsbmv_debug__Float64x_SOURCES=../common/Rsbmv.debug.cpp $(debug_srcdepends)
Rsbmv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rsbmv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rscal_debug__Float64x_ref_SOURCES=../common/Rscal.debug.cpp $(debug_srcdepends)
Rscal_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rscal_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rscal_debug__Float64x_SOURCES=../common/Rscal.debug.cpp $(debug_srcdepends)
Rscal_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rscal_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspmv_debug__Float64x_ref_SOURCES=../common/Rspmv.debug.cpp $(debug_srcdepends)
Rspmv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rspmv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspmv_debug__Float64x_SOURCES=../common/Rspmv.debug.cpp $(debug_srcdepends)
Rspmv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rspmv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr_debug__Float64x_ref_SOURCES=../common/Rspr.debug.cpp $(debug_srcdepends)
Rspr_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rspr_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr_debug__Float64x_SOURCES=../common/Rspr.debug.cpp $(debug_srcdepends)
Rspr_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rspr_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr2_debug__Float64x_ref_SOURCES=../common/Rspr2.debug.cpp $(debug_srcdepends)
Rspr2_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rspr2_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr2_debug__Float64x_SOURCES=../common/Rspr2.debug.cpp $(debug_srcdepends)
Rspr2_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rspr2_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rswap_debug__Float64x_ref_SOURCES=../common/Rswap.debug.cpp $(debug_srcdepends)
Rswap_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rswap_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rswap_debug__Float64x_SOURCES=../common/Rswap.debug.cpp $(debug_srcdepends)
Rswap_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rswap_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymm_debug__Float64x_ref_SOURCES=../common/Rsymm.debug.cpp $(debug_srcdepends)
Rsymm_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rsymm_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymm_debug__Float64x_SOURCES=../common/Rsymm.debug.cpp $(debug_srcdepends)
Rsymm_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rsymm_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymv_debug__Float64x_ref_SOURCES=../common/Rsymv.debug.cpp $(debug_srcdepends)
Rsymv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rsymv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymv_debug__Float64x_SOURCES=../common/Rsymv.debug.cpp $(debug_srcdepends)
Rsymv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rsymv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr_debug__Float64x_ref_SOURCES=../common/Rsyr.debug.cpp $(debug_srcdepends)
Rsyr_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr_debug__Float64x_SOURCES=../common/Rsyr.debug.cpp $(debug_srcdepends)
Rsyr_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rsyr_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2_debug__Float64x_ref_SOURCES=../common/Rsyr2.debug.cpp $(debug_srcdepends)
Rsyr2_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr2_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2_debug__Float64x_SOURCES=../common/Rsyr2.debug.cpp $(debug_srcdepends)
Rsyr2_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rsyr2_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2k_debug__Float64x_ref_SOURCES=../common/Rsyr2k.debug.cpp $(debug_srcdepends)
Rsyr2k_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr2k_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2k_debug__Float64x_SOURCES=../common/Rsyr2k.debug.cpp $(debug_srcdepends)
Rsyr2k_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rsyr2k_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyrk_debug__Float64x_ref_SOURCES=../common/Rsyrk.debug.cpp $(debug_srcdepends)
Rsyrk_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rsyrk_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyrk_debug__Float64x_SOURCES=../common/Rsyrk.debug.cpp $(debug_srcdepends)
Rsyrk_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rsyrk_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbmv_debug__Float64x_ref_SOURCES=../common/Rtbmv.debug.cpp $(debug_srcdepends)
Rtbmv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rtbmv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbmv_debug__Float64x_SOURCES=../common/Rtbmv.debug.cpp $(debug_srcdepends)
Rtbmv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rtbmv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbsv_debug__Float64x_ref_SOURCES=../common/Rtbsv.debug.cpp $(debug_srcdepends)
Rtbsv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rtbsv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbsv_debug__Float64x_SOURCES=../common/Rtbsv.debug.cpp $(debug_srcdepends)
Rtbsv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rtbsv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpmv_debug__Float64x_ref_SOURCES=../common/Rtpmv.debug.cpp $(debug_srcdepends)
Rtpmv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rtpmv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpmv_debug__Float64x_SOURCES=../common/Rtpmv.debug.cpp $(debug_srcdepends)
Rtpmv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rtpmv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpsv_debug__Float64x_ref_SOURCES=../common/Rtpsv.debug.cpp $(debug_srcdepends)
Rtpsv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rtpsv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpsv_debug__Float64x_SOURCES=../common/Rtpsv.debug.cpp $(debug_srcdepends)
Rtpsv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rtpsv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmm_debug__Float64x_ref_SOURCES=../common/Rtrmm.debug.cpp $(debug_srcdepends)
Rtrmm_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rtrmm_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmm_debug__Float64x_SOURCES=../common/Rtrmm.debug.cpp $(debug_srcdepends)
Rtrmm_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rtrmm_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmv_debug__Float64x_ref_SOURCES=../common/Rtrmv.debug.cpp $(debug_srcdepends)
Rtrmv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rtrmv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmv_debug__Float64x_SOURCES=../common/Rtrmv.debug.cpp $(debug_srcdepends)
Rtrmv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rtrmv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsm_debug__Float64x_ref_SOURCES=../common/Rtrsm.debug.cpp $(debug_srcdepends)
Rtrsm_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rtrsm_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsm_debug__Float64x_SOURCES=../common/Rtrsm.debug.cpp $(debug_srcdepends)
Rtrsm_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rtrsm_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsv_debug__Float64x_ref_SOURCES=../common/Rtrsv.debug.cpp $(debug_srcdepends)
Rtrsv_debug__Float64x_ref_CXXFLAGS=$(debug_cxxflags)
Rtrsv_debug__Float64x_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsv_debug__Float64x_SOURCES=../common/Rtrsv.debug.cpp $(debug_srcdepends)
Rtrsv_debug__Float64x_CXXFLAGS=$(debug_cxxflags)
Rtrsv_debug__Float64x_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)
