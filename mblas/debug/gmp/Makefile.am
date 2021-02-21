#CTBSV, CTPSV might have bugs or intrinsically inaccurate
#Rnrm2 should be checked more
check_PROGRAMS = $(mblas_gmp_test_PROGRAMS) 
TESTS = $(check_PROGRAMS)
mblas_gmp_test_PROGRAMS = \
arithmetic.debug_gmp \
complex.debug_gmp \
iCamax.debug.gmp_ref iCamax.debug.gmp   \
iRamax.debug.gmp_ref iRamax.debug.gmp   \
Mlsame.debug.gmp_ref Mlsame.debug.gmp   \
Mxerbla.debug.gmp_ref Mxerbla.debug.gmp \
CRrot.debug.gmp_ref  CRrot.debug.gmp    \
CRscal.debug.gmp_ref CRscal.debug.gmp   \
Caxpy.debug.gmp_ref  Caxpy.debug.gmp    \
Ccopy.debug.gmp_ref  Ccopy.debug.gmp    \
Cdotc.debug.gmp_ref  Cdotc.debug.gmp	\
Cdotu.debug.gmp_ref  Cdotu.debug.gmp	\
Cgbmv.debug.gmp_ref  Cgbmv.debug.gmp	\
Cgemm.debug.gmp_ref  Cgemm.debug.gmp	\
Cgemv.debug.gmp_ref  Cgemv.debug.gmp	\
Cgerc.debug.gmp_ref  Cgerc.debug.gmp	\
Cgeru.debug.gmp_ref  Cgeru.debug.gmp    \
Chbmv.debug.gmp_ref  Chbmv.debug.gmp	\
Chemm.debug.gmp_ref  Chemm.debug.gmp	\
Chemv.debug.gmp_ref  Chemv.debug.gmp	\
Cher.debug.gmp_ref   Cher.debug.gmp 	\
Cher2.debug.gmp_ref  Cher2.debug.gmp	\
Cher2k.debug.gmp_ref Cher2k.debug.gmp	\
Cherk.debug.gmp_ref  Cherk.debug.gmp	\
Chpmv.debug.gmp_ref  Chpmv.debug.gmp	\
Chpr.debug.gmp_ref   Chpr.debug.gmp	\
Chpr2.debug.gmp_ref  Chpr2.debug.gmp	\
Crotg.debug.gmp_ref  Crotg.debug.gmp	\
Cscal.debug.gmp_ref  Cscal.debug.gmp	\
Cswap.debug.gmp_ref  Cswap.debug.gmp	\
Csymm.debug.gmp_ref  Csymm.debug.gmp    \
Csyr2k.debug.gmp_ref Csyr2k.debug.gmp	\
Csyrk.debug.gmp_ref  Csyrk.debug.gmp	\
Ctbmv.debug.gmp_ref  Ctbmv.debug.gmp	\
Ctbsv.debug.gmp_ref  Ctbsv.debug.gmp	\
Ctpmv.debug.gmp_ref  Ctpmv.debug.gmp	\
Ctpsv.debug.gmp_ref  Ctpsv.debug.gmp	\
Ctrmm.debug.gmp_ref  Ctrmm.debug.gmp	\
Ctrmv.debug.gmp_ref  Ctrmv.debug.gmp	\
Ctrsm.debug.gmp_ref  Ctrsm.debug.gmp	\
Ctrsv.debug.gmp_ref  Ctrsv.debug.gmp	\
RCabs1.debug.gmp_ref RCabs1.debug.gmp	\
RCasum.debug.gmp_ref RCasum.debug.gmp	\
RCnrm2.debug.gmp_ref RCnrm2.debug.gmp	\
Rasum.debug.gmp_ref  Rasum.debug.gmp    \
Raxpy.debug.gmp_ref  Raxpy.debug.gmp	\
Rcopy.debug.gmp_ref  Rcopy.debug.gmp	\
Rdot.debug.gmp_ref   Rdot.debug.gmp	\
Rgbmv.debug.gmp_ref  Rgbmv.debug.gmp	\
Rgemm.debug.gmp_ref  Rgemm.debug.gmp	\
Rgemv.debug.gmp_ref  Rgemv.debug.gmp	\
Rger.debug.gmp_ref   Rger.debug.gmp	\
Rnrm2.debug.gmp_ref  Rnrm2.debug.gmp	\
Rrot.debug.gmp_ref   Rrot.debug.gmp	\
Rrotg.debug.gmp_ref  Rrotg.debug.gmp	\
Rrotm.debug.gmp_ref  Rrotm.debug.gmp	\
Rsbmv.debug.gmp_ref  Rsbmv.debug.gmp	\
Rscal.debug.gmp_ref  Rscal.debug.gmp	\
Rspmv.debug.gmp_ref  Rspmv.debug.gmp    \
Rspr.debug.gmp_ref   Rspr.debug.gmp	\
Rspr2.debug.gmp_ref  Rspr2.debug.gmp	\
Rswap.debug.gmp_ref  Rswap.debug.gmp	\
Rsymm.debug.gmp_ref  Rsymm.debug.gmp	\
Rsymv.debug.gmp_ref  Rsymv.debug.gmp	\
Rsyr.debug.gmp_ref   Rsyr.debug.gmp	\
Rsyr2.debug.gmp_ref  Rsyr2.debug.gmp	\
Rsyr2k.debug.gmp_ref Rsyr2k.debug.gmp	\
Rsyrk.debug.gmp_ref  Rsyrk.debug.gmp	\
Rtbmv.debug.gmp_ref  Rtbmv.debug.gmp	\
Rtbsv.debug.gmp_ref  Rtbsv.debug.gmp	\
Rtpmv.debug.gmp_ref  Rtpmv.debug.gmp	\
Rtpsv.debug.gmp_ref  Rtpsv.debug.gmp	\
Rtrmm.debug.gmp_ref  Rtrmm.debug.gmp    \
Rtrmv.debug.gmp_ref  Rtrmv.debug.gmp	\
Rtrsm.debug.gmp_ref  Rtrsm.debug.gmp	\
Rtrsv.debug.gmp_ref  Rtrsv.debug.gmp

mblas_gmp_testdir = $(prefix)/share/mpack/testing/gmp
ref_libdepends = -L$(top_builddir)/mblas/reference -lmblas_gmp_ref -lmblas_mpfr_ref
opt_libdepends = -L$(top_builddir)/mblas/optimized/gmp -lmblas_gmp -L$(top_builddir)/mblas/reference -lmblas_mpfr_ref
mplibs=-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp
refblas=-L$(top_builddir)/external/i/LAPACK/lib -lblas $(FCLIBS)

debug_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPACK_BUILD_WITH_GMP___

if !IS_MACOSX
debug_srcdepends = ../common/xerbla.cpp ../common/mpack.debug.cpp ../common/Mxerbla.override.cpp
endif

if IS_MACOSX
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPACK_BUILD_WITH_GMP___
libmxerbla_override_la_SOURCES =  ../common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS= -Wl,-flat_namespace,-undefined,dynamic_lookup
debug_srcdepends = ../common/xerbla.cpp ../common/mpack.debug.cpp
mplibs+=-Wl,-flat_namespace,-undefined,dynamic_lookup
TESTS_ENVIRONMENT=DYLD_INSERT_LIBRARIES=.libs/libmxerbla_override.0.dylib
endif

##
arithmetic_debug_gmp_SOURCES=../common/arithmetic.debug.cpp $(debug_srcdepends)
arithmetic_debug_gmp_CXXFLAGS=$(debug_cxxflags)
arithmetic_debug_gmp_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

complex_debug_gmp_SOURCES=../common/complex.debug.cpp $(debug_srcdepends)
complex_debug_gmp_CXXFLAGS=$(debug_cxxflags)
complex_debug_gmp_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_debug_gmp_ref_SOURCES=../common/iCamax.debug.cpp $(debug_srcdepends)
iCamax_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
iCamax_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_debug_gmp_SOURCES=../common/iCamax.debug.cpp $(debug_srcdepends)
iCamax_debug_gmp_CXXFLAGS=$(debug_cxxflags)
iCamax_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iRamax_debug_gmp_ref_SOURCES=../common/iRamax.debug.cpp $(debug_srcdepends)
iRamax_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
iRamax_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iRamax_debug_gmp_SOURCES=../common/iRamax.debug.cpp $(debug_srcdepends)
iRamax_debug_gmp_CXXFLAGS=$(debug_cxxflags)
iRamax_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mlsame_debug_gmp_ref_SOURCES=../common/Mlsame.debug.cpp $(debug_srcdepends)
Mlsame_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Mlsame_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mlsame_debug_gmp_SOURCES=../common/Mlsame.debug.cpp $(debug_srcdepends)
Mlsame_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Mlsame_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mxerbla_debug_gmp_ref_SOURCES=../common/Mxerbla.debug.cpp $(debug_srcdepends)
Mxerbla_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Mxerbla_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mxerbla_debug_gmp_SOURCES=../common/Mxerbla.debug.cpp $(debug_srcdepends)
Mxerbla_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Mxerbla_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRrot_debug_gmp_ref_SOURCES=../common/CRrot.debug.cpp $(debug_srcdepends)
CRrot_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
CRrot_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRrot_debug_gmp_SOURCES=../common/CRrot.debug.cpp $(debug_srcdepends)
CRrot_debug_gmp_CXXFLAGS=$(debug_cxxflags)
CRrot_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRscal_debug_gmp_ref_SOURCES=../common/CRscal.debug.cpp $(debug_srcdepends)
CRscal_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
CRscal_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRscal_debug_gmp_SOURCES=../common/CRscal.debug.cpp $(debug_srcdepends)
CRscal_debug_gmp_CXXFLAGS=$(debug_cxxflags)
CRscal_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_debug_gmp_SOURCES=../common/Caxpy.debug.cpp $(debug_srcdepends)
Caxpy_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Caxpy_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_debug_gmp_ref_SOURCES=../common/Caxpy.debug.cpp $(debug_srcdepends)
Caxpy_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Caxpy_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ccopy_debug_gmp_SOURCES=../common/Ccopy.debug.cpp $(debug_srcdepends)
Ccopy_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Ccopy_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ccopy_debug_gmp_ref_SOURCES=../common/Ccopy.debug.cpp $(debug_srcdepends)
Ccopy_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Ccopy_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotc_debug_gmp_SOURCES=../common/Cdotc.debug.cpp $(debug_srcdepends)
Cdotc_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Cdotc_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotc_debug_gmp_ref_SOURCES=../common/Cdotc.debug.cpp $(debug_srcdepends)
Cdotc_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Cdotc_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotu_debug_gmp_SOURCES=../common/Cdotu.debug.cpp $(debug_srcdepends)
Cdotu_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Cdotu_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotu_debug_gmp_ref_SOURCES=../common/Cdotu.debug.cpp $(debug_srcdepends)
Cdotu_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Cdotu_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgbmv_debug_gmp_SOURCES=../common/Cgbmv.debug.cpp $(debug_srcdepends)
Cgbmv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Cgbmv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgbmv_debug_gmp_ref_SOURCES=../common/Cgbmv.debug.cpp $(debug_srcdepends)
Cgbmv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Cgbmv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemm_debug_gmp_SOURCES=../common/Cgemm.debug.cpp $(debug_srcdepends)
Cgemm_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Cgemm_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemm_debug_gmp_ref_SOURCES=../common/Cgemm.debug.cpp $(debug_srcdepends)
Cgemm_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Cgemm_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemv_debug_gmp_SOURCES=../common/Cgemv.debug.cpp $(debug_srcdepends)
Cgemv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Cgemv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemv_debug_gmp_ref_SOURCES=../common/Cgemv.debug.cpp $(debug_srcdepends)
Cgemv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Cgemv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgerc_debug_gmp_SOURCES=../common/Cgerc.debug.cpp $(debug_srcdepends)
Cgerc_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Cgerc_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgerc_debug_gmp_ref_SOURCES=../common/Cgerc.debug.cpp $(debug_srcdepends)
Cgerc_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Cgerc_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgeru_debug_gmp_SOURCES=../common/Cgeru.debug.cpp $(debug_srcdepends)
Cgeru_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Cgeru_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgeru_debug_gmp_ref_SOURCES=../common/Cgeru.debug.cpp $(debug_srcdepends)
Cgeru_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Cgeru_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chbmv_debug_gmp_SOURCES=../common/Chbmv.debug.cpp $(debug_srcdepends)
Chbmv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Chbmv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chbmv_debug_gmp_ref_SOURCES=../common/Chbmv.debug.cpp $(debug_srcdepends)
Chbmv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Chbmv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemm_debug_gmp_SOURCES=../common/Chemm.debug.cpp $(debug_srcdepends)
Chemm_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Chemm_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemm_debug_gmp_ref_SOURCES=../common/Chemm.debug.cpp $(debug_srcdepends)
Chemm_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Chemm_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemv_debug_gmp_SOURCES=../common/Chemv.debug.cpp $(debug_srcdepends)
Chemv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Chemv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemv_debug_gmp_ref_SOURCES=../common/Chemv.debug.cpp $(debug_srcdepends)
Chemv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Chemv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher_debug_gmp_SOURCES=../common/Cher.debug.cpp $(debug_srcdepends)
Cher_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Cher_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher_debug_gmp_ref_SOURCES=../common/Cher.debug.cpp $(debug_srcdepends)
Cher_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Cher_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2_debug_gmp_SOURCES=../common/Cher2.debug.cpp $(debug_srcdepends)
Cher2_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Cher2_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2_debug_gmp_ref_SOURCES=../common/Cher2.debug.cpp $(debug_srcdepends)
Cher2_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Cher2_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2k_debug_gmp_SOURCES=../common/Cher2k.debug.cpp $(debug_srcdepends)
Cher2k_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Cher2k_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2k_debug_gmp_ref_SOURCES=../common/Cher2k.debug.cpp $(debug_srcdepends)
Cher2k_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Cher2k_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cherk_debug_gmp_SOURCES=../common/Cherk.debug.cpp $(debug_srcdepends)
Cherk_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Cherk_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cherk_debug_gmp_ref_SOURCES=../common/Cherk.debug.cpp $(debug_srcdepends)
Cherk_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Cherk_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpmv_debug_gmp_SOURCES=../common/Chpmv.debug.cpp $(debug_srcdepends)
Chpmv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Chpmv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpmv_debug_gmp_ref_SOURCES=../common/Chpmv.debug.cpp $(debug_srcdepends)
Chpmv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Chpmv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr_debug_gmp_SOURCES=../common/Chpr.debug.cpp $(debug_srcdepends)
Chpr_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Chpr_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr_debug_gmp_ref_SOURCES=../common/Chpr.debug.cpp $(debug_srcdepends)
Chpr_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Chpr_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr2_debug_gmp_SOURCES=../common/Chpr2.debug.cpp $(debug_srcdepends)
Chpr2_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Chpr2_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr2_debug_gmp_ref_SOURCES=../common/Chpr2.debug.cpp $(debug_srcdepends)
Chpr2_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Chpr2_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Crotg_debug_gmp_SOURCES=../common/Crotg.debug.cpp $(debug_srcdepends)
Crotg_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Crotg_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Crotg_debug_gmp_ref_SOURCES=../common/Crotg.debug.cpp $(debug_srcdepends)
Crotg_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Crotg_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cscal_debug_gmp_SOURCES=../common/Cscal.debug.cpp $(debug_srcdepends)
Cscal_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Cscal_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cscal_debug_gmp_ref_SOURCES=../common/Cscal.debug.cpp $(debug_srcdepends)
Cscal_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Cscal_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cswap_debug_gmp_SOURCES=../common/Cswap.debug.cpp $(debug_srcdepends)
Cswap_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Cswap_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cswap_debug_gmp_ref_SOURCES=../common/Cswap.debug.cpp $(debug_srcdepends)
Cswap_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Cswap_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csymm_debug_gmp_SOURCES=../common/Csymm.debug.cpp $(debug_srcdepends)
Csymm_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Csymm_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csymm_debug_gmp_ref_SOURCES=../common/Csymm.debug.cpp $(debug_srcdepends)
Csymm_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Csymm_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyr2k_debug_gmp_SOURCES=../common/Csyr2k.debug.cpp $(debug_srcdepends)
Csyr2k_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Csyr2k_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyr2k_debug_gmp_ref_SOURCES=../common/Csyr2k.debug.cpp $(debug_srcdepends)
Csyr2k_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Csyr2k_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyrk_debug_gmp_SOURCES=../common/Csyrk.debug.cpp $(debug_srcdepends)
Csyrk_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Csyrk_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyrk_debug_gmp_ref_SOURCES=../common/Csyrk.debug.cpp $(debug_srcdepends)
Csyrk_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Csyrk_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbmv_debug_gmp_SOURCES=../common/Ctbmv.debug.cpp $(debug_srcdepends)
Ctbmv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Ctbmv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbmv_debug_gmp_ref_SOURCES=../common/Ctbmv.debug.cpp $(debug_srcdepends)
Ctbmv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Ctbmv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbsv_debug_gmp_SOURCES=../common/Ctbsv.debug.cpp $(debug_srcdepends)
Ctbsv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Ctbsv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbsv_debug_gmp_ref_SOURCES=../common/Ctbsv.debug.cpp $(debug_srcdepends)
Ctbsv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Ctbsv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpmv_debug_gmp_SOURCES=../common/Ctpmv.debug.cpp $(debug_srcdepends)
Ctpmv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Ctpmv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpmv_debug_gmp_ref_SOURCES=../common/Ctpmv.debug.cpp $(debug_srcdepends)
Ctpmv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Ctpmv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpsv_debug_gmp_SOURCES=../common/Ctpsv.debug.cpp $(debug_srcdepends)
Ctpsv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Ctpsv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpsv_debug_gmp_ref_SOURCES=../common/Ctpsv.debug.cpp $(debug_srcdepends)
Ctpsv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Ctpsv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmm_debug_gmp_SOURCES=../common/Ctrmm.debug.cpp $(debug_srcdepends)
Ctrmm_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Ctrmm_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmm_debug_gmp_ref_SOURCES=../common/Ctrmm.debug.cpp $(debug_srcdepends)
Ctrmm_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Ctrmm_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmv_debug_gmp_SOURCES=../common/Ctrmv.debug.cpp $(debug_srcdepends)
Ctrmv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Ctrmv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmv_debug_gmp_ref_SOURCES=../common/Ctrmv.debug.cpp $(debug_srcdepends)
Ctrmv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Ctrmv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsm_debug_gmp_SOURCES=../common/Ctrsm.debug.cpp $(debug_srcdepends)
Ctrsm_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Ctrsm_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsm_debug_gmp_ref_SOURCES=../common/Ctrsm.debug.cpp $(debug_srcdepends)
Ctrsm_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Ctrsm_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsv_debug_gmp_SOURCES=../common/Ctrsv.debug.cpp $(debug_srcdepends)
Ctrsv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Ctrsv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsv_debug_gmp_ref_SOURCES=../common/Ctrsv.debug.cpp $(debug_srcdepends)
Ctrsv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Ctrsv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_debug_gmp_ref_SOURCES=../common/RCabs1.debug.cpp $(debug_srcdepends)
RCabs1_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
RCabs1_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_debug_gmp_SOURCES=../common/RCabs1.debug.cpp $(debug_srcdepends)
RCabs1_debug_gmp_CXXFLAGS=$(debug_cxxflags)
RCabs1_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCasum_debug_gmp_ref_SOURCES=../common/RCasum.debug.cpp $(debug_srcdepends)
RCasum_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
RCasum_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCasum_debug_gmp_SOURCES=../common/RCasum.debug.cpp $(debug_srcdepends)
RCasum_debug_gmp_CXXFLAGS=$(debug_cxxflags)
RCasum_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCnrm2_debug_gmp_ref_SOURCES=../common/RCnrm2.debug.cpp $(debug_srcdepends)
RCnrm2_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
RCnrm2_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCnrm2_debug_gmp_SOURCES=../common/RCnrm2.debug.cpp $(debug_srcdepends)
RCnrm2_debug_gmp_CXXFLAGS=$(debug_cxxflags)
RCnrm2_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rasum_debug_gmp_ref_SOURCES=../common/Rasum.debug.cpp $(debug_srcdepends)
Rasum_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rasum_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rasum_debug_gmp_SOURCES=../common/Rasum.debug.cpp $(debug_srcdepends)
Rasum_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rasum_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Raxpy_debug_gmp_ref_SOURCES=../common/Raxpy.debug.cpp $(debug_srcdepends)
Raxpy_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Raxpy_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Raxpy_debug_gmp_SOURCES=../common/Raxpy.debug.cpp $(debug_srcdepends)
Raxpy_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Raxpy_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rcopy_debug_gmp_ref_SOURCES=../common/Rcopy.debug.cpp $(debug_srcdepends)
Rcopy_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rcopy_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rcopy_debug_gmp_SOURCES=../common/Rcopy.debug.cpp $(debug_srcdepends)
Rcopy_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rcopy_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rdot_debug_gmp_ref_SOURCES=../common/Rdot.debug.cpp $(debug_srcdepends)
Rdot_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rdot_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rdot_debug_gmp_SOURCES=../common/Rdot.debug.cpp $(debug_srcdepends)
Rdot_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rdot_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgbmv_debug_gmp_ref_SOURCES=../common/Rgbmv.debug.cpp $(debug_srcdepends)
Rgbmv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rgbmv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgbmv_debug_gmp_SOURCES=../common/Rgbmv.debug.cpp $(debug_srcdepends)
Rgbmv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rgbmv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemm_debug_gmp_ref_SOURCES=../common/Rgemm.debug.cpp $(debug_srcdepends)
Rgemm_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rgemm_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemm_debug_gmp_SOURCES=../common/Rgemm.debug.cpp $(debug_srcdepends)
Rgemm_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rgemm_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemv_debug_gmp_ref_SOURCES=../common/Rgemv.debug.cpp $(debug_srcdepends)
Rgemv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rgemv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemv_debug_gmp_SOURCES=../common/Rgemv.debug.cpp $(debug_srcdepends)
Rgemv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rgemv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rger_debug_gmp_ref_SOURCES=../common/Rger.debug.cpp $(debug_srcdepends)
Rger_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rger_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rger_debug_gmp_SOURCES=../common/Rger.debug.cpp $(debug_srcdepends)
Rger_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rger_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rnrm2_debug_gmp_ref_SOURCES=../common/Rnrm2.debug.cpp $(debug_srcdepends)
Rnrm2_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rnrm2_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rnrm2_debug_gmp_SOURCES=../common/Rnrm2.debug.cpp $(debug_srcdepends)
Rnrm2_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rnrm2_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrot_debug_gmp_ref_SOURCES=../common/Rrot.debug.cpp $(debug_srcdepends)
Rrot_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rrot_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrot_debug_gmp_SOURCES=../common/Rrot.debug.cpp $(debug_srcdepends)
Rrot_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rrot_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotg_debug_gmp_ref_SOURCES=../common/Rrotg.debug.cpp $(debug_srcdepends)
Rrotg_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rrotg_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotg_debug_gmp_SOURCES=../common/Rrotg.debug.cpp $(debug_srcdepends)
Rrotg_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rrotg_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotm_debug_gmp_ref_SOURCES=../common/Rrotm.debug.cpp $(debug_srcdepends)
Rrotm_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rrotm_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotm_debug_gmp_SOURCES=../common/Rrotm.debug.cpp $(debug_srcdepends)
Rrotm_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rrotm_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsbmv_debug_gmp_ref_SOURCES=../common/Rsbmv.debug.cpp $(debug_srcdepends)
Rsbmv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rsbmv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsbmv_debug_gmp_SOURCES=../common/Rsbmv.debug.cpp $(debug_srcdepends)
Rsbmv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rsbmv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rscal_debug_gmp_ref_SOURCES=../common/Rscal.debug.cpp $(debug_srcdepends)
Rscal_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rscal_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rscal_debug_gmp_SOURCES=../common/Rscal.debug.cpp $(debug_srcdepends)
Rscal_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rscal_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspmv_debug_gmp_ref_SOURCES=../common/Rspmv.debug.cpp $(debug_srcdepends)
Rspmv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rspmv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspmv_debug_gmp_SOURCES=../common/Rspmv.debug.cpp $(debug_srcdepends)
Rspmv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rspmv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr_debug_gmp_ref_SOURCES=../common/Rspr.debug.cpp $(debug_srcdepends)
Rspr_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rspr_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr_debug_gmp_SOURCES=../common/Rspr.debug.cpp $(debug_srcdepends)
Rspr_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rspr_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr2_debug_gmp_ref_SOURCES=../common/Rspr2.debug.cpp $(debug_srcdepends)
Rspr2_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rspr2_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr2_debug_gmp_SOURCES=../common/Rspr2.debug.cpp $(debug_srcdepends)
Rspr2_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rspr2_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rswap_debug_gmp_ref_SOURCES=../common/Rswap.debug.cpp $(debug_srcdepends)
Rswap_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rswap_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rswap_debug_gmp_SOURCES=../common/Rswap.debug.cpp $(debug_srcdepends)
Rswap_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rswap_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymm_debug_gmp_ref_SOURCES=../common/Rsymm.debug.cpp $(debug_srcdepends)
Rsymm_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rsymm_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymm_debug_gmp_SOURCES=../common/Rsymm.debug.cpp $(debug_srcdepends)
Rsymm_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rsymm_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymv_debug_gmp_ref_SOURCES=../common/Rsymv.debug.cpp $(debug_srcdepends)
Rsymv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rsymv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymv_debug_gmp_SOURCES=../common/Rsymv.debug.cpp $(debug_srcdepends)
Rsymv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rsymv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr_debug_gmp_ref_SOURCES=../common/Rsyr.debug.cpp $(debug_srcdepends)
Rsyr_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr_debug_gmp_SOURCES=../common/Rsyr.debug.cpp $(debug_srcdepends)
Rsyr_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rsyr_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2_debug_gmp_ref_SOURCES=../common/Rsyr2.debug.cpp $(debug_srcdepends)
Rsyr2_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr2_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2_debug_gmp_SOURCES=../common/Rsyr2.debug.cpp $(debug_srcdepends)
Rsyr2_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rsyr2_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2k_debug_gmp_ref_SOURCES=../common/Rsyr2k.debug.cpp $(debug_srcdepends)
Rsyr2k_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr2k_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2k_debug_gmp_SOURCES=../common/Rsyr2k.debug.cpp $(debug_srcdepends)
Rsyr2k_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rsyr2k_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyrk_debug_gmp_ref_SOURCES=../common/Rsyrk.debug.cpp $(debug_srcdepends)
Rsyrk_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rsyrk_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyrk_debug_gmp_SOURCES=../common/Rsyrk.debug.cpp $(debug_srcdepends)
Rsyrk_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rsyrk_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbmv_debug_gmp_ref_SOURCES=../common/Rtbmv.debug.cpp $(debug_srcdepends)
Rtbmv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rtbmv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbmv_debug_gmp_SOURCES=../common/Rtbmv.debug.cpp $(debug_srcdepends)
Rtbmv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rtbmv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbsv_debug_gmp_ref_SOURCES=../common/Rtbsv.debug.cpp $(debug_srcdepends)
Rtbsv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rtbsv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbsv_debug_gmp_SOURCES=../common/Rtbsv.debug.cpp $(debug_srcdepends)
Rtbsv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rtbsv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpmv_debug_gmp_ref_SOURCES=../common/Rtpmv.debug.cpp $(debug_srcdepends)
Rtpmv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rtpmv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpmv_debug_gmp_SOURCES=../common/Rtpmv.debug.cpp $(debug_srcdepends)
Rtpmv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rtpmv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpsv_debug_gmp_ref_SOURCES=../common/Rtpsv.debug.cpp $(debug_srcdepends)
Rtpsv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rtpsv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpsv_debug_gmp_SOURCES=../common/Rtpsv.debug.cpp $(debug_srcdepends)
Rtpsv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rtpsv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmm_debug_gmp_ref_SOURCES=../common/Rtrmm.debug.cpp $(debug_srcdepends)
Rtrmm_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rtrmm_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmm_debug_gmp_SOURCES=../common/Rtrmm.debug.cpp $(debug_srcdepends)
Rtrmm_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rtrmm_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmv_debug_gmp_ref_SOURCES=../common/Rtrmv.debug.cpp $(debug_srcdepends)
Rtrmv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rtrmv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmv_debug_gmp_SOURCES=../common/Rtrmv.debug.cpp $(debug_srcdepends)
Rtrmv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rtrmv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsm_debug_gmp_ref_SOURCES=../common/Rtrsm.debug.cpp $(debug_srcdepends)
Rtrsm_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rtrsm_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsm_debug_gmp_SOURCES=../common/Rtrsm.debug.cpp $(debug_srcdepends)
Rtrsm_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rtrsm_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsv_debug_gmp_ref_SOURCES=../common/Rtrsv.debug.cpp $(debug_srcdepends)
Rtrsv_debug_gmp_ref_CXXFLAGS=$(debug_cxxflags)
Rtrsv_debug_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsv_debug_gmp_SOURCES=../common/Rtrsv.debug.cpp $(debug_srcdepends)
Rtrsv_debug_gmp_CXXFLAGS=$(debug_cxxflags)
Rtrsv_debug_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)
