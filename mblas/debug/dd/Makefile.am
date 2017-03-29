#CTBSV, CTPSV might have bugs or intrinsically inaccurate
#Rnrm2 should be checked more
check_PROGRAMS = $(mblas_dd_test_PROGRAMS) 
TESTS = $(check_PROGRAMS)
mblas_dd_test_PROGRAMS = \
arithmetic.debug_dd \
complex.debug_dd \
iCamax.debug.dd_ref iCamax.debug.dd   \
iRamax.debug.dd_ref iRamax.debug.dd   \
Mlsame.debug.dd_ref Mlsame.debug.dd   \
Mxerbla.debug.dd_ref Mxerbla.debug.dd \
CRrot.debug.dd_ref  CRrot.debug.dd    \
CRscal.debug.dd_ref CRscal.debug.dd   \
Caxpy.debug.dd_ref  Caxpy.debug.dd    \
Ccopy.debug.dd_ref  Ccopy.debug.dd    \
Cdotc.debug.dd_ref  Cdotc.debug.dd	\
Cdotu.debug.dd_ref  Cdotu.debug.dd	\
Cgbmv.debug.dd_ref  Cgbmv.debug.dd	\
Cgemm.debug.dd_ref  Cgemm.debug.dd	\
Cgemv.debug.dd_ref  Cgemv.debug.dd	\
Cgerc.debug.dd_ref  Cgerc.debug.dd	\
Cgeru.debug.dd_ref  Cgeru.debug.dd    \
Chbmv.debug.dd_ref  Chbmv.debug.dd	\
Chemm.debug.dd_ref  Chemm.debug.dd	\
Chemv.debug.dd_ref  Chemv.debug.dd	\
Cher.debug.dd_ref   Cher.debug.dd 	\
Cher2.debug.dd_ref  Cher2.debug.dd	\
Cher2k.debug.dd_ref Cher2k.debug.dd	\
Cherk.debug.dd_ref  Cherk.debug.dd	\
Chpmv.debug.dd_ref  Chpmv.debug.dd	\
Chpr.debug.dd_ref   Chpr.debug.dd	\
Chpr2.debug.dd_ref  Chpr2.debug.dd	\
Crotg.debug.dd_ref  Crotg.debug.dd	\
Cscal.debug.dd_ref  Cscal.debug.dd	\
Cswap.debug.dd_ref  Cswap.debug.dd	\
Csymm.debug.dd_ref  Csymm.debug.dd    \
Csyr2k.debug.dd_ref Csyr2k.debug.dd	\
Csyrk.debug.dd_ref  Csyrk.debug.dd	\
Ctbmv.debug.dd_ref  Ctbmv.debug.dd	\
Ctbsv.debug.dd_ref  Ctbsv.debug.dd	\
Ctpmv.debug.dd_ref  Ctpmv.debug.dd	\
Ctpsv.debug.dd_ref  Ctpsv.debug.dd	\
Ctrmm.debug.dd_ref  Ctrmm.debug.dd	\
Ctrmv.debug.dd_ref  Ctrmv.debug.dd	\
Ctrsm.debug.dd_ref  Ctrsm.debug.dd	\
Ctrsv.debug.dd_ref  Ctrsv.debug.dd	\
RCabs1.debug.dd_ref RCabs1.debug.dd	\
RCasum.debug.dd_ref RCasum.debug.dd	\
RCnrm2.debug.dd_ref RCnrm2.debug.dd	\
Rasum.debug.dd_ref  Rasum.debug.dd    \
Raxpy.debug.dd_ref  Raxpy.debug.dd	\
Rcopy.debug.dd_ref  Rcopy.debug.dd	\
Rdot.debug.dd_ref   Rdot.debug.dd	\
Rgbmv.debug.dd_ref  Rgbmv.debug.dd	\
Rgemm.debug.dd_ref  Rgemm.debug.dd	\
Rgemv.debug.dd_ref  Rgemv.debug.dd	\
Rger.debug.dd_ref   Rger.debug.dd	\
Rnrm2.debug.dd_ref  Rnrm2.debug.dd	\
Rrot.debug.dd_ref   Rrot.debug.dd	\
Rrotg.debug.dd_ref  Rrotg.debug.dd	\
Rrotm.debug.dd_ref  Rrotm.debug.dd	\
Rsbmv.debug.dd_ref  Rsbmv.debug.dd	\
Rscal.debug.dd_ref  Rscal.debug.dd	\
Rspmv.debug.dd_ref  Rspmv.debug.dd    \
Rspr.debug.dd_ref   Rspr.debug.dd	\
Rspr2.debug.dd_ref  Rspr2.debug.dd	\
Rswap.debug.dd_ref  Rswap.debug.dd	\
Rsymm.debug.dd_ref  Rsymm.debug.dd	\
Rsymv.debug.dd_ref  Rsymv.debug.dd	\
Rsyr.debug.dd_ref   Rsyr.debug.dd	\
Rsyr2.debug.dd_ref  Rsyr2.debug.dd	\
Rsyr2k.debug.dd_ref Rsyr2k.debug.dd	\
Rsyrk.debug.dd_ref  Rsyrk.debug.dd	\
Rtbmv.debug.dd_ref  Rtbmv.debug.dd	\
Rtbsv.debug.dd_ref  Rtbsv.debug.dd	\
Rtpmv.debug.dd_ref  Rtpmv.debug.dd	\
Rtpsv.debug.dd_ref  Rtpsv.debug.dd	\
Rtrmm.debug.dd_ref  Rtrmm.debug.dd    \
Rtrmv.debug.dd_ref  Rtrmv.debug.dd	\
Rtrsm.debug.dd_ref  Rtrsm.debug.dd	\
Rtrsv.debug.dd_ref  Rtrsv.debug.dd

mblas_dd_testdir = $(prefix)/share/mpack/testing/dd
ref_libdepends = -L$(top_builddir)/mblas/reference -L$(top_builddir)/mblas/reference -lmblas_dd_ref -lmblas_mpfr_ref
opt_libdepends = -L$(top_builddir)/mblas/reference -L$(top_builddir)/mblas/optimized/dd -lmblas_dd -lmblas_mpfr_ref
mplibs=-L$(top_srcdir)/$(MPC_LIBDIR) -L$(top_srcdir)/$(MPFR_LIBDIR) -L$(top_srcdir)/$(GMP_LIBDIR) -L$(top_builddir)/mpfrc++ -lmpfrcxx -lmpc -lmpfr -lgmpxx -lgmp -L$(top_srcdir)/$(QD_LIBDIR) -lqd
refblas=-L$(top_builddir)/external/i/LAPACK/lib -lblas $(FCLIBS)

debug_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/$(GMP_INCLUDEDIR) -I$(top_srcdir)/$(MPFR_INCLUDEDIR) -I$(top_srcdir)/$(MPC_INCLUDEDIR) -I$(top_srcdir)/$(QD_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPACK_BUILD_WITH_DD___

if !IS_MACOSX
debug_srcdepends = ../common/xerbla.cpp ../common/mpack.debug.cpp ../common/Mxerbla.override.cpp
endif

if IS_MACOSX
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPACK_BUILD_WITH_DD___
libmxerbla_override_la_SOURCES =  ../common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS= -Wl,-flat_namespace,-undefined,dynamic_lookup
debug_srcdepends = ../common/xerbla.cpp ../common/mpack.debug.cpp
mplibs+=-Wl,-flat_namespace,-undefined,dynamic_lookup
TESTS_ENVIRONMENT=DYLD_INSERT_LIBRARIES=.libs/libmxerbla_override.0.dylib
endif

##
arithmetic_debug_dd_SOURCES=../common/arithmetic.debug.cpp $(debug_srcdepends)
arithmetic_debug_dd_CXXFLAGS=$(debug_cxxflags)
arithmetic_debug_dd_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

complex_debug_dd_SOURCES=../common/complex.debug.cpp $(debug_srcdepends)
complex_debug_dd_CXXFLAGS=$(debug_cxxflags)
complex_debug_dd_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_debug_dd_ref_SOURCES=../common/iCamax.debug.cpp $(debug_srcdepends)
iCamax_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
iCamax_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_debug_dd_SOURCES=../common/iCamax.debug.cpp $(debug_srcdepends)
iCamax_debug_dd_CXXFLAGS=$(debug_cxxflags)
iCamax_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iRamax_debug_dd_ref_SOURCES=../common/iRamax.debug.cpp $(debug_srcdepends)
iRamax_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
iRamax_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iRamax_debug_dd_SOURCES=../common/iRamax.debug.cpp $(debug_srcdepends)
iRamax_debug_dd_CXXFLAGS=$(debug_cxxflags)
iRamax_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mlsame_debug_dd_ref_SOURCES=../common/Mlsame.debug.cpp $(debug_srcdepends)
Mlsame_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Mlsame_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mlsame_debug_dd_SOURCES=../common/Mlsame.debug.cpp $(debug_srcdepends)
Mlsame_debug_dd_CXXFLAGS=$(debug_cxxflags)
Mlsame_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mxerbla_debug_dd_ref_SOURCES=../common/Mxerbla.debug.cpp $(debug_srcdepends)
Mxerbla_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Mxerbla_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mxerbla_debug_dd_SOURCES=../common/Mxerbla.debug.cpp $(debug_srcdepends)
Mxerbla_debug_dd_CXXFLAGS=$(debug_cxxflags)
Mxerbla_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRrot_debug_dd_ref_SOURCES=../common/CRrot.debug.cpp $(debug_srcdepends)
CRrot_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
CRrot_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRrot_debug_dd_SOURCES=../common/CRrot.debug.cpp $(debug_srcdepends)
CRrot_debug_dd_CXXFLAGS=$(debug_cxxflags)
CRrot_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRscal_debug_dd_ref_SOURCES=../common/CRscal.debug.cpp $(debug_srcdepends)
CRscal_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
CRscal_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRscal_debug_dd_SOURCES=../common/CRscal.debug.cpp $(debug_srcdepends)
CRscal_debug_dd_CXXFLAGS=$(debug_cxxflags)
CRscal_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_debug_dd_SOURCES=../common/Caxpy.debug.cpp $(debug_srcdepends)
Caxpy_debug_dd_CXXFLAGS=$(debug_cxxflags)
Caxpy_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_debug_dd_ref_SOURCES=../common/Caxpy.debug.cpp $(debug_srcdepends)
Caxpy_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Caxpy_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ccopy_debug_dd_SOURCES=../common/Ccopy.debug.cpp $(debug_srcdepends)
Ccopy_debug_dd_CXXFLAGS=$(debug_cxxflags)
Ccopy_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ccopy_debug_dd_ref_SOURCES=../common/Ccopy.debug.cpp $(debug_srcdepends)
Ccopy_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Ccopy_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotc_debug_dd_SOURCES=../common/Cdotc.debug.cpp $(debug_srcdepends)
Cdotc_debug_dd_CXXFLAGS=$(debug_cxxflags)
Cdotc_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotc_debug_dd_ref_SOURCES=../common/Cdotc.debug.cpp $(debug_srcdepends)
Cdotc_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Cdotc_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotu_debug_dd_SOURCES=../common/Cdotu.debug.cpp $(debug_srcdepends)
Cdotu_debug_dd_CXXFLAGS=$(debug_cxxflags)
Cdotu_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotu_debug_dd_ref_SOURCES=../common/Cdotu.debug.cpp $(debug_srcdepends)
Cdotu_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Cdotu_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgbmv_debug_dd_SOURCES=../common/Cgbmv.debug.cpp $(debug_srcdepends)
Cgbmv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Cgbmv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgbmv_debug_dd_ref_SOURCES=../common/Cgbmv.debug.cpp $(debug_srcdepends)
Cgbmv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Cgbmv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemm_debug_dd_SOURCES=../common/Cgemm.debug.cpp $(debug_srcdepends)
Cgemm_debug_dd_CXXFLAGS=$(debug_cxxflags)
Cgemm_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemm_debug_dd_ref_SOURCES=../common/Cgemm.debug.cpp $(debug_srcdepends)
Cgemm_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Cgemm_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemv_debug_dd_SOURCES=../common/Cgemv.debug.cpp $(debug_srcdepends)
Cgemv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Cgemv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemv_debug_dd_ref_SOURCES=../common/Cgemv.debug.cpp $(debug_srcdepends)
Cgemv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Cgemv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgerc_debug_dd_SOURCES=../common/Cgerc.debug.cpp $(debug_srcdepends)
Cgerc_debug_dd_CXXFLAGS=$(debug_cxxflags)
Cgerc_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgerc_debug_dd_ref_SOURCES=../common/Cgerc.debug.cpp $(debug_srcdepends)
Cgerc_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Cgerc_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgeru_debug_dd_SOURCES=../common/Cgeru.debug.cpp $(debug_srcdepends)
Cgeru_debug_dd_CXXFLAGS=$(debug_cxxflags)
Cgeru_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgeru_debug_dd_ref_SOURCES=../common/Cgeru.debug.cpp $(debug_srcdepends)
Cgeru_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Cgeru_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chbmv_debug_dd_SOURCES=../common/Chbmv.debug.cpp $(debug_srcdepends)
Chbmv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Chbmv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chbmv_debug_dd_ref_SOURCES=../common/Chbmv.debug.cpp $(debug_srcdepends)
Chbmv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Chbmv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemm_debug_dd_SOURCES=../common/Chemm.debug.cpp $(debug_srcdepends)
Chemm_debug_dd_CXXFLAGS=$(debug_cxxflags)
Chemm_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemm_debug_dd_ref_SOURCES=../common/Chemm.debug.cpp $(debug_srcdepends)
Chemm_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Chemm_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemv_debug_dd_SOURCES=../common/Chemv.debug.cpp $(debug_srcdepends)
Chemv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Chemv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemv_debug_dd_ref_SOURCES=../common/Chemv.debug.cpp $(debug_srcdepends)
Chemv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Chemv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher_debug_dd_SOURCES=../common/Cher.debug.cpp $(debug_srcdepends)
Cher_debug_dd_CXXFLAGS=$(debug_cxxflags)
Cher_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher_debug_dd_ref_SOURCES=../common/Cher.debug.cpp $(debug_srcdepends)
Cher_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Cher_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2_debug_dd_SOURCES=../common/Cher2.debug.cpp $(debug_srcdepends)
Cher2_debug_dd_CXXFLAGS=$(debug_cxxflags)
Cher2_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2_debug_dd_ref_SOURCES=../common/Cher2.debug.cpp $(debug_srcdepends)
Cher2_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Cher2_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2k_debug_dd_SOURCES=../common/Cher2k.debug.cpp $(debug_srcdepends)
Cher2k_debug_dd_CXXFLAGS=$(debug_cxxflags)
Cher2k_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2k_debug_dd_ref_SOURCES=../common/Cher2k.debug.cpp $(debug_srcdepends)
Cher2k_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Cher2k_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cherk_debug_dd_SOURCES=../common/Cherk.debug.cpp $(debug_srcdepends)
Cherk_debug_dd_CXXFLAGS=$(debug_cxxflags)
Cherk_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cherk_debug_dd_ref_SOURCES=../common/Cherk.debug.cpp $(debug_srcdepends)
Cherk_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Cherk_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpmv_debug_dd_SOURCES=../common/Chpmv.debug.cpp $(debug_srcdepends)
Chpmv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Chpmv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpmv_debug_dd_ref_SOURCES=../common/Chpmv.debug.cpp $(debug_srcdepends)
Chpmv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Chpmv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr_debug_dd_SOURCES=../common/Chpr.debug.cpp $(debug_srcdepends)
Chpr_debug_dd_CXXFLAGS=$(debug_cxxflags)
Chpr_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr_debug_dd_ref_SOURCES=../common/Chpr.debug.cpp $(debug_srcdepends)
Chpr_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Chpr_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr2_debug_dd_SOURCES=../common/Chpr2.debug.cpp $(debug_srcdepends)
Chpr2_debug_dd_CXXFLAGS=$(debug_cxxflags)
Chpr2_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr2_debug_dd_ref_SOURCES=../common/Chpr2.debug.cpp $(debug_srcdepends)
Chpr2_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Chpr2_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Crotg_debug_dd_SOURCES=../common/Crotg.debug.cpp $(debug_srcdepends)
Crotg_debug_dd_CXXFLAGS=$(debug_cxxflags)
Crotg_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Crotg_debug_dd_ref_SOURCES=../common/Crotg.debug.cpp $(debug_srcdepends)
Crotg_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Crotg_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cscal_debug_dd_SOURCES=../common/Cscal.debug.cpp $(debug_srcdepends)
Cscal_debug_dd_CXXFLAGS=$(debug_cxxflags)
Cscal_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cscal_debug_dd_ref_SOURCES=../common/Cscal.debug.cpp $(debug_srcdepends)
Cscal_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Cscal_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cswap_debug_dd_SOURCES=../common/Cswap.debug.cpp $(debug_srcdepends)
Cswap_debug_dd_CXXFLAGS=$(debug_cxxflags)
Cswap_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cswap_debug_dd_ref_SOURCES=../common/Cswap.debug.cpp $(debug_srcdepends)
Cswap_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Cswap_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csymm_debug_dd_SOURCES=../common/Csymm.debug.cpp $(debug_srcdepends)
Csymm_debug_dd_CXXFLAGS=$(debug_cxxflags)
Csymm_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csymm_debug_dd_ref_SOURCES=../common/Csymm.debug.cpp $(debug_srcdepends)
Csymm_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Csymm_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyr2k_debug_dd_SOURCES=../common/Csyr2k.debug.cpp $(debug_srcdepends)
Csyr2k_debug_dd_CXXFLAGS=$(debug_cxxflags)
Csyr2k_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyr2k_debug_dd_ref_SOURCES=../common/Csyr2k.debug.cpp $(debug_srcdepends)
Csyr2k_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Csyr2k_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyrk_debug_dd_SOURCES=../common/Csyrk.debug.cpp $(debug_srcdepends)
Csyrk_debug_dd_CXXFLAGS=$(debug_cxxflags)
Csyrk_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyrk_debug_dd_ref_SOURCES=../common/Csyrk.debug.cpp $(debug_srcdepends)
Csyrk_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Csyrk_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbmv_debug_dd_SOURCES=../common/Ctbmv.debug.cpp $(debug_srcdepends)
Ctbmv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Ctbmv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbmv_debug_dd_ref_SOURCES=../common/Ctbmv.debug.cpp $(debug_srcdepends)
Ctbmv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Ctbmv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbsv_debug_dd_SOURCES=../common/Ctbsv.debug.cpp $(debug_srcdepends)
Ctbsv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Ctbsv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbsv_debug_dd_ref_SOURCES=../common/Ctbsv.debug.cpp $(debug_srcdepends)
Ctbsv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Ctbsv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpmv_debug_dd_SOURCES=../common/Ctpmv.debug.cpp $(debug_srcdepends)
Ctpmv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Ctpmv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpmv_debug_dd_ref_SOURCES=../common/Ctpmv.debug.cpp $(debug_srcdepends)
Ctpmv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Ctpmv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpsv_debug_dd_SOURCES=../common/Ctpsv.debug.cpp $(debug_srcdepends)
Ctpsv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Ctpsv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpsv_debug_dd_ref_SOURCES=../common/Ctpsv.debug.cpp $(debug_srcdepends)
Ctpsv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Ctpsv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmm_debug_dd_SOURCES=../common/Ctrmm.debug.cpp $(debug_srcdepends)
Ctrmm_debug_dd_CXXFLAGS=$(debug_cxxflags)
Ctrmm_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmm_debug_dd_ref_SOURCES=../common/Ctrmm.debug.cpp $(debug_srcdepends)
Ctrmm_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Ctrmm_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmv_debug_dd_SOURCES=../common/Ctrmv.debug.cpp $(debug_srcdepends)
Ctrmv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Ctrmv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmv_debug_dd_ref_SOURCES=../common/Ctrmv.debug.cpp $(debug_srcdepends)
Ctrmv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Ctrmv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsm_debug_dd_SOURCES=../common/Ctrsm.debug.cpp $(debug_srcdepends)
Ctrsm_debug_dd_CXXFLAGS=$(debug_cxxflags)
Ctrsm_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsm_debug_dd_ref_SOURCES=../common/Ctrsm.debug.cpp $(debug_srcdepends)
Ctrsm_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Ctrsm_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsv_debug_dd_SOURCES=../common/Ctrsv.debug.cpp $(debug_srcdepends)
Ctrsv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Ctrsv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsv_debug_dd_ref_SOURCES=../common/Ctrsv.debug.cpp $(debug_srcdepends)
Ctrsv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Ctrsv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_debug_dd_ref_SOURCES=../common/RCabs1.debug.cpp $(debug_srcdepends)
RCabs1_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
RCabs1_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_debug_dd_SOURCES=../common/RCabs1.debug.cpp $(debug_srcdepends)
RCabs1_debug_dd_CXXFLAGS=$(debug_cxxflags)
RCabs1_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCasum_debug_dd_ref_SOURCES=../common/RCasum.debug.cpp $(debug_srcdepends)
RCasum_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
RCasum_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCasum_debug_dd_SOURCES=../common/RCasum.debug.cpp $(debug_srcdepends)
RCasum_debug_dd_CXXFLAGS=$(debug_cxxflags)
RCasum_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCnrm2_debug_dd_ref_SOURCES=../common/RCnrm2.debug.cpp $(debug_srcdepends)
RCnrm2_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
RCnrm2_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCnrm2_debug_dd_SOURCES=../common/RCnrm2.debug.cpp $(debug_srcdepends)
RCnrm2_debug_dd_CXXFLAGS=$(debug_cxxflags)
RCnrm2_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rasum_debug_dd_ref_SOURCES=../common/Rasum.debug.cpp $(debug_srcdepends)
Rasum_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rasum_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rasum_debug_dd_SOURCES=../common/Rasum.debug.cpp $(debug_srcdepends)
Rasum_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rasum_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Raxpy_debug_dd_ref_SOURCES=../common/Raxpy.debug.cpp $(debug_srcdepends)
Raxpy_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Raxpy_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Raxpy_debug_dd_SOURCES=../common/Raxpy.debug.cpp $(debug_srcdepends)
Raxpy_debug_dd_CXXFLAGS=$(debug_cxxflags)
Raxpy_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rcopy_debug_dd_ref_SOURCES=../common/Rcopy.debug.cpp $(debug_srcdepends)
Rcopy_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rcopy_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rcopy_debug_dd_SOURCES=../common/Rcopy.debug.cpp $(debug_srcdepends)
Rcopy_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rcopy_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rdot_debug_dd_ref_SOURCES=../common/Rdot.debug.cpp $(debug_srcdepends)
Rdot_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rdot_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rdot_debug_dd_SOURCES=../common/Rdot.debug.cpp $(debug_srcdepends)
Rdot_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rdot_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgbmv_debug_dd_ref_SOURCES=../common/Rgbmv.debug.cpp $(debug_srcdepends)
Rgbmv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rgbmv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgbmv_debug_dd_SOURCES=../common/Rgbmv.debug.cpp $(debug_srcdepends)
Rgbmv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rgbmv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemm_debug_dd_ref_SOURCES=../common/Rgemm.debug.cpp $(debug_srcdepends)
Rgemm_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rgemm_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemm_debug_dd_SOURCES=../common/Rgemm.debug.cpp $(debug_srcdepends)
Rgemm_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rgemm_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemv_debug_dd_ref_SOURCES=../common/Rgemv.debug.cpp $(debug_srcdepends)
Rgemv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rgemv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemv_debug_dd_SOURCES=../common/Rgemv.debug.cpp $(debug_srcdepends)
Rgemv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rgemv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rger_debug_dd_ref_SOURCES=../common/Rger.debug.cpp $(debug_srcdepends)
Rger_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rger_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rger_debug_dd_SOURCES=../common/Rger.debug.cpp $(debug_srcdepends)
Rger_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rger_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rnrm2_debug_dd_ref_SOURCES=../common/Rnrm2.debug.cpp $(debug_srcdepends)
Rnrm2_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rnrm2_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rnrm2_debug_dd_SOURCES=../common/Rnrm2.debug.cpp $(debug_srcdepends)
Rnrm2_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rnrm2_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrot_debug_dd_ref_SOURCES=../common/Rrot.debug.cpp $(debug_srcdepends)
Rrot_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rrot_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrot_debug_dd_SOURCES=../common/Rrot.debug.cpp $(debug_srcdepends)
Rrot_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rrot_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotg_debug_dd_ref_SOURCES=../common/Rrotg.debug.cpp $(debug_srcdepends)
Rrotg_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rrotg_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotg_debug_dd_SOURCES=../common/Rrotg.debug.cpp $(debug_srcdepends)
Rrotg_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rrotg_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotm_debug_dd_ref_SOURCES=../common/Rrotm.debug.cpp $(debug_srcdepends)
Rrotm_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rrotm_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotm_debug_dd_SOURCES=../common/Rrotm.debug.cpp $(debug_srcdepends)
Rrotm_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rrotm_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsbmv_debug_dd_ref_SOURCES=../common/Rsbmv.debug.cpp $(debug_srcdepends)
Rsbmv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rsbmv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsbmv_debug_dd_SOURCES=../common/Rsbmv.debug.cpp $(debug_srcdepends)
Rsbmv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rsbmv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rscal_debug_dd_ref_SOURCES=../common/Rscal.debug.cpp $(debug_srcdepends)
Rscal_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rscal_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rscal_debug_dd_SOURCES=../common/Rscal.debug.cpp $(debug_srcdepends)
Rscal_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rscal_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspmv_debug_dd_ref_SOURCES=../common/Rspmv.debug.cpp $(debug_srcdepends)
Rspmv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rspmv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspmv_debug_dd_SOURCES=../common/Rspmv.debug.cpp $(debug_srcdepends)
Rspmv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rspmv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr_debug_dd_ref_SOURCES=../common/Rspr.debug.cpp $(debug_srcdepends)
Rspr_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rspr_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr_debug_dd_SOURCES=../common/Rspr.debug.cpp $(debug_srcdepends)
Rspr_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rspr_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr2_debug_dd_ref_SOURCES=../common/Rspr2.debug.cpp $(debug_srcdepends)
Rspr2_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rspr2_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr2_debug_dd_SOURCES=../common/Rspr2.debug.cpp $(debug_srcdepends)
Rspr2_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rspr2_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rswap_debug_dd_ref_SOURCES=../common/Rswap.debug.cpp $(debug_srcdepends)
Rswap_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rswap_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rswap_debug_dd_SOURCES=../common/Rswap.debug.cpp $(debug_srcdepends)
Rswap_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rswap_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymm_debug_dd_ref_SOURCES=../common/Rsymm.debug.cpp $(debug_srcdepends)
Rsymm_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rsymm_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymm_debug_dd_SOURCES=../common/Rsymm.debug.cpp $(debug_srcdepends)
Rsymm_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rsymm_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymv_debug_dd_ref_SOURCES=../common/Rsymv.debug.cpp $(debug_srcdepends)
Rsymv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rsymv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymv_debug_dd_SOURCES=../common/Rsymv.debug.cpp $(debug_srcdepends)
Rsymv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rsymv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr_debug_dd_ref_SOURCES=../common/Rsyr.debug.cpp $(debug_srcdepends)
Rsyr_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr_debug_dd_SOURCES=../common/Rsyr.debug.cpp $(debug_srcdepends)
Rsyr_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rsyr_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2_debug_dd_ref_SOURCES=../common/Rsyr2.debug.cpp $(debug_srcdepends)
Rsyr2_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr2_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2_debug_dd_SOURCES=../common/Rsyr2.debug.cpp $(debug_srcdepends)
Rsyr2_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rsyr2_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2k_debug_dd_ref_SOURCES=../common/Rsyr2k.debug.cpp $(debug_srcdepends)
Rsyr2k_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr2k_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2k_debug_dd_SOURCES=../common/Rsyr2k.debug.cpp $(debug_srcdepends)
Rsyr2k_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rsyr2k_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyrk_debug_dd_ref_SOURCES=../common/Rsyrk.debug.cpp $(debug_srcdepends)
Rsyrk_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rsyrk_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyrk_debug_dd_SOURCES=../common/Rsyrk.debug.cpp $(debug_srcdepends)
Rsyrk_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rsyrk_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbmv_debug_dd_ref_SOURCES=../common/Rtbmv.debug.cpp $(debug_srcdepends)
Rtbmv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rtbmv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbmv_debug_dd_SOURCES=../common/Rtbmv.debug.cpp $(debug_srcdepends)
Rtbmv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rtbmv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbsv_debug_dd_ref_SOURCES=../common/Rtbsv.debug.cpp $(debug_srcdepends)
Rtbsv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rtbsv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbsv_debug_dd_SOURCES=../common/Rtbsv.debug.cpp $(debug_srcdepends)
Rtbsv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rtbsv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpmv_debug_dd_ref_SOURCES=../common/Rtpmv.debug.cpp $(debug_srcdepends)
Rtpmv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rtpmv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpmv_debug_dd_SOURCES=../common/Rtpmv.debug.cpp $(debug_srcdepends)
Rtpmv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rtpmv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpsv_debug_dd_ref_SOURCES=../common/Rtpsv.debug.cpp $(debug_srcdepends)
Rtpsv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rtpsv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpsv_debug_dd_SOURCES=../common/Rtpsv.debug.cpp $(debug_srcdepends)
Rtpsv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rtpsv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmm_debug_dd_ref_SOURCES=../common/Rtrmm.debug.cpp $(debug_srcdepends)
Rtrmm_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rtrmm_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmm_debug_dd_SOURCES=../common/Rtrmm.debug.cpp $(debug_srcdepends)
Rtrmm_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rtrmm_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmv_debug_dd_ref_SOURCES=../common/Rtrmv.debug.cpp $(debug_srcdepends)
Rtrmv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rtrmv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmv_debug_dd_SOURCES=../common/Rtrmv.debug.cpp $(debug_srcdepends)
Rtrmv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rtrmv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsm_debug_dd_ref_SOURCES=../common/Rtrsm.debug.cpp $(debug_srcdepends)
Rtrsm_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rtrsm_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsm_debug_dd_SOURCES=../common/Rtrsm.debug.cpp $(debug_srcdepends)
Rtrsm_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rtrsm_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsv_debug_dd_ref_SOURCES=../common/Rtrsv.debug.cpp $(debug_srcdepends)
Rtrsv_debug_dd_ref_CXXFLAGS=$(debug_cxxflags)
Rtrsv_debug_dd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsv_debug_dd_SOURCES=../common/Rtrsv.debug.cpp $(debug_srcdepends)
Rtrsv_debug_dd_CXXFLAGS=$(debug_cxxflags)
Rtrsv_debug_dd_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

