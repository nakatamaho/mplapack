#CTBSV, CTPSV might have bugs or intrinsically inaccurate
#Rnrm2 should be checked more
check_PROGRAMS = $(mblas_mpfr_test_PROGRAMS) 
TESTS = $(check_PROGRAMS)
mblas_mpfr_test_PROGRAMS = \
arithmetic.debug_mpfr \
complex.debug_mpfr \
iCamax.debug.mpfr_ref iCamax.debug.mpfr   \
iRamax.debug.mpfr_ref iRamax.debug.mpfr   \
Mlsame.debug.mpfr_ref Mlsame.debug.mpfr   \
Mxerbla.debug.mpfr_ref Mxerbla.debug.mpfr \
CRrot.debug.mpfr_ref  CRrot.debug.mpfr    \
CRscal.debug.mpfr_ref CRscal.debug.mpfr   \
Caxpy.debug.mpfr_ref  Caxpy.debug.mpfr    \
Ccopy.debug.mpfr_ref  Ccopy.debug.mpfr    \
Cdotc.debug.mpfr_ref  Cdotc.debug.mpfr	\
Cdotu.debug.mpfr_ref  Cdotu.debug.mpfr	\
Cgbmv.debug.mpfr_ref  Cgbmv.debug.mpfr	\
Cgemm.debug.mpfr_ref  Cgemm.debug.mpfr	\
Cgemv.debug.mpfr_ref  Cgemv.debug.mpfr	\
Cgerc.debug.mpfr_ref  Cgerc.debug.mpfr	\
Cgeru.debug.mpfr_ref  Cgeru.debug.mpfr    \
Chbmv.debug.mpfr_ref  Chbmv.debug.mpfr	\
Chemm.debug.mpfr_ref  Chemm.debug.mpfr	\
Chemv.debug.mpfr_ref  Chemv.debug.mpfr	\
Cher.debug.mpfr_ref   Cher.debug.mpfr 	\
Cher2.debug.mpfr_ref  Cher2.debug.mpfr	\
Cher2k.debug.mpfr_ref Cher2k.debug.mpfr	\
Cherk.debug.mpfr_ref  Cherk.debug.mpfr	\
Chpmv.debug.mpfr_ref  Chpmv.debug.mpfr	\
Chpr.debug.mpfr_ref   Chpr.debug.mpfr	\
Chpr2.debug.mpfr_ref  Chpr2.debug.mpfr	\
Crotg.debug.mpfr_ref  Crotg.debug.mpfr	\
Cscal.debug.mpfr_ref  Cscal.debug.mpfr	\
Cswap.debug.mpfr_ref  Cswap.debug.mpfr	\
Csymm.debug.mpfr_ref  Csymm.debug.mpfr    \
Csyr2k.debug.mpfr_ref Csyr2k.debug.mpfr	\
Csyrk.debug.mpfr_ref  Csyrk.debug.mpfr	\
Ctbmv.debug.mpfr_ref  Ctbmv.debug.mpfr	\
Ctbsv.debug.mpfr_ref  Ctbsv.debug.mpfr	\
Ctpmv.debug.mpfr_ref  Ctpmv.debug.mpfr	\
Ctpsv.debug.mpfr_ref  Ctpsv.debug.mpfr	\
Ctrmm.debug.mpfr_ref  Ctrmm.debug.mpfr	\
Ctrmv.debug.mpfr_ref  Ctrmv.debug.mpfr	\
Ctrsm.debug.mpfr_ref  Ctrsm.debug.mpfr	\
Ctrsv.debug.mpfr_ref  Ctrsv.debug.mpfr	\
RCabs1.debug.mpfr_ref RCabs1.debug.mpfr	\
RCasum.debug.mpfr_ref RCasum.debug.mpfr	\
RCnrm2.debug.mpfr_ref RCnrm2.debug.mpfr	\
Rasum.debug.mpfr_ref  Rasum.debug.mpfr    \
Raxpy.debug.mpfr_ref  Raxpy.debug.mpfr	\
Rcopy.debug.mpfr_ref  Rcopy.debug.mpfr	\
Rdot.debug.mpfr_ref   Rdot.debug.mpfr	\
Rgbmv.debug.mpfr_ref  Rgbmv.debug.mpfr	\
Rgemm.debug.mpfr_ref  Rgemm.debug.mpfr	\
Rgemv.debug.mpfr_ref  Rgemv.debug.mpfr	\
Rger.debug.mpfr_ref   Rger.debug.mpfr	\
Rnrm2.debug.mpfr_ref  Rnrm2.debug.mpfr	\
Rrot.debug.mpfr_ref   Rrot.debug.mpfr	\
Rrotg.debug.mpfr_ref  Rrotg.debug.mpfr	\
Rrotm.debug.mpfr_ref  Rrotm.debug.mpfr	\
Rsbmv.debug.mpfr_ref  Rsbmv.debug.mpfr	\
Rscal.debug.mpfr_ref  Rscal.debug.mpfr	\
Rspmv.debug.mpfr_ref  Rspmv.debug.mpfr    \
Rspr.debug.mpfr_ref   Rspr.debug.mpfr	\
Rspr2.debug.mpfr_ref  Rspr2.debug.mpfr	\
Rswap.debug.mpfr_ref  Rswap.debug.mpfr	\
Rsymm.debug.mpfr_ref  Rsymm.debug.mpfr	\
Rsymv.debug.mpfr_ref  Rsymv.debug.mpfr	\
Rsyr.debug.mpfr_ref   Rsyr.debug.mpfr	\
Rsyr2.debug.mpfr_ref  Rsyr2.debug.mpfr	\
Rsyr2k.debug.mpfr_ref Rsyr2k.debug.mpfr	\
Rsyrk.debug.mpfr_ref  Rsyrk.debug.mpfr	\
Rtbmv.debug.mpfr_ref  Rtbmv.debug.mpfr	\
Rtbsv.debug.mpfr_ref  Rtbsv.debug.mpfr	\
Rtpmv.debug.mpfr_ref  Rtpmv.debug.mpfr	\
Rtpsv.debug.mpfr_ref  Rtpsv.debug.mpfr	\
Rtrmm.debug.mpfr_ref  Rtrmm.debug.mpfr    \
Rtrmv.debug.mpfr_ref  Rtrmv.debug.mpfr	\
Rtrsm.debug.mpfr_ref  Rtrsm.debug.mpfr	\
Rtrsv.debug.mpfr_ref  Rtrsv.debug.mpfr

mblas_mpfr_testdir = $(prefix)/share/mpack/testing/mpfr
ref_libdepends = -L$(top_builddir)/mblas/reference -lmblas_mpfr_ref
opt_libdepends = -L$(top_builddir)/mblas/optimized/mpfr -lmblas_mpfr
mplibs=-L$(top_builddir)/mpfrc++ -lmpfrcxx -L$(top_srcdir)/$(MPC_LIBDIR) -L$(top_srcdir)/$(MPFR_LIBDIR) -L$(top_srcdir)/$(GMP_LIBDIR) -lmpfr -lmpc -lgmp
refblas=-L$(top_builddir)/external/i/LAPACK/lib -lblas $(FCLIBS)

debug_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(top_srcdir)/$(GMP_INCLUDEDIR) -I$(top_srcdir)/$(MPFR_INCLUDEDIR) -I$(top_srcdir)/$(MPC_INCLUDEDIR) -D___MPACK_BUILD_WITH_MPFR___

if !IS_MACOSX
debug_srcdepends = ../common/xerbla.cpp ../common/mpack.debug.cpp ../common/Mxerbla.override.cpp
endif

if IS_MACOSX
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPACK_BUILD_WITH_MPFR___
libmxerbla_override_la_SOURCES =  ../common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS= -Wl,-flat_namespace,-undefined,dynamic_lookup
debug_srcdepends = ../common/xerbla.cpp ../common/mpack.debug.cpp
mplibs+=-Wl,-flat_namespace,-undefined,dynamic_lookup
TESTS_ENVIRONMENT=DYLD_INSERT_LIBRARIES=.libs/libmxerbla_override.0.dylib
endif

##
arithmetic_debug_mpfr_SOURCES=../common/arithmetic.debug.cpp $(debug_srcdepends)
arithmetic_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
arithmetic_debug_mpfr_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

complex_debug_mpfr_SOURCES=../common/complex.debug.cpp $(debug_srcdepends)
complex_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
complex_debug_mpfr_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_debug_mpfr_ref_SOURCES=../common/iCamax.debug.cpp $(debug_srcdepends)
iCamax_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
iCamax_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_debug_mpfr_SOURCES=../common/iCamax.debug.cpp $(debug_srcdepends)
iCamax_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
iCamax_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iRamax_debug_mpfr_ref_SOURCES=../common/iRamax.debug.cpp $(debug_srcdepends)
iRamax_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
iRamax_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iRamax_debug_mpfr_SOURCES=../common/iRamax.debug.cpp $(debug_srcdepends)
iRamax_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
iRamax_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mlsame_debug_mpfr_ref_SOURCES=../common/Mlsame.debug.cpp $(debug_srcdepends)
Mlsame_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Mlsame_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mlsame_debug_mpfr_SOURCES=../common/Mlsame.debug.cpp $(debug_srcdepends)
Mlsame_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Mlsame_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mxerbla_debug_mpfr_ref_SOURCES=../common/Mxerbla.debug.cpp $(debug_srcdepends)
Mxerbla_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Mxerbla_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mxerbla_debug_mpfr_SOURCES=../common/Mxerbla.debug.cpp $(debug_srcdepends)
Mxerbla_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Mxerbla_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRrot_debug_mpfr_ref_SOURCES=../common/CRrot.debug.cpp $(debug_srcdepends)
CRrot_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
CRrot_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRrot_debug_mpfr_SOURCES=../common/CRrot.debug.cpp $(debug_srcdepends)
CRrot_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
CRrot_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRscal_debug_mpfr_ref_SOURCES=../common/CRscal.debug.cpp $(debug_srcdepends)
CRscal_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
CRscal_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRscal_debug_mpfr_SOURCES=../common/CRscal.debug.cpp $(debug_srcdepends)
CRscal_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
CRscal_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_debug_mpfr_SOURCES=../common/Caxpy.debug.cpp $(debug_srcdepends)
Caxpy_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Caxpy_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_debug_mpfr_ref_SOURCES=../common/Caxpy.debug.cpp $(debug_srcdepends)
Caxpy_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Caxpy_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ccopy_debug_mpfr_SOURCES=../common/Ccopy.debug.cpp $(debug_srcdepends)
Ccopy_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Ccopy_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ccopy_debug_mpfr_ref_SOURCES=../common/Ccopy.debug.cpp $(debug_srcdepends)
Ccopy_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Ccopy_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotc_debug_mpfr_SOURCES=../common/Cdotc.debug.cpp $(debug_srcdepends)
Cdotc_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Cdotc_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotc_debug_mpfr_ref_SOURCES=../common/Cdotc.debug.cpp $(debug_srcdepends)
Cdotc_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Cdotc_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotu_debug_mpfr_SOURCES=../common/Cdotu.debug.cpp $(debug_srcdepends)
Cdotu_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Cdotu_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotu_debug_mpfr_ref_SOURCES=../common/Cdotu.debug.cpp $(debug_srcdepends)
Cdotu_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Cdotu_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgbmv_debug_mpfr_SOURCES=../common/Cgbmv.debug.cpp $(debug_srcdepends)
Cgbmv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Cgbmv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgbmv_debug_mpfr_ref_SOURCES=../common/Cgbmv.debug.cpp $(debug_srcdepends)
Cgbmv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Cgbmv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemm_debug_mpfr_SOURCES=../common/Cgemm.debug.cpp $(debug_srcdepends)
Cgemm_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Cgemm_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemm_debug_mpfr_ref_SOURCES=../common/Cgemm.debug.cpp $(debug_srcdepends)
Cgemm_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Cgemm_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemv_debug_mpfr_SOURCES=../common/Cgemv.debug.cpp $(debug_srcdepends)
Cgemv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Cgemv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemv_debug_mpfr_ref_SOURCES=../common/Cgemv.debug.cpp $(debug_srcdepends)
Cgemv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Cgemv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgerc_debug_mpfr_SOURCES=../common/Cgerc.debug.cpp $(debug_srcdepends)
Cgerc_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Cgerc_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgerc_debug_mpfr_ref_SOURCES=../common/Cgerc.debug.cpp $(debug_srcdepends)
Cgerc_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Cgerc_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgeru_debug_mpfr_SOURCES=../common/Cgeru.debug.cpp $(debug_srcdepends)
Cgeru_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Cgeru_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgeru_debug_mpfr_ref_SOURCES=../common/Cgeru.debug.cpp $(debug_srcdepends)
Cgeru_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Cgeru_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chbmv_debug_mpfr_SOURCES=../common/Chbmv.debug.cpp $(debug_srcdepends)
Chbmv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Chbmv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chbmv_debug_mpfr_ref_SOURCES=../common/Chbmv.debug.cpp $(debug_srcdepends)
Chbmv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Chbmv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemm_debug_mpfr_SOURCES=../common/Chemm.debug.cpp $(debug_srcdepends)
Chemm_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Chemm_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemm_debug_mpfr_ref_SOURCES=../common/Chemm.debug.cpp $(debug_srcdepends)
Chemm_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Chemm_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemv_debug_mpfr_SOURCES=../common/Chemv.debug.cpp $(debug_srcdepends)
Chemv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Chemv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemv_debug_mpfr_ref_SOURCES=../common/Chemv.debug.cpp $(debug_srcdepends)
Chemv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Chemv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher_debug_mpfr_SOURCES=../common/Cher.debug.cpp $(debug_srcdepends)
Cher_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Cher_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher_debug_mpfr_ref_SOURCES=../common/Cher.debug.cpp $(debug_srcdepends)
Cher_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Cher_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2_debug_mpfr_SOURCES=../common/Cher2.debug.cpp $(debug_srcdepends)
Cher2_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Cher2_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2_debug_mpfr_ref_SOURCES=../common/Cher2.debug.cpp $(debug_srcdepends)
Cher2_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Cher2_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2k_debug_mpfr_SOURCES=../common/Cher2k.debug.cpp $(debug_srcdepends)
Cher2k_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Cher2k_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2k_debug_mpfr_ref_SOURCES=../common/Cher2k.debug.cpp $(debug_srcdepends)
Cher2k_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Cher2k_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cherk_debug_mpfr_SOURCES=../common/Cherk.debug.cpp $(debug_srcdepends)
Cherk_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Cherk_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cherk_debug_mpfr_ref_SOURCES=../common/Cherk.debug.cpp $(debug_srcdepends)
Cherk_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Cherk_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpmv_debug_mpfr_SOURCES=../common/Chpmv.debug.cpp $(debug_srcdepends)
Chpmv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Chpmv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpmv_debug_mpfr_ref_SOURCES=../common/Chpmv.debug.cpp $(debug_srcdepends)
Chpmv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Chpmv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr_debug_mpfr_SOURCES=../common/Chpr.debug.cpp $(debug_srcdepends)
Chpr_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Chpr_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr_debug_mpfr_ref_SOURCES=../common/Chpr.debug.cpp $(debug_srcdepends)
Chpr_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Chpr_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr2_debug_mpfr_SOURCES=../common/Chpr2.debug.cpp $(debug_srcdepends)
Chpr2_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Chpr2_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr2_debug_mpfr_ref_SOURCES=../common/Chpr2.debug.cpp $(debug_srcdepends)
Chpr2_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Chpr2_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Crotg_debug_mpfr_SOURCES=../common/Crotg.debug.cpp $(debug_srcdepends)
Crotg_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Crotg_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Crotg_debug_mpfr_ref_SOURCES=../common/Crotg.debug.cpp $(debug_srcdepends)
Crotg_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Crotg_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cscal_debug_mpfr_SOURCES=../common/Cscal.debug.cpp $(debug_srcdepends)
Cscal_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Cscal_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cscal_debug_mpfr_ref_SOURCES=../common/Cscal.debug.cpp $(debug_srcdepends)
Cscal_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Cscal_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cswap_debug_mpfr_SOURCES=../common/Cswap.debug.cpp $(debug_srcdepends)
Cswap_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Cswap_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cswap_debug_mpfr_ref_SOURCES=../common/Cswap.debug.cpp $(debug_srcdepends)
Cswap_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Cswap_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csymm_debug_mpfr_SOURCES=../common/Csymm.debug.cpp $(debug_srcdepends)
Csymm_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Csymm_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csymm_debug_mpfr_ref_SOURCES=../common/Csymm.debug.cpp $(debug_srcdepends)
Csymm_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Csymm_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyr2k_debug_mpfr_SOURCES=../common/Csyr2k.debug.cpp $(debug_srcdepends)
Csyr2k_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Csyr2k_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyr2k_debug_mpfr_ref_SOURCES=../common/Csyr2k.debug.cpp $(debug_srcdepends)
Csyr2k_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Csyr2k_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyrk_debug_mpfr_SOURCES=../common/Csyrk.debug.cpp $(debug_srcdepends)
Csyrk_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Csyrk_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyrk_debug_mpfr_ref_SOURCES=../common/Csyrk.debug.cpp $(debug_srcdepends)
Csyrk_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Csyrk_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbmv_debug_mpfr_SOURCES=../common/Ctbmv.debug.cpp $(debug_srcdepends)
Ctbmv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Ctbmv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbmv_debug_mpfr_ref_SOURCES=../common/Ctbmv.debug.cpp $(debug_srcdepends)
Ctbmv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Ctbmv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbsv_debug_mpfr_SOURCES=../common/Ctbsv.debug.cpp $(debug_srcdepends)
Ctbsv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Ctbsv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbsv_debug_mpfr_ref_SOURCES=../common/Ctbsv.debug.cpp $(debug_srcdepends)
Ctbsv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Ctbsv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpmv_debug_mpfr_SOURCES=../common/Ctpmv.debug.cpp $(debug_srcdepends)
Ctpmv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Ctpmv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpmv_debug_mpfr_ref_SOURCES=../common/Ctpmv.debug.cpp $(debug_srcdepends)
Ctpmv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Ctpmv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpsv_debug_mpfr_SOURCES=../common/Ctpsv.debug.cpp $(debug_srcdepends)
Ctpsv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Ctpsv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpsv_debug_mpfr_ref_SOURCES=../common/Ctpsv.debug.cpp $(debug_srcdepends)
Ctpsv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Ctpsv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmm_debug_mpfr_SOURCES=../common/Ctrmm.debug.cpp $(debug_srcdepends)
Ctrmm_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Ctrmm_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmm_debug_mpfr_ref_SOURCES=../common/Ctrmm.debug.cpp $(debug_srcdepends)
Ctrmm_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Ctrmm_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmv_debug_mpfr_SOURCES=../common/Ctrmv.debug.cpp $(debug_srcdepends)
Ctrmv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Ctrmv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmv_debug_mpfr_ref_SOURCES=../common/Ctrmv.debug.cpp $(debug_srcdepends)
Ctrmv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Ctrmv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsm_debug_mpfr_SOURCES=../common/Ctrsm.debug.cpp $(debug_srcdepends)
Ctrsm_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Ctrsm_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsm_debug_mpfr_ref_SOURCES=../common/Ctrsm.debug.cpp $(debug_srcdepends)
Ctrsm_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Ctrsm_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsv_debug_mpfr_SOURCES=../common/Ctrsv.debug.cpp $(debug_srcdepends)
Ctrsv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Ctrsv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsv_debug_mpfr_ref_SOURCES=../common/Ctrsv.debug.cpp $(debug_srcdepends)
Ctrsv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Ctrsv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_debug_mpfr_ref_SOURCES=../common/RCabs1.debug.cpp $(debug_srcdepends)
RCabs1_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
RCabs1_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_debug_mpfr_SOURCES=../common/RCabs1.debug.cpp $(debug_srcdepends)
RCabs1_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
RCabs1_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCasum_debug_mpfr_ref_SOURCES=../common/RCasum.debug.cpp $(debug_srcdepends)
RCasum_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
RCasum_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCasum_debug_mpfr_SOURCES=../common/RCasum.debug.cpp $(debug_srcdepends)
RCasum_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
RCasum_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCnrm2_debug_mpfr_ref_SOURCES=../common/RCnrm2.debug.cpp $(debug_srcdepends)
RCnrm2_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
RCnrm2_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCnrm2_debug_mpfr_SOURCES=../common/RCnrm2.debug.cpp $(debug_srcdepends)
RCnrm2_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
RCnrm2_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rasum_debug_mpfr_ref_SOURCES=../common/Rasum.debug.cpp $(debug_srcdepends)
Rasum_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rasum_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rasum_debug_mpfr_SOURCES=../common/Rasum.debug.cpp $(debug_srcdepends)
Rasum_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rasum_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Raxpy_debug_mpfr_ref_SOURCES=../common/Raxpy.debug.cpp $(debug_srcdepends)
Raxpy_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Raxpy_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Raxpy_debug_mpfr_SOURCES=../common/Raxpy.debug.cpp $(debug_srcdepends)
Raxpy_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Raxpy_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rcopy_debug_mpfr_ref_SOURCES=../common/Rcopy.debug.cpp $(debug_srcdepends)
Rcopy_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rcopy_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rcopy_debug_mpfr_SOURCES=../common/Rcopy.debug.cpp $(debug_srcdepends)
Rcopy_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rcopy_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rdot_debug_mpfr_ref_SOURCES=../common/Rdot.debug.cpp $(debug_srcdepends)
Rdot_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rdot_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rdot_debug_mpfr_SOURCES=../common/Rdot.debug.cpp $(debug_srcdepends)
Rdot_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rdot_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgbmv_debug_mpfr_ref_SOURCES=../common/Rgbmv.debug.cpp $(debug_srcdepends)
Rgbmv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rgbmv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgbmv_debug_mpfr_SOURCES=../common/Rgbmv.debug.cpp $(debug_srcdepends)
Rgbmv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rgbmv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemm_debug_mpfr_ref_SOURCES=../common/Rgemm.debug.cpp $(debug_srcdepends)
Rgemm_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rgemm_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemm_debug_mpfr_SOURCES=../common/Rgemm.debug.cpp $(debug_srcdepends)
Rgemm_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rgemm_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemv_debug_mpfr_ref_SOURCES=../common/Rgemv.debug.cpp $(debug_srcdepends)
Rgemv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rgemv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemv_debug_mpfr_SOURCES=../common/Rgemv.debug.cpp $(debug_srcdepends)
Rgemv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rgemv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rger_debug_mpfr_ref_SOURCES=../common/Rger.debug.cpp $(debug_srcdepends)
Rger_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rger_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rger_debug_mpfr_SOURCES=../common/Rger.debug.cpp $(debug_srcdepends)
Rger_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rger_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rnrm2_debug_mpfr_ref_SOURCES=../common/Rnrm2.debug.cpp $(debug_srcdepends)
Rnrm2_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rnrm2_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rnrm2_debug_mpfr_SOURCES=../common/Rnrm2.debug.cpp $(debug_srcdepends)
Rnrm2_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rnrm2_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrot_debug_mpfr_ref_SOURCES=../common/Rrot.debug.cpp $(debug_srcdepends)
Rrot_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rrot_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrot_debug_mpfr_SOURCES=../common/Rrot.debug.cpp $(debug_srcdepends)
Rrot_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rrot_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotg_debug_mpfr_ref_SOURCES=../common/Rrotg.debug.cpp $(debug_srcdepends)
Rrotg_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rrotg_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotg_debug_mpfr_SOURCES=../common/Rrotg.debug.cpp $(debug_srcdepends)
Rrotg_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rrotg_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotm_debug_mpfr_ref_SOURCES=../common/Rrotm.debug.cpp $(debug_srcdepends)
Rrotm_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rrotm_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotm_debug_mpfr_SOURCES=../common/Rrotm.debug.cpp $(debug_srcdepends)
Rrotm_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rrotm_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsbmv_debug_mpfr_ref_SOURCES=../common/Rsbmv.debug.cpp $(debug_srcdepends)
Rsbmv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rsbmv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsbmv_debug_mpfr_SOURCES=../common/Rsbmv.debug.cpp $(debug_srcdepends)
Rsbmv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rsbmv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rscal_debug_mpfr_ref_SOURCES=../common/Rscal.debug.cpp $(debug_srcdepends)
Rscal_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rscal_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rscal_debug_mpfr_SOURCES=../common/Rscal.debug.cpp $(debug_srcdepends)
Rscal_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rscal_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspmv_debug_mpfr_ref_SOURCES=../common/Rspmv.debug.cpp $(debug_srcdepends)
Rspmv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rspmv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspmv_debug_mpfr_SOURCES=../common/Rspmv.debug.cpp $(debug_srcdepends)
Rspmv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rspmv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr_debug_mpfr_ref_SOURCES=../common/Rspr.debug.cpp $(debug_srcdepends)
Rspr_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rspr_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr_debug_mpfr_SOURCES=../common/Rspr.debug.cpp $(debug_srcdepends)
Rspr_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rspr_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr2_debug_mpfr_ref_SOURCES=../common/Rspr2.debug.cpp $(debug_srcdepends)
Rspr2_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rspr2_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr2_debug_mpfr_SOURCES=../common/Rspr2.debug.cpp $(debug_srcdepends)
Rspr2_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rspr2_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rswap_debug_mpfr_ref_SOURCES=../common/Rswap.debug.cpp $(debug_srcdepends)
Rswap_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rswap_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rswap_debug_mpfr_SOURCES=../common/Rswap.debug.cpp $(debug_srcdepends)
Rswap_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rswap_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymm_debug_mpfr_ref_SOURCES=../common/Rsymm.debug.cpp $(debug_srcdepends)
Rsymm_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rsymm_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymm_debug_mpfr_SOURCES=../common/Rsymm.debug.cpp $(debug_srcdepends)
Rsymm_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rsymm_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymv_debug_mpfr_ref_SOURCES=../common/Rsymv.debug.cpp $(debug_srcdepends)
Rsymv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rsymv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymv_debug_mpfr_SOURCES=../common/Rsymv.debug.cpp $(debug_srcdepends)
Rsymv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rsymv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr_debug_mpfr_ref_SOURCES=../common/Rsyr.debug.cpp $(debug_srcdepends)
Rsyr_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr_debug_mpfr_SOURCES=../common/Rsyr.debug.cpp $(debug_srcdepends)
Rsyr_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rsyr_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2_debug_mpfr_ref_SOURCES=../common/Rsyr2.debug.cpp $(debug_srcdepends)
Rsyr2_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr2_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2_debug_mpfr_SOURCES=../common/Rsyr2.debug.cpp $(debug_srcdepends)
Rsyr2_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rsyr2_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2k_debug_mpfr_ref_SOURCES=../common/Rsyr2k.debug.cpp $(debug_srcdepends)
Rsyr2k_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr2k_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2k_debug_mpfr_SOURCES=../common/Rsyr2k.debug.cpp $(debug_srcdepends)
Rsyr2k_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rsyr2k_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyrk_debug_mpfr_ref_SOURCES=../common/Rsyrk.debug.cpp $(debug_srcdepends)
Rsyrk_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rsyrk_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyrk_debug_mpfr_SOURCES=../common/Rsyrk.debug.cpp $(debug_srcdepends)
Rsyrk_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rsyrk_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbmv_debug_mpfr_ref_SOURCES=../common/Rtbmv.debug.cpp $(debug_srcdepends)
Rtbmv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rtbmv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbmv_debug_mpfr_SOURCES=../common/Rtbmv.debug.cpp $(debug_srcdepends)
Rtbmv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rtbmv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbsv_debug_mpfr_ref_SOURCES=../common/Rtbsv.debug.cpp $(debug_srcdepends)
Rtbsv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rtbsv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbsv_debug_mpfr_SOURCES=../common/Rtbsv.debug.cpp $(debug_srcdepends)
Rtbsv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rtbsv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpmv_debug_mpfr_ref_SOURCES=../common/Rtpmv.debug.cpp $(debug_srcdepends)
Rtpmv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rtpmv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpmv_debug_mpfr_SOURCES=../common/Rtpmv.debug.cpp $(debug_srcdepends)
Rtpmv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rtpmv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpsv_debug_mpfr_ref_SOURCES=../common/Rtpsv.debug.cpp $(debug_srcdepends)
Rtpsv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rtpsv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpsv_debug_mpfr_SOURCES=../common/Rtpsv.debug.cpp $(debug_srcdepends)
Rtpsv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rtpsv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmm_debug_mpfr_ref_SOURCES=../common/Rtrmm.debug.cpp $(debug_srcdepends)
Rtrmm_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rtrmm_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmm_debug_mpfr_SOURCES=../common/Rtrmm.debug.cpp $(debug_srcdepends)
Rtrmm_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rtrmm_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmv_debug_mpfr_ref_SOURCES=../common/Rtrmv.debug.cpp $(debug_srcdepends)
Rtrmv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rtrmv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmv_debug_mpfr_SOURCES=../common/Rtrmv.debug.cpp $(debug_srcdepends)
Rtrmv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rtrmv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsm_debug_mpfr_ref_SOURCES=../common/Rtrsm.debug.cpp $(debug_srcdepends)
Rtrsm_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rtrsm_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsm_debug_mpfr_SOURCES=../common/Rtrsm.debug.cpp $(debug_srcdepends)
Rtrsm_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rtrsm_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsv_debug_mpfr_ref_SOURCES=../common/Rtrsv.debug.cpp $(debug_srcdepends)
Rtrsv_debug_mpfr_ref_CXXFLAGS=$(debug_cxxflags)
Rtrsv_debug_mpfr_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsv_debug_mpfr_SOURCES=../common/Rtrsv.debug.cpp $(debug_srcdepends)
Rtrsv_debug_mpfr_CXXFLAGS=$(debug_cxxflags)
Rtrsv_debug_mpfr_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

