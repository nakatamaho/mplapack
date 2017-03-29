#CTBSV, CTPSV might have bugs or intrinsically inaccurate
#Rnrm2 should be checked more
check_PROGRAMS = $(mblas___float128_test_PROGRAMS) 
TESTS = $(check_PROGRAMS)
mblas___float128_test_PROGRAMS = \
arithmetic.debug___float128 \
complex.debug___float128 \
iCamax.debug.__float128_ref iCamax.debug.__float128   \
iRamax.debug.__float128_ref iRamax.debug.__float128   \
Mlsame.debug.__float128_ref Mlsame.debug.__float128   \
Mxerbla.debug.__float128_ref Mxerbla.debug.__float128 \
CRrot.debug.__float128_ref  CRrot.debug.__float128    \
CRscal.debug.__float128_ref CRscal.debug.__float128   \
Caxpy.debug.__float128_ref  Caxpy.debug.__float128    \
Ccopy.debug.__float128_ref  Ccopy.debug.__float128    \
Cdotc.debug.__float128_ref  Cdotc.debug.__float128	\
Cdotu.debug.__float128_ref  Cdotu.debug.__float128	\
Cgbmv.debug.__float128_ref  Cgbmv.debug.__float128	\
Cgemm.debug.__float128_ref  Cgemm.debug.__float128	\
Cgemv.debug.__float128_ref  Cgemv.debug.__float128	\
Cgerc.debug.__float128_ref  Cgerc.debug.__float128	\
Cgeru.debug.__float128_ref  Cgeru.debug.__float128    \
Chbmv.debug.__float128_ref  Chbmv.debug.__float128	\
Chemm.debug.__float128_ref  Chemm.debug.__float128	\
Chemv.debug.__float128_ref  Chemv.debug.__float128	\
Cher.debug.__float128_ref   Cher.debug.__float128 	\
Cher2.debug.__float128_ref  Cher2.debug.__float128	\
Cher2k.debug.__float128_ref Cher2k.debug.__float128	\
Cherk.debug.__float128_ref  Cherk.debug.__float128	\
Chpmv.debug.__float128_ref  Chpmv.debug.__float128	\
Chpr.debug.__float128_ref   Chpr.debug.__float128	\
Chpr2.debug.__float128_ref  Chpr2.debug.__float128	\
Crotg.debug.__float128_ref  Crotg.debug.__float128	\
Cscal.debug.__float128_ref  Cscal.debug.__float128	\
Cswap.debug.__float128_ref  Cswap.debug.__float128	\
Csymm.debug.__float128_ref  Csymm.debug.__float128    \
Csyr2k.debug.__float128_ref Csyr2k.debug.__float128	\
Csyrk.debug.__float128_ref  Csyrk.debug.__float128	\
Ctbmv.debug.__float128_ref  Ctbmv.debug.__float128	\
Ctbsv.debug.__float128_ref  Ctbsv.debug.__float128	\
Ctpmv.debug.__float128_ref  Ctpmv.debug.__float128	\
Ctpsv.debug.__float128_ref  Ctpsv.debug.__float128	\
Ctrmm.debug.__float128_ref  Ctrmm.debug.__float128	\
Ctrmv.debug.__float128_ref  Ctrmv.debug.__float128	\
Ctrsm.debug.__float128_ref  Ctrsm.debug.__float128	\
Ctrsv.debug.__float128_ref  Ctrsv.debug.__float128	\
RCabs1.debug.__float128_ref RCabs1.debug.__float128	\
RCasum.debug.__float128_ref RCasum.debug.__float128	\
RCnrm2.debug.__float128_ref RCnrm2.debug.__float128	\
Rasum.debug.__float128_ref  Rasum.debug.__float128    \
Raxpy.debug.__float128_ref  Raxpy.debug.__float128	\
Rcopy.debug.__float128_ref  Rcopy.debug.__float128	\
Rdot.debug.__float128_ref   Rdot.debug.__float128	\
Rgbmv.debug.__float128_ref  Rgbmv.debug.__float128	\
Rgemm.debug.__float128_ref  Rgemm.debug.__float128	\
Rgemv.debug.__float128_ref  Rgemv.debug.__float128	\
Rger.debug.__float128_ref   Rger.debug.__float128	\
Rnrm2.debug.__float128_ref  Rnrm2.debug.__float128	\
Rrot.debug.__float128_ref   Rrot.debug.__float128	\
Rrotg.debug.__float128_ref  Rrotg.debug.__float128	\
Rrotm.debug.__float128_ref  Rrotm.debug.__float128	\
Rsbmv.debug.__float128_ref  Rsbmv.debug.__float128	\
Rscal.debug.__float128_ref  Rscal.debug.__float128	\
Rspmv.debug.__float128_ref  Rspmv.debug.__float128    \
Rspr.debug.__float128_ref   Rspr.debug.__float128	\
Rspr2.debug.__float128_ref  Rspr2.debug.__float128	\
Rswap.debug.__float128_ref  Rswap.debug.__float128	\
Rsymm.debug.__float128_ref  Rsymm.debug.__float128	\
Rsymv.debug.__float128_ref  Rsymv.debug.__float128	\
Rsyr.debug.__float128_ref   Rsyr.debug.__float128	\
Rsyr2.debug.__float128_ref  Rsyr2.debug.__float128	\
Rsyr2k.debug.__float128_ref Rsyr2k.debug.__float128	\
Rsyrk.debug.__float128_ref  Rsyrk.debug.__float128	\
Rtbmv.debug.__float128_ref  Rtbmv.debug.__float128	\
Rtbsv.debug.__float128_ref  Rtbsv.debug.__float128	\
Rtpmv.debug.__float128_ref  Rtpmv.debug.__float128	\
Rtpsv.debug.__float128_ref  Rtpsv.debug.__float128	\
Rtrmm.debug.__float128_ref  Rtrmm.debug.__float128    \
Rtrmv.debug.__float128_ref  Rtrmv.debug.__float128	\
Rtrsm.debug.__float128_ref  Rtrsm.debug.__float128	\
Rtrsv.debug.__float128_ref  Rtrsv.debug.__float128

mblas___float128_testdir = $(prefix)/share/mpack/testing/__float128
ref_libdepends = -L$(top_builddir)/mblas/reference -lmblas___float128_ref -lmblas_mpfr_ref
opt_libdepends = -L$(top_builddir)/mblas/optimized/__float128 -lmblas___float128 -L$(top_builddir)/mblas/reference -lmblas_mpfr_ref
mplibs=-L$(top_builddir)/mpfrc++ -lmpfrcxx -L$(top_srcdir)/$(MPC_LIBDIR) -L$(top_srcdir)/$(MPFR_LIBDIR) -L$(top_srcdir)/$(GMP_LIBDIR) -lmpfr -lmpc -lgmp -lquadmath
refblas=-L$(top_builddir)/external/i/LAPACK/lib -lblas $(FCLIBS)

debug_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/$(GMP_INCLUDEDIR) -I$(top_srcdir)/$(MPFR_INCLUDEDIR) -I$(top_srcdir)/$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPACK_BUILD_WITH___FLOAT128___

if !IS_MACOSX
debug_srcdepends = ../common/xerbla.cpp ../common/mpack.debug.cpp ../common/Mxerbla.override.cpp
endif

if IS_MACOSX
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPACK_BUILD_WITH___FLOAT128___
libmxerbla_override_la_SOURCES =  ../common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS= -Wl,-flat_namespace,-undefined,dynamic_lookup
debug_srcdepends = ../common/xerbla.cpp ../common/mpack.debug.cpp
mplibs+=-Wl,-flat_namespace,-undefined,dynamic_lookup
TESTS_ENVIRONMENT=DYLD_INSERT_LIBRARIES=.libs/libmxerbla_override.0.dylib
endif

##
arithmetic_debug___float128_SOURCES=../common/arithmetic.debug.cpp $(debug_srcdepends)
arithmetic_debug___float128_CXXFLAGS=$(debug_cxxflags)
arithmetic_debug___float128_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

complex_debug___float128_SOURCES=../common/complex.debug.cpp $(debug_srcdepends)
complex_debug___float128_CXXFLAGS=$(debug_cxxflags)
complex_debug___float128_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_debug___float128_ref_SOURCES=../common/iCamax.debug.cpp $(debug_srcdepends)
iCamax_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
iCamax_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_debug___float128_SOURCES=../common/iCamax.debug.cpp $(debug_srcdepends)
iCamax_debug___float128_CXXFLAGS=$(debug_cxxflags)
iCamax_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iRamax_debug___float128_ref_SOURCES=../common/iRamax.debug.cpp $(debug_srcdepends)
iRamax_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
iRamax_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iRamax_debug___float128_SOURCES=../common/iRamax.debug.cpp $(debug_srcdepends)
iRamax_debug___float128_CXXFLAGS=$(debug_cxxflags)
iRamax_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mlsame_debug___float128_ref_SOURCES=../common/Mlsame.debug.cpp $(debug_srcdepends)
Mlsame_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Mlsame_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mlsame_debug___float128_SOURCES=../common/Mlsame.debug.cpp $(debug_srcdepends)
Mlsame_debug___float128_CXXFLAGS=$(debug_cxxflags)
Mlsame_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mxerbla_debug___float128_ref_SOURCES=../common/Mxerbla.debug.cpp $(debug_srcdepends)
Mxerbla_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Mxerbla_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mxerbla_debug___float128_SOURCES=../common/Mxerbla.debug.cpp $(debug_srcdepends)
Mxerbla_debug___float128_CXXFLAGS=$(debug_cxxflags)
Mxerbla_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRrot_debug___float128_ref_SOURCES=../common/CRrot.debug.cpp $(debug_srcdepends)
CRrot_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
CRrot_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRrot_debug___float128_SOURCES=../common/CRrot.debug.cpp $(debug_srcdepends)
CRrot_debug___float128_CXXFLAGS=$(debug_cxxflags)
CRrot_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRscal_debug___float128_ref_SOURCES=../common/CRscal.debug.cpp $(debug_srcdepends)
CRscal_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
CRscal_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRscal_debug___float128_SOURCES=../common/CRscal.debug.cpp $(debug_srcdepends)
CRscal_debug___float128_CXXFLAGS=$(debug_cxxflags)
CRscal_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_debug___float128_SOURCES=../common/Caxpy.debug.cpp $(debug_srcdepends)
Caxpy_debug___float128_CXXFLAGS=$(debug_cxxflags)
Caxpy_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_debug___float128_ref_SOURCES=../common/Caxpy.debug.cpp $(debug_srcdepends)
Caxpy_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Caxpy_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ccopy_debug___float128_SOURCES=../common/Ccopy.debug.cpp $(debug_srcdepends)
Ccopy_debug___float128_CXXFLAGS=$(debug_cxxflags)
Ccopy_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ccopy_debug___float128_ref_SOURCES=../common/Ccopy.debug.cpp $(debug_srcdepends)
Ccopy_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Ccopy_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotc_debug___float128_SOURCES=../common/Cdotc.debug.cpp $(debug_srcdepends)
Cdotc_debug___float128_CXXFLAGS=$(debug_cxxflags)
Cdotc_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotc_debug___float128_ref_SOURCES=../common/Cdotc.debug.cpp $(debug_srcdepends)
Cdotc_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Cdotc_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotu_debug___float128_SOURCES=../common/Cdotu.debug.cpp $(debug_srcdepends)
Cdotu_debug___float128_CXXFLAGS=$(debug_cxxflags)
Cdotu_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotu_debug___float128_ref_SOURCES=../common/Cdotu.debug.cpp $(debug_srcdepends)
Cdotu_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Cdotu_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgbmv_debug___float128_SOURCES=../common/Cgbmv.debug.cpp $(debug_srcdepends)
Cgbmv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Cgbmv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgbmv_debug___float128_ref_SOURCES=../common/Cgbmv.debug.cpp $(debug_srcdepends)
Cgbmv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Cgbmv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemm_debug___float128_SOURCES=../common/Cgemm.debug.cpp $(debug_srcdepends)
Cgemm_debug___float128_CXXFLAGS=$(debug_cxxflags)
Cgemm_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemm_debug___float128_ref_SOURCES=../common/Cgemm.debug.cpp $(debug_srcdepends)
Cgemm_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Cgemm_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemv_debug___float128_SOURCES=../common/Cgemv.debug.cpp $(debug_srcdepends)
Cgemv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Cgemv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemv_debug___float128_ref_SOURCES=../common/Cgemv.debug.cpp $(debug_srcdepends)
Cgemv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Cgemv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgerc_debug___float128_SOURCES=../common/Cgerc.debug.cpp $(debug_srcdepends)
Cgerc_debug___float128_CXXFLAGS=$(debug_cxxflags)
Cgerc_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgerc_debug___float128_ref_SOURCES=../common/Cgerc.debug.cpp $(debug_srcdepends)
Cgerc_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Cgerc_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgeru_debug___float128_SOURCES=../common/Cgeru.debug.cpp $(debug_srcdepends)
Cgeru_debug___float128_CXXFLAGS=$(debug_cxxflags)
Cgeru_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgeru_debug___float128_ref_SOURCES=../common/Cgeru.debug.cpp $(debug_srcdepends)
Cgeru_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Cgeru_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chbmv_debug___float128_SOURCES=../common/Chbmv.debug.cpp $(debug_srcdepends)
Chbmv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Chbmv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chbmv_debug___float128_ref_SOURCES=../common/Chbmv.debug.cpp $(debug_srcdepends)
Chbmv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Chbmv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemm_debug___float128_SOURCES=../common/Chemm.debug.cpp $(debug_srcdepends)
Chemm_debug___float128_CXXFLAGS=$(debug_cxxflags)
Chemm_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemm_debug___float128_ref_SOURCES=../common/Chemm.debug.cpp $(debug_srcdepends)
Chemm_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Chemm_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemv_debug___float128_SOURCES=../common/Chemv.debug.cpp $(debug_srcdepends)
Chemv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Chemv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemv_debug___float128_ref_SOURCES=../common/Chemv.debug.cpp $(debug_srcdepends)
Chemv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Chemv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher_debug___float128_SOURCES=../common/Cher.debug.cpp $(debug_srcdepends)
Cher_debug___float128_CXXFLAGS=$(debug_cxxflags)
Cher_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher_debug___float128_ref_SOURCES=../common/Cher.debug.cpp $(debug_srcdepends)
Cher_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Cher_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2_debug___float128_SOURCES=../common/Cher2.debug.cpp $(debug_srcdepends)
Cher2_debug___float128_CXXFLAGS=$(debug_cxxflags)
Cher2_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2_debug___float128_ref_SOURCES=../common/Cher2.debug.cpp $(debug_srcdepends)
Cher2_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Cher2_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2k_debug___float128_SOURCES=../common/Cher2k.debug.cpp $(debug_srcdepends)
Cher2k_debug___float128_CXXFLAGS=$(debug_cxxflags)
Cher2k_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2k_debug___float128_ref_SOURCES=../common/Cher2k.debug.cpp $(debug_srcdepends)
Cher2k_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Cher2k_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cherk_debug___float128_SOURCES=../common/Cherk.debug.cpp $(debug_srcdepends)
Cherk_debug___float128_CXXFLAGS=$(debug_cxxflags)
Cherk_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cherk_debug___float128_ref_SOURCES=../common/Cherk.debug.cpp $(debug_srcdepends)
Cherk_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Cherk_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpmv_debug___float128_SOURCES=../common/Chpmv.debug.cpp $(debug_srcdepends)
Chpmv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Chpmv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpmv_debug___float128_ref_SOURCES=../common/Chpmv.debug.cpp $(debug_srcdepends)
Chpmv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Chpmv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr_debug___float128_SOURCES=../common/Chpr.debug.cpp $(debug_srcdepends)
Chpr_debug___float128_CXXFLAGS=$(debug_cxxflags)
Chpr_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr_debug___float128_ref_SOURCES=../common/Chpr.debug.cpp $(debug_srcdepends)
Chpr_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Chpr_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr2_debug___float128_SOURCES=../common/Chpr2.debug.cpp $(debug_srcdepends)
Chpr2_debug___float128_CXXFLAGS=$(debug_cxxflags)
Chpr2_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr2_debug___float128_ref_SOURCES=../common/Chpr2.debug.cpp $(debug_srcdepends)
Chpr2_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Chpr2_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Crotg_debug___float128_SOURCES=../common/Crotg.debug.cpp $(debug_srcdepends)
Crotg_debug___float128_CXXFLAGS=$(debug_cxxflags)
Crotg_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Crotg_debug___float128_ref_SOURCES=../common/Crotg.debug.cpp $(debug_srcdepends)
Crotg_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Crotg_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cscal_debug___float128_SOURCES=../common/Cscal.debug.cpp $(debug_srcdepends)
Cscal_debug___float128_CXXFLAGS=$(debug_cxxflags)
Cscal_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cscal_debug___float128_ref_SOURCES=../common/Cscal.debug.cpp $(debug_srcdepends)
Cscal_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Cscal_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cswap_debug___float128_SOURCES=../common/Cswap.debug.cpp $(debug_srcdepends)
Cswap_debug___float128_CXXFLAGS=$(debug_cxxflags)
Cswap_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cswap_debug___float128_ref_SOURCES=../common/Cswap.debug.cpp $(debug_srcdepends)
Cswap_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Cswap_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csymm_debug___float128_SOURCES=../common/Csymm.debug.cpp $(debug_srcdepends)
Csymm_debug___float128_CXXFLAGS=$(debug_cxxflags)
Csymm_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csymm_debug___float128_ref_SOURCES=../common/Csymm.debug.cpp $(debug_srcdepends)
Csymm_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Csymm_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyr2k_debug___float128_SOURCES=../common/Csyr2k.debug.cpp $(debug_srcdepends)
Csyr2k_debug___float128_CXXFLAGS=$(debug_cxxflags)
Csyr2k_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyr2k_debug___float128_ref_SOURCES=../common/Csyr2k.debug.cpp $(debug_srcdepends)
Csyr2k_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Csyr2k_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyrk_debug___float128_SOURCES=../common/Csyrk.debug.cpp $(debug_srcdepends)
Csyrk_debug___float128_CXXFLAGS=$(debug_cxxflags)
Csyrk_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyrk_debug___float128_ref_SOURCES=../common/Csyrk.debug.cpp $(debug_srcdepends)
Csyrk_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Csyrk_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbmv_debug___float128_SOURCES=../common/Ctbmv.debug.cpp $(debug_srcdepends)
Ctbmv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Ctbmv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbmv_debug___float128_ref_SOURCES=../common/Ctbmv.debug.cpp $(debug_srcdepends)
Ctbmv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctbmv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbsv_debug___float128_SOURCES=../common/Ctbsv.debug.cpp $(debug_srcdepends)
Ctbsv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Ctbsv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbsv_debug___float128_ref_SOURCES=../common/Ctbsv.debug.cpp $(debug_srcdepends)
Ctbsv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctbsv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpmv_debug___float128_SOURCES=../common/Ctpmv.debug.cpp $(debug_srcdepends)
Ctpmv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Ctpmv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpmv_debug___float128_ref_SOURCES=../common/Ctpmv.debug.cpp $(debug_srcdepends)
Ctpmv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctpmv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpsv_debug___float128_SOURCES=../common/Ctpsv.debug.cpp $(debug_srcdepends)
Ctpsv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Ctpsv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpsv_debug___float128_ref_SOURCES=../common/Ctpsv.debug.cpp $(debug_srcdepends)
Ctpsv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctpsv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmm_debug___float128_SOURCES=../common/Ctrmm.debug.cpp $(debug_srcdepends)
Ctrmm_debug___float128_CXXFLAGS=$(debug_cxxflags)
Ctrmm_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmm_debug___float128_ref_SOURCES=../common/Ctrmm.debug.cpp $(debug_srcdepends)
Ctrmm_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctrmm_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmv_debug___float128_SOURCES=../common/Ctrmv.debug.cpp $(debug_srcdepends)
Ctrmv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Ctrmv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmv_debug___float128_ref_SOURCES=../common/Ctrmv.debug.cpp $(debug_srcdepends)
Ctrmv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctrmv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsm_debug___float128_SOURCES=../common/Ctrsm.debug.cpp $(debug_srcdepends)
Ctrsm_debug___float128_CXXFLAGS=$(debug_cxxflags)
Ctrsm_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsm_debug___float128_ref_SOURCES=../common/Ctrsm.debug.cpp $(debug_srcdepends)
Ctrsm_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctrsm_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsv_debug___float128_SOURCES=../common/Ctrsv.debug.cpp $(debug_srcdepends)
Ctrsv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Ctrsv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsv_debug___float128_ref_SOURCES=../common/Ctrsv.debug.cpp $(debug_srcdepends)
Ctrsv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctrsv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_debug___float128_ref_SOURCES=../common/RCabs1.debug.cpp $(debug_srcdepends)
RCabs1_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
RCabs1_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_debug___float128_SOURCES=../common/RCabs1.debug.cpp $(debug_srcdepends)
RCabs1_debug___float128_CXXFLAGS=$(debug_cxxflags)
RCabs1_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCasum_debug___float128_ref_SOURCES=../common/RCasum.debug.cpp $(debug_srcdepends)
RCasum_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
RCasum_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCasum_debug___float128_SOURCES=../common/RCasum.debug.cpp $(debug_srcdepends)
RCasum_debug___float128_CXXFLAGS=$(debug_cxxflags)
RCasum_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCnrm2_debug___float128_ref_SOURCES=../common/RCnrm2.debug.cpp $(debug_srcdepends)
RCnrm2_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
RCnrm2_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCnrm2_debug___float128_SOURCES=../common/RCnrm2.debug.cpp $(debug_srcdepends)
RCnrm2_debug___float128_CXXFLAGS=$(debug_cxxflags)
RCnrm2_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rasum_debug___float128_ref_SOURCES=../common/Rasum.debug.cpp $(debug_srcdepends)
Rasum_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rasum_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rasum_debug___float128_SOURCES=../common/Rasum.debug.cpp $(debug_srcdepends)
Rasum_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rasum_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Raxpy_debug___float128_ref_SOURCES=../common/Raxpy.debug.cpp $(debug_srcdepends)
Raxpy_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Raxpy_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Raxpy_debug___float128_SOURCES=../common/Raxpy.debug.cpp $(debug_srcdepends)
Raxpy_debug___float128_CXXFLAGS=$(debug_cxxflags)
Raxpy_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rcopy_debug___float128_ref_SOURCES=../common/Rcopy.debug.cpp $(debug_srcdepends)
Rcopy_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rcopy_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rcopy_debug___float128_SOURCES=../common/Rcopy.debug.cpp $(debug_srcdepends)
Rcopy_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rcopy_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rdot_debug___float128_ref_SOURCES=../common/Rdot.debug.cpp $(debug_srcdepends)
Rdot_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rdot_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rdot_debug___float128_SOURCES=../common/Rdot.debug.cpp $(debug_srcdepends)
Rdot_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rdot_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgbmv_debug___float128_ref_SOURCES=../common/Rgbmv.debug.cpp $(debug_srcdepends)
Rgbmv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rgbmv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgbmv_debug___float128_SOURCES=../common/Rgbmv.debug.cpp $(debug_srcdepends)
Rgbmv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rgbmv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemm_debug___float128_ref_SOURCES=../common/Rgemm.debug.cpp $(debug_srcdepends)
Rgemm_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rgemm_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemm_debug___float128_SOURCES=../common/Rgemm.debug.cpp $(debug_srcdepends)
Rgemm_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rgemm_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemv_debug___float128_ref_SOURCES=../common/Rgemv.debug.cpp $(debug_srcdepends)
Rgemv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rgemv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemv_debug___float128_SOURCES=../common/Rgemv.debug.cpp $(debug_srcdepends)
Rgemv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rgemv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rger_debug___float128_ref_SOURCES=../common/Rger.debug.cpp $(debug_srcdepends)
Rger_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rger_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rger_debug___float128_SOURCES=../common/Rger.debug.cpp $(debug_srcdepends)
Rger_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rger_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rnrm2_debug___float128_ref_SOURCES=../common/Rnrm2.debug.cpp $(debug_srcdepends)
Rnrm2_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rnrm2_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rnrm2_debug___float128_SOURCES=../common/Rnrm2.debug.cpp $(debug_srcdepends)
Rnrm2_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rnrm2_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrot_debug___float128_ref_SOURCES=../common/Rrot.debug.cpp $(debug_srcdepends)
Rrot_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rrot_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrot_debug___float128_SOURCES=../common/Rrot.debug.cpp $(debug_srcdepends)
Rrot_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rrot_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotg_debug___float128_ref_SOURCES=../common/Rrotg.debug.cpp $(debug_srcdepends)
Rrotg_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rrotg_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotg_debug___float128_SOURCES=../common/Rrotg.debug.cpp $(debug_srcdepends)
Rrotg_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rrotg_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotm_debug___float128_ref_SOURCES=../common/Rrotm.debug.cpp $(debug_srcdepends)
Rrotm_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rrotm_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotm_debug___float128_SOURCES=../common/Rrotm.debug.cpp $(debug_srcdepends)
Rrotm_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rrotm_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsbmv_debug___float128_ref_SOURCES=../common/Rsbmv.debug.cpp $(debug_srcdepends)
Rsbmv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rsbmv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsbmv_debug___float128_SOURCES=../common/Rsbmv.debug.cpp $(debug_srcdepends)
Rsbmv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rsbmv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rscal_debug___float128_ref_SOURCES=../common/Rscal.debug.cpp $(debug_srcdepends)
Rscal_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rscal_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rscal_debug___float128_SOURCES=../common/Rscal.debug.cpp $(debug_srcdepends)
Rscal_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rscal_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspmv_debug___float128_ref_SOURCES=../common/Rspmv.debug.cpp $(debug_srcdepends)
Rspmv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rspmv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspmv_debug___float128_SOURCES=../common/Rspmv.debug.cpp $(debug_srcdepends)
Rspmv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rspmv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr_debug___float128_ref_SOURCES=../common/Rspr.debug.cpp $(debug_srcdepends)
Rspr_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rspr_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr_debug___float128_SOURCES=../common/Rspr.debug.cpp $(debug_srcdepends)
Rspr_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rspr_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr2_debug___float128_ref_SOURCES=../common/Rspr2.debug.cpp $(debug_srcdepends)
Rspr2_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rspr2_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr2_debug___float128_SOURCES=../common/Rspr2.debug.cpp $(debug_srcdepends)
Rspr2_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rspr2_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rswap_debug___float128_ref_SOURCES=../common/Rswap.debug.cpp $(debug_srcdepends)
Rswap_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rswap_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rswap_debug___float128_SOURCES=../common/Rswap.debug.cpp $(debug_srcdepends)
Rswap_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rswap_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymm_debug___float128_ref_SOURCES=../common/Rsymm.debug.cpp $(debug_srcdepends)
Rsymm_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rsymm_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymm_debug___float128_SOURCES=../common/Rsymm.debug.cpp $(debug_srcdepends)
Rsymm_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rsymm_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymv_debug___float128_ref_SOURCES=../common/Rsymv.debug.cpp $(debug_srcdepends)
Rsymv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rsymv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymv_debug___float128_SOURCES=../common/Rsymv.debug.cpp $(debug_srcdepends)
Rsymv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rsymv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr_debug___float128_ref_SOURCES=../common/Rsyr.debug.cpp $(debug_srcdepends)
Rsyr_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr_debug___float128_SOURCES=../common/Rsyr.debug.cpp $(debug_srcdepends)
Rsyr_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rsyr_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2_debug___float128_ref_SOURCES=../common/Rsyr2.debug.cpp $(debug_srcdepends)
Rsyr2_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr2_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2_debug___float128_SOURCES=../common/Rsyr2.debug.cpp $(debug_srcdepends)
Rsyr2_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rsyr2_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2k_debug___float128_ref_SOURCES=../common/Rsyr2k.debug.cpp $(debug_srcdepends)
Rsyr2k_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr2k_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2k_debug___float128_SOURCES=../common/Rsyr2k.debug.cpp $(debug_srcdepends)
Rsyr2k_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rsyr2k_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyrk_debug___float128_ref_SOURCES=../common/Rsyrk.debug.cpp $(debug_srcdepends)
Rsyrk_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rsyrk_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyrk_debug___float128_SOURCES=../common/Rsyrk.debug.cpp $(debug_srcdepends)
Rsyrk_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rsyrk_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbmv_debug___float128_ref_SOURCES=../common/Rtbmv.debug.cpp $(debug_srcdepends)
Rtbmv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtbmv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbmv_debug___float128_SOURCES=../common/Rtbmv.debug.cpp $(debug_srcdepends)
Rtbmv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rtbmv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbsv_debug___float128_ref_SOURCES=../common/Rtbsv.debug.cpp $(debug_srcdepends)
Rtbsv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtbsv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbsv_debug___float128_SOURCES=../common/Rtbsv.debug.cpp $(debug_srcdepends)
Rtbsv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rtbsv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpmv_debug___float128_ref_SOURCES=../common/Rtpmv.debug.cpp $(debug_srcdepends)
Rtpmv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtpmv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpmv_debug___float128_SOURCES=../common/Rtpmv.debug.cpp $(debug_srcdepends)
Rtpmv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rtpmv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpsv_debug___float128_ref_SOURCES=../common/Rtpsv.debug.cpp $(debug_srcdepends)
Rtpsv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtpsv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpsv_debug___float128_SOURCES=../common/Rtpsv.debug.cpp $(debug_srcdepends)
Rtpsv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rtpsv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmm_debug___float128_ref_SOURCES=../common/Rtrmm.debug.cpp $(debug_srcdepends)
Rtrmm_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtrmm_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmm_debug___float128_SOURCES=../common/Rtrmm.debug.cpp $(debug_srcdepends)
Rtrmm_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rtrmm_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmv_debug___float128_ref_SOURCES=../common/Rtrmv.debug.cpp $(debug_srcdepends)
Rtrmv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtrmv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmv_debug___float128_SOURCES=../common/Rtrmv.debug.cpp $(debug_srcdepends)
Rtrmv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rtrmv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsm_debug___float128_ref_SOURCES=../common/Rtrsm.debug.cpp $(debug_srcdepends)
Rtrsm_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtrsm_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsm_debug___float128_SOURCES=../common/Rtrsm.debug.cpp $(debug_srcdepends)
Rtrsm_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rtrsm_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsv_debug___float128_ref_SOURCES=../common/Rtrsv.debug.cpp $(debug_srcdepends)
Rtrsv_debug___float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtrsv_debug___float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsv_debug___float128_SOURCES=../common/Rtrsv.debug.cpp $(debug_srcdepends)
Rtrsv_debug___float128_CXXFLAGS=$(debug_cxxflags)
Rtrsv_debug___float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)
