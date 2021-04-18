#CTBSV, CTPSV might have bugs or intrinsically inaccurate
#Rnrm2 should be checked more
check_PROGRAMS = $(mpblas__Float128_test_PROGRAMS) 
TESTS = $(check_PROGRAMS)
mpblas__Float128_test_PROGRAMS = \
arithmetic.debug__Float128 \
complex.debug__Float128 \
iCamax.debug._Float128_ref iCamax.debug._Float128   \
iRamax.debug._Float128_ref iRamax.debug._Float128   \
Mlsame.debug._Float128_ref Mlsame.debug._Float128   \
Mxerbla.debug._Float128_ref Mxerbla.debug._Float128 \
CRrot.debug._Float128_ref  CRrot.debug._Float128    \
CRscal.debug._Float128_ref CRscal.debug._Float128   \
Caxpy.debug._Float128_ref  Caxpy.debug._Float128    \
Ccopy.debug._Float128_ref  Ccopy.debug._Float128    \
Cdotc.debug._Float128_ref  Cdotc.debug._Float128	\
Cdotu.debug._Float128_ref  Cdotu.debug._Float128	\
Cgbmv.debug._Float128_ref  Cgbmv.debug._Float128	\
Cgemm.debug._Float128_ref  Cgemm.debug._Float128	\
Cgemv.debug._Float128_ref  Cgemv.debug._Float128	\
Cgerc.debug._Float128_ref  Cgerc.debug._Float128	\
Cgeru.debug._Float128_ref  Cgeru.debug._Float128    \
Chbmv.debug._Float128_ref  Chbmv.debug._Float128	\
Chemm.debug._Float128_ref  Chemm.debug._Float128	\
Chemv.debug._Float128_ref  Chemv.debug._Float128	\
Cher.debug._Float128_ref   Cher.debug._Float128 	\
Cher2.debug._Float128_ref  Cher2.debug._Float128	\
Cher2k.debug._Float128_ref Cher2k.debug._Float128	\
Cherk.debug._Float128_ref  Cherk.debug._Float128	\
Chpmv.debug._Float128_ref  Chpmv.debug._Float128	\
Chpr.debug._Float128_ref   Chpr.debug._Float128	\
Chpr2.debug._Float128_ref  Chpr2.debug._Float128	\
Crotg.debug._Float128_ref  Crotg.debug._Float128	\
Cscal.debug._Float128_ref  Cscal.debug._Float128	\
Cswap.debug._Float128_ref  Cswap.debug._Float128	\
Csymm.debug._Float128_ref  Csymm.debug._Float128    \
Csyr2k.debug._Float128_ref Csyr2k.debug._Float128	\
Csyrk.debug._Float128_ref  Csyrk.debug._Float128	\
Ctbmv.debug._Float128_ref  Ctbmv.debug._Float128	\
Ctbsv.debug._Float128_ref  Ctbsv.debug._Float128	\
Ctpmv.debug._Float128_ref  Ctpmv.debug._Float128	\
Ctpsv.debug._Float128_ref  Ctpsv.debug._Float128	\
Ctrmm.debug._Float128_ref  Ctrmm.debug._Float128	\
Ctrmv.debug._Float128_ref  Ctrmv.debug._Float128	\
Ctrsm.debug._Float128_ref  Ctrsm.debug._Float128	\
Ctrsv.debug._Float128_ref  Ctrsv.debug._Float128	\
RCabs1.debug._Float128_ref RCabs1.debug._Float128	\
RCasum.debug._Float128_ref RCasum.debug._Float128	\
RCnrm2.debug._Float128_ref RCnrm2.debug._Float128	\
Rasum.debug._Float128_ref  Rasum.debug._Float128    \
Raxpy.debug._Float128_ref  Raxpy.debug._Float128	\
Rcopy.debug._Float128_ref  Rcopy.debug._Float128	\
Rdot.debug._Float128_ref   Rdot.debug._Float128	\
Rgbmv.debug._Float128_ref  Rgbmv.debug._Float128	\
Rgemm.debug._Float128_ref  Rgemm.debug._Float128	\
Rgemv.debug._Float128_ref  Rgemv.debug._Float128	\
Rger.debug._Float128_ref   Rger.debug._Float128	\
Rnrm2.debug._Float128_ref  Rnrm2.debug._Float128	\
Rrot.debug._Float128_ref   Rrot.debug._Float128	\
Rrotg.debug._Float128_ref  Rrotg.debug._Float128	\
Rrotm.debug._Float128_ref  Rrotm.debug._Float128	\
Rsbmv.debug._Float128_ref  Rsbmv.debug._Float128	\
Rscal.debug._Float128_ref  Rscal.debug._Float128	\
Rspmv.debug._Float128_ref  Rspmv.debug._Float128    \
Rspr.debug._Float128_ref   Rspr.debug._Float128	\
Rspr2.debug._Float128_ref  Rspr2.debug._Float128	\
Rswap.debug._Float128_ref  Rswap.debug._Float128	\
Rsymm.debug._Float128_ref  Rsymm.debug._Float128	\
Rsymv.debug._Float128_ref  Rsymv.debug._Float128	\
Rsyr.debug._Float128_ref   Rsyr.debug._Float128	\
Rsyr2.debug._Float128_ref  Rsyr2.debug._Float128	\
Rsyr2k.debug._Float128_ref Rsyr2k.debug._Float128	\
Rsyrk.debug._Float128_ref  Rsyrk.debug._Float128	\
Rtbmv.debug._Float128_ref  Rtbmv.debug._Float128	\
Rtbsv.debug._Float128_ref  Rtbsv.debug._Float128	\
Rtpmv.debug._Float128_ref  Rtpmv.debug._Float128	\
Rtpsv.debug._Float128_ref  Rtpsv.debug._Float128	\
Rtrmm.debug._Float128_ref  Rtrmm.debug._Float128    \
Rtrmv.debug._Float128_ref  Rtrmv.debug._Float128	\
Rtrsm.debug._Float128_ref  Rtrsm.debug._Float128	\
Rtrsv.debug._Float128_ref  Rtrsv.debug._Float128

mpblas__Float128_testdir = $(prefix)/share/mplapack/debug/mpblas/_Float128
ref_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas__Float128_ref -lmpblas_mpfr_ref
opt_libdepends = -L$(top_builddir)/mpblas/optimized/_Float128 -lmpblas__Float128 -L$(top_builddir)/mpblas/reference -lmpblas_mpfr_ref
mplibs=-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp
refblas=-L$(top_builddir)/external/i/LAPACK/lib -lblas $(FCLIBS)

debug_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH__FLOAT128___

if !IS_MACOSX
debug_srcdepends = ../common/xerbla.cpp ../common/mplapack.debug.cpp ../common/Mxerbla.override.cpp
endif

if IS_MACOSX
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT128___
libmxerbla_override_la_SOURCES =  ../common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS= -Wl,-flat_namespace,-undefined,dynamic_lookup
debug_srcdepends = ../common/xerbla.cpp ../common/mplapack.debug.cpp
mplibs+=-Wl,-flat_namespace,-undefined,dynamic_lookup
TESTS_ENVIRONMENT=DYLD_INSERT_LIBRARIES=.libs/libmxerbla_override.0.dylib
endif

##
arithmetic_debug__Float128_SOURCES=../common/arithmetic.debug.cpp $(debug_srcdepends)
arithmetic_debug__Float128_CXXFLAGS=$(debug_cxxflags)
arithmetic_debug__Float128_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

complex_debug__Float128_SOURCES=../common/complex.debug.cpp $(debug_srcdepends)
complex_debug__Float128_CXXFLAGS=$(debug_cxxflags)
complex_debug__Float128_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_debug__Float128_ref_SOURCES=../common/iCamax.debug.cpp $(debug_srcdepends)
iCamax_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
iCamax_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_debug__Float128_SOURCES=../common/iCamax.debug.cpp $(debug_srcdepends)
iCamax_debug__Float128_CXXFLAGS=$(debug_cxxflags)
iCamax_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iRamax_debug__Float128_ref_SOURCES=../common/iRamax.debug.cpp $(debug_srcdepends)
iRamax_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
iRamax_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iRamax_debug__Float128_SOURCES=../common/iRamax.debug.cpp $(debug_srcdepends)
iRamax_debug__Float128_CXXFLAGS=$(debug_cxxflags)
iRamax_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mlsame_debug__Float128_ref_SOURCES=../common/Mlsame.debug.cpp $(debug_srcdepends)
Mlsame_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Mlsame_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mlsame_debug__Float128_SOURCES=../common/Mlsame.debug.cpp $(debug_srcdepends)
Mlsame_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Mlsame_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mxerbla_debug__Float128_ref_SOURCES=../common/Mxerbla.debug.cpp $(debug_srcdepends)
Mxerbla_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Mxerbla_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mxerbla_debug__Float128_SOURCES=../common/Mxerbla.debug.cpp $(debug_srcdepends)
Mxerbla_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Mxerbla_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRrot_debug__Float128_ref_SOURCES=../common/CRrot.debug.cpp $(debug_srcdepends)
CRrot_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
CRrot_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRrot_debug__Float128_SOURCES=../common/CRrot.debug.cpp $(debug_srcdepends)
CRrot_debug__Float128_CXXFLAGS=$(debug_cxxflags)
CRrot_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRscal_debug__Float128_ref_SOURCES=../common/CRscal.debug.cpp $(debug_srcdepends)
CRscal_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
CRscal_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRscal_debug__Float128_SOURCES=../common/CRscal.debug.cpp $(debug_srcdepends)
CRscal_debug__Float128_CXXFLAGS=$(debug_cxxflags)
CRscal_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_debug__Float128_SOURCES=../common/Caxpy.debug.cpp $(debug_srcdepends)
Caxpy_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Caxpy_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_debug__Float128_ref_SOURCES=../common/Caxpy.debug.cpp $(debug_srcdepends)
Caxpy_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Caxpy_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ccopy_debug__Float128_SOURCES=../common/Ccopy.debug.cpp $(debug_srcdepends)
Ccopy_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Ccopy_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ccopy_debug__Float128_ref_SOURCES=../common/Ccopy.debug.cpp $(debug_srcdepends)
Ccopy_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Ccopy_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotc_debug__Float128_SOURCES=../common/Cdotc.debug.cpp $(debug_srcdepends)
Cdotc_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Cdotc_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotc_debug__Float128_ref_SOURCES=../common/Cdotc.debug.cpp $(debug_srcdepends)
Cdotc_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Cdotc_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotu_debug__Float128_SOURCES=../common/Cdotu.debug.cpp $(debug_srcdepends)
Cdotu_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Cdotu_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotu_debug__Float128_ref_SOURCES=../common/Cdotu.debug.cpp $(debug_srcdepends)
Cdotu_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Cdotu_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgbmv_debug__Float128_SOURCES=../common/Cgbmv.debug.cpp $(debug_srcdepends)
Cgbmv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Cgbmv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgbmv_debug__Float128_ref_SOURCES=../common/Cgbmv.debug.cpp $(debug_srcdepends)
Cgbmv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Cgbmv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemm_debug__Float128_SOURCES=../common/Cgemm.debug.cpp $(debug_srcdepends)
Cgemm_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Cgemm_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemm_debug__Float128_ref_SOURCES=../common/Cgemm.debug.cpp $(debug_srcdepends)
Cgemm_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Cgemm_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemv_debug__Float128_SOURCES=../common/Cgemv.debug.cpp $(debug_srcdepends)
Cgemv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Cgemv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemv_debug__Float128_ref_SOURCES=../common/Cgemv.debug.cpp $(debug_srcdepends)
Cgemv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Cgemv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgerc_debug__Float128_SOURCES=../common/Cgerc.debug.cpp $(debug_srcdepends)
Cgerc_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Cgerc_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgerc_debug__Float128_ref_SOURCES=../common/Cgerc.debug.cpp $(debug_srcdepends)
Cgerc_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Cgerc_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgeru_debug__Float128_SOURCES=../common/Cgeru.debug.cpp $(debug_srcdepends)
Cgeru_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Cgeru_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgeru_debug__Float128_ref_SOURCES=../common/Cgeru.debug.cpp $(debug_srcdepends)
Cgeru_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Cgeru_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chbmv_debug__Float128_SOURCES=../common/Chbmv.debug.cpp $(debug_srcdepends)
Chbmv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Chbmv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chbmv_debug__Float128_ref_SOURCES=../common/Chbmv.debug.cpp $(debug_srcdepends)
Chbmv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Chbmv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemm_debug__Float128_SOURCES=../common/Chemm.debug.cpp $(debug_srcdepends)
Chemm_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Chemm_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemm_debug__Float128_ref_SOURCES=../common/Chemm.debug.cpp $(debug_srcdepends)
Chemm_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Chemm_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemv_debug__Float128_SOURCES=../common/Chemv.debug.cpp $(debug_srcdepends)
Chemv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Chemv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemv_debug__Float128_ref_SOURCES=../common/Chemv.debug.cpp $(debug_srcdepends)
Chemv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Chemv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher_debug__Float128_SOURCES=../common/Cher.debug.cpp $(debug_srcdepends)
Cher_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Cher_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher_debug__Float128_ref_SOURCES=../common/Cher.debug.cpp $(debug_srcdepends)
Cher_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Cher_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2_debug__Float128_SOURCES=../common/Cher2.debug.cpp $(debug_srcdepends)
Cher2_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Cher2_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2_debug__Float128_ref_SOURCES=../common/Cher2.debug.cpp $(debug_srcdepends)
Cher2_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Cher2_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2k_debug__Float128_SOURCES=../common/Cher2k.debug.cpp $(debug_srcdepends)
Cher2k_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Cher2k_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2k_debug__Float128_ref_SOURCES=../common/Cher2k.debug.cpp $(debug_srcdepends)
Cher2k_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Cher2k_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cherk_debug__Float128_SOURCES=../common/Cherk.debug.cpp $(debug_srcdepends)
Cherk_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Cherk_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cherk_debug__Float128_ref_SOURCES=../common/Cherk.debug.cpp $(debug_srcdepends)
Cherk_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Cherk_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpmv_debug__Float128_SOURCES=../common/Chpmv.debug.cpp $(debug_srcdepends)
Chpmv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Chpmv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpmv_debug__Float128_ref_SOURCES=../common/Chpmv.debug.cpp $(debug_srcdepends)
Chpmv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Chpmv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr_debug__Float128_SOURCES=../common/Chpr.debug.cpp $(debug_srcdepends)
Chpr_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Chpr_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr_debug__Float128_ref_SOURCES=../common/Chpr.debug.cpp $(debug_srcdepends)
Chpr_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Chpr_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr2_debug__Float128_SOURCES=../common/Chpr2.debug.cpp $(debug_srcdepends)
Chpr2_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Chpr2_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr2_debug__Float128_ref_SOURCES=../common/Chpr2.debug.cpp $(debug_srcdepends)
Chpr2_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Chpr2_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Crotg_debug__Float128_SOURCES=../common/Crotg.debug.cpp $(debug_srcdepends)
Crotg_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Crotg_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Crotg_debug__Float128_ref_SOURCES=../common/Crotg.debug.cpp $(debug_srcdepends)
Crotg_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Crotg_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cscal_debug__Float128_SOURCES=../common/Cscal.debug.cpp $(debug_srcdepends)
Cscal_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Cscal_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cscal_debug__Float128_ref_SOURCES=../common/Cscal.debug.cpp $(debug_srcdepends)
Cscal_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Cscal_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cswap_debug__Float128_SOURCES=../common/Cswap.debug.cpp $(debug_srcdepends)
Cswap_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Cswap_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cswap_debug__Float128_ref_SOURCES=../common/Cswap.debug.cpp $(debug_srcdepends)
Cswap_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Cswap_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csymm_debug__Float128_SOURCES=../common/Csymm.debug.cpp $(debug_srcdepends)
Csymm_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Csymm_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csymm_debug__Float128_ref_SOURCES=../common/Csymm.debug.cpp $(debug_srcdepends)
Csymm_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Csymm_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyr2k_debug__Float128_SOURCES=../common/Csyr2k.debug.cpp $(debug_srcdepends)
Csyr2k_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Csyr2k_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyr2k_debug__Float128_ref_SOURCES=../common/Csyr2k.debug.cpp $(debug_srcdepends)
Csyr2k_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Csyr2k_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyrk_debug__Float128_SOURCES=../common/Csyrk.debug.cpp $(debug_srcdepends)
Csyrk_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Csyrk_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyrk_debug__Float128_ref_SOURCES=../common/Csyrk.debug.cpp $(debug_srcdepends)
Csyrk_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Csyrk_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbmv_debug__Float128_SOURCES=../common/Ctbmv.debug.cpp $(debug_srcdepends)
Ctbmv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Ctbmv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbmv_debug__Float128_ref_SOURCES=../common/Ctbmv.debug.cpp $(debug_srcdepends)
Ctbmv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctbmv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbsv_debug__Float128_SOURCES=../common/Ctbsv.debug.cpp $(debug_srcdepends)
Ctbsv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Ctbsv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbsv_debug__Float128_ref_SOURCES=../common/Ctbsv.debug.cpp $(debug_srcdepends)
Ctbsv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctbsv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpmv_debug__Float128_SOURCES=../common/Ctpmv.debug.cpp $(debug_srcdepends)
Ctpmv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Ctpmv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpmv_debug__Float128_ref_SOURCES=../common/Ctpmv.debug.cpp $(debug_srcdepends)
Ctpmv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctpmv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpsv_debug__Float128_SOURCES=../common/Ctpsv.debug.cpp $(debug_srcdepends)
Ctpsv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Ctpsv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpsv_debug__Float128_ref_SOURCES=../common/Ctpsv.debug.cpp $(debug_srcdepends)
Ctpsv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctpsv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmm_debug__Float128_SOURCES=../common/Ctrmm.debug.cpp $(debug_srcdepends)
Ctrmm_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Ctrmm_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmm_debug__Float128_ref_SOURCES=../common/Ctrmm.debug.cpp $(debug_srcdepends)
Ctrmm_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctrmm_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmv_debug__Float128_SOURCES=../common/Ctrmv.debug.cpp $(debug_srcdepends)
Ctrmv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Ctrmv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmv_debug__Float128_ref_SOURCES=../common/Ctrmv.debug.cpp $(debug_srcdepends)
Ctrmv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctrmv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsm_debug__Float128_SOURCES=../common/Ctrsm.debug.cpp $(debug_srcdepends)
Ctrsm_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Ctrsm_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsm_debug__Float128_ref_SOURCES=../common/Ctrsm.debug.cpp $(debug_srcdepends)
Ctrsm_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctrsm_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsv_debug__Float128_SOURCES=../common/Ctrsv.debug.cpp $(debug_srcdepends)
Ctrsv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Ctrsv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsv_debug__Float128_ref_SOURCES=../common/Ctrsv.debug.cpp $(debug_srcdepends)
Ctrsv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Ctrsv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_debug__Float128_ref_SOURCES=../common/RCabs1.debug.cpp $(debug_srcdepends)
RCabs1_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
RCabs1_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_debug__Float128_SOURCES=../common/RCabs1.debug.cpp $(debug_srcdepends)
RCabs1_debug__Float128_CXXFLAGS=$(debug_cxxflags)
RCabs1_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCasum_debug__Float128_ref_SOURCES=../common/RCasum.debug.cpp $(debug_srcdepends)
RCasum_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
RCasum_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCasum_debug__Float128_SOURCES=../common/RCasum.debug.cpp $(debug_srcdepends)
RCasum_debug__Float128_CXXFLAGS=$(debug_cxxflags)
RCasum_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCnrm2_debug__Float128_ref_SOURCES=../common/RCnrm2.debug.cpp $(debug_srcdepends)
RCnrm2_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
RCnrm2_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCnrm2_debug__Float128_SOURCES=../common/RCnrm2.debug.cpp $(debug_srcdepends)
RCnrm2_debug__Float128_CXXFLAGS=$(debug_cxxflags)
RCnrm2_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rasum_debug__Float128_ref_SOURCES=../common/Rasum.debug.cpp $(debug_srcdepends)
Rasum_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rasum_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rasum_debug__Float128_SOURCES=../common/Rasum.debug.cpp $(debug_srcdepends)
Rasum_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rasum_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Raxpy_debug__Float128_ref_SOURCES=../common/Raxpy.debug.cpp $(debug_srcdepends)
Raxpy_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Raxpy_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Raxpy_debug__Float128_SOURCES=../common/Raxpy.debug.cpp $(debug_srcdepends)
Raxpy_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Raxpy_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rcopy_debug__Float128_ref_SOURCES=../common/Rcopy.debug.cpp $(debug_srcdepends)
Rcopy_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rcopy_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rcopy_debug__Float128_SOURCES=../common/Rcopy.debug.cpp $(debug_srcdepends)
Rcopy_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rcopy_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rdot_debug__Float128_ref_SOURCES=../common/Rdot.debug.cpp $(debug_srcdepends)
Rdot_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rdot_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rdot_debug__Float128_SOURCES=../common/Rdot.debug.cpp $(debug_srcdepends)
Rdot_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rdot_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgbmv_debug__Float128_ref_SOURCES=../common/Rgbmv.debug.cpp $(debug_srcdepends)
Rgbmv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rgbmv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgbmv_debug__Float128_SOURCES=../common/Rgbmv.debug.cpp $(debug_srcdepends)
Rgbmv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rgbmv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemm_debug__Float128_ref_SOURCES=../common/Rgemm.debug.cpp $(debug_srcdepends)
Rgemm_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rgemm_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemm_debug__Float128_SOURCES=../common/Rgemm.debug.cpp $(debug_srcdepends)
Rgemm_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rgemm_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemv_debug__Float128_ref_SOURCES=../common/Rgemv.debug.cpp $(debug_srcdepends)
Rgemv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rgemv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemv_debug__Float128_SOURCES=../common/Rgemv.debug.cpp $(debug_srcdepends)
Rgemv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rgemv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rger_debug__Float128_ref_SOURCES=../common/Rger.debug.cpp $(debug_srcdepends)
Rger_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rger_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rger_debug__Float128_SOURCES=../common/Rger.debug.cpp $(debug_srcdepends)
Rger_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rger_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rnrm2_debug__Float128_ref_SOURCES=../common/Rnrm2.debug.cpp $(debug_srcdepends)
Rnrm2_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rnrm2_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rnrm2_debug__Float128_SOURCES=../common/Rnrm2.debug.cpp $(debug_srcdepends)
Rnrm2_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rnrm2_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrot_debug__Float128_ref_SOURCES=../common/Rrot.debug.cpp $(debug_srcdepends)
Rrot_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rrot_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrot_debug__Float128_SOURCES=../common/Rrot.debug.cpp $(debug_srcdepends)
Rrot_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rrot_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotg_debug__Float128_ref_SOURCES=../common/Rrotg.debug.cpp $(debug_srcdepends)
Rrotg_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rrotg_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotg_debug__Float128_SOURCES=../common/Rrotg.debug.cpp $(debug_srcdepends)
Rrotg_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rrotg_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotm_debug__Float128_ref_SOURCES=../common/Rrotm.debug.cpp $(debug_srcdepends)
Rrotm_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rrotm_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotm_debug__Float128_SOURCES=../common/Rrotm.debug.cpp $(debug_srcdepends)
Rrotm_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rrotm_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsbmv_debug__Float128_ref_SOURCES=../common/Rsbmv.debug.cpp $(debug_srcdepends)
Rsbmv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rsbmv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsbmv_debug__Float128_SOURCES=../common/Rsbmv.debug.cpp $(debug_srcdepends)
Rsbmv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rsbmv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rscal_debug__Float128_ref_SOURCES=../common/Rscal.debug.cpp $(debug_srcdepends)
Rscal_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rscal_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rscal_debug__Float128_SOURCES=../common/Rscal.debug.cpp $(debug_srcdepends)
Rscal_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rscal_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspmv_debug__Float128_ref_SOURCES=../common/Rspmv.debug.cpp $(debug_srcdepends)
Rspmv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rspmv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspmv_debug__Float128_SOURCES=../common/Rspmv.debug.cpp $(debug_srcdepends)
Rspmv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rspmv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr_debug__Float128_ref_SOURCES=../common/Rspr.debug.cpp $(debug_srcdepends)
Rspr_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rspr_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr_debug__Float128_SOURCES=../common/Rspr.debug.cpp $(debug_srcdepends)
Rspr_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rspr_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr2_debug__Float128_ref_SOURCES=../common/Rspr2.debug.cpp $(debug_srcdepends)
Rspr2_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rspr2_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr2_debug__Float128_SOURCES=../common/Rspr2.debug.cpp $(debug_srcdepends)
Rspr2_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rspr2_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rswap_debug__Float128_ref_SOURCES=../common/Rswap.debug.cpp $(debug_srcdepends)
Rswap_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rswap_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rswap_debug__Float128_SOURCES=../common/Rswap.debug.cpp $(debug_srcdepends)
Rswap_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rswap_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymm_debug__Float128_ref_SOURCES=../common/Rsymm.debug.cpp $(debug_srcdepends)
Rsymm_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rsymm_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymm_debug__Float128_SOURCES=../common/Rsymm.debug.cpp $(debug_srcdepends)
Rsymm_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rsymm_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymv_debug__Float128_ref_SOURCES=../common/Rsymv.debug.cpp $(debug_srcdepends)
Rsymv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rsymv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymv_debug__Float128_SOURCES=../common/Rsymv.debug.cpp $(debug_srcdepends)
Rsymv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rsymv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr_debug__Float128_ref_SOURCES=../common/Rsyr.debug.cpp $(debug_srcdepends)
Rsyr_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr_debug__Float128_SOURCES=../common/Rsyr.debug.cpp $(debug_srcdepends)
Rsyr_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rsyr_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2_debug__Float128_ref_SOURCES=../common/Rsyr2.debug.cpp $(debug_srcdepends)
Rsyr2_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr2_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2_debug__Float128_SOURCES=../common/Rsyr2.debug.cpp $(debug_srcdepends)
Rsyr2_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rsyr2_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2k_debug__Float128_ref_SOURCES=../common/Rsyr2k.debug.cpp $(debug_srcdepends)
Rsyr2k_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr2k_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2k_debug__Float128_SOURCES=../common/Rsyr2k.debug.cpp $(debug_srcdepends)
Rsyr2k_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rsyr2k_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyrk_debug__Float128_ref_SOURCES=../common/Rsyrk.debug.cpp $(debug_srcdepends)
Rsyrk_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rsyrk_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyrk_debug__Float128_SOURCES=../common/Rsyrk.debug.cpp $(debug_srcdepends)
Rsyrk_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rsyrk_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbmv_debug__Float128_ref_SOURCES=../common/Rtbmv.debug.cpp $(debug_srcdepends)
Rtbmv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtbmv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbmv_debug__Float128_SOURCES=../common/Rtbmv.debug.cpp $(debug_srcdepends)
Rtbmv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rtbmv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbsv_debug__Float128_ref_SOURCES=../common/Rtbsv.debug.cpp $(debug_srcdepends)
Rtbsv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtbsv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbsv_debug__Float128_SOURCES=../common/Rtbsv.debug.cpp $(debug_srcdepends)
Rtbsv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rtbsv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpmv_debug__Float128_ref_SOURCES=../common/Rtpmv.debug.cpp $(debug_srcdepends)
Rtpmv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtpmv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpmv_debug__Float128_SOURCES=../common/Rtpmv.debug.cpp $(debug_srcdepends)
Rtpmv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rtpmv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpsv_debug__Float128_ref_SOURCES=../common/Rtpsv.debug.cpp $(debug_srcdepends)
Rtpsv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtpsv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpsv_debug__Float128_SOURCES=../common/Rtpsv.debug.cpp $(debug_srcdepends)
Rtpsv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rtpsv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmm_debug__Float128_ref_SOURCES=../common/Rtrmm.debug.cpp $(debug_srcdepends)
Rtrmm_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtrmm_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmm_debug__Float128_SOURCES=../common/Rtrmm.debug.cpp $(debug_srcdepends)
Rtrmm_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rtrmm_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmv_debug__Float128_ref_SOURCES=../common/Rtrmv.debug.cpp $(debug_srcdepends)
Rtrmv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtrmv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmv_debug__Float128_SOURCES=../common/Rtrmv.debug.cpp $(debug_srcdepends)
Rtrmv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rtrmv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsm_debug__Float128_ref_SOURCES=../common/Rtrsm.debug.cpp $(debug_srcdepends)
Rtrsm_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtrsm_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsm_debug__Float128_SOURCES=../common/Rtrsm.debug.cpp $(debug_srcdepends)
Rtrsm_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rtrsm_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsv_debug__Float128_ref_SOURCES=../common/Rtrsv.debug.cpp $(debug_srcdepends)
Rtrsv_debug__Float128_ref_CXXFLAGS=$(debug_cxxflags)
Rtrsv_debug__Float128_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsv_debug__Float128_SOURCES=../common/Rtrsv.debug.cpp $(debug_srcdepends)
Rtrsv_debug__Float128_CXXFLAGS=$(debug_cxxflags)
Rtrsv_debug__Float128_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)
