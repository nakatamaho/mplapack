#CTBSV, CTPSV might have bugs or intrinsically inaccurate
#Rnrm2 should be checked more
check_PROGRAMS = $(mpblas_gmp_test_PROGRAMS) 
TESTS = $(check_PROGRAMS)
mpblas_gmp_test_PROGRAMS = \
arithmetic.test_gmp \
complex.test_gmp \
iCamax.test.gmp_ref iCamax.test.gmp   \
iRamax.test.gmp_ref iRamax.test.gmp   \
Mlsame.test.gmp_ref Mlsame.test.gmp   \
Mxerbla.test.gmp_ref Mxerbla.test.gmp \
CRrot.test.gmp_ref  CRrot.test.gmp    \
CRscal.test.gmp_ref CRscal.test.gmp   \
Caxpy.test.gmp_ref  Caxpy.test.gmp    \
Ccopy.test.gmp_ref  Ccopy.test.gmp    \
Cdotc.test.gmp_ref  Cdotc.test.gmp	\
Cdotu.test.gmp_ref  Cdotu.test.gmp	\
Cgbmv.test.gmp_ref  Cgbmv.test.gmp	\
Cgemm.test.gmp_ref  Cgemm.test.gmp	\
Cgemv.test.gmp_ref  Cgemv.test.gmp	\
Cgerc.test.gmp_ref  Cgerc.test.gmp	\
Cgeru.test.gmp_ref  Cgeru.test.gmp    \
Chbmv.test.gmp_ref  Chbmv.test.gmp	\
Chemm.test.gmp_ref  Chemm.test.gmp	\
Chemv.test.gmp_ref  Chemv.test.gmp	\
Cher.test.gmp_ref   Cher.test.gmp 	\
Cher2.test.gmp_ref  Cher2.test.gmp	\
Cher2k.test.gmp_ref Cher2k.test.gmp	\
Cherk.test.gmp_ref  Cherk.test.gmp	\
Chpmv.test.gmp_ref  Chpmv.test.gmp	\
Chpr.test.gmp_ref   Chpr.test.gmp	\
Chpr2.test.gmp_ref  Chpr2.test.gmp	\
Crotg.test.gmp_ref  Crotg.test.gmp	\
Cscal.test.gmp_ref  Cscal.test.gmp	\
Cswap.test.gmp_ref  Cswap.test.gmp	\
Csymm.test.gmp_ref  Csymm.test.gmp    \
Csyr2k.test.gmp_ref Csyr2k.test.gmp	\
Csyrk.test.gmp_ref  Csyrk.test.gmp	\
Ctbmv.test.gmp_ref  Ctbmv.test.gmp	\
Ctbsv.test.gmp_ref  Ctbsv.test.gmp	\
Ctpmv.test.gmp_ref  Ctpmv.test.gmp	\
Ctpsv.test.gmp_ref  Ctpsv.test.gmp	\
Ctrmm.test.gmp_ref  Ctrmm.test.gmp	\
Ctrmv.test.gmp_ref  Ctrmv.test.gmp	\
Ctrsm.test.gmp_ref  Ctrsm.test.gmp	\
Ctrsv.test.gmp_ref  Ctrsv.test.gmp	\
RCabs1.test.gmp_ref RCabs1.test.gmp	\
RCasum.test.gmp_ref RCasum.test.gmp	\
RCnrm2.test.gmp_ref RCnrm2.test.gmp	\
Rasum.test.gmp_ref  Rasum.test.gmp    \
Raxpy.test.gmp_ref  Raxpy.test.gmp	\
Rcopy.test.gmp_ref  Rcopy.test.gmp	\
Rdot.test.gmp_ref   Rdot.test.gmp	\
Rgbmv.test.gmp_ref  Rgbmv.test.gmp	\
Rgemm.test.gmp_ref  Rgemm.test.gmp	\
Rgemv.test.gmp_ref  Rgemv.test.gmp	\
Rger.test.gmp_ref   Rger.test.gmp	\
Rnrm2.test.gmp_ref  Rnrm2.test.gmp	\
Rrot.test.gmp_ref   Rrot.test.gmp	\
Rrotg.test.gmp_ref  Rrotg.test.gmp	\
Rrotm.test.gmp_ref  Rrotm.test.gmp	\
Rsbmv.test.gmp_ref  Rsbmv.test.gmp	\
Rscal.test.gmp_ref  Rscal.test.gmp	\
Rspmv.test.gmp_ref  Rspmv.test.gmp    \
Rspr.test.gmp_ref   Rspr.test.gmp	\
Rspr2.test.gmp_ref  Rspr2.test.gmp	\
Rswap.test.gmp_ref  Rswap.test.gmp	\
Rsymm.test.gmp_ref  Rsymm.test.gmp	\
Rsymv.test.gmp_ref  Rsymv.test.gmp	\
Rsyr.test.gmp_ref   Rsyr.test.gmp	\
Rsyr2.test.gmp_ref  Rsyr2.test.gmp	\
Rsyr2k.test.gmp_ref Rsyr2k.test.gmp	\
Rsyrk.test.gmp_ref  Rsyrk.test.gmp	\
Rtbmv.test.gmp_ref  Rtbmv.test.gmp	\
Rtbsv.test.gmp_ref  Rtbsv.test.gmp	\
Rtpmv.test.gmp_ref  Rtpmv.test.gmp	\
Rtpsv.test.gmp_ref  Rtpsv.test.gmp	\
Rtrmm.test.gmp_ref  Rtrmm.test.gmp    \
Rtrmv.test.gmp_ref  Rtrmv.test.gmp	\
Rtrsm.test.gmp_ref  Rtrsm.test.gmp	\
Rtrsv.test.gmp_ref  Rtrsv.test.gmp

mpblas_gmp_testdir = $(prefix)/share/mplapack/test/mpblas/gmp
ref_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_gmp_ref -lmpblas_mpfr_ref
opt_libdepends = -L$(top_builddir)/mpblas/optimized/gmp -lmpblas_gmp -L$(top_builddir)/mpblas/reference -lmpblas_mpfr_ref
mplibs=-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp
refblas=-L$(top_builddir)/external/i/LAPACK/lib -lblas $(FCLIBS)

test_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_GMP___

if !IS_MACOSX
test_srcdepends = ../common/xerbla.cpp ../common/mplapack.test.cpp ../common/Mxerbla.override.cpp
endif

if IS_MACOSX
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_GMP___
libmxerbla_override_la_SOURCES =  ../common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS= -Wl,-flat_namespace,-undefined,dynamic_lookup
test_srcdepends = ../common/xerbla.cpp ../common/mplapack.test.cpp
mplibs+=-Wl,-flat_namespace,-undefined,dynamic_lookup
TESTS_ENVIRONMENT=DYLD_INSERT_LIBRARIES=.libs/libmxerbla_override.0.dylib
endif

##
arithmetic_test_gmp_SOURCES=../common/arithmetic.test.cpp $(test_srcdepends)
arithmetic_test_gmp_CXXFLAGS=$(test_cxxflags)
arithmetic_test_gmp_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

complex_test_gmp_SOURCES=../common/complex.test.cpp $(test_srcdepends)
complex_test_gmp_CXXFLAGS=$(test_cxxflags)
complex_test_gmp_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_test_gmp_ref_SOURCES=../common/iCamax.test.cpp $(test_srcdepends)
iCamax_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
iCamax_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_test_gmp_SOURCES=../common/iCamax.test.cpp $(test_srcdepends)
iCamax_test_gmp_CXXFLAGS=$(test_cxxflags)
iCamax_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iRamax_test_gmp_ref_SOURCES=../common/iRamax.test.cpp $(test_srcdepends)
iRamax_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
iRamax_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iRamax_test_gmp_SOURCES=../common/iRamax.test.cpp $(test_srcdepends)
iRamax_test_gmp_CXXFLAGS=$(test_cxxflags)
iRamax_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mlsame_test_gmp_ref_SOURCES=../common/Mlsame.test.cpp $(test_srcdepends)
Mlsame_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Mlsame_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mlsame_test_gmp_SOURCES=../common/Mlsame.test.cpp $(test_srcdepends)
Mlsame_test_gmp_CXXFLAGS=$(test_cxxflags)
Mlsame_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mxerbla_test_gmp_ref_SOURCES=../common/Mxerbla.test.cpp $(test_srcdepends)
Mxerbla_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Mxerbla_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mxerbla_test_gmp_SOURCES=../common/Mxerbla.test.cpp $(test_srcdepends)
Mxerbla_test_gmp_CXXFLAGS=$(test_cxxflags)
Mxerbla_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRrot_test_gmp_ref_SOURCES=../common/CRrot.test.cpp $(test_srcdepends)
CRrot_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
CRrot_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRrot_test_gmp_SOURCES=../common/CRrot.test.cpp $(test_srcdepends)
CRrot_test_gmp_CXXFLAGS=$(test_cxxflags)
CRrot_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRscal_test_gmp_ref_SOURCES=../common/CRscal.test.cpp $(test_srcdepends)
CRscal_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
CRscal_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRscal_test_gmp_SOURCES=../common/CRscal.test.cpp $(test_srcdepends)
CRscal_test_gmp_CXXFLAGS=$(test_cxxflags)
CRscal_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_test_gmp_SOURCES=../common/Caxpy.test.cpp $(test_srcdepends)
Caxpy_test_gmp_CXXFLAGS=$(test_cxxflags)
Caxpy_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_test_gmp_ref_SOURCES=../common/Caxpy.test.cpp $(test_srcdepends)
Caxpy_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Caxpy_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ccopy_test_gmp_SOURCES=../common/Ccopy.test.cpp $(test_srcdepends)
Ccopy_test_gmp_CXXFLAGS=$(test_cxxflags)
Ccopy_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ccopy_test_gmp_ref_SOURCES=../common/Ccopy.test.cpp $(test_srcdepends)
Ccopy_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Ccopy_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotc_test_gmp_SOURCES=../common/Cdotc.test.cpp $(test_srcdepends)
Cdotc_test_gmp_CXXFLAGS=$(test_cxxflags)
Cdotc_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotc_test_gmp_ref_SOURCES=../common/Cdotc.test.cpp $(test_srcdepends)
Cdotc_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Cdotc_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotu_test_gmp_SOURCES=../common/Cdotu.test.cpp $(test_srcdepends)
Cdotu_test_gmp_CXXFLAGS=$(test_cxxflags)
Cdotu_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotu_test_gmp_ref_SOURCES=../common/Cdotu.test.cpp $(test_srcdepends)
Cdotu_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Cdotu_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgbmv_test_gmp_SOURCES=../common/Cgbmv.test.cpp $(test_srcdepends)
Cgbmv_test_gmp_CXXFLAGS=$(test_cxxflags)
Cgbmv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgbmv_test_gmp_ref_SOURCES=../common/Cgbmv.test.cpp $(test_srcdepends)
Cgbmv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Cgbmv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemm_test_gmp_SOURCES=../common/Cgemm.test.cpp $(test_srcdepends)
Cgemm_test_gmp_CXXFLAGS=$(test_cxxflags)
Cgemm_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemm_test_gmp_ref_SOURCES=../common/Cgemm.test.cpp $(test_srcdepends)
Cgemm_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Cgemm_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemv_test_gmp_SOURCES=../common/Cgemv.test.cpp $(test_srcdepends)
Cgemv_test_gmp_CXXFLAGS=$(test_cxxflags)
Cgemv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemv_test_gmp_ref_SOURCES=../common/Cgemv.test.cpp $(test_srcdepends)
Cgemv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Cgemv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgerc_test_gmp_SOURCES=../common/Cgerc.test.cpp $(test_srcdepends)
Cgerc_test_gmp_CXXFLAGS=$(test_cxxflags)
Cgerc_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgerc_test_gmp_ref_SOURCES=../common/Cgerc.test.cpp $(test_srcdepends)
Cgerc_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Cgerc_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgeru_test_gmp_SOURCES=../common/Cgeru.test.cpp $(test_srcdepends)
Cgeru_test_gmp_CXXFLAGS=$(test_cxxflags)
Cgeru_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgeru_test_gmp_ref_SOURCES=../common/Cgeru.test.cpp $(test_srcdepends)
Cgeru_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Cgeru_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chbmv_test_gmp_SOURCES=../common/Chbmv.test.cpp $(test_srcdepends)
Chbmv_test_gmp_CXXFLAGS=$(test_cxxflags)
Chbmv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chbmv_test_gmp_ref_SOURCES=../common/Chbmv.test.cpp $(test_srcdepends)
Chbmv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Chbmv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemm_test_gmp_SOURCES=../common/Chemm.test.cpp $(test_srcdepends)
Chemm_test_gmp_CXXFLAGS=$(test_cxxflags)
Chemm_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemm_test_gmp_ref_SOURCES=../common/Chemm.test.cpp $(test_srcdepends)
Chemm_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Chemm_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemv_test_gmp_SOURCES=../common/Chemv.test.cpp $(test_srcdepends)
Chemv_test_gmp_CXXFLAGS=$(test_cxxflags)
Chemv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemv_test_gmp_ref_SOURCES=../common/Chemv.test.cpp $(test_srcdepends)
Chemv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Chemv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher_test_gmp_SOURCES=../common/Cher.test.cpp $(test_srcdepends)
Cher_test_gmp_CXXFLAGS=$(test_cxxflags)
Cher_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher_test_gmp_ref_SOURCES=../common/Cher.test.cpp $(test_srcdepends)
Cher_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Cher_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2_test_gmp_SOURCES=../common/Cher2.test.cpp $(test_srcdepends)
Cher2_test_gmp_CXXFLAGS=$(test_cxxflags)
Cher2_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2_test_gmp_ref_SOURCES=../common/Cher2.test.cpp $(test_srcdepends)
Cher2_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Cher2_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2k_test_gmp_SOURCES=../common/Cher2k.test.cpp $(test_srcdepends)
Cher2k_test_gmp_CXXFLAGS=$(test_cxxflags)
Cher2k_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2k_test_gmp_ref_SOURCES=../common/Cher2k.test.cpp $(test_srcdepends)
Cher2k_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Cher2k_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cherk_test_gmp_SOURCES=../common/Cherk.test.cpp $(test_srcdepends)
Cherk_test_gmp_CXXFLAGS=$(test_cxxflags)
Cherk_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cherk_test_gmp_ref_SOURCES=../common/Cherk.test.cpp $(test_srcdepends)
Cherk_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Cherk_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpmv_test_gmp_SOURCES=../common/Chpmv.test.cpp $(test_srcdepends)
Chpmv_test_gmp_CXXFLAGS=$(test_cxxflags)
Chpmv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpmv_test_gmp_ref_SOURCES=../common/Chpmv.test.cpp $(test_srcdepends)
Chpmv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Chpmv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr_test_gmp_SOURCES=../common/Chpr.test.cpp $(test_srcdepends)
Chpr_test_gmp_CXXFLAGS=$(test_cxxflags)
Chpr_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr_test_gmp_ref_SOURCES=../common/Chpr.test.cpp $(test_srcdepends)
Chpr_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Chpr_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr2_test_gmp_SOURCES=../common/Chpr2.test.cpp $(test_srcdepends)
Chpr2_test_gmp_CXXFLAGS=$(test_cxxflags)
Chpr2_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr2_test_gmp_ref_SOURCES=../common/Chpr2.test.cpp $(test_srcdepends)
Chpr2_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Chpr2_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Crotg_test_gmp_SOURCES=../common/Crotg.test.cpp $(test_srcdepends)
Crotg_test_gmp_CXXFLAGS=$(test_cxxflags)
Crotg_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Crotg_test_gmp_ref_SOURCES=../common/Crotg.test.cpp $(test_srcdepends)
Crotg_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Crotg_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cscal_test_gmp_SOURCES=../common/Cscal.test.cpp $(test_srcdepends)
Cscal_test_gmp_CXXFLAGS=$(test_cxxflags)
Cscal_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cscal_test_gmp_ref_SOURCES=../common/Cscal.test.cpp $(test_srcdepends)
Cscal_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Cscal_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cswap_test_gmp_SOURCES=../common/Cswap.test.cpp $(test_srcdepends)
Cswap_test_gmp_CXXFLAGS=$(test_cxxflags)
Cswap_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cswap_test_gmp_ref_SOURCES=../common/Cswap.test.cpp $(test_srcdepends)
Cswap_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Cswap_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csymm_test_gmp_SOURCES=../common/Csymm.test.cpp $(test_srcdepends)
Csymm_test_gmp_CXXFLAGS=$(test_cxxflags)
Csymm_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csymm_test_gmp_ref_SOURCES=../common/Csymm.test.cpp $(test_srcdepends)
Csymm_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Csymm_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyr2k_test_gmp_SOURCES=../common/Csyr2k.test.cpp $(test_srcdepends)
Csyr2k_test_gmp_CXXFLAGS=$(test_cxxflags)
Csyr2k_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyr2k_test_gmp_ref_SOURCES=../common/Csyr2k.test.cpp $(test_srcdepends)
Csyr2k_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Csyr2k_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyrk_test_gmp_SOURCES=../common/Csyrk.test.cpp $(test_srcdepends)
Csyrk_test_gmp_CXXFLAGS=$(test_cxxflags)
Csyrk_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyrk_test_gmp_ref_SOURCES=../common/Csyrk.test.cpp $(test_srcdepends)
Csyrk_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Csyrk_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbmv_test_gmp_SOURCES=../common/Ctbmv.test.cpp $(test_srcdepends)
Ctbmv_test_gmp_CXXFLAGS=$(test_cxxflags)
Ctbmv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbmv_test_gmp_ref_SOURCES=../common/Ctbmv.test.cpp $(test_srcdepends)
Ctbmv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Ctbmv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbsv_test_gmp_SOURCES=../common/Ctbsv.test.cpp $(test_srcdepends)
Ctbsv_test_gmp_CXXFLAGS=$(test_cxxflags)
Ctbsv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbsv_test_gmp_ref_SOURCES=../common/Ctbsv.test.cpp $(test_srcdepends)
Ctbsv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Ctbsv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpmv_test_gmp_SOURCES=../common/Ctpmv.test.cpp $(test_srcdepends)
Ctpmv_test_gmp_CXXFLAGS=$(test_cxxflags)
Ctpmv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpmv_test_gmp_ref_SOURCES=../common/Ctpmv.test.cpp $(test_srcdepends)
Ctpmv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Ctpmv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpsv_test_gmp_SOURCES=../common/Ctpsv.test.cpp $(test_srcdepends)
Ctpsv_test_gmp_CXXFLAGS=$(test_cxxflags)
Ctpsv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpsv_test_gmp_ref_SOURCES=../common/Ctpsv.test.cpp $(test_srcdepends)
Ctpsv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Ctpsv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmm_test_gmp_SOURCES=../common/Ctrmm.test.cpp $(test_srcdepends)
Ctrmm_test_gmp_CXXFLAGS=$(test_cxxflags)
Ctrmm_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmm_test_gmp_ref_SOURCES=../common/Ctrmm.test.cpp $(test_srcdepends)
Ctrmm_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Ctrmm_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmv_test_gmp_SOURCES=../common/Ctrmv.test.cpp $(test_srcdepends)
Ctrmv_test_gmp_CXXFLAGS=$(test_cxxflags)
Ctrmv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmv_test_gmp_ref_SOURCES=../common/Ctrmv.test.cpp $(test_srcdepends)
Ctrmv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Ctrmv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsm_test_gmp_SOURCES=../common/Ctrsm.test.cpp $(test_srcdepends)
Ctrsm_test_gmp_CXXFLAGS=$(test_cxxflags)
Ctrsm_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsm_test_gmp_ref_SOURCES=../common/Ctrsm.test.cpp $(test_srcdepends)
Ctrsm_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Ctrsm_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsv_test_gmp_SOURCES=../common/Ctrsv.test.cpp $(test_srcdepends)
Ctrsv_test_gmp_CXXFLAGS=$(test_cxxflags)
Ctrsv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsv_test_gmp_ref_SOURCES=../common/Ctrsv.test.cpp $(test_srcdepends)
Ctrsv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Ctrsv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_test_gmp_ref_SOURCES=../common/RCabs1.test.cpp $(test_srcdepends)
RCabs1_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
RCabs1_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_test_gmp_SOURCES=../common/RCabs1.test.cpp $(test_srcdepends)
RCabs1_test_gmp_CXXFLAGS=$(test_cxxflags)
RCabs1_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCasum_test_gmp_ref_SOURCES=../common/RCasum.test.cpp $(test_srcdepends)
RCasum_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
RCasum_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCasum_test_gmp_SOURCES=../common/RCasum.test.cpp $(test_srcdepends)
RCasum_test_gmp_CXXFLAGS=$(test_cxxflags)
RCasum_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCnrm2_test_gmp_ref_SOURCES=../common/RCnrm2.test.cpp $(test_srcdepends)
RCnrm2_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
RCnrm2_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCnrm2_test_gmp_SOURCES=../common/RCnrm2.test.cpp $(test_srcdepends)
RCnrm2_test_gmp_CXXFLAGS=$(test_cxxflags)
RCnrm2_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rasum_test_gmp_ref_SOURCES=../common/Rasum.test.cpp $(test_srcdepends)
Rasum_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rasum_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rasum_test_gmp_SOURCES=../common/Rasum.test.cpp $(test_srcdepends)
Rasum_test_gmp_CXXFLAGS=$(test_cxxflags)
Rasum_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Raxpy_test_gmp_ref_SOURCES=../common/Raxpy.test.cpp $(test_srcdepends)
Raxpy_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Raxpy_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Raxpy_test_gmp_SOURCES=../common/Raxpy.test.cpp $(test_srcdepends)
Raxpy_test_gmp_CXXFLAGS=$(test_cxxflags)
Raxpy_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rcopy_test_gmp_ref_SOURCES=../common/Rcopy.test.cpp $(test_srcdepends)
Rcopy_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rcopy_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rcopy_test_gmp_SOURCES=../common/Rcopy.test.cpp $(test_srcdepends)
Rcopy_test_gmp_CXXFLAGS=$(test_cxxflags)
Rcopy_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rdot_test_gmp_ref_SOURCES=../common/Rdot.test.cpp $(test_srcdepends)
Rdot_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rdot_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rdot_test_gmp_SOURCES=../common/Rdot.test.cpp $(test_srcdepends)
Rdot_test_gmp_CXXFLAGS=$(test_cxxflags)
Rdot_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgbmv_test_gmp_ref_SOURCES=../common/Rgbmv.test.cpp $(test_srcdepends)
Rgbmv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rgbmv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgbmv_test_gmp_SOURCES=../common/Rgbmv.test.cpp $(test_srcdepends)
Rgbmv_test_gmp_CXXFLAGS=$(test_cxxflags)
Rgbmv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemm_test_gmp_ref_SOURCES=../common/Rgemm.test.cpp $(test_srcdepends)
Rgemm_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rgemm_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemm_test_gmp_SOURCES=../common/Rgemm.test.cpp $(test_srcdepends)
Rgemm_test_gmp_CXXFLAGS=$(test_cxxflags)
Rgemm_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemv_test_gmp_ref_SOURCES=../common/Rgemv.test.cpp $(test_srcdepends)
Rgemv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rgemv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemv_test_gmp_SOURCES=../common/Rgemv.test.cpp $(test_srcdepends)
Rgemv_test_gmp_CXXFLAGS=$(test_cxxflags)
Rgemv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rger_test_gmp_ref_SOURCES=../common/Rger.test.cpp $(test_srcdepends)
Rger_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rger_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rger_test_gmp_SOURCES=../common/Rger.test.cpp $(test_srcdepends)
Rger_test_gmp_CXXFLAGS=$(test_cxxflags)
Rger_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rnrm2_test_gmp_ref_SOURCES=../common/Rnrm2.test.cpp $(test_srcdepends)
Rnrm2_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rnrm2_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rnrm2_test_gmp_SOURCES=../common/Rnrm2.test.cpp $(test_srcdepends)
Rnrm2_test_gmp_CXXFLAGS=$(test_cxxflags)
Rnrm2_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrot_test_gmp_ref_SOURCES=../common/Rrot.test.cpp $(test_srcdepends)
Rrot_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rrot_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrot_test_gmp_SOURCES=../common/Rrot.test.cpp $(test_srcdepends)
Rrot_test_gmp_CXXFLAGS=$(test_cxxflags)
Rrot_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotg_test_gmp_ref_SOURCES=../common/Rrotg.test.cpp $(test_srcdepends)
Rrotg_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rrotg_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotg_test_gmp_SOURCES=../common/Rrotg.test.cpp $(test_srcdepends)
Rrotg_test_gmp_CXXFLAGS=$(test_cxxflags)
Rrotg_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotm_test_gmp_ref_SOURCES=../common/Rrotm.test.cpp $(test_srcdepends)
Rrotm_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rrotm_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotm_test_gmp_SOURCES=../common/Rrotm.test.cpp $(test_srcdepends)
Rrotm_test_gmp_CXXFLAGS=$(test_cxxflags)
Rrotm_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsbmv_test_gmp_ref_SOURCES=../common/Rsbmv.test.cpp $(test_srcdepends)
Rsbmv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rsbmv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsbmv_test_gmp_SOURCES=../common/Rsbmv.test.cpp $(test_srcdepends)
Rsbmv_test_gmp_CXXFLAGS=$(test_cxxflags)
Rsbmv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rscal_test_gmp_ref_SOURCES=../common/Rscal.test.cpp $(test_srcdepends)
Rscal_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rscal_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rscal_test_gmp_SOURCES=../common/Rscal.test.cpp $(test_srcdepends)
Rscal_test_gmp_CXXFLAGS=$(test_cxxflags)
Rscal_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspmv_test_gmp_ref_SOURCES=../common/Rspmv.test.cpp $(test_srcdepends)
Rspmv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rspmv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspmv_test_gmp_SOURCES=../common/Rspmv.test.cpp $(test_srcdepends)
Rspmv_test_gmp_CXXFLAGS=$(test_cxxflags)
Rspmv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr_test_gmp_ref_SOURCES=../common/Rspr.test.cpp $(test_srcdepends)
Rspr_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rspr_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr_test_gmp_SOURCES=../common/Rspr.test.cpp $(test_srcdepends)
Rspr_test_gmp_CXXFLAGS=$(test_cxxflags)
Rspr_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr2_test_gmp_ref_SOURCES=../common/Rspr2.test.cpp $(test_srcdepends)
Rspr2_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rspr2_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr2_test_gmp_SOURCES=../common/Rspr2.test.cpp $(test_srcdepends)
Rspr2_test_gmp_CXXFLAGS=$(test_cxxflags)
Rspr2_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rswap_test_gmp_ref_SOURCES=../common/Rswap.test.cpp $(test_srcdepends)
Rswap_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rswap_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rswap_test_gmp_SOURCES=../common/Rswap.test.cpp $(test_srcdepends)
Rswap_test_gmp_CXXFLAGS=$(test_cxxflags)
Rswap_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymm_test_gmp_ref_SOURCES=../common/Rsymm.test.cpp $(test_srcdepends)
Rsymm_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rsymm_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymm_test_gmp_SOURCES=../common/Rsymm.test.cpp $(test_srcdepends)
Rsymm_test_gmp_CXXFLAGS=$(test_cxxflags)
Rsymm_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymv_test_gmp_ref_SOURCES=../common/Rsymv.test.cpp $(test_srcdepends)
Rsymv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rsymv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymv_test_gmp_SOURCES=../common/Rsymv.test.cpp $(test_srcdepends)
Rsymv_test_gmp_CXXFLAGS=$(test_cxxflags)
Rsymv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr_test_gmp_ref_SOURCES=../common/Rsyr.test.cpp $(test_srcdepends)
Rsyr_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rsyr_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr_test_gmp_SOURCES=../common/Rsyr.test.cpp $(test_srcdepends)
Rsyr_test_gmp_CXXFLAGS=$(test_cxxflags)
Rsyr_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2_test_gmp_ref_SOURCES=../common/Rsyr2.test.cpp $(test_srcdepends)
Rsyr2_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rsyr2_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2_test_gmp_SOURCES=../common/Rsyr2.test.cpp $(test_srcdepends)
Rsyr2_test_gmp_CXXFLAGS=$(test_cxxflags)
Rsyr2_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2k_test_gmp_ref_SOURCES=../common/Rsyr2k.test.cpp $(test_srcdepends)
Rsyr2k_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rsyr2k_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2k_test_gmp_SOURCES=../common/Rsyr2k.test.cpp $(test_srcdepends)
Rsyr2k_test_gmp_CXXFLAGS=$(test_cxxflags)
Rsyr2k_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyrk_test_gmp_ref_SOURCES=../common/Rsyrk.test.cpp $(test_srcdepends)
Rsyrk_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rsyrk_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyrk_test_gmp_SOURCES=../common/Rsyrk.test.cpp $(test_srcdepends)
Rsyrk_test_gmp_CXXFLAGS=$(test_cxxflags)
Rsyrk_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbmv_test_gmp_ref_SOURCES=../common/Rtbmv.test.cpp $(test_srcdepends)
Rtbmv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rtbmv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbmv_test_gmp_SOURCES=../common/Rtbmv.test.cpp $(test_srcdepends)
Rtbmv_test_gmp_CXXFLAGS=$(test_cxxflags)
Rtbmv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbsv_test_gmp_ref_SOURCES=../common/Rtbsv.test.cpp $(test_srcdepends)
Rtbsv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rtbsv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbsv_test_gmp_SOURCES=../common/Rtbsv.test.cpp $(test_srcdepends)
Rtbsv_test_gmp_CXXFLAGS=$(test_cxxflags)
Rtbsv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpmv_test_gmp_ref_SOURCES=../common/Rtpmv.test.cpp $(test_srcdepends)
Rtpmv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rtpmv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpmv_test_gmp_SOURCES=../common/Rtpmv.test.cpp $(test_srcdepends)
Rtpmv_test_gmp_CXXFLAGS=$(test_cxxflags)
Rtpmv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpsv_test_gmp_ref_SOURCES=../common/Rtpsv.test.cpp $(test_srcdepends)
Rtpsv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rtpsv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpsv_test_gmp_SOURCES=../common/Rtpsv.test.cpp $(test_srcdepends)
Rtpsv_test_gmp_CXXFLAGS=$(test_cxxflags)
Rtpsv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmm_test_gmp_ref_SOURCES=../common/Rtrmm.test.cpp $(test_srcdepends)
Rtrmm_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rtrmm_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmm_test_gmp_SOURCES=../common/Rtrmm.test.cpp $(test_srcdepends)
Rtrmm_test_gmp_CXXFLAGS=$(test_cxxflags)
Rtrmm_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmv_test_gmp_ref_SOURCES=../common/Rtrmv.test.cpp $(test_srcdepends)
Rtrmv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rtrmv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmv_test_gmp_SOURCES=../common/Rtrmv.test.cpp $(test_srcdepends)
Rtrmv_test_gmp_CXXFLAGS=$(test_cxxflags)
Rtrmv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsm_test_gmp_ref_SOURCES=../common/Rtrsm.test.cpp $(test_srcdepends)
Rtrsm_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rtrsm_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsm_test_gmp_SOURCES=../common/Rtrsm.test.cpp $(test_srcdepends)
Rtrsm_test_gmp_CXXFLAGS=$(test_cxxflags)
Rtrsm_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsv_test_gmp_ref_SOURCES=../common/Rtrsv.test.cpp $(test_srcdepends)
Rtrsv_test_gmp_ref_CXXFLAGS=$(test_cxxflags)
Rtrsv_test_gmp_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsv_test_gmp_SOURCES=../common/Rtrsv.test.cpp $(test_srcdepends)
Rtrsv_test_gmp_CXXFLAGS=$(test_cxxflags)
Rtrsv_test_gmp_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)
