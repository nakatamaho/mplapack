#CTBSV, CTPSV might have bugs or intrinsically inaccurate
#Rnrm2 should be checked more
check_PROGRAMS = $(mpblas__Float128_test_PROGRAMS) 
TESTS = $(check_PROGRAMS)
mpblas__Float128_test_PROGRAMS = \
arithmetic.test__Float128 \
complex.test__Float128 \
iCamax.test._Float128_opt iCamax.test._Float128   \
iRamax.test._Float128_opt iRamax.test._Float128   \
Mlsame.test._Float128_opt Mlsame.test._Float128   \
Mxerbla.test._Float128_opt Mxerbla.test._Float128 \
CRrot.test._Float128_opt  CRrot.test._Float128    \
CRscal.test._Float128_opt CRscal.test._Float128   \
Caxpy.test._Float128_opt  Caxpy.test._Float128    \
Ccopy.test._Float128_opt  Ccopy.test._Float128    \
Cdotc.test._Float128_opt  Cdotc.test._Float128	\
Cdotu.test._Float128_opt  Cdotu.test._Float128	\
Cgbmv.test._Float128_opt  Cgbmv.test._Float128	\
Cgemm.test._Float128_opt  Cgemm.test._Float128	\
Cgemv.test._Float128_opt  Cgemv.test._Float128	\
Cgerc.test._Float128_opt  Cgerc.test._Float128	\
Cgeru.test._Float128_opt  Cgeru.test._Float128    \
Chbmv.test._Float128_opt  Chbmv.test._Float128	\
Chemm.test._Float128_opt  Chemm.test._Float128	\
Chemv.test._Float128_opt  Chemv.test._Float128	\
Cher.test._Float128_opt   Cher.test._Float128 	\
Cher2.test._Float128_opt  Cher2.test._Float128	\
Cher2k.test._Float128_opt Cher2k.test._Float128	\
Cherk.test._Float128_opt  Cherk.test._Float128	\
Chpmv.test._Float128_opt  Chpmv.test._Float128	\
Chpr.test._Float128_opt   Chpr.test._Float128	\
Chpr2.test._Float128_opt  Chpr2.test._Float128	\
Crotg.test._Float128_opt  Crotg.test._Float128	\
Cscal.test._Float128_opt  Cscal.test._Float128	\
Cswap.test._Float128_opt  Cswap.test._Float128	\
Csymm.test._Float128_opt  Csymm.test._Float128    \
Csyr2k.test._Float128_opt Csyr2k.test._Float128	\
Csyrk.test._Float128_opt  Csyrk.test._Float128	\
Ctbmv.test._Float128_opt  Ctbmv.test._Float128	\
Ctbsv.test._Float128_opt  Ctbsv.test._Float128	\
Ctpmv.test._Float128_opt  Ctpmv.test._Float128	\
Ctpsv.test._Float128_opt  Ctpsv.test._Float128	\
Ctrmm.test._Float128_opt  Ctrmm.test._Float128	\
Ctrmv.test._Float128_opt  Ctrmv.test._Float128	\
Ctrsm.test._Float128_opt  Ctrsm.test._Float128	\
Ctrsv.test._Float128_opt  Ctrsv.test._Float128	\
RCabs1.test._Float128_opt RCabs1.test._Float128	\
RCasum.test._Float128_opt RCasum.test._Float128	\
RCnrm2.test._Float128_opt RCnrm2.test._Float128	\
Rasum.test._Float128_opt  Rasum.test._Float128    \
Raxpy.test._Float128_opt  Raxpy.test._Float128	\
Rcopy.test._Float128_opt  Rcopy.test._Float128	\
Rdot.test._Float128_opt   Rdot.test._Float128	\
Rgbmv.test._Float128_opt  Rgbmv.test._Float128	\
Rgemm.test._Float128_opt  Rgemm.test._Float128	\
Rgemv.test._Float128_opt  Rgemv.test._Float128	\
Rger.test._Float128_opt   Rger.test._Float128	\
Rnrm2.test._Float128_opt  Rnrm2.test._Float128	\
Rrot.test._Float128_opt   Rrot.test._Float128	\
Rrotg.test._Float128_opt  Rrotg.test._Float128	\
Rrotm.test._Float128_opt  Rrotm.test._Float128	\
Rsbmv.test._Float128_opt  Rsbmv.test._Float128	\
Rscal.test._Float128_opt  Rscal.test._Float128	\
Rspmv.test._Float128_opt  Rspmv.test._Float128    \
Rspr.test._Float128_opt   Rspr.test._Float128	\
Rspr2.test._Float128_opt  Rspr2.test._Float128	\
Rswap.test._Float128_opt  Rswap.test._Float128	\
Rsymm.test._Float128_opt  Rsymm.test._Float128	\
Rsymv.test._Float128_opt  Rsymv.test._Float128	\
Rsyr.test._Float128_opt   Rsyr.test._Float128	\
Rsyr2.test._Float128_opt  Rsyr2.test._Float128	\
Rsyr2k.test._Float128_opt Rsyr2k.test._Float128	\
Rsyrk.test._Float128_opt  Rsyrk.test._Float128	\
Rtbmv.test._Float128_opt  Rtbmv.test._Float128	\
Rtbsv.test._Float128_opt  Rtbsv.test._Float128	\
Rtpmv.test._Float128_opt  Rtpmv.test._Float128	\
Rtpsv.test._Float128_opt  Rtpsv.test._Float128	\
Rtrmm.test._Float128_opt  Rtrmm.test._Float128    \
Rtrmv.test._Float128_opt  Rtrmv.test._Float128	\
Rtrsm.test._Float128_opt  Rtrsm.test._Float128	\
Rtrsv.test._Float128_opt  Rtrsv.test._Float128

mpblas__Float128_testdir = $(prefix)/lib/$(target)/mplapack/test/compare/_Float128

libdepends = -L$(top_builddir)/mpblas/reference -lmpblas__Float128 -lmpblas_mpfr
opt_libdepends = -L$(top_builddir)/mpblas/optimized/_Float128 -lmpblas__Float128_opt -L$(top_builddir)/mpblas/reference -lmpblas_mpfr
mplibs=-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp
refblas=-L$(top_builddir)/external/i/LAPACK/lib -lblas $(FCLIBS)

test_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH__FLOAT128___ -D___MPLAPACK_INTERNAL___

if !IS_MACOS
test_srcdepends = ../common/xerbla.cpp ../common/mplapack.test.cpp ../common/Mxerbla.override.cpp
endif

if IS_MACOS
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT128___
libmxerbla_override_la_SOURCES =  ../common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS= -Wl,-flat_namespace,-undefined,dynamic_lookup
test_srcdepends = ../common/xerbla.cpp ../common/mplapack.test.cpp
mplibs+=-Wl,-flat_namespace,-undefined,dynamic_lookup
TESTS_ENVIRONMENT=DYLD_INSERT_LIBRARIES=.libs/libmxerbla_override.0.dylib
endif

##
arithmetic_test__Float128_SOURCES=../common/arithmetic.test.cpp $(test_srcdepends)
arithmetic_test__Float128_CXXFLAGS=$(test_cxxflags)
arithmetic_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

complex_test__Float128_SOURCES=../common/complex.test.cpp $(test_srcdepends)
complex_test__Float128_CXXFLAGS=$(test_cxxflags)
complex_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

iCamax_test__Float128_opt_SOURCES=../common/iCamax.test.cpp $(test_srcdepends)
iCamax_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
iCamax_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iCamax_test__Float128_SOURCES=../common/iCamax.test.cpp $(test_srcdepends)
iCamax_test__Float128_CXXFLAGS=$(test_cxxflags)
iCamax_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

iRamax_test__Float128_opt_SOURCES=../common/iRamax.test.cpp $(test_srcdepends)
iRamax_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
iRamax_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iRamax_test__Float128_SOURCES=../common/iRamax.test.cpp $(test_srcdepends)
iRamax_test__Float128_CXXFLAGS=$(test_cxxflags)
iRamax_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Mlsame_test__Float128_opt_SOURCES=../common/Mlsame.test.cpp $(test_srcdepends)
Mlsame_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Mlsame_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mlsame_test__Float128_SOURCES=../common/Mlsame.test.cpp $(test_srcdepends)
Mlsame_test__Float128_CXXFLAGS=$(test_cxxflags)
Mlsame_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Mxerbla_test__Float128_opt_SOURCES=../common/Mxerbla.test.cpp $(test_srcdepends)
Mxerbla_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Mxerbla_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mxerbla_test__Float128_SOURCES=../common/Mxerbla.test.cpp $(test_srcdepends)
Mxerbla_test__Float128_CXXFLAGS=$(test_cxxflags)
Mxerbla_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

CRrot_test__Float128_opt_SOURCES=../common/CRrot.test.cpp $(test_srcdepends)
CRrot_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
CRrot_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRrot_test__Float128_SOURCES=../common/CRrot.test.cpp $(test_srcdepends)
CRrot_test__Float128_CXXFLAGS=$(test_cxxflags)
CRrot_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

CRscal_test__Float128_opt_SOURCES=../common/CRscal.test.cpp $(test_srcdepends)
CRscal_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
CRscal_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRscal_test__Float128_SOURCES=../common/CRscal.test.cpp $(test_srcdepends)
CRscal_test__Float128_CXXFLAGS=$(test_cxxflags)
CRscal_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Caxpy_test__Float128_SOURCES=../common/Caxpy.test.cpp $(test_srcdepends)
Caxpy_test__Float128_CXXFLAGS=$(test_cxxflags)
Caxpy_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Caxpy_test__Float128_opt_SOURCES=../common/Caxpy.test.cpp $(test_srcdepends)
Caxpy_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Caxpy_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ccopy_test__Float128_SOURCES=../common/Ccopy.test.cpp $(test_srcdepends)
Ccopy_test__Float128_CXXFLAGS=$(test_cxxflags)
Ccopy_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ccopy_test__Float128_opt_SOURCES=../common/Ccopy.test.cpp $(test_srcdepends)
Ccopy_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Ccopy_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotc_test__Float128_SOURCES=../common/Cdotc.test.cpp $(test_srcdepends)
Cdotc_test__Float128_CXXFLAGS=$(test_cxxflags)
Cdotc_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cdotc_test__Float128_opt_SOURCES=../common/Cdotc.test.cpp $(test_srcdepends)
Cdotc_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Cdotc_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotu_test__Float128_SOURCES=../common/Cdotu.test.cpp $(test_srcdepends)
Cdotu_test__Float128_CXXFLAGS=$(test_cxxflags)
Cdotu_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cdotu_test__Float128_opt_SOURCES=../common/Cdotu.test.cpp $(test_srcdepends)
Cdotu_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Cdotu_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgbmv_test__Float128_SOURCES=../common/Cgbmv.test.cpp $(test_srcdepends)
Cgbmv_test__Float128_CXXFLAGS=$(test_cxxflags)
Cgbmv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgbmv_test__Float128_opt_SOURCES=../common/Cgbmv.test.cpp $(test_srcdepends)
Cgbmv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Cgbmv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemm_test__Float128_SOURCES=../common/Cgemm.test.cpp $(test_srcdepends)
Cgemm_test__Float128_CXXFLAGS=$(test_cxxflags)
Cgemm_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgemm_test__Float128_opt_SOURCES=../common/Cgemm.test.cpp $(test_srcdepends)
Cgemm_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Cgemm_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemv_test__Float128_SOURCES=../common/Cgemv.test.cpp $(test_srcdepends)
Cgemv_test__Float128_CXXFLAGS=$(test_cxxflags)
Cgemv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgemv_test__Float128_opt_SOURCES=../common/Cgemv.test.cpp $(test_srcdepends)
Cgemv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Cgemv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgerc_test__Float128_SOURCES=../common/Cgerc.test.cpp $(test_srcdepends)
Cgerc_test__Float128_CXXFLAGS=$(test_cxxflags)
Cgerc_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgerc_test__Float128_opt_SOURCES=../common/Cgerc.test.cpp $(test_srcdepends)
Cgerc_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Cgerc_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgeru_test__Float128_SOURCES=../common/Cgeru.test.cpp $(test_srcdepends)
Cgeru_test__Float128_CXXFLAGS=$(test_cxxflags)
Cgeru_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgeru_test__Float128_opt_SOURCES=../common/Cgeru.test.cpp $(test_srcdepends)
Cgeru_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Cgeru_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chbmv_test__Float128_SOURCES=../common/Chbmv.test.cpp $(test_srcdepends)
Chbmv_test__Float128_CXXFLAGS=$(test_cxxflags)
Chbmv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chbmv_test__Float128_opt_SOURCES=../common/Chbmv.test.cpp $(test_srcdepends)
Chbmv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Chbmv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemm_test__Float128_SOURCES=../common/Chemm.test.cpp $(test_srcdepends)
Chemm_test__Float128_CXXFLAGS=$(test_cxxflags)
Chemm_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chemm_test__Float128_opt_SOURCES=../common/Chemm.test.cpp $(test_srcdepends)
Chemm_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Chemm_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemv_test__Float128_SOURCES=../common/Chemv.test.cpp $(test_srcdepends)
Chemv_test__Float128_CXXFLAGS=$(test_cxxflags)
Chemv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chemv_test__Float128_opt_SOURCES=../common/Chemv.test.cpp $(test_srcdepends)
Chemv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Chemv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher_test__Float128_SOURCES=../common/Cher.test.cpp $(test_srcdepends)
Cher_test__Float128_CXXFLAGS=$(test_cxxflags)
Cher_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher_test__Float128_opt_SOURCES=../common/Cher.test.cpp $(test_srcdepends)
Cher_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Cher_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2_test__Float128_SOURCES=../common/Cher2.test.cpp $(test_srcdepends)
Cher2_test__Float128_CXXFLAGS=$(test_cxxflags)
Cher2_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher2_test__Float128_opt_SOURCES=../common/Cher2.test.cpp $(test_srcdepends)
Cher2_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Cher2_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2k_test__Float128_SOURCES=../common/Cher2k.test.cpp $(test_srcdepends)
Cher2k_test__Float128_CXXFLAGS=$(test_cxxflags)
Cher2k_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher2k_test__Float128_opt_SOURCES=../common/Cher2k.test.cpp $(test_srcdepends)
Cher2k_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Cher2k_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cherk_test__Float128_SOURCES=../common/Cherk.test.cpp $(test_srcdepends)
Cherk_test__Float128_CXXFLAGS=$(test_cxxflags)
Cherk_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cherk_test__Float128_opt_SOURCES=../common/Cherk.test.cpp $(test_srcdepends)
Cherk_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Cherk_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpmv_test__Float128_SOURCES=../common/Chpmv.test.cpp $(test_srcdepends)
Chpmv_test__Float128_CXXFLAGS=$(test_cxxflags)
Chpmv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpmv_test__Float128_opt_SOURCES=../common/Chpmv.test.cpp $(test_srcdepends)
Chpmv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Chpmv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr_test__Float128_SOURCES=../common/Chpr.test.cpp $(test_srcdepends)
Chpr_test__Float128_CXXFLAGS=$(test_cxxflags)
Chpr_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpr_test__Float128_opt_SOURCES=../common/Chpr.test.cpp $(test_srcdepends)
Chpr_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Chpr_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr2_test__Float128_SOURCES=../common/Chpr2.test.cpp $(test_srcdepends)
Chpr2_test__Float128_CXXFLAGS=$(test_cxxflags)
Chpr2_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpr2_test__Float128_opt_SOURCES=../common/Chpr2.test.cpp $(test_srcdepends)
Chpr2_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Chpr2_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Crotg_test__Float128_SOURCES=../common/Crotg.test.cpp $(test_srcdepends)
Crotg_test__Float128_CXXFLAGS=$(test_cxxflags)
Crotg_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Crotg_test__Float128_opt_SOURCES=../common/Crotg.test.cpp $(test_srcdepends)
Crotg_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Crotg_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cscal_test__Float128_SOURCES=../common/Cscal.test.cpp $(test_srcdepends)
Cscal_test__Float128_CXXFLAGS=$(test_cxxflags)
Cscal_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cscal_test__Float128_opt_SOURCES=../common/Cscal.test.cpp $(test_srcdepends)
Cscal_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Cscal_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cswap_test__Float128_SOURCES=../common/Cswap.test.cpp $(test_srcdepends)
Cswap_test__Float128_CXXFLAGS=$(test_cxxflags)
Cswap_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cswap_test__Float128_opt_SOURCES=../common/Cswap.test.cpp $(test_srcdepends)
Cswap_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Cswap_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csymm_test__Float128_SOURCES=../common/Csymm.test.cpp $(test_srcdepends)
Csymm_test__Float128_CXXFLAGS=$(test_cxxflags)
Csymm_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csymm_test__Float128_opt_SOURCES=../common/Csymm.test.cpp $(test_srcdepends)
Csymm_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Csymm_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyr2k_test__Float128_SOURCES=../common/Csyr2k.test.cpp $(test_srcdepends)
Csyr2k_test__Float128_CXXFLAGS=$(test_cxxflags)
Csyr2k_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csyr2k_test__Float128_opt_SOURCES=../common/Csyr2k.test.cpp $(test_srcdepends)
Csyr2k_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Csyr2k_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyrk_test__Float128_SOURCES=../common/Csyrk.test.cpp $(test_srcdepends)
Csyrk_test__Float128_CXXFLAGS=$(test_cxxflags)
Csyrk_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csyrk_test__Float128_opt_SOURCES=../common/Csyrk.test.cpp $(test_srcdepends)
Csyrk_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Csyrk_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbmv_test__Float128_SOURCES=../common/Ctbmv.test.cpp $(test_srcdepends)
Ctbmv_test__Float128_CXXFLAGS=$(test_cxxflags)
Ctbmv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctbmv_test__Float128_opt_SOURCES=../common/Ctbmv.test.cpp $(test_srcdepends)
Ctbmv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Ctbmv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbsv_test__Float128_SOURCES=../common/Ctbsv.test.cpp $(test_srcdepends)
Ctbsv_test__Float128_CXXFLAGS=$(test_cxxflags)
Ctbsv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctbsv_test__Float128_opt_SOURCES=../common/Ctbsv.test.cpp $(test_srcdepends)
Ctbsv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Ctbsv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpmv_test__Float128_SOURCES=../common/Ctpmv.test.cpp $(test_srcdepends)
Ctpmv_test__Float128_CXXFLAGS=$(test_cxxflags)
Ctpmv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctpmv_test__Float128_opt_SOURCES=../common/Ctpmv.test.cpp $(test_srcdepends)
Ctpmv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Ctpmv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpsv_test__Float128_SOURCES=../common/Ctpsv.test.cpp $(test_srcdepends)
Ctpsv_test__Float128_CXXFLAGS=$(test_cxxflags)
Ctpsv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctpsv_test__Float128_opt_SOURCES=../common/Ctpsv.test.cpp $(test_srcdepends)
Ctpsv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Ctpsv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmm_test__Float128_SOURCES=../common/Ctrmm.test.cpp $(test_srcdepends)
Ctrmm_test__Float128_CXXFLAGS=$(test_cxxflags)
Ctrmm_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrmm_test__Float128_opt_SOURCES=../common/Ctrmm.test.cpp $(test_srcdepends)
Ctrmm_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Ctrmm_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmv_test__Float128_SOURCES=../common/Ctrmv.test.cpp $(test_srcdepends)
Ctrmv_test__Float128_CXXFLAGS=$(test_cxxflags)
Ctrmv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrmv_test__Float128_opt_SOURCES=../common/Ctrmv.test.cpp $(test_srcdepends)
Ctrmv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Ctrmv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsm_test__Float128_SOURCES=../common/Ctrsm.test.cpp $(test_srcdepends)
Ctrsm_test__Float128_CXXFLAGS=$(test_cxxflags)
Ctrsm_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrsm_test__Float128_opt_SOURCES=../common/Ctrsm.test.cpp $(test_srcdepends)
Ctrsm_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Ctrsm_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsv_test__Float128_SOURCES=../common/Ctrsv.test.cpp $(test_srcdepends)
Ctrsv_test__Float128_CXXFLAGS=$(test_cxxflags)
Ctrsv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrsv_test__Float128_opt_SOURCES=../common/Ctrsv.test.cpp $(test_srcdepends)
Ctrsv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Ctrsv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCabs1_test__Float128_opt_SOURCES=../common/RCabs1.test.cpp $(test_srcdepends)
RCabs1_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
RCabs1_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCabs1_test__Float128_SOURCES=../common/RCabs1.test.cpp $(test_srcdepends)
RCabs1_test__Float128_CXXFLAGS=$(test_cxxflags)
RCabs1_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

RCasum_test__Float128_opt_SOURCES=../common/RCasum.test.cpp $(test_srcdepends)
RCasum_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
RCasum_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCasum_test__Float128_SOURCES=../common/RCasum.test.cpp $(test_srcdepends)
RCasum_test__Float128_CXXFLAGS=$(test_cxxflags)
RCasum_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

RCnrm2_test__Float128_opt_SOURCES=../common/RCnrm2.test.cpp $(test_srcdepends)
RCnrm2_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
RCnrm2_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCnrm2_test__Float128_SOURCES=../common/RCnrm2.test.cpp $(test_srcdepends)
RCnrm2_test__Float128_CXXFLAGS=$(test_cxxflags)
RCnrm2_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rasum_test__Float128_opt_SOURCES=../common/Rasum.test.cpp $(test_srcdepends)
Rasum_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rasum_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rasum_test__Float128_SOURCES=../common/Rasum.test.cpp $(test_srcdepends)
Rasum_test__Float128_CXXFLAGS=$(test_cxxflags)
Rasum_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Raxpy_test__Float128_opt_SOURCES=../common/Raxpy.test.cpp $(test_srcdepends)
Raxpy_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Raxpy_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Raxpy_test__Float128_SOURCES=../common/Raxpy.test.cpp $(test_srcdepends)
Raxpy_test__Float128_CXXFLAGS=$(test_cxxflags)
Raxpy_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rcopy_test__Float128_opt_SOURCES=../common/Rcopy.test.cpp $(test_srcdepends)
Rcopy_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rcopy_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rcopy_test__Float128_SOURCES=../common/Rcopy.test.cpp $(test_srcdepends)
Rcopy_test__Float128_CXXFLAGS=$(test_cxxflags)
Rcopy_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rdot_test__Float128_opt_SOURCES=../common/Rdot.test.cpp $(test_srcdepends)
Rdot_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rdot_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rdot_test__Float128_SOURCES=../common/Rdot.test.cpp $(test_srcdepends)
Rdot_test__Float128_CXXFLAGS=$(test_cxxflags)
Rdot_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgbmv_test__Float128_opt_SOURCES=../common/Rgbmv.test.cpp $(test_srcdepends)
Rgbmv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rgbmv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgbmv_test__Float128_SOURCES=../common/Rgbmv.test.cpp $(test_srcdepends)
Rgbmv_test__Float128_CXXFLAGS=$(test_cxxflags)
Rgbmv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgemm_test__Float128_opt_SOURCES=../common/Rgemm.test.cpp $(test_srcdepends)
Rgemm_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rgemm_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemm_test__Float128_SOURCES=../common/Rgemm.test.cpp $(test_srcdepends)
Rgemm_test__Float128_CXXFLAGS=$(test_cxxflags)
Rgemm_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgemv_test__Float128_opt_SOURCES=../common/Rgemv.test.cpp $(test_srcdepends)
Rgemv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rgemv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemv_test__Float128_SOURCES=../common/Rgemv.test.cpp $(test_srcdepends)
Rgemv_test__Float128_CXXFLAGS=$(test_cxxflags)
Rgemv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rger_test__Float128_opt_SOURCES=../common/Rger.test.cpp $(test_srcdepends)
Rger_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rger_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rger_test__Float128_SOURCES=../common/Rger.test.cpp $(test_srcdepends)
Rger_test__Float128_CXXFLAGS=$(test_cxxflags)
Rger_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rnrm2_test__Float128_opt_SOURCES=../common/Rnrm2.test.cpp $(test_srcdepends)
Rnrm2_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rnrm2_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rnrm2_test__Float128_SOURCES=../common/Rnrm2.test.cpp $(test_srcdepends)
Rnrm2_test__Float128_CXXFLAGS=$(test_cxxflags)
Rnrm2_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrot_test__Float128_opt_SOURCES=../common/Rrot.test.cpp $(test_srcdepends)
Rrot_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rrot_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrot_test__Float128_SOURCES=../common/Rrot.test.cpp $(test_srcdepends)
Rrot_test__Float128_CXXFLAGS=$(test_cxxflags)
Rrot_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrotg_test__Float128_opt_SOURCES=../common/Rrotg.test.cpp $(test_srcdepends)
Rrotg_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rrotg_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotg_test__Float128_SOURCES=../common/Rrotg.test.cpp $(test_srcdepends)
Rrotg_test__Float128_CXXFLAGS=$(test_cxxflags)
Rrotg_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrotm_test__Float128_opt_SOURCES=../common/Rrotm.test.cpp $(test_srcdepends)
Rrotm_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rrotm_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotm_test__Float128_SOURCES=../common/Rrotm.test.cpp $(test_srcdepends)
Rrotm_test__Float128_CXXFLAGS=$(test_cxxflags)
Rrotm_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsbmv_test__Float128_opt_SOURCES=../common/Rsbmv.test.cpp $(test_srcdepends)
Rsbmv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rsbmv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsbmv_test__Float128_SOURCES=../common/Rsbmv.test.cpp $(test_srcdepends)
Rsbmv_test__Float128_CXXFLAGS=$(test_cxxflags)
Rsbmv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rscal_test__Float128_opt_SOURCES=../common/Rscal.test.cpp $(test_srcdepends)
Rscal_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rscal_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rscal_test__Float128_SOURCES=../common/Rscal.test.cpp $(test_srcdepends)
Rscal_test__Float128_CXXFLAGS=$(test_cxxflags)
Rscal_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspmv_test__Float128_opt_SOURCES=../common/Rspmv.test.cpp $(test_srcdepends)
Rspmv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rspmv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspmv_test__Float128_SOURCES=../common/Rspmv.test.cpp $(test_srcdepends)
Rspmv_test__Float128_CXXFLAGS=$(test_cxxflags)
Rspmv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspr_test__Float128_opt_SOURCES=../common/Rspr.test.cpp $(test_srcdepends)
Rspr_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rspr_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr_test__Float128_SOURCES=../common/Rspr.test.cpp $(test_srcdepends)
Rspr_test__Float128_CXXFLAGS=$(test_cxxflags)
Rspr_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspr2_test__Float128_opt_SOURCES=../common/Rspr2.test.cpp $(test_srcdepends)
Rspr2_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rspr2_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr2_test__Float128_SOURCES=../common/Rspr2.test.cpp $(test_srcdepends)
Rspr2_test__Float128_CXXFLAGS=$(test_cxxflags)
Rspr2_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rswap_test__Float128_opt_SOURCES=../common/Rswap.test.cpp $(test_srcdepends)
Rswap_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rswap_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rswap_test__Float128_SOURCES=../common/Rswap.test.cpp $(test_srcdepends)
Rswap_test__Float128_CXXFLAGS=$(test_cxxflags)
Rswap_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsymm_test__Float128_opt_SOURCES=../common/Rsymm.test.cpp $(test_srcdepends)
Rsymm_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rsymm_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymm_test__Float128_SOURCES=../common/Rsymm.test.cpp $(test_srcdepends)
Rsymm_test__Float128_CXXFLAGS=$(test_cxxflags)
Rsymm_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsymv_test__Float128_opt_SOURCES=../common/Rsymv.test.cpp $(test_srcdepends)
Rsymv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rsymv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymv_test__Float128_SOURCES=../common/Rsymv.test.cpp $(test_srcdepends)
Rsymv_test__Float128_CXXFLAGS=$(test_cxxflags)
Rsymv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr_test__Float128_opt_SOURCES=../common/Rsyr.test.cpp $(test_srcdepends)
Rsyr_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rsyr_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr_test__Float128_SOURCES=../common/Rsyr.test.cpp $(test_srcdepends)
Rsyr_test__Float128_CXXFLAGS=$(test_cxxflags)
Rsyr_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr2_test__Float128_opt_SOURCES=../common/Rsyr2.test.cpp $(test_srcdepends)
Rsyr2_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rsyr2_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2_test__Float128_SOURCES=../common/Rsyr2.test.cpp $(test_srcdepends)
Rsyr2_test__Float128_CXXFLAGS=$(test_cxxflags)
Rsyr2_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr2k_test__Float128_opt_SOURCES=../common/Rsyr2k.test.cpp $(test_srcdepends)
Rsyr2k_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rsyr2k_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2k_test__Float128_SOURCES=../common/Rsyr2k.test.cpp $(test_srcdepends)
Rsyr2k_test__Float128_CXXFLAGS=$(test_cxxflags)
Rsyr2k_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyrk_test__Float128_opt_SOURCES=../common/Rsyrk.test.cpp $(test_srcdepends)
Rsyrk_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rsyrk_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyrk_test__Float128_SOURCES=../common/Rsyrk.test.cpp $(test_srcdepends)
Rsyrk_test__Float128_CXXFLAGS=$(test_cxxflags)
Rsyrk_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtbmv_test__Float128_opt_SOURCES=../common/Rtbmv.test.cpp $(test_srcdepends)
Rtbmv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rtbmv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbmv_test__Float128_SOURCES=../common/Rtbmv.test.cpp $(test_srcdepends)
Rtbmv_test__Float128_CXXFLAGS=$(test_cxxflags)
Rtbmv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtbsv_test__Float128_opt_SOURCES=../common/Rtbsv.test.cpp $(test_srcdepends)
Rtbsv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rtbsv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbsv_test__Float128_SOURCES=../common/Rtbsv.test.cpp $(test_srcdepends)
Rtbsv_test__Float128_CXXFLAGS=$(test_cxxflags)
Rtbsv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtpmv_test__Float128_opt_SOURCES=../common/Rtpmv.test.cpp $(test_srcdepends)
Rtpmv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rtpmv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpmv_test__Float128_SOURCES=../common/Rtpmv.test.cpp $(test_srcdepends)
Rtpmv_test__Float128_CXXFLAGS=$(test_cxxflags)
Rtpmv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtpsv_test__Float128_opt_SOURCES=../common/Rtpsv.test.cpp $(test_srcdepends)
Rtpsv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rtpsv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpsv_test__Float128_SOURCES=../common/Rtpsv.test.cpp $(test_srcdepends)
Rtpsv_test__Float128_CXXFLAGS=$(test_cxxflags)
Rtpsv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrmm_test__Float128_opt_SOURCES=../common/Rtrmm.test.cpp $(test_srcdepends)
Rtrmm_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rtrmm_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmm_test__Float128_SOURCES=../common/Rtrmm.test.cpp $(test_srcdepends)
Rtrmm_test__Float128_CXXFLAGS=$(test_cxxflags)
Rtrmm_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrmv_test__Float128_opt_SOURCES=../common/Rtrmv.test.cpp $(test_srcdepends)
Rtrmv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rtrmv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmv_test__Float128_SOURCES=../common/Rtrmv.test.cpp $(test_srcdepends)
Rtrmv_test__Float128_CXXFLAGS=$(test_cxxflags)
Rtrmv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrsm_test__Float128_opt_SOURCES=../common/Rtrsm.test.cpp $(test_srcdepends)
Rtrsm_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rtrsm_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsm_test__Float128_SOURCES=../common/Rtrsm.test.cpp $(test_srcdepends)
Rtrsm_test__Float128_CXXFLAGS=$(test_cxxflags)
Rtrsm_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrsv_test__Float128_opt_SOURCES=../common/Rtrsv.test.cpp $(test_srcdepends)
Rtrsv_test__Float128_opt_CXXFLAGS=$(test_cxxflags)
Rtrsv_test__Float128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsv_test__Float128_SOURCES=../common/Rtrsv.test.cpp $(test_srcdepends)
Rtrsv_test__Float128_CXXFLAGS=$(test_cxxflags)
Rtrsv_test__Float128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)
