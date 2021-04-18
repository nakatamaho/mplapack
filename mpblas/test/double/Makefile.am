#CTBSV, CTPSV might have bugs or intrinsically inaccurate
#Rnrm2 should be checked more
check_PROGRAMS = $(mpblas_double_test_PROGRAMS) 
TESTS = $(check_PROGRAMS)
mpblas_double_test_PROGRAMS = \
arithmetic.test_double \
complex.test_double \
iCamax.test.double_ref iCamax.test.double   \
iRamax.test.double_ref iRamax.test.double   \
Mlsame.test.double_ref Mlsame.test.double   \
Mxerbla.test.double_ref Mxerbla.test.double \
CRrot.test.double_ref  CRrot.test.double    \
CRscal.test.double_ref CRscal.test.double   \
Caxpy.test.double_ref  Caxpy.test.double    \
Ccopy.test.double_ref  Ccopy.test.double    \
Cdotc.test.double_ref  Cdotc.test.double	\
Cdotu.test.double_ref  Cdotu.test.double	\
Cgbmv.test.double_ref  Cgbmv.test.double	\
Cgemm.test.double_ref  Cgemm.test.double	\
Cgemv.test.double_ref  Cgemv.test.double	\
Cgerc.test.double_ref  Cgerc.test.double	\
Cgeru.test.double_ref  Cgeru.test.double    \
Chbmv.test.double_ref  Chbmv.test.double	\
Chemm.test.double_ref  Chemm.test.double	\
Chemv.test.double_ref  Chemv.test.double	\
Cher.test.double_ref   Cher.test.double 	\
Cher2.test.double_ref  Cher2.test.double	\
Cher2k.test.double_ref Cher2k.test.double	\
Cherk.test.double_ref  Cherk.test.double	\
Chpmv.test.double_ref  Chpmv.test.double	\
Chpr.test.double_ref   Chpr.test.double	\
Chpr2.test.double_ref  Chpr2.test.double	\
Crotg.test.double_ref  Crotg.test.double	\
Cscal.test.double_ref  Cscal.test.double	\
Cswap.test.double_ref  Cswap.test.double	\
Csymm.test.double_ref  Csymm.test.double    \
Csyr2k.test.double_ref Csyr2k.test.double	\
Csyrk.test.double_ref  Csyrk.test.double	\
Ctbmv.test.double_ref  Ctbmv.test.double	\
Ctbsv.test.double_ref  Ctbsv.test.double	\
Ctpmv.test.double_ref  Ctpmv.test.double	\
Ctpsv.test.double_ref  Ctpsv.test.double	\
Ctrmm.test.double_ref  Ctrmm.test.double	\
Ctrmv.test.double_ref  Ctrmv.test.double	\
Ctrsm.test.double_ref  Ctrsm.test.double	\
Ctrsv.test.double_ref  Ctrsv.test.double	\
RCabs1.test.double_ref RCabs1.test.double	\
RCasum.test.double_ref RCasum.test.double	\
RCnrm2.test.double_ref RCnrm2.test.double	\
Rasum.test.double_ref  Rasum.test.double    \
Raxpy.test.double_ref  Raxpy.test.double	\
Rcopy.test.double_ref  Rcopy.test.double	\
Rdot.test.double_ref   Rdot.test.double	\
Rgbmv.test.double_ref  Rgbmv.test.double	\
Rgemm.test.double_ref  Rgemm.test.double	\
Rgemv.test.double_ref  Rgemv.test.double	\
Rger.test.double_ref   Rger.test.double	\
Rnrm2.test.double_ref  Rnrm2.test.double	\
Rrot.test.double_ref   Rrot.test.double	\
Rrotg.test.double_ref  Rrotg.test.double	\
Rrotm.test.double_ref  Rrotm.test.double	\
Rsbmv.test.double_ref  Rsbmv.test.double	\
Rscal.test.double_ref  Rscal.test.double	\
Rspmv.test.double_ref  Rspmv.test.double    \
Rspr.test.double_ref   Rspr.test.double	\
Rspr2.test.double_ref  Rspr2.test.double	\
Rswap.test.double_ref  Rswap.test.double	\
Rsymm.test.double_ref  Rsymm.test.double	\
Rsymv.test.double_ref  Rsymv.test.double	\
Rsyr.test.double_ref   Rsyr.test.double	\
Rsyr2.test.double_ref  Rsyr2.test.double	\
Rsyr2k.test.double_ref Rsyr2k.test.double	\
Rsyrk.test.double_ref  Rsyrk.test.double	\
Rtbmv.test.double_ref  Rtbmv.test.double	\
Rtbsv.test.double_ref  Rtbsv.test.double	\
Rtpmv.test.double_ref  Rtpmv.test.double	\
Rtpsv.test.double_ref  Rtpsv.test.double	\
Rtrmm.test.double_ref  Rtrmm.test.double    \
Rtrmv.test.double_ref  Rtrmv.test.double	\
Rtrsm.test.double_ref  Rtrsm.test.double	\
Rtrsv.test.double_ref  Rtrsv.test.double

mpblas_double_testdir = $(prefix)/share/mplapack/test/mpblas/double
ref_libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_double_ref -lmpblas_mpfr_ref
opt_libdepends = -L$(top_builddir)/mpblas/optimized/double -lmpblas_double -L$(top_builddir)/mpblas/reference -lmpblas_mpfr_ref
mplibs=-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp
refblas=-L$(top_builddir)/external/i/LAPACK/lib -lblas $(FCLIBS)

test_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_DOUBLE___

if !IS_MACOSX
test_srcdepends = ../common/xerbla.cpp ../common/mplapack.test.cpp ../common/Mxerbla.override.cpp
endif

if IS_MACOSX
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_DOUBLE___
libmxerbla_override_la_SOURCES =  ../common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS= -Wl,-flat_namespace,-undefined,dynamic_lookup
test_srcdepends = ../common/xerbla.cpp ../common/mplapack.test.cpp
mplibs+=-Wl,-flat_namespace,-undefined,dynamic_lookup
TESTS_ENVIRONMENT=DYLD_INSERT_LIBRARIES=.libs/libmxerbla_override.0.dylib
endif

##
arithmetic_test_double_SOURCES=../common/arithmetic.test.cpp $(test_srcdepends)
arithmetic_test_double_CXXFLAGS=$(test_cxxflags)
arithmetic_test_double_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

complex_test_double_SOURCES=../common/complex.test.cpp $(test_srcdepends)
complex_test_double_CXXFLAGS=$(test_cxxflags)
complex_test_double_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_test_double_ref_SOURCES=../common/iCamax.test.cpp $(test_srcdepends)
iCamax_test_double_ref_CXXFLAGS=$(test_cxxflags)
iCamax_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iCamax_test_double_SOURCES=../common/iCamax.test.cpp $(test_srcdepends)
iCamax_test_double_CXXFLAGS=$(test_cxxflags)
iCamax_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iRamax_test_double_ref_SOURCES=../common/iRamax.test.cpp $(test_srcdepends)
iRamax_test_double_ref_CXXFLAGS=$(test_cxxflags)
iRamax_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

iRamax_test_double_SOURCES=../common/iRamax.test.cpp $(test_srcdepends)
iRamax_test_double_CXXFLAGS=$(test_cxxflags)
iRamax_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mlsame_test_double_ref_SOURCES=../common/Mlsame.test.cpp $(test_srcdepends)
Mlsame_test_double_ref_CXXFLAGS=$(test_cxxflags)
Mlsame_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mlsame_test_double_SOURCES=../common/Mlsame.test.cpp $(test_srcdepends)
Mlsame_test_double_CXXFLAGS=$(test_cxxflags)
Mlsame_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mxerbla_test_double_ref_SOURCES=../common/Mxerbla.test.cpp $(test_srcdepends)
Mxerbla_test_double_ref_CXXFLAGS=$(test_cxxflags)
Mxerbla_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Mxerbla_test_double_SOURCES=../common/Mxerbla.test.cpp $(test_srcdepends)
Mxerbla_test_double_CXXFLAGS=$(test_cxxflags)
Mxerbla_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRrot_test_double_ref_SOURCES=../common/CRrot.test.cpp $(test_srcdepends)
CRrot_test_double_ref_CXXFLAGS=$(test_cxxflags)
CRrot_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRrot_test_double_SOURCES=../common/CRrot.test.cpp $(test_srcdepends)
CRrot_test_double_CXXFLAGS=$(test_cxxflags)
CRrot_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRscal_test_double_ref_SOURCES=../common/CRscal.test.cpp $(test_srcdepends)
CRscal_test_double_ref_CXXFLAGS=$(test_cxxflags)
CRscal_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

CRscal_test_double_SOURCES=../common/CRscal.test.cpp $(test_srcdepends)
CRscal_test_double_CXXFLAGS=$(test_cxxflags)
CRscal_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_test_double_SOURCES=../common/Caxpy.test.cpp $(test_srcdepends)
Caxpy_test_double_CXXFLAGS=$(test_cxxflags)
Caxpy_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Caxpy_test_double_ref_SOURCES=../common/Caxpy.test.cpp $(test_srcdepends)
Caxpy_test_double_ref_CXXFLAGS=$(test_cxxflags)
Caxpy_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ccopy_test_double_SOURCES=../common/Ccopy.test.cpp $(test_srcdepends)
Ccopy_test_double_CXXFLAGS=$(test_cxxflags)
Ccopy_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ccopy_test_double_ref_SOURCES=../common/Ccopy.test.cpp $(test_srcdepends)
Ccopy_test_double_ref_CXXFLAGS=$(test_cxxflags)
Ccopy_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotc_test_double_SOURCES=../common/Cdotc.test.cpp $(test_srcdepends)
Cdotc_test_double_CXXFLAGS=$(test_cxxflags)
Cdotc_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotc_test_double_ref_SOURCES=../common/Cdotc.test.cpp $(test_srcdepends)
Cdotc_test_double_ref_CXXFLAGS=$(test_cxxflags)
Cdotc_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cdotu_test_double_SOURCES=../common/Cdotu.test.cpp $(test_srcdepends)
Cdotu_test_double_CXXFLAGS=$(test_cxxflags)
Cdotu_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotu_test_double_ref_SOURCES=../common/Cdotu.test.cpp $(test_srcdepends)
Cdotu_test_double_ref_CXXFLAGS=$(test_cxxflags)
Cdotu_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgbmv_test_double_SOURCES=../common/Cgbmv.test.cpp $(test_srcdepends)
Cgbmv_test_double_CXXFLAGS=$(test_cxxflags)
Cgbmv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgbmv_test_double_ref_SOURCES=../common/Cgbmv.test.cpp $(test_srcdepends)
Cgbmv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Cgbmv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemm_test_double_SOURCES=../common/Cgemm.test.cpp $(test_srcdepends)
Cgemm_test_double_CXXFLAGS=$(test_cxxflags)
Cgemm_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemm_test_double_ref_SOURCES=../common/Cgemm.test.cpp $(test_srcdepends)
Cgemm_test_double_ref_CXXFLAGS=$(test_cxxflags)
Cgemm_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgemv_test_double_SOURCES=../common/Cgemv.test.cpp $(test_srcdepends)
Cgemv_test_double_CXXFLAGS=$(test_cxxflags)
Cgemv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemv_test_double_ref_SOURCES=../common/Cgemv.test.cpp $(test_srcdepends)
Cgemv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Cgemv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgerc_test_double_SOURCES=../common/Cgerc.test.cpp $(test_srcdepends)
Cgerc_test_double_CXXFLAGS=$(test_cxxflags)
Cgerc_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgerc_test_double_ref_SOURCES=../common/Cgerc.test.cpp $(test_srcdepends)
Cgerc_test_double_ref_CXXFLAGS=$(test_cxxflags)
Cgerc_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cgeru_test_double_SOURCES=../common/Cgeru.test.cpp $(test_srcdepends)
Cgeru_test_double_CXXFLAGS=$(test_cxxflags)
Cgeru_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgeru_test_double_ref_SOURCES=../common/Cgeru.test.cpp $(test_srcdepends)
Cgeru_test_double_ref_CXXFLAGS=$(test_cxxflags)
Cgeru_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chbmv_test_double_SOURCES=../common/Chbmv.test.cpp $(test_srcdepends)
Chbmv_test_double_CXXFLAGS=$(test_cxxflags)
Chbmv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chbmv_test_double_ref_SOURCES=../common/Chbmv.test.cpp $(test_srcdepends)
Chbmv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Chbmv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemm_test_double_SOURCES=../common/Chemm.test.cpp $(test_srcdepends)
Chemm_test_double_CXXFLAGS=$(test_cxxflags)
Chemm_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemm_test_double_ref_SOURCES=../common/Chemm.test.cpp $(test_srcdepends)
Chemm_test_double_ref_CXXFLAGS=$(test_cxxflags)
Chemm_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chemv_test_double_SOURCES=../common/Chemv.test.cpp $(test_srcdepends)
Chemv_test_double_CXXFLAGS=$(test_cxxflags)
Chemv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemv_test_double_ref_SOURCES=../common/Chemv.test.cpp $(test_srcdepends)
Chemv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Chemv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher_test_double_SOURCES=../common/Cher.test.cpp $(test_srcdepends)
Cher_test_double_CXXFLAGS=$(test_cxxflags)
Cher_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher_test_double_ref_SOURCES=../common/Cher.test.cpp $(test_srcdepends)
Cher_test_double_ref_CXXFLAGS=$(test_cxxflags)
Cher_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2_test_double_SOURCES=../common/Cher2.test.cpp $(test_srcdepends)
Cher2_test_double_CXXFLAGS=$(test_cxxflags)
Cher2_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2_test_double_ref_SOURCES=../common/Cher2.test.cpp $(test_srcdepends)
Cher2_test_double_ref_CXXFLAGS=$(test_cxxflags)
Cher2_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cher2k_test_double_SOURCES=../common/Cher2k.test.cpp $(test_srcdepends)
Cher2k_test_double_CXXFLAGS=$(test_cxxflags)
Cher2k_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2k_test_double_ref_SOURCES=../common/Cher2k.test.cpp $(test_srcdepends)
Cher2k_test_double_ref_CXXFLAGS=$(test_cxxflags)
Cher2k_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cherk_test_double_SOURCES=../common/Cherk.test.cpp $(test_srcdepends)
Cherk_test_double_CXXFLAGS=$(test_cxxflags)
Cherk_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cherk_test_double_ref_SOURCES=../common/Cherk.test.cpp $(test_srcdepends)
Cherk_test_double_ref_CXXFLAGS=$(test_cxxflags)
Cherk_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpmv_test_double_SOURCES=../common/Chpmv.test.cpp $(test_srcdepends)
Chpmv_test_double_CXXFLAGS=$(test_cxxflags)
Chpmv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpmv_test_double_ref_SOURCES=../common/Chpmv.test.cpp $(test_srcdepends)
Chpmv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Chpmv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr_test_double_SOURCES=../common/Chpr.test.cpp $(test_srcdepends)
Chpr_test_double_CXXFLAGS=$(test_cxxflags)
Chpr_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr_test_double_ref_SOURCES=../common/Chpr.test.cpp $(test_srcdepends)
Chpr_test_double_ref_CXXFLAGS=$(test_cxxflags)
Chpr_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Chpr2_test_double_SOURCES=../common/Chpr2.test.cpp $(test_srcdepends)
Chpr2_test_double_CXXFLAGS=$(test_cxxflags)
Chpr2_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr2_test_double_ref_SOURCES=../common/Chpr2.test.cpp $(test_srcdepends)
Chpr2_test_double_ref_CXXFLAGS=$(test_cxxflags)
Chpr2_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Crotg_test_double_SOURCES=../common/Crotg.test.cpp $(test_srcdepends)
Crotg_test_double_CXXFLAGS=$(test_cxxflags)
Crotg_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Crotg_test_double_ref_SOURCES=../common/Crotg.test.cpp $(test_srcdepends)
Crotg_test_double_ref_CXXFLAGS=$(test_cxxflags)
Crotg_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cscal_test_double_SOURCES=../common/Cscal.test.cpp $(test_srcdepends)
Cscal_test_double_CXXFLAGS=$(test_cxxflags)
Cscal_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cscal_test_double_ref_SOURCES=../common/Cscal.test.cpp $(test_srcdepends)
Cscal_test_double_ref_CXXFLAGS=$(test_cxxflags)
Cscal_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Cswap_test_double_SOURCES=../common/Cswap.test.cpp $(test_srcdepends)
Cswap_test_double_CXXFLAGS=$(test_cxxflags)
Cswap_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cswap_test_double_ref_SOURCES=../common/Cswap.test.cpp $(test_srcdepends)
Cswap_test_double_ref_CXXFLAGS=$(test_cxxflags)
Cswap_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csymm_test_double_SOURCES=../common/Csymm.test.cpp $(test_srcdepends)
Csymm_test_double_CXXFLAGS=$(test_cxxflags)
Csymm_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csymm_test_double_ref_SOURCES=../common/Csymm.test.cpp $(test_srcdepends)
Csymm_test_double_ref_CXXFLAGS=$(test_cxxflags)
Csymm_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyr2k_test_double_SOURCES=../common/Csyr2k.test.cpp $(test_srcdepends)
Csyr2k_test_double_CXXFLAGS=$(test_cxxflags)
Csyr2k_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyr2k_test_double_ref_SOURCES=../common/Csyr2k.test.cpp $(test_srcdepends)
Csyr2k_test_double_ref_CXXFLAGS=$(test_cxxflags)
Csyr2k_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Csyrk_test_double_SOURCES=../common/Csyrk.test.cpp $(test_srcdepends)
Csyrk_test_double_CXXFLAGS=$(test_cxxflags)
Csyrk_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyrk_test_double_ref_SOURCES=../common/Csyrk.test.cpp $(test_srcdepends)
Csyrk_test_double_ref_CXXFLAGS=$(test_cxxflags)
Csyrk_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbmv_test_double_SOURCES=../common/Ctbmv.test.cpp $(test_srcdepends)
Ctbmv_test_double_CXXFLAGS=$(test_cxxflags)
Ctbmv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbmv_test_double_ref_SOURCES=../common/Ctbmv.test.cpp $(test_srcdepends)
Ctbmv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Ctbmv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctbsv_test_double_SOURCES=../common/Ctbsv.test.cpp $(test_srcdepends)
Ctbsv_test_double_CXXFLAGS=$(test_cxxflags)
Ctbsv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbsv_test_double_ref_SOURCES=../common/Ctbsv.test.cpp $(test_srcdepends)
Ctbsv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Ctbsv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpmv_test_double_SOURCES=../common/Ctpmv.test.cpp $(test_srcdepends)
Ctpmv_test_double_CXXFLAGS=$(test_cxxflags)
Ctpmv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpmv_test_double_ref_SOURCES=../common/Ctpmv.test.cpp $(test_srcdepends)
Ctpmv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Ctpmv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctpsv_test_double_SOURCES=../common/Ctpsv.test.cpp $(test_srcdepends)
Ctpsv_test_double_CXXFLAGS=$(test_cxxflags)
Ctpsv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpsv_test_double_ref_SOURCES=../common/Ctpsv.test.cpp $(test_srcdepends)
Ctpsv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Ctpsv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmm_test_double_SOURCES=../common/Ctrmm.test.cpp $(test_srcdepends)
Ctrmm_test_double_CXXFLAGS=$(test_cxxflags)
Ctrmm_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmm_test_double_ref_SOURCES=../common/Ctrmm.test.cpp $(test_srcdepends)
Ctrmm_test_double_ref_CXXFLAGS=$(test_cxxflags)
Ctrmm_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrmv_test_double_SOURCES=../common/Ctrmv.test.cpp $(test_srcdepends)
Ctrmv_test_double_CXXFLAGS=$(test_cxxflags)
Ctrmv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmv_test_double_ref_SOURCES=../common/Ctrmv.test.cpp $(test_srcdepends)
Ctrmv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Ctrmv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsm_test_double_SOURCES=../common/Ctrsm.test.cpp $(test_srcdepends)
Ctrsm_test_double_CXXFLAGS=$(test_cxxflags)
Ctrsm_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsm_test_double_ref_SOURCES=../common/Ctrsm.test.cpp $(test_srcdepends)
Ctrsm_test_double_ref_CXXFLAGS=$(test_cxxflags)
Ctrsm_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Ctrsv_test_double_SOURCES=../common/Ctrsv.test.cpp $(test_srcdepends)
Ctrsv_test_double_CXXFLAGS=$(test_cxxflags)
Ctrsv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsv_test_double_ref_SOURCES=../common/Ctrsv.test.cpp $(test_srcdepends)
Ctrsv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Ctrsv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_test_double_ref_SOURCES=../common/RCabs1.test.cpp $(test_srcdepends)
RCabs1_test_double_ref_CXXFLAGS=$(test_cxxflags)
RCabs1_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCabs1_test_double_SOURCES=../common/RCabs1.test.cpp $(test_srcdepends)
RCabs1_test_double_CXXFLAGS=$(test_cxxflags)
RCabs1_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCasum_test_double_ref_SOURCES=../common/RCasum.test.cpp $(test_srcdepends)
RCasum_test_double_ref_CXXFLAGS=$(test_cxxflags)
RCasum_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCasum_test_double_SOURCES=../common/RCasum.test.cpp $(test_srcdepends)
RCasum_test_double_CXXFLAGS=$(test_cxxflags)
RCasum_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCnrm2_test_double_ref_SOURCES=../common/RCnrm2.test.cpp $(test_srcdepends)
RCnrm2_test_double_ref_CXXFLAGS=$(test_cxxflags)
RCnrm2_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

RCnrm2_test_double_SOURCES=../common/RCnrm2.test.cpp $(test_srcdepends)
RCnrm2_test_double_CXXFLAGS=$(test_cxxflags)
RCnrm2_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rasum_test_double_ref_SOURCES=../common/Rasum.test.cpp $(test_srcdepends)
Rasum_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rasum_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rasum_test_double_SOURCES=../common/Rasum.test.cpp $(test_srcdepends)
Rasum_test_double_CXXFLAGS=$(test_cxxflags)
Rasum_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Raxpy_test_double_ref_SOURCES=../common/Raxpy.test.cpp $(test_srcdepends)
Raxpy_test_double_ref_CXXFLAGS=$(test_cxxflags)
Raxpy_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Raxpy_test_double_SOURCES=../common/Raxpy.test.cpp $(test_srcdepends)
Raxpy_test_double_CXXFLAGS=$(test_cxxflags)
Raxpy_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rcopy_test_double_ref_SOURCES=../common/Rcopy.test.cpp $(test_srcdepends)
Rcopy_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rcopy_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rcopy_test_double_SOURCES=../common/Rcopy.test.cpp $(test_srcdepends)
Rcopy_test_double_CXXFLAGS=$(test_cxxflags)
Rcopy_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rdot_test_double_ref_SOURCES=../common/Rdot.test.cpp $(test_srcdepends)
Rdot_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rdot_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rdot_test_double_SOURCES=../common/Rdot.test.cpp $(test_srcdepends)
Rdot_test_double_CXXFLAGS=$(test_cxxflags)
Rdot_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgbmv_test_double_ref_SOURCES=../common/Rgbmv.test.cpp $(test_srcdepends)
Rgbmv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rgbmv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgbmv_test_double_SOURCES=../common/Rgbmv.test.cpp $(test_srcdepends)
Rgbmv_test_double_CXXFLAGS=$(test_cxxflags)
Rgbmv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemm_test_double_ref_SOURCES=../common/Rgemm.test.cpp $(test_srcdepends)
Rgemm_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rgemm_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemm_test_double_SOURCES=../common/Rgemm.test.cpp $(test_srcdepends)
Rgemm_test_double_CXXFLAGS=$(test_cxxflags)
Rgemm_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemv_test_double_ref_SOURCES=../common/Rgemv.test.cpp $(test_srcdepends)
Rgemv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rgemv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rgemv_test_double_SOURCES=../common/Rgemv.test.cpp $(test_srcdepends)
Rgemv_test_double_CXXFLAGS=$(test_cxxflags)
Rgemv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rger_test_double_ref_SOURCES=../common/Rger.test.cpp $(test_srcdepends)
Rger_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rger_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rger_test_double_SOURCES=../common/Rger.test.cpp $(test_srcdepends)
Rger_test_double_CXXFLAGS=$(test_cxxflags)
Rger_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rnrm2_test_double_ref_SOURCES=../common/Rnrm2.test.cpp $(test_srcdepends)
Rnrm2_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rnrm2_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rnrm2_test_double_SOURCES=../common/Rnrm2.test.cpp $(test_srcdepends)
Rnrm2_test_double_CXXFLAGS=$(test_cxxflags)
Rnrm2_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrot_test_double_ref_SOURCES=../common/Rrot.test.cpp $(test_srcdepends)
Rrot_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rrot_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrot_test_double_SOURCES=../common/Rrot.test.cpp $(test_srcdepends)
Rrot_test_double_CXXFLAGS=$(test_cxxflags)
Rrot_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotg_test_double_ref_SOURCES=../common/Rrotg.test.cpp $(test_srcdepends)
Rrotg_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rrotg_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotg_test_double_SOURCES=../common/Rrotg.test.cpp $(test_srcdepends)
Rrotg_test_double_CXXFLAGS=$(test_cxxflags)
Rrotg_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotm_test_double_ref_SOURCES=../common/Rrotm.test.cpp $(test_srcdepends)
Rrotm_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rrotm_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rrotm_test_double_SOURCES=../common/Rrotm.test.cpp $(test_srcdepends)
Rrotm_test_double_CXXFLAGS=$(test_cxxflags)
Rrotm_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsbmv_test_double_ref_SOURCES=../common/Rsbmv.test.cpp $(test_srcdepends)
Rsbmv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rsbmv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsbmv_test_double_SOURCES=../common/Rsbmv.test.cpp $(test_srcdepends)
Rsbmv_test_double_CXXFLAGS=$(test_cxxflags)
Rsbmv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rscal_test_double_ref_SOURCES=../common/Rscal.test.cpp $(test_srcdepends)
Rscal_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rscal_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rscal_test_double_SOURCES=../common/Rscal.test.cpp $(test_srcdepends)
Rscal_test_double_CXXFLAGS=$(test_cxxflags)
Rscal_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspmv_test_double_ref_SOURCES=../common/Rspmv.test.cpp $(test_srcdepends)
Rspmv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rspmv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspmv_test_double_SOURCES=../common/Rspmv.test.cpp $(test_srcdepends)
Rspmv_test_double_CXXFLAGS=$(test_cxxflags)
Rspmv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr_test_double_ref_SOURCES=../common/Rspr.test.cpp $(test_srcdepends)
Rspr_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rspr_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr_test_double_SOURCES=../common/Rspr.test.cpp $(test_srcdepends)
Rspr_test_double_CXXFLAGS=$(test_cxxflags)
Rspr_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr2_test_double_ref_SOURCES=../common/Rspr2.test.cpp $(test_srcdepends)
Rspr2_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rspr2_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rspr2_test_double_SOURCES=../common/Rspr2.test.cpp $(test_srcdepends)
Rspr2_test_double_CXXFLAGS=$(test_cxxflags)
Rspr2_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rswap_test_double_ref_SOURCES=../common/Rswap.test.cpp $(test_srcdepends)
Rswap_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rswap_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rswap_test_double_SOURCES=../common/Rswap.test.cpp $(test_srcdepends)
Rswap_test_double_CXXFLAGS=$(test_cxxflags)
Rswap_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymm_test_double_ref_SOURCES=../common/Rsymm.test.cpp $(test_srcdepends)
Rsymm_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rsymm_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymm_test_double_SOURCES=../common/Rsymm.test.cpp $(test_srcdepends)
Rsymm_test_double_CXXFLAGS=$(test_cxxflags)
Rsymm_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymv_test_double_ref_SOURCES=../common/Rsymv.test.cpp $(test_srcdepends)
Rsymv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rsymv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsymv_test_double_SOURCES=../common/Rsymv.test.cpp $(test_srcdepends)
Rsymv_test_double_CXXFLAGS=$(test_cxxflags)
Rsymv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr_test_double_ref_SOURCES=../common/Rsyr.test.cpp $(test_srcdepends)
Rsyr_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rsyr_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr_test_double_SOURCES=../common/Rsyr.test.cpp $(test_srcdepends)
Rsyr_test_double_CXXFLAGS=$(test_cxxflags)
Rsyr_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2_test_double_ref_SOURCES=../common/Rsyr2.test.cpp $(test_srcdepends)
Rsyr2_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rsyr2_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2_test_double_SOURCES=../common/Rsyr2.test.cpp $(test_srcdepends)
Rsyr2_test_double_CXXFLAGS=$(test_cxxflags)
Rsyr2_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2k_test_double_ref_SOURCES=../common/Rsyr2k.test.cpp $(test_srcdepends)
Rsyr2k_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rsyr2k_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyr2k_test_double_SOURCES=../common/Rsyr2k.test.cpp $(test_srcdepends)
Rsyr2k_test_double_CXXFLAGS=$(test_cxxflags)
Rsyr2k_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyrk_test_double_ref_SOURCES=../common/Rsyrk.test.cpp $(test_srcdepends)
Rsyrk_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rsyrk_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rsyrk_test_double_SOURCES=../common/Rsyrk.test.cpp $(test_srcdepends)
Rsyrk_test_double_CXXFLAGS=$(test_cxxflags)
Rsyrk_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbmv_test_double_ref_SOURCES=../common/Rtbmv.test.cpp $(test_srcdepends)
Rtbmv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rtbmv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbmv_test_double_SOURCES=../common/Rtbmv.test.cpp $(test_srcdepends)
Rtbmv_test_double_CXXFLAGS=$(test_cxxflags)
Rtbmv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbsv_test_double_ref_SOURCES=../common/Rtbsv.test.cpp $(test_srcdepends)
Rtbsv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rtbsv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtbsv_test_double_SOURCES=../common/Rtbsv.test.cpp $(test_srcdepends)
Rtbsv_test_double_CXXFLAGS=$(test_cxxflags)
Rtbsv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpmv_test_double_ref_SOURCES=../common/Rtpmv.test.cpp $(test_srcdepends)
Rtpmv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rtpmv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpmv_test_double_SOURCES=../common/Rtpmv.test.cpp $(test_srcdepends)
Rtpmv_test_double_CXXFLAGS=$(test_cxxflags)
Rtpmv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpsv_test_double_ref_SOURCES=../common/Rtpsv.test.cpp $(test_srcdepends)
Rtpsv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rtpsv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtpsv_test_double_SOURCES=../common/Rtpsv.test.cpp $(test_srcdepends)
Rtpsv_test_double_CXXFLAGS=$(test_cxxflags)
Rtpsv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmm_test_double_ref_SOURCES=../common/Rtrmm.test.cpp $(test_srcdepends)
Rtrmm_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rtrmm_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmm_test_double_SOURCES=../common/Rtrmm.test.cpp $(test_srcdepends)
Rtrmm_test_double_CXXFLAGS=$(test_cxxflags)
Rtrmm_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmv_test_double_ref_SOURCES=../common/Rtrmv.test.cpp $(test_srcdepends)
Rtrmv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rtrmv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrmv_test_double_SOURCES=../common/Rtrmv.test.cpp $(test_srcdepends)
Rtrmv_test_double_CXXFLAGS=$(test_cxxflags)
Rtrmv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsm_test_double_ref_SOURCES=../common/Rtrsm.test.cpp $(test_srcdepends)
Rtrsm_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rtrsm_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsm_test_double_SOURCES=../common/Rtrsm.test.cpp $(test_srcdepends)
Rtrsm_test_double_CXXFLAGS=$(test_cxxflags)
Rtrsm_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsv_test_double_ref_SOURCES=../common/Rtrsv.test.cpp $(test_srcdepends)
Rtrsv_test_double_ref_CXXFLAGS=$(test_cxxflags)
Rtrsv_test_double_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(refblas)

Rtrsv_test_double_SOURCES=../common/Rtrsv.test.cpp $(test_srcdepends)
Rtrsv_test_double_CXXFLAGS=$(test_cxxflags)
Rtrsv_test_double_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

