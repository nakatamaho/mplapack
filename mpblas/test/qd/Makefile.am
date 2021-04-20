#CTBSV, CTPSV might have bugs or intrinsically inaccurate
#Rnrm2 should be checked more
check_PROGRAMS = $(mpblas_qd_test_PROGRAMS) 
TESTS = $(check_PROGRAMS)
mpblas_qd_test_PROGRAMS = \
arithmetic.test_qd \
complex.test_qd \
iCamax.test.qd_opt iCamax.test.qd   \
iRamax.test.qd_opt iRamax.test.qd   \
Mlsame.test.qd_opt Mlsame.test.qd   \
Mxerbla.test.qd_opt Mxerbla.test.qd \
CRrot.test.qd_opt  CRrot.test.qd    \
CRscal.test.qd_opt CRscal.test.qd   \
Caxpy.test.qd_opt  Caxpy.test.qd    \
Ccopy.test.qd_opt  Ccopy.test.qd    \
Cdotc.test.qd_opt  Cdotc.test.qd	\
Cdotu.test.qd_opt  Cdotu.test.qd	\
Cgbmv.test.qd_opt  Cgbmv.test.qd	\
Cgemm.test.qd_opt  Cgemm.test.qd	\
Cgemv.test.qd_opt  Cgemv.test.qd	\
Cgerc.test.qd_opt  Cgerc.test.qd	\
Cgeru.test.qd_opt  Cgeru.test.qd    \
Chbmv.test.qd_opt  Chbmv.test.qd	\
Chemm.test.qd_opt  Chemm.test.qd	\
Chemv.test.qd_opt  Chemv.test.qd	\
Cher.test.qd_opt   Cher.test.qd 	\
Cher2.test.qd_opt  Cher2.test.qd	\
Cher2k.test.qd_opt Cher2k.test.qd	\
Cherk.test.qd_opt  Cherk.test.qd	\
Chpmv.test.qd_opt  Chpmv.test.qd	\
Chpr.test.qd_opt   Chpr.test.qd	\
Chpr2.test.qd_opt  Chpr2.test.qd	\
Crotg.test.qd_opt  Crotg.test.qd	\
Cscal.test.qd_opt  Cscal.test.qd	\
Cswap.test.qd_opt  Cswap.test.qd	\
Csymm.test.qd_opt  Csymm.test.qd    \
Csyr2k.test.qd_opt Csyr2k.test.qd	\
Csyrk.test.qd_opt  Csyrk.test.qd	\
Ctbmv.test.qd_opt  Ctbmv.test.qd	\
Ctbsv.test.qd_opt  Ctbsv.test.qd	\
Ctpmv.test.qd_opt  Ctpmv.test.qd	\
Ctpsv.test.qd_opt  Ctpsv.test.qd	\
Ctrmm.test.qd_opt  Ctrmm.test.qd	\
Ctrmv.test.qd_opt  Ctrmv.test.qd	\
Ctrsm.test.qd_opt  Ctrsm.test.qd	\
Ctrsv.test.qd_opt  Ctrsv.test.qd	\
RCabs1.test.qd_opt RCabs1.test.qd	\
RCasum.test.qd_opt RCasum.test.qd	\
RCnrm2.test.qd_opt RCnrm2.test.qd	\
Rasum.test.qd_opt  Rasum.test.qd    \
Raxpy.test.qd_opt  Raxpy.test.qd	\
Rcopy.test.qd_opt  Rcopy.test.qd	\
Rdot.test.qd_opt   Rdot.test.qd	\
Rgbmv.test.qd_opt  Rgbmv.test.qd	\
Rgemm.test.qd_opt  Rgemm.test.qd	\
Rgemv.test.qd_opt  Rgemv.test.qd	\
Rger.test.qd_opt   Rger.test.qd	\
Rnrm2.test.qd_opt  Rnrm2.test.qd	\
Rrot.test.qd_opt   Rrot.test.qd	\
Rrotg.test.qd_opt  Rrotg.test.qd	\
Rrotm.test.qd_opt  Rrotm.test.qd	\
Rsbmv.test.qd_opt  Rsbmv.test.qd	\
Rscal.test.qd_opt  Rscal.test.qd	\
Rspmv.test.qd_opt  Rspmv.test.qd    \
Rspr.test.qd_opt   Rspr.test.qd	\
Rspr2.test.qd_opt  Rspr2.test.qd	\
Rswap.test.qd_opt  Rswap.test.qd	\
Rsymm.test.qd_opt  Rsymm.test.qd	\
Rsymv.test.qd_opt  Rsymv.test.qd	\
Rsyr.test.qd_opt   Rsyr.test.qd	\
Rsyr2.test.qd_opt  Rsyr2.test.qd	\
Rsyr2k.test.qd_opt Rsyr2k.test.qd	\
Rsyrk.test.qd_opt  Rsyrk.test.qd	\
Rtbmv.test.qd_opt  Rtbmv.test.qd	\
Rtbsv.test.qd_opt  Rtbsv.test.qd	\
Rtpmv.test.qd_opt  Rtpmv.test.qd	\
Rtpsv.test.qd_opt  Rtpsv.test.qd	\
Rtrmm.test.qd_opt  Rtrmm.test.qd    \
Rtrmv.test.qd_opt  Rtrmv.test.qd	\
Rtrsm.test.qd_opt  Rtrsm.test.qd	\
Rtrsv.test.qd_opt  Rtrsv.test.qd

mpblas_qd_testdir = $(prefix)/share/mplapack/test/mpblas/qd
libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_qd -lmpblas_mpfr
opt_libdepends = -L$(top_builddir)/mpblas/optimized/qd -lmpblas_qd_opt -L$(top_builddir)/mpblas/reference -lmpblas_mpfr
mplibs=-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpc -lmpfr -lgmpxx -lgmp -L$(QD_LIBDIR) -lqd
refblas=-L$(top_builddir)/external/i/LAPACK/lib -lblas $(FCLIBS)

test_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___

if !IS_MACOSX
test_srcdepends = ../common/xerbla.cpp ../common/mplapack.test.cpp ../common/Mxerbla.override.cpp
endif

if IS_MACOSX
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_QD___
libmxerbla_override_la_SOURCES =  ../common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS= -Wl,-flat_namespace,-undefined,dynamic_lookup
test_srcdepends = ../common/xerbla.cpp ../common/mplapack.test.cpp
mplibs+=-Wl,-flat_namespace,-undefined,dynamic_lookup
TESTS_ENVIRONMENT=DYLD_INSERT_LIBRARIES=.libs/libmxerbla_override.0.dylib
endif

##
arithmetic_test_qd_SOURCES=../common/arithmetic.test.cpp $(test_srcdepends)
arithmetic_test_qd_CXXFLAGS=$(test_cxxflags)
arithmetic_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

complex_test_qd_SOURCES=../common/complex.test.cpp $(test_srcdepends)
complex_test_qd_CXXFLAGS=$(test_cxxflags)
complex_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

iCamax_test_qd_opt_SOURCES=../common/iCamax.test.cpp $(test_srcdepends)
iCamax_test_qd_opt_CXXFLAGS=$(test_cxxflags)
iCamax_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iCamax_test_qd_SOURCES=../common/iCamax.test.cpp $(test_srcdepends)
iCamax_test_qd_CXXFLAGS=$(test_cxxflags)
iCamax_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

iRamax_test_qd_opt_SOURCES=../common/iRamax.test.cpp $(test_srcdepends)
iRamax_test_qd_opt_CXXFLAGS=$(test_cxxflags)
iRamax_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iRamax_test_qd_SOURCES=../common/iRamax.test.cpp $(test_srcdepends)
iRamax_test_qd_CXXFLAGS=$(test_cxxflags)
iRamax_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Mlsame_test_qd_opt_SOURCES=../common/Mlsame.test.cpp $(test_srcdepends)
Mlsame_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Mlsame_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mlsame_test_qd_SOURCES=../common/Mlsame.test.cpp $(test_srcdepends)
Mlsame_test_qd_CXXFLAGS=$(test_cxxflags)
Mlsame_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Mxerbla_test_qd_opt_SOURCES=../common/Mxerbla.test.cpp $(test_srcdepends)
Mxerbla_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Mxerbla_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mxerbla_test_qd_SOURCES=../common/Mxerbla.test.cpp $(test_srcdepends)
Mxerbla_test_qd_CXXFLAGS=$(test_cxxflags)
Mxerbla_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

CRrot_test_qd_opt_SOURCES=../common/CRrot.test.cpp $(test_srcdepends)
CRrot_test_qd_opt_CXXFLAGS=$(test_cxxflags)
CRrot_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRrot_test_qd_SOURCES=../common/CRrot.test.cpp $(test_srcdepends)
CRrot_test_qd_CXXFLAGS=$(test_cxxflags)
CRrot_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

CRscal_test_qd_opt_SOURCES=../common/CRscal.test.cpp $(test_srcdepends)
CRscal_test_qd_opt_CXXFLAGS=$(test_cxxflags)
CRscal_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRscal_test_qd_SOURCES=../common/CRscal.test.cpp $(test_srcdepends)
CRscal_test_qd_CXXFLAGS=$(test_cxxflags)
CRscal_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Caxpy_test_qd_SOURCES=../common/Caxpy.test.cpp $(test_srcdepends)
Caxpy_test_qd_CXXFLAGS=$(test_cxxflags)
Caxpy_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Caxpy_test_qd_opt_SOURCES=../common/Caxpy.test.cpp $(test_srcdepends)
Caxpy_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Caxpy_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ccopy_test_qd_SOURCES=../common/Ccopy.test.cpp $(test_srcdepends)
Ccopy_test_qd_CXXFLAGS=$(test_cxxflags)
Ccopy_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ccopy_test_qd_opt_SOURCES=../common/Ccopy.test.cpp $(test_srcdepends)
Ccopy_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Ccopy_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotc_test_qd_SOURCES=../common/Cdotc.test.cpp $(test_srcdepends)
Cdotc_test_qd_CXXFLAGS=$(test_cxxflags)
Cdotc_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cdotc_test_qd_opt_SOURCES=../common/Cdotc.test.cpp $(test_srcdepends)
Cdotc_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Cdotc_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotu_test_qd_SOURCES=../common/Cdotu.test.cpp $(test_srcdepends)
Cdotu_test_qd_CXXFLAGS=$(test_cxxflags)
Cdotu_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cdotu_test_qd_opt_SOURCES=../common/Cdotu.test.cpp $(test_srcdepends)
Cdotu_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Cdotu_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgbmv_test_qd_SOURCES=../common/Cgbmv.test.cpp $(test_srcdepends)
Cgbmv_test_qd_CXXFLAGS=$(test_cxxflags)
Cgbmv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgbmv_test_qd_opt_SOURCES=../common/Cgbmv.test.cpp $(test_srcdepends)
Cgbmv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Cgbmv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemm_test_qd_SOURCES=../common/Cgemm.test.cpp $(test_srcdepends)
Cgemm_test_qd_CXXFLAGS=$(test_cxxflags)
Cgemm_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgemm_test_qd_opt_SOURCES=../common/Cgemm.test.cpp $(test_srcdepends)
Cgemm_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Cgemm_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemv_test_qd_SOURCES=../common/Cgemv.test.cpp $(test_srcdepends)
Cgemv_test_qd_CXXFLAGS=$(test_cxxflags)
Cgemv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgemv_test_qd_opt_SOURCES=../common/Cgemv.test.cpp $(test_srcdepends)
Cgemv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Cgemv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgerc_test_qd_SOURCES=../common/Cgerc.test.cpp $(test_srcdepends)
Cgerc_test_qd_CXXFLAGS=$(test_cxxflags)
Cgerc_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgerc_test_qd_opt_SOURCES=../common/Cgerc.test.cpp $(test_srcdepends)
Cgerc_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Cgerc_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgeru_test_qd_SOURCES=../common/Cgeru.test.cpp $(test_srcdepends)
Cgeru_test_qd_CXXFLAGS=$(test_cxxflags)
Cgeru_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgeru_test_qd_opt_SOURCES=../common/Cgeru.test.cpp $(test_srcdepends)
Cgeru_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Cgeru_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chbmv_test_qd_SOURCES=../common/Chbmv.test.cpp $(test_srcdepends)
Chbmv_test_qd_CXXFLAGS=$(test_cxxflags)
Chbmv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chbmv_test_qd_opt_SOURCES=../common/Chbmv.test.cpp $(test_srcdepends)
Chbmv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Chbmv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemm_test_qd_SOURCES=../common/Chemm.test.cpp $(test_srcdepends)
Chemm_test_qd_CXXFLAGS=$(test_cxxflags)
Chemm_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chemm_test_qd_opt_SOURCES=../common/Chemm.test.cpp $(test_srcdepends)
Chemm_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Chemm_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemv_test_qd_SOURCES=../common/Chemv.test.cpp $(test_srcdepends)
Chemv_test_qd_CXXFLAGS=$(test_cxxflags)
Chemv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chemv_test_qd_opt_SOURCES=../common/Chemv.test.cpp $(test_srcdepends)
Chemv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Chemv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher_test_qd_SOURCES=../common/Cher.test.cpp $(test_srcdepends)
Cher_test_qd_CXXFLAGS=$(test_cxxflags)
Cher_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher_test_qd_opt_SOURCES=../common/Cher.test.cpp $(test_srcdepends)
Cher_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Cher_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2_test_qd_SOURCES=../common/Cher2.test.cpp $(test_srcdepends)
Cher2_test_qd_CXXFLAGS=$(test_cxxflags)
Cher2_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher2_test_qd_opt_SOURCES=../common/Cher2.test.cpp $(test_srcdepends)
Cher2_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Cher2_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2k_test_qd_SOURCES=../common/Cher2k.test.cpp $(test_srcdepends)
Cher2k_test_qd_CXXFLAGS=$(test_cxxflags)
Cher2k_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher2k_test_qd_opt_SOURCES=../common/Cher2k.test.cpp $(test_srcdepends)
Cher2k_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Cher2k_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cherk_test_qd_SOURCES=../common/Cherk.test.cpp $(test_srcdepends)
Cherk_test_qd_CXXFLAGS=$(test_cxxflags)
Cherk_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cherk_test_qd_opt_SOURCES=../common/Cherk.test.cpp $(test_srcdepends)
Cherk_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Cherk_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpmv_test_qd_SOURCES=../common/Chpmv.test.cpp $(test_srcdepends)
Chpmv_test_qd_CXXFLAGS=$(test_cxxflags)
Chpmv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpmv_test_qd_opt_SOURCES=../common/Chpmv.test.cpp $(test_srcdepends)
Chpmv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Chpmv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr_test_qd_SOURCES=../common/Chpr.test.cpp $(test_srcdepends)
Chpr_test_qd_CXXFLAGS=$(test_cxxflags)
Chpr_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpr_test_qd_opt_SOURCES=../common/Chpr.test.cpp $(test_srcdepends)
Chpr_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Chpr_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr2_test_qd_SOURCES=../common/Chpr2.test.cpp $(test_srcdepends)
Chpr2_test_qd_CXXFLAGS=$(test_cxxflags)
Chpr2_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpr2_test_qd_opt_SOURCES=../common/Chpr2.test.cpp $(test_srcdepends)
Chpr2_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Chpr2_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Crotg_test_qd_SOURCES=../common/Crotg.test.cpp $(test_srcdepends)
Crotg_test_qd_CXXFLAGS=$(test_cxxflags)
Crotg_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Crotg_test_qd_opt_SOURCES=../common/Crotg.test.cpp $(test_srcdepends)
Crotg_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Crotg_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cscal_test_qd_SOURCES=../common/Cscal.test.cpp $(test_srcdepends)
Cscal_test_qd_CXXFLAGS=$(test_cxxflags)
Cscal_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cscal_test_qd_opt_SOURCES=../common/Cscal.test.cpp $(test_srcdepends)
Cscal_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Cscal_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cswap_test_qd_SOURCES=../common/Cswap.test.cpp $(test_srcdepends)
Cswap_test_qd_CXXFLAGS=$(test_cxxflags)
Cswap_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cswap_test_qd_opt_SOURCES=../common/Cswap.test.cpp $(test_srcdepends)
Cswap_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Cswap_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csymm_test_qd_SOURCES=../common/Csymm.test.cpp $(test_srcdepends)
Csymm_test_qd_CXXFLAGS=$(test_cxxflags)
Csymm_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csymm_test_qd_opt_SOURCES=../common/Csymm.test.cpp $(test_srcdepends)
Csymm_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Csymm_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyr2k_test_qd_SOURCES=../common/Csyr2k.test.cpp $(test_srcdepends)
Csyr2k_test_qd_CXXFLAGS=$(test_cxxflags)
Csyr2k_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csyr2k_test_qd_opt_SOURCES=../common/Csyr2k.test.cpp $(test_srcdepends)
Csyr2k_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Csyr2k_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyrk_test_qd_SOURCES=../common/Csyrk.test.cpp $(test_srcdepends)
Csyrk_test_qd_CXXFLAGS=$(test_cxxflags)
Csyrk_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csyrk_test_qd_opt_SOURCES=../common/Csyrk.test.cpp $(test_srcdepends)
Csyrk_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Csyrk_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbmv_test_qd_SOURCES=../common/Ctbmv.test.cpp $(test_srcdepends)
Ctbmv_test_qd_CXXFLAGS=$(test_cxxflags)
Ctbmv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctbmv_test_qd_opt_SOURCES=../common/Ctbmv.test.cpp $(test_srcdepends)
Ctbmv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Ctbmv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbsv_test_qd_SOURCES=../common/Ctbsv.test.cpp $(test_srcdepends)
Ctbsv_test_qd_CXXFLAGS=$(test_cxxflags)
Ctbsv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctbsv_test_qd_opt_SOURCES=../common/Ctbsv.test.cpp $(test_srcdepends)
Ctbsv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Ctbsv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpmv_test_qd_SOURCES=../common/Ctpmv.test.cpp $(test_srcdepends)
Ctpmv_test_qd_CXXFLAGS=$(test_cxxflags)
Ctpmv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctpmv_test_qd_opt_SOURCES=../common/Ctpmv.test.cpp $(test_srcdepends)
Ctpmv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Ctpmv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpsv_test_qd_SOURCES=../common/Ctpsv.test.cpp $(test_srcdepends)
Ctpsv_test_qd_CXXFLAGS=$(test_cxxflags)
Ctpsv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctpsv_test_qd_opt_SOURCES=../common/Ctpsv.test.cpp $(test_srcdepends)
Ctpsv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Ctpsv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmm_test_qd_SOURCES=../common/Ctrmm.test.cpp $(test_srcdepends)
Ctrmm_test_qd_CXXFLAGS=$(test_cxxflags)
Ctrmm_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrmm_test_qd_opt_SOURCES=../common/Ctrmm.test.cpp $(test_srcdepends)
Ctrmm_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Ctrmm_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmv_test_qd_SOURCES=../common/Ctrmv.test.cpp $(test_srcdepends)
Ctrmv_test_qd_CXXFLAGS=$(test_cxxflags)
Ctrmv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrmv_test_qd_opt_SOURCES=../common/Ctrmv.test.cpp $(test_srcdepends)
Ctrmv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Ctrmv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsm_test_qd_SOURCES=../common/Ctrsm.test.cpp $(test_srcdepends)
Ctrsm_test_qd_CXXFLAGS=$(test_cxxflags)
Ctrsm_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrsm_test_qd_opt_SOURCES=../common/Ctrsm.test.cpp $(test_srcdepends)
Ctrsm_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Ctrsm_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsv_test_qd_SOURCES=../common/Ctrsv.test.cpp $(test_srcdepends)
Ctrsv_test_qd_CXXFLAGS=$(test_cxxflags)
Ctrsv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrsv_test_qd_opt_SOURCES=../common/Ctrsv.test.cpp $(test_srcdepends)
Ctrsv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Ctrsv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCabs1_test_qd_opt_SOURCES=../common/RCabs1.test.cpp $(test_srcdepends)
RCabs1_test_qd_opt_CXXFLAGS=$(test_cxxflags)
RCabs1_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCabs1_test_qd_SOURCES=../common/RCabs1.test.cpp $(test_srcdepends)
RCabs1_test_qd_CXXFLAGS=$(test_cxxflags)
RCabs1_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

RCasum_test_qd_opt_SOURCES=../common/RCasum.test.cpp $(test_srcdepends)
RCasum_test_qd_opt_CXXFLAGS=$(test_cxxflags)
RCasum_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCasum_test_qd_SOURCES=../common/RCasum.test.cpp $(test_srcdepends)
RCasum_test_qd_CXXFLAGS=$(test_cxxflags)
RCasum_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

RCnrm2_test_qd_opt_SOURCES=../common/RCnrm2.test.cpp $(test_srcdepends)
RCnrm2_test_qd_opt_CXXFLAGS=$(test_cxxflags)
RCnrm2_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCnrm2_test_qd_SOURCES=../common/RCnrm2.test.cpp $(test_srcdepends)
RCnrm2_test_qd_CXXFLAGS=$(test_cxxflags)
RCnrm2_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rasum_test_qd_opt_SOURCES=../common/Rasum.test.cpp $(test_srcdepends)
Rasum_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rasum_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rasum_test_qd_SOURCES=../common/Rasum.test.cpp $(test_srcdepends)
Rasum_test_qd_CXXFLAGS=$(test_cxxflags)
Rasum_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Raxpy_test_qd_opt_SOURCES=../common/Raxpy.test.cpp $(test_srcdepends)
Raxpy_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Raxpy_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Raxpy_test_qd_SOURCES=../common/Raxpy.test.cpp $(test_srcdepends)
Raxpy_test_qd_CXXFLAGS=$(test_cxxflags)
Raxpy_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rcopy_test_qd_opt_SOURCES=../common/Rcopy.test.cpp $(test_srcdepends)
Rcopy_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rcopy_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rcopy_test_qd_SOURCES=../common/Rcopy.test.cpp $(test_srcdepends)
Rcopy_test_qd_CXXFLAGS=$(test_cxxflags)
Rcopy_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rdot_test_qd_opt_SOURCES=../common/Rdot.test.cpp $(test_srcdepends)
Rdot_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rdot_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rdot_test_qd_SOURCES=../common/Rdot.test.cpp $(test_srcdepends)
Rdot_test_qd_CXXFLAGS=$(test_cxxflags)
Rdot_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgbmv_test_qd_opt_SOURCES=../common/Rgbmv.test.cpp $(test_srcdepends)
Rgbmv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rgbmv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgbmv_test_qd_SOURCES=../common/Rgbmv.test.cpp $(test_srcdepends)
Rgbmv_test_qd_CXXFLAGS=$(test_cxxflags)
Rgbmv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgemm_test_qd_opt_SOURCES=../common/Rgemm.test.cpp $(test_srcdepends)
Rgemm_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rgemm_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemm_test_qd_SOURCES=../common/Rgemm.test.cpp $(test_srcdepends)
Rgemm_test_qd_CXXFLAGS=$(test_cxxflags)
Rgemm_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgemv_test_qd_opt_SOURCES=../common/Rgemv.test.cpp $(test_srcdepends)
Rgemv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rgemv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemv_test_qd_SOURCES=../common/Rgemv.test.cpp $(test_srcdepends)
Rgemv_test_qd_CXXFLAGS=$(test_cxxflags)
Rgemv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rger_test_qd_opt_SOURCES=../common/Rger.test.cpp $(test_srcdepends)
Rger_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rger_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rger_test_qd_SOURCES=../common/Rger.test.cpp $(test_srcdepends)
Rger_test_qd_CXXFLAGS=$(test_cxxflags)
Rger_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rnrm2_test_qd_opt_SOURCES=../common/Rnrm2.test.cpp $(test_srcdepends)
Rnrm2_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rnrm2_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rnrm2_test_qd_SOURCES=../common/Rnrm2.test.cpp $(test_srcdepends)
Rnrm2_test_qd_CXXFLAGS=$(test_cxxflags)
Rnrm2_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrot_test_qd_opt_SOURCES=../common/Rrot.test.cpp $(test_srcdepends)
Rrot_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rrot_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrot_test_qd_SOURCES=../common/Rrot.test.cpp $(test_srcdepends)
Rrot_test_qd_CXXFLAGS=$(test_cxxflags)
Rrot_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrotg_test_qd_opt_SOURCES=../common/Rrotg.test.cpp $(test_srcdepends)
Rrotg_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rrotg_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotg_test_qd_SOURCES=../common/Rrotg.test.cpp $(test_srcdepends)
Rrotg_test_qd_CXXFLAGS=$(test_cxxflags)
Rrotg_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrotm_test_qd_opt_SOURCES=../common/Rrotm.test.cpp $(test_srcdepends)
Rrotm_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rrotm_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotm_test_qd_SOURCES=../common/Rrotm.test.cpp $(test_srcdepends)
Rrotm_test_qd_CXXFLAGS=$(test_cxxflags)
Rrotm_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsbmv_test_qd_opt_SOURCES=../common/Rsbmv.test.cpp $(test_srcdepends)
Rsbmv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rsbmv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsbmv_test_qd_SOURCES=../common/Rsbmv.test.cpp $(test_srcdepends)
Rsbmv_test_qd_CXXFLAGS=$(test_cxxflags)
Rsbmv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rscal_test_qd_opt_SOURCES=../common/Rscal.test.cpp $(test_srcdepends)
Rscal_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rscal_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rscal_test_qd_SOURCES=../common/Rscal.test.cpp $(test_srcdepends)
Rscal_test_qd_CXXFLAGS=$(test_cxxflags)
Rscal_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspmv_test_qd_opt_SOURCES=../common/Rspmv.test.cpp $(test_srcdepends)
Rspmv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rspmv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspmv_test_qd_SOURCES=../common/Rspmv.test.cpp $(test_srcdepends)
Rspmv_test_qd_CXXFLAGS=$(test_cxxflags)
Rspmv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspr_test_qd_opt_SOURCES=../common/Rspr.test.cpp $(test_srcdepends)
Rspr_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rspr_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr_test_qd_SOURCES=../common/Rspr.test.cpp $(test_srcdepends)
Rspr_test_qd_CXXFLAGS=$(test_cxxflags)
Rspr_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspr2_test_qd_opt_SOURCES=../common/Rspr2.test.cpp $(test_srcdepends)
Rspr2_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rspr2_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr2_test_qd_SOURCES=../common/Rspr2.test.cpp $(test_srcdepends)
Rspr2_test_qd_CXXFLAGS=$(test_cxxflags)
Rspr2_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rswap_test_qd_opt_SOURCES=../common/Rswap.test.cpp $(test_srcdepends)
Rswap_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rswap_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rswap_test_qd_SOURCES=../common/Rswap.test.cpp $(test_srcdepends)
Rswap_test_qd_CXXFLAGS=$(test_cxxflags)
Rswap_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsymm_test_qd_opt_SOURCES=../common/Rsymm.test.cpp $(test_srcdepends)
Rsymm_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rsymm_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymm_test_qd_SOURCES=../common/Rsymm.test.cpp $(test_srcdepends)
Rsymm_test_qd_CXXFLAGS=$(test_cxxflags)
Rsymm_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsymv_test_qd_opt_SOURCES=../common/Rsymv.test.cpp $(test_srcdepends)
Rsymv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rsymv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymv_test_qd_SOURCES=../common/Rsymv.test.cpp $(test_srcdepends)
Rsymv_test_qd_CXXFLAGS=$(test_cxxflags)
Rsymv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr_test_qd_opt_SOURCES=../common/Rsyr.test.cpp $(test_srcdepends)
Rsyr_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rsyr_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr_test_qd_SOURCES=../common/Rsyr.test.cpp $(test_srcdepends)
Rsyr_test_qd_CXXFLAGS=$(test_cxxflags)
Rsyr_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr2_test_qd_opt_SOURCES=../common/Rsyr2.test.cpp $(test_srcdepends)
Rsyr2_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rsyr2_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2_test_qd_SOURCES=../common/Rsyr2.test.cpp $(test_srcdepends)
Rsyr2_test_qd_CXXFLAGS=$(test_cxxflags)
Rsyr2_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr2k_test_qd_opt_SOURCES=../common/Rsyr2k.test.cpp $(test_srcdepends)
Rsyr2k_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rsyr2k_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2k_test_qd_SOURCES=../common/Rsyr2k.test.cpp $(test_srcdepends)
Rsyr2k_test_qd_CXXFLAGS=$(test_cxxflags)
Rsyr2k_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyrk_test_qd_opt_SOURCES=../common/Rsyrk.test.cpp $(test_srcdepends)
Rsyrk_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rsyrk_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyrk_test_qd_SOURCES=../common/Rsyrk.test.cpp $(test_srcdepends)
Rsyrk_test_qd_CXXFLAGS=$(test_cxxflags)
Rsyrk_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtbmv_test_qd_opt_SOURCES=../common/Rtbmv.test.cpp $(test_srcdepends)
Rtbmv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rtbmv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbmv_test_qd_SOURCES=../common/Rtbmv.test.cpp $(test_srcdepends)
Rtbmv_test_qd_CXXFLAGS=$(test_cxxflags)
Rtbmv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtbsv_test_qd_opt_SOURCES=../common/Rtbsv.test.cpp $(test_srcdepends)
Rtbsv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rtbsv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbsv_test_qd_SOURCES=../common/Rtbsv.test.cpp $(test_srcdepends)
Rtbsv_test_qd_CXXFLAGS=$(test_cxxflags)
Rtbsv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtpmv_test_qd_opt_SOURCES=../common/Rtpmv.test.cpp $(test_srcdepends)
Rtpmv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rtpmv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpmv_test_qd_SOURCES=../common/Rtpmv.test.cpp $(test_srcdepends)
Rtpmv_test_qd_CXXFLAGS=$(test_cxxflags)
Rtpmv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtpsv_test_qd_opt_SOURCES=../common/Rtpsv.test.cpp $(test_srcdepends)
Rtpsv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rtpsv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpsv_test_qd_SOURCES=../common/Rtpsv.test.cpp $(test_srcdepends)
Rtpsv_test_qd_CXXFLAGS=$(test_cxxflags)
Rtpsv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrmm_test_qd_opt_SOURCES=../common/Rtrmm.test.cpp $(test_srcdepends)
Rtrmm_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rtrmm_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmm_test_qd_SOURCES=../common/Rtrmm.test.cpp $(test_srcdepends)
Rtrmm_test_qd_CXXFLAGS=$(test_cxxflags)
Rtrmm_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrmv_test_qd_opt_SOURCES=../common/Rtrmv.test.cpp $(test_srcdepends)
Rtrmv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rtrmv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmv_test_qd_SOURCES=../common/Rtrmv.test.cpp $(test_srcdepends)
Rtrmv_test_qd_CXXFLAGS=$(test_cxxflags)
Rtrmv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrsm_test_qd_opt_SOURCES=../common/Rtrsm.test.cpp $(test_srcdepends)
Rtrsm_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rtrsm_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsm_test_qd_SOURCES=../common/Rtrsm.test.cpp $(test_srcdepends)
Rtrsm_test_qd_CXXFLAGS=$(test_cxxflags)
Rtrsm_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrsv_test_qd_opt_SOURCES=../common/Rtrsv.test.cpp $(test_srcdepends)
Rtrsv_test_qd_opt_CXXFLAGS=$(test_cxxflags)
Rtrsv_test_qd_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsv_test_qd_SOURCES=../common/Rtrsv.test.cpp $(test_srcdepends)
Rtrsv_test_qd_CXXFLAGS=$(test_cxxflags)
Rtrsv_test_qd_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

