#CTBSV, CTPSV might have bugs or intrinsically inaccurate
#Rnrm2 should be checked more
check_PROGRAMS = $(mpblas_mpfr_test_PROGRAMS) 
TESTS = $(check_PROGRAMS)
mpblas_mpfr_test_PROGRAMS = \
arithmetic.test_mpfr \
complex.test_mpfr \
iCamax.test.mpfr_opt iCamax.test.mpfr   \
iRamax.test.mpfr_opt iRamax.test.mpfr   \
Mlsame.test.mpfr_opt Mlsame.test.mpfr   \
Mxerbla.test.mpfr_opt Mxerbla.test.mpfr \
CRrot.test.mpfr_opt  CRrot.test.mpfr    \
CRscal.test.mpfr_opt CRscal.test.mpfr   \
Caxpy.test.mpfr_opt  Caxpy.test.mpfr    \
Ccopy.test.mpfr_opt  Ccopy.test.mpfr    \
Cdotc.test.mpfr_opt  Cdotc.test.mpfr	\
Cdotu.test.mpfr_opt  Cdotu.test.mpfr	\
Cgbmv.test.mpfr_opt  Cgbmv.test.mpfr	\
Cgemm.test.mpfr_opt  Cgemm.test.mpfr	\
Cgemv.test.mpfr_opt  Cgemv.test.mpfr	\
Cgerc.test.mpfr_opt  Cgerc.test.mpfr	\
Cgeru.test.mpfr_opt  Cgeru.test.mpfr    \
Chbmv.test.mpfr_opt  Chbmv.test.mpfr	\
Chemm.test.mpfr_opt  Chemm.test.mpfr	\
Chemv.test.mpfr_opt  Chemv.test.mpfr	\
Cher.test.mpfr_opt   Cher.test.mpfr 	\
Cher2.test.mpfr_opt  Cher2.test.mpfr	\
Cher2k.test.mpfr_opt Cher2k.test.mpfr	\
Cherk.test.mpfr_opt  Cherk.test.mpfr	\
Chpmv.test.mpfr_opt  Chpmv.test.mpfr	\
Chpr.test.mpfr_opt   Chpr.test.mpfr	\
Chpr2.test.mpfr_opt  Chpr2.test.mpfr	\
Crotg.test.mpfr_opt  Crotg.test.mpfr	\
Cscal.test.mpfr_opt  Cscal.test.mpfr	\
Cswap.test.mpfr_opt  Cswap.test.mpfr	\
Csymm.test.mpfr_opt  Csymm.test.mpfr    \
Csyr2k.test.mpfr_opt Csyr2k.test.mpfr	\
Csyrk.test.mpfr_opt  Csyrk.test.mpfr	\
Ctbmv.test.mpfr_opt  Ctbmv.test.mpfr	\
Ctbsv.test.mpfr_opt  Ctbsv.test.mpfr	\
Ctpmv.test.mpfr_opt  Ctpmv.test.mpfr	\
Ctpsv.test.mpfr_opt  Ctpsv.test.mpfr	\
Ctrmm.test.mpfr_opt  Ctrmm.test.mpfr	\
Ctrmv.test.mpfr_opt  Ctrmv.test.mpfr	\
Ctrsm.test.mpfr_opt  Ctrsm.test.mpfr	\
Ctrsv.test.mpfr_opt  Ctrsv.test.mpfr	\
RCabs1.test.mpfr_opt RCabs1.test.mpfr	\
RCasum.test.mpfr_opt RCasum.test.mpfr	\
RCnrm2.test.mpfr_opt RCnrm2.test.mpfr	\
Rasum.test.mpfr_opt  Rasum.test.mpfr    \
Raxpy.test.mpfr_opt  Raxpy.test.mpfr	\
Rcopy.test.mpfr_opt  Rcopy.test.mpfr	\
Rdot.test.mpfr_opt   Rdot.test.mpfr	\
Rgbmv.test.mpfr_opt  Rgbmv.test.mpfr	\
Rgemm.test.mpfr_opt  Rgemm.test.mpfr	\
Rgemv.test.mpfr_opt  Rgemv.test.mpfr	\
Rger.test.mpfr_opt   Rger.test.mpfr	\
Rnrm2.test.mpfr_opt  Rnrm2.test.mpfr	\
Rrot.test.mpfr_opt   Rrot.test.mpfr	\
Rrotg.test.mpfr_opt  Rrotg.test.mpfr	\
Rrotm.test.mpfr_opt  Rrotm.test.mpfr	\
Rsbmv.test.mpfr_opt  Rsbmv.test.mpfr	\
Rscal.test.mpfr_opt  Rscal.test.mpfr	\
Rspmv.test.mpfr_opt  Rspmv.test.mpfr    \
Rspr.test.mpfr_opt   Rspr.test.mpfr	\
Rspr2.test.mpfr_opt  Rspr2.test.mpfr	\
Rswap.test.mpfr_opt  Rswap.test.mpfr	\
Rsymm.test.mpfr_opt  Rsymm.test.mpfr	\
Rsymv.test.mpfr_opt  Rsymv.test.mpfr	\
Rsyr.test.mpfr_opt   Rsyr.test.mpfr	\
Rsyr2.test.mpfr_opt  Rsyr2.test.mpfr	\
Rsyr2k.test.mpfr_opt Rsyr2k.test.mpfr	\
Rsyrk.test.mpfr_opt  Rsyrk.test.mpfr	\
Rtbmv.test.mpfr_opt  Rtbmv.test.mpfr	\
Rtbsv.test.mpfr_opt  Rtbsv.test.mpfr	\
Rtpmv.test.mpfr_opt  Rtpmv.test.mpfr	\
Rtpsv.test.mpfr_opt  Rtpsv.test.mpfr	\
Rtrmm.test.mpfr_opt  Rtrmm.test.mpfr    \
Rtrmv.test.mpfr_opt  Rtrmv.test.mpfr	\
Rtrsm.test.mpfr_opt  Rtrsm.test.mpfr	\
Rtrsv.test.mpfr_opt  Rtrsv.test.mpfr

mpblas_mpfr_testdir = $(prefix)/lib/$(target)/mplapack/test/compare/mpfr

mplibs=-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp
refblas=-L$(top_builddir)/external/i/LAPACK/lib -lblas $(FCLIBS)

test_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_MPFR___ -D___MPLAPACK_INTERNAL___

if !IS_MACOS
test_srcdepends = ../common/xerbla.cpp ../common/mplapack.test.cpp ../common/Mxerbla.override.cpp
endif

if IS_MACOS
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_MPFR___
libmxerbla_override_la_SOURCES =  ../common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS= -Wl,-flat_namespace,-undefined,dynamic_lookup
test_srcdepends = ../common/xerbla.cpp ../common/mplapack.test.cpp
mplibs+=-Wl,-flat_namespace,-undefined,dynamic_lookup
TESTS_ENVIRONMENT=DYLD_INSERT_LIBRARIES=.libs/libmxerbla_override.0.dylib
endif

if IS_MINGW
$(check_PROGRAMS): libmxerbla_override.la
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_MPFR___
libmxerbla_override_la_SOURCES =  ../common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS = -no-undefined
endif

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,.libs/libmxerbla_override.a,$(top_builddir)/mpblas/reference/.libs/libmpblas_mpfr.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -lmpblas_mpfr
endif

if IS_MINGW
opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,.libs/libmxerbla_override.a,$(top_builddir)/mpblas/optimized/mpfr/.libs/libmpblas_mpfr_opt.a,--no-whole-archive
else
opt_libdepends = -L$(top_builddir)/mpblas/optimized/mpfr -lmpblas_mpfr_opt -L$(top_builddir)/mpblas/reference -lmpblas_mpfr
endif

##
arithmetic_test_mpfr_SOURCES=../common/arithmetic.test.cpp $(test_srcdepends)
arithmetic_test_mpfr_CXXFLAGS=$(test_cxxflags)
arithmetic_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

complex_test_mpfr_SOURCES=../common/complex.test.cpp $(test_srcdepends)
complex_test_mpfr_CXXFLAGS=$(test_cxxflags)
complex_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

iCamax_test_mpfr_opt_SOURCES=../common/iCamax.test.cpp $(test_srcdepends)
iCamax_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
iCamax_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iCamax_test_mpfr_SOURCES=../common/iCamax.test.cpp $(test_srcdepends)
iCamax_test_mpfr_CXXFLAGS=$(test_cxxflags)
iCamax_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

iRamax_test_mpfr_opt_SOURCES=../common/iRamax.test.cpp $(test_srcdepends)
iRamax_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
iRamax_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iRamax_test_mpfr_SOURCES=../common/iRamax.test.cpp $(test_srcdepends)
iRamax_test_mpfr_CXXFLAGS=$(test_cxxflags)
iRamax_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Mlsame_test_mpfr_opt_SOURCES=../common/Mlsame.test.cpp $(test_srcdepends)
Mlsame_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Mlsame_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mlsame_test_mpfr_SOURCES=../common/Mlsame.test.cpp $(test_srcdepends)
Mlsame_test_mpfr_CXXFLAGS=$(test_cxxflags)
Mlsame_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Mxerbla_test_mpfr_opt_SOURCES=../common/Mxerbla.test.cpp $(test_srcdepends)
Mxerbla_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Mxerbla_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mxerbla_test_mpfr_SOURCES=../common/Mxerbla.test.cpp $(test_srcdepends)
Mxerbla_test_mpfr_CXXFLAGS=$(test_cxxflags)
Mxerbla_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

CRrot_test_mpfr_opt_SOURCES=../common/CRrot.test.cpp $(test_srcdepends)
CRrot_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
CRrot_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRrot_test_mpfr_SOURCES=../common/CRrot.test.cpp $(test_srcdepends)
CRrot_test_mpfr_CXXFLAGS=$(test_cxxflags)
CRrot_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

CRscal_test_mpfr_opt_SOURCES=../common/CRscal.test.cpp $(test_srcdepends)
CRscal_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
CRscal_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRscal_test_mpfr_SOURCES=../common/CRscal.test.cpp $(test_srcdepends)
CRscal_test_mpfr_CXXFLAGS=$(test_cxxflags)
CRscal_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Caxpy_test_mpfr_SOURCES=../common/Caxpy.test.cpp $(test_srcdepends)
Caxpy_test_mpfr_CXXFLAGS=$(test_cxxflags)
Caxpy_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Caxpy_test_mpfr_opt_SOURCES=../common/Caxpy.test.cpp $(test_srcdepends)
Caxpy_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Caxpy_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ccopy_test_mpfr_SOURCES=../common/Ccopy.test.cpp $(test_srcdepends)
Ccopy_test_mpfr_CXXFLAGS=$(test_cxxflags)
Ccopy_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ccopy_test_mpfr_opt_SOURCES=../common/Ccopy.test.cpp $(test_srcdepends)
Ccopy_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Ccopy_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotc_test_mpfr_SOURCES=../common/Cdotc.test.cpp $(test_srcdepends)
Cdotc_test_mpfr_CXXFLAGS=$(test_cxxflags)
Cdotc_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cdotc_test_mpfr_opt_SOURCES=../common/Cdotc.test.cpp $(test_srcdepends)
Cdotc_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Cdotc_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotu_test_mpfr_SOURCES=../common/Cdotu.test.cpp $(test_srcdepends)
Cdotu_test_mpfr_CXXFLAGS=$(test_cxxflags)
Cdotu_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cdotu_test_mpfr_opt_SOURCES=../common/Cdotu.test.cpp $(test_srcdepends)
Cdotu_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Cdotu_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgbmv_test_mpfr_SOURCES=../common/Cgbmv.test.cpp $(test_srcdepends)
Cgbmv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Cgbmv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgbmv_test_mpfr_opt_SOURCES=../common/Cgbmv.test.cpp $(test_srcdepends)
Cgbmv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Cgbmv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemm_test_mpfr_SOURCES=../common/Cgemm.test.cpp $(test_srcdepends)
Cgemm_test_mpfr_CXXFLAGS=$(test_cxxflags)
Cgemm_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgemm_test_mpfr_opt_SOURCES=../common/Cgemm.test.cpp $(test_srcdepends)
Cgemm_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Cgemm_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemv_test_mpfr_SOURCES=../common/Cgemv.test.cpp $(test_srcdepends)
Cgemv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Cgemv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgemv_test_mpfr_opt_SOURCES=../common/Cgemv.test.cpp $(test_srcdepends)
Cgemv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Cgemv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgerc_test_mpfr_SOURCES=../common/Cgerc.test.cpp $(test_srcdepends)
Cgerc_test_mpfr_CXXFLAGS=$(test_cxxflags)
Cgerc_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgerc_test_mpfr_opt_SOURCES=../common/Cgerc.test.cpp $(test_srcdepends)
Cgerc_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Cgerc_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgeru_test_mpfr_SOURCES=../common/Cgeru.test.cpp $(test_srcdepends)
Cgeru_test_mpfr_CXXFLAGS=$(test_cxxflags)
Cgeru_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgeru_test_mpfr_opt_SOURCES=../common/Cgeru.test.cpp $(test_srcdepends)
Cgeru_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Cgeru_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chbmv_test_mpfr_SOURCES=../common/Chbmv.test.cpp $(test_srcdepends)
Chbmv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Chbmv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chbmv_test_mpfr_opt_SOURCES=../common/Chbmv.test.cpp $(test_srcdepends)
Chbmv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Chbmv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemm_test_mpfr_SOURCES=../common/Chemm.test.cpp $(test_srcdepends)
Chemm_test_mpfr_CXXFLAGS=$(test_cxxflags)
Chemm_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chemm_test_mpfr_opt_SOURCES=../common/Chemm.test.cpp $(test_srcdepends)
Chemm_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Chemm_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemv_test_mpfr_SOURCES=../common/Chemv.test.cpp $(test_srcdepends)
Chemv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Chemv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chemv_test_mpfr_opt_SOURCES=../common/Chemv.test.cpp $(test_srcdepends)
Chemv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Chemv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher_test_mpfr_SOURCES=../common/Cher.test.cpp $(test_srcdepends)
Cher_test_mpfr_CXXFLAGS=$(test_cxxflags)
Cher_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher_test_mpfr_opt_SOURCES=../common/Cher.test.cpp $(test_srcdepends)
Cher_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Cher_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2_test_mpfr_SOURCES=../common/Cher2.test.cpp $(test_srcdepends)
Cher2_test_mpfr_CXXFLAGS=$(test_cxxflags)
Cher2_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher2_test_mpfr_opt_SOURCES=../common/Cher2.test.cpp $(test_srcdepends)
Cher2_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Cher2_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2k_test_mpfr_SOURCES=../common/Cher2k.test.cpp $(test_srcdepends)
Cher2k_test_mpfr_CXXFLAGS=$(test_cxxflags)
Cher2k_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher2k_test_mpfr_opt_SOURCES=../common/Cher2k.test.cpp $(test_srcdepends)
Cher2k_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Cher2k_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cherk_test_mpfr_SOURCES=../common/Cherk.test.cpp $(test_srcdepends)
Cherk_test_mpfr_CXXFLAGS=$(test_cxxflags)
Cherk_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cherk_test_mpfr_opt_SOURCES=../common/Cherk.test.cpp $(test_srcdepends)
Cherk_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Cherk_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpmv_test_mpfr_SOURCES=../common/Chpmv.test.cpp $(test_srcdepends)
Chpmv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Chpmv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpmv_test_mpfr_opt_SOURCES=../common/Chpmv.test.cpp $(test_srcdepends)
Chpmv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Chpmv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr_test_mpfr_SOURCES=../common/Chpr.test.cpp $(test_srcdepends)
Chpr_test_mpfr_CXXFLAGS=$(test_cxxflags)
Chpr_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpr_test_mpfr_opt_SOURCES=../common/Chpr.test.cpp $(test_srcdepends)
Chpr_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Chpr_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr2_test_mpfr_SOURCES=../common/Chpr2.test.cpp $(test_srcdepends)
Chpr2_test_mpfr_CXXFLAGS=$(test_cxxflags)
Chpr2_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpr2_test_mpfr_opt_SOURCES=../common/Chpr2.test.cpp $(test_srcdepends)
Chpr2_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Chpr2_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Crotg_test_mpfr_SOURCES=../common/Crotg.test.cpp $(test_srcdepends)
Crotg_test_mpfr_CXXFLAGS=$(test_cxxflags)
Crotg_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Crotg_test_mpfr_opt_SOURCES=../common/Crotg.test.cpp $(test_srcdepends)
Crotg_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Crotg_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cscal_test_mpfr_SOURCES=../common/Cscal.test.cpp $(test_srcdepends)
Cscal_test_mpfr_CXXFLAGS=$(test_cxxflags)
Cscal_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cscal_test_mpfr_opt_SOURCES=../common/Cscal.test.cpp $(test_srcdepends)
Cscal_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Cscal_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cswap_test_mpfr_SOURCES=../common/Cswap.test.cpp $(test_srcdepends)
Cswap_test_mpfr_CXXFLAGS=$(test_cxxflags)
Cswap_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cswap_test_mpfr_opt_SOURCES=../common/Cswap.test.cpp $(test_srcdepends)
Cswap_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Cswap_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csymm_test_mpfr_SOURCES=../common/Csymm.test.cpp $(test_srcdepends)
Csymm_test_mpfr_CXXFLAGS=$(test_cxxflags)
Csymm_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csymm_test_mpfr_opt_SOURCES=../common/Csymm.test.cpp $(test_srcdepends)
Csymm_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Csymm_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyr2k_test_mpfr_SOURCES=../common/Csyr2k.test.cpp $(test_srcdepends)
Csyr2k_test_mpfr_CXXFLAGS=$(test_cxxflags)
Csyr2k_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csyr2k_test_mpfr_opt_SOURCES=../common/Csyr2k.test.cpp $(test_srcdepends)
Csyr2k_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Csyr2k_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyrk_test_mpfr_SOURCES=../common/Csyrk.test.cpp $(test_srcdepends)
Csyrk_test_mpfr_CXXFLAGS=$(test_cxxflags)
Csyrk_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csyrk_test_mpfr_opt_SOURCES=../common/Csyrk.test.cpp $(test_srcdepends)
Csyrk_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Csyrk_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbmv_test_mpfr_SOURCES=../common/Ctbmv.test.cpp $(test_srcdepends)
Ctbmv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Ctbmv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctbmv_test_mpfr_opt_SOURCES=../common/Ctbmv.test.cpp $(test_srcdepends)
Ctbmv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Ctbmv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbsv_test_mpfr_SOURCES=../common/Ctbsv.test.cpp $(test_srcdepends)
Ctbsv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Ctbsv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctbsv_test_mpfr_opt_SOURCES=../common/Ctbsv.test.cpp $(test_srcdepends)
Ctbsv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Ctbsv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpmv_test_mpfr_SOURCES=../common/Ctpmv.test.cpp $(test_srcdepends)
Ctpmv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Ctpmv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctpmv_test_mpfr_opt_SOURCES=../common/Ctpmv.test.cpp $(test_srcdepends)
Ctpmv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Ctpmv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpsv_test_mpfr_SOURCES=../common/Ctpsv.test.cpp $(test_srcdepends)
Ctpsv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Ctpsv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctpsv_test_mpfr_opt_SOURCES=../common/Ctpsv.test.cpp $(test_srcdepends)
Ctpsv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Ctpsv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmm_test_mpfr_SOURCES=../common/Ctrmm.test.cpp $(test_srcdepends)
Ctrmm_test_mpfr_CXXFLAGS=$(test_cxxflags)
Ctrmm_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrmm_test_mpfr_opt_SOURCES=../common/Ctrmm.test.cpp $(test_srcdepends)
Ctrmm_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Ctrmm_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmv_test_mpfr_SOURCES=../common/Ctrmv.test.cpp $(test_srcdepends)
Ctrmv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Ctrmv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrmv_test_mpfr_opt_SOURCES=../common/Ctrmv.test.cpp $(test_srcdepends)
Ctrmv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Ctrmv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsm_test_mpfr_SOURCES=../common/Ctrsm.test.cpp $(test_srcdepends)
Ctrsm_test_mpfr_CXXFLAGS=$(test_cxxflags)
Ctrsm_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrsm_test_mpfr_opt_SOURCES=../common/Ctrsm.test.cpp $(test_srcdepends)
Ctrsm_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Ctrsm_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsv_test_mpfr_SOURCES=../common/Ctrsv.test.cpp $(test_srcdepends)
Ctrsv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Ctrsv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrsv_test_mpfr_opt_SOURCES=../common/Ctrsv.test.cpp $(test_srcdepends)
Ctrsv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Ctrsv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCabs1_test_mpfr_opt_SOURCES=../common/RCabs1.test.cpp $(test_srcdepends)
RCabs1_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
RCabs1_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCabs1_test_mpfr_SOURCES=../common/RCabs1.test.cpp $(test_srcdepends)
RCabs1_test_mpfr_CXXFLAGS=$(test_cxxflags)
RCabs1_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

RCasum_test_mpfr_opt_SOURCES=../common/RCasum.test.cpp $(test_srcdepends)
RCasum_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
RCasum_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCasum_test_mpfr_SOURCES=../common/RCasum.test.cpp $(test_srcdepends)
RCasum_test_mpfr_CXXFLAGS=$(test_cxxflags)
RCasum_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

RCnrm2_test_mpfr_opt_SOURCES=../common/RCnrm2.test.cpp $(test_srcdepends)
RCnrm2_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
RCnrm2_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCnrm2_test_mpfr_SOURCES=../common/RCnrm2.test.cpp $(test_srcdepends)
RCnrm2_test_mpfr_CXXFLAGS=$(test_cxxflags)
RCnrm2_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rasum_test_mpfr_opt_SOURCES=../common/Rasum.test.cpp $(test_srcdepends)
Rasum_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rasum_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rasum_test_mpfr_SOURCES=../common/Rasum.test.cpp $(test_srcdepends)
Rasum_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rasum_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Raxpy_test_mpfr_opt_SOURCES=../common/Raxpy.test.cpp $(test_srcdepends)
Raxpy_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Raxpy_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Raxpy_test_mpfr_SOURCES=../common/Raxpy.test.cpp $(test_srcdepends)
Raxpy_test_mpfr_CXXFLAGS=$(test_cxxflags)
Raxpy_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rcopy_test_mpfr_opt_SOURCES=../common/Rcopy.test.cpp $(test_srcdepends)
Rcopy_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rcopy_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rcopy_test_mpfr_SOURCES=../common/Rcopy.test.cpp $(test_srcdepends)
Rcopy_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rcopy_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rdot_test_mpfr_opt_SOURCES=../common/Rdot.test.cpp $(test_srcdepends)
Rdot_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rdot_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rdot_test_mpfr_SOURCES=../common/Rdot.test.cpp $(test_srcdepends)
Rdot_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rdot_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgbmv_test_mpfr_opt_SOURCES=../common/Rgbmv.test.cpp $(test_srcdepends)
Rgbmv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rgbmv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgbmv_test_mpfr_SOURCES=../common/Rgbmv.test.cpp $(test_srcdepends)
Rgbmv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rgbmv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgemm_test_mpfr_opt_SOURCES=../common/Rgemm.test.cpp $(test_srcdepends)
Rgemm_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rgemm_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemm_test_mpfr_SOURCES=../common/Rgemm.test.cpp $(test_srcdepends)
Rgemm_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rgemm_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgemv_test_mpfr_opt_SOURCES=../common/Rgemv.test.cpp $(test_srcdepends)
Rgemv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rgemv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemv_test_mpfr_SOURCES=../common/Rgemv.test.cpp $(test_srcdepends)
Rgemv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rgemv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rger_test_mpfr_opt_SOURCES=../common/Rger.test.cpp $(test_srcdepends)
Rger_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rger_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rger_test_mpfr_SOURCES=../common/Rger.test.cpp $(test_srcdepends)
Rger_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rger_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rnrm2_test_mpfr_opt_SOURCES=../common/Rnrm2.test.cpp $(test_srcdepends)
Rnrm2_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rnrm2_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rnrm2_test_mpfr_SOURCES=../common/Rnrm2.test.cpp $(test_srcdepends)
Rnrm2_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rnrm2_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrot_test_mpfr_opt_SOURCES=../common/Rrot.test.cpp $(test_srcdepends)
Rrot_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rrot_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrot_test_mpfr_SOURCES=../common/Rrot.test.cpp $(test_srcdepends)
Rrot_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rrot_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrotg_test_mpfr_opt_SOURCES=../common/Rrotg.test.cpp $(test_srcdepends)
Rrotg_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rrotg_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotg_test_mpfr_SOURCES=../common/Rrotg.test.cpp $(test_srcdepends)
Rrotg_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rrotg_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrotm_test_mpfr_opt_SOURCES=../common/Rrotm.test.cpp $(test_srcdepends)
Rrotm_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rrotm_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotm_test_mpfr_SOURCES=../common/Rrotm.test.cpp $(test_srcdepends)
Rrotm_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rrotm_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsbmv_test_mpfr_opt_SOURCES=../common/Rsbmv.test.cpp $(test_srcdepends)
Rsbmv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rsbmv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsbmv_test_mpfr_SOURCES=../common/Rsbmv.test.cpp $(test_srcdepends)
Rsbmv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rsbmv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rscal_test_mpfr_opt_SOURCES=../common/Rscal.test.cpp $(test_srcdepends)
Rscal_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rscal_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rscal_test_mpfr_SOURCES=../common/Rscal.test.cpp $(test_srcdepends)
Rscal_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rscal_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspmv_test_mpfr_opt_SOURCES=../common/Rspmv.test.cpp $(test_srcdepends)
Rspmv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rspmv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspmv_test_mpfr_SOURCES=../common/Rspmv.test.cpp $(test_srcdepends)
Rspmv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rspmv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspr_test_mpfr_opt_SOURCES=../common/Rspr.test.cpp $(test_srcdepends)
Rspr_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rspr_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr_test_mpfr_SOURCES=../common/Rspr.test.cpp $(test_srcdepends)
Rspr_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rspr_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspr2_test_mpfr_opt_SOURCES=../common/Rspr2.test.cpp $(test_srcdepends)
Rspr2_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rspr2_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr2_test_mpfr_SOURCES=../common/Rspr2.test.cpp $(test_srcdepends)
Rspr2_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rspr2_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rswap_test_mpfr_opt_SOURCES=../common/Rswap.test.cpp $(test_srcdepends)
Rswap_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rswap_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rswap_test_mpfr_SOURCES=../common/Rswap.test.cpp $(test_srcdepends)
Rswap_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rswap_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsymm_test_mpfr_opt_SOURCES=../common/Rsymm.test.cpp $(test_srcdepends)
Rsymm_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rsymm_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymm_test_mpfr_SOURCES=../common/Rsymm.test.cpp $(test_srcdepends)
Rsymm_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rsymm_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsymv_test_mpfr_opt_SOURCES=../common/Rsymv.test.cpp $(test_srcdepends)
Rsymv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rsymv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymv_test_mpfr_SOURCES=../common/Rsymv.test.cpp $(test_srcdepends)
Rsymv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rsymv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr_test_mpfr_opt_SOURCES=../common/Rsyr.test.cpp $(test_srcdepends)
Rsyr_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rsyr_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr_test_mpfr_SOURCES=../common/Rsyr.test.cpp $(test_srcdepends)
Rsyr_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rsyr_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr2_test_mpfr_opt_SOURCES=../common/Rsyr2.test.cpp $(test_srcdepends)
Rsyr2_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rsyr2_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2_test_mpfr_SOURCES=../common/Rsyr2.test.cpp $(test_srcdepends)
Rsyr2_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rsyr2_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr2k_test_mpfr_opt_SOURCES=../common/Rsyr2k.test.cpp $(test_srcdepends)
Rsyr2k_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rsyr2k_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2k_test_mpfr_SOURCES=../common/Rsyr2k.test.cpp $(test_srcdepends)
Rsyr2k_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rsyr2k_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyrk_test_mpfr_opt_SOURCES=../common/Rsyrk.test.cpp $(test_srcdepends)
Rsyrk_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rsyrk_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyrk_test_mpfr_SOURCES=../common/Rsyrk.test.cpp $(test_srcdepends)
Rsyrk_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rsyrk_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtbmv_test_mpfr_opt_SOURCES=../common/Rtbmv.test.cpp $(test_srcdepends)
Rtbmv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rtbmv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbmv_test_mpfr_SOURCES=../common/Rtbmv.test.cpp $(test_srcdepends)
Rtbmv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rtbmv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtbsv_test_mpfr_opt_SOURCES=../common/Rtbsv.test.cpp $(test_srcdepends)
Rtbsv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rtbsv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbsv_test_mpfr_SOURCES=../common/Rtbsv.test.cpp $(test_srcdepends)
Rtbsv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rtbsv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtpmv_test_mpfr_opt_SOURCES=../common/Rtpmv.test.cpp $(test_srcdepends)
Rtpmv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rtpmv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpmv_test_mpfr_SOURCES=../common/Rtpmv.test.cpp $(test_srcdepends)
Rtpmv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rtpmv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtpsv_test_mpfr_opt_SOURCES=../common/Rtpsv.test.cpp $(test_srcdepends)
Rtpsv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rtpsv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpsv_test_mpfr_SOURCES=../common/Rtpsv.test.cpp $(test_srcdepends)
Rtpsv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rtpsv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrmm_test_mpfr_opt_SOURCES=../common/Rtrmm.test.cpp $(test_srcdepends)
Rtrmm_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rtrmm_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmm_test_mpfr_SOURCES=../common/Rtrmm.test.cpp $(test_srcdepends)
Rtrmm_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rtrmm_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrmv_test_mpfr_opt_SOURCES=../common/Rtrmv.test.cpp $(test_srcdepends)
Rtrmv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rtrmv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmv_test_mpfr_SOURCES=../common/Rtrmv.test.cpp $(test_srcdepends)
Rtrmv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rtrmv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrsm_test_mpfr_opt_SOURCES=../common/Rtrsm.test.cpp $(test_srcdepends)
Rtrsm_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rtrsm_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsm_test_mpfr_SOURCES=../common/Rtrsm.test.cpp $(test_srcdepends)
Rtrsm_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rtrsm_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrsv_test_mpfr_opt_SOURCES=../common/Rtrsv.test.cpp $(test_srcdepends)
Rtrsv_test_mpfr_opt_CXXFLAGS=$(test_cxxflags)
Rtrsv_test_mpfr_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsv_test_mpfr_SOURCES=../common/Rtrsv.test.cpp $(test_srcdepends)
Rtrsv_test_mpfr_CXXFLAGS=$(test_cxxflags)
Rtrsv_test_mpfr_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mpblas_mpfr_testdir) $(prefix)
endif
