#CTBSV, CTPSV might have bugs or intrinsically inaccurate
#Rnrm2 should be checked more
check_PROGRAMS = $(mpblas_binary128_test_PROGRAMS)
mpblas_binary128_smoke_TESTS = arithmetic.test_binary128 complex.test_binary128
mpblas_binary128_smoke_WRAPPERS = $(mpblas_binary128_smoke_TESTS:%=run_%.sh)
mpblas_binary128_test_WRAPPERS = $(mpblas_binary128_test_PROGRAMS:%=run_%.sh)
TESTS = $(mpblas_binary128_smoke_WRAPPERS)
TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = $(SHELL)
AM_TESTS_ENVIRONMENT = \
	LOG_COMPILER='$(LOG_COMPILER)'; export LOG_COMPILER; \
	EXEEXT='$(EXEEXT)'; export EXEEXT;
CLEANFILES = $(mpblas_binary128_smoke_WRAPPERS) $(mpblas_binary128_test_WRAPPERS)
.PHONY: check-all
check-all: $(check_PROGRAMS) $(mpblas_binary128_test_WRAPPERS)
	$(MAKE) $(AM_MAKEFLAGS) check TESTS="$(mpblas_binary128_test_WRAPPERS)"

run_%.sh: Makefile
	$(AM_V_GEN)printf '%s\n' '#!/bin/sh' > $@
	$(AM_V_at)printf '%s\n' '$${LOG_COMPILER} ./$*$${EXEEXT}' >> $@
	$(AM_V_at)chmod +x $@
mpblas_binary128_test_PROGRAMS = \
arithmetic.test_binary128 \
complex.test_binary128 \
iCamax.test.binary128_opt iCamax.test.binary128   \
iRamax.test.binary128_opt iRamax.test.binary128   \
Mlsame.test.binary128_opt Mlsame.test.binary128   \
Mxerbla.test.binary128_opt Mxerbla.test.binary128 \
CRrot.test.binary128_opt  CRrot.test.binary128    \
CRscal.test.binary128_opt CRscal.test.binary128   \
Caxpy.test.binary128_opt  Caxpy.test.binary128    \
Ccopy.test.binary128_opt  Ccopy.test.binary128    \
Cdotc.test.binary128_opt  Cdotc.test.binary128	\
Cdotu.test.binary128_opt  Cdotu.test.binary128	\
Cgbmv.test.binary128_opt  Cgbmv.test.binary128	\
Cgemm.test.binary128_opt  Cgemm.test.binary128	\
Cgemv.test.binary128_opt  Cgemv.test.binary128	\
Cgerc.test.binary128_opt  Cgerc.test.binary128	\
Cgeru.test.binary128_opt  Cgeru.test.binary128    \
Chbmv.test.binary128_opt  Chbmv.test.binary128	\
Chemm.test.binary128_opt  Chemm.test.binary128	\
Chemv.test.binary128_opt  Chemv.test.binary128	\
Cher.test.binary128_opt   Cher.test.binary128 	\
Cher2.test.binary128_opt  Cher2.test.binary128	\
Cher2k.test.binary128_opt Cher2k.test.binary128	\
Cherk.test.binary128_opt  Cherk.test.binary128	\
Chpmv.test.binary128_opt  Chpmv.test.binary128	\
Chpr.test.binary128_opt   Chpr.test.binary128	\
Chpr2.test.binary128_opt  Chpr2.test.binary128	\
Crotg.test.binary128_opt  Crotg.test.binary128	\
Cscal.test.binary128_opt  Cscal.test.binary128	\
Cswap.test.binary128_opt  Cswap.test.binary128	\
Csymm.test.binary128_opt  Csymm.test.binary128    \
Csyr2k.test.binary128_opt Csyr2k.test.binary128	\
Csyrk.test.binary128_opt  Csyrk.test.binary128	\
Ctbmv.test.binary128_opt  Ctbmv.test.binary128	\
Ctbsv.test.binary128_opt  Ctbsv.test.binary128	\
Ctpmv.test.binary128_opt  Ctpmv.test.binary128	\
Ctpsv.test.binary128_opt  Ctpsv.test.binary128	\
Ctrmm.test.binary128_opt  Ctrmm.test.binary128	\
Ctrmv.test.binary128_opt  Ctrmv.test.binary128	\
Ctrsm.test.binary128_opt  Ctrsm.test.binary128	\
Ctrsv.test.binary128_opt  Ctrsv.test.binary128	\
RCabs1.test.binary128_opt RCabs1.test.binary128	\
RCasum.test.binary128_opt RCasum.test.binary128	\
RCnrm2.test.binary128_opt RCnrm2.test.binary128	\
Rasum.test.binary128_opt  Rasum.test.binary128    \
Raxpy.test.binary128_opt  Raxpy.test.binary128	\
Rcopy.test.binary128_opt  Rcopy.test.binary128	\
Rdot.test.binary128_opt   Rdot.test.binary128	\
Rgbmv.test.binary128_opt  Rgbmv.test.binary128	\
Rgemm.test.binary128_opt  Rgemm.test.binary128	\
Rgemv.test.binary128_opt  Rgemv.test.binary128	\
Rger.test.binary128_opt   Rger.test.binary128	\
Rnrm2.test.binary128_opt  Rnrm2.test.binary128	\
Rrot.test.binary128_opt   Rrot.test.binary128	\
Rrotg.test.binary128_opt  Rrotg.test.binary128	\
Rrotm.test.binary128_opt  Rrotm.test.binary128	\
Rsbmv.test.binary128_opt  Rsbmv.test.binary128	\
Rscal.test.binary128_opt  Rscal.test.binary128	\
Rspmv.test.binary128_opt  Rspmv.test.binary128    \
Rspr.test.binary128_opt   Rspr.test.binary128	\
Rspr2.test.binary128_opt  Rspr2.test.binary128	\
Rswap.test.binary128_opt  Rswap.test.binary128	\
Rsymm.test.binary128_opt  Rsymm.test.binary128	\
Rsymv.test.binary128_opt  Rsymv.test.binary128	\
Rsyr.test.binary128_opt   Rsyr.test.binary128	\
Rsyr2.test.binary128_opt  Rsyr2.test.binary128	\
Rsyr2k.test.binary128_opt Rsyr2k.test.binary128	\
Rsyrk.test.binary128_opt  Rsyrk.test.binary128	\
Rtbmv.test.binary128_opt  Rtbmv.test.binary128	\
Rtbsv.test.binary128_opt  Rtbsv.test.binary128	\
Rtpmv.test.binary128_opt  Rtpmv.test.binary128	\
Rtpsv.test.binary128_opt  Rtpsv.test.binary128	\
Rtrmm.test.binary128_opt  Rtrmm.test.binary128    \
Rtrmv.test.binary128_opt  Rtrmv.test.binary128	\
Rtrsm.test.binary128_opt  Rtrsm.test.binary128	\
Rtrsv.test.binary128_opt  Rtrsv.test.binary128

mpblas_binary128_testdir = $(prefix)/lib/$(target)/mplapack/test/compare/binary128
mplibs=-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp
refblas=-L$(top_builddir)/external/i/LAPACK/lib -lblas $(FCLIBS)
test_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_BINARY128___ -D___MPLAPACK_INTERNAL___
test_srcdepends = ../common/xerbla.cpp ../common/mplapack.test.cpp ../common/Mxerbla.override.cpp

if IS_MACOS
mplibs+=-Wl,-flat_namespace,-undefined,dynamic_lookup
endif

if IS_MINGW
$(check_PROGRAMS): libmxerbla_override.la
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_BINARY128___
libmxerbla_override_la_SOURCES = ../common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS = -no-undefined
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,.libs/libmxerbla_override.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_binary128.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_mpfr.a,--no-whole-archive
opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,.libs/libmxerbla_override.a,$(top_builddir)/mpblas/optimized/binary128/.libs/libmpblas_binary128_opt.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_mpfr.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -lmplapack_mpfr
opt_libdepends = -L$(top_builddir)/mpblas/optimized/binary128 -lmpblas_binary128_opt -L$(top_builddir)/mplapack/reference -lmplapack_mpfr
endif

##
arithmetic_test_binary128_SOURCES=../common/arithmetic.test.cpp $(test_srcdepends)
arithmetic_test_binary128_CXXFLAGS=$(test_cxxflags)
arithmetic_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

complex_test_binary128_SOURCES=../common/complex.test.cpp $(test_srcdepends)
complex_test_binary128_CXXFLAGS=$(test_cxxflags)
complex_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

iCamax_test_binary128_opt_SOURCES=../common/iCamax.test.cpp $(test_srcdepends)
iCamax_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
iCamax_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iCamax_test_binary128_SOURCES=../common/iCamax.test.cpp $(test_srcdepends)
iCamax_test_binary128_CXXFLAGS=$(test_cxxflags)
iCamax_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

iRamax_test_binary128_opt_SOURCES=../common/iRamax.test.cpp $(test_srcdepends)
iRamax_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
iRamax_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iRamax_test_binary128_SOURCES=../common/iRamax.test.cpp $(test_srcdepends)
iRamax_test_binary128_CXXFLAGS=$(test_cxxflags)
iRamax_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Mlsame_test_binary128_opt_SOURCES=../common/Mlsame.test.cpp $(test_srcdepends)
Mlsame_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Mlsame_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mlsame_test_binary128_SOURCES=../common/Mlsame.test.cpp $(test_srcdepends)
Mlsame_test_binary128_CXXFLAGS=$(test_cxxflags)
Mlsame_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Mxerbla_test_binary128_opt_SOURCES=../common/Mxerbla.test.cpp $(test_srcdepends)
Mxerbla_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Mxerbla_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mxerbla_test_binary128_SOURCES=../common/Mxerbla.test.cpp $(test_srcdepends)
Mxerbla_test_binary128_CXXFLAGS=$(test_cxxflags)
Mxerbla_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

CRrot_test_binary128_opt_SOURCES=../common/CRrot.test.cpp $(test_srcdepends)
CRrot_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
CRrot_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRrot_test_binary128_SOURCES=../common/CRrot.test.cpp $(test_srcdepends)
CRrot_test_binary128_CXXFLAGS=$(test_cxxflags)
CRrot_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

CRscal_test_binary128_opt_SOURCES=../common/CRscal.test.cpp $(test_srcdepends)
CRscal_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
CRscal_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRscal_test_binary128_SOURCES=../common/CRscal.test.cpp $(test_srcdepends)
CRscal_test_binary128_CXXFLAGS=$(test_cxxflags)
CRscal_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Caxpy_test_binary128_SOURCES=../common/Caxpy.test.cpp $(test_srcdepends)
Caxpy_test_binary128_CXXFLAGS=$(test_cxxflags)
Caxpy_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Caxpy_test_binary128_opt_SOURCES=../common/Caxpy.test.cpp $(test_srcdepends)
Caxpy_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Caxpy_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ccopy_test_binary128_SOURCES=../common/Ccopy.test.cpp $(test_srcdepends)
Ccopy_test_binary128_CXXFLAGS=$(test_cxxflags)
Ccopy_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ccopy_test_binary128_opt_SOURCES=../common/Ccopy.test.cpp $(test_srcdepends)
Ccopy_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Ccopy_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotc_test_binary128_SOURCES=../common/Cdotc.test.cpp $(test_srcdepends)
Cdotc_test_binary128_CXXFLAGS=$(test_cxxflags)
Cdotc_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cdotc_test_binary128_opt_SOURCES=../common/Cdotc.test.cpp $(test_srcdepends)
Cdotc_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Cdotc_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotu_test_binary128_SOURCES=../common/Cdotu.test.cpp $(test_srcdepends)
Cdotu_test_binary128_CXXFLAGS=$(test_cxxflags)
Cdotu_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cdotu_test_binary128_opt_SOURCES=../common/Cdotu.test.cpp $(test_srcdepends)
Cdotu_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Cdotu_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgbmv_test_binary128_SOURCES=../common/Cgbmv.test.cpp $(test_srcdepends)
Cgbmv_test_binary128_CXXFLAGS=$(test_cxxflags)
Cgbmv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgbmv_test_binary128_opt_SOURCES=../common/Cgbmv.test.cpp $(test_srcdepends)
Cgbmv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Cgbmv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemm_test_binary128_SOURCES=../common/Cgemm.test.cpp $(test_srcdepends)
Cgemm_test_binary128_CXXFLAGS=$(test_cxxflags)
Cgemm_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgemm_test_binary128_opt_SOURCES=../common/Cgemm.test.cpp $(test_srcdepends)
Cgemm_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Cgemm_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemv_test_binary128_SOURCES=../common/Cgemv.test.cpp $(test_srcdepends)
Cgemv_test_binary128_CXXFLAGS=$(test_cxxflags)
Cgemv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgemv_test_binary128_opt_SOURCES=../common/Cgemv.test.cpp $(test_srcdepends)
Cgemv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Cgemv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgerc_test_binary128_SOURCES=../common/Cgerc.test.cpp $(test_srcdepends)
Cgerc_test_binary128_CXXFLAGS=$(test_cxxflags)
Cgerc_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgerc_test_binary128_opt_SOURCES=../common/Cgerc.test.cpp $(test_srcdepends)
Cgerc_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Cgerc_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgeru_test_binary128_SOURCES=../common/Cgeru.test.cpp $(test_srcdepends)
Cgeru_test_binary128_CXXFLAGS=$(test_cxxflags)
Cgeru_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgeru_test_binary128_opt_SOURCES=../common/Cgeru.test.cpp $(test_srcdepends)
Cgeru_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Cgeru_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chbmv_test_binary128_SOURCES=../common/Chbmv.test.cpp $(test_srcdepends)
Chbmv_test_binary128_CXXFLAGS=$(test_cxxflags)
Chbmv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chbmv_test_binary128_opt_SOURCES=../common/Chbmv.test.cpp $(test_srcdepends)
Chbmv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Chbmv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemm_test_binary128_SOURCES=../common/Chemm.test.cpp $(test_srcdepends)
Chemm_test_binary128_CXXFLAGS=$(test_cxxflags)
Chemm_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chemm_test_binary128_opt_SOURCES=../common/Chemm.test.cpp $(test_srcdepends)
Chemm_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Chemm_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemv_test_binary128_SOURCES=../common/Chemv.test.cpp $(test_srcdepends)
Chemv_test_binary128_CXXFLAGS=$(test_cxxflags)
Chemv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chemv_test_binary128_opt_SOURCES=../common/Chemv.test.cpp $(test_srcdepends)
Chemv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Chemv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher_test_binary128_SOURCES=../common/Cher.test.cpp $(test_srcdepends)
Cher_test_binary128_CXXFLAGS=$(test_cxxflags)
Cher_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher_test_binary128_opt_SOURCES=../common/Cher.test.cpp $(test_srcdepends)
Cher_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Cher_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2_test_binary128_SOURCES=../common/Cher2.test.cpp $(test_srcdepends)
Cher2_test_binary128_CXXFLAGS=$(test_cxxflags)
Cher2_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher2_test_binary128_opt_SOURCES=../common/Cher2.test.cpp $(test_srcdepends)
Cher2_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Cher2_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2k_test_binary128_SOURCES=../common/Cher2k.test.cpp $(test_srcdepends)
Cher2k_test_binary128_CXXFLAGS=$(test_cxxflags)
Cher2k_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher2k_test_binary128_opt_SOURCES=../common/Cher2k.test.cpp $(test_srcdepends)
Cher2k_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Cher2k_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cherk_test_binary128_SOURCES=../common/Cherk.test.cpp $(test_srcdepends)
Cherk_test_binary128_CXXFLAGS=$(test_cxxflags)
Cherk_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cherk_test_binary128_opt_SOURCES=../common/Cherk.test.cpp $(test_srcdepends)
Cherk_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Cherk_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpmv_test_binary128_SOURCES=../common/Chpmv.test.cpp $(test_srcdepends)
Chpmv_test_binary128_CXXFLAGS=$(test_cxxflags)
Chpmv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpmv_test_binary128_opt_SOURCES=../common/Chpmv.test.cpp $(test_srcdepends)
Chpmv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Chpmv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr_test_binary128_SOURCES=../common/Chpr.test.cpp $(test_srcdepends)
Chpr_test_binary128_CXXFLAGS=$(test_cxxflags)
Chpr_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpr_test_binary128_opt_SOURCES=../common/Chpr.test.cpp $(test_srcdepends)
Chpr_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Chpr_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr2_test_binary128_SOURCES=../common/Chpr2.test.cpp $(test_srcdepends)
Chpr2_test_binary128_CXXFLAGS=$(test_cxxflags)
Chpr2_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpr2_test_binary128_opt_SOURCES=../common/Chpr2.test.cpp $(test_srcdepends)
Chpr2_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Chpr2_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Crotg_test_binary128_SOURCES=../common/Crotg.test.cpp $(test_srcdepends)
Crotg_test_binary128_CXXFLAGS=$(test_cxxflags)
Crotg_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Crotg_test_binary128_opt_SOURCES=../common/Crotg.test.cpp $(test_srcdepends)
Crotg_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Crotg_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cscal_test_binary128_SOURCES=../common/Cscal.test.cpp $(test_srcdepends)
Cscal_test_binary128_CXXFLAGS=$(test_cxxflags)
Cscal_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cscal_test_binary128_opt_SOURCES=../common/Cscal.test.cpp $(test_srcdepends)
Cscal_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Cscal_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cswap_test_binary128_SOURCES=../common/Cswap.test.cpp $(test_srcdepends)
Cswap_test_binary128_CXXFLAGS=$(test_cxxflags)
Cswap_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cswap_test_binary128_opt_SOURCES=../common/Cswap.test.cpp $(test_srcdepends)
Cswap_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Cswap_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csymm_test_binary128_SOURCES=../common/Csymm.test.cpp $(test_srcdepends)
Csymm_test_binary128_CXXFLAGS=$(test_cxxflags)
Csymm_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csymm_test_binary128_opt_SOURCES=../common/Csymm.test.cpp $(test_srcdepends)
Csymm_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Csymm_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyr2k_test_binary128_SOURCES=../common/Csyr2k.test.cpp $(test_srcdepends)
Csyr2k_test_binary128_CXXFLAGS=$(test_cxxflags)
Csyr2k_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csyr2k_test_binary128_opt_SOURCES=../common/Csyr2k.test.cpp $(test_srcdepends)
Csyr2k_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Csyr2k_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyrk_test_binary128_SOURCES=../common/Csyrk.test.cpp $(test_srcdepends)
Csyrk_test_binary128_CXXFLAGS=$(test_cxxflags)
Csyrk_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csyrk_test_binary128_opt_SOURCES=../common/Csyrk.test.cpp $(test_srcdepends)
Csyrk_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Csyrk_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbmv_test_binary128_SOURCES=../common/Ctbmv.test.cpp $(test_srcdepends)
Ctbmv_test_binary128_CXXFLAGS=$(test_cxxflags)
Ctbmv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctbmv_test_binary128_opt_SOURCES=../common/Ctbmv.test.cpp $(test_srcdepends)
Ctbmv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Ctbmv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbsv_test_binary128_SOURCES=../common/Ctbsv.test.cpp $(test_srcdepends)
Ctbsv_test_binary128_CXXFLAGS=$(test_cxxflags)
Ctbsv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctbsv_test_binary128_opt_SOURCES=../common/Ctbsv.test.cpp $(test_srcdepends)
Ctbsv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Ctbsv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpmv_test_binary128_SOURCES=../common/Ctpmv.test.cpp $(test_srcdepends)
Ctpmv_test_binary128_CXXFLAGS=$(test_cxxflags)
Ctpmv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctpmv_test_binary128_opt_SOURCES=../common/Ctpmv.test.cpp $(test_srcdepends)
Ctpmv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Ctpmv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpsv_test_binary128_SOURCES=../common/Ctpsv.test.cpp $(test_srcdepends)
Ctpsv_test_binary128_CXXFLAGS=$(test_cxxflags)
Ctpsv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctpsv_test_binary128_opt_SOURCES=../common/Ctpsv.test.cpp $(test_srcdepends)
Ctpsv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Ctpsv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmm_test_binary128_SOURCES=../common/Ctrmm.test.cpp $(test_srcdepends)
Ctrmm_test_binary128_CXXFLAGS=$(test_cxxflags)
Ctrmm_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrmm_test_binary128_opt_SOURCES=../common/Ctrmm.test.cpp $(test_srcdepends)
Ctrmm_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Ctrmm_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmv_test_binary128_SOURCES=../common/Ctrmv.test.cpp $(test_srcdepends)
Ctrmv_test_binary128_CXXFLAGS=$(test_cxxflags)
Ctrmv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrmv_test_binary128_opt_SOURCES=../common/Ctrmv.test.cpp $(test_srcdepends)
Ctrmv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Ctrmv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsm_test_binary128_SOURCES=../common/Ctrsm.test.cpp $(test_srcdepends)
Ctrsm_test_binary128_CXXFLAGS=$(test_cxxflags)
Ctrsm_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrsm_test_binary128_opt_SOURCES=../common/Ctrsm.test.cpp $(test_srcdepends)
Ctrsm_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Ctrsm_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsv_test_binary128_SOURCES=../common/Ctrsv.test.cpp $(test_srcdepends)
Ctrsv_test_binary128_CXXFLAGS=$(test_cxxflags)
Ctrsv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrsv_test_binary128_opt_SOURCES=../common/Ctrsv.test.cpp $(test_srcdepends)
Ctrsv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Ctrsv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCabs1_test_binary128_opt_SOURCES=../common/RCabs1.test.cpp $(test_srcdepends)
RCabs1_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
RCabs1_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCabs1_test_binary128_SOURCES=../common/RCabs1.test.cpp $(test_srcdepends)
RCabs1_test_binary128_CXXFLAGS=$(test_cxxflags)
RCabs1_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

RCasum_test_binary128_opt_SOURCES=../common/RCasum.test.cpp $(test_srcdepends)
RCasum_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
RCasum_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCasum_test_binary128_SOURCES=../common/RCasum.test.cpp $(test_srcdepends)
RCasum_test_binary128_CXXFLAGS=$(test_cxxflags)
RCasum_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

RCnrm2_test_binary128_opt_SOURCES=../common/RCnrm2.test.cpp $(test_srcdepends)
RCnrm2_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
RCnrm2_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCnrm2_test_binary128_SOURCES=../common/RCnrm2.test.cpp $(test_srcdepends)
RCnrm2_test_binary128_CXXFLAGS=$(test_cxxflags)
RCnrm2_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rasum_test_binary128_opt_SOURCES=../common/Rasum.test.cpp $(test_srcdepends)
Rasum_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rasum_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rasum_test_binary128_SOURCES=../common/Rasum.test.cpp $(test_srcdepends)
Rasum_test_binary128_CXXFLAGS=$(test_cxxflags)
Rasum_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Raxpy_test_binary128_opt_SOURCES=../common/Raxpy.test.cpp $(test_srcdepends)
Raxpy_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Raxpy_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Raxpy_test_binary128_SOURCES=../common/Raxpy.test.cpp $(test_srcdepends)
Raxpy_test_binary128_CXXFLAGS=$(test_cxxflags)
Raxpy_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rcopy_test_binary128_opt_SOURCES=../common/Rcopy.test.cpp $(test_srcdepends)
Rcopy_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rcopy_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rcopy_test_binary128_SOURCES=../common/Rcopy.test.cpp $(test_srcdepends)
Rcopy_test_binary128_CXXFLAGS=$(test_cxxflags)
Rcopy_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rdot_test_binary128_opt_SOURCES=../common/Rdot.test.cpp $(test_srcdepends)
Rdot_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rdot_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rdot_test_binary128_SOURCES=../common/Rdot.test.cpp $(test_srcdepends)
Rdot_test_binary128_CXXFLAGS=$(test_cxxflags)
Rdot_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgbmv_test_binary128_opt_SOURCES=../common/Rgbmv.test.cpp $(test_srcdepends)
Rgbmv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rgbmv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgbmv_test_binary128_SOURCES=../common/Rgbmv.test.cpp $(test_srcdepends)
Rgbmv_test_binary128_CXXFLAGS=$(test_cxxflags)
Rgbmv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgemm_test_binary128_opt_SOURCES=../common/Rgemm.test.cpp $(test_srcdepends)
Rgemm_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rgemm_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemm_test_binary128_SOURCES=../common/Rgemm.test.cpp $(test_srcdepends)
Rgemm_test_binary128_CXXFLAGS=$(test_cxxflags)
Rgemm_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgemv_test_binary128_opt_SOURCES=../common/Rgemv.test.cpp $(test_srcdepends)
Rgemv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rgemv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemv_test_binary128_SOURCES=../common/Rgemv.test.cpp $(test_srcdepends)
Rgemv_test_binary128_CXXFLAGS=$(test_cxxflags)
Rgemv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rger_test_binary128_opt_SOURCES=../common/Rger.test.cpp $(test_srcdepends)
Rger_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rger_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rger_test_binary128_SOURCES=../common/Rger.test.cpp $(test_srcdepends)
Rger_test_binary128_CXXFLAGS=$(test_cxxflags)
Rger_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rnrm2_test_binary128_opt_SOURCES=../common/Rnrm2.test.cpp $(test_srcdepends)
Rnrm2_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rnrm2_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rnrm2_test_binary128_SOURCES=../common/Rnrm2.test.cpp $(test_srcdepends)
Rnrm2_test_binary128_CXXFLAGS=$(test_cxxflags)
Rnrm2_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrot_test_binary128_opt_SOURCES=../common/Rrot.test.cpp $(test_srcdepends)
Rrot_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rrot_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrot_test_binary128_SOURCES=../common/Rrot.test.cpp $(test_srcdepends)
Rrot_test_binary128_CXXFLAGS=$(test_cxxflags)
Rrot_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrotg_test_binary128_opt_SOURCES=../common/Rrotg.test.cpp $(test_srcdepends)
Rrotg_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rrotg_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotg_test_binary128_SOURCES=../common/Rrotg.test.cpp $(test_srcdepends)
Rrotg_test_binary128_CXXFLAGS=$(test_cxxflags)
Rrotg_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrotm_test_binary128_opt_SOURCES=../common/Rrotm.test.cpp $(test_srcdepends)
Rrotm_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rrotm_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotm_test_binary128_SOURCES=../common/Rrotm.test.cpp $(test_srcdepends)
Rrotm_test_binary128_CXXFLAGS=$(test_cxxflags)
Rrotm_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsbmv_test_binary128_opt_SOURCES=../common/Rsbmv.test.cpp $(test_srcdepends)
Rsbmv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rsbmv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsbmv_test_binary128_SOURCES=../common/Rsbmv.test.cpp $(test_srcdepends)
Rsbmv_test_binary128_CXXFLAGS=$(test_cxxflags)
Rsbmv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rscal_test_binary128_opt_SOURCES=../common/Rscal.test.cpp $(test_srcdepends)
Rscal_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rscal_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rscal_test_binary128_SOURCES=../common/Rscal.test.cpp $(test_srcdepends)
Rscal_test_binary128_CXXFLAGS=$(test_cxxflags)
Rscal_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspmv_test_binary128_opt_SOURCES=../common/Rspmv.test.cpp $(test_srcdepends)
Rspmv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rspmv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspmv_test_binary128_SOURCES=../common/Rspmv.test.cpp $(test_srcdepends)
Rspmv_test_binary128_CXXFLAGS=$(test_cxxflags)
Rspmv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspr_test_binary128_opt_SOURCES=../common/Rspr.test.cpp $(test_srcdepends)
Rspr_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rspr_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr_test_binary128_SOURCES=../common/Rspr.test.cpp $(test_srcdepends)
Rspr_test_binary128_CXXFLAGS=$(test_cxxflags)
Rspr_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspr2_test_binary128_opt_SOURCES=../common/Rspr2.test.cpp $(test_srcdepends)
Rspr2_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rspr2_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr2_test_binary128_SOURCES=../common/Rspr2.test.cpp $(test_srcdepends)
Rspr2_test_binary128_CXXFLAGS=$(test_cxxflags)
Rspr2_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rswap_test_binary128_opt_SOURCES=../common/Rswap.test.cpp $(test_srcdepends)
Rswap_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rswap_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rswap_test_binary128_SOURCES=../common/Rswap.test.cpp $(test_srcdepends)
Rswap_test_binary128_CXXFLAGS=$(test_cxxflags)
Rswap_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsymm_test_binary128_opt_SOURCES=../common/Rsymm.test.cpp $(test_srcdepends)
Rsymm_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rsymm_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymm_test_binary128_SOURCES=../common/Rsymm.test.cpp $(test_srcdepends)
Rsymm_test_binary128_CXXFLAGS=$(test_cxxflags)
Rsymm_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsymv_test_binary128_opt_SOURCES=../common/Rsymv.test.cpp $(test_srcdepends)
Rsymv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rsymv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymv_test_binary128_SOURCES=../common/Rsymv.test.cpp $(test_srcdepends)
Rsymv_test_binary128_CXXFLAGS=$(test_cxxflags)
Rsymv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr_test_binary128_opt_SOURCES=../common/Rsyr.test.cpp $(test_srcdepends)
Rsyr_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rsyr_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr_test_binary128_SOURCES=../common/Rsyr.test.cpp $(test_srcdepends)
Rsyr_test_binary128_CXXFLAGS=$(test_cxxflags)
Rsyr_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr2_test_binary128_opt_SOURCES=../common/Rsyr2.test.cpp $(test_srcdepends)
Rsyr2_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rsyr2_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2_test_binary128_SOURCES=../common/Rsyr2.test.cpp $(test_srcdepends)
Rsyr2_test_binary128_CXXFLAGS=$(test_cxxflags)
Rsyr2_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr2k_test_binary128_opt_SOURCES=../common/Rsyr2k.test.cpp $(test_srcdepends)
Rsyr2k_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rsyr2k_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2k_test_binary128_SOURCES=../common/Rsyr2k.test.cpp $(test_srcdepends)
Rsyr2k_test_binary128_CXXFLAGS=$(test_cxxflags)
Rsyr2k_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyrk_test_binary128_opt_SOURCES=../common/Rsyrk.test.cpp $(test_srcdepends)
Rsyrk_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rsyrk_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyrk_test_binary128_SOURCES=../common/Rsyrk.test.cpp $(test_srcdepends)
Rsyrk_test_binary128_CXXFLAGS=$(test_cxxflags)
Rsyrk_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtbmv_test_binary128_opt_SOURCES=../common/Rtbmv.test.cpp $(test_srcdepends)
Rtbmv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rtbmv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbmv_test_binary128_SOURCES=../common/Rtbmv.test.cpp $(test_srcdepends)
Rtbmv_test_binary128_CXXFLAGS=$(test_cxxflags)
Rtbmv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtbsv_test_binary128_opt_SOURCES=../common/Rtbsv.test.cpp $(test_srcdepends)
Rtbsv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rtbsv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbsv_test_binary128_SOURCES=../common/Rtbsv.test.cpp $(test_srcdepends)
Rtbsv_test_binary128_CXXFLAGS=$(test_cxxflags)
Rtbsv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtpmv_test_binary128_opt_SOURCES=../common/Rtpmv.test.cpp $(test_srcdepends)
Rtpmv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rtpmv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpmv_test_binary128_SOURCES=../common/Rtpmv.test.cpp $(test_srcdepends)
Rtpmv_test_binary128_CXXFLAGS=$(test_cxxflags)
Rtpmv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtpsv_test_binary128_opt_SOURCES=../common/Rtpsv.test.cpp $(test_srcdepends)
Rtpsv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rtpsv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpsv_test_binary128_SOURCES=../common/Rtpsv.test.cpp $(test_srcdepends)
Rtpsv_test_binary128_CXXFLAGS=$(test_cxxflags)
Rtpsv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrmm_test_binary128_opt_SOURCES=../common/Rtrmm.test.cpp $(test_srcdepends)
Rtrmm_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rtrmm_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmm_test_binary128_SOURCES=../common/Rtrmm.test.cpp $(test_srcdepends)
Rtrmm_test_binary128_CXXFLAGS=$(test_cxxflags)
Rtrmm_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrmv_test_binary128_opt_SOURCES=../common/Rtrmv.test.cpp $(test_srcdepends)
Rtrmv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rtrmv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmv_test_binary128_SOURCES=../common/Rtrmv.test.cpp $(test_srcdepends)
Rtrmv_test_binary128_CXXFLAGS=$(test_cxxflags)
Rtrmv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrsm_test_binary128_opt_SOURCES=../common/Rtrsm.test.cpp $(test_srcdepends)
Rtrsm_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rtrsm_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsm_test_binary128_SOURCES=../common/Rtrsm.test.cpp $(test_srcdepends)
Rtrsm_test_binary128_CXXFLAGS=$(test_cxxflags)
Rtrsm_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrsv_test_binary128_opt_SOURCES=../common/Rtrsv.test.cpp $(test_srcdepends)
Rtrsv_test_binary128_opt_CXXFLAGS=$(test_cxxflags)
Rtrsv_test_binary128_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsv_test_binary128_SOURCES=../common/Rtrsv.test.cpp $(test_srcdepends)
Rtrsv_test_binary128_CXXFLAGS=$(test_cxxflags)
Rtrsv_test_binary128_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mpblas_binary128_testdir) $(prefix)
endif
