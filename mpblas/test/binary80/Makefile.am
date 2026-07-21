#CTBSV, CTPSV might have bugs or intrinsically inaccurate
#Rnrm2 should be checked more
check_PROGRAMS = $(mpblas_binary80_test_PROGRAMS) 
mpblas_binary80_smoke_TESTS = arithmetic.test_binary80 complex.test_binary80
mpblas_binary80_smoke_WRAPPERS = $(mpblas_binary80_smoke_TESTS:%=run_%.sh)
mpblas_binary80_test_WRAPPERS = $(mpblas_binary80_test_PROGRAMS:%=run_%.sh)
TESTS = $(mpblas_binary80_smoke_WRAPPERS)
TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = $(SHELL)
AM_TESTS_ENVIRONMENT = \
	LOG_COMPILER='$(LOG_COMPILER)'; export LOG_COMPILER; \
	EXEEXT='$(EXEEXT)'; export EXEEXT;
CLEANFILES = $(mpblas_binary80_smoke_WRAPPERS) $(mpblas_binary80_test_WRAPPERS)
.PHONY: check-all
check-all: $(check_PROGRAMS) $(mpblas_binary80_test_WRAPPERS)
	$(MAKE) $(AM_MAKEFLAGS) check TESTS="$(mpblas_binary80_test_WRAPPERS)"

run_%.sh: Makefile
	$(AM_V_GEN)printf '%s\n' '#!/bin/sh' > $@
	$(AM_V_at)printf '%s\n' '$${LOG_COMPILER} ./$*$${EXEEXT}' >> $@
	$(AM_V_at)chmod +x $@
mpblas_binary80_test_PROGRAMS = \
arithmetic.test_binary80 \
complex.test_binary80 \
iCamax.test.binary80_opt iCamax.test.binary80   \
iRamax.test.binary80_opt iRamax.test.binary80   \
Mlsame.test.binary80_opt Mlsame.test.binary80   \
Mxerbla.test.binary80_opt Mxerbla.test.binary80 \
CRrot.test.binary80_opt  CRrot.test.binary80    \
CRscal.test.binary80_opt CRscal.test.binary80   \
Caxpy.test.binary80_opt  Caxpy.test.binary80    \
Ccopy.test.binary80_opt  Ccopy.test.binary80    \
Cdotc.test.binary80_opt  Cdotc.test.binary80	\
Cdotu.test.binary80_opt  Cdotu.test.binary80	\
Cgbmv.test.binary80_opt  Cgbmv.test.binary80	\
Cgemm.test.binary80_opt  Cgemm.test.binary80	\
Cgemv.test.binary80_opt  Cgemv.test.binary80	\
Cgerc.test.binary80_opt  Cgerc.test.binary80	\
Cgeru.test.binary80_opt  Cgeru.test.binary80    \
Chbmv.test.binary80_opt  Chbmv.test.binary80	\
Chemm.test.binary80_opt  Chemm.test.binary80	\
Chemv.test.binary80_opt  Chemv.test.binary80	\
Cher.test.binary80_opt   Cher.test.binary80 	\
Cher2.test.binary80_opt  Cher2.test.binary80	\
Cher2k.test.binary80_opt Cher2k.test.binary80	\
Cherk.test.binary80_opt  Cherk.test.binary80	\
Chpmv.test.binary80_opt  Chpmv.test.binary80	\
Chpr.test.binary80_opt   Chpr.test.binary80	\
Chpr2.test.binary80_opt  Chpr2.test.binary80	\
Crotg.test.binary80_opt  Crotg.test.binary80	\
Cscal.test.binary80_opt  Cscal.test.binary80	\
Cswap.test.binary80_opt  Cswap.test.binary80	\
Csymm.test.binary80_opt  Csymm.test.binary80    \
Csyr2k.test.binary80_opt Csyr2k.test.binary80	\
Csyrk.test.binary80_opt  Csyrk.test.binary80	\
Ctbmv.test.binary80_opt  Ctbmv.test.binary80	\
Ctbsv.test.binary80_opt  Ctbsv.test.binary80	\
Ctpmv.test.binary80_opt  Ctpmv.test.binary80	\
Ctpsv.test.binary80_opt  Ctpsv.test.binary80	\
Ctrmm.test.binary80_opt  Ctrmm.test.binary80	\
Ctrmv.test.binary80_opt  Ctrmv.test.binary80	\
Ctrsm.test.binary80_opt  Ctrsm.test.binary80	\
Ctrsv.test.binary80_opt  Ctrsv.test.binary80	\
RCabs1.test.binary80_opt RCabs1.test.binary80	\
RCasum.test.binary80_opt RCasum.test.binary80	\
RCnrm2.test.binary80_opt RCnrm2.test.binary80	\
Rasum.test.binary80_opt  Rasum.test.binary80    \
Raxpy.test.binary80_opt  Raxpy.test.binary80	\
Rcopy.test.binary80_opt  Rcopy.test.binary80	\
Rdot.test.binary80_opt   Rdot.test.binary80	\
Rgbmv.test.binary80_opt  Rgbmv.test.binary80	\
Rgemm.test.binary80_opt  Rgemm.test.binary80	\
Rgemv.test.binary80_opt  Rgemv.test.binary80	\
Rger.test.binary80_opt   Rger.test.binary80	\
Rnrm2.test.binary80_opt  Rnrm2.test.binary80	\
Rrot.test.binary80_opt   Rrot.test.binary80	\
Rrotg.test.binary80_opt  Rrotg.test.binary80	\
Rrotm.test.binary80_opt  Rrotm.test.binary80	\
Rsbmv.test.binary80_opt  Rsbmv.test.binary80	\
Rscal.test.binary80_opt  Rscal.test.binary80	\
Rspmv.test.binary80_opt  Rspmv.test.binary80    \
Rspr.test.binary80_opt   Rspr.test.binary80	\
Rspr2.test.binary80_opt  Rspr2.test.binary80	\
Rswap.test.binary80_opt  Rswap.test.binary80	\
Rsymm.test.binary80_opt  Rsymm.test.binary80	\
Rsymv.test.binary80_opt  Rsymv.test.binary80	\
Rsyr.test.binary80_opt   Rsyr.test.binary80	\
Rsyr2.test.binary80_opt  Rsyr2.test.binary80	\
Rsyr2k.test.binary80_opt Rsyr2k.test.binary80	\
Rsyrk.test.binary80_opt  Rsyrk.test.binary80	\
Rtbmv.test.binary80_opt  Rtbmv.test.binary80	\
Rtbsv.test.binary80_opt  Rtbsv.test.binary80	\
Rtpmv.test.binary80_opt  Rtpmv.test.binary80	\
Rtpsv.test.binary80_opt  Rtpsv.test.binary80	\
Rtrmm.test.binary80_opt  Rtrmm.test.binary80    \
Rtrmv.test.binary80_opt  Rtrmv.test.binary80	\
Rtrsm.test.binary80_opt  Rtrsm.test.binary80	\
Rtrsv.test.binary80_opt  Rtrsv.test.binary80

mpblas_binary80_testdir = $(prefix)/lib/$(target)/mplapack/test/compare/binary80
mplibs=-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp
refblas=-L$(top_builddir)/external/i/LAPACK/lib -lblas $(FCLIBS)
test_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_BINARY80___ -D___MPLAPACK_INTERNAL___
test_srcdepends = ../common/xerbla.cpp ../common/mplapack.test.cpp ../common/Mxerbla.override.cpp

if IS_MACOS
mplibs+=-Wl,-flat_namespace,-undefined,dynamic_lookup
endif

if IS_MINGW
$(check_PROGRAMS): libmxerbla_override.la
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_BINARY80___
libmxerbla_override_la_SOURCES = ../common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS = -no-undefined
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,.libs/libmxerbla_override.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_binary80.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_mpfr.a,--no-whole-archive
opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,.libs/libmxerbla_override.a,$(top_builddir)/mpblas/optimized/binary80/.libs/libmpblas_binary80_opt.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_mpfr.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_binary80 -lmplapack_mpfr
opt_libdepends = -L$(top_builddir)/mpblas/optimized/binary80 -lmpblas_binary80_opt -L$(top_builddir)/mplapack/reference -lmplapack_mpfr
endif

##
arithmetic_test_binary80_SOURCES=../common/arithmetic.test.cpp $(test_srcdepends)
arithmetic_test_binary80_CXXFLAGS=$(test_cxxflags)
arithmetic_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

complex_test_binary80_SOURCES=../common/complex.test.cpp $(test_srcdepends)
complex_test_binary80_CXXFLAGS=$(test_cxxflags)
complex_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

iCamax_test_binary80_opt_SOURCES=../common/iCamax.test.cpp $(test_srcdepends)
iCamax_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
iCamax_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iCamax_test_binary80_SOURCES=../common/iCamax.test.cpp $(test_srcdepends)
iCamax_test_binary80_CXXFLAGS=$(test_cxxflags)
iCamax_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

iRamax_test_binary80_opt_SOURCES=../common/iRamax.test.cpp $(test_srcdepends)
iRamax_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
iRamax_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

iRamax_test_binary80_SOURCES=../common/iRamax.test.cpp $(test_srcdepends)
iRamax_test_binary80_CXXFLAGS=$(test_cxxflags)
iRamax_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Mlsame_test_binary80_opt_SOURCES=../common/Mlsame.test.cpp $(test_srcdepends)
Mlsame_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Mlsame_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mlsame_test_binary80_SOURCES=../common/Mlsame.test.cpp $(test_srcdepends)
Mlsame_test_binary80_CXXFLAGS=$(test_cxxflags)
Mlsame_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Mxerbla_test_binary80_opt_SOURCES=../common/Mxerbla.test.cpp $(test_srcdepends)
Mxerbla_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Mxerbla_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Mxerbla_test_binary80_SOURCES=../common/Mxerbla.test.cpp $(test_srcdepends)
Mxerbla_test_binary80_CXXFLAGS=$(test_cxxflags)
Mxerbla_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

CRrot_test_binary80_opt_SOURCES=../common/CRrot.test.cpp $(test_srcdepends)
CRrot_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
CRrot_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRrot_test_binary80_SOURCES=../common/CRrot.test.cpp $(test_srcdepends)
CRrot_test_binary80_CXXFLAGS=$(test_cxxflags)
CRrot_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

CRscal_test_binary80_opt_SOURCES=../common/CRscal.test.cpp $(test_srcdepends)
CRscal_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
CRscal_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

CRscal_test_binary80_SOURCES=../common/CRscal.test.cpp $(test_srcdepends)
CRscal_test_binary80_CXXFLAGS=$(test_cxxflags)
CRscal_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Caxpy_test_binary80_SOURCES=../common/Caxpy.test.cpp $(test_srcdepends)
Caxpy_test_binary80_CXXFLAGS=$(test_cxxflags)
Caxpy_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Caxpy_test_binary80_opt_SOURCES=../common/Caxpy.test.cpp $(test_srcdepends)
Caxpy_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Caxpy_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ccopy_test_binary80_SOURCES=../common/Ccopy.test.cpp $(test_srcdepends)
Ccopy_test_binary80_CXXFLAGS=$(test_cxxflags)
Ccopy_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ccopy_test_binary80_opt_SOURCES=../common/Ccopy.test.cpp $(test_srcdepends)
Ccopy_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Ccopy_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotc_test_binary80_SOURCES=../common/Cdotc.test.cpp $(test_srcdepends)
Cdotc_test_binary80_CXXFLAGS=$(test_cxxflags)
Cdotc_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cdotc_test_binary80_opt_SOURCES=../common/Cdotc.test.cpp $(test_srcdepends)
Cdotc_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Cdotc_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cdotu_test_binary80_SOURCES=../common/Cdotu.test.cpp $(test_srcdepends)
Cdotu_test_binary80_CXXFLAGS=$(test_cxxflags)
Cdotu_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cdotu_test_binary80_opt_SOURCES=../common/Cdotu.test.cpp $(test_srcdepends)
Cdotu_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Cdotu_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgbmv_test_binary80_SOURCES=../common/Cgbmv.test.cpp $(test_srcdepends)
Cgbmv_test_binary80_CXXFLAGS=$(test_cxxflags)
Cgbmv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgbmv_test_binary80_opt_SOURCES=../common/Cgbmv.test.cpp $(test_srcdepends)
Cgbmv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Cgbmv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemm_test_binary80_SOURCES=../common/Cgemm.test.cpp $(test_srcdepends)
Cgemm_test_binary80_CXXFLAGS=$(test_cxxflags)
Cgemm_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgemm_test_binary80_opt_SOURCES=../common/Cgemm.test.cpp $(test_srcdepends)
Cgemm_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Cgemm_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgemv_test_binary80_SOURCES=../common/Cgemv.test.cpp $(test_srcdepends)
Cgemv_test_binary80_CXXFLAGS=$(test_cxxflags)
Cgemv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgemv_test_binary80_opt_SOURCES=../common/Cgemv.test.cpp $(test_srcdepends)
Cgemv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Cgemv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgerc_test_binary80_SOURCES=../common/Cgerc.test.cpp $(test_srcdepends)
Cgerc_test_binary80_CXXFLAGS=$(test_cxxflags)
Cgerc_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgerc_test_binary80_opt_SOURCES=../common/Cgerc.test.cpp $(test_srcdepends)
Cgerc_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Cgerc_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cgeru_test_binary80_SOURCES=../common/Cgeru.test.cpp $(test_srcdepends)
Cgeru_test_binary80_CXXFLAGS=$(test_cxxflags)
Cgeru_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cgeru_test_binary80_opt_SOURCES=../common/Cgeru.test.cpp $(test_srcdepends)
Cgeru_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Cgeru_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chbmv_test_binary80_SOURCES=../common/Chbmv.test.cpp $(test_srcdepends)
Chbmv_test_binary80_CXXFLAGS=$(test_cxxflags)
Chbmv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chbmv_test_binary80_opt_SOURCES=../common/Chbmv.test.cpp $(test_srcdepends)
Chbmv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Chbmv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemm_test_binary80_SOURCES=../common/Chemm.test.cpp $(test_srcdepends)
Chemm_test_binary80_CXXFLAGS=$(test_cxxflags)
Chemm_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chemm_test_binary80_opt_SOURCES=../common/Chemm.test.cpp $(test_srcdepends)
Chemm_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Chemm_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chemv_test_binary80_SOURCES=../common/Chemv.test.cpp $(test_srcdepends)
Chemv_test_binary80_CXXFLAGS=$(test_cxxflags)
Chemv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chemv_test_binary80_opt_SOURCES=../common/Chemv.test.cpp $(test_srcdepends)
Chemv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Chemv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher_test_binary80_SOURCES=../common/Cher.test.cpp $(test_srcdepends)
Cher_test_binary80_CXXFLAGS=$(test_cxxflags)
Cher_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher_test_binary80_opt_SOURCES=../common/Cher.test.cpp $(test_srcdepends)
Cher_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Cher_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2_test_binary80_SOURCES=../common/Cher2.test.cpp $(test_srcdepends)
Cher2_test_binary80_CXXFLAGS=$(test_cxxflags)
Cher2_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher2_test_binary80_opt_SOURCES=../common/Cher2.test.cpp $(test_srcdepends)
Cher2_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Cher2_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cher2k_test_binary80_SOURCES=../common/Cher2k.test.cpp $(test_srcdepends)
Cher2k_test_binary80_CXXFLAGS=$(test_cxxflags)
Cher2k_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cher2k_test_binary80_opt_SOURCES=../common/Cher2k.test.cpp $(test_srcdepends)
Cher2k_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Cher2k_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cherk_test_binary80_SOURCES=../common/Cherk.test.cpp $(test_srcdepends)
Cherk_test_binary80_CXXFLAGS=$(test_cxxflags)
Cherk_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cherk_test_binary80_opt_SOURCES=../common/Cherk.test.cpp $(test_srcdepends)
Cherk_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Cherk_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpmv_test_binary80_SOURCES=../common/Chpmv.test.cpp $(test_srcdepends)
Chpmv_test_binary80_CXXFLAGS=$(test_cxxflags)
Chpmv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpmv_test_binary80_opt_SOURCES=../common/Chpmv.test.cpp $(test_srcdepends)
Chpmv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Chpmv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr_test_binary80_SOURCES=../common/Chpr.test.cpp $(test_srcdepends)
Chpr_test_binary80_CXXFLAGS=$(test_cxxflags)
Chpr_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpr_test_binary80_opt_SOURCES=../common/Chpr.test.cpp $(test_srcdepends)
Chpr_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Chpr_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Chpr2_test_binary80_SOURCES=../common/Chpr2.test.cpp $(test_srcdepends)
Chpr2_test_binary80_CXXFLAGS=$(test_cxxflags)
Chpr2_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Chpr2_test_binary80_opt_SOURCES=../common/Chpr2.test.cpp $(test_srcdepends)
Chpr2_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Chpr2_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Crotg_test_binary80_SOURCES=../common/Crotg.test.cpp $(test_srcdepends)
Crotg_test_binary80_CXXFLAGS=$(test_cxxflags)
Crotg_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Crotg_test_binary80_opt_SOURCES=../common/Crotg.test.cpp $(test_srcdepends)
Crotg_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Crotg_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cscal_test_binary80_SOURCES=../common/Cscal.test.cpp $(test_srcdepends)
Cscal_test_binary80_CXXFLAGS=$(test_cxxflags)
Cscal_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cscal_test_binary80_opt_SOURCES=../common/Cscal.test.cpp $(test_srcdepends)
Cscal_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Cscal_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Cswap_test_binary80_SOURCES=../common/Cswap.test.cpp $(test_srcdepends)
Cswap_test_binary80_CXXFLAGS=$(test_cxxflags)
Cswap_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Cswap_test_binary80_opt_SOURCES=../common/Cswap.test.cpp $(test_srcdepends)
Cswap_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Cswap_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csymm_test_binary80_SOURCES=../common/Csymm.test.cpp $(test_srcdepends)
Csymm_test_binary80_CXXFLAGS=$(test_cxxflags)
Csymm_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csymm_test_binary80_opt_SOURCES=../common/Csymm.test.cpp $(test_srcdepends)
Csymm_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Csymm_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyr2k_test_binary80_SOURCES=../common/Csyr2k.test.cpp $(test_srcdepends)
Csyr2k_test_binary80_CXXFLAGS=$(test_cxxflags)
Csyr2k_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csyr2k_test_binary80_opt_SOURCES=../common/Csyr2k.test.cpp $(test_srcdepends)
Csyr2k_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Csyr2k_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Csyrk_test_binary80_SOURCES=../common/Csyrk.test.cpp $(test_srcdepends)
Csyrk_test_binary80_CXXFLAGS=$(test_cxxflags)
Csyrk_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Csyrk_test_binary80_opt_SOURCES=../common/Csyrk.test.cpp $(test_srcdepends)
Csyrk_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Csyrk_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbmv_test_binary80_SOURCES=../common/Ctbmv.test.cpp $(test_srcdepends)
Ctbmv_test_binary80_CXXFLAGS=$(test_cxxflags)
Ctbmv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctbmv_test_binary80_opt_SOURCES=../common/Ctbmv.test.cpp $(test_srcdepends)
Ctbmv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Ctbmv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctbsv_test_binary80_SOURCES=../common/Ctbsv.test.cpp $(test_srcdepends)
Ctbsv_test_binary80_CXXFLAGS=$(test_cxxflags)
Ctbsv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctbsv_test_binary80_opt_SOURCES=../common/Ctbsv.test.cpp $(test_srcdepends)
Ctbsv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Ctbsv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpmv_test_binary80_SOURCES=../common/Ctpmv.test.cpp $(test_srcdepends)
Ctpmv_test_binary80_CXXFLAGS=$(test_cxxflags)
Ctpmv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctpmv_test_binary80_opt_SOURCES=../common/Ctpmv.test.cpp $(test_srcdepends)
Ctpmv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Ctpmv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctpsv_test_binary80_SOURCES=../common/Ctpsv.test.cpp $(test_srcdepends)
Ctpsv_test_binary80_CXXFLAGS=$(test_cxxflags)
Ctpsv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctpsv_test_binary80_opt_SOURCES=../common/Ctpsv.test.cpp $(test_srcdepends)
Ctpsv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Ctpsv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmm_test_binary80_SOURCES=../common/Ctrmm.test.cpp $(test_srcdepends)
Ctrmm_test_binary80_CXXFLAGS=$(test_cxxflags)
Ctrmm_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrmm_test_binary80_opt_SOURCES=../common/Ctrmm.test.cpp $(test_srcdepends)
Ctrmm_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Ctrmm_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrmv_test_binary80_SOURCES=../common/Ctrmv.test.cpp $(test_srcdepends)
Ctrmv_test_binary80_CXXFLAGS=$(test_cxxflags)
Ctrmv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrmv_test_binary80_opt_SOURCES=../common/Ctrmv.test.cpp $(test_srcdepends)
Ctrmv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Ctrmv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsm_test_binary80_SOURCES=../common/Ctrsm.test.cpp $(test_srcdepends)
Ctrsm_test_binary80_CXXFLAGS=$(test_cxxflags)
Ctrsm_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrsm_test_binary80_opt_SOURCES=../common/Ctrsm.test.cpp $(test_srcdepends)
Ctrsm_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Ctrsm_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Ctrsv_test_binary80_SOURCES=../common/Ctrsv.test.cpp $(test_srcdepends)
Ctrsv_test_binary80_CXXFLAGS=$(test_cxxflags)
Ctrsv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Ctrsv_test_binary80_opt_SOURCES=../common/Ctrsv.test.cpp $(test_srcdepends)
Ctrsv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Ctrsv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCabs1_test_binary80_opt_SOURCES=../common/RCabs1.test.cpp $(test_srcdepends)
RCabs1_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
RCabs1_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCabs1_test_binary80_SOURCES=../common/RCabs1.test.cpp $(test_srcdepends)
RCabs1_test_binary80_CXXFLAGS=$(test_cxxflags)
RCabs1_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

RCasum_test_binary80_opt_SOURCES=../common/RCasum.test.cpp $(test_srcdepends)
RCasum_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
RCasum_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCasum_test_binary80_SOURCES=../common/RCasum.test.cpp $(test_srcdepends)
RCasum_test_binary80_CXXFLAGS=$(test_cxxflags)
RCasum_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

RCnrm2_test_binary80_opt_SOURCES=../common/RCnrm2.test.cpp $(test_srcdepends)
RCnrm2_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
RCnrm2_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

RCnrm2_test_binary80_SOURCES=../common/RCnrm2.test.cpp $(test_srcdepends)
RCnrm2_test_binary80_CXXFLAGS=$(test_cxxflags)
RCnrm2_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rasum_test_binary80_opt_SOURCES=../common/Rasum.test.cpp $(test_srcdepends)
Rasum_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rasum_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rasum_test_binary80_SOURCES=../common/Rasum.test.cpp $(test_srcdepends)
Rasum_test_binary80_CXXFLAGS=$(test_cxxflags)
Rasum_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Raxpy_test_binary80_opt_SOURCES=../common/Raxpy.test.cpp $(test_srcdepends)
Raxpy_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Raxpy_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Raxpy_test_binary80_SOURCES=../common/Raxpy.test.cpp $(test_srcdepends)
Raxpy_test_binary80_CXXFLAGS=$(test_cxxflags)
Raxpy_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rcopy_test_binary80_opt_SOURCES=../common/Rcopy.test.cpp $(test_srcdepends)
Rcopy_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rcopy_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rcopy_test_binary80_SOURCES=../common/Rcopy.test.cpp $(test_srcdepends)
Rcopy_test_binary80_CXXFLAGS=$(test_cxxflags)
Rcopy_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rdot_test_binary80_opt_SOURCES=../common/Rdot.test.cpp $(test_srcdepends)
Rdot_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rdot_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rdot_test_binary80_SOURCES=../common/Rdot.test.cpp $(test_srcdepends)
Rdot_test_binary80_CXXFLAGS=$(test_cxxflags)
Rdot_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgbmv_test_binary80_opt_SOURCES=../common/Rgbmv.test.cpp $(test_srcdepends)
Rgbmv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rgbmv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgbmv_test_binary80_SOURCES=../common/Rgbmv.test.cpp $(test_srcdepends)
Rgbmv_test_binary80_CXXFLAGS=$(test_cxxflags)
Rgbmv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgemm_test_binary80_opt_SOURCES=../common/Rgemm.test.cpp $(test_srcdepends)
Rgemm_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rgemm_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemm_test_binary80_SOURCES=../common/Rgemm.test.cpp $(test_srcdepends)
Rgemm_test_binary80_CXXFLAGS=$(test_cxxflags)
Rgemm_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rgemv_test_binary80_opt_SOURCES=../common/Rgemv.test.cpp $(test_srcdepends)
Rgemv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rgemv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rgemv_test_binary80_SOURCES=../common/Rgemv.test.cpp $(test_srcdepends)
Rgemv_test_binary80_CXXFLAGS=$(test_cxxflags)
Rgemv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rger_test_binary80_opt_SOURCES=../common/Rger.test.cpp $(test_srcdepends)
Rger_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rger_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rger_test_binary80_SOURCES=../common/Rger.test.cpp $(test_srcdepends)
Rger_test_binary80_CXXFLAGS=$(test_cxxflags)
Rger_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rnrm2_test_binary80_opt_SOURCES=../common/Rnrm2.test.cpp $(test_srcdepends)
Rnrm2_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rnrm2_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rnrm2_test_binary80_SOURCES=../common/Rnrm2.test.cpp $(test_srcdepends)
Rnrm2_test_binary80_CXXFLAGS=$(test_cxxflags)
Rnrm2_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrot_test_binary80_opt_SOURCES=../common/Rrot.test.cpp $(test_srcdepends)
Rrot_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rrot_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrot_test_binary80_SOURCES=../common/Rrot.test.cpp $(test_srcdepends)
Rrot_test_binary80_CXXFLAGS=$(test_cxxflags)
Rrot_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrotg_test_binary80_opt_SOURCES=../common/Rrotg.test.cpp $(test_srcdepends)
Rrotg_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rrotg_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotg_test_binary80_SOURCES=../common/Rrotg.test.cpp $(test_srcdepends)
Rrotg_test_binary80_CXXFLAGS=$(test_cxxflags)
Rrotg_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rrotm_test_binary80_opt_SOURCES=../common/Rrotm.test.cpp $(test_srcdepends)
Rrotm_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rrotm_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rrotm_test_binary80_SOURCES=../common/Rrotm.test.cpp $(test_srcdepends)
Rrotm_test_binary80_CXXFLAGS=$(test_cxxflags)
Rrotm_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsbmv_test_binary80_opt_SOURCES=../common/Rsbmv.test.cpp $(test_srcdepends)
Rsbmv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rsbmv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsbmv_test_binary80_SOURCES=../common/Rsbmv.test.cpp $(test_srcdepends)
Rsbmv_test_binary80_CXXFLAGS=$(test_cxxflags)
Rsbmv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rscal_test_binary80_opt_SOURCES=../common/Rscal.test.cpp $(test_srcdepends)
Rscal_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rscal_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rscal_test_binary80_SOURCES=../common/Rscal.test.cpp $(test_srcdepends)
Rscal_test_binary80_CXXFLAGS=$(test_cxxflags)
Rscal_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspmv_test_binary80_opt_SOURCES=../common/Rspmv.test.cpp $(test_srcdepends)
Rspmv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rspmv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspmv_test_binary80_SOURCES=../common/Rspmv.test.cpp $(test_srcdepends)
Rspmv_test_binary80_CXXFLAGS=$(test_cxxflags)
Rspmv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspr_test_binary80_opt_SOURCES=../common/Rspr.test.cpp $(test_srcdepends)
Rspr_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rspr_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr_test_binary80_SOURCES=../common/Rspr.test.cpp $(test_srcdepends)
Rspr_test_binary80_CXXFLAGS=$(test_cxxflags)
Rspr_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rspr2_test_binary80_opt_SOURCES=../common/Rspr2.test.cpp $(test_srcdepends)
Rspr2_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rspr2_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rspr2_test_binary80_SOURCES=../common/Rspr2.test.cpp $(test_srcdepends)
Rspr2_test_binary80_CXXFLAGS=$(test_cxxflags)
Rspr2_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rswap_test_binary80_opt_SOURCES=../common/Rswap.test.cpp $(test_srcdepends)
Rswap_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rswap_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rswap_test_binary80_SOURCES=../common/Rswap.test.cpp $(test_srcdepends)
Rswap_test_binary80_CXXFLAGS=$(test_cxxflags)
Rswap_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsymm_test_binary80_opt_SOURCES=../common/Rsymm.test.cpp $(test_srcdepends)
Rsymm_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rsymm_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymm_test_binary80_SOURCES=../common/Rsymm.test.cpp $(test_srcdepends)
Rsymm_test_binary80_CXXFLAGS=$(test_cxxflags)
Rsymm_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsymv_test_binary80_opt_SOURCES=../common/Rsymv.test.cpp $(test_srcdepends)
Rsymv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rsymv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsymv_test_binary80_SOURCES=../common/Rsymv.test.cpp $(test_srcdepends)
Rsymv_test_binary80_CXXFLAGS=$(test_cxxflags)
Rsymv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr_test_binary80_opt_SOURCES=../common/Rsyr.test.cpp $(test_srcdepends)
Rsyr_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rsyr_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr_test_binary80_SOURCES=../common/Rsyr.test.cpp $(test_srcdepends)
Rsyr_test_binary80_CXXFLAGS=$(test_cxxflags)
Rsyr_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr2_test_binary80_opt_SOURCES=../common/Rsyr2.test.cpp $(test_srcdepends)
Rsyr2_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rsyr2_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2_test_binary80_SOURCES=../common/Rsyr2.test.cpp $(test_srcdepends)
Rsyr2_test_binary80_CXXFLAGS=$(test_cxxflags)
Rsyr2_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyr2k_test_binary80_opt_SOURCES=../common/Rsyr2k.test.cpp $(test_srcdepends)
Rsyr2k_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rsyr2k_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyr2k_test_binary80_SOURCES=../common/Rsyr2k.test.cpp $(test_srcdepends)
Rsyr2k_test_binary80_CXXFLAGS=$(test_cxxflags)
Rsyr2k_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rsyrk_test_binary80_opt_SOURCES=../common/Rsyrk.test.cpp $(test_srcdepends)
Rsyrk_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rsyrk_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rsyrk_test_binary80_SOURCES=../common/Rsyrk.test.cpp $(test_srcdepends)
Rsyrk_test_binary80_CXXFLAGS=$(test_cxxflags)
Rsyrk_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtbmv_test_binary80_opt_SOURCES=../common/Rtbmv.test.cpp $(test_srcdepends)
Rtbmv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rtbmv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbmv_test_binary80_SOURCES=../common/Rtbmv.test.cpp $(test_srcdepends)
Rtbmv_test_binary80_CXXFLAGS=$(test_cxxflags)
Rtbmv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtbsv_test_binary80_opt_SOURCES=../common/Rtbsv.test.cpp $(test_srcdepends)
Rtbsv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rtbsv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtbsv_test_binary80_SOURCES=../common/Rtbsv.test.cpp $(test_srcdepends)
Rtbsv_test_binary80_CXXFLAGS=$(test_cxxflags)
Rtbsv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtpmv_test_binary80_opt_SOURCES=../common/Rtpmv.test.cpp $(test_srcdepends)
Rtpmv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rtpmv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpmv_test_binary80_SOURCES=../common/Rtpmv.test.cpp $(test_srcdepends)
Rtpmv_test_binary80_CXXFLAGS=$(test_cxxflags)
Rtpmv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtpsv_test_binary80_opt_SOURCES=../common/Rtpsv.test.cpp $(test_srcdepends)
Rtpsv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rtpsv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtpsv_test_binary80_SOURCES=../common/Rtpsv.test.cpp $(test_srcdepends)
Rtpsv_test_binary80_CXXFLAGS=$(test_cxxflags)
Rtpsv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrmm_test_binary80_opt_SOURCES=../common/Rtrmm.test.cpp $(test_srcdepends)
Rtrmm_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rtrmm_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmm_test_binary80_SOURCES=../common/Rtrmm.test.cpp $(test_srcdepends)
Rtrmm_test_binary80_CXXFLAGS=$(test_cxxflags)
Rtrmm_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrmv_test_binary80_opt_SOURCES=../common/Rtrmv.test.cpp $(test_srcdepends)
Rtrmv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rtrmv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrmv_test_binary80_SOURCES=../common/Rtrmv.test.cpp $(test_srcdepends)
Rtrmv_test_binary80_CXXFLAGS=$(test_cxxflags)
Rtrmv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrsm_test_binary80_opt_SOURCES=../common/Rtrsm.test.cpp $(test_srcdepends)
Rtrsm_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rtrsm_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsm_test_binary80_SOURCES=../common/Rtrsm.test.cpp $(test_srcdepends)
Rtrsm_test_binary80_CXXFLAGS=$(test_cxxflags)
Rtrsm_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

Rtrsv_test_binary80_opt_SOURCES=../common/Rtrsv.test.cpp $(test_srcdepends)
Rtrsv_test_binary80_opt_CXXFLAGS=$(test_cxxflags)
Rtrsv_test_binary80_opt_LDFLAGS=$(opt_libdepends) $(mplibs) $(refblas)

Rtrsv_test_binary80_SOURCES=../common/Rtrsv.test.cpp $(test_srcdepends)
Rtrsv_test_binary80_CXXFLAGS=$(test_cxxflags)
Rtrsv_test_binary80_LDFLAGS=$(libdepends) $(mplibs) $(refblas)

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mpblas_binary80_testdir) $(prefix)
endif
