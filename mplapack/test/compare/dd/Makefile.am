check_PROGRAMS = $(mplapack_dd_test_PROGRAMS)
TESTS = run_Rlaruv_test.sh run_Rlamch_test.sh
CLEANFILES = run_Rlaruv_test.sh run_Rlamch_test.sh Rlaruv.txt Rlamch.txt

TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = $(SHELL)

AM_TESTS_ENVIRONMENT = \
    LOG_COMPILER='$(LOG_COMPILER)'; export LOG_COMPILER; \
    EXEEXT='$(EXEEXT)'; export EXEEXT; \
    srcdir='$(srcdir)'; export srcdir;

run_Rlaruv_test.sh: Makefile
	@echo '#!/bin/sh' > $@
	@echo 'NUM_DIFF_PY="$(top_srcdir)/misc/num_diff.py"' >> $@
	@echo '$${LOG_COMPILER} ./Rlaruv.test_dd$${EXEEXT} || true' >> $@
	@echo 'if python3 "$${NUM_DIFF_PY}" "$${srcdir}/Rlaruv_reference.txt" Rlaruv.txt --tol 1e-30; then' >> $@
	@echo '  echo "PASS: Rlaruv"; exit 0' >> $@
	@echo 'else' >> $@
	@echo '  echo "FAIL: Rlaruv"; exit 1' >> $@
	@echo 'fi' >> $@
	@chmod +x $@

run_Rlamch_test.sh: Makefile
	@echo '#!/bin/sh' > $@
	@echo 'NUM_DIFF_PY="$(top_srcdir)/misc/num_diff.py"' >> $@
	@echo '$${LOG_COMPILER} ./Rlamch.test_dd$${EXEEXT} || true' >> $@
	@echo 'if python3 "$${NUM_DIFF_PY}" "$${srcdir}/Rlamch_reference.txt" Rlamch.txt --tol 1e-30; then' >> $@
	@echo '  echo "PASS: Rlamch"; exit 0' >> $@
	@echo 'else' >> $@
	@echo '  echo "FAIL: Rlamch"; exit 1' >> $@
	@echo 'fi' >> $@
	@chmod +x $@

EXTRA_DIST = Rlaruv_reference.txt Rlamch_reference.txt

.PHONY: check-all
check-all: $(check_PROGRAMS)
	$(MAKE) $(AM_MAKEFLAGS) check TESTS="$(mplapack_dd_test_PROGRAMS)"

mplapack_dd_test_PROGRAMS = \
Cgelq2.test_dd \
Cgeqr2.test_dd \
Cgeqrf.test_dd \
Cgesv.test_dd \
Cgetf2.test_dd \
Cgetrf.test_dd \
Cgetri.test_dd \
Cgetrs.test_dd \
Cheev.test_dd \
Chetd2.test_dd \
Chetrd.test_dd \
Clacgv.test_dd \
Clacrm.test_dd \
Clacrt.test_dd \
Cladiv.test_dd \
Claesy.test_dd \
Claev2.test_dd \
Clahef.test_dd \
Clanhe.test_dd \
Clanht.test_dd \
Clansy.test_dd \
Clarfb.test_dd \
Clarfg.test_dd \
Clarf.test_dd \
Clarft.test_dd \
Clartg.test_dd \
Clascl.test_dd \
Claset.test_dd \
Clasr.test_dd \
Classq.test_dd \
Claswp.test_dd \
Clasyf.test_dd \
Clatrd.test_dd \
Cpotf2.test_dd \
Crot.test_dd \
Cspmv.test_dd \
Cspr.test_dd \
Csteqr.test_dd \
Csymv.test_dd \
Csyr.test_dd \
Ctrti2.test_dd \
Ctrtri.test_dd \
Ctrtrs.test_dd \
Cung2l.test_dd \
Cung2r.test_dd \
Cungql.test_dd \
Cungqr.test_dd \
Cungtr.test_dd \
Cunmqr.test_dd \
iCmax1.test_dd \
iMlaenv.test_dd \
Mlsamen.test_dd \
Mutils.test_dd \
RCsum1.test_dd \
Rgecon.test_dd \
Rgeequ.test_dd \
Rgeqlf.test_dd \
Rgesv.test_dd \
Rgetf2.test_dd \
Rgetrf.test_dd \
Rgetri.test_dd \
Rgetrs.test_dd \
Rladiv.test_dd \
Rlae2.test_dd \
Rlaev2.test_dd \
Rlamch.test_dd \
Rlange.test_dd \
Rlanst.test_dd \
Rlansy.test_dd \
Rlapy2.test_dd \
Rlapy3.test_dd \
Rlarfb.test_dd \
Rlarfg.test_dd \
Rlarf.test_dd \
Rlarft.test_dd \
Rlartg.test_dd \
Rlaruv.test_dd \
Rlascl.test_dd \
Rlaset.test_dd \
Rlasr.test_dd \
Rlasrt.test_dd \
Rlassq.test_dd \
Rlaswp.test_dd \
Rlasyf.test_dd \
Rlatrd.test_dd \
Rlatrs.test_dd \
Rlauu2.test_dd \
Rlauum.test_dd \
Rnaninf.test_dd \
Rorg2l.test_dd \
Rorg2r.test_dd \
Rorgql.test_dd \
Rorgqr.test_dd \
Rorgtr.test_dd \
Rpocon.test_dd \
Rposv.test_dd \
Rpotf2.test_dd \
Rpotri.test_dd \
Rpotrs.test_dd \
Rsteqr.test_dd \
Rsterf.test_dd \
Rsyev.test_dd \
Rsytd2.test_dd \
Rsytrd.test_dd \
Rtrti2.test_dd \
Rtrtri.test_dd \
Rtrtrs.test_dd 


mplapack_dd_testdir = $(prefix)/lib/$(target)/mplapack/test/compare/dd

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_dd_testdir) $(prefix)
endif

mplibs=-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(QD_LIBDIR) -lmpc -lmpfr -lgmpxx -lgmp -lqd
refblas=-L$(top_builddir)/external/i/LAPACK/lib -lblas $(FCLIBS)
reflapack=-L$(top_builddir)/external/i/LAPACK/lib -llapack $(FCLIBS)

test_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(QD_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_DD___ -D___MPLAPACK_INTERNAL___ $(DD_CXXFLAGS)
AM_FFLAGS = -x f77-cpp-input -g

test_srcdepends = ../../../../mpblas/test/common/xerbla.cpp ../../../../mpblas/test/common/mplapack.test.cpp

if IS_MACOS
mplibs+=-Wl,-flat_namespace,-undefined,dynamic_lookup
endif

if IS_MINGW
$(check_PROGRAMS): libmxerbla_override.la
lib_LTLIBRARIES = libmxerbla_override.la
libmxerbla_override_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_DD___
libmxerbla_override_la_SOURCES = ../../../../mpblas/test/common/Mxerbla.override.cpp
libmxerbla_override_la_LDFLAGS = -no-undefined

test_libdepends = $(mplibs) $(reflapack) $(refblas) -Wl,--allow-multiple-definition
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,.libs/libmxerbla_override.a,$(top_builddir)/mpblas/reference/.libs/libmpblas_dd.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_dd.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_mpfr.a,$(top_builddir)/mpblas/reference/.libs/libmpblas_mpfr.a,--no-whole-archive
else
test_srcdepends += ../../../../mpblas/test/common/Mxerbla.override.cpp
test_libdepends = $(mplibs) $(reflapack) $(refblas)
libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_dd -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_dd -lmpblas_mpfr
endif

Mutils_test_dd_SOURCES  = ../common/Mutils.test.cpp $(test_srcdepends) ../common/ilaenv.cpp
Mutils_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Mutils_test_dd_CXXFLAGS = $(test_cxxflags)

Cgelq2_test_dd_SOURCES  = ../common/Cgelq2.test.cpp $(test_srcdepends)
Cgelq2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cgelq2_test_dd_CXXFLAGS = $(test_cxxflags)

Cgeqr2_test_dd_SOURCES  = ../common/Cgeqr2.test.cpp $(test_srcdepends)
Cgeqr2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cgeqr2_test_dd_CXXFLAGS = $(test_cxxflags)

Cgeqrf_test_dd_SOURCES  = ../common/Cgeqrf.test.cpp $(test_srcdepends)
Cgeqrf_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cgeqrf_test_dd_CXXFLAGS = $(test_cxxflags)

Cgesv_test_dd_SOURCES  = ../common/Cgesv.test.cpp $(test_srcdepends)
Cgesv_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cgesv_test_dd_CXXFLAGS = $(test_cxxflags)

Cgetf2_test_dd_SOURCES  = ../common/Cgetf2.test.cpp $(test_srcdepends)
Cgetf2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cgetf2_test_dd_CXXFLAGS = $(test_cxxflags)

Cgetrf_test_dd_SOURCES  = ../common/Cgetrf.test.cpp $(test_srcdepends)
Cgetrf_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cgetrf_test_dd_CXXFLAGS = $(test_cxxflags)

Cgetri_test_dd_SOURCES  = ../common/Cgetri.test.cpp $(test_srcdepends)
Cgetri_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cgetri_test_dd_CXXFLAGS = $(test_cxxflags)

Cgetrs_test_dd_SOURCES  = ../common/Cgetrs.test.cpp $(test_srcdepends)
Cgetrs_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cgetrs_test_dd_CXXFLAGS = $(test_cxxflags)

Cheev_test_dd_SOURCES  = ../common/Cheev.test.cpp $(test_srcdepends)
Cheev_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cheev_test_dd_CXXFLAGS = $(test_cxxflags)

Chetd2_test_dd_SOURCES  = ../common/Chetd2.test.cpp $(test_srcdepends)
Chetd2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Chetd2_test_dd_CXXFLAGS = $(test_cxxflags)

Chetrd_test_dd_SOURCES  = ../common/Chetrd.test.cpp $(test_srcdepends)
Chetrd_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Chetrd_test_dd_CXXFLAGS = $(test_cxxflags)

Clacgv_test_dd_SOURCES  = ../common/Clacgv.test.cpp $(test_srcdepends)
Clacgv_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clacgv_test_dd_CXXFLAGS = $(test_cxxflags)

Clacrm_test_dd_SOURCES  = ../common/Clacrm.test.cpp $(test_srcdepends)
Clacrm_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clacrm_test_dd_CXXFLAGS = $(test_cxxflags)

Clacrt_test_dd_SOURCES  = ../common/Clacrt.test.cpp $(test_srcdepends)
Clacrt_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clacrt_test_dd_CXXFLAGS = $(test_cxxflags)

Cladiv_test_dd_SOURCES  = ../common/Cladiv.test.cpp $(test_srcdepends)
Cladiv_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cladiv_test_dd_CXXFLAGS = $(test_cxxflags)

Claesy_test_dd_SOURCES  = ../common/Claesy.test.cpp $(test_srcdepends)
Claesy_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Claesy_test_dd_CXXFLAGS = $(test_cxxflags)

Claev2_test_dd_SOURCES  = ../common/Claev2.test.cpp $(test_srcdepends)
Claev2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Claev2_test_dd_CXXFLAGS = $(test_cxxflags)

Clahef_test_dd_SOURCES  = ../common/Clahef.test.cpp $(test_srcdepends)
Clahef_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clahef_test_dd_CXXFLAGS = $(test_cxxflags)

Clanhe_test_dd_SOURCES  = ../common/Clanhe.test.cpp $(test_srcdepends)
Clanhe_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clanhe_test_dd_CXXFLAGS = $(test_cxxflags)

Clanht_test_dd_SOURCES  = ../common/Clanht.test.cpp $(test_srcdepends)
Clanht_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clanht_test_dd_CXXFLAGS = $(test_cxxflags)

Clansy_test_dd_SOURCES  = ../common/Clansy.test.cpp $(test_srcdepends)
Clansy_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clansy_test_dd_CXXFLAGS = $(test_cxxflags)

Clarfb_test_dd_SOURCES  = ../common/Clarfb.test.cpp $(test_srcdepends)
Clarfb_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clarfb_test_dd_CXXFLAGS = $(test_cxxflags)

Clarfg_test_dd_SOURCES  = ../common/Clarfg.test.cpp $(test_srcdepends)
Clarfg_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clarfg_test_dd_CXXFLAGS = $(test_cxxflags)

Clarf_test_dd_SOURCES  = ../common/Clarf.test.cpp $(test_srcdepends)
Clarf_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clarf_test_dd_CXXFLAGS = $(test_cxxflags)

Clarft_test_dd_SOURCES  = ../common/Clarft.test.cpp $(test_srcdepends)
Clarft_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clarft_test_dd_CXXFLAGS = $(test_cxxflags)

Clartg_test_dd_SOURCES  = ../common/Clartg.test.cpp $(test_srcdepends)
Clartg_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clartg_test_dd_CXXFLAGS = $(test_cxxflags)

Clascl_test_dd_SOURCES  = ../common/Clascl.test.cpp $(test_srcdepends)
Clascl_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clascl_test_dd_CXXFLAGS = $(test_cxxflags)

Claset_test_dd_SOURCES  = ../common/Claset.test.cpp $(test_srcdepends)
Claset_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Claset_test_dd_CXXFLAGS = $(test_cxxflags)

Clasr_test_dd_SOURCES  = ../common/Clasr.test.cpp $(test_srcdepends)
Clasr_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clasr_test_dd_CXXFLAGS = $(test_cxxflags)

Classq_test_dd_SOURCES  = ../common/Classq.test.cpp $(test_srcdepends)
Classq_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Classq_test_dd_CXXFLAGS = $(test_cxxflags)

Claswp_test_dd_SOURCES  = ../common/Claswp.test.cpp $(test_srcdepends)
Claswp_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Claswp_test_dd_CXXFLAGS = $(test_cxxflags)

Clasyf_test_dd_SOURCES  = ../common/Clasyf.test.cpp $(test_srcdepends)
Clasyf_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clasyf_test_dd_CXXFLAGS = $(test_cxxflags)

Clatrd_test_dd_SOURCES  = ../common/Clatrd.test.cpp $(test_srcdepends)
Clatrd_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Clatrd_test_dd_CXXFLAGS = $(test_cxxflags)

Cpotf2_test_dd_SOURCES  = ../common/Cpotf2.test.cpp $(test_srcdepends)
Cpotf2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cpotf2_test_dd_CXXFLAGS = $(test_cxxflags)

Crot_test_dd_SOURCES  = ../common/Crot.test.cpp $(test_srcdepends)
Crot_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Crot_test_dd_CXXFLAGS = $(test_cxxflags)

Cspmv_test_dd_SOURCES  = ../common/Cspmv.test.cpp $(test_srcdepends)
Cspmv_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cspmv_test_dd_CXXFLAGS = $(test_cxxflags)

Cspr_test_dd_SOURCES  = ../common/Cspr.test.cpp $(test_srcdepends)
Cspr_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cspr_test_dd_CXXFLAGS = $(test_cxxflags)

Csteqr_test_dd_SOURCES  = ../common/Csteqr.test.cpp $(test_srcdepends)
Csteqr_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Csteqr_test_dd_CXXFLAGS = $(test_cxxflags)

Csymv_test_dd_SOURCES  = ../common/Csymv.test.cpp $(test_srcdepends)
Csymv_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Csymv_test_dd_CXXFLAGS = $(test_cxxflags)

Csyr_test_dd_SOURCES  = ../common/Csyr.test.cpp $(test_srcdepends)
Csyr_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Csyr_test_dd_CXXFLAGS = $(test_cxxflags)

Ctrti2_test_dd_SOURCES  = ../common/Ctrti2.test.cpp $(test_srcdepends)
Ctrti2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Ctrti2_test_dd_CXXFLAGS = $(test_cxxflags)

Ctrtri_test_dd_SOURCES  = ../common/Ctrtri.test.cpp $(test_srcdepends)
Ctrtri_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Ctrtri_test_dd_CXXFLAGS = $(test_cxxflags)

Ctrtrs_test_dd_SOURCES  = ../common/Ctrtrs.test.cpp $(test_srcdepends)
Ctrtrs_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Ctrtrs_test_dd_CXXFLAGS = $(test_cxxflags)

Cung2l_test_dd_SOURCES  = ../common/Cung2l.test.cpp $(test_srcdepends)
Cung2l_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cung2l_test_dd_CXXFLAGS = $(test_cxxflags)

Cung2r_test_dd_SOURCES  = ../common/Cung2r.test.cpp $(test_srcdepends)
Cung2r_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cung2r_test_dd_CXXFLAGS = $(test_cxxflags)

Cungql_test_dd_SOURCES  = ../common/Cungql.test.cpp $(test_srcdepends)
Cungql_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cungql_test_dd_CXXFLAGS = $(test_cxxflags)

Cungqr_test_dd_SOURCES  = ../common/Cungqr.test.cpp $(test_srcdepends)
Cungqr_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cungqr_test_dd_CXXFLAGS = $(test_cxxflags)

Cungtr_test_dd_SOURCES  = ../common/Cungtr.test.cpp $(test_srcdepends)
Cungtr_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cungtr_test_dd_CXXFLAGS = $(test_cxxflags)

Cunmqr_test_dd_SOURCES  = ../common/Cunmqr.test.cpp $(test_srcdepends)
Cunmqr_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Cunmqr_test_dd_CXXFLAGS = $(test_cxxflags)

iCmax1_test_dd_SOURCES  = ../common/iCmax1.test.cpp $(test_srcdepends)
iCmax1_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
iCmax1_test_dd_CXXFLAGS = $(test_cxxflags)

iMlaenv_test_dd_SOURCES  = ../common/iMlaenv.test.cpp $(test_srcdepends)
iMlaenv_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
iMlaenv_test_dd_CXXFLAGS = $(test_cxxflags)

Mlsamen_test_dd_SOURCES  = ../common/Mlsamen.test.cpp $(test_srcdepends)
Mlsamen_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Mlsamen_test_dd_CXXFLAGS = $(test_cxxflags)

RCsum1_test_dd_SOURCES  = ../common/RCsum1.test.cpp $(test_srcdepends)
RCsum1_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
RCsum1_test_dd_CXXFLAGS = $(test_cxxflags)

Rgecon_test_dd_SOURCES  = ../common/Rgecon.test.cpp $(test_srcdepends)
Rgecon_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rgecon_test_dd_CXXFLAGS = $(test_cxxflags)

Rgeequ_test_dd_SOURCES  = ../common/Rgeequ.test.cpp $(test_srcdepends)
Rgeequ_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rgeequ_test_dd_CXXFLAGS = $(test_cxxflags)

Rgeqlf_test_dd_SOURCES  = ../common/Rgeqlf.test.cpp $(test_srcdepends)
Rgeqlf_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rgeqlf_test_dd_CXXFLAGS = $(test_cxxflags)

Rgesv_test_dd_SOURCES  = ../common/Rgesv.test.cpp $(test_srcdepends)
Rgesv_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rgesv_test_dd_CXXFLAGS = $(test_cxxflags)

Rgetf2_test_dd_SOURCES  = ../common/Rgetf2.test.cpp $(test_srcdepends)
Rgetf2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rgetf2_test_dd_CXXFLAGS = $(test_cxxflags)

Rgetrf_test_dd_SOURCES  = ../common/Rgetrf.test.cpp $(test_srcdepends)
Rgetrf_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rgetrf_test_dd_CXXFLAGS = $(test_cxxflags)

Rgetri_test_dd_SOURCES  = ../common/Rgetri.test.cpp $(test_srcdepends)
Rgetri_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rgetri_test_dd_CXXFLAGS = $(test_cxxflags)

Rgetrs_test_dd_SOURCES  = ../common/Rgetrs.test.cpp $(test_srcdepends)
Rgetrs_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rgetrs_test_dd_CXXFLAGS = $(test_cxxflags)

Rladiv_test_dd_SOURCES  = ../common/Rladiv.test.cpp $(test_srcdepends)
Rladiv_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rladiv_test_dd_CXXFLAGS = $(test_cxxflags)

Rlae2_test_dd_SOURCES  = ../common/Rlae2.test.cpp $(test_srcdepends)
Rlae2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlae2_test_dd_CXXFLAGS = $(test_cxxflags)

Rlaev2_test_dd_SOURCES  = ../common/Rlaev2.test.cpp $(test_srcdepends)
Rlaev2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlaev2_test_dd_CXXFLAGS = $(test_cxxflags)

Rlamch_test_dd_SOURCES  = ../common/Rlamch.test.cpp $(test_srcdepends)
Rlamch_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlamch_test_dd_CXXFLAGS = $(test_cxxflags)

Rlange_test_dd_SOURCES  = ../common/Rlange.test.cpp $(test_srcdepends)
Rlange_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlange_test_dd_CXXFLAGS = $(test_cxxflags)

Rlanst_test_dd_SOURCES  = ../common/Rlanst.test.cpp $(test_srcdepends)
Rlanst_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlanst_test_dd_CXXFLAGS = $(test_cxxflags)

Rlansy_test_dd_SOURCES  = ../common/Rlansy.test.cpp $(test_srcdepends)
Rlansy_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlansy_test_dd_CXXFLAGS = $(test_cxxflags)

Rlapy2_test_dd_SOURCES  = ../common/Rlapy2.test.cpp $(test_srcdepends)
Rlapy2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlapy2_test_dd_CXXFLAGS = $(test_cxxflags)

Rlapy3_test_dd_SOURCES  = ../common/Rlapy3.test.cpp $(test_srcdepends)
Rlapy3_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlapy3_test_dd_CXXFLAGS = $(test_cxxflags)

Rlarfb_test_dd_SOURCES  = ../common/Rlarfb.test.cpp $(test_srcdepends)
Rlarfb_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlarfb_test_dd_CXXFLAGS = $(test_cxxflags)

Rlarfg_test_dd_SOURCES  = ../common/Rlarfg.test.cpp $(test_srcdepends)
Rlarfg_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlarfg_test_dd_CXXFLAGS = $(test_cxxflags)

Rlarf_test_dd_SOURCES  = ../common/Rlarf.test.cpp $(test_srcdepends)
Rlarf_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlarf_test_dd_CXXFLAGS = $(test_cxxflags)

Rlarft_test_dd_SOURCES  = ../common/Rlarft.test.cpp $(test_srcdepends)
Rlarft_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlarft_test_dd_CXXFLAGS = $(test_cxxflags)

Rlartg_test_dd_SOURCES  = ../common/Rlartg.test.cpp $(test_srcdepends)
Rlartg_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlartg_test_dd_CXXFLAGS = $(test_cxxflags)

Rlaruv_test_dd_SOURCES  = ../common/Rlaruv.test.cpp $(test_srcdepends)
Rlaruv_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlaruv_test_dd_CXXFLAGS = $(test_cxxflags)

Rlascl_test_dd_SOURCES  = ../common/Rlascl.test.cpp $(test_srcdepends)
Rlascl_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlascl_test_dd_CXXFLAGS = $(test_cxxflags)

Rlaset_test_dd_SOURCES  = ../common/Rlaset.test.cpp $(test_srcdepends)
Rlaset_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlaset_test_dd_CXXFLAGS = $(test_cxxflags)

Rlasr_test_dd_SOURCES  = ../common/Rlasr.test.cpp $(test_srcdepends)
Rlasr_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlasr_test_dd_CXXFLAGS = $(test_cxxflags)

Rlasrt_test_dd_SOURCES  = ../common/Rlasrt.test.cpp $(test_srcdepends)
Rlasrt_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlasrt_test_dd_CXXFLAGS = $(test_cxxflags)

Rlassq_test_dd_SOURCES  = ../common/Rlassq.test.cpp $(test_srcdepends)
Rlassq_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlassq_test_dd_CXXFLAGS = $(test_cxxflags)

Rlaswp_test_dd_SOURCES  = ../common/Rlaswp.test.cpp $(test_srcdepends)
Rlaswp_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlaswp_test_dd_CXXFLAGS = $(test_cxxflags)

Rlasyf_test_dd_SOURCES  = ../common/Rlasyf.test.cpp $(test_srcdepends)
Rlasyf_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlasyf_test_dd_CXXFLAGS = $(test_cxxflags)

Rlatrd_test_dd_SOURCES  = ../common/Rlatrd.test.cpp $(test_srcdepends)
Rlatrd_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlatrd_test_dd_CXXFLAGS = $(test_cxxflags)

Rlatrs_test_dd_SOURCES  = ../common/Rlatrs.test.cpp $(test_srcdepends)
Rlatrs_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlatrs_test_dd_CXXFLAGS = $(test_cxxflags)

Rlauu2_test_dd_SOURCES  = ../common/Rlauu2.test.cpp $(test_srcdepends)
Rlauu2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlauu2_test_dd_CXXFLAGS = $(test_cxxflags)

Rlauum_test_dd_SOURCES  = ../common/Rlauum.test.cpp $(test_srcdepends)
Rlauum_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rlauum_test_dd_CXXFLAGS = $(test_cxxflags)

Rnaninf_test_dd_SOURCES  = ../common/Rnaninf.test.cpp $(test_srcdepends)
Rnaninf_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rnaninf_test_dd_CXXFLAGS = $(test_cxxflags)

Rorg2l_test_dd_SOURCES  = ../common/Rorg2l.test.cpp $(test_srcdepends)
Rorg2l_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rorg2l_test_dd_CXXFLAGS = $(test_cxxflags)

Rorg2r_test_dd_SOURCES  = ../common/Rorg2r.test.cpp $(test_srcdepends)
Rorg2r_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rorg2r_test_dd_CXXFLAGS = $(test_cxxflags)

Rorgql_test_dd_SOURCES  = ../common/Rorgql.test.cpp $(test_srcdepends)
Rorgql_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rorgql_test_dd_CXXFLAGS = $(test_cxxflags)

Rorgqr_test_dd_SOURCES  = ../common/Rorgqr.test.cpp $(test_srcdepends)
Rorgqr_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rorgqr_test_dd_CXXFLAGS = $(test_cxxflags)

Rorgtr_test_dd_SOURCES  = ../common/Rorgtr.test.cpp $(test_srcdepends)
Rorgtr_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rorgtr_test_dd_CXXFLAGS = $(test_cxxflags)

Rpocon_test_dd_SOURCES  = ../common/Rpocon.test.cpp $(test_srcdepends)
Rpocon_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rpocon_test_dd_CXXFLAGS = $(test_cxxflags)

Rposv_test_dd_SOURCES  = ../common/Rposv.test.cpp $(test_srcdepends)
Rposv_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rposv_test_dd_CXXFLAGS = $(test_cxxflags)

Rpotf2_test_dd_SOURCES  = ../common/Rpotf2.test.cpp $(test_srcdepends)
Rpotf2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rpotf2_test_dd_CXXFLAGS = $(test_cxxflags)

Rpotri_test_dd_SOURCES  = ../common/Rpotri.test.cpp $(test_srcdepends)
Rpotri_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rpotri_test_dd_CXXFLAGS = $(test_cxxflags)

Rpotrs_test_dd_SOURCES  = ../common/Rpotrs.test.cpp $(test_srcdepends)
Rpotrs_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rpotrs_test_dd_CXXFLAGS = $(test_cxxflags)

Rsteqr_test_dd_SOURCES  = ../common/Rsteqr.test.cpp $(test_srcdepends)
Rsteqr_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rsteqr_test_dd_CXXFLAGS = $(test_cxxflags)

Rsterf_test_dd_SOURCES  = ../common/Rsterf.test.cpp $(test_srcdepends)
Rsterf_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rsterf_test_dd_CXXFLAGS = $(test_cxxflags)

Rsyev_test_dd_SOURCES  = ../common/Rsyev.test.cpp $(test_srcdepends)
Rsyev_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rsyev_test_dd_CXXFLAGS = $(test_cxxflags)

Rsytd2_test_dd_SOURCES  = ../common/Rsytd2.test.cpp $(test_srcdepends)
Rsytd2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rsytd2_test_dd_CXXFLAGS = $(test_cxxflags)

Rsytrd_test_dd_SOURCES  = ../common/Rsytrd.test.cpp $(test_srcdepends)
Rsytrd_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rsytrd_test_dd_CXXFLAGS = $(test_cxxflags)

Rtrti2_test_dd_SOURCES  = ../common/Rtrti2.test.cpp $(test_srcdepends)
Rtrti2_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rtrti2_test_dd_CXXFLAGS = $(test_cxxflags)

Rtrtri_test_dd_SOURCES  = ../common/Rtrtri.test.cpp $(test_srcdepends)
Rtrtri_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rtrtri_test_dd_CXXFLAGS = $(test_cxxflags)

Rtrtrs_test_dd_SOURCES  = ../common/Rtrtrs.test.cpp $(test_srcdepends)
Rtrtrs_test_dd_LDFLAGS  = $(libdepends) $(test_libdepends)
Rtrtrs_test_dd_CXXFLAGS = $(test_cxxflags)



