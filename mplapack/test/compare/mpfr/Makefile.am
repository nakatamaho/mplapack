MPFR_COMPARE_TEST_ENV = MPFRXX_DEFAULT_EMAX=32768 MPFRXX_DEFAULT_EMIN=-32765
backend_tests_environment = \
    MPFRXX_DEFAULT_EMAX=32768; export MPFRXX_DEFAULT_EMAX; \
    MPFRXX_DEFAULT_EMIN=-32765; export MPFRXX_DEFAULT_EMIN;

run_Rlaruv_test.sh: Makefile
	@echo '#!/bin/sh' > $@
	@echo '$(MPFR_COMPARE_TEST_ENV) $${LOG_COMPILER} ./Rlaruv.test_mpfr$${EXEEXT} || true' >> $@
	@echo 'if diff -u $${srcdir}/Rlaruv_reference.txt Rlaruv.txt > /dev/null 2>&1; then' >> $@
	@echo '  echo "PASS: Rlaruv"; exit 0' >> $@
	@echo 'else' >> $@
	@echo '  echo "FAIL: Rlaruv"; diff $${srcdir}/Rlaruv_reference.txt Rlaruv.txt; exit 1' >> $@
	@echo 'fi' >> $@
	@chmod +x $@

run_Rlamch_test.sh: Makefile
	@echo '#!/bin/sh' > $@
	@echo '$(MPFR_COMPARE_TEST_ENV) $${LOG_COMPILER} ./Rlamch.test_mpfr$${EXEEXT} || true' >> $@
	@echo 'if diff -u $${srcdir}/Rlamch_reference.txt Rlamch.txt > /dev/null 2>&1; then' >> $@
	@echo '  echo "PASS: Rlamch"; exit 0' >> $@
	@echo 'else' >> $@
	@echo '  echo "FAIL: Rlamch"; diff $${srcdir}/Rlamch_reference.txt Rlamch.txt; exit 1' >> $@
	@echo 'fi' >> $@
	@chmod +x $@

EXTRA_DIST = Rlaruv_reference.txt Rlamch_reference.txt

backend_define = $(GMPFRXX_MKII_CPPFLAGS) -DMPLAPACK_BUILD_WITH_MPFR
backend_cppflags = -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) $(backend_define) -DMPLAPACK_INTERNAL
backend_cxxflags =
backend_reference_archives = $(top_builddir)/mplapack/reference/.libs/libmplapack_mpfr.a
backend_reference_libs = -lmplapack_mpfr
backend_whole_archive_libs =
backend_link_libs =
backend_testdir = $(prefix)/lib/$(target)/mplapack/test/compare/mpfr
mplibs = -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpc -lmpfr -lgmp

compare_tests = \
    Cgelq2.test_mpfr \
    Cgeqr2.test_mpfr \
    Cgeqrf.test_mpfr \
    Cgesv.test_mpfr \
    Cgetf2.test_mpfr \
    Cgetrf.test_mpfr \
    Cgetri.test_mpfr \
    Cgetrs.test_mpfr \
    Cheev.test_mpfr \
    Chetd2.test_mpfr \
    Chetrd.test_mpfr \
    Clacgv.test_mpfr \
    Clacrm.test_mpfr \
    Clacrt.test_mpfr \
    Cladiv.test_mpfr \
    Claesy.test_mpfr \
    Claev2.test_mpfr \
    Clahef.test_mpfr \
    Clanhe.test_mpfr \
    Clanht.test_mpfr \
    Clansy.test_mpfr \
    Clarf.test_mpfr \
    Clarfb.test_mpfr \
    Clarfg.test_mpfr \
    Clarft.test_mpfr \
    Clartg.test_mpfr \
    Clascl.test_mpfr \
    Claset.test_mpfr \
    Clasr.test_mpfr \
    Classq.test_mpfr \
    Claswp.test_mpfr \
    Clasyf.test_mpfr \
    Clatrd.test_mpfr \
    Cpotf2.test_mpfr \
    Crot.test_mpfr \
    Cspmv.test_mpfr \
    Cspr.test_mpfr \
    Csteqr.test_mpfr \
    Csymv.test_mpfr \
    Csyr.test_mpfr \
    Ctrti2.test_mpfr \
    Ctrtri.test_mpfr \
    Ctrtrs.test_mpfr \
    Cung2l.test_mpfr \
    Cung2r.test_mpfr \
    Cungql.test_mpfr \
    Cungqr.test_mpfr \
    Cungtr.test_mpfr \
    Cunmqr.test_mpfr \
    Mlsamen.test_mpfr \
    Mutils.test_mpfr \
    RCsum1.test_mpfr \
    Rgecon.test_mpfr \
    Rgeequ.test_mpfr \
    Rgeqlf.test_mpfr \
    Rgesv.test_mpfr \
    Rgetf2.test_mpfr \
    Rgetrf.test_mpfr \
    Rgetri.test_mpfr \
    Rgetrs.test_mpfr \
    Rladiv.test_mpfr \
    Rlae2.test_mpfr \
    Rlaev2.test_mpfr \
    Rlamch.test_mpfr \
    Rlange.test_mpfr \
    Rlanst.test_mpfr \
    Rlansy.test_mpfr \
    Rlapy2.test_mpfr \
    Rlapy3.test_mpfr \
    Rlarf.test_mpfr \
    Rlarfb.test_mpfr \
    Rlarfg.test_mpfr \
    Rlarft.test_mpfr \
    Rlartg.test_mpfr \
    Rlaruv.test_mpfr \
    Rlascl.test_mpfr \
    Rlaset.test_mpfr \
    Rlasr.test_mpfr \
    Rlasrt.test_mpfr \
    Rlassq.test_mpfr \
    Rlaswp.test_mpfr \
    Rlasyf.test_mpfr \
    Rlatrd.test_mpfr \
    Rlatrs.test_mpfr \
    Rlauu2.test_mpfr \
    Rlauum.test_mpfr \
    Rnaninf.test_mpfr \
    Rorg2l.test_mpfr \
    Rorg2r.test_mpfr \
    Rorgql.test_mpfr \
    Rorgqr.test_mpfr \
    Rorgtr.test_mpfr \
    Rpocon.test_mpfr \
    Rposv.test_mpfr \
    Rpotf2.test_mpfr \
    Rpotri.test_mpfr \
    Rpotrs.test_mpfr \
    Rsteqr.test_mpfr \
    Rsterf.test_mpfr \
    Rsyev.test_mpfr \
    Rsytd2.test_mpfr \
    Rsytrd.test_mpfr \
    Rtrti2.test_mpfr \
    Rtrtri.test_mpfr \
    Rtrtrs.test_mpfr \
    iCmax1.test_mpfr \
    iMlaenv.test_mpfr

Cgelq2_test_mpfr_SOURCES = ../common/Cgelq2.test.cpp $(test_srcdepends)
Cgelq2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cgeqr2_test_mpfr_SOURCES = ../common/Cgeqr2.test.cpp $(test_srcdepends)
Cgeqr2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cgeqrf_test_mpfr_SOURCES = ../common/Cgeqrf.test.cpp $(test_srcdepends)
Cgeqrf_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cgesv_test_mpfr_SOURCES = ../common/Cgesv.test.cpp $(test_srcdepends)
Cgesv_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cgetf2_test_mpfr_SOURCES = ../common/Cgetf2.test.cpp $(test_srcdepends)
Cgetf2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cgetrf_test_mpfr_SOURCES = ../common/Cgetrf.test.cpp $(test_srcdepends)
Cgetrf_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cgetri_test_mpfr_SOURCES = ../common/Cgetri.test.cpp $(test_srcdepends)
Cgetri_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cgetrs_test_mpfr_SOURCES = ../common/Cgetrs.test.cpp $(test_srcdepends)
Cgetrs_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cheev_test_mpfr_SOURCES = ../common/Cheev.test.cpp $(test_srcdepends)
Cheev_test_mpfr_CXXFLAGS = $(test_cxxflags)
Chetd2_test_mpfr_SOURCES = ../common/Chetd2.test.cpp $(test_srcdepends)
Chetd2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Chetrd_test_mpfr_SOURCES = ../common/Chetrd.test.cpp $(test_srcdepends)
Chetrd_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clacgv_test_mpfr_SOURCES = ../common/Clacgv.test.cpp $(test_srcdepends)
Clacgv_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clacrm_test_mpfr_SOURCES = ../common/Clacrm.test.cpp $(test_srcdepends)
Clacrm_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clacrt_test_mpfr_SOURCES = ../common/Clacrt.test.cpp $(test_srcdepends)
Clacrt_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cladiv_test_mpfr_SOURCES = ../common/Cladiv.test.cpp $(test_srcdepends)
Cladiv_test_mpfr_CXXFLAGS = $(test_cxxflags)
Claesy_test_mpfr_SOURCES = ../common/Claesy.test.cpp $(test_srcdepends)
Claesy_test_mpfr_CXXFLAGS = $(test_cxxflags)
Claev2_test_mpfr_SOURCES = ../common/Claev2.test.cpp $(test_srcdepends)
Claev2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clahef_test_mpfr_SOURCES = ../common/Clahef.test.cpp $(test_srcdepends)
Clahef_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clanhe_test_mpfr_SOURCES = ../common/Clanhe.test.cpp $(test_srcdepends)
Clanhe_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clanht_test_mpfr_SOURCES = ../common/Clanht.test.cpp $(test_srcdepends)
Clanht_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clansy_test_mpfr_SOURCES = ../common/Clansy.test.cpp $(test_srcdepends)
Clansy_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clarf_test_mpfr_SOURCES = ../common/Clarf.test.cpp $(test_srcdepends)
Clarf_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clarfb_test_mpfr_SOURCES = ../common/Clarfb.test.cpp $(test_srcdepends)
Clarfb_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clarfg_test_mpfr_SOURCES = ../common/Clarfg.test.cpp $(test_srcdepends)
Clarfg_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clarft_test_mpfr_SOURCES = ../common/Clarft.test.cpp $(test_srcdepends)
Clarft_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clartg_test_mpfr_SOURCES = ../common/Clartg.test.cpp $(test_srcdepends)
Clartg_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clascl_test_mpfr_SOURCES = ../common/Clascl.test.cpp $(test_srcdepends)
Clascl_test_mpfr_CXXFLAGS = $(test_cxxflags)
Claset_test_mpfr_SOURCES = ../common/Claset.test.cpp $(test_srcdepends)
Claset_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clasr_test_mpfr_SOURCES = ../common/Clasr.test.cpp $(test_srcdepends)
Clasr_test_mpfr_CXXFLAGS = $(test_cxxflags)
Classq_test_mpfr_SOURCES = ../common/Classq.test.cpp $(test_srcdepends)
Classq_test_mpfr_CXXFLAGS = $(test_cxxflags)
Claswp_test_mpfr_SOURCES = ../common/Claswp.test.cpp $(test_srcdepends)
Claswp_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clasyf_test_mpfr_SOURCES = ../common/Clasyf.test.cpp $(test_srcdepends)
Clasyf_test_mpfr_CXXFLAGS = $(test_cxxflags)
Clatrd_test_mpfr_SOURCES = ../common/Clatrd.test.cpp $(test_srcdepends)
Clatrd_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cpotf2_test_mpfr_SOURCES = ../common/Cpotf2.test.cpp $(test_srcdepends)
Cpotf2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Crot_test_mpfr_SOURCES = ../common/Crot.test.cpp $(test_srcdepends)
Crot_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cspmv_test_mpfr_SOURCES = ../common/Cspmv.test.cpp $(test_srcdepends)
Cspmv_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cspr_test_mpfr_SOURCES = ../common/Cspr.test.cpp $(test_srcdepends)
Cspr_test_mpfr_CXXFLAGS = $(test_cxxflags)
Csteqr_test_mpfr_SOURCES = ../common/Csteqr.test.cpp $(test_srcdepends)
Csteqr_test_mpfr_CXXFLAGS = $(test_cxxflags)
Csymv_test_mpfr_SOURCES = ../common/Csymv.test.cpp $(test_srcdepends)
Csymv_test_mpfr_CXXFLAGS = $(test_cxxflags)
Csyr_test_mpfr_SOURCES = ../common/Csyr.test.cpp $(test_srcdepends)
Csyr_test_mpfr_CXXFLAGS = $(test_cxxflags)
Ctrti2_test_mpfr_SOURCES = ../common/Ctrti2.test.cpp $(test_srcdepends)
Ctrti2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Ctrtri_test_mpfr_SOURCES = ../common/Ctrtri.test.cpp $(test_srcdepends)
Ctrtri_test_mpfr_CXXFLAGS = $(test_cxxflags)
Ctrtrs_test_mpfr_SOURCES = ../common/Ctrtrs.test.cpp $(test_srcdepends)
Ctrtrs_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cung2l_test_mpfr_SOURCES = ../common/Cung2l.test.cpp $(test_srcdepends)
Cung2l_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cung2r_test_mpfr_SOURCES = ../common/Cung2r.test.cpp $(test_srcdepends)
Cung2r_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cungql_test_mpfr_SOURCES = ../common/Cungql.test.cpp $(test_srcdepends)
Cungql_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cungqr_test_mpfr_SOURCES = ../common/Cungqr.test.cpp $(test_srcdepends)
Cungqr_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cungtr_test_mpfr_SOURCES = ../common/Cungtr.test.cpp $(test_srcdepends)
Cungtr_test_mpfr_CXXFLAGS = $(test_cxxflags)
Cunmqr_test_mpfr_SOURCES = ../common/Cunmqr.test.cpp $(test_srcdepends)
Cunmqr_test_mpfr_CXXFLAGS = $(test_cxxflags)
Mlsamen_test_mpfr_SOURCES = ../common/Mlsamen.test.cpp $(test_srcdepends)
Mlsamen_test_mpfr_CXXFLAGS = $(test_cxxflags)
Mutils_test_mpfr_SOURCES = ../common/Mutils.test.cpp $(test_srcdepends) ../common/ilaenv.cpp
Mutils_test_mpfr_CXXFLAGS = $(test_cxxflags)
RCsum1_test_mpfr_SOURCES = ../common/RCsum1.test.cpp $(test_srcdepends)
RCsum1_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rgecon_test_mpfr_SOURCES = ../common/Rgecon.test.cpp $(test_srcdepends)
Rgecon_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rgeequ_test_mpfr_SOURCES = ../common/Rgeequ.test.cpp $(test_srcdepends)
Rgeequ_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rgeqlf_test_mpfr_SOURCES = ../common/Rgeqlf.test.cpp $(test_srcdepends)
Rgeqlf_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rgesv_test_mpfr_SOURCES = ../common/Rgesv.test.cpp $(test_srcdepends)
Rgesv_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rgetf2_test_mpfr_SOURCES = ../common/Rgetf2.test.cpp $(test_srcdepends)
Rgetf2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rgetrf_test_mpfr_SOURCES = ../common/Rgetrf.test.cpp $(test_srcdepends)
Rgetrf_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rgetri_test_mpfr_SOURCES = ../common/Rgetri.test.cpp $(test_srcdepends)
Rgetri_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rgetrs_test_mpfr_SOURCES = ../common/Rgetrs.test.cpp $(test_srcdepends)
Rgetrs_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rladiv_test_mpfr_SOURCES = ../common/Rladiv.test.cpp $(test_srcdepends)
Rladiv_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlae2_test_mpfr_SOURCES = ../common/Rlae2.test.cpp $(test_srcdepends)
Rlae2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlaev2_test_mpfr_SOURCES = ../common/Rlaev2.test.cpp $(test_srcdepends)
Rlaev2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlamch_test_mpfr_SOURCES = ../common/Rlamch.test.cpp $(test_srcdepends)
Rlamch_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlange_test_mpfr_SOURCES = ../common/Rlange.test.cpp $(test_srcdepends)
Rlange_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlanst_test_mpfr_SOURCES = ../common/Rlanst.test.cpp $(test_srcdepends)
Rlanst_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlansy_test_mpfr_SOURCES = ../common/Rlansy.test.cpp $(test_srcdepends)
Rlansy_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlapy2_test_mpfr_SOURCES = ../common/Rlapy2.test.cpp $(test_srcdepends)
Rlapy2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlapy3_test_mpfr_SOURCES = ../common/Rlapy3.test.cpp $(test_srcdepends)
Rlapy3_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlarf_test_mpfr_SOURCES = ../common/Rlarf.test.cpp $(test_srcdepends)
Rlarf_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlarfb_test_mpfr_SOURCES = ../common/Rlarfb.test.cpp $(test_srcdepends)
Rlarfb_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlarfg_test_mpfr_SOURCES = ../common/Rlarfg.test.cpp $(test_srcdepends)
Rlarfg_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlarft_test_mpfr_SOURCES = ../common/Rlarft.test.cpp $(test_srcdepends)
Rlarft_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlartg_test_mpfr_SOURCES = ../common/Rlartg.test.cpp $(test_srcdepends)
Rlartg_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlaruv_test_mpfr_SOURCES = ../common/Rlaruv.test.cpp $(test_srcdepends)
Rlaruv_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlascl_test_mpfr_SOURCES = ../common/Rlascl.test.cpp $(test_srcdepends)
Rlascl_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlaset_test_mpfr_SOURCES = ../common/Rlaset.test.cpp $(test_srcdepends)
Rlaset_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlasr_test_mpfr_SOURCES = ../common/Rlasr.test.cpp $(test_srcdepends)
Rlasr_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlasrt_test_mpfr_SOURCES = ../common/Rlasrt.test.cpp $(test_srcdepends)
Rlasrt_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlassq_test_mpfr_SOURCES = ../common/Rlassq.test.cpp $(test_srcdepends)
Rlassq_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlaswp_test_mpfr_SOURCES = ../common/Rlaswp.test.cpp $(test_srcdepends)
Rlaswp_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlasyf_test_mpfr_SOURCES = ../common/Rlasyf.test.cpp $(test_srcdepends)
Rlasyf_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlatrd_test_mpfr_SOURCES = ../common/Rlatrd.test.cpp $(test_srcdepends)
Rlatrd_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlatrs_test_mpfr_SOURCES = ../common/Rlatrs.test.cpp $(test_srcdepends)
Rlatrs_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlauu2_test_mpfr_SOURCES = ../common/Rlauu2.test.cpp $(test_srcdepends)
Rlauu2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rlauum_test_mpfr_SOURCES = ../common/Rlauum.test.cpp $(test_srcdepends)
Rlauum_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rnaninf_test_mpfr_SOURCES = ../common/Rnaninf.test.cpp $(test_srcdepends)
Rnaninf_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rorg2l_test_mpfr_SOURCES = ../common/Rorg2l.test.cpp $(test_srcdepends)
Rorg2l_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rorg2r_test_mpfr_SOURCES = ../common/Rorg2r.test.cpp $(test_srcdepends)
Rorg2r_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rorgql_test_mpfr_SOURCES = ../common/Rorgql.test.cpp $(test_srcdepends)
Rorgql_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rorgqr_test_mpfr_SOURCES = ../common/Rorgqr.test.cpp $(test_srcdepends)
Rorgqr_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rorgtr_test_mpfr_SOURCES = ../common/Rorgtr.test.cpp $(test_srcdepends)
Rorgtr_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rpocon_test_mpfr_SOURCES = ../common/Rpocon.test.cpp $(test_srcdepends)
Rpocon_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rposv_test_mpfr_SOURCES = ../common/Rposv.test.cpp $(test_srcdepends)
Rposv_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rpotf2_test_mpfr_SOURCES = ../common/Rpotf2.test.cpp $(test_srcdepends)
Rpotf2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rpotri_test_mpfr_SOURCES = ../common/Rpotri.test.cpp $(test_srcdepends)
Rpotri_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rpotrs_test_mpfr_SOURCES = ../common/Rpotrs.test.cpp $(test_srcdepends)
Rpotrs_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rsteqr_test_mpfr_SOURCES = ../common/Rsteqr.test.cpp $(test_srcdepends)
Rsteqr_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rsterf_test_mpfr_SOURCES = ../common/Rsterf.test.cpp $(test_srcdepends)
Rsterf_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rsyev_test_mpfr_SOURCES = ../common/Rsyev.test.cpp $(test_srcdepends)
Rsyev_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rsytd2_test_mpfr_SOURCES = ../common/Rsytd2.test.cpp $(test_srcdepends)
Rsytd2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rsytrd_test_mpfr_SOURCES = ../common/Rsytrd.test.cpp $(test_srcdepends)
Rsytrd_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rtrti2_test_mpfr_SOURCES = ../common/Rtrti2.test.cpp $(test_srcdepends)
Rtrti2_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rtrtri_test_mpfr_SOURCES = ../common/Rtrtri.test.cpp $(test_srcdepends)
Rtrtri_test_mpfr_CXXFLAGS = $(test_cxxflags)
Rtrtrs_test_mpfr_SOURCES = ../common/Rtrtrs.test.cpp $(test_srcdepends)
Rtrtrs_test_mpfr_CXXFLAGS = $(test_cxxflags)
iCmax1_test_mpfr_SOURCES = ../common/iCmax1.test.cpp $(test_srcdepends)
iCmax1_test_mpfr_CXXFLAGS = $(test_cxxflags)
iMlaenv_test_mpfr_SOURCES = ../common/iMlaenv.test.cpp $(test_srcdepends)
iMlaenv_test_mpfr_CXXFLAGS = $(test_cxxflags)

include $(top_srcdir)/mplapack/test/compare/Makefile.backend.am
include $(top_srcdir)/mplapack/test/compare/Makefile.tests.am
