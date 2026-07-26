run_Rlaruv_test.sh: Makefile
	@echo '#!/bin/sh' > $@
	@echo '$${LOG_COMPILER} ./Rlaruv.test_binary128$${EXEEXT} || true' >> $@
	@echo 'if diff -u $${srcdir}/Rlaruv_reference.txt Rlaruv.txt > /dev/null 2>&1; then' >> $@
	@echo '  echo "PASS: Rlaruv"; exit 0' >> $@
	@echo 'else' >> $@
	@echo '  echo "FAIL: Rlaruv"; diff $${srcdir}/Rlaruv_reference.txt Rlaruv.txt; exit 1' >> $@
	@echo 'fi' >> $@
	@chmod +x $@

run_Rlamch_test.sh: Makefile
	@echo '#!/bin/sh' > $@
	@echo '$${LOG_COMPILER} ./Rlamch.test_binary128$${EXEEXT} || true' >> $@
	@echo 'if diff -u $${srcdir}/Rlamch_reference.txt Rlamch.txt > /dev/null 2>&1; then' >> $@
	@echo '  echo "PASS: Rlamch"; exit 0' >> $@
	@echo 'else' >> $@
	@echo '  echo "FAIL: Rlamch"; diff $${srcdir}/Rlamch_reference.txt Rlamch.txt; exit 1' >> $@
	@echo 'fi' >> $@
	@chmod +x $@

EXTRA_DIST = Rlaruv_reference.txt Rlamch_reference.txt

backend_define = -D___MPLAPACK_BUILD_WITH_BINARY128___
backend_cppflags = -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) $(backend_define) -D___MPLAPACK_INTERNAL___
backend_cxxflags =
backend_reference_archives = $(top_builddir)/mplapack/reference/.libs/libmplapack_binary128.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_mpfr.a
backend_reference_libs = -lmplapack_binary128 -lmplapack_mpfr
backend_whole_archive_libs = -lquadmath
backend_link_libs =
backend_testdir = $(prefix)/lib/$(target)/mplapack/test/compare/binary128
mplibs = -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -lmpc -lmpfr -lgmp

compare_tests = \
    Cgelq2.test_binary128 \
    Cgeqr2.test_binary128 \
    Cgeqrf.test_binary128 \
    Cgesv.test_binary128 \
    Cgetf2.test_binary128 \
    Cgetrf.test_binary128 \
    Cgetri.test_binary128 \
    Cgetrs.test_binary128 \
    Cheev.test_binary128 \
    Chetd2.test_binary128 \
    Chetrd.test_binary128 \
    Clacgv.test_binary128 \
    Clacrm.test_binary128 \
    Clacrt.test_binary128 \
    Cladiv.test_binary128 \
    Claesy.test_binary128 \
    Claev2.test_binary128 \
    Clahef.test_binary128 \
    Clanhe.test_binary128 \
    Clanht.test_binary128 \
    Clansy.test_binary128 \
    Clarf.test_binary128 \
    Clarfb.test_binary128 \
    Clarfg.test_binary128 \
    Clarft.test_binary128 \
    Clartg.test_binary128 \
    Clascl.test_binary128 \
    Claset.test_binary128 \
    Clasr.test_binary128 \
    Classq.test_binary128 \
    Claswp.test_binary128 \
    Clasyf.test_binary128 \
    Clatrd.test_binary128 \
    Cpotf2.test_binary128 \
    Crot.test_binary128 \
    Cspmv.test_binary128 \
    Cspr.test_binary128 \
    Csteqr.test_binary128 \
    Csymv.test_binary128 \
    Csyr.test_binary128 \
    Ctrti2.test_binary128 \
    Ctrtri.test_binary128 \
    Ctrtrs.test_binary128 \
    Cung2l.test_binary128 \
    Cung2r.test_binary128 \
    Cungql.test_binary128 \
    Cungqr.test_binary128 \
    Cungtr.test_binary128 \
    Cunmqr.test_binary128 \
    Mlsamen.test_binary128 \
    Mutils.test_binary128 \
    RCsum1.test_binary128 \
    Rgecon.test_binary128 \
    Rgeequ.test_binary128 \
    Rgeqlf.test_binary128 \
    Rgesv.test_binary128 \
    Rgetf2.test_binary128 \
    Rgetrf.test_binary128 \
    Rgetri.test_binary128 \
    Rgetrs.test_binary128 \
    Rladiv.test_binary128 \
    Rlae2.test_binary128 \
    Rlaev2.test_binary128 \
    Rlamch.test_binary128 \
    Rlange.test_binary128 \
    Rlanst.test_binary128 \
    Rlansy.test_binary128 \
    Rlapy2.test_binary128 \
    Rlapy3.test_binary128 \
    Rlarf.test_binary128 \
    Rlarfb.test_binary128 \
    Rlarfg.test_binary128 \
    Rlarft.test_binary128 \
    Rlartg.test_binary128 \
    Rlaruv.test_binary128 \
    Rlascl.test_binary128 \
    Rlaset.test_binary128 \
    Rlasr.test_binary128 \
    Rlasrt.test_binary128 \
    Rlassq.test_binary128 \
    Rlaswp.test_binary128 \
    Rlasyf.test_binary128 \
    Rlatrd.test_binary128 \
    Rlatrs.test_binary128 \
    Rlauu2.test_binary128 \
    Rlauum.test_binary128 \
    Rnaninf.test_binary128 \
    Rorg2l.test_binary128 \
    Rorg2r.test_binary128 \
    Rorgql.test_binary128 \
    Rorgqr.test_binary128 \
    Rorgtr.test_binary128 \
    Rpocon.test_binary128 \
    Rposv.test_binary128 \
    Rpotf2.test_binary128 \
    Rpotri.test_binary128 \
    Rpotrs.test_binary128 \
    Rsteqr.test_binary128 \
    Rsterf.test_binary128 \
    Rsyev.test_binary128 \
    Rsytd2.test_binary128 \
    Rsytrd.test_binary128 \
    Rtrti2.test_binary128 \
    Rtrtri.test_binary128 \
    Rtrtrs.test_binary128 \
    iCmax1.test_binary128 \
    iMlaenv.test_binary128

Cgelq2_test_binary128_SOURCES = ../common/Cgelq2.test.cpp $(test_srcdepends)
Cgelq2_test_binary128_CXXFLAGS = $(test_cxxflags)
Cgeqr2_test_binary128_SOURCES = ../common/Cgeqr2.test.cpp $(test_srcdepends)
Cgeqr2_test_binary128_CXXFLAGS = $(test_cxxflags)
Cgeqrf_test_binary128_SOURCES = ../common/Cgeqrf.test.cpp $(test_srcdepends)
Cgeqrf_test_binary128_CXXFLAGS = $(test_cxxflags)
Cgesv_test_binary128_SOURCES = ../common/Cgesv.test.cpp $(test_srcdepends)
Cgesv_test_binary128_CXXFLAGS = $(test_cxxflags)
Cgetf2_test_binary128_SOURCES = ../common/Cgetf2.test.cpp $(test_srcdepends)
Cgetf2_test_binary128_CXXFLAGS = $(test_cxxflags)
Cgetrf_test_binary128_SOURCES = ../common/Cgetrf.test.cpp $(test_srcdepends)
Cgetrf_test_binary128_CXXFLAGS = $(test_cxxflags)
Cgetri_test_binary128_SOURCES = ../common/Cgetri.test.cpp $(test_srcdepends)
Cgetri_test_binary128_CXXFLAGS = $(test_cxxflags)
Cgetrs_test_binary128_SOURCES = ../common/Cgetrs.test.cpp $(test_srcdepends)
Cgetrs_test_binary128_CXXFLAGS = $(test_cxxflags)
Cheev_test_binary128_SOURCES = ../common/Cheev.test.cpp $(test_srcdepends)
Cheev_test_binary128_CXXFLAGS = $(test_cxxflags)
Chetd2_test_binary128_SOURCES = ../common/Chetd2.test.cpp $(test_srcdepends)
Chetd2_test_binary128_CXXFLAGS = $(test_cxxflags)
Chetrd_test_binary128_SOURCES = ../common/Chetrd.test.cpp $(test_srcdepends)
Chetrd_test_binary128_CXXFLAGS = $(test_cxxflags)
Clacgv_test_binary128_SOURCES = ../common/Clacgv.test.cpp $(test_srcdepends)
Clacgv_test_binary128_CXXFLAGS = $(test_cxxflags)
Clacrm_test_binary128_SOURCES = ../common/Clacrm.test.cpp $(test_srcdepends)
Clacrm_test_binary128_CXXFLAGS = $(test_cxxflags)
Clacrt_test_binary128_SOURCES = ../common/Clacrt.test.cpp $(test_srcdepends)
Clacrt_test_binary128_CXXFLAGS = $(test_cxxflags)
Cladiv_test_binary128_SOURCES = ../common/Cladiv.test.cpp $(test_srcdepends)
Cladiv_test_binary128_CXXFLAGS = $(test_cxxflags)
Claesy_test_binary128_SOURCES = ../common/Claesy.test.cpp $(test_srcdepends)
Claesy_test_binary128_CXXFLAGS = $(test_cxxflags)
Claev2_test_binary128_SOURCES = ../common/Claev2.test.cpp $(test_srcdepends)
Claev2_test_binary128_CXXFLAGS = $(test_cxxflags)
Clahef_test_binary128_SOURCES = ../common/Clahef.test.cpp $(test_srcdepends)
Clahef_test_binary128_CXXFLAGS = $(test_cxxflags)
Clanhe_test_binary128_SOURCES = ../common/Clanhe.test.cpp $(test_srcdepends)
Clanhe_test_binary128_CXXFLAGS = $(test_cxxflags)
Clanht_test_binary128_SOURCES = ../common/Clanht.test.cpp $(test_srcdepends)
Clanht_test_binary128_CXXFLAGS = $(test_cxxflags)
Clansy_test_binary128_SOURCES = ../common/Clansy.test.cpp $(test_srcdepends)
Clansy_test_binary128_CXXFLAGS = $(test_cxxflags)
Clarf_test_binary128_SOURCES = ../common/Clarf.test.cpp $(test_srcdepends)
Clarf_test_binary128_CXXFLAGS = $(test_cxxflags)
Clarfb_test_binary128_SOURCES = ../common/Clarfb.test.cpp $(test_srcdepends)
Clarfb_test_binary128_CXXFLAGS = $(test_cxxflags)
Clarfg_test_binary128_SOURCES = ../common/Clarfg.test.cpp $(test_srcdepends)
Clarfg_test_binary128_CXXFLAGS = $(test_cxxflags)
Clarft_test_binary128_SOURCES = ../common/Clarft.test.cpp $(test_srcdepends)
Clarft_test_binary128_CXXFLAGS = $(test_cxxflags)
Clartg_test_binary128_SOURCES = ../common/Clartg.test.cpp $(test_srcdepends)
Clartg_test_binary128_CXXFLAGS = $(test_cxxflags)
Clascl_test_binary128_SOURCES = ../common/Clascl.test.cpp $(test_srcdepends)
Clascl_test_binary128_CXXFLAGS = $(test_cxxflags)
Claset_test_binary128_SOURCES = ../common/Claset.test.cpp $(test_srcdepends)
Claset_test_binary128_CXXFLAGS = $(test_cxxflags)
Clasr_test_binary128_SOURCES = ../common/Clasr.test.cpp $(test_srcdepends)
Clasr_test_binary128_CXXFLAGS = $(test_cxxflags)
Classq_test_binary128_SOURCES = ../common/Classq.test.cpp $(test_srcdepends)
Classq_test_binary128_CXXFLAGS = $(test_cxxflags)
Claswp_test_binary128_SOURCES = ../common/Claswp.test.cpp $(test_srcdepends)
Claswp_test_binary128_CXXFLAGS = $(test_cxxflags)
Clasyf_test_binary128_SOURCES = ../common/Clasyf.test.cpp $(test_srcdepends)
Clasyf_test_binary128_CXXFLAGS = $(test_cxxflags)
Clatrd_test_binary128_SOURCES = ../common/Clatrd.test.cpp $(test_srcdepends)
Clatrd_test_binary128_CXXFLAGS = $(test_cxxflags)
Cpotf2_test_binary128_SOURCES = ../common/Cpotf2.test.cpp $(test_srcdepends)
Cpotf2_test_binary128_CXXFLAGS = $(test_cxxflags)
Crot_test_binary128_SOURCES = ../common/Crot.test.cpp $(test_srcdepends)
Crot_test_binary128_CXXFLAGS = $(test_cxxflags)
Cspmv_test_binary128_SOURCES = ../common/Cspmv.test.cpp $(test_srcdepends)
Cspmv_test_binary128_CXXFLAGS = $(test_cxxflags)
Cspr_test_binary128_SOURCES = ../common/Cspr.test.cpp $(test_srcdepends)
Cspr_test_binary128_CXXFLAGS = $(test_cxxflags)
Csteqr_test_binary128_SOURCES = ../common/Csteqr.test.cpp $(test_srcdepends)
Csteqr_test_binary128_CXXFLAGS = $(test_cxxflags)
Csymv_test_binary128_SOURCES = ../common/Csymv.test.cpp $(test_srcdepends)
Csymv_test_binary128_CXXFLAGS = $(test_cxxflags)
Csyr_test_binary128_SOURCES = ../common/Csyr.test.cpp $(test_srcdepends)
Csyr_test_binary128_CXXFLAGS = $(test_cxxflags)
Ctrti2_test_binary128_SOURCES = ../common/Ctrti2.test.cpp $(test_srcdepends)
Ctrti2_test_binary128_CXXFLAGS = $(test_cxxflags)
Ctrtri_test_binary128_SOURCES = ../common/Ctrtri.test.cpp $(test_srcdepends)
Ctrtri_test_binary128_CXXFLAGS = $(test_cxxflags)
Ctrtrs_test_binary128_SOURCES = ../common/Ctrtrs.test.cpp $(test_srcdepends)
Ctrtrs_test_binary128_CXXFLAGS = $(test_cxxflags)
Cung2l_test_binary128_SOURCES = ../common/Cung2l.test.cpp $(test_srcdepends)
Cung2l_test_binary128_CXXFLAGS = $(test_cxxflags)
Cung2r_test_binary128_SOURCES = ../common/Cung2r.test.cpp $(test_srcdepends)
Cung2r_test_binary128_CXXFLAGS = $(test_cxxflags)
Cungql_test_binary128_SOURCES = ../common/Cungql.test.cpp $(test_srcdepends)
Cungql_test_binary128_CXXFLAGS = $(test_cxxflags)
Cungqr_test_binary128_SOURCES = ../common/Cungqr.test.cpp $(test_srcdepends)
Cungqr_test_binary128_CXXFLAGS = $(test_cxxflags)
Cungtr_test_binary128_SOURCES = ../common/Cungtr.test.cpp $(test_srcdepends)
Cungtr_test_binary128_CXXFLAGS = $(test_cxxflags)
Cunmqr_test_binary128_SOURCES = ../common/Cunmqr.test.cpp $(test_srcdepends)
Cunmqr_test_binary128_CXXFLAGS = $(test_cxxflags)
Mlsamen_test_binary128_SOURCES = ../common/Mlsamen.test.cpp $(test_srcdepends)
Mlsamen_test_binary128_CXXFLAGS = $(test_cxxflags)
Mutils_test_binary128_SOURCES = ../common/Mutils.test.cpp $(test_srcdepends) ../common/ilaenv.cpp
Mutils_test_binary128_CXXFLAGS = $(test_cxxflags)
RCsum1_test_binary128_SOURCES = ../common/RCsum1.test.cpp $(test_srcdepends)
RCsum1_test_binary128_CXXFLAGS = $(test_cxxflags)
Rgecon_test_binary128_SOURCES = ../common/Rgecon.test.cpp $(test_srcdepends)
Rgecon_test_binary128_CXXFLAGS = $(test_cxxflags)
Rgeequ_test_binary128_SOURCES = ../common/Rgeequ.test.cpp $(test_srcdepends)
Rgeequ_test_binary128_CXXFLAGS = $(test_cxxflags)
Rgeqlf_test_binary128_SOURCES = ../common/Rgeqlf.test.cpp $(test_srcdepends)
Rgeqlf_test_binary128_CXXFLAGS = $(test_cxxflags)
Rgesv_test_binary128_SOURCES = ../common/Rgesv.test.cpp $(test_srcdepends)
Rgesv_test_binary128_CXXFLAGS = $(test_cxxflags)
Rgetf2_test_binary128_SOURCES = ../common/Rgetf2.test.cpp $(test_srcdepends)
Rgetf2_test_binary128_CXXFLAGS = $(test_cxxflags)
Rgetrf_test_binary128_SOURCES = ../common/Rgetrf.test.cpp $(test_srcdepends)
Rgetrf_test_binary128_CXXFLAGS = $(test_cxxflags)
Rgetri_test_binary128_SOURCES = ../common/Rgetri.test.cpp $(test_srcdepends)
Rgetri_test_binary128_CXXFLAGS = $(test_cxxflags)
Rgetrs_test_binary128_SOURCES = ../common/Rgetrs.test.cpp $(test_srcdepends)
Rgetrs_test_binary128_CXXFLAGS = $(test_cxxflags)
Rladiv_test_binary128_SOURCES = ../common/Rladiv.test.cpp $(test_srcdepends)
Rladiv_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlae2_test_binary128_SOURCES = ../common/Rlae2.test.cpp $(test_srcdepends)
Rlae2_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlaev2_test_binary128_SOURCES = ../common/Rlaev2.test.cpp $(test_srcdepends)
Rlaev2_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlamch_test_binary128_SOURCES = ../common/Rlamch.test.cpp $(test_srcdepends)
Rlamch_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlange_test_binary128_SOURCES = ../common/Rlange.test.cpp $(test_srcdepends)
Rlange_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlanst_test_binary128_SOURCES = ../common/Rlanst.test.cpp $(test_srcdepends)
Rlanst_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlansy_test_binary128_SOURCES = ../common/Rlansy.test.cpp $(test_srcdepends)
Rlansy_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlapy2_test_binary128_SOURCES = ../common/Rlapy2.test.cpp $(test_srcdepends)
Rlapy2_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlapy3_test_binary128_SOURCES = ../common/Rlapy3.test.cpp $(test_srcdepends)
Rlapy3_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlarf_test_binary128_SOURCES = ../common/Rlarf.test.cpp $(test_srcdepends)
Rlarf_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlarfb_test_binary128_SOURCES = ../common/Rlarfb.test.cpp $(test_srcdepends)
Rlarfb_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlarfg_test_binary128_SOURCES = ../common/Rlarfg.test.cpp $(test_srcdepends)
Rlarfg_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlarft_test_binary128_SOURCES = ../common/Rlarft.test.cpp $(test_srcdepends)
Rlarft_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlartg_test_binary128_SOURCES = ../common/Rlartg.test.cpp $(test_srcdepends)
Rlartg_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlaruv_test_binary128_SOURCES = ../common/Rlaruv.test.cpp $(test_srcdepends)
Rlaruv_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlascl_test_binary128_SOURCES = ../common/Rlascl.test.cpp $(test_srcdepends)
Rlascl_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlaset_test_binary128_SOURCES = ../common/Rlaset.test.cpp $(test_srcdepends)
Rlaset_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlasr_test_binary128_SOURCES = ../common/Rlasr.test.cpp $(test_srcdepends)
Rlasr_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlasrt_test_binary128_SOURCES = ../common/Rlasrt.test.cpp $(test_srcdepends)
Rlasrt_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlassq_test_binary128_SOURCES = ../common/Rlassq.test.cpp $(test_srcdepends)
Rlassq_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlaswp_test_binary128_SOURCES = ../common/Rlaswp.test.cpp $(test_srcdepends)
Rlaswp_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlasyf_test_binary128_SOURCES = ../common/Rlasyf.test.cpp $(test_srcdepends)
Rlasyf_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlatrd_test_binary128_SOURCES = ../common/Rlatrd.test.cpp $(test_srcdepends)
Rlatrd_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlatrs_test_binary128_SOURCES = ../common/Rlatrs.test.cpp $(test_srcdepends)
Rlatrs_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlauu2_test_binary128_SOURCES = ../common/Rlauu2.test.cpp $(test_srcdepends)
Rlauu2_test_binary128_CXXFLAGS = $(test_cxxflags)
Rlauum_test_binary128_SOURCES = ../common/Rlauum.test.cpp $(test_srcdepends)
Rlauum_test_binary128_CXXFLAGS = $(test_cxxflags)
Rnaninf_test_binary128_SOURCES = ../common/Rnaninf.test.cpp $(test_srcdepends)
Rnaninf_test_binary128_CXXFLAGS = $(test_cxxflags)
Rorg2l_test_binary128_SOURCES = ../common/Rorg2l.test.cpp $(test_srcdepends)
Rorg2l_test_binary128_CXXFLAGS = $(test_cxxflags)
Rorg2r_test_binary128_SOURCES = ../common/Rorg2r.test.cpp $(test_srcdepends)
Rorg2r_test_binary128_CXXFLAGS = $(test_cxxflags)
Rorgql_test_binary128_SOURCES = ../common/Rorgql.test.cpp $(test_srcdepends)
Rorgql_test_binary128_CXXFLAGS = $(test_cxxflags)
Rorgqr_test_binary128_SOURCES = ../common/Rorgqr.test.cpp $(test_srcdepends)
Rorgqr_test_binary128_CXXFLAGS = $(test_cxxflags)
Rorgtr_test_binary128_SOURCES = ../common/Rorgtr.test.cpp $(test_srcdepends)
Rorgtr_test_binary128_CXXFLAGS = $(test_cxxflags)
Rpocon_test_binary128_SOURCES = ../common/Rpocon.test.cpp $(test_srcdepends)
Rpocon_test_binary128_CXXFLAGS = $(test_cxxflags)
Rposv_test_binary128_SOURCES = ../common/Rposv.test.cpp $(test_srcdepends)
Rposv_test_binary128_CXXFLAGS = $(test_cxxflags)
Rpotf2_test_binary128_SOURCES = ../common/Rpotf2.test.cpp $(test_srcdepends)
Rpotf2_test_binary128_CXXFLAGS = $(test_cxxflags)
Rpotri_test_binary128_SOURCES = ../common/Rpotri.test.cpp $(test_srcdepends)
Rpotri_test_binary128_CXXFLAGS = $(test_cxxflags)
Rpotrs_test_binary128_SOURCES = ../common/Rpotrs.test.cpp $(test_srcdepends)
Rpotrs_test_binary128_CXXFLAGS = $(test_cxxflags)
Rsteqr_test_binary128_SOURCES = ../common/Rsteqr.test.cpp $(test_srcdepends)
Rsteqr_test_binary128_CXXFLAGS = $(test_cxxflags)
Rsterf_test_binary128_SOURCES = ../common/Rsterf.test.cpp $(test_srcdepends)
Rsterf_test_binary128_CXXFLAGS = $(test_cxxflags)
Rsyev_test_binary128_SOURCES = ../common/Rsyev.test.cpp $(test_srcdepends)
Rsyev_test_binary128_CXXFLAGS = $(test_cxxflags)
Rsytd2_test_binary128_SOURCES = ../common/Rsytd2.test.cpp $(test_srcdepends)
Rsytd2_test_binary128_CXXFLAGS = $(test_cxxflags)
Rsytrd_test_binary128_SOURCES = ../common/Rsytrd.test.cpp $(test_srcdepends)
Rsytrd_test_binary128_CXXFLAGS = $(test_cxxflags)
Rtrti2_test_binary128_SOURCES = ../common/Rtrti2.test.cpp $(test_srcdepends)
Rtrti2_test_binary128_CXXFLAGS = $(test_cxxflags)
Rtrtri_test_binary128_SOURCES = ../common/Rtrtri.test.cpp $(test_srcdepends)
Rtrtri_test_binary128_CXXFLAGS = $(test_cxxflags)
Rtrtrs_test_binary128_SOURCES = ../common/Rtrtrs.test.cpp $(test_srcdepends)
Rtrtrs_test_binary128_CXXFLAGS = $(test_cxxflags)
iCmax1_test_binary128_SOURCES = ../common/iCmax1.test.cpp $(test_srcdepends)
iCmax1_test_binary128_CXXFLAGS = $(test_cxxflags)
iMlaenv_test_binary128_SOURCES = ../common/iMlaenv.test.cpp $(test_srcdepends)
iMlaenv_test_binary128_CXXFLAGS = $(test_cxxflags)

include $(top_srcdir)/mplapack/test/compare/Makefile.backend.am
include $(top_srcdir)/mplapack/test/compare/Makefile.tests.am
