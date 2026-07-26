run_Rlaruv_test.sh: Makefile
	@echo '#!/bin/sh' > $@
	@echo 'NUM_DIFF_PY="$(top_srcdir)/misc/num_diff.py"' >> $@
	@echo '$${LOG_COMPILER} ./Rlaruv.test_qd$${EXEEXT} || true' >> $@
	@echo 'if python3 "$${NUM_DIFF_PY}" "$${srcdir}/Rlaruv_reference.txt" Rlaruv.txt --tol 1e-60; then' >> $@
	@echo '  echo "PASS: Rlaruv"; exit 0' >> $@
	@echo 'else' >> $@
	@echo '  echo "FAIL: Rlaruv"; exit 1' >> $@
	@echo 'fi' >> $@
	@chmod +x $@

run_Rlamch_test.sh: Makefile
	@echo '#!/bin/sh' > $@
	@echo 'NUM_DIFF_PY="$(top_srcdir)/misc/num_diff.py"' >> $@
	@echo '$${LOG_COMPILER} ./Rlamch.test_qd$${EXEEXT} || true' >> $@
	@echo 'if python3 "$${NUM_DIFF_PY}" "$${srcdir}/Rlamch_reference.txt" Rlamch.txt --tol 1e-60; then' >> $@
	@echo '  echo "PASS: Rlamch"; exit 0' >> $@
	@echo 'else' >> $@
	@echo '  echo "FAIL: Rlamch"; exit 1' >> $@
	@echo 'fi' >> $@
	@chmod +x $@

EXTRA_DIST = Rlaruv_reference.txt Rlamch_reference.txt

backend_define = -D___MPLAPACK_BUILD_WITH_QD___
backend_cppflags = -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(QD_INCLUDEDIR) $(backend_define) -D___MPLAPACK_INTERNAL___
backend_cxxflags =
backend_reference_archives = $(top_builddir)/mplapack/reference/.libs/libmplapack_qd.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_mpfr.a
backend_reference_libs = -lmplapack_qd -lmplapack_mpfr
backend_whole_archive_libs =
backend_link_libs =
backend_testdir = $(prefix)/lib/$(target)/mplapack/test/compare/qd
mplibs = -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(QD_LIBDIR) -lmpc -lmpfr -lgmp -lqd

compare_tests = \
    Cgelq2.test_qd \
    Cgeqr2.test_qd \
    Cgeqrf.test_qd \
    Cgesv.test_qd \
    Cgetf2.test_qd \
    Cgetrf.test_qd \
    Cgetri.test_qd \
    Cgetrs.test_qd \
    Cheev.test_qd \
    Chetd2.test_qd \
    Chetrd.test_qd \
    Clacgv.test_qd \
    Clacrm.test_qd \
    Clacrt.test_qd \
    Cladiv.test_qd \
    Claesy.test_qd \
    Claev2.test_qd \
    Clahef.test_qd \
    Clanhe.test_qd \
    Clanht.test_qd \
    Clansy.test_qd \
    Clarf.test_qd \
    Clarfb.test_qd \
    Clarfg.test_qd \
    Clarft.test_qd \
    Clartg.test_qd \
    Clascl.test_qd \
    Claset.test_qd \
    Clasr.test_qd \
    Classq.test_qd \
    Claswp.test_qd \
    Clasyf.test_qd \
    Clatrd.test_qd \
    Cpotf2.test_qd \
    Crot.test_qd \
    Cspmv.test_qd \
    Cspr.test_qd \
    Csteqr.test_qd \
    Csymv.test_qd \
    Csyr.test_qd \
    Ctrti2.test_qd \
    Ctrtri.test_qd \
    Ctrtrs.test_qd \
    Cung2l.test_qd \
    Cung2r.test_qd \
    Cungql.test_qd \
    Cungqr.test_qd \
    Cungtr.test_qd \
    Cunmqr.test_qd \
    Mlsamen.test_qd \
    Mutils.test_qd \
    RCsum1.test_qd \
    Rgecon.test_qd \
    Rgeequ.test_qd \
    Rgeqlf.test_qd \
    Rgesv.test_qd \
    Rgetf2.test_qd \
    Rgetrf.test_qd \
    Rgetri.test_qd \
    Rgetrs.test_qd \
    Rladiv.test_qd \
    Rlae2.test_qd \
    Rlaev2.test_qd \
    Rlamch.test_qd \
    Rlange.test_qd \
    Rlanst.test_qd \
    Rlansy.test_qd \
    Rlapy2.test_qd \
    Rlapy3.test_qd \
    Rlarf.test_qd \
    Rlarfb.test_qd \
    Rlarfg.test_qd \
    Rlarft.test_qd \
    Rlartg.test_qd \
    Rlaruv.test_qd \
    Rlascl.test_qd \
    Rlaset.test_qd \
    Rlasr.test_qd \
    Rlasrt.test_qd \
    Rlassq.test_qd \
    Rlaswp.test_qd \
    Rlasyf.test_qd \
    Rlatrd.test_qd \
    Rlatrs.test_qd \
    Rlauu2.test_qd \
    Rlauum.test_qd \
    Rnaninf.test_qd \
    Rorg2l.test_qd \
    Rorg2r.test_qd \
    Rorgql.test_qd \
    Rorgqr.test_qd \
    Rorgtr.test_qd \
    Rpocon.test_qd \
    Rposv.test_qd \
    Rpotf2.test_qd \
    Rpotri.test_qd \
    Rpotrs.test_qd \
    Rsteqr.test_qd \
    Rsterf.test_qd \
    Rsyev.test_qd \
    Rsytd2.test_qd \
    Rsytrd.test_qd \
    Rtrti2.test_qd \
    Rtrtri.test_qd \
    Rtrtrs.test_qd \
    iCmax1.test_qd \
    iMlaenv.test_qd

Cgelq2_test_qd_SOURCES = ../common/Cgelq2.test.cpp $(test_srcdepends)
Cgelq2_test_qd_CXXFLAGS = $(test_cxxflags)
Cgeqr2_test_qd_SOURCES = ../common/Cgeqr2.test.cpp $(test_srcdepends)
Cgeqr2_test_qd_CXXFLAGS = $(test_cxxflags)
Cgeqrf_test_qd_SOURCES = ../common/Cgeqrf.test.cpp $(test_srcdepends)
Cgeqrf_test_qd_CXXFLAGS = $(test_cxxflags)
Cgesv_test_qd_SOURCES = ../common/Cgesv.test.cpp $(test_srcdepends)
Cgesv_test_qd_CXXFLAGS = $(test_cxxflags)
Cgetf2_test_qd_SOURCES = ../common/Cgetf2.test.cpp $(test_srcdepends)
Cgetf2_test_qd_CXXFLAGS = $(test_cxxflags)
Cgetrf_test_qd_SOURCES = ../common/Cgetrf.test.cpp $(test_srcdepends)
Cgetrf_test_qd_CXXFLAGS = $(test_cxxflags)
Cgetri_test_qd_SOURCES = ../common/Cgetri.test.cpp $(test_srcdepends)
Cgetri_test_qd_CXXFLAGS = $(test_cxxflags)
Cgetrs_test_qd_SOURCES = ../common/Cgetrs.test.cpp $(test_srcdepends)
Cgetrs_test_qd_CXXFLAGS = $(test_cxxflags)
Cheev_test_qd_SOURCES = ../common/Cheev.test.cpp $(test_srcdepends)
Cheev_test_qd_CXXFLAGS = $(test_cxxflags)
Chetd2_test_qd_SOURCES = ../common/Chetd2.test.cpp $(test_srcdepends)
Chetd2_test_qd_CXXFLAGS = $(test_cxxflags)
Chetrd_test_qd_SOURCES = ../common/Chetrd.test.cpp $(test_srcdepends)
Chetrd_test_qd_CXXFLAGS = $(test_cxxflags)
Clacgv_test_qd_SOURCES = ../common/Clacgv.test.cpp $(test_srcdepends)
Clacgv_test_qd_CXXFLAGS = $(test_cxxflags)
Clacrm_test_qd_SOURCES = ../common/Clacrm.test.cpp $(test_srcdepends)
Clacrm_test_qd_CXXFLAGS = $(test_cxxflags)
Clacrt_test_qd_SOURCES = ../common/Clacrt.test.cpp $(test_srcdepends)
Clacrt_test_qd_CXXFLAGS = $(test_cxxflags)
Cladiv_test_qd_SOURCES = ../common/Cladiv.test.cpp $(test_srcdepends)
Cladiv_test_qd_CXXFLAGS = $(test_cxxflags)
Claesy_test_qd_SOURCES = ../common/Claesy.test.cpp $(test_srcdepends)
Claesy_test_qd_CXXFLAGS = $(test_cxxflags)
Claev2_test_qd_SOURCES = ../common/Claev2.test.cpp $(test_srcdepends)
Claev2_test_qd_CXXFLAGS = $(test_cxxflags)
Clahef_test_qd_SOURCES = ../common/Clahef.test.cpp $(test_srcdepends)
Clahef_test_qd_CXXFLAGS = $(test_cxxflags)
Clanhe_test_qd_SOURCES = ../common/Clanhe.test.cpp $(test_srcdepends)
Clanhe_test_qd_CXXFLAGS = $(test_cxxflags)
Clanht_test_qd_SOURCES = ../common/Clanht.test.cpp $(test_srcdepends)
Clanht_test_qd_CXXFLAGS = $(test_cxxflags)
Clansy_test_qd_SOURCES = ../common/Clansy.test.cpp $(test_srcdepends)
Clansy_test_qd_CXXFLAGS = $(test_cxxflags)
Clarf_test_qd_SOURCES = ../common/Clarf.test.cpp $(test_srcdepends)
Clarf_test_qd_CXXFLAGS = $(test_cxxflags)
Clarfb_test_qd_SOURCES = ../common/Clarfb.test.cpp $(test_srcdepends)
Clarfb_test_qd_CXXFLAGS = $(test_cxxflags)
Clarfg_test_qd_SOURCES = ../common/Clarfg.test.cpp $(test_srcdepends)
Clarfg_test_qd_CXXFLAGS = $(test_cxxflags)
Clarft_test_qd_SOURCES = ../common/Clarft.test.cpp $(test_srcdepends)
Clarft_test_qd_CXXFLAGS = $(test_cxxflags)
Clartg_test_qd_SOURCES = ../common/Clartg.test.cpp $(test_srcdepends)
Clartg_test_qd_CXXFLAGS = $(test_cxxflags)
Clascl_test_qd_SOURCES = ../common/Clascl.test.cpp $(test_srcdepends)
Clascl_test_qd_CXXFLAGS = $(test_cxxflags)
Claset_test_qd_SOURCES = ../common/Claset.test.cpp $(test_srcdepends)
Claset_test_qd_CXXFLAGS = $(test_cxxflags)
Clasr_test_qd_SOURCES = ../common/Clasr.test.cpp $(test_srcdepends)
Clasr_test_qd_CXXFLAGS = $(test_cxxflags)
Classq_test_qd_SOURCES = ../common/Classq.test.cpp $(test_srcdepends)
Classq_test_qd_CXXFLAGS = $(test_cxxflags)
Claswp_test_qd_SOURCES = ../common/Claswp.test.cpp $(test_srcdepends)
Claswp_test_qd_CXXFLAGS = $(test_cxxflags)
Clasyf_test_qd_SOURCES = ../common/Clasyf.test.cpp $(test_srcdepends)
Clasyf_test_qd_CXXFLAGS = $(test_cxxflags)
Clatrd_test_qd_SOURCES = ../common/Clatrd.test.cpp $(test_srcdepends)
Clatrd_test_qd_CXXFLAGS = $(test_cxxflags)
Cpotf2_test_qd_SOURCES = ../common/Cpotf2.test.cpp $(test_srcdepends)
Cpotf2_test_qd_CXXFLAGS = $(test_cxxflags)
Crot_test_qd_SOURCES = ../common/Crot.test.cpp $(test_srcdepends)
Crot_test_qd_CXXFLAGS = $(test_cxxflags)
Cspmv_test_qd_SOURCES = ../common/Cspmv.test.cpp $(test_srcdepends)
Cspmv_test_qd_CXXFLAGS = $(test_cxxflags)
Cspr_test_qd_SOURCES = ../common/Cspr.test.cpp $(test_srcdepends)
Cspr_test_qd_CXXFLAGS = $(test_cxxflags)
Csteqr_test_qd_SOURCES = ../common/Csteqr.test.cpp $(test_srcdepends)
Csteqr_test_qd_CXXFLAGS = $(test_cxxflags)
Csymv_test_qd_SOURCES = ../common/Csymv.test.cpp $(test_srcdepends)
Csymv_test_qd_CXXFLAGS = $(test_cxxflags)
Csyr_test_qd_SOURCES = ../common/Csyr.test.cpp $(test_srcdepends)
Csyr_test_qd_CXXFLAGS = $(test_cxxflags)
Ctrti2_test_qd_SOURCES = ../common/Ctrti2.test.cpp $(test_srcdepends)
Ctrti2_test_qd_CXXFLAGS = $(test_cxxflags)
Ctrtri_test_qd_SOURCES = ../common/Ctrtri.test.cpp $(test_srcdepends)
Ctrtri_test_qd_CXXFLAGS = $(test_cxxflags)
Ctrtrs_test_qd_SOURCES = ../common/Ctrtrs.test.cpp $(test_srcdepends)
Ctrtrs_test_qd_CXXFLAGS = $(test_cxxflags)
Cung2l_test_qd_SOURCES = ../common/Cung2l.test.cpp $(test_srcdepends)
Cung2l_test_qd_CXXFLAGS = $(test_cxxflags)
Cung2r_test_qd_SOURCES = ../common/Cung2r.test.cpp $(test_srcdepends)
Cung2r_test_qd_CXXFLAGS = $(test_cxxflags)
Cungql_test_qd_SOURCES = ../common/Cungql.test.cpp $(test_srcdepends)
Cungql_test_qd_CXXFLAGS = $(test_cxxflags)
Cungqr_test_qd_SOURCES = ../common/Cungqr.test.cpp $(test_srcdepends)
Cungqr_test_qd_CXXFLAGS = $(test_cxxflags)
Cungtr_test_qd_SOURCES = ../common/Cungtr.test.cpp $(test_srcdepends)
Cungtr_test_qd_CXXFLAGS = $(test_cxxflags)
Cunmqr_test_qd_SOURCES = ../common/Cunmqr.test.cpp $(test_srcdepends)
Cunmqr_test_qd_CXXFLAGS = $(test_cxxflags)
Mlsamen_test_qd_SOURCES = ../common/Mlsamen.test.cpp $(test_srcdepends)
Mlsamen_test_qd_CXXFLAGS = $(test_cxxflags)
Mutils_test_qd_SOURCES = ../common/Mutils.test.cpp $(test_srcdepends) ../common/ilaenv.cpp
Mutils_test_qd_CXXFLAGS = $(test_cxxflags)
RCsum1_test_qd_SOURCES = ../common/RCsum1.test.cpp $(test_srcdepends)
RCsum1_test_qd_CXXFLAGS = $(test_cxxflags)
Rgecon_test_qd_SOURCES = ../common/Rgecon.test.cpp $(test_srcdepends)
Rgecon_test_qd_CXXFLAGS = $(test_cxxflags)
Rgeequ_test_qd_SOURCES = ../common/Rgeequ.test.cpp $(test_srcdepends)
Rgeequ_test_qd_CXXFLAGS = $(test_cxxflags)
Rgeqlf_test_qd_SOURCES = ../common/Rgeqlf.test.cpp $(test_srcdepends)
Rgeqlf_test_qd_CXXFLAGS = $(test_cxxflags)
Rgesv_test_qd_SOURCES = ../common/Rgesv.test.cpp $(test_srcdepends)
Rgesv_test_qd_CXXFLAGS = $(test_cxxflags)
Rgetf2_test_qd_SOURCES = ../common/Rgetf2.test.cpp $(test_srcdepends)
Rgetf2_test_qd_CXXFLAGS = $(test_cxxflags)
Rgetrf_test_qd_SOURCES = ../common/Rgetrf.test.cpp $(test_srcdepends)
Rgetrf_test_qd_CXXFLAGS = $(test_cxxflags)
Rgetri_test_qd_SOURCES = ../common/Rgetri.test.cpp $(test_srcdepends)
Rgetri_test_qd_CXXFLAGS = $(test_cxxflags)
Rgetrs_test_qd_SOURCES = ../common/Rgetrs.test.cpp $(test_srcdepends)
Rgetrs_test_qd_CXXFLAGS = $(test_cxxflags)
Rladiv_test_qd_SOURCES = ../common/Rladiv.test.cpp $(test_srcdepends)
Rladiv_test_qd_CXXFLAGS = $(test_cxxflags)
Rlae2_test_qd_SOURCES = ../common/Rlae2.test.cpp $(test_srcdepends)
Rlae2_test_qd_CXXFLAGS = $(test_cxxflags)
Rlaev2_test_qd_SOURCES = ../common/Rlaev2.test.cpp $(test_srcdepends)
Rlaev2_test_qd_CXXFLAGS = $(test_cxxflags)
Rlamch_test_qd_SOURCES = ../common/Rlamch.test.cpp $(test_srcdepends)
Rlamch_test_qd_CXXFLAGS = $(test_cxxflags)
Rlange_test_qd_SOURCES = ../common/Rlange.test.cpp $(test_srcdepends)
Rlange_test_qd_CXXFLAGS = $(test_cxxflags)
Rlanst_test_qd_SOURCES = ../common/Rlanst.test.cpp $(test_srcdepends)
Rlanst_test_qd_CXXFLAGS = $(test_cxxflags)
Rlansy_test_qd_SOURCES = ../common/Rlansy.test.cpp $(test_srcdepends)
Rlansy_test_qd_CXXFLAGS = $(test_cxxflags)
Rlapy2_test_qd_SOURCES = ../common/Rlapy2.test.cpp $(test_srcdepends)
Rlapy2_test_qd_CXXFLAGS = $(test_cxxflags)
Rlapy3_test_qd_SOURCES = ../common/Rlapy3.test.cpp $(test_srcdepends)
Rlapy3_test_qd_CXXFLAGS = $(test_cxxflags)
Rlarf_test_qd_SOURCES = ../common/Rlarf.test.cpp $(test_srcdepends)
Rlarf_test_qd_CXXFLAGS = $(test_cxxflags)
Rlarfb_test_qd_SOURCES = ../common/Rlarfb.test.cpp $(test_srcdepends)
Rlarfb_test_qd_CXXFLAGS = $(test_cxxflags)
Rlarfg_test_qd_SOURCES = ../common/Rlarfg.test.cpp $(test_srcdepends)
Rlarfg_test_qd_CXXFLAGS = $(test_cxxflags)
Rlarft_test_qd_SOURCES = ../common/Rlarft.test.cpp $(test_srcdepends)
Rlarft_test_qd_CXXFLAGS = $(test_cxxflags)
Rlartg_test_qd_SOURCES = ../common/Rlartg.test.cpp $(test_srcdepends)
Rlartg_test_qd_CXXFLAGS = $(test_cxxflags)
Rlaruv_test_qd_SOURCES = ../common/Rlaruv.test.cpp $(test_srcdepends)
Rlaruv_test_qd_CXXFLAGS = $(test_cxxflags)
Rlascl_test_qd_SOURCES = ../common/Rlascl.test.cpp $(test_srcdepends)
Rlascl_test_qd_CXXFLAGS = $(test_cxxflags)
Rlaset_test_qd_SOURCES = ../common/Rlaset.test.cpp $(test_srcdepends)
Rlaset_test_qd_CXXFLAGS = $(test_cxxflags)
Rlasr_test_qd_SOURCES = ../common/Rlasr.test.cpp $(test_srcdepends)
Rlasr_test_qd_CXXFLAGS = $(test_cxxflags)
Rlasrt_test_qd_SOURCES = ../common/Rlasrt.test.cpp $(test_srcdepends)
Rlasrt_test_qd_CXXFLAGS = $(test_cxxflags)
Rlassq_test_qd_SOURCES = ../common/Rlassq.test.cpp $(test_srcdepends)
Rlassq_test_qd_CXXFLAGS = $(test_cxxflags)
Rlaswp_test_qd_SOURCES = ../common/Rlaswp.test.cpp $(test_srcdepends)
Rlaswp_test_qd_CXXFLAGS = $(test_cxxflags)
Rlasyf_test_qd_SOURCES = ../common/Rlasyf.test.cpp $(test_srcdepends)
Rlasyf_test_qd_CXXFLAGS = $(test_cxxflags)
Rlatrd_test_qd_SOURCES = ../common/Rlatrd.test.cpp $(test_srcdepends)
Rlatrd_test_qd_CXXFLAGS = $(test_cxxflags)
Rlatrs_test_qd_SOURCES = ../common/Rlatrs.test.cpp $(test_srcdepends)
Rlatrs_test_qd_CXXFLAGS = $(test_cxxflags)
Rlauu2_test_qd_SOURCES = ../common/Rlauu2.test.cpp $(test_srcdepends)
Rlauu2_test_qd_CXXFLAGS = $(test_cxxflags)
Rlauum_test_qd_SOURCES = ../common/Rlauum.test.cpp $(test_srcdepends)
Rlauum_test_qd_CXXFLAGS = $(test_cxxflags)
Rnaninf_test_qd_SOURCES = ../common/Rnaninf.test.cpp $(test_srcdepends)
Rnaninf_test_qd_CXXFLAGS = $(test_cxxflags)
Rorg2l_test_qd_SOURCES = ../common/Rorg2l.test.cpp $(test_srcdepends)
Rorg2l_test_qd_CXXFLAGS = $(test_cxxflags)
Rorg2r_test_qd_SOURCES = ../common/Rorg2r.test.cpp $(test_srcdepends)
Rorg2r_test_qd_CXXFLAGS = $(test_cxxflags)
Rorgql_test_qd_SOURCES = ../common/Rorgql.test.cpp $(test_srcdepends)
Rorgql_test_qd_CXXFLAGS = $(test_cxxflags)
Rorgqr_test_qd_SOURCES = ../common/Rorgqr.test.cpp $(test_srcdepends)
Rorgqr_test_qd_CXXFLAGS = $(test_cxxflags)
Rorgtr_test_qd_SOURCES = ../common/Rorgtr.test.cpp $(test_srcdepends)
Rorgtr_test_qd_CXXFLAGS = $(test_cxxflags)
Rpocon_test_qd_SOURCES = ../common/Rpocon.test.cpp $(test_srcdepends)
Rpocon_test_qd_CXXFLAGS = $(test_cxxflags)
Rposv_test_qd_SOURCES = ../common/Rposv.test.cpp $(test_srcdepends)
Rposv_test_qd_CXXFLAGS = $(test_cxxflags)
Rpotf2_test_qd_SOURCES = ../common/Rpotf2.test.cpp $(test_srcdepends)
Rpotf2_test_qd_CXXFLAGS = $(test_cxxflags)
Rpotri_test_qd_SOURCES = ../common/Rpotri.test.cpp $(test_srcdepends)
Rpotri_test_qd_CXXFLAGS = $(test_cxxflags)
Rpotrs_test_qd_SOURCES = ../common/Rpotrs.test.cpp $(test_srcdepends)
Rpotrs_test_qd_CXXFLAGS = $(test_cxxflags)
Rsteqr_test_qd_SOURCES = ../common/Rsteqr.test.cpp $(test_srcdepends)
Rsteqr_test_qd_CXXFLAGS = $(test_cxxflags)
Rsterf_test_qd_SOURCES = ../common/Rsterf.test.cpp $(test_srcdepends)
Rsterf_test_qd_CXXFLAGS = $(test_cxxflags)
Rsyev_test_qd_SOURCES = ../common/Rsyev.test.cpp $(test_srcdepends)
Rsyev_test_qd_CXXFLAGS = $(test_cxxflags)
Rsytd2_test_qd_SOURCES = ../common/Rsytd2.test.cpp $(test_srcdepends)
Rsytd2_test_qd_CXXFLAGS = $(test_cxxflags)
Rsytrd_test_qd_SOURCES = ../common/Rsytrd.test.cpp $(test_srcdepends)
Rsytrd_test_qd_CXXFLAGS = $(test_cxxflags)
Rtrti2_test_qd_SOURCES = ../common/Rtrti2.test.cpp $(test_srcdepends)
Rtrti2_test_qd_CXXFLAGS = $(test_cxxflags)
Rtrtri_test_qd_SOURCES = ../common/Rtrtri.test.cpp $(test_srcdepends)
Rtrtri_test_qd_CXXFLAGS = $(test_cxxflags)
Rtrtrs_test_qd_SOURCES = ../common/Rtrtrs.test.cpp $(test_srcdepends)
Rtrtrs_test_qd_CXXFLAGS = $(test_cxxflags)
iCmax1_test_qd_SOURCES = ../common/iCmax1.test.cpp $(test_srcdepends)
iCmax1_test_qd_CXXFLAGS = $(test_cxxflags)
iMlaenv_test_qd_SOURCES = ../common/iMlaenv.test.cpp $(test_srcdepends)
iMlaenv_test_qd_CXXFLAGS = $(test_cxxflags)

include $(top_srcdir)/mplapack/test/compare/Makefile.backend.am
include $(top_srcdir)/mplapack/test/compare/Makefile.tests.am
