run_Rlaruv_test.sh: Makefile
	@echo '#!/bin/sh' > $@
	@echo '$${LOG_COMPILER} ./Rlaruv.test_gmp$${EXEEXT} || true' >> $@
	@echo 'if diff -u $${srcdir}/Rlaruv_reference.txt Rlaruv.txt > /dev/null 2>&1; then' >> $@
	@echo '  echo "PASS: Rlaruv"; exit 0' >> $@
	@echo 'else' >> $@
	@echo '  echo "FAIL: Rlaruv"; diff $${srcdir}/Rlaruv_reference.txt Rlaruv.txt; exit 1' >> $@
	@echo 'fi' >> $@
	@chmod +x $@

run_Rlamch_test.sh: Makefile
	@echo '#!/bin/sh' > $@
	@echo '$${LOG_COMPILER} ./Rlamch.test_gmp$${EXEEXT} || true' >> $@
	@echo 'if diff -u $${srcdir}/Rlamch_reference@ABI_BITS@.txt Rlamch.txt > /dev/null 2>&1; then' >> $@
	@echo '  echo "PASS: Rlamch (@ABI_BITS@-bit)"; exit 0' >> $@
	@echo 'else' >> $@
	@echo '  echo "FAIL: Rlamch (@ABI_BITS@-bit)"; diff $${srcdir}/Rlamch_reference@ABI_BITS@.txt Rlamch.txt; exit 1' >> $@
	@echo 'fi' >> $@
	@chmod +x $@

EXTRA_DIST = Rlamch_reference32.txt Rlamch_reference64.txt Rlaruv_reference.txt

backend_define = $(GMPFRXX_MKII_CPPFLAGS) -DMPLAPACK_BUILD_WITH_GMP
backend_cppflags = -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) $(backend_define) -DMPLAPACK_INTERNAL
backend_cxxflags =
backend_reference_archives = $(top_builddir)/mplapack/reference/.libs/libmplapack_gmp.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_mpfr.a
backend_reference_libs = -lmplapack_gmp -lmplapack_mpfr
backend_whole_archive_libs =
backend_link_libs =
backend_testdir = $(prefix)/lib/$(target)/mplapack/test/compare/gmp
mplibs = -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -lmpc -lmpfr -lgmp

compare_tests = \
    Cgelq2.test_gmp \
    Cgeqr2.test_gmp \
    Cgeqrf.test_gmp \
    Cgesv.test_gmp \
    Cgetf2.test_gmp \
    Cgetrf.test_gmp \
    Cgetri.test_gmp \
    Cgetrs.test_gmp \
    Cheev.test_gmp \
    Chetd2.test_gmp \
    Chetrd.test_gmp \
    Clacgv.test_gmp \
    Clacrm.test_gmp \
    Clacrt.test_gmp \
    Cladiv.test_gmp \
    Claesy.test_gmp \
    Claev2.test_gmp \
    Clahef.test_gmp \
    Clanhe.test_gmp \
    Clanht.test_gmp \
    Clansy.test_gmp \
    Clarf.test_gmp \
    Clarfb.test_gmp \
    Clarfg.test_gmp \
    Clarft.test_gmp \
    Clartg.test_gmp \
    Clascl.test_gmp \
    Claset.test_gmp \
    Clasr.test_gmp \
    Classq.test_gmp \
    Claswp.test_gmp \
    Clasyf.test_gmp \
    Clatrd.test_gmp \
    Cpotf2.test_gmp \
    Crot.test_gmp \
    Cspmv.test_gmp \
    Cspr.test_gmp \
    Csteqr.test_gmp \
    Csymv.test_gmp \
    Csyr.test_gmp \
    Ctrti2.test_gmp \
    Ctrtri.test_gmp \
    Ctrtrs.test_gmp \
    Cung2l.test_gmp \
    Cung2r.test_gmp \
    Cungql.test_gmp \
    Cungqr.test_gmp \
    Cungtr.test_gmp \
    Cunmqr.test_gmp \
    Mlsamen.test_gmp \
    Mutils.test_gmp \
    RCsum1.test_gmp \
    Rgecon.test_gmp \
    Rgeequ.test_gmp \
    Rgeqlf.test_gmp \
    Rgesv.test_gmp \
    Rgetf2.test_gmp \
    Rgetrf.test_gmp \
    Rgetri.test_gmp \
    Rgetrs.test_gmp \
    Rladiv.test_gmp \
    Rlae2.test_gmp \
    Rlaev2.test_gmp \
    Rlamch.test_gmp \
    Rlange.test_gmp \
    Rlanst.test_gmp \
    Rlansy.test_gmp \
    Rlapy2.test_gmp \
    Rlapy3.test_gmp \
    Rlarf.test_gmp \
    Rlarfb.test_gmp \
    Rlarfg.test_gmp \
    Rlarft.test_gmp \
    Rlartg.test_gmp \
    Rlaruv.test_gmp \
    Rlascl.test_gmp \
    Rlaset.test_gmp \
    Rlasr.test_gmp \
    Rlasrt.test_gmp \
    Rlassq.test_gmp \
    Rlaswp.test_gmp \
    Rlasyf.test_gmp \
    Rlatrd.test_gmp \
    Rlatrs.test_gmp \
    Rlauu2.test_gmp \
    Rlauum.test_gmp \
    Rnaninf.test_gmp \
    Rorg2l.test_gmp \
    Rorg2r.test_gmp \
    Rorgql.test_gmp \
    Rorgqr.test_gmp \
    Rorgtr.test_gmp \
    Rpocon.test_gmp \
    Rposv.test_gmp \
    Rpotf2.test_gmp \
    Rpotri.test_gmp \
    Rpotrs.test_gmp \
    Rsteqr.test_gmp \
    Rsterf.test_gmp \
    Rsyev.test_gmp \
    Rsytd2.test_gmp \
    Rsytrd.test_gmp \
    Rtrti2.test_gmp \
    Rtrtri.test_gmp \
    Rtrtrs.test_gmp \
    iCmax1.test_gmp \
    iMlaenv.test_gmp

Cgelq2_test_gmp_SOURCES = ../common/Cgelq2.test.cpp $(test_srcdepends)
Cgelq2_test_gmp_CXXFLAGS = $(test_cxxflags)
Cgeqr2_test_gmp_SOURCES = ../common/Cgeqr2.test.cpp $(test_srcdepends)
Cgeqr2_test_gmp_CXXFLAGS = $(test_cxxflags)
Cgeqrf_test_gmp_SOURCES = ../common/Cgeqrf.test.cpp $(test_srcdepends)
Cgeqrf_test_gmp_CXXFLAGS = $(test_cxxflags)
Cgesv_test_gmp_SOURCES = ../common/Cgesv.test.cpp $(test_srcdepends)
Cgesv_test_gmp_CXXFLAGS = $(test_cxxflags)
Cgetf2_test_gmp_SOURCES = ../common/Cgetf2.test.cpp $(test_srcdepends)
Cgetf2_test_gmp_CXXFLAGS = $(test_cxxflags)
Cgetrf_test_gmp_SOURCES = ../common/Cgetrf.test.cpp $(test_srcdepends)
Cgetrf_test_gmp_CXXFLAGS = $(test_cxxflags)
Cgetri_test_gmp_SOURCES = ../common/Cgetri.test.cpp $(test_srcdepends)
Cgetri_test_gmp_CXXFLAGS = $(test_cxxflags)
Cgetrs_test_gmp_SOURCES = ../common/Cgetrs.test.cpp $(test_srcdepends)
Cgetrs_test_gmp_CXXFLAGS = $(test_cxxflags)
Cheev_test_gmp_SOURCES = ../common/Cheev.test.cpp $(test_srcdepends)
Cheev_test_gmp_CXXFLAGS = $(test_cxxflags)
Chetd2_test_gmp_SOURCES = ../common/Chetd2.test.cpp $(test_srcdepends)
Chetd2_test_gmp_CXXFLAGS = $(test_cxxflags)
Chetrd_test_gmp_SOURCES = ../common/Chetrd.test.cpp $(test_srcdepends)
Chetrd_test_gmp_CXXFLAGS = $(test_cxxflags)
Clacgv_test_gmp_SOURCES = ../common/Clacgv.test.cpp $(test_srcdepends)
Clacgv_test_gmp_CXXFLAGS = $(test_cxxflags)
Clacrm_test_gmp_SOURCES = ../common/Clacrm.test.cpp $(test_srcdepends)
Clacrm_test_gmp_CXXFLAGS = $(test_cxxflags)
Clacrt_test_gmp_SOURCES = ../common/Clacrt.test.cpp $(test_srcdepends)
Clacrt_test_gmp_CXXFLAGS = $(test_cxxflags)
Cladiv_test_gmp_SOURCES = ../common/Cladiv.test.cpp $(test_srcdepends)
Cladiv_test_gmp_CXXFLAGS = $(test_cxxflags)
Claesy_test_gmp_SOURCES = ../common/Claesy.test.cpp $(test_srcdepends)
Claesy_test_gmp_CXXFLAGS = $(test_cxxflags)
Claev2_test_gmp_SOURCES = ../common/Claev2.test.cpp $(test_srcdepends)
Claev2_test_gmp_CXXFLAGS = $(test_cxxflags)
Clahef_test_gmp_SOURCES = ../common/Clahef.test.cpp $(test_srcdepends)
Clahef_test_gmp_CXXFLAGS = $(test_cxxflags)
Clanhe_test_gmp_SOURCES = ../common/Clanhe.test.cpp $(test_srcdepends)
Clanhe_test_gmp_CXXFLAGS = $(test_cxxflags)
Clanht_test_gmp_SOURCES = ../common/Clanht.test.cpp $(test_srcdepends)
Clanht_test_gmp_CXXFLAGS = $(test_cxxflags)
Clansy_test_gmp_SOURCES = ../common/Clansy.test.cpp $(test_srcdepends)
Clansy_test_gmp_CXXFLAGS = $(test_cxxflags)
Clarf_test_gmp_SOURCES = ../common/Clarf.test.cpp $(test_srcdepends)
Clarf_test_gmp_CXXFLAGS = $(test_cxxflags)
Clarfb_test_gmp_SOURCES = ../common/Clarfb.test.cpp $(test_srcdepends)
Clarfb_test_gmp_CXXFLAGS = $(test_cxxflags)
Clarfg_test_gmp_SOURCES = ../common/Clarfg.test.cpp $(test_srcdepends)
Clarfg_test_gmp_CXXFLAGS = $(test_cxxflags)
Clarft_test_gmp_SOURCES = ../common/Clarft.test.cpp $(test_srcdepends)
Clarft_test_gmp_CXXFLAGS = $(test_cxxflags)
Clartg_test_gmp_SOURCES = ../common/Clartg.test.cpp $(test_srcdepends)
Clartg_test_gmp_CXXFLAGS = $(test_cxxflags)
Clascl_test_gmp_SOURCES = ../common/Clascl.test.cpp $(test_srcdepends)
Clascl_test_gmp_CXXFLAGS = $(test_cxxflags)
Claset_test_gmp_SOURCES = ../common/Claset.test.cpp $(test_srcdepends)
Claset_test_gmp_CXXFLAGS = $(test_cxxflags)
Clasr_test_gmp_SOURCES = ../common/Clasr.test.cpp $(test_srcdepends)
Clasr_test_gmp_CXXFLAGS = $(test_cxxflags)
Classq_test_gmp_SOURCES = ../common/Classq.test.cpp $(test_srcdepends)
Classq_test_gmp_CXXFLAGS = $(test_cxxflags)
Claswp_test_gmp_SOURCES = ../common/Claswp.test.cpp $(test_srcdepends)
Claswp_test_gmp_CXXFLAGS = $(test_cxxflags)
Clasyf_test_gmp_SOURCES = ../common/Clasyf.test.cpp $(test_srcdepends)
Clasyf_test_gmp_CXXFLAGS = $(test_cxxflags)
Clatrd_test_gmp_SOURCES = ../common/Clatrd.test.cpp $(test_srcdepends)
Clatrd_test_gmp_CXXFLAGS = $(test_cxxflags)
Cpotf2_test_gmp_SOURCES = ../common/Cpotf2.test.cpp $(test_srcdepends)
Cpotf2_test_gmp_CXXFLAGS = $(test_cxxflags)
Crot_test_gmp_SOURCES = ../common/Crot.test.cpp $(test_srcdepends)
Crot_test_gmp_CXXFLAGS = $(test_cxxflags)
Cspmv_test_gmp_SOURCES = ../common/Cspmv.test.cpp $(test_srcdepends)
Cspmv_test_gmp_CXXFLAGS = $(test_cxxflags)
Cspr_test_gmp_SOURCES = ../common/Cspr.test.cpp $(test_srcdepends)
Cspr_test_gmp_CXXFLAGS = $(test_cxxflags)
Csteqr_test_gmp_SOURCES = ../common/Csteqr.test.cpp $(test_srcdepends)
Csteqr_test_gmp_CXXFLAGS = $(test_cxxflags)
Csymv_test_gmp_SOURCES = ../common/Csymv.test.cpp $(test_srcdepends)
Csymv_test_gmp_CXXFLAGS = $(test_cxxflags)
Csyr_test_gmp_SOURCES = ../common/Csyr.test.cpp $(test_srcdepends)
Csyr_test_gmp_CXXFLAGS = $(test_cxxflags)
Ctrti2_test_gmp_SOURCES = ../common/Ctrti2.test.cpp $(test_srcdepends)
Ctrti2_test_gmp_CXXFLAGS = $(test_cxxflags)
Ctrtri_test_gmp_SOURCES = ../common/Ctrtri.test.cpp $(test_srcdepends)
Ctrtri_test_gmp_CXXFLAGS = $(test_cxxflags)
Ctrtrs_test_gmp_SOURCES = ../common/Ctrtrs.test.cpp $(test_srcdepends)
Ctrtrs_test_gmp_CXXFLAGS = $(test_cxxflags)
Cung2l_test_gmp_SOURCES = ../common/Cung2l.test.cpp $(test_srcdepends)
Cung2l_test_gmp_CXXFLAGS = $(test_cxxflags)
Cung2r_test_gmp_SOURCES = ../common/Cung2r.test.cpp $(test_srcdepends)
Cung2r_test_gmp_CXXFLAGS = $(test_cxxflags)
Cungql_test_gmp_SOURCES = ../common/Cungql.test.cpp $(test_srcdepends)
Cungql_test_gmp_CXXFLAGS = $(test_cxxflags)
Cungqr_test_gmp_SOURCES = ../common/Cungqr.test.cpp $(test_srcdepends)
Cungqr_test_gmp_CXXFLAGS = $(test_cxxflags)
Cungtr_test_gmp_SOURCES = ../common/Cungtr.test.cpp $(test_srcdepends)
Cungtr_test_gmp_CXXFLAGS = $(test_cxxflags)
Cunmqr_test_gmp_SOURCES = ../common/Cunmqr.test.cpp $(test_srcdepends)
Cunmqr_test_gmp_CXXFLAGS = $(test_cxxflags)
Mlsamen_test_gmp_SOURCES = ../common/Mlsamen.test.cpp $(test_srcdepends)
Mlsamen_test_gmp_CXXFLAGS = $(test_cxxflags)
Mutils_test_gmp_SOURCES = ../common/Mutils.test.cpp $(test_srcdepends) ../common/ilaenv.cpp
Mutils_test_gmp_CXXFLAGS = $(test_cxxflags)
RCsum1_test_gmp_SOURCES = ../common/RCsum1.test.cpp $(test_srcdepends)
RCsum1_test_gmp_CXXFLAGS = $(test_cxxflags)
Rgecon_test_gmp_SOURCES = ../common/Rgecon.test.cpp $(test_srcdepends)
Rgecon_test_gmp_CXXFLAGS = $(test_cxxflags)
Rgeequ_test_gmp_SOURCES = ../common/Rgeequ.test.cpp $(test_srcdepends)
Rgeequ_test_gmp_CXXFLAGS = $(test_cxxflags)
Rgeqlf_test_gmp_SOURCES = ../common/Rgeqlf.test.cpp $(test_srcdepends)
Rgeqlf_test_gmp_CXXFLAGS = $(test_cxxflags)
Rgesv_test_gmp_SOURCES = ../common/Rgesv.test.cpp $(test_srcdepends)
Rgesv_test_gmp_CXXFLAGS = $(test_cxxflags)
Rgetf2_test_gmp_SOURCES = ../common/Rgetf2.test.cpp $(test_srcdepends)
Rgetf2_test_gmp_CXXFLAGS = $(test_cxxflags)
Rgetrf_test_gmp_SOURCES = ../common/Rgetrf.test.cpp $(test_srcdepends)
Rgetrf_test_gmp_CXXFLAGS = $(test_cxxflags)
Rgetri_test_gmp_SOURCES = ../common/Rgetri.test.cpp $(test_srcdepends)
Rgetri_test_gmp_CXXFLAGS = $(test_cxxflags)
Rgetrs_test_gmp_SOURCES = ../common/Rgetrs.test.cpp $(test_srcdepends)
Rgetrs_test_gmp_CXXFLAGS = $(test_cxxflags)
Rladiv_test_gmp_SOURCES = ../common/Rladiv.test.cpp $(test_srcdepends)
Rladiv_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlae2_test_gmp_SOURCES = ../common/Rlae2.test.cpp $(test_srcdepends)
Rlae2_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlaev2_test_gmp_SOURCES = ../common/Rlaev2.test.cpp $(test_srcdepends)
Rlaev2_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlamch_test_gmp_SOURCES = ../common/Rlamch.test.cpp $(test_srcdepends)
Rlamch_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlange_test_gmp_SOURCES = ../common/Rlange.test.cpp $(test_srcdepends)
Rlange_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlanst_test_gmp_SOURCES = ../common/Rlanst.test.cpp $(test_srcdepends)
Rlanst_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlansy_test_gmp_SOURCES = ../common/Rlansy.test.cpp $(test_srcdepends)
Rlansy_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlapy2_test_gmp_SOURCES = ../common/Rlapy2.test.cpp $(test_srcdepends)
Rlapy2_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlapy3_test_gmp_SOURCES = ../common/Rlapy3.test.cpp $(test_srcdepends)
Rlapy3_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlarf_test_gmp_SOURCES = ../common/Rlarf.test.cpp $(test_srcdepends)
Rlarf_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlarfb_test_gmp_SOURCES = ../common/Rlarfb.test.cpp $(test_srcdepends)
Rlarfb_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlarfg_test_gmp_SOURCES = ../common/Rlarfg.test.cpp $(test_srcdepends)
Rlarfg_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlarft_test_gmp_SOURCES = ../common/Rlarft.test.cpp $(test_srcdepends)
Rlarft_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlartg_test_gmp_SOURCES = ../common/Rlartg.test.cpp $(test_srcdepends)
Rlartg_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlaruv_test_gmp_SOURCES = ../common/Rlaruv.test.cpp $(test_srcdepends)
Rlaruv_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlascl_test_gmp_SOURCES = ../common/Rlascl.test.cpp $(test_srcdepends)
Rlascl_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlaset_test_gmp_SOURCES = ../common/Rlaset.test.cpp $(test_srcdepends)
Rlaset_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlasr_test_gmp_SOURCES = ../common/Rlasr.test.cpp $(test_srcdepends)
Rlasr_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlasrt_test_gmp_SOURCES = ../common/Rlasrt.test.cpp $(test_srcdepends)
Rlasrt_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlassq_test_gmp_SOURCES = ../common/Rlassq.test.cpp $(test_srcdepends)
Rlassq_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlaswp_test_gmp_SOURCES = ../common/Rlaswp.test.cpp $(test_srcdepends)
Rlaswp_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlasyf_test_gmp_SOURCES = ../common/Rlasyf.test.cpp $(test_srcdepends)
Rlasyf_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlatrd_test_gmp_SOURCES = ../common/Rlatrd.test.cpp $(test_srcdepends)
Rlatrd_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlatrs_test_gmp_SOURCES = ../common/Rlatrs.test.cpp $(test_srcdepends)
Rlatrs_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlauu2_test_gmp_SOURCES = ../common/Rlauu2.test.cpp $(test_srcdepends)
Rlauu2_test_gmp_CXXFLAGS = $(test_cxxflags)
Rlauum_test_gmp_SOURCES = ../common/Rlauum.test.cpp $(test_srcdepends)
Rlauum_test_gmp_CXXFLAGS = $(test_cxxflags)
Rnaninf_test_gmp_SOURCES = ../common/Rnaninf.test.cpp $(test_srcdepends)
Rnaninf_test_gmp_CXXFLAGS = $(test_cxxflags)
Rorg2l_test_gmp_SOURCES = ../common/Rorg2l.test.cpp $(test_srcdepends)
Rorg2l_test_gmp_CXXFLAGS = $(test_cxxflags)
Rorg2r_test_gmp_SOURCES = ../common/Rorg2r.test.cpp $(test_srcdepends)
Rorg2r_test_gmp_CXXFLAGS = $(test_cxxflags)
Rorgql_test_gmp_SOURCES = ../common/Rorgql.test.cpp $(test_srcdepends)
Rorgql_test_gmp_CXXFLAGS = $(test_cxxflags)
Rorgqr_test_gmp_SOURCES = ../common/Rorgqr.test.cpp $(test_srcdepends)
Rorgqr_test_gmp_CXXFLAGS = $(test_cxxflags)
Rorgtr_test_gmp_SOURCES = ../common/Rorgtr.test.cpp $(test_srcdepends)
Rorgtr_test_gmp_CXXFLAGS = $(test_cxxflags)
Rpocon_test_gmp_SOURCES = ../common/Rpocon.test.cpp $(test_srcdepends)
Rpocon_test_gmp_CXXFLAGS = $(test_cxxflags)
Rposv_test_gmp_SOURCES = ../common/Rposv.test.cpp $(test_srcdepends)
Rposv_test_gmp_CXXFLAGS = $(test_cxxflags)
Rpotf2_test_gmp_SOURCES = ../common/Rpotf2.test.cpp $(test_srcdepends)
Rpotf2_test_gmp_CXXFLAGS = $(test_cxxflags)
Rpotri_test_gmp_SOURCES = ../common/Rpotri.test.cpp $(test_srcdepends)
Rpotri_test_gmp_CXXFLAGS = $(test_cxxflags)
Rpotrs_test_gmp_SOURCES = ../common/Rpotrs.test.cpp $(test_srcdepends)
Rpotrs_test_gmp_CXXFLAGS = $(test_cxxflags)
Rsteqr_test_gmp_SOURCES = ../common/Rsteqr.test.cpp $(test_srcdepends)
Rsteqr_test_gmp_CXXFLAGS = $(test_cxxflags)
Rsterf_test_gmp_SOURCES = ../common/Rsterf.test.cpp $(test_srcdepends)
Rsterf_test_gmp_CXXFLAGS = $(test_cxxflags)
Rsyev_test_gmp_SOURCES = ../common/Rsyev.test.cpp $(test_srcdepends)
Rsyev_test_gmp_CXXFLAGS = $(test_cxxflags)
Rsytd2_test_gmp_SOURCES = ../common/Rsytd2.test.cpp $(test_srcdepends)
Rsytd2_test_gmp_CXXFLAGS = $(test_cxxflags)
Rsytrd_test_gmp_SOURCES = ../common/Rsytrd.test.cpp $(test_srcdepends)
Rsytrd_test_gmp_CXXFLAGS = $(test_cxxflags)
Rtrti2_test_gmp_SOURCES = ../common/Rtrti2.test.cpp $(test_srcdepends)
Rtrti2_test_gmp_CXXFLAGS = $(test_cxxflags)
Rtrtri_test_gmp_SOURCES = ../common/Rtrtri.test.cpp $(test_srcdepends)
Rtrtri_test_gmp_CXXFLAGS = $(test_cxxflags)
Rtrtrs_test_gmp_SOURCES = ../common/Rtrtrs.test.cpp $(test_srcdepends)
Rtrtrs_test_gmp_CXXFLAGS = $(test_cxxflags)
iCmax1_test_gmp_SOURCES = ../common/iCmax1.test.cpp $(test_srcdepends)
iCmax1_test_gmp_CXXFLAGS = $(test_cxxflags)
iMlaenv_test_gmp_SOURCES = ../common/iMlaenv.test.cpp $(test_srcdepends)
iMlaenv_test_gmp_CXXFLAGS = $(test_cxxflags)

include $(top_srcdir)/mplapack/test/compare/Makefile.backend.am
include $(top_srcdir)/mplapack/test/compare/Makefile.tests.am
