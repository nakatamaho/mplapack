lib_LTLIBRARIES = 

if ENABLE_MPFR
lib_LTLIBRARIES += libmplapack_mpfr_ref.la
endif

if ENABLE_GMP
lib_LTLIBRARIES += libmplapack_gmp_ref.la
endif

if ENABLE_QD
lib_LTLIBRARIES += libmplapack_qd_ref.la
endif

if ENABLE_DD
lib_LTLIBRARIES += libmplapack_dd_ref.la
endif

if ENABLE_DOUBLE
lib_LTLIBRARIES += libmplapack_double_ref.la
endif

if ENABLE__FLOAT64X
lib_LTLIBRARIES += libmplapack__Float64x_ref.la
endif

if ENABLE__FLOAT128
lib_LTLIBRARIES += libmplapack__Float128_ref.la
endif

if ENABLE_GMP
libmplapack_gmp_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___
libmplapack_gmp_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack_gmp_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_MPFR
libmplapack_mpfr_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___
libmplapack_mpfr_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack_mpfr_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_QD
libmplapack_qd_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___
libmplapack_qd_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack_qd_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DD
libmplapack_dd_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___
libmplapack_dd_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack_dd_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DOUBLE
libmplapack_double_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_DOUBLE___
libmplapack_double_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack_double_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE__FLOAT64X
libmplapack__Float64x_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT64X___
libmplapack__Float64x_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack__Float64x_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE__FLOAT128
libmplapack__Float128_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT128___
libmplapack__Float128_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack__Float128_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

MPLAPACK_SOURCES = \
Rgebrd.cpp \
Rgesvd.cpp \
Rormbr.cpp \
Rormql.cpp \
Rormqr.cpp \
Rorml2.cpp \
Rorm2r.cpp \
Rorglq.cpp \
Rorgl2.cpp \
Rgeqrf.cpp \
Rgelqf.cpp \
Rgeqr2.cpp \
Chetri.cpp \
Chetri_rook.cpp \
Clarfy.cpp \
Rgelq2.cpp \
Rlabrd.cpp \
Rgebd2.cpp \
Rbdsqr.cpp \
iMlazlr.cpp \
iMlazlc.cpp \
iMladlr.cpp \
iMladlc.cpp \
Claev2.cpp \
Clansy.cpp \
Clanhe.cpp \
Cladiv.cpp \
Clarfg.cpp \
Clartg.cpp \
Cpotf2.cpp \
Clascl.cpp \
Chetd2.cpp \
Clatrd.cpp \
Rsytrd.cpp \
Chetrd.cpp \
Rlarf.cpp \
Clarf.cpp \
Cung2l.cpp \
Cung2r.cpp \
Clarft.cpp \
Clarfb.cpp \
Cungql.cpp \
Cheev.cpp \
Clacrm.cpp \
Rtrtri.cpp \
Ctrtri.cpp \
Cgetri.cpp \
Clasyf.cpp \
Clahef.cpp \
Claesy.cpp \
Cdrscl.cpp \
Cgebak.cpp \
Cgelq.cpp \
Cgelqt.cpp \
Cgemlqt.cpp \
Cgemqrt.cpp \
Cgeqr.cpp \
Cgeqrt.cpp \
Cgesv.cpp \
Cgetc2.cpp \
Cgetf2.cpp \
Cgetrf.cpp \
Cgetrf2.cpp \
Cgetrs.cpp \
Cggbak.cpp \
Cgtts2.cpp \
Chegst.cpp \
Chesv.cpp \
Chesv_rook.cpp \
Cheswapr.cpp \
Chetri2.cpp \
Chetri_3.cpp \
Chetrs.cpp \
Chetrs2.cpp \
Chetrs_3.cpp \
Chetrs_aa.cpp \
Chetrs_rook.cpp \
Chfrk.cpp \
Chpsv.cpp \
Cla_wwaddw.cpp \
Clacgv.cpp \
Clacp2.cpp \
Clacpy.cpp \
Clacrt.cpp \
Clagtm.cpp \
Clangt.cpp \
Clanht.cpp \
Clapmr.cpp \
Clapmt.cpp \
Claqgb.cpp \
Claqge.cpp \
Claqhb.cpp \
Claqhe.cpp \
Claqhp.cpp \
Claqsb.cpp \
Claqsp.cpp \
Claqsy.cpp \
Clar2v.cpp \
Clarfb_gett.cpp \
Clarfx.cpp \
Clarscl2.cpp \
Clartv.cpp \
Clarz.cpp \
Clarzt.cpp \
Clascl2.cpp \
Claset.cpp \
Clasr.cpp \
Classq.cpp \
Claswlq.cpp \
Claswp.cpp \
Clatsqr.cpp \
Claunhr_col_getrfnp.cpp \
Clauum.cpp \
Cpbequ.cpp \
Cpbsv.cpp \
Cpbtrs.cpp \
Cpoequ.cpp \
Cposv.cpp \
Cpotrf.cpp \
Cpotrf2.cpp \
Cpotri.cpp \
Cpotrs.cpp \
Cppequ.cpp \
Cppsv.cpp \
Cpptrs.cpp \
Cptcon.cpp \
Cptsv.cpp \
Cpttrf.cpp \
Cpttrs.cpp \
Cptts2.cpp \
Crot.cpp \
Cspmv.cpp \
Cspr.cpp \
Cspsv.cpp \
Csyconv.cpp \
Csyconvf.cpp \
Csyconvf_rook.cpp \
Csymv.cpp \
Csyr.cpp \
Csyswapr.cpp \
Csytri2.cpp \
Csytri_3.cpp \
Csytrs.cpp \
Csytrs2.cpp \
Csytrs_3.cpp \
Csytrs_aa.cpp \
Csytrs_rook.cpp \
Ctbtrs.cpp \
Ctgexc.cpp \
Ctplqt.cpp \
Ctpmlqt.cpp \
Ctpmqrt.cpp \
Ctpqrt.cpp \
Ctprfb.cpp \
Ctptrs.cpp \
Ctpttr.cpp \
Ctrti2.cpp \
Ctrtrs.cpp \
Ctrttp.cpp \
Cunbdb5.cpp \
Cunbdb6.cpp \
Cunghr.cpp \
Cungl2.cpp \
Cungqr.cpp \
Cungr2.cpp \
Cungtr.cpp \
Cungtsqr.cpp \
Cunhr_col.cpp \
Cunm2l.cpp \
Cunm2r.cpp \
Cunml2.cpp \
Cunmr2.cpp \
Cunmr3.cpp \
Mutils.cpp \
RCsum1.cpp \
Rcombssq.cpp \
Rgebak.cpp \
Rgelq.cpp \
Rgelqt.cpp \
Rgemlqt.cpp \
Rgemqrt.cpp \
Rgeqr.cpp \
Rgeqrt.cpp \
Rgesv.cpp \
Rgetc2.cpp \
Rgetf2.cpp \
Rgetrf.cpp \
Rgetrf2.cpp \
Rgetri.cpp \
Rgetrs.cpp \
Rggbak.cpp \
Rgtsv.cpp \
Rgttrf.cpp \
Rgtts2.cpp \
Risnan.cpp \
Rla_gerpvgrw.cpp \
Rla_wwaddw.cpp \
Rlabad.cpp \
Rlacpy.cpp \
Rladiv.cpp \
Rlae2.cpp \
Rlaed5.cpp \
Rlaev2.cpp \
Rlagtf.cpp \
Rlagtm.cpp \
Rlaic1.cpp \
Rlaisnan.cpp \
Rlamch.cpp \
Rlamrg.cpp \
Rlamswlq.cpp \
Rlaneg.cpp \
Rlangt.cpp \
Rlanst.cpp \
Rlansy.cpp \
Rlanv2.cpp \
Rlaorhr_col_getrfnp.cpp \
Rlaorhr_col_getrfnp2.cpp \
Rlapmr.cpp \
Rlapmt.cpp \
Rlapy2.cpp \
Rlapy3.cpp \
Rlaqgb.cpp \
Rlaqge.cpp \
Rlaqr1.cpp \
Rlaqsb.cpp \
Rlaqsp.cpp \
Rlaqsy.cpp \
Rlar1v.cpp \
Rlar2v.cpp \
Rlarfb.cpp \
Rlarfb_gett.cpp \
Rlarfg.cpp \
Rlarfgp.cpp \
Rlarft.cpp \
Rlarfx.cpp \
Rlarfy.cpp \
Rlargv.cpp \
Rlarra.cpp \
Rlarrc.cpp \
Rlarrj.cpp \
Rlarrr.cpp \
Rlarscl2.cpp \
Rlartg.cpp \
Rlartgp.cpp \
Rlartgs.cpp \
Rlartv.cpp \
Rlarz.cpp \
Rlarzt.cpp \
Rlas2.cpp \
Rlascl.cpp \
Rlascl2.cpp \
Rlasd5.cpp \
Rlaset.cpp \
Rlasq1.cpp \
Rlasq4.cpp \
Rlasq6.cpp \
Rlasr.cpp \
Rlasrt.cpp \
Rlassq.cpp \
Rlasv2.cpp \
Rlaswlq.cpp \
Rlaswp.cpp \
Rlasyf.cpp \
Rlasyf_aa.cpp \
Rlasyf_rk.cpp \
Rlasyf_rook.cpp \
Rlatrd.cpp \
Rlatsqr.cpp \
Rlauu2.cpp \
Rlauum.cpp \
Rorbdb5.cpp \
Rorbdb6.cpp \
Rorg2l.cpp \
Rorg2r.cpp \
Rorgbr.cpp \
Rorghr.cpp \
Rorgql.cpp \
Rorgqr.cpp \
Rorgtr.cpp \
Rorhr_col.cpp \
Rpbequ.cpp \
Rpbsv.cpp \
Rpbtrs.cpp \
Rpoequ.cpp \
Rpoequb.cpp \
Rposv.cpp \
Rposvx.cpp \
Rpotf2.cpp \
Rpotrf.cpp \
Rpotrf2.cpp \
Rpotri.cpp \
Rpotrs.cpp \
Rppequ.cpp \
Rppsv.cpp \
Rpptrs.cpp \
Rptcon.cpp \
Rptsv.cpp \
Rptsvx.cpp \
Rpttrf.cpp \
Rpttrs.cpp \
Rptts2.cpp \
Rrscl.cpp \
Rsfrk.cpp \
Rspsv.cpp \
Rsteqr.cpp \
Rsterf.cpp \
Rstev.cpp \
Rstevd.cpp \
Rsyconv.cpp \
Rsyconvf.cpp \
Rsyconvf_rook.cpp \
Rsyequb.cpp \
Rsyev.cpp \
Rsyevd.cpp \
Rsygs2.cpp \
Rsygst.cpp \
Rsysv.cpp \
Rsysv_aa.cpp \
Rsysv_aa_2stage.cpp \
Rsysv_rk.cpp \
Rsysv_rook.cpp \
Rsysvx.cpp \
Rsyswapr.cpp \
Rsytd2.cpp \
Rsytf2.cpp \
Rsytf2_rk.cpp \
Rsytf2_rook.cpp \
Rsytri.cpp \
Rsytri2.cpp \
Rsytri_3.cpp \
Rsytri_rook.cpp \
Rsytrs.cpp \
Rsytrs2.cpp \
Rsytrs_3.cpp \
Rsytrs_aa.cpp \
Rsytrs_aa_2stage.cpp \
Rsytrs_rook.cpp \
Rtbtrs.cpp \
Rtgexc.cpp \
Rtplqt.cpp \
Rtpmlqt.cpp \
Rtpmqrt.cpp \
Rtpqrt.cpp \
Rtprfb.cpp \
Rtptrs.cpp \
Rtpttr.cpp \
Rtrexc.cpp \
Rtrti2.cpp \
Rtrtrs.cpp \
Rtrttp.cpp \
iCmax1.cpp \
Csteqr.cpp \
iMlaenv.cpp 
