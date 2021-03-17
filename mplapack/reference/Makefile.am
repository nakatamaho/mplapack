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
lib_LTLIBRARIES += libmplapack_longdouble_ref.la
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
libmplapack_longdouble_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT64X___
libmplapack_longdouble_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack_longdouble_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE__FLOAT128
libmplapack__Float128_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT128___
libmplapack__Float128_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack__Float128_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

#done...90 (do not count Mulits.cpp)
MPLAPACK_SOURCES = Mutils.cpp \
iMlaenv.cpp Rlae2.cpp  Rlamch.cpp Rlaev2.cpp Claev2.cpp\
Rlassq.cpp Classq.cpp Rlanst.cpp Clanht.cpp Rlansy.cpp \
Clansy.cpp Clanhe.cpp Rlapy2.cpp Rlarfg.cpp Rlapy3.cpp \
Rladiv.cpp Cladiv.cpp Clarfg.cpp Rlartg.cpp Clartg.cpp \
Rlaset.cpp Claset.cpp Rlasr.cpp  Clasr.cpp  Rpotf2.cpp \
Clacgv.cpp Cpotf2.cpp Rlascl.cpp Clascl.cpp Rlasrt.cpp \
Rsteqr.cpp Csteqr.cpp Rsterf.cpp Rsytd2.cpp Chetd2.cpp \
Rlatrd.cpp Clatrd.cpp Rsytrd.cpp Chetrd.cpp Rlarf.cpp  \
Clarf.cpp  Rorg2l.cpp Cung2l.cpp Rorg2r.cpp Cung2r.cpp \
Rlarft.cpp Clarft.cpp Rlarfb.cpp Clarfb.cpp Rorgqr.cpp \
Cungqr.cpp Rorgql.cpp Cungql.cpp Rorgtr.cpp Cungtr.cpp \
Rsyev.cpp  Cheev.cpp  Rpotrf.cpp Cpotrf.cpp Clacrm.cpp \
Rtrti2.cpp Ctrti2.cpp Rtrtri.cpp Ctrtri.cpp Rgetf2.cpp \
Cgetf2.cpp Rlaswp.cpp Claswp.cpp Rgetrf.cpp Cgetrf.cpp \
Rgetri.cpp Cgetri.cpp Rgetrs.cpp Cgetrs.cpp Rgesv.cpp  \
Cgesv.cpp  Rtrtrs.cpp Ctrtrs.cpp Rlasyf.cpp Clasyf.cpp \
Clahef.cpp Clacrt.cpp Claesy.cpp Crot.cpp   Cspmv.cpp  \
Cspr.cpp   Csymv.cpp  Csyr.cpp   iCmax1.cpp RCsum1.cpp

#QA work in progress...4
MPLAPACK_SOURCES += Rlacn2.cpp Clacn2.cpp Rlatrs.cpp Clatrs.cpp

#nothing done yet, auxiliary routines...221
MPLAPACK_SOURCES += \
Rgbtf2.cpp Cgbtf2.cpp Rgebd2.cpp Cgebd2.cpp Rgehd2.cpp \
Cgehd2.cpp Rgelq2.cpp Cgelq2.cpp Rgeql2.cpp Cgeql2.cpp \
Rgeqr2.cpp Cgeqr2.cpp Rgerq2.cpp Cgerq2.cpp Rgesc2.cpp \
Cgesc2.cpp Rgetc2.cpp Cgetc2.cpp Rgtts2.cpp Cgtts2.cpp \
Rlabrd.cpp Clabrd.cpp Rlacpy.cpp Clacpy.cpp Clarcm.cpp \
Rlaebz.cpp Rlaed0.cpp Claed0.cpp Rlaed1.cpp Rlaed2.cpp \
Rlaed3.cpp Rlaed4.cpp Rlaed5.cpp Rlaed6.cpp Rlaed7.cpp \
Claed7.cpp Rlaed8.cpp Claed8.cpp Rlaed9.cpp Rlaeda.cpp \
Rlaein.cpp Claein.cpp Rlaexc.cpp Rlag2.cpp  Rlags2.cpp \
Clags2.cpp Rlagtf.cpp Rlagtm.cpp Clagtm.cpp Rlagts.cpp \
Rlagv2.cpp Rlahqr.cpp Clahqr.cpp Rlahr2.cpp Clahr2.cpp \
Rlahrd.cpp Clahrd.cpp Rlaic1.cpp Claic1.cpp Rlaln2.cpp \
Rlals0.cpp Clals0.cpp Rlalsa.cpp Clalsa.cpp Rlalsd.cpp \
Clalsd.cpp Rlamrg.cpp Rlangb.cpp Clangb.cpp Rlange.cpp \
Clange.cpp Rlangt.cpp Clangt.cpp Rlanhs.cpp Clanhs.cpp \
Rlansb.cpp Clansb.cpp Clanhb.cpp Rlansp.cpp Clansp.cpp \
Clanhp.cpp Rlantb.cpp Clantb.cpp Rlantp.cpp Clantp.cpp \
Rlantr.cpp Clantr.cpp Rlanv2.cpp Rlapll.cpp Clapll.cpp \
Rlapmt.cpp Clapmt.cpp Rlaqgb.cpp Claqgb.cpp Rlaqge.cpp \
Claqge.cpp Rlaqp2.cpp Claqp2.cpp Rlaqps.cpp Claqps.cpp \
Rlaqsb.cpp Claqsb.cpp Claqhb.cpp Rlaqsp.cpp Claqsp.cpp \
Claqhp.cpp Rlaqsy.cpp Claqsy.cpp Claqhe.cpp Rlaqtr.cpp \
Rlar1v.cpp Clar1v.cpp Rlar2v.cpp Clar2v.cpp Rlarfx.cpp \
Clarfx.cpp Rlargv.cpp Clargv.cpp Rlarnv.cpp Clarnv.cpp \
Rlaneg.cpp Rlarra.cpp Rlarrb.cpp Rlarrc.cpp Rlarrd.cpp \
Rlarre.cpp Rlarrf.cpp Rlarrj.cpp Rlarrk.cpp Rlarrr.cpp \
Rlartv.cpp Clartv.cpp Rlarrv.cpp Clarrv.cpp Rlaruv.cpp \
Rlarz.cpp  Clarz.cpp  Rlarzb.cpp Clarzb.cpp Rlarzt.cpp \
Clarzt.cpp Rlas2.cpp  Rlasd0.cpp Rlasd1.cpp Rlasd2.cpp \
Rlasd3.cpp Rlasd4.cpp Rlasd5.cpp Rlasd6.cpp Rlasd7.cpp \
Rlasd8.cpp Rlasda.cpp Rlasdq.cpp Rlasdt.cpp Rlasq1.cpp \
Rlasq2.cpp Rlasq3.cpp Rlasq4.cpp Rlasq5.cpp Rlasq6.cpp \
Rlasv2.cpp Rlasy2.cpp Rlatbs.cpp Clatbs.cpp Rlatdf.cpp \
Clatdf.cpp Rlatps.cpp Clatps.cpp Rlauu2.cpp Clauu2.cpp \
Rlauum.cpp Clauum.cpp Rlatrz.cpp Clatrz.cpp Rorgl2.cpp \
Cungl2.cpp Rorgr2.cpp Cungr2.cpp Rorm2l.cpp Cunm2l.cpp \
Rorm2r.cpp Cunm2r.cpp Rorml2.cpp Cunml2.cpp Rormr2.cpp \
Cunmr2.cpp Rormr3.cpp Cunmr3.cpp Rpbtf2.cpp Cpbtf2.cpp \
Rptts2.cpp Cptts2.cpp Rrscl.cpp  CRrscl.cpp Rsygs2.cpp \
Chegs2.cpp Rsytf2.cpp Csytf2.cpp Chetf2.cpp Rtgex2.cpp \
Ctgex2.cpp Rtgsy2.cpp Ctgsy2.cpp Rlazq3.cpp Rlazq4.cpp \
Rlaqr0.cpp Claqr0.cpp Rlaqr1.cpp Claqr1.cpp Rlaqr2.cpp \
Claqr2.cpp Rlaqr3.cpp Claqr3.cpp Rlaqr4.cpp Claqr4.cpp \
Rlaqr5.cpp Claqr5.cpp Rlatzm.cpp Clatzm.cpp Clacp2.cpp \
Risnan.cpp

#not yet done, driver and computational routines...351
MPLAPACK_SOURCES += \
Rbdsdc.cpp Rbdsqr.cpp Cbdsqr.cpp Rdisna.cpp Rgbbrd.cpp \
Cgbbrd.cpp Rgbcon.cpp Cgbcon.cpp Rgbequ.cpp Cgbequ.cpp \
Rgbrfs.cpp Cgbrfs.cpp Rgbsv.cpp  Cgbsv.cpp  Rgbsvx.cpp \
Cgbsvx.cpp Rgbtrf.cpp Cgbtrf.cpp Rgbtrs.cpp Cgbtrs.cpp \
Rgebak.cpp Cgebak.cpp Rgebal.cpp Cgebal.cpp Rgebrd.cpp \
Cgebrd.cpp Rgecon.cpp Cgecon.cpp Rgeequ.cpp Cgeequ.cpp \
Rgees.cpp  Cgees.cpp  Rgeesx.cpp Cgeesx.cpp Rgeev.cpp  \
Cgeev.cpp  Rgeevx.cpp Cgeevx.cpp Cgehrd.cpp Rgehrd.cpp \
Rgelqf.cpp Cgelqf.cpp Rgels.cpp  Cgels.cpp  Rgelsx.cpp \
Cgelsx.cpp Rgelsd.cpp Cgelsd.cpp Rgelss.cpp Cgelss.cpp \
Rgelsy.cpp Cgelsy.cpp Rgeqlf.cpp Cgeqlf.cpp Rgeqpf.cpp \
Cgeqpf.cpp Rgeqp3.cpp Cgeqp3.cpp Rgeqrf.cpp Cgeqrf.cpp \
Rgerfs.cpp Cgerfs.cpp Rgerqf.cpp Cgerqf.cpp Rgesdd.cpp \
Cgesdd.cpp Rgesvd.cpp Cgesvd.cpp Rgesvx.cpp Cgesvx.cpp \
Rggbak.cpp Cggbak.cpp Rggbal.cpp Cggbal.cpp Rgges.cpp  \
Cgges.cpp  Rgegs.cpp  Cgegs.cpp  Rggesx.cpp Cggesx.cpp \
Rggev.cpp  Cggev.cpp  Rgegv.cpp  Cgegv.cpp  Rggevx.cpp \
Cggevx.cpp Rggglm.cpp Cggglm.cpp Rgghrd.cpp Cgghrd.cpp \
Rgglse.cpp Cgglse.cpp Rggqrf.cpp Cggqrf.cpp Rggrqf.cpp \
Cggrqf.cpp Rggsvd.cpp Cggsvd.cpp Rggsvp.cpp Cggsvp.cpp \
Rgtcon.cpp Cgtcon.cpp Rgtrfs.cpp Cgtrfs.cpp Rgtsv.cpp  \
Cgtsv.cpp  Rgtsvx.cpp Cgtsvx.cpp Rgttrf.cpp Cgttrf.cpp \
Rgttrs.cpp Cgttrs.cpp Rhgeqz.cpp Chgeqz.cpp Rhsein.cpp \
Chsein.cpp Rhseqr.cpp Chseqr.cpp Ropgtr.cpp Cupgtr.cpp \
Ropmtr.cpp Cupmtr.cpp Rorgbr.cpp Cungbr.cpp Rorghr.cpp \
Cunghr.cpp Rorglq.cpp Cunglq.cpp Rorgrq.cpp Cungrq.cpp \
Rormbr.cpp Cunmbr.cpp Rormhr.cpp Cunmhr.cpp Rormlq.cpp \
Cunmlq.cpp Rormql.cpp Cunmql.cpp Rormqr.cpp Cunmqr.cpp \
Rormrq.cpp Cunmrq.cpp Rormrz.cpp Cunmrz.cpp Rormtr.cpp \
Cunmtr.cpp Rpbcon.cpp Cpbcon.cpp Rpbequ.cpp Cpbequ.cpp \
Rpbrfs.cpp Cpbrfs.cpp Rpbstf.cpp Cpbstf.cpp Rpbsv.cpp  \
Cpbsv.cpp  Rpbsvx.cpp Cpbsvx.cpp Rpbtrf.cpp Cpbtrf.cpp \
Rpbtrs.cpp Cpbtrs.cpp Rpocon.cpp Cpocon.cpp Rpoequ.cpp \
Cpoequ.cpp Rporfs.cpp Cporfs.cpp Rposv.cpp  Cposv.cpp  \
Rposvx.cpp Cposvx.cpp Rpotri.cpp Cpotri.cpp Rpotrs.cpp \
Cpotrs.cpp Rppcon.cpp Cppcon.cpp Rppequ.cpp Cppequ.cpp \
Rpptrf.cpp Cpptrf.cpp Rpptri.cpp Cpptri.cpp Rpptrs.cpp \
Cpptrs.cpp Rpprfs.cpp Cpprfs.cpp Rppsv.cpp  Cppsv.cpp  \
Rppsvx.cpp Cppsvx.cpp Rptcon.cpp Cptcon.cpp Rpteqr.cpp \
Cpteqr.cpp Rptrfs.cpp Cptrfs.cpp Rptsv.cpp  Cptsv.cpp  \
Rptsvx.cpp Cptsvx.cpp Rpttrf.cpp Cpttrf.cpp Rpttrs.cpp \
Cpttrs.cpp Rsbev.cpp  Chbev.cpp  Rsbevd.cpp Chbevd.cpp \
Rsbevx.cpp Chbevx.cpp Rsbgst.cpp Chbgst.cpp Rsbgv.cpp  \
Chbgv.cpp  Rsbgvd.cpp Chbgvd.cpp Rsbgvx.cpp Chbgvx.cpp \
Rsbtrd.cpp Chbtrd.cpp Rspcon.cpp Cspcon.cpp Chpcon.cpp \
Rspev.cpp  Chpev.cpp  Rspevd.cpp Chpevd.cpp Rspevx.cpp \
Chpevx.cpp Rspgst.cpp Chpgst.cpp Rspgv.cpp  Chpgv.cpp  \
Rspgvd.cpp Chpgvd.cpp Rspgvx.cpp Chpgvx.cpp Rsprfs.cpp \
Csprfs.cpp Chprfs.cpp Rspsv.cpp  Cspsv.cpp  Chpsv.cpp  \
Rspsvx.cpp Cspsvx.cpp Chpsvx.cpp Rsptrd.cpp Chptrd.cpp \
Rsptrf.cpp Csptrf.cpp Chptrf.cpp Rsptri.cpp Csptri.cpp \
Chptri.cpp Rsptrs.cpp Csptrs.cpp Chptrs.cpp Rstebz.cpp \
Rstedc.cpp Cstedc.cpp Rstegr.cpp Cstegr.cpp Rstein.cpp \
Cstein.cpp Rstev.cpp  Rstevd.cpp Rstevr.cpp Rstevx.cpp \
Rstemr.cpp Cstemr.cpp Rsycon.cpp Csycon.cpp Checon.cpp \
Rsyevd.cpp Cheevd.cpp Rsyevr.cpp Cheevr.cpp Rsyevx.cpp \
Cheevx.cpp Rsygst.cpp Chegst.cpp Rsygv.cpp  Chegv.cpp  \
Rsygvd.cpp Chegvd.cpp Rsygvx.cpp Chegvx.cpp Rsyrfs.cpp \
Csyrfs.cpp Cherfs.cpp Rsysv.cpp  Csysv.cpp  Chesv.cpp  \
Rsysvx.cpp Csysvx.cpp Chesvx.cpp Rsytrf.cpp Csytrf.cpp \
Chetrf.cpp Rsytri.cpp Csytri.cpp Chetri.cpp Rsytrs.cpp \
Csytrs.cpp Chetrs.cpp Rtbcon.cpp Ctbcon.cpp Rtbrfs.cpp \
Ctbrfs.cpp Rtbtrs.cpp Ctbtrs.cpp Rtgevc.cpp Ctgevc.cpp \
Rtgexc.cpp Ctgexc.cpp Rtgsen.cpp Ctgsen.cpp Rtgsja.cpp \
Ctgsja.cpp Rtgsna.cpp Ctgsna.cpp Rtgsyl.cpp Ctgsyl.cpp \
Rtpcon.cpp Ctpcon.cpp Rtprfs.cpp Ctprfs.cpp Rtptri.cpp \
Ctptri.cpp Rtptrs.cpp Ctptrs.cpp Rtrcon.cpp Ctrcon.cpp \
Rtrevc.cpp Ctrevc.cpp Rtrexc.cpp Ctrexc.cpp Rtrrfs.cpp \
Ctrrfs.cpp Rtrsen.cpp Ctrsen.cpp Rtrsna.cpp Ctrsna.cpp \
Rtrsyl.cpp Ctrsyl.cpp Rtzrzf.cpp Ctzrzf.cpp Rtzrqf.cpp \
Ctzrqf.cpp 

# # of MPLAPACK routines = 666.
# LAPACK 3.1.1 d*f, z*f, izmax1, ilaenv, dlamch = 674
# not implemented 8
# 674 - 8 = 666
#not to be implemtented
#1.Rlag2s.cpp : convert double prec. to single prec.
#2.Clag2c.cpp : convert double prec. to single prec.
#3.Rsgesv.cpp : single prec. i/o, internal double prec.
#4.Czgesv.cpp : single prec. i/o, internal double prec.
#5.Rlabad.cpp : Clay specific?
#6.Rlaisnan.cpp : use Risnan.cpp without optimization (?)
#7.Rlacon.cpp: use Rlacn2.
#8:Clacon.cpp: use Clacn2.

#will be removed 
#Rgegs.cpp:  deprecated by Rgges.cpp
#Rgegv.cpp:  deprecated by Rggev.cpp
#Rgelsx.cpp: deprecated by Rgelsy.cpp.
#Rgeqpf.cpp: deprecated by Rgeqp3.cpp
#Rlahrd.cpp: deprecated
#Rlatzm.cpp: deprecated by Rormrz.cpp
#Rtzrqf.cpp: deprecated by Rtzrzf.cpp
#Cgegs.cpp:  deprecated by Cgges.cpp
#Cgegv.cpp:  deprecated by Cggev.cpp
#Cgelsx.cpp: deprecated by Cgelsy.cpp
#Cgeqpf.cpp: deprecated by Cgeqp3.cpp
#Clahrd.cpp: deprecated
#Clatzm.cpp: deprecated by Cunmrz.cpp
#Ctzrqf.cpp: deprecated by Ctzrzf.cpp
