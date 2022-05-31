SUBDIRS = . 
lib_LTLIBRARIES = 

mplapack_lin_DATA = Ctest.in Ctest_rfp.in Rtest.in Rtest_rfp.in
mplapack_lindir = $(prefix)/lib/$(target)/mplapack/test/lin

if ENABLE_MPFR
lib_LTLIBRARIES += liblin_mpfr.la liblin_override_mpfr.la
SUBDIRS += mpfr
endif

if ENABLE_GMP
lib_LTLIBRARIES += liblin_gmp.la liblin_override_gmp.la
SUBDIRS += gmp
endif

if ENABLE_DD
lib_LTLIBRARIES += liblin_dd.la liblin_override_dd.la
SUBDIRS += dd
endif

if ENABLE_QD
lib_LTLIBRARIES += liblin_qd.la liblin_override_qd.la
SUBDIRS += qd
endif

if ENABLE_DOUBLE
lib_LTLIBRARIES += liblin_double.la liblin_override_double.la
SUBDIRS += double
endif

if ENABLE__FLOAT64X
lib_LTLIBRARIES += liblin__Float64x.la liblin_override__Float64x.la
SUBDIRS += _Float64x
endif

if ENABLE__FLOAT128
lib_LTLIBRARIES += liblin__Float128.la liblin_override__Float128.la
SUBDIRS += _Float128
endif

if ENABLE_MPFR
liblin_override_mpfr_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -I$(top_srcdir)/mplapack/test/fable
liblin_override_mpfr_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp
liblin_mpfr_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -D___MPLAPACK_INTERNAL___
liblin_mpfr_la_SOURCES = $(LIN_SOURCES)
if IS_MACOSX
liblin_mpfr_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_GMP
liblin_override_gmp_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___ -I$(top_srcdir)/mplapack/test/fable
liblin_override_gmp_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp
liblin_gmp_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___ -D___MPLAPACK_INTERNAL___
liblin_gmp_la_SOURCES = $(LIN_SOURCES)
if IS_MACOSX
liblin_gmp_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_QD
liblin_override_qd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___ -I$(top_srcdir)/mplapack/test/fable
liblin_override_qd_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp
liblin_qd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___ -D___MPLAPACK_INTERNAL___
liblin_qd_la_SOURCES = $(LIN_SOURCES)
if IS_MACOSX
liblin_qd_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DD
liblin_override_dd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___ -I$(top_srcdir)/mplapack/test/fable
liblin_override_dd_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp
liblin_dd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___ -D___MPLAPACK_INTERNAL___
liblin_dd_la_SOURCES = $(LIN_SOURCES)
if IS_MACOSX
liblin_dd_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DOUBLE
liblin_override_double_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DOUBLE___ -I$(top_srcdir)/mplapack/test/fable
liblin_override_double_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp
liblin_double_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -D___MPLAPACK_BUILD_WITH_DOUBLE___ -D___MPLAPACK_INTERNAL___
liblin_double_la_SOURCES = $(LIN_SOURCES)
if IS_MACOSX
liblin_double_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE__FLOAT64X
liblin_override__Float64x_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT64X___ -I$(top_srcdir)/mplapack/test/fable
liblin_override__Float64x_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp
liblin__Float64x_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -D___MPLAPACK_BUILD_WITH__FLOAT64X___ -D___MPLAPACK_INTERNAL___
liblin__Float64x_la_SOURCES = $(LIN_SOURCES)
if IS_MACOSX
liblin__Float64x_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE__FLOAT128
liblin_override__Float128_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT128___ -I$(top_srcdir)/mplapack/test/fable
liblin_override__Float128_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp
liblin__Float128_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -D___MPLAPACK_BUILD_WITH__FLOAT128___ -D___MPLAPACK_INTERNAL___
liblin__Float128_la_SOURCES = $(LIN_SOURCES)
if IS_MACOSX
liblin__Float128_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

LIN_SOURCES = \
common/xlaenv.cpp \
common/Rdrvpp.cpp \
common/Cdrvpp.cpp \
common/Alaesm.cpp \
common/chkxer.cpp \
common/Claipd.cpp \
common/Alasum.cpp \
common/Rgennd.cpp \
common/Alasvm.cpp \
common/Cgennd.cpp \
common/Rget06.cpp \
common/Rgeqrs.cpp \
common/Cgeqrs.cpp \
common/Rgeqls.cpp \
common/Cgeqls.cpp \
common/Rget03.cpp \
common/Rqrt11.cpp \
common/Rgelqs.cpp \
common/Rgerqs.cpp \
common/Cqrt11.cpp \
common/Rptt02.cpp \
common/Rget04.cpp \
common/Cgelqs.cpp \
common/Cgerqs.cpp \
common/Cptt02.cpp \
common/Cpot02.cpp \
common/Rpot02.cpp \
common/Rrzt02.cpp \
common/Rppt02.cpp \
common/Crzt02.cpp \
common/Rgtt02.cpp \
common/Cget03.cpp \
common/Cgtt02.cpp \
common/Ctpt06.cpp \
common/Csyt02.cpp \
common/Rpbt02.cpp \
common/Rlaord.cpp \
common/Rpot06.cpp \
common/Cspt02.cpp \
common/Cppt02.cpp \
common/Rqrt13.cpp \
common/Ctrt06.cpp \
common/Cqrt13.cpp \
common/Cpbt02.cpp \
common/Rerrqp.cpp \
common/Ctbt06.cpp \
common/Rget02.cpp \
common/Rptt01.cpp \
common/Cptt01.cpp \
common/Rtpt02.cpp \
common/Rtpt06.cpp \
common/Rerrtz.cpp \
common/Cget04.cpp \
common/Ctpt02.cpp \
common/Rtrt06.cpp \
common/Rtrt02.cpp \
common/Rrzt01.cpp \
common/Rtbt06.cpp \
common/Rget08.cpp \
common/Rtbt02.cpp \
common/Crzt01.cpp \
common/Ctrt02.cpp \
common/Cerrqp.cpp \
common/Cget02.cpp \
common/Rqrt16.cpp \
common/Ctbt02.cpp \
common/Cpot06.cpp \
common/Cerrtz.cpp \
common/Rtplqt.cpp \
common/Cqrt16.cpp \
common/Rgbt02.cpp \
common/Rqpt01.cpp \
common/Cqpt01.cpp \
common/Cget08.cpp \
common/Rqrt02.cpp \
common/Rspt01.cpp \
common/Rlqt02.cpp \
common/Rpot03.cpp \
common/Cgbt02.cpp \
common/Rerrps.cpp \
common/Rsyt01.cpp \
common/Cerrps.cpp \
common/Rtrt01.cpp \
common/Ctrt01.cpp \
common/Rsyt01_rook.cpp \
common/Cspt01.cpp \
common/Cqrt02.cpp \
common/Clqt02.cpp \
common/Cpot03.cpp \
common/Csyt01.cpp \
common/Rsyt01_aa.cpp \
common/Csyt01_aa.cpp \
common/Csyt03.cpp \
common/Csyt01_rook.cpp \
common/Rqrt01.cpp \
common/Rqrt01p.cpp \
common/Rlqt01.cpp \
common/Rtpt01.cpp \
common/Ctpt01.cpp \
common/Cqrt01.cpp \
common/Cqrt01p.cpp \
common/Rlahilb.cpp \
common/Clqt01.cpp \
common/Rsyt01_3.cpp \
common/Clatb5.cpp \
common/Rlatb5.cpp \
common/Rerrorhr_col.cpp \
common/Cerrunhr_col.cpp \
common/Cchklqt.cpp \
common/Rchklqt.cpp \
common/Rchkqrt.cpp \
common/Cchkqrt.cpp \
common/Csyt01_3.cpp \
common/Rppt01.cpp \
common/Rlaptm.cpp \
common/Rqrt17.cpp \
common/Rget01.cpp \
common/Ctrt03.cpp \
common/Cqrt17.cpp \
common/Rtrt03.cpp \
common/Cget01.cpp \
common/Chet01.cpp \
common/Chet01_rook.cpp \
common/Rrqt02.cpp \
common/Cchkqrtp.cpp \
common/Rqlt02.cpp \
common/Ctbt03.cpp \
common/Ctpt03.cpp \
common/Cchklqtp.cpp \
common/Rchklqtp.cpp \
common/Rchkqrtp.cpp \
common/Crqt02.cpp \
common/Rtbt03.cpp \
common/Cqlt02.cpp \
common/Rtpt03.cpp \
common/Rpot01.cpp \
common/Chpt01.cpp \
common/Rrqt01.cpp \
common/Rqlt01.cpp \
common/Chet01_3.cpp \
common/Crqt01.cpp \
common/Cqlt01.cpp \
common/Cqrt12.cpp \
common/Rget07.cpp \
common/Rptt05.cpp \
common/Rgbt01.cpp \
common/Rqrt12.cpp \
common/Rgbt05.cpp \
common/Rqrt03.cpp \
common/Rlqt03.cpp \
common/Rppt03.cpp \
common/Cgbt01.cpp \
common/Cqrt03.cpp \
common/Rpbt01.cpp \
common/Rgtt01.cpp \
common/Cgtt01.cpp \
common/Cptt05.cpp \
common/Cpot01.cpp \
common/Rpot05.cpp \
common/Cget07.cpp \
common/Cgbt05.cpp \
common/Rppt05.cpp \
common/Cppt03.cpp \
common/Rerrlqt.cpp \
common/Rpbt05.cpp \
common/Cppt01.cpp \
common/Cerrlqt.cpp \
common/Rrqt03.cpp \
common/Rqlt03.cpp \
common/Cpot05.cpp \
common/Crqt03.cpp \
common/Rqrt14.cpp \
common/Cqlt03.cpp \
common/Cqrt14.cpp \
common/Cpbt01.cpp \
common/Cppt05.cpp \
common/Rerrqrt.cpp \
common/Cpbt05.cpp \
common/Cerrqrt.cpp \
common/Rgtt05.cpp \
common/Alareq.cpp \
common/Rqrt15.cpp \
common/Rtrt05.cpp \
common/Rerrlqtp.cpp \
common/Cqrt15.cpp \
common/Cgtt05.cpp \
common/Cchktsqr.cpp \
common/Rchktsqr.cpp \
common/Cerrlqtp.cpp \
common/Clahilb.cpp \
common/Rtbt05.cpp \
common/Rerrls.cpp \
common/Rtpt05.cpp \
common/Rerrqrtp.cpp \
common/Ctrt05.cpp \
common/Cerrqrtp.cpp \
common/Cerrls.cpp \
common/Rpst01.cpp \
common/Claptm.cpp \
common/Ctbt05.cpp \
common/Ctpt05.cpp \
common/Clatsp.cpp \
common/Rerrgt.cpp \
common/Cspt03.cpp \
common/Cerrtsqr.cpp \
common/Rlarhs.cpp \
common/Rerrtsqr.cpp \
common/Cerrgt.cpp \
common/Csbmv.cpp \
common/Cpst01.cpp \
common/Rdrvrf2.cpp \
common/Cdrvrf2.cpp \
common/Clatsy.cpp \
common/Rchktz.cpp \
common/Cchktz.cpp \
common/Rerrrfp.cpp \
common/Cerrrfp.cpp \
common/Clarhs.cpp \
common/Rdrvrf1.cpp \
common/Cdrvrf1.cpp \
common/Rerrlq.cpp \
common/Rerrql.cpp \
common/Rerrrq.cpp \
common/Cerrlq.cpp \
common/Cerrql.cpp \
common/Cerrrq.cpp \
common/Cchkunhr_col.cpp \
common/Rchkorhr_col.cpp \
common/Rlqt05.cpp \
common/Rqrt05.cpp \
common/Clqt05.cpp \
common/Cdrvls.cpp \
common/Rdrvls.cpp \
common/icopy.cpp \
common/Rtsqr01.cpp \
common/Cqrt05.cpp \
common/Clqt04.cpp \
common/Cqrt04.cpp \
common/Rlqt04.cpp \
common/Rqrt04.cpp \
common/Cunhr_col01.cpp \
common/Cunhr_col02.cpp \
common/Rorhr_col01.cpp \
common/Rorhr_col02.cpp \
common/Cchkps.cpp \
common/Cchkq3.cpp \
common/Rchkps.cpp \
common/Rchkql.cpp \
common/Cchklq.cpp \
common/Cchkql.cpp \
common/Cchkrq.cpp \
common/Rchkrq.cpp \
common/Ctsqr01.cpp \
common/Rchksy_rook.cpp \
common/Rchksy_rk.cpp \
common/Cchkhe_rk.cpp \
common/Rchkq3.cpp \
common/Rerrqr.cpp \
common/Cerrqr.cpp \
common/Rchklq.cpp \
common/Rchkqr.cpp \
common/Cchkqr.cpp \
common/Chet01_aa.cpp \
common/Rdrvrf4.cpp \
common/Cdrvrf4.cpp \
common/Rlatb4.cpp \
common/Rerrge.cpp \
common/Clatb4.cpp \
common/Cerrge.cpp \
common/Rdrvsy_aa.cpp \
common/Rerrpo.cpp \
common/Cdrvsy_aa.cpp \
common/Rdrvrf3.cpp \
common/Cdrvhe_aa.cpp \
common/Rchkpp.cpp \
common/Cerrpo.cpp \
common/Cchkhe_rook.cpp\
common/Cchksy_rook.cpp\
common/Rdrvsy_aa_2stage.cpp \
common/Cchkpp.cpp \
common/Cdrvsy_aa_2stage.cpp \
common/Cdrvhe_aa_2stage.cpp \
common/Alahd.cpp \
common/Alaerh.cpp \
common/Cdrvrf3.cpp \
common/Rchkpo.cpp \
common/Rchkpt.cpp \
common/Rdrvrfp.cpp \
common/Cchkpo.cpp \
common/Cerrtr.cpp \
common/Rerrtr.cpp \
common/Cdrvrfp.cpp \
common/Cdrvhe_rook.cpp \
common/Rlavsy.cpp \
common/Rdrvsy_rook.cpp \
common/Cchkpt.cpp \
common/Clavsy.cpp \
common/Rlavsp.cpp \
common/Cdrvhe_rk.cpp \
common/Rdrvsy_rk.cpp \
common/Rchksp.cpp \
common/Clavsp.cpp \
common/Rchkeq.cpp \
common/Rchktp.cpp \
common/Cchktp.cpp \
common/Rlattp.cpp \
common/Clattp.cpp \
common/Aladhd.cpp \
common/Clqt03.cpp \
common/Rchkgt.cpp \
common/Cchkgt.cpp \
common/Rdrvgb.cpp \
common/Cdrvgb.cpp \
common/Cchkhp.cpp \
common/Cdrvsy_rook.cpp \
common/Cdrvsy_rk.cpp \
common/Rchksy_aa_2stage.cpp \
common/Rchksy_aa.cpp \
common/Cchksy_aa_2stage.cpp \
common/Cchksy_aa.cpp \
common/Cchksp.cpp \
common/Clavhp.cpp \
common/Cchkhe_aa.cpp \
common/Clavhe.cpp \
common/Cchkeq.cpp \
common/Cchkhe_aa_2stage.cpp \
common/Rlavsy_rook.cpp \
common/Clavsy_rook.cpp \
common/Rchkge.cpp \
common/Rdrvpt.cpp \
common/Cchkge.cpp \
common/Cdrvpt.cpp \
common/Rchktr.cpp \
common/Cchktr.cpp \
common/Rchkgb.cpp \
common/Rlattr.cpp \
common/Clattr.cpp \
common/Rdrvpb.cpp \
common/Cdrvpb.cpp \
common/Rdrvge.cpp \
common/Cdrvge.cpp \
common/Cchkgb.cpp \
common/Rerrvx.cpp \
common/Cerrvx.cpp \
common/Cchksy_rk.cpp \
common/Clavhe_rook.cpp \
common/Rlattb.cpp \
common/Clattb.cpp \
common/Rchkpb.cpp \
common/Cchkpb.cpp \
common/Rerrsy.cpp \
common/Cerrhe.cpp \
common/Cchktb.cpp \
common/Rchktb.cpp \
common/Cchkhe.cpp \
common/Rchksy.cpp \
common/Cchksy.cpp \
common/Cerrsy.cpp \
common/Rdrvsp.cpp \
common/Cdrvsp.cpp \
common/Rdrvgt.cpp \
common/Cdrvgt.cpp \
common/Rdrvsy.cpp \
common/Cdrvhp.cpp \
common/Cdrvhe.cpp \
common/Cdrvsy.cpp \
common/Rdrvpo.cpp \
common/Cdrvpo.cpp 
