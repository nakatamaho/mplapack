SUBDIRS = .
lib_LTLIBRARIES = 

mplapack_eig_DATA = Cbak.in Cbal.in Cbal_double.in Cbb.in Cec.in Ced.in Cgbak.in Cgbal.in Cgd.in Cgg.in \
Csb.in Csg.in Rbak.in Rbal.in Rbal_double.in Rbb.in Rec.in Red.in Rgbak.in Rgbal.in Rgd.in Rgg.in Rsb.in Rsg.in \
csd.in glm.in gqr.in gsv.in lse.in nep.in se2.in sep.in svd.in
mplapack_eigdir = $(prefix)/lib/$(target)/mplapack/test/eig

if ENABLE_MPFR
lib_LTLIBRARIES += libeig_mpfr.la libeig_override_mpfr.la
SUBDIRS += mpfr
endif

if ENABLE_GMP
lib_LTLIBRARIES += libeig_gmp.la libeig_override_gmp.la
SUBDIRS += gmp
endif

if ENABLE_DD
lib_LTLIBRARIES += libeig_dd.la libeig_override_dd.la
SUBDIRS += dd
endif

if ENABLE_QD
lib_LTLIBRARIES += libeig_qd.la libeig_override_qd.la
SUBDIRS += qd
endif

if ENABLE_DOUBLE
lib_LTLIBRARIES += libeig_double.la libeig_override_double.la
SUBDIRS += double
endif

if ENABLE__FLOAT64X
lib_LTLIBRARIES += libeig__Float64x.la libeig_override__Float64x.la
SUBDIRS += _Float64x
endif

if ENABLE__FLOAT128
lib_LTLIBRARIES += libeig__Float128.la libeig_override__Float128.la
SUBDIRS += _Float128
endif

if ENABLE_MPFR
libeig_mpfr_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -I$(LAPACK_INCLUDEDIR) -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -D___MPLAPACK_INTERNAL___
libeig_mpfr_la_SOURCES = $(EIG_SOURCES)
libeig_override_mpfr_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -I$(top_srcdir)/mplapack/test/fable
libeig_override_mpfr_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp
if IS_MACOSX
libeig_mpfr_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_GMP
libeig_gmp_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -I$(LAPACK_INCLUDEDIR) -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___ -D___MPLAPACK_INTERNAL___
libeig_gmp_la_SOURCES = $(EIG_SOURCES)
libeig_override_gmp_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(top_srcdir)/mplapack/test/fable -D___MPLAPACK_BUILD_WITH_GMP___
libeig_override_gmp_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp
if IS_MACOSX
libeig_gmp_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_QD
libeig_qd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -I$(LAPACK_INCLUDEDIR) -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___ -D___MPLAPACK_INTERNAL___
libeig_qd_la_SOURCES = $(EIG_SOURCES)
libeig_override_qd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___
libeig_override_qd_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp
if IS_MACOSX
libeig_qd_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DD
libeig_dd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -I$(LAPACK_INCLUDEDIR) -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___ -D___MPLAPACK_INTERNAL___
libeig_dd_la_SOURCES = $(EIG_SOURCES)
libeig_override_dd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___
libeig_override_dd_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp
if IS_MACOSX
libeig_dd_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DOUBLE
libeig_double_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -I$(LAPACK_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DOUBLE___ -D___MPLAPACK_INTERNAL___
libeig_double_la_SOURCES = $(EIG_SOURCES)
libeig_override_double_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -D___MPLAPACK_BUILD_WITH_DOUBLE___
libeig_override_double_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp
if IS_MACOSX
libeig_double_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE__FLOAT64X
libeig__Float64x_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -I$(LAPACK_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH__FLOAT64X___ -D___MPLAPACK_INTERNAL___
libeig__Float64x_la_SOURCES = $(EIG_SOURCES)
libeig_override__Float64x_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -D___MPLAPACK_BUILD_WITH__FLOAT64X___
libeig_override__Float64x_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp
if IS_MACOSX
libeig__Float64x_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE__FLOAT128
libeig__Float128_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -I$(LAPACK_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH__FLOAT128___ -D___MPLAPACK_INTERNAL___
libeig__Float128_la_SOURCES = $(EIG_SOURCES)
libeig_override__Float128_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mplapack/test/fable -D___MPLAPACK_BUILD_WITH__FLOAT128___
libeig_override__Float128_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp
if IS_MACOSX
libeig__Float128_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

EIG_SOURCES = \
common/Cslect.cpp \
common/Rslect.cpp \
common/Rlctsx.cpp \
common/Rlacsg.cpp \
common/Clacsg.cpp \
common/Rlasum.cpp \
common/Rlctes.cpp \
common/Alasmg.cpp \
common/Alasum.cpp \
common/Alasvm.cpp \
common/Rlarfy.cpp \
common/Clarfy.cpp \
common/Clctes.cpp \
common/Rsxt1.cpp \
common/Clctsx.cpp \
common/Rget10.cpp \
common/Cget10.cpp \
common/Rlsets.cpp \
common/Clsets.cpp \
common/Rbdt02.cpp \
common/Cbdt02.cpp \
common/Rbdt05.cpp \
common/Rget02.cpp \
common/Cget02.cpp \
common/Cbdt05.cpp \
common/Rglmts.cpp \
common/Rhst01.cpp \
common/Cglmts.cpp \
common/Chst01.cpp \
common/Rstect.cpp \
common/Rsgt01.cpp \
common/Rget54.cpp \
common/Rsyt22.cpp \
common/Csgt01.cpp \
common/Cget54.cpp \
common/Rsvdct.cpp \
common/Chet22.cpp \
common/Rstt21.cpp \
common/Rlafts.cpp \
common/Cstt21.cpp \
common/Rort01.cpp \
common/Rsvdch.cpp \
common/Cunt01.cpp \
common/Rbdt04.cpp \
common/Cchkec.cpp \
common/Rstech.cpp \
common/Rort03.cpp \
common/Rget51.cpp \
common/Cget51.cpp \
common/Rstt22.cpp \
common/Rsbt21.cpp \
common/Cstt22.cpp \
common/Chbt21.cpp \
common/Cget52.cpp \
common/Rbdt03.cpp \
common/Cbdt03.cpp \
common/Rget53.cpp \
common/Cunt03.cpp \
common/Rckgsv.cpp \
common/Cget22.cpp \
common/Rgqrts.cpp \
common/Rgrqts.cpp \
common/Cckgsv.cpp \
common/Rbdt01.cpp \
common/Cgqrts.cpp \
common/Cgrqts.cpp \
common/Cbdt01.cpp \
common/Rckglm.cpp \
common/Cckglm.cpp \
common/Rlatb9.cpp \
common/Rckcsd.cpp \
common/Cckcsd.cpp \
common/Rcklse.cpp \
common/Ccklse.cpp \
common/Rchkec.cpp \
common/Rgsvts3.cpp \
common/Csbmv.cpp \
common/Rlarhs.cpp \
common/Cgsvts3.cpp \
common/Rget22.cpp \
common/Rget52.cpp \
common/Rget33.cpp \
common/Rsyt21.cpp \
common/Chet21.cpp \
common/Rspt21.cpp \
common/Clarhs.cpp \
common/Cerrec.cpp \
common/Chpt21.cpp \
common/Cerrbd.cpp \
common/Rerrec.cpp \
common/Cckgqr.cpp \
common/Rckgqr.cpp \
common/Clatm4.cpp \
common/Alahdg.cpp \
common/Rlatm4.cpp \
common/Rerrbd.cpp \
common/Rcsdts.cpp \
common/Rget39.cpp \
common/Ccsdts.cpp \
common/Cerrhs.cpp \
common/Rerrhs.cpp \
common/Rchkbb.cpp \
common/Cchkbb.cpp \
common/Rget32.cpp \
common/Rget34.cpp \
common/Cget23.cpp \
common/Cerred.cpp \
common/Rget23.cpp \
common/Rerred.cpp \
common/Rchksb.cpp \
common/Cchkhb.cpp \
common/Rget31.cpp \
common/Cget24.cpp \
common/Cchkbd.cpp \
common/Cdrves.cpp \
common/Cdrvev.cpp \
common/Rdrvev.cpp \
common/Rchksb2stg.cpp \
common/Rget24.cpp \
common/Cchkhb2stg.cpp \
common/Rdrgev.cpp \
common/Rdrgev3.cpp \
common/Cdrgev.cpp \
common/Cdrgev3.cpp \
common/Rdrves.cpp \
common/Rlahd2.cpp \
common/Cdrges.cpp \
common/Cdrges3.cpp \
common/Cchkhs.cpp \
common/Rchkhs.cpp \
common/Rdrges.cpp \
common/Rdrges3.cpp \
common/Rchkgg.cpp \
common/Cchkgg.cpp \
common/Rdrvbd.cpp \
common/Rchkbd.cpp \
common/Rdrvsg.cpp \
common/Cdrvbd.cpp \
common/Cdrvsg.cpp \
common/Cerrgg.cpp \
common/Cerrst.cpp \
common/Rerrgg.cpp \
common/Rerrst.cpp \
common/Rdrvsg2stg.cpp \
common/Cdrvsg2stg.cpp \
common/Rchkst.cpp \
common/Cchkst.cpp \
common/Rchkst2stg.cpp \
common/Cchkst2stg.cpp \
common/Cdrvst.cpp \
common/Cdrvst2stg.cpp \
common/Rdrvst.cpp \
common/Rdrvst2stg.cpp \
common/Cchkbk.cpp \
common/Rchkbk.cpp \
common/Cget36.cpp \
common/Rchkgl.cpp \
common/Cchkgl.cpp \
common/Rget36.cpp \
common/Rchkgk.cpp \
common/Cchkgk.cpp \
common/Cget35.cpp \
common/Rget35.cpp \
common/Rget38.cpp \
common/Cget38.cpp \
common/Cget37.cpp \
common/Rget37.cpp \
common/Cdrgvx.cpp \
common/Rget40.cpp \
common/Rdrgvx.cpp \
common/Rchkbl.cpp \
common/Cchkbl.cpp \
common/Rdrvvx.cpp \
common/Cdrvvx.cpp \
common/Cdrvsx.cpp \
common/Rdrvsx.cpp \
common/Cdrgsx.cpp \
common/Rdrgsx.cpp \
common/Alareq.cpp \
common/Alarqg.cpp \
common/iMlaenv.cpp \
common/chkxer.cpp \
common/xlaenv.cpp 
