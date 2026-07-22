SUBDIRS = .
lib_LTLIBRARIES = 
# -----------------------------------------------------------------------------
# Parallel 'make -j check' across independent backend subdirectories.
#
# Automake's default recursive rule (check-recursive) walks SUBDIRS serially
# even when GNU make is invoked with -j.  Backend test suites are independent,
# so we expose each backend check as a separate make target.  GNU make can
# then schedule them in parallel while the jobserver is correctly inherited
# by the sub-makes via '+$(MAKE)'.
#
# NOTE: We override only 'check-recursive'.  All other recursive targets
# (all/install/clean/dist/...) keep Automake's standard behavior.
# -----------------------------------------------------------------------------
CHECK_BACKENDS =

mplapack_eig_DATA = Cbak.in Cbal.in Cbb.in Cec.in Ced.in Cgbak.in Cgbal.in Cgd.in Cgg.in \
Csb.in Csg.in Rbak.in Rbal.in Rbb.in Rec.in Red.in Rgbak.in Rgbal.in Rgd.in Rgg.in Rsb.in Rsg.in \
csd.in glm.in gqr.in gsv.in lse.in nep.in se2.in se2_opt.in sep.in svd.in
mplapack_eigdir = $(prefix)/lib/$(target)/mplapack/test/eig

if ENABLE_MPFR
lib_LTLIBRARIES += libeig_mpfr.la libeig_override_mpfr.la
SUBDIRS += mpfr
CHECK_BACKENDS += mpfr
endif

if ENABLE_GMP
lib_LTLIBRARIES += libeig_gmp.la libeig_override_gmp.la
SUBDIRS += gmp
CHECK_BACKENDS += gmp
endif

if ENABLE_DD
lib_LTLIBRARIES += libeig_dd.la libeig_override_dd.la
SUBDIRS += dd
CHECK_BACKENDS += dd
endif

if ENABLE_QD
lib_LTLIBRARIES += libeig_qd.la libeig_override_qd.la
SUBDIRS += qd
CHECK_BACKENDS += qd
endif

if ENABLE_DOUBLE
lib_LTLIBRARIES += libeig_double.la libeig_override_double.la
SUBDIRS += double
CHECK_BACKENDS += double
endif

if ENABLE_BINARY80
lib_LTLIBRARIES += libeig_binary80.la libeig_override_binary80.la
SUBDIRS += binary80
CHECK_BACKENDS += binary80
endif

if ENABLE_BINARY128
lib_LTLIBRARIES += libeig_binary128.la libeig_override_binary128.la
SUBDIRS += binary128
CHECK_BACKENDS += binary128
endif

# -----------------------------------------------------------------------------
# Parallel backend checks (GNU make jobserver aware)
# -----------------------------------------------------------------------------
CHECK_BACKEND_TARGETS = $(patsubst %,check-backend-%,$(CHECK_BACKENDS))

.PHONY: check-parallel

# Run local checks first (check-am builds this directory and runs any local
# tests), then run each backend check in parallel.
check-parallel: check-am $(CHECK_BACKEND_TARGETS)

# Each backend is an independent recursive make. The order-only prerequisite
# ensures that 'check-am' (and thus 'all-am') has completed before any backend
# starts, without forcing rebuilds when the directory list changes.
# NOTE: Do not declare check-backend-* targets as phony.
# GNU make does not apply pattern rules to phony targets.
check-backend-%: | check-am
	+$(MAKE) $(AM_MAKEFLAGS) -C $* check

# Override Automake's default serial recursive driver for 'make check' in this
# directory. Other recursive targets (all/install/clean/dist/...) keep the
# standard Automake behavior.
check-recursive: check-parallel
	@:

if ENABLE_MPFR
libeig_mpfr_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(LAPACK_INCLUDEDIR) -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -D___MPLAPACK_INTERNAL___
libeig_mpfr_la_SOURCES = $(EIG_SOURCES)
libeig_mpfr_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
libeig_override_mpfr_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -I$(top_srcdir)/fable
libeig_override_mpfr_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
libeig_override_mpfr_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
libeig_mpfr_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
libeig_override_mpfr_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_GMP
libeig_gmp_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(LAPACK_INCLUDEDIR) -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___ -D___MPLAPACK_INTERNAL___
libeig_gmp_la_SOURCES = $(EIG_SOURCES)
libeig_gmp_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
libeig_override_gmp_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(top_srcdir)/fable -D___MPLAPACK_BUILD_WITH_GMP___
libeig_override_gmp_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
libeig_override_gmp_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
libeig_gmp_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
libeig_override_gmp_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_QD
libeig_qd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(LAPACK_INCLUDEDIR) -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___ @qd_extra_cxxflags@ -D___MPLAPACK_INTERNAL___
libeig_qd_la_SOURCES = $(EIG_SOURCES)
libeig_qd_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
libeig_override_qd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___ @qd_extra_cxxflags@
libeig_override_qd_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
libeig_override_qd_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
libeig_qd_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
libeig_override_qd_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DD
libeig_dd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(LAPACK_INCLUDEDIR) -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___ @qd_extra_cxxflags@ -D___MPLAPACK_INTERNAL___
libeig_dd_la_SOURCES = $(EIG_SOURCES)
libeig_dd_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
libeig_override_dd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___ @qd_extra_cxxflags@
libeig_override_dd_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
libeig_override_dd_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
libeig_dd_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
libeig_override_dd_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DOUBLE
libeig_double_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(LAPACK_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DOUBLE___ -D___MPLAPACK_INTERNAL___
libeig_double_la_SOURCES = $(EIG_SOURCES)
libeig_double_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
libeig_override_double_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -D___MPLAPACK_BUILD_WITH_DOUBLE___
libeig_override_double_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
libeig_override_double_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
libeig_double_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
libeig_override_double_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_BINARY80
libeig_binary80_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(LAPACK_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_BINARY80___ -D___MPLAPACK_INTERNAL___
libeig_binary80_la_SOURCES = $(EIG_SOURCES)
libeig_binary80_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
libeig_override_binary80_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -D___MPLAPACK_BUILD_WITH_BINARY80___
libeig_override_binary80_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
libeig_override_binary80_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
libeig_binary80_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
libeig_override_binary80_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_BINARY128
libeig_binary128_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(LAPACK_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_BINARY128___ -D___MPLAPACK_INTERNAL___
libeig_binary128_la_SOURCES = $(EIG_SOURCES)
libeig_binary128_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
libeig_override_binary128_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -D___MPLAPACK_BUILD_WITH_BINARY128___
libeig_override_binary128_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
libeig_override_binary128_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
libeig_binary128_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
libeig_override_binary128_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

EIG_SOURCES = \
common/Alahdg.cpp \
common/Alareq.cpp \
common/Alarqg.cpp \
common/Alasmg.cpp \
common/Alasum.cpp \
common/Alasvm.cpp \
common/Cbdt01.cpp \
common/Cbdt02.cpp \
common/Cbdt03.cpp \
common/Cbdt05.cpp \
common/Cchkbb.cpp \
common/Cchkbd.cpp \
common/Cchkbk.cpp \
common/Cchkbl.cpp \
common/Cchkec.cpp \
common/Cchkgg.cpp \
common/Cchkgk.cpp \
common/Cchkgl.cpp \
common/Cchkhb.cpp \
common/Cchkhb2stg.cpp \
common/Cchkhs.cpp \
common/Cchkst.cpp \
common/Cchkst2stg.cpp \
common/Cckcsd.cpp \
common/Cckglm.cpp \
common/Cckgqr.cpp \
common/Cckgsv.cpp \
common/Ccklse.cpp \
common/Ccsdts.cpp \
common/Cdrges.cpp \
common/Cdrges3.cpp \
common/Cdrgev.cpp \
common/Cdrgev3.cpp \
common/Cdrgsx.cpp \
common/Cdrgvx.cpp \
common/Cdrvbd.cpp \
common/Cdrves.cpp \
common/Cdrvev.cpp \
common/Cdrvsg.cpp \
common/Cdrvsg2stg.cpp \
common/Cdrvst.cpp \
common/Cdrvst2stg.cpp \
common/Cdrvsx.cpp \
common/Cdrvvx.cpp \
common/Cerrbd.cpp \
common/Cerrec.cpp \
common/Cerred.cpp \
common/Cerrgg.cpp \
common/Cerrhs.cpp \
common/Cerrst.cpp \
common/Cget02.cpp \
common/Cget10.cpp \
common/Cget22.cpp \
common/Cget23.cpp \
common/Cget24.cpp \
common/Cget35.cpp \
common/Cget36.cpp \
common/Cget37.cpp \
common/Cget38.cpp \
common/Cget51.cpp \
common/Cget52.cpp \
common/Cget54.cpp \
common/Cglmts.cpp \
common/Cgqrts.cpp \
common/Cgrqts.cpp \
common/Cgsvts3.cpp \
common/Chbt21.cpp \
common/Chet21.cpp \
common/Chet22.cpp \
common/Chkxer.cpp \
common/Chpt21.cpp \
common/Chst01.cpp \
common/Clarfy.cpp \
common/Clarhs.cpp \
common/Clatm4.cpp \
common/Clctes.cpp \
common/Clctsx.cpp \
common/Clsets.cpp \
common/Csbmv.cpp \
common/Csgt01.cpp \
common/Cslect.cpp \
common/Cstt21.cpp \
common/Cstt22.cpp \
common/Csyl01.cpp \
common/Cunt01.cpp \
common/Cunt03.cpp \
common/Rbdt01.cpp \
common/Rbdt02.cpp \
common/Rbdt03.cpp \
common/Rbdt04.cpp \
common/Rbdt05.cpp \
common/Rchkbb.cpp \
common/Rchkbd.cpp \
common/Rchkbk.cpp \
common/Rchkbl.cpp \
common/Rchkec.cpp \
common/Rchkgg.cpp \
common/Rchkgk.cpp \
common/Rchkgl.cpp \
common/Rchkhs.cpp \
common/Rchksb.cpp \
common/Rchksb2stg.cpp \
common/Rchkst.cpp \
common/Rchkst2stg.cpp \
common/Rckcsd.cpp \
common/Rckglm.cpp \
common/Rckgqr.cpp \
common/Rckgsv.cpp \
common/Rcklse.cpp \
common/Rcsdts.cpp \
common/Rdrges.cpp \
common/Rdrges3.cpp \
common/Rdrgev.cpp \
common/Rdrgev3.cpp \
common/Rdrgsx.cpp \
common/Rdrgvx.cpp \
common/Rdrvbd.cpp \
common/Rdrves.cpp \
common/Rdrvev.cpp \
common/Rdrvsg.cpp \
common/Rdrvsg2stg.cpp \
common/Rdrvst.cpp \
common/Rdrvst2stg.cpp \
common/Rdrvsx.cpp \
common/Rdrvvx.cpp \
common/Rerrbd.cpp \
common/Rerrec.cpp \
common/Rerred.cpp \
common/Rerrgg.cpp \
common/Rerrhs.cpp \
common/Rerrst.cpp \
common/Rget02.cpp \
common/Rget10.cpp \
common/Rget22.cpp \
common/Rget23.cpp \
common/Rget24.cpp \
common/Rget31.cpp \
common/Rget32.cpp \
common/Rget33.cpp \
common/Rget34.cpp \
common/Rget35.cpp \
common/Rget36.cpp \
common/Rget37.cpp \
common/Rget38.cpp \
common/Rget39.cpp \
common/Rget40.cpp \
common/Rget51.cpp \
common/Rget52.cpp \
common/Rget53.cpp \
common/Rget54.cpp \
common/Rglmts.cpp \
common/Rgqrts.cpp \
common/Rgrqts.cpp \
common/Rgsvts3.cpp \
common/Rhst01.cpp \
common/Rlafts.cpp \
common/Rlahd2.cpp \
common/Rlarfy.cpp \
common/Rlarhs.cpp \
common/Rlasum.cpp \
common/Rlatb9.cpp \
common/Rlatm4.cpp \
common/Rlctes.cpp \
common/Rlctsx.cpp \
common/Rlsets.cpp \
common/Rort01.cpp \
common/Rort03.cpp \
common/Rsbt21.cpp \
common/Rsgt01.cpp \
common/Rslect.cpp \
common/Rspt21.cpp \
common/Rstech.cpp \
common/Rstect.cpp \
common/Rstt21.cpp \
common/Rstt22.cpp \
common/Rsvdch.cpp \
common/Rsvdct.cpp \
common/Rsxt1.cpp \
common/Rsyl01.cpp \
common/Rsyt21.cpp \
common/Rsyt22.cpp \
common/iMlaenv.cpp

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_eigdir) $(prefix)
endif
