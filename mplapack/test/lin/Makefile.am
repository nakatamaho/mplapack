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

mplapack_lin_DATA = Ctest.in Ctest_rfp.in Rtest.in Rtest_rfp.in
mplapack_lindir = $(prefix)/lib/$(target)/mplapack/test/lin

if ENABLE_MPFR
lib_LTLIBRARIES += liblin_mpfr.la liblin_override_mpfr.la liblin_mpfr_opt.la liblin_override_mpfr_opt.la
SUBDIRS += mpfr
CHECK_BACKENDS += mpfr
endif

if ENABLE_GMP
lib_LTLIBRARIES += liblin_gmp.la liblin_override_gmp.la liblin_gmp_opt.la liblin_override_gmp_opt.la
SUBDIRS += gmp
CHECK_BACKENDS += gmp
endif

if ENABLE_DD
lib_LTLIBRARIES += liblin_dd.la liblin_override_dd.la liblin_dd_opt.la liblin_override_dd_opt.la
SUBDIRS += dd
CHECK_BACKENDS += dd
endif

if ENABLE_QD
lib_LTLIBRARIES += liblin_qd.la liblin_override_qd.la liblin_qd_opt.la liblin_override_qd_opt.la
SUBDIRS += qd
CHECK_BACKENDS += qd
endif

if ENABLE_DOUBLE
lib_LTLIBRARIES += liblin_double.la liblin_override_double.la liblin_double_opt.la liblin_override_double_opt.la
SUBDIRS += double
CHECK_BACKENDS += double
endif

if ENABLE_BINARY80
lib_LTLIBRARIES += liblin_binary80.la liblin_override_binary80.la liblin_binary80_opt.la liblin_override_binary80_opt.la
SUBDIRS += binary80
CHECK_BACKENDS += binary80
endif

if ENABLE_BINARY128
lib_LTLIBRARIES += liblin_binary128.la liblin_override_binary128.la liblin_binary128_opt.la liblin_override_binary128_opt.la
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
liblin_override_mpfr_la_CPPFLAGS = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_MPFR___
liblin_override_mpfr_opt_la_CPPFLAGS = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_MPFR___
liblin_override_mpfr_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
liblin_override_mpfr_opt_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
liblin_override_mpfr_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_override_mpfr_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_mpfr_la_CPPFLAGS = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_MPFR___ -D___MPLAPACK_INTERNAL___
liblin_mpfr_opt_la_CPPFLAGS = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_MPFR___ -D___MPLAPACK_INTERNAL___
liblin_mpfr_la_SOURCES = $(LIN_SOURCES)
liblin_mpfr_opt_la_SOURCES = $(LIN_SOURCES)
liblin_mpfr_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_mpfr_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
liblin_mpfr_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_mpfr_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_override_mpfr_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_override_mpfr_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_GMP
liblin_gmp_la_CPPFLAGS = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___ -D___MPLAPACK_INTERNAL___
liblin_gmp_opt_la_CPPFLAGS = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___ -D___MPLAPACK_INTERNAL___
liblin_gmp_la_SOURCES = $(LIN_SOURCES)
liblin_gmp_opt_la_SOURCES = $(LIN_SOURCES)
liblin_gmp_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_gmp_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_override_gmp_la_CPPFLAGS = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___
liblin_override_gmp_opt_la_CPPFLAGS = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___
liblin_override_gmp_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
liblin_override_gmp_opt_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
liblin_override_gmp_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_override_gmp_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
liblin_gmp_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_gmp_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_override_gmp_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_override_gmp_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_QD
liblin_qd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___ @qd_extra_cxxflags@ -D___MPLAPACK_INTERNAL___
liblin_qd_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___ @qd_extra_cxxflags@ -D___MPLAPACK_INTERNAL___
liblin_qd_la_SOURCES = $(LIN_SOURCES)
liblin_qd_opt_la_SOURCES = $(LIN_SOURCES)
liblin_qd_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_qd_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_override_qd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___ @qd_extra_cxxflags@
liblin_override_qd_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___ @qd_extra_cxxflags@
liblin_override_qd_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
liblin_override_qd_opt_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
liblin_override_qd_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_override_qd_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
liblin_qd_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_qd_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_override_qd_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_override_qd_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DD
liblin_dd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___ @qd_extra_cxxflags@ -D___MPLAPACK_INTERNAL___
liblin_dd_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___ @qd_extra_cxxflags@ -D___MPLAPACK_INTERNAL___
liblin_dd_la_SOURCES = $(LIN_SOURCES)
liblin_dd_opt_la_SOURCES = $(LIN_SOURCES)
liblin_dd_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_dd_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_override_dd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___ @qd_extra_cxxflags@
liblin_override_dd_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___ @qd_extra_cxxflags@
liblin_override_dd_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
liblin_override_dd_opt_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
liblin_override_dd_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_override_dd_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
liblin_dd_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_dd_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_override_dd_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_override_dd_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DOUBLE
liblin_double_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -D___MPLAPACK_BUILD_WITH_DOUBLE___ -D___MPLAPACK_INTERNAL___
liblin_double_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -D___MPLAPACK_BUILD_WITH_DOUBLE___ -D___MPLAPACK_INTERNAL___
liblin_double_la_SOURCES = $(LIN_SOURCES)
liblin_double_opt_la_SOURCES = $(LIN_SOURCES)
liblin_double_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_double_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_override_double_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DOUBLE___
liblin_override_double_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DOUBLE___
liblin_override_double_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
liblin_override_double_opt_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
liblin_override_double_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_override_double_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
liblin_double_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_double_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_override_double_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_override_double_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_BINARY80
liblin_binary80_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -D___MPLAPACK_BUILD_WITH_BINARY80___ -D___MPLAPACK_INTERNAL___
liblin_binary80_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -D___MPLAPACK_BUILD_WITH_BINARY80___ -D___MPLAPACK_INTERNAL___
liblin_binary80_la_SOURCES = $(LIN_SOURCES)
liblin_binary80_opt_la_SOURCES = $(LIN_SOURCES)
liblin_binary80_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_binary80_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_override_binary80_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -D___MPLAPACK_BUILD_WITH_BINARY80___
liblin_override_binary80_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -D___MPLAPACK_BUILD_WITH_BINARY80___
liblin_override_binary80_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
liblin_override_binary80_opt_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
liblin_override_binary80_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_override_binary80_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
liblin_binary80_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_binary80_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_override_binary80_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_override_binary80_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_BINARY128
liblin_binary128_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -D___MPLAPACK_BUILD_WITH_BINARY128___ -D___MPLAPACK_INTERNAL___
liblin_binary128_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -D___MPLAPACK_BUILD_WITH_BINARY128___ -D___MPLAPACK_INTERNAL___
liblin_binary128_la_SOURCES = $(LIN_SOURCES)
liblin_binary128_opt_la_SOURCES = $(LIN_SOURCES)
liblin_binary128_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_binary128_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_override_binary128_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -D___MPLAPACK_BUILD_WITH_BINARY128___
liblin_override_binary128_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -D___MPLAPACK_BUILD_WITH_BINARY128___
liblin_override_binary128_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
liblin_override_binary128_opt_la_SOURCES = common/Mxerbla.cpp common/iMlaenv.cpp common/Mxlaenv.cpp
liblin_override_binary128_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
liblin_override_binary128_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
liblin_binary128_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_binary128_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_override_binary128_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
liblin_override_binary128_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

LIN_SOURCES = \
common/Aladhd.cpp \
common/Alaerh.cpp \
common/Alaesm.cpp \
common/Alahd.cpp \
common/Alareq.cpp \
common/Alasum.cpp \
common/Alasvm.cpp \
common/Cchkeq.cpp \
common/Cchkgb.cpp \
common/Cchkge.cpp \
common/Cchkgt.cpp \
common/Cchkhe.cpp \
common/Cchkhe_aa.cpp \
common/Cchkhe_aa_2stage.cpp \
common/Cchkhe_rk.cpp \
common/Cchkhe_rook.cpp \
common/Cchkhp.cpp \
common/Cchklq.cpp \
common/Cchklqt.cpp \
common/Cchklqtp.cpp \
common/Cchkpb.cpp \
common/Cchkpo.cpp \
common/Cchkpp.cpp \
common/Cchkps.cpp \
common/Cchkpt.cpp \
common/Cchkq3.cpp \
common/Cchkql.cpp \
common/Cchkqp3rk.cpp \
common/Cchkqr.cpp \
common/Cchkqrt.cpp \
common/Cchkqrtp.cpp \
common/Cchkrq.cpp \
common/Cchksp.cpp \
common/Cchksy.cpp \
common/Cchksy_aa.cpp \
common/Cchksy_aa_2stage.cpp \
common/Cchksy_rk.cpp \
common/Cchksy_rook.cpp \
common/Cchktb.cpp \
common/Cchktp.cpp \
common/Cchktr.cpp \
common/Cchktsqr.cpp \
common/Cchktz.cpp \
common/Cchkunhr_col.cpp \
common/Cdrvgb.cpp \
common/Cdrvge.cpp \
common/Cdrvgt.cpp \
common/Cdrvhe.cpp \
common/Cdrvhe_aa.cpp \
common/Cdrvhe_aa_2stage.cpp \
common/Cdrvhe_rk.cpp \
common/Cdrvhe_rook.cpp \
common/Cdrvhp.cpp \
common/Cdrvls.cpp \
common/Cdrvpb.cpp \
common/Cdrvpo.cpp \
common/Cdrvpp.cpp \
common/Cdrvpt.cpp \
common/Cdrvrf1.cpp \
common/Cdrvrf2.cpp \
common/Cdrvrf3.cpp \
common/Cdrvrf4.cpp \
common/Cdrvrfp.cpp \
common/Cdrvsp.cpp \
common/Cdrvsy.cpp \
common/Cdrvsy_aa.cpp \
common/Cdrvsy_aa_2stage.cpp \
common/Cdrvsy_rk.cpp \
common/Cdrvsy_rook.cpp \
common/Cerrge.cpp \
common/Cerrgt.cpp \
common/Cerrhe.cpp \
common/Cerrlq.cpp \
common/Cerrlqt.cpp \
common/Cerrlqtp.cpp \
common/Cerrls.cpp \
common/Cerrpo.cpp \
common/Cerrps.cpp \
common/Cerrql.cpp \
common/Cerrqp.cpp \
common/Cerrqr.cpp \
common/Cerrqrt.cpp \
common/Cerrqrtp.cpp \
common/Cerrrfp.cpp \
common/Cerrrq.cpp \
common/Cerrsy.cpp \
common/Cerrtr.cpp \
common/Cerrtsqr.cpp \
common/Cerrtz.cpp \
common/Cerrunhr_col.cpp \
common/Cerrvx.cpp \
common/Cgbt01.cpp \
common/Cgbt02.cpp \
common/Cgbt05.cpp \
common/Cgennd.cpp \
common/Cgeqls.cpp \
common/Cgerqs.cpp \
common/Cget01.cpp \
common/Cget02.cpp \
common/Cget03.cpp \
common/Cget04.cpp \
common/Cget07.cpp \
common/Cget08.cpp \
common/Cgtt01.cpp \
common/Cgtt02.cpp \
common/Cgtt05.cpp \
common/Chet01.cpp \
common/Chet01_3.cpp \
common/Chet01_aa.cpp \
common/Chet01_rook.cpp \
common/Chkxer.cpp \
common/Chpt01.cpp \
common/Clahilb.cpp \
common/Claipd.cpp \
common/Claptm.cpp \
common/Clarhs.cpp \
common/Clatb4.cpp \
common/Clatb5.cpp \
common/Clatsp.cpp \
common/Clatsy.cpp \
common/Clattb.cpp \
common/Clattp.cpp \
common/Clattr.cpp \
common/Clavhe.cpp \
common/Clavhe_rook.cpp \
common/Clavhp.cpp \
common/Clavsp.cpp \
common/Clavsy.cpp \
common/Clavsy_rook.cpp \
common/Clqt01.cpp \
common/Clqt02.cpp \
common/Clqt03.cpp \
common/Clqt04.cpp \
common/Clqt05.cpp \
common/Cpbt01.cpp \
common/Cpbt02.cpp \
common/Cpbt05.cpp \
common/Cpot01.cpp \
common/Cpot02.cpp \
common/Cpot03.cpp \
common/Cpot05.cpp \
common/Cpot06.cpp \
common/Cppt01.cpp \
common/Cppt02.cpp \
common/Cppt03.cpp \
common/Cppt05.cpp \
common/Cpst01.cpp \
common/Cptt01.cpp \
common/Cptt02.cpp \
common/Cptt05.cpp \
common/Cqlt01.cpp \
common/Cqlt02.cpp \
common/Cqlt03.cpp \
common/Cqpt01.cpp \
common/Cqrt01.cpp \
common/Cqrt01p.cpp \
common/Cqrt02.cpp \
common/Cqrt03.cpp \
common/Cqrt04.cpp \
common/Cqrt05.cpp \
common/Cqrt11.cpp \
common/Cqrt12.cpp \
common/Cqrt13.cpp \
common/Cqrt14.cpp \
common/Cqrt15.cpp \
common/Cqrt16.cpp \
common/Cqrt17.cpp \
common/Crqt01.cpp \
common/Crqt02.cpp \
common/Crqt03.cpp \
common/Crzt01.cpp \
common/Crzt02.cpp \
common/Csbmv.cpp \
common/Cspt01.cpp \
common/Cspt02.cpp \
common/Cspt03.cpp \
common/Csyt01.cpp \
common/Csyt01_3.cpp \
common/Csyt01_aa.cpp \
common/Csyt01_rook.cpp \
common/Csyt02.cpp \
common/Csyt03.cpp \
common/Ctbt02.cpp \
common/Ctbt03.cpp \
common/Ctbt05.cpp \
common/Ctbt06.cpp \
common/Ctpt01.cpp \
common/Ctpt02.cpp \
common/Ctpt03.cpp \
common/Ctpt05.cpp \
common/Ctpt06.cpp \
common/Ctrt01.cpp \
common/Ctrt02.cpp \
common/Ctrt03.cpp \
common/Ctrt05.cpp \
common/Ctrt06.cpp \
common/Ctsqr01.cpp \
common/Cunhr_col01.cpp \
common/Cunhr_col02.cpp \
common/Rchkeq.cpp \
common/Rchkgb.cpp \
common/Rchkge.cpp \
common/Rchkgt.cpp \
common/Rchklq.cpp \
common/Rchklqt.cpp \
common/Rchklqtp.cpp \
common/Rchkorhr_col.cpp \
common/Rchkpb.cpp \
common/Rchkpo.cpp \
common/Rchkpp.cpp \
common/Rchkps.cpp \
common/Rchkpt.cpp \
common/Rchkq3.cpp \
common/Rchkql.cpp \
common/Rchkqp3rk.cpp \
common/Rchkqr.cpp \
common/Rchkqrt.cpp \
common/Rchkqrtp.cpp \
common/Rchkrq.cpp \
common/Rchksp.cpp \
common/Rchksy.cpp \
common/Rchksy_aa.cpp \
common/Rchksy_aa_2stage.cpp \
common/Rchksy_rk.cpp \
common/Rchksy_rook.cpp \
common/Rchktb.cpp \
common/Rchktp.cpp \
common/Rchktr.cpp \
common/Rchktsqr.cpp \
common/Rchktz.cpp \
common/Rdrvgb.cpp \
common/Rdrvge.cpp \
common/Rdrvgt.cpp \
common/Rdrvls.cpp \
common/Rdrvpb.cpp \
common/Rdrvpo.cpp \
common/Rdrvpp.cpp \
common/Rdrvpt.cpp \
common/Rdrvrf1.cpp \
common/Rdrvrf2.cpp \
common/Rdrvrf3.cpp \
common/Rdrvrf4.cpp \
common/Rdrvrfp.cpp \
common/Rdrvsp.cpp \
common/Rdrvsy.cpp \
common/Rdrvsy_aa.cpp \
common/Rdrvsy_aa_2stage.cpp \
common/Rdrvsy_rk.cpp \
common/Rdrvsy_rook.cpp \
common/Rerrge.cpp \
common/Rerrgt.cpp \
common/Rerrlq.cpp \
common/Rerrlqt.cpp \
common/Rerrlqtp.cpp \
common/Rerrls.cpp \
common/Rerrorhr_col.cpp \
common/Rerrpo.cpp \
common/Rerrps.cpp \
common/Rerrql.cpp \
common/Rerrqp.cpp \
common/Rerrqr.cpp \
common/Rerrqrt.cpp \
common/Rerrqrtp.cpp \
common/Rerrrfp.cpp \
common/Rerrrq.cpp \
common/Rerrsy.cpp \
common/Rerrtr.cpp \
common/Rerrtsqr.cpp \
common/Rerrtz.cpp \
common/Rerrvx.cpp \
common/Rgbt01.cpp \
common/Rgbt02.cpp \
common/Rgbt05.cpp \
common/Rgennd.cpp \
common/Rgeqls.cpp \
common/Rgerqs.cpp \
common/Rget01.cpp \
common/Rget02.cpp \
common/Rget03.cpp \
common/Rget04.cpp \
common/Rget06.cpp \
common/Rget07.cpp \
common/Rget08.cpp \
common/Rgtt01.cpp \
common/Rgtt02.cpp \
common/Rgtt05.cpp \
common/Rlahilb.cpp \
common/Rlaord.cpp \
common/Rlaptm.cpp \
common/Rlarhs.cpp \
common/Rlatb4.cpp \
common/Rlatb5.cpp \
common/Rlattb.cpp \
common/Rlattp.cpp \
common/Rlattr.cpp \
common/Rlavsp.cpp \
common/Rlavsy.cpp \
common/Rlavsy_rook.cpp \
common/Rlqt01.cpp \
common/Rlqt02.cpp \
common/Rlqt03.cpp \
common/Rlqt04.cpp \
common/Rlqt05.cpp \
common/Rorhr_col01.cpp \
common/Rorhr_col02.cpp \
common/Rpbt01.cpp \
common/Rpbt02.cpp \
common/Rpbt05.cpp \
common/Rpot01.cpp \
common/Rpot02.cpp \
common/Rpot03.cpp \
common/Rpot05.cpp \
common/Rpot06.cpp \
common/Rppt01.cpp \
common/Rppt02.cpp \
common/Rppt03.cpp \
common/Rppt05.cpp \
common/Rpst01.cpp \
common/Rptt01.cpp \
common/Rptt02.cpp \
common/Rptt05.cpp \
common/Rqlt01.cpp \
common/Rqlt02.cpp \
common/Rqlt03.cpp \
common/Rqpt01.cpp \
common/Rqrt01.cpp \
common/Rqrt01p.cpp \
common/Rqrt02.cpp \
common/Rqrt03.cpp \
common/Rqrt04.cpp \
common/Rqrt05.cpp \
common/Rqrt11.cpp \
common/Rqrt12.cpp \
common/Rqrt13.cpp \
common/Rqrt14.cpp \
common/Rqrt15.cpp \
common/Rqrt16.cpp \
common/Rqrt17.cpp \
common/Rrqt01.cpp \
common/Rrqt02.cpp \
common/Rrqt03.cpp \
common/Rrzt01.cpp \
common/Rrzt02.cpp \
common/Rspt01.cpp \
common/Rsyt01.cpp \
common/Rsyt01_3.cpp \
common/Rsyt01_aa.cpp \
common/Rsyt01_rook.cpp \
common/Rtbt02.cpp \
common/Rtbt03.cpp \
common/Rtbt05.cpp \
common/Rtbt06.cpp \
common/Rtplqt.cpp \
common/Rtpt01.cpp \
common/Rtpt02.cpp \
common/Rtpt03.cpp \
common/Rtpt05.cpp \
common/Rtpt06.cpp \
common/Rtrt01.cpp \
common/Rtrt02.cpp \
common/Rtrt03.cpp \
common/Rtrt05.cpp \
common/Rtrt06.cpp \
common/Rtsqr01.cpp \
common/iMlaenv.cpp \
common/icopy.cpp

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_eigdir) $(prefix)
endif
