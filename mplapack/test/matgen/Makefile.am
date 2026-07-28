lib_LTLIBRARIES =

if ENABLE_MPFR
lib_LTLIBRARIES += libmatgen_mpfr.la libmatgen_mpfr_opt.la
endif

if ENABLE_GMP
lib_LTLIBRARIES += libmatgen_gmp.la libmatgen_gmp_opt.la
endif

if ENABLE_DD
lib_LTLIBRARIES += libmatgen_dd.la libmatgen_dd_opt.la
endif

if ENABLE_QD
lib_LTLIBRARIES += libmatgen_qd.la libmatgen_qd_opt.la
endif

if ENABLE_DOUBLE
lib_LTLIBRARIES += libmatgen_double.la libmatgen_double_opt.la
endif

if ENABLE_BINARY80
lib_LTLIBRARIES += libmatgen_binary80.la libmatgen_binary80_opt.la
endif

if ENABLE_BINARY128
lib_LTLIBRARIES += libmatgen_binary128.la libmatgen_binary128_opt.la
endif


if ENABLE_MPFR
libmatgen_mpfr_la_CPPFLAGS = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_MPFR -DMPLAPACK_INTERNAL
libmatgen_mpfr_opt_la_CPPFLAGS = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_MPFR -DMPLAPACK_INTERNAL
libmatgen_mpfr_la_SOURCES = $(MATGEN_SOURCES)
libmatgen_mpfr_opt_la_SOURCES = $(MATGEN_SOURCES)
libmatgen_mpfr_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
libmatgen_mpfr_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
libmatgen_mpfr_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
libmatgen_mpfr_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_GMP
libmatgen_gmp_la_CPPFLAGS = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(GMP_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_GMP -DMPLAPACK_INTERNAL
libmatgen_gmp_opt_la_CPPFLAGS = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(GMP_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_GMP -DMPLAPACK_INTERNAL
libmatgen_gmp_la_SOURCES = $(MATGEN_SOURCES)
libmatgen_gmp_opt_la_SOURCES = $(MATGEN_SOURCES)
libmatgen_gmp_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
libmatgen_gmp_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
libmatgen_gmp_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
libmatgen_gmp_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_QD
libmatgen_qd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_QD @qd_extra_cxxflags@ -DMPLAPACK_INTERNAL
libmatgen_qd_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_QD @qd_extra_cxxflags@ -DMPLAPACK_INTERNAL
libmatgen_qd_la_SOURCES = $(MATGEN_SOURCES)
libmatgen_qd_opt_la_SOURCES = $(MATGEN_SOURCES)
libmatgen_qd_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
libmatgen_qd_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
libmatgen_qd_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
libmatgen_qd_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DD
libmatgen_dd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_DD @qd_extra_cxxflags@ -DMPLAPACK_INTERNAL
libmatgen_dd_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(QD_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_DD @qd_extra_cxxflags@ -DMPLAPACK_INTERNAL
libmatgen_dd_la_SOURCES = $(MATGEN_SOURCES)
libmatgen_dd_opt_la_SOURCES = $(MATGEN_SOURCES)
libmatgen_dd_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
libmatgen_dd_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
libmatgen_dd_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
libmatgen_dd_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DOUBLE
libmatgen_double_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -DMPLAPACK_BUILD_WITH_DOUBLE -DMPLAPACK_INTERNAL
libmatgen_double_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -DMPLAPACK_BUILD_WITH_DOUBLE -DMPLAPACK_INTERNAL
libmatgen_double_la_SOURCES = $(MATGEN_SOURCES)
libmatgen_double_opt_la_SOURCES = $(MATGEN_SOURCES)
libmatgen_double_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
libmatgen_double_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
libmatgen_double_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
libmatgen_double_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_BINARY80
libmatgen_binary80_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -DMPLAPACK_BUILD_WITH_BINARY80 -DMPLAPACK_INTERNAL
libmatgen_binary80_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -DMPLAPACK_BUILD_WITH_BINARY80 -DMPLAPACK_INTERNAL
libmatgen_binary80_la_SOURCES = $(MATGEN_SOURCES)
libmatgen_binary80_opt_la_SOURCES = $(MATGEN_SOURCES)
libmatgen_binary80_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
libmatgen_binary80_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
libmatgen_binary80_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
libmatgen_binary80_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_BINARY128
libmatgen_binary128_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -DMPLAPACK_BUILD_WITH_BINARY128 -DMPLAPACK_INTERNAL
libmatgen_binary128_opt_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -DMPLAPACK_BUILD_WITH_BINARY128 -DMPLAPACK_INTERNAL
libmatgen_binary128_la_SOURCES = $(MATGEN_SOURCES)
libmatgen_binary128_opt_la_SOURCES = $(MATGEN_SOURCES)
libmatgen_binary128_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
libmatgen_binary128_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO)
if IS_MACOS
libmatgen_binary128_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
libmatgen_binary128_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

MATGEN_SOURCES = \
Clagge.cpp \
Claghe.cpp \
Clagsy.cpp \
Clahilb.cpp \
Clakf2.cpp \
Clarge.cpp \
Clarnd.cpp \
Claror.cpp \
Clarot.cpp \
Clatm1.cpp \
Clatm2.cpp \
Clatm3.cpp \
Clatm5.cpp \
Clatm6.cpp \
Clatme.cpp \
Clatmr.cpp \
Clatms.cpp \
Clatmt.cpp \
Rlagge.cpp \
Rlagsy.cpp \
Rlahilb.cpp \
Rlakf2.cpp \
Rlaran.cpp \
Rlarge.cpp \
Rlarnd.cpp \
Rlaror.cpp \
Rlarot.cpp \
Rlatm1.cpp \
Rlatm2.cpp \
Rlatm3.cpp \
Rlatm5.cpp \
Rlatm6.cpp \
Rlatm7.cpp \
Rlatme.cpp \
Rlatmr.cpp \
Rlatms.cpp \
Rlatmt.cpp
