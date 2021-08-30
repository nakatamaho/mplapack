lib_LTLIBRARIES =

if ENABLE_MPFR
lib_LTLIBRARIES += libmatgen_mpfr.la
endif

if ENABLE_GMP
lib_LTLIBRARIES += libmatgen_gmp.la
endif

if ENABLE_DD
lib_LTLIBRARIES += libmatgen_dd.la
endif

if ENABLE_QD
lib_LTLIBRARIES += libmatgen_qd.la
endif

if ENABLE_DOUBLE
lib_LTLIBRARIES += libmatgen_double.la
endif

if ENABLE__FLOAT64X
lib_LTLIBRARIES += libmatgen__Float64x.la
endif

if ENABLE__FLOAT128
lib_LTLIBRARIES += libmatgen__Float128.la
endif


if ENABLE_MPFR
libmatgen_mpfr_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -D___MPLAPACK_INTERNAL___
libmatgen_mpfr_la_SOURCES = $(MATGEN_SOURCES)
if IS_MACOSX
libmatgen_mpfr_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_GMP
libmatgen_gmp_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___ -D___MPLAPACK_INTERNAL___
libmatgen_gmp_la_SOURCES = $(MATGEN_SOURCES)
if IS_MACOSX
libmatgen_gmp_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_QD
libmatgen_qd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___ -D___MPLAPACK_INTERNAL___
libmatgen_qd_la_SOURCES = $(MATGEN_SOURCES)
if IS_MACOSX
libmatgen_qd_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DD
libmatgen_dd_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___ -D___MPLAPACK_INTERNAL___
libmatgen_dd_la_SOURCES = $(MATGEN_SOURCES)
if IS_MACOSX
libmatgen_dd_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DOUBLE
libmatgen_double_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_DOUBLE___ -D___MPLAPACK_INTERNAL___
libmatgen_double_la_SOURCES = $(MATGEN_SOURCES)
if IS_MACOSX
libmatgen_double_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE__FLOAT64X
libmatgen__Float64x_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT64X___ -D___MPLAPACK_INTERNAL___
libmatgen__Float64x_la_SOURCES = $(MATGEN_SOURCES)
if IS_MACOSX
libmatgen__Float64x_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE__FLOAT128
libmatgen__Float128_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT128___ -D___MPLAPACK_INTERNAL___
libmatgen__Float128_la_SOURCES = $(MATGEN_SOURCES)
if IS_MACOSX
libmatgen__Float128_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
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
