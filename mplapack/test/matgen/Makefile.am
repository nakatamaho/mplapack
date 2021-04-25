lib_LTLIBRARIES = libmatgen_mpfr.la

libmatgen_mpfr_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___
libmatgen_mpfr_la_SOURCES = $(MATGEN_SOURCES)
if IS_MACOSX
libmatgen_mpfr_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif

MATGEN_SOURCES = \
Clagge.cpp \
Claghe.cpp \
Clagsy.cpp \
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
