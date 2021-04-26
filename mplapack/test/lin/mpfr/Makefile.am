lib_LTLIBRARIES = liblin_mpfr.la

liblin_mpfr_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___
liblin_mpfr_la_SOURCES = $(LIN_SOURCES)
if IS_MACOSX
liblin_mpfr_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif

LIN_SOURCES = \
../common/Rlqt01.cpp
