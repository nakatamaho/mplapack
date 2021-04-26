lib_LTLIBRARIES = liblin_mpfr.la

liblin_mpfr_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___
liblin_mpfr_la_SOURCES = $(LIN_SOURCES)

LIN_SOURCES = \
../common/Rlqt01.cpp \
../common/Rlatb4.cpp \
../common/Rchklq.cpp

libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr
mplibs=-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpc -lmpfr -lgmpxx -lgmp -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___

check_PROGRAMS = $(mplapack_mpfr_lintest_PROGRAMS)
TESTS = $(check_PROGRAMS)
mplapack_mpfr_lintestdir = $(prefix)/share/mplapack/test/mplapack/lin/mpfr
mplapack_mpfr_lintest_PROGRAMS = \
Rchkaa_mpfr

Rchkaa_mpfr_SOURCES  = Rchkaa.cpp
Rchkaa_mpfr_LDFLAGS  = $(libdepends) $(mplibs)
Rchkaa_mpfr_CXXFLAGS = $(cxxflags)

