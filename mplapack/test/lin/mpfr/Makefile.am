lib_LTLIBRARIES = liblin_mpfr.la

liblin_mpfr_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -I$(top_srcdir)/mplapack/test/fable
liblin_mpfr_la_SOURCES = $(LIN_SOURCES)

LIN_SOURCES = \
../common/Rerrge.cpp \
../common/Rchkge.cpp \
../common/chkxer.cpp \
../common/Alasum.cpp \
../common/Alareq.cpp \
../common/Aladhd.cpp \
../common/Alahd.cpp \
../common/Alaesm.cpp \
../common/Alaerh.cpp \
../common/Rget01.cpp \
../common/Rget02.cpp \
../common/Rget03.cpp \
../common/Rget04.cpp \
../common/Rget06.cpp \
../common/Rget07.cpp \
../common/Rlatb4.cpp \
../common/Rlarhs.cpp \
../common/mplapack.test.cpp

libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L. -llin_mpfr -L$(top_builddir)/mplapack/test/matgen -lmatgen_mpfr 
mplibs=-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpc -lmpfr -lgmpxx -lgmp -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -I$(top_srcdir)/mplapack/test/fable

check_PROGRAMS = $(mplapack_mpfr_lintest_PROGRAMS)
TESTS = $(check_PROGRAMS)
mplapack_mpfr_lintestdir = $(prefix)/share/mplapack/test/mplapack/lin/mpfr
mplapack_mpfr_lintest_PROGRAMS = \
xlintstR_mpfr

xlintstR_mpfr_SOURCES  = ../common/Rchkaa.cpp ../common/mplapack.test.cpp
xlintstR_mpfr_LDFLAGS  = $(libdepends) $(mplibs)
xlintstR_mpfr_CXXFLAGS = $(cxxflags)
