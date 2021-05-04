lib_LTLIBRARIES = libeig_mpfr.la

libeig_mpfr_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -I$(top_srcdir)/mplapack/test/fable
libeig_mpfr_la_SOURCES = $(EIG_SOURCES)

EIG_SOURCES = \
../common/Rerrge.cpp \
../common/Rchkge.cpp \
../common/Rdrvge.cpp \
../common/Rchkgb.cpp \
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
../common/Rgbt01.cpp \
../common/Rgbt02.cpp \
../common/Rgbt05.cpp

libdepends = -L. -leig_mpfr -L$(top_builddir)/mplapack/test/matgen -lmatgen_mpfr -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr 
mplibs=-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpc -lmpfr -lgmpxx -lgmp -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -I$(top_srcdir)/mplapack/test/fable

check_PROGRAMS = $(mplapack_mpfr_eigtest_PROGRAMS)
TESTS = $(check_PROGRAMS)
mplapack_mpfr_eigtestdir = $(prefix)/share/mplapack/test/mplapack/eig/mpfr
mplapack_mpfr_eigtest_PROGRAMS = \
xeigtstR_mpfr

xeigtstR_mpfr_SOURCES  = ../common/Rchkaa.cpp $(top_builddir)/mpblas/test/common/mplapack.test.cpp
xeigtstR_mpfr_LDFLAGS  = $(libdepends) $(mplibs)
xeigtstR_mpfr_CXXFLAGS = $(cxxflags)
