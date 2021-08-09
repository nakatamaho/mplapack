bin_PROGRAMS = xlintstR_mpfr \
xlintstC_mpfr \
xlintstrfR_mpfr \
xlintstrfC_mpfr

mplapack_mpfr_lintestdir = $(prefix)/share/mplapack/test/mplapack/lin/mpfr

$(bin_PROGRAMS): ../liblin_mpfr.la ../liblin_override_mpfr.la

libdepends = -L$(top_builddir)/mplapack/test/lin -llin_mpfr -llin_override_mpfr -L$(top_builddir)/mplapack/test/matgen -lmatgen_mpfr -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr
mplibs = -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpc -lmpfr -lgmpxx -lgmp -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/mplapack/test/fable

xlintstR_mpfr_SOURCES  = ../common/Rchkaa.cpp
xlintstR_mpfr_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_mpfr_CXXFLAGS = $(cxxflags) 

xlintstC_mpfr_SOURCES  = ../common/Cchkaa.cpp
xlintstC_mpfr_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_mpfr_CXXFLAGS = $(cxxflags)

xlintstrfR_mpfr_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_mpfr_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_mpfr_CXXFLAGS = $(cxxflags)

xlintstrfC_mpfr_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_mpfr_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_mpfr_CXXFLAGS = $(cxxflags)

test: Rtest.out Ctest.out Rtest_rfp.out Ctest_rfp.out

clean-local:
	rm -f *.out *~ \#*

Rtest.out: ../Rtest.in xlintstR_mpfr
	@echo Testing Real MPLAPACK linear equation routines
	./xlintstR_mpfr < ../Rtest.in > $@ 2>&1

Ctest.out: ../Ctest.in xlintstC_mpfr
	@echo Testing Complex MPLAPACK linear equation routines
	./xlintstC_mpfr < ../Ctest.in > $@ 2>&1

Rtest_rfp.out: ../Rtest_rfp.in xlintstrfR_mpfr
	@echo Testing Real MPLAPACK RFP prototype linear equation routines
	./xlintstrfR_mpfr < ../Rtest_rfp.in > $@ 2>&1

Ctest_rfp.out: ../Ctest_rfp.in xlintstrfC_mpfr
	@echo Testing Complex MPLAPACK RFP prototype linear equation routines
	./xlintstrfC_mpfr < ../Ctest_rfp.in > $@ 2>&1
