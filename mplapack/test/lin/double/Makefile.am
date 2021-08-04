bin_PROGRAMS = xlintstR_double \
xlintstC_double \
xlintstrfR_double \
xlintstrfC_double

mplapack_double_lintestdir = $(prefix)/share/mplapack/test/mplapack/lin/double

$(bin_PROGRAMS): ../liblin_double.la ../liblin_override_double.la

libdepends = -L.. -llin_double -llin_override_double -L$(top_builddir)/mplapack/test/matgen -lmatgen_double -L$(top_builddir)/mplapack/reference -lmplapack_double -L$(top_builddir)/mpblas/reference -lmpblas_double
mplibs = -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_DOUBLE___ -I$(top_srcdir)/mplapack/test/fable

xlintstR_double_SOURCES  = ../common/Rchkaa.cpp
xlintstR_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_double_CXXFLAGS = $(cxxflags) 

xlintstC_double_SOURCES  = ../common/Cchkaa.cpp
xlintstC_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_double_CXXFLAGS = $(cxxflags)

xlintstrfR_double_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_double_CXXFLAGS = $(cxxflags)

xlintstrfC_double_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_double_CXXFLAGS = $(cxxflags)

do-tests: Rtest.out Ctest.out Rtest_rfp.out Ctest_rfp.out

clean-local:
	rm -f *.out *~ \#*

Rtest.out: Rtest.in xlintstR_double
	@echo Testing Real MPLAPACK linear equation routines
	./xlintstR_double < Rtest.in > $@ 2>&1

Ctest.out: Ctest.in xlintstC_double
	@echo Testing Complex MPLAPACK linear equation routines
	./xlintstC_double < Ctest.in > $@ 2>&1

Rtest_rfp.out: Rtest_rfp.in 
	@echo Testing Real MPLAPACK RFP prototype linear equation routines
	./xlintstrfR_double < Rtest_rfp.in > $@ 2>&1

Ctest_rfp.out: Ctest_rfp.in
	@echo Testing Complex MPLAPACK RFP prototype linear equation routines
	./xlintstrfC_double < Ctest_rfp.in > $@ 2>&1
