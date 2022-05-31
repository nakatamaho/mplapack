check_PROGRAMS = $(mplapack__Float64x_test_PROGRAMS)

mplapack__Float64x_test_PROGRAMS = xlintstR__Float64x xlintstC__Float64x xlintstrfR__Float64x xlintstrfC__Float64x
mplapack__Float64x_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(check_PROGRAMS): ../liblin__Float64x.la ../liblin_override__Float64x.la

libdepends = -L.. -llin__Float64x -llin_override__Float64x -L$(top_builddir)/mplapack/test/matgen -lmatgen__Float64x -L$(top_builddir)/mplapack/reference -lmplapack__Float64x -L$(top_builddir)/mpblas/reference -lmpblas__Float64x
mplibs = -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT64X___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/mplapack/test/fable

xlintstR__Float64x_SOURCES  = ../common/Rchkaa.cpp
xlintstR__Float64x_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR__Float64x_CXXFLAGS = $(cxxflags) 

xlintstC__Float64x_SOURCES  = ../common/Cchkaa.cpp
xlintstC__Float64x_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC__Float64x_CXXFLAGS = $(cxxflags)

xlintstrfR__Float64x_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR__Float64x_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR__Float64x_CXXFLAGS = $(cxxflags)

xlintstrfC__Float64x_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC__Float64x_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC__Float64x_CXXFLAGS = $(cxxflags)

check: Rtest.out Ctest.out Rtest_rfp.out Ctest_rfp.out

clean-local:
	rm -f *.out *~ \#*

Rtest.out: ../Rtest.in xlintstR__Float64x
	@echo Testing Real MPLAPACK linear equation routines
	./xlintstR__Float64x < ../Rtest.in > $@ 2>&1

Ctest.out: ../Ctest.in xlintstC__Float64x
	@echo Testing Complex MPLAPACK linear equation routines
	./xlintstC__Float64x < ../Ctest.in > $@ 2>&1

Rtest_rfp.out: ../Rtest_rfp.in xlintstrfR__Float64x
	@echo Testing Real MPLAPACK RFP prototype linear equation routines
	./xlintstrfR__Float64x < ../Rtest_rfp.in > $@ 2>&1

Ctest_rfp.out: ../Ctest_rfp.in xlintstrfC__Float64x
	@echo Testing Complex MPLAPACK RFP prototype linear equation routines
	./xlintstrfC__Float64x < ../Ctest_rfp.in > $@ 2>&1
