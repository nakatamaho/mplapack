check_PROGRAMS = $(mplapack_dd_test_PROGRAMS)

mplapack_dd_test_PROGRAMS = xlintstR_dd xlintstC_dd xlintstrfR_dd xlintstrfC_dd
mplapack_dd_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(check_PROGRAMS): ../liblin_dd.la ../liblin_override_dd.la

libdepends = -L$(top_builddir)/mplapack/test/lin -llin_dd -llin_override_dd -L$(top_builddir)/mplapack/test/matgen -lmatgen_dd -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(top_builddir)/mpblas/reference -lmpblas_dd
mplibs = -L$(MPC_LIBDIR) -L$(QD_LIBDIR) -lqd -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/mplapack/test/fable

xlintstR_dd_SOURCES  = ../common/Rchkaa.cpp
xlintstR_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_dd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_dd_CXXFLAGS = $(cxxflags) 

xlintstC_dd_SOURCES  = ../common/Cchkaa.cpp
xlintstC_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_dd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_dd_CXXFLAGS = $(cxxflags)

xlintstrfR_dd_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_dd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_dd_CXXFLAGS = $(cxxflags)

xlintstrfC_dd_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_dd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_dd_CXXFLAGS = $(cxxflags)

check: Rtest.out Ctest.out Rtest_rfp.out Ctest_rfp.out

clean-local:
	rm -f *.out *~ \#*

Rtest.out: ../Rtest.in xlintstR_dd
	@echo Testing Real MPLAPACK linear equation routines
	./xlintstR_dd < ../Rtest.in > $@ 2>&1

Ctest.out: ../Ctest.in xlintstC_dd
	@echo Testing Complex MPLAPACK linear equation routines
	./xlintstC_dd < ../Ctest.in > $@ 2>&1

Rtest_rfp.out: ../Rtest_rfp.in ./xlintstrfR_dd
	@echo Testing Real MPLAPACK RFP prototype linear equation routines
	./xlintstrfR_dd < ../Rtest_rfp.in > $@ 2>&1

Ctest_rfp.out: ../Ctest_rfp.in ./xlintstrfC_dd
	@echo Testing Complex MPLAPACK RFP prototype linear equation routines
	./xlintstrfC_dd < ../Ctest_rfp.in > $@ 2>&1
