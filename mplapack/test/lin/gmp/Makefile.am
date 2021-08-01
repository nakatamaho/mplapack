bin_PROGRAMS = xlintstR_gmp \
xlintstC_gmp \
xlintstrfR_gmp \
xlintstrfC_gmp

mplapack_gmp_lintestdir = $(prefix)/share/mplapack/test/mplapack/lin/gmp

$(bin_PROGRAMS): ../liblin_gmp.la ../liblin_override.la

libdepends = -L$(top_builddir)/mplapack/test/lin -llin_gmp -llin_override -L$(top_builddir)/mplapack/test/matgen -lmatgen_gmp -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp
mplibs = -L$(MPC_LIBDIR) -L$(GMP_LIBDIR) -L$(GMP_LIBDIR) -lmpc -lgmp -lgmpxx -lgmp -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(GMP_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/gmpc++ -D___MPLAPACK_BUILD_WITH_GMP___ -I$(top_srcdir)/mplapack/test/fable

xlintstR_gmp_SOURCES  = ../common/Rchkaa.cpp ../../../../mpblas/test/common/mplapack.test.cpp 
xlintstR_gmp_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_gmp_CXXFLAGS = $(cxxflags) 

xlintstC_gmp_SOURCES  = ../common/Cchkaa.cpp ../../../../mpblas/test/common/mplapack.test.cpp
xlintstC_gmp_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_gmp_CXXFLAGS = $(cxxflags)

xlintstrfR_gmp_SOURCES  = ../common/Rchkrfp.cpp ../../../../mpblas/test/common/mplapack.test.cpp
xlintstrfR_gmp_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_gmp_CXXFLAGS = $(cxxflags)

xlintstrfC_gmp_SOURCES  = ../common/Cchkrfp.cpp ../../../../mpblas/test/common/mplapack.test.cpp
xlintstrfC_gmp_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_gmp_CXXFLAGS = $(cxxflags)

test: Rtest.out Ctest.out Rtest_rfp.out Ctest_rfp.out

clean-local:
	rm -f *.out *~ \#*

Rtest.out: Rtest.in xlintstR_gmp
	@echo Testing Real MPLAPACK linear equation routines
	./xlintstR_gmp < Rtest.in > $@ 2>&1

Ctest.out: Ctest.in xlintstC_gmp
	@echo Testing Complex MPLAPACK linear equation routines
	./xlintstC_gmp < Ctest.in > $@ 2>&1

Rtest_rfp.out: Rtest_rfp.in 
	@echo Testing Real MPLAPACK RFP prototype linear equation routines
	./xlintstrfR_gmp < Rtest_rfp.in > $@ 2>&1

Ctest_rfp.out: Ctest_rfp.in
	@echo Testing Complex MPLAPACK RFP prototype linear equation routines
	./xlintstrfC_gmp < Ctest_rfp.in > $@ 2>&1