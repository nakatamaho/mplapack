check_PROGRAMS = $(mplapack_gmp_test_PROGRAMS)

mplapack_gmp_test_PROGRAMS = xlintstR_gmp xlintstC_gmp xlintstrfR_gmp xlintstrfC_gmp
mplapack_gmp_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(check_PROGRAMS): ../liblin_gmp.la ../liblin_override_gmp.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_gmp.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_gmp.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_gmp.a,$(top_builddir)/mpblas/reference/.libs/libmpblas_gmp.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_gmp.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/lin -llin_gmp -llin_override_gmp -L$(top_builddir)/mplapack/test/matgen -lmatgen_gmp -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp
endif

mplibs = -L$(GMP_LIBDIR) -lgmpxx -lgmp -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/mplapack/test/fable

xlintstR_gmp_SOURCES  = ../common/Rchkaa.cpp
xlintstR_gmp_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_gmp_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_gmp_CXXFLAGS = $(cxxflags) 

xlintstC_gmp_SOURCES  = ../common/Cchkaa.cpp
xlintstC_gmp_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_gmp_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_gmp_CXXFLAGS = $(cxxflags)

xlintstrfR_gmp_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_gmp_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_gmp_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_gmp_CXXFLAGS = $(cxxflags)

xlintstrfC_gmp_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_gmp_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_gmp_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_gmp_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapack_gmp_testdir) $(prefix)
endif

check: Rtest.out Ctest.out Rtest_rfp.out Ctest_rfp.out

clean-local:
	rm -f *.out *~ \#*

Rtest.out: ../Rtest.in xlintstR_gmp
	@echo Testing Real MPLAPACK linear equation routines
	./xlintstR_gmp < ../Rtest.in > $@ 2>&1

Ctest.out: ../Ctest.in xlintstC_gmp
	@echo Testing Complex MPLAPACK linear equation routines
	./xlintstC_gmp < ../Ctest.in > $@ 2>&1

Rtest_rfp.out: ../Rtest_rfp.in xlintstrfR_gmp
	@echo Testing Real MPLAPACK RFP prototype linear equation routines
	./xlintstrfR_gmp < ../Rtest_rfp.in > $@ 2>&1

Ctest_rfp.out: ../Ctest_rfp.in xlintstrfC_gmp
	@echo Testing Complex MPLAPACK RFP prototype linear equation routines
	./xlintstrfC_gmp < ../Ctest_rfp.in > $@ 2>&1
