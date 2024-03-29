check_PROGRAMS = $(mplapack_qd_test_PROGRAMS)

mplapack_qd_test_PROGRAMS = xlintstR_qd xlintstC_qd xlintstrfR_qd xlintstrfC_qd
mplapack_qd_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(check_PROGRAMS): ../liblin_qd.la ../liblin_override_qd.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_qd.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_qd.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_qd.a,$(top_builddir)/mpblas/reference/.libs/libmpblas_qd.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_qd.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/lin -llin_qd -llin_override_qd -L$(top_builddir)/mplapack/test/matgen -lmatgen_qd -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(top_builddir)/mpblas/reference -lmpblas_qd
endif

mplibs = -L$(QD_LIBDIR) -lqd -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/mplapack/test/fable

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapack_qd_testdir) $(prefix)
endif

xlintstR_qd_SOURCES  = ../common/Rchkaa.cpp
xlintstR_qd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_qd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_qd_CXXFLAGS = $(cxxflags) 

xlintstC_qd_SOURCES  = ../common/Cchkaa.cpp
xlintstC_qd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_qd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_qd_CXXFLAGS = $(cxxflags)

xlintstrfR_qd_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_qd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_qd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_qd_CXXFLAGS = $(cxxflags)

xlintstrfC_qd_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_qd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_qd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_qd_CXXFLAGS = $(cxxflags)

check: Rtest.out Ctest.out Rtest_rfp.out Ctest_rfp.out

clean-local:
	rm -f *.out *~ \#*

Rtest.out: ../Rtest.in xlintstR_qd
	@echo Testing Real MPLAPACK linear equation routines
	./xlintstR_qd < ../Rtest.in > $@ 2>&1

Ctest.out: ../Ctest.in xlintstC_qd
	@echo Testing Complex MPLAPACK linear equation routines
	./xlintstC_qd < ../Ctest.in > $@ 2>&1

Rtest_rfp.out: ../Rtest_rfp.in 
	@echo Testing Real MPLAPACK RFP prototype linear equation routines
	./xlintstrfR_qd < ../Rtest_rfp.in > $@ 2>&1

Ctest_rfp.out: ../Ctest_rfp.in
	@echo Testing Complex MPLAPACK RFP prototype linear equation routines
	./xlintstrfC_qd < ../Ctest_rfp.in > $@ 2>&1
