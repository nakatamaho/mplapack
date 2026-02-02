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
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/fable

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

# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <test>.out)
# -----------------------------------------------------------------------------

GMP_LIN_TESTS = Rtest Ctest Rtest_rfp Ctest_rfp

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
# -----------------------------------------------------------------------------

VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended (last assignment wins)
GMP_EXTRA_ENV ?=

# -----------------------------------------------------------------------------
# Cross-compilation support (wine for MinGW)
# -----------------------------------------------------------------------------
if IS_MINGW
WINE = wine
else
WINE =
endif

# -----------------------------------------------------------------------------
# Hook into Automake's "make check" (parallel-capable)
# -----------------------------------------------------------------------------

# GNU make: allow deferred expansion of prerequisites (Automake may reorder vars/rules)
.SECONDEXPANSION:

# Force test execution even if *.out already exists
.PHONY: gmp-check FORCE
FORCE:

# All output logs: <test>.out
GMP_LIN_OUTS = $(foreach t,$(GMP_LIN_TESTS),$(t).out)
CLEANFILES = $(GMP_LIN_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(GMP_LIN_OUTS): FORCE

check-local: gmp-check

.PRECIOUS: $(GMP_LIN_OUTS)

gmp-check: $$(GMP_LIN_OUTS)
	@echo "GMP lin tests done."

clean-local:
	rm -f *.out *~ \#*

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

Rtest.out: ../Rtest.in xlintstR_gmp$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK linear equation routines"; \
	set -e; { \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xlintstR_gmp$(EXEEXT) < ../Rtest.in; \
	} > $@ 2>&1

Ctest.out: ../Ctest.in xlintstC_gmp$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK linear equation routines"; \
	set -e; { \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xlintstC_gmp$(EXEEXT) < ../Ctest.in; \
	} > $@ 2>&1

Rtest_rfp.out: ../Rtest_rfp.in xlintstrfR_gmp$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfR_gmp$(EXEEXT) < ../Rtest_rfp.in; \
	} > $@ 2>&1

Ctest_rfp.out: ../Ctest_rfp.in xlintstrfC_gmp$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfC_gmp$(EXEEXT) < ../Ctest_rfp.in; \
	} > $@ 2>&1

