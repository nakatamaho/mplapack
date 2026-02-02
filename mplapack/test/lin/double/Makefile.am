check_PROGRAMS = $(mplapack_double_test_PROGRAMS)

mplapack_double_test_PROGRAMS = xlintstR_double xlintstC_double xlintstrfR_double xlintstrfC_double
mplapack_double_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(check_PROGRAMS): ../liblin_double.la ../liblin_override_double.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_double.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_double.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_double.a,$(top_builddir)/mpblas/reference/.libs/libmpblas_double.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_double.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/lin -llin_double -llin_override_double -L$(top_builddir)/mplapack/test/matgen -lmatgen_double -L$(top_builddir)/mplapack/reference -lmplapack_double -L$(top_builddir)/mpblas/reference -lmpblas_double
endif

mplibs = -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_DOUBLE___ -I$(top_srcdir)/fable -D___MPLAPACK_INTERNAL___

xlintstR_double_SOURCES  = ../common/Rchkaa.cpp
xlintstR_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_double_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_double_CXXFLAGS = $(cxxflags) 

xlintstC_double_SOURCES  = ../common/Cchkaa.cpp
xlintstC_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_double_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_double_CXXFLAGS = $(cxxflags)

xlintstrfR_double_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_double_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_double_CXXFLAGS = $(cxxflags)

xlintstrfC_double_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_double_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_double_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapack_double_testdir) $(prefix)
endif

# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <test>.out)
# -----------------------------------------------------------------------------

DOUBLE_LIN_TESTS = Rtest Ctest Rtest_rfp Ctest_rfp

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
# -----------------------------------------------------------------------------

VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended (last assignment wins)
DOUBLE_EXTRA_ENV ?=

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
.PHONY: double-check FORCE
FORCE:

# All output logs: <test>.out
DOUBLE_LIN_OUTS = $(foreach t,$(DOUBLE_LIN_TESTS),$(t).out)
CLEANFILES = $(DOUBLE_LIN_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(DOUBLE_LIN_OUTS): FORCE

check-local: double-check

.PRECIOUS: $(DOUBLE_LIN_OUTS)

double-check: $$(DOUBLE_LIN_OUTS)
	@echo "double lin tests done."

clean-local:
	rm -f *.out *~ \#*

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

Rtest.out: ../Rtest.in xlintstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK linear equation routines"; \
	set -e; { \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xlintstR_double$(EXEEXT) < ../Rtest.in; \
	} > $@ 2>&1

Ctest.out: ../Ctest.in xlintstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK linear equation routines"; \
	set -e; { \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xlintstC_double$(EXEEXT) < ../Ctest.in; \
	} > $@ 2>&1

Rtest_rfp.out: ../Rtest_rfp.in xlintstrfR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfR_double$(EXEEXT) < ../Rtest_rfp.in; \
	} > $@ 2>&1

Ctest_rfp.out: ../Ctest_rfp.in xlintstrfC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfC_double$(EXEEXT) < ../Ctest_rfp.in; \
	} > $@ 2>&1
