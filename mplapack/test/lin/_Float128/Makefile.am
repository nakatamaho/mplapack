check_PROGRAMS = $(mplapack__Float128_test_PROGRAMS)

mplapack__Float128_test_PROGRAMS = xlintstR__Float128 xlintstC__Float128 xlintstrfR__Float128 xlintstrfC__Float128
mplapack__Float128_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(check_PROGRAMS): ../liblin__Float128.la ../liblin_override__Float128.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin__Float128.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override__Float128.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen__Float128.a,$(top_builddir)/mpblas/reference/.libs/libmpblas__Float128.a,$(top_builddir)/mplapack/reference/.libs/libmplapack__Float128.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/lin -llin__Float128 -llin_override__Float128 -L$(top_builddir)/mplapack/test/matgen -lmatgen__Float128 -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128
endif

if WANT_QUADMATH
mplibs = -lm -lquadmath
else
mplibs = -lm
endif

cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT128___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/fable

xlintstR__Float128_SOURCES  = ../common/Rchkaa.cpp
xlintstR__Float128_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR__Float128_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR__Float128_CXXFLAGS = $(cxxflags)

xlintstC__Float128_SOURCES  = ../common/Cchkaa.cpp
xlintstC__Float128_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC__Float128_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC__Float128_CXXFLAGS = $(cxxflags)

xlintstrfR__Float128_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR__Float128_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR__Float128_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR__Float128_CXXFLAGS = $(cxxflags)

xlintstrfC__Float128_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC__Float128_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC__Float128_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC__Float128_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapack__Float128_testdir) $(prefix)
endif

# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <test>.out)
# -----------------------------------------------------------------------------

F128_LIN_TESTS = Rtest Ctest Rtest_rfp Ctest_rfp

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
# -----------------------------------------------------------------------------

VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended (last assignment wins)
F128_EXTRA_ENV ?=

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
.PHONY: f128-check FORCE
FORCE:

# All output logs: <test>.out
F128_LIN_OUTS = $(foreach t,$(F128_LIN_TESTS),$(t).out)
CLEANFILES = $(F128_LIN_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(F128_LIN_OUTS): FORCE

check-local: f128-check

.PRECIOUS: $(F128_LIN_OUTS)

f128-check: $$(F128_LIN_OUTS)
	@echo "_Float128 lin tests done."

clean-local:
	rm -f *.out *~ \#*

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

Rtest.out: ../Rtest.in xlintstR__Float128$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK linear equation routines"; \
	set -e; { \
		$(F128_EXTRA_ENV) $$wrapper $(WINE) ./xlintstR__Float128$(EXEEXT) < ../Rtest.in; \
	} > $@ 2>&1

Ctest.out: ../Ctest.in xlintstC__Float128$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK linear equation routines"; \
	set -e; { \
		$(F128_EXTRA_ENV) $$wrapper $(WINE) ./xlintstC__Float128$(EXEEXT) < ../Ctest.in; \
	} > $@ 2>&1

Rtest_rfp.out: ../Rtest_rfp.in xlintstrfR__Float128$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(F128_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfR__Float128$(EXEEXT) < ../Rtest_rfp.in; \
	} > $@ 2>&1

Ctest_rfp.out: ../Ctest_rfp.in xlintstrfC__Float128$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(F128_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfC__Float128$(EXEEXT) < ../Ctest_rfp.in; \
	} > $@ 2>&1
