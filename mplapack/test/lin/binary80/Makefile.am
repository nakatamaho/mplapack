check_PROGRAMS = $(mplapack__Float64x_test_PROGRAMS)

mplapack__Float64x_test_PROGRAMS = xlintstR_binary80 xlintstC_binary80 xlintstrfR_binary80 xlintstrfC_binary80
mplapack__Float64x_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(check_PROGRAMS): ../liblin_binary80.la ../liblin_override_binary80.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_binary80.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_binary80.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_binary80.a,$(top_builddir)/mpblas/reference/.libs/libmpblas_binary80.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_binary80.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/lin -llin_binary80 -llin_override_binary80 -L$(top_builddir)/mplapack/test/matgen -lmatgen_binary80 -L$(top_builddir)/mplapack/reference -lmplapack_binary80 -L$(top_builddir)/mpblas/reference -lmpblas_binary80
endif

mplibs = -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_BINARY80___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/fable

xlintstR__Float64x_SOURCES  = ../common/Rchkaa.cpp
xlintstR__Float64x_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR__Float64x_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR__Float64x_CXXFLAGS = $(cxxflags)

xlintstC__Float64x_SOURCES  = ../common/Cchkaa.cpp
xlintstC__Float64x_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC__Float64x_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC__Float64x_CXXFLAGS = $(cxxflags)

xlintstrfR__Float64x_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR__Float64x_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR__Float64x_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR__Float64x_CXXFLAGS = $(cxxflags)

xlintstrfC__Float64x_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC__Float64x_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC__Float64x_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC__Float64x_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapack__Float64x_testdir) $(prefix)
endif

# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <test>.out)
# -----------------------------------------------------------------------------

F64X_LIN_TESTS = Rtest Ctest Rtest_rfp Ctest_rfp

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
# -----------------------------------------------------------------------------

VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended (last assignment wins)
F64X_EXTRA_ENV ?=

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
.PHONY: f64x-check FORCE
FORCE:

# All output logs: <test>.out
F64X_LIN_OUTS = $(foreach t,$(F64X_LIN_TESTS),$(t).out)
CLEANFILES = $(F64X_LIN_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(F64X_LIN_OUTS): FORCE

check-local: f64x-check

.PRECIOUS: $(F64X_LIN_OUTS)

f64x-check: $$(F64X_LIN_OUTS)
	@echo "_Float64x lin tests done."

clean-local:
	rm -f *.out *~ \#*

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

Rtest.out: ../Rtest.in xlintstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK linear equation routines"; \
	set -e; { \
		$(F64X_EXTRA_ENV) $$wrapper $(WINE) ./xlintstR_binary80$(EXEEXT) < ../Rtest.in; \
	} > $@ 2>&1

Ctest.out: ../Ctest.in xlintstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK linear equation routines"; \
	set -e; { \
		$(F64X_EXTRA_ENV) $$wrapper $(WINE) ./xlintstC_binary80$(EXEEXT) < ../Ctest.in; \
	} > $@ 2>&1

Rtest_rfp.out: ../Rtest_rfp.in xlintstrfR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(F64X_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfR_binary80$(EXEEXT) < ../Rtest_rfp.in; \
	} > $@ 2>&1

Ctest_rfp.out: ../Ctest_rfp.in xlintstrfC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(F64X_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfC_binary80$(EXEEXT) < ../Ctest_rfp.in; \
	} > $@ 2>&1
