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
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/fable

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

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapack_qd_testdir) $(prefix)
endif

# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <test>.out)
# -----------------------------------------------------------------------------

QD_LIN_TESTS = Rtest Ctest Rtest_rfp Ctest_rfp

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
# -----------------------------------------------------------------------------

VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended (last assignment wins)
QD_EXTRA_ENV ?=

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
.PHONY: qd-check FORCE
FORCE:

# All output logs: <test>.out
QD_LIN_OUTS = $(foreach t,$(QD_LIN_TESTS),$(t).out)
CLEANFILES = $(QD_LIN_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(QD_LIN_OUTS): FORCE

check-local: qd-check

.PRECIOUS: $(QD_LIN_OUTS)

qd-check: $$(QD_LIN_OUTS)
	@echo "QD lin tests done."

clean-local:
	rm -f *.out *~ \#*

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

Rtest.out: ../Rtest.in xlintstR_qd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK linear equation routines"; \
	set -e; { \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstR_qd$(EXEEXT) < ../Rtest.in; \
	} > $@ 2>&1

Ctest.out: ../Ctest.in xlintstC_qd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK linear equation routines"; \
	set -e; { \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstC_qd$(EXEEXT) < ../Ctest.in; \
	} > $@ 2>&1

Rtest_rfp.out: ../Rtest_rfp.in xlintstrfR_qd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfR_qd$(EXEEXT) < ../Rtest_rfp.in; \
	} > $@ 2>&1

Ctest_rfp.out: ../Ctest_rfp.in xlintstrfC_qd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfC_qd$(EXEEXT) < ../Ctest_rfp.in; \
	} > $@ 2>&1
