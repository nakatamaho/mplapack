check_PROGRAMS = $(mplapack_double_test_PROGRAMS) $(mplapack_double_opt_test_PROGRAMS)

mplapack_double_test_PROGRAMS = xlintstR_double xlintstC_double xlintstrfR_double xlintstrfC_double
mplapack_double_opt_test_PROGRAMS = xlintstR_double_opt xlintstC_double_opt xlintstrfR_double_opt xlintstrfC_double_opt
mplapack_double_testdir = $(prefix)/lib/$(target)/mplapack/test/lin
mplapack_double_opt_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(mplapack_double_test_PROGRAMS): ../liblin_double.la ../liblin_override_double.la
$(mplapack_double_opt_test_PROGRAMS): ../liblin_double_opt.la ../liblin_override_double_opt.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_double.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_double.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_double.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_double.a,--no-whole-archive
opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_double_opt.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_double_opt.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_double_opt.a,$(top_builddir)/mpblas/optimized/double/.libs/libmplapack_double_opt.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/lin -llin_double -llin_override_double -L$(top_builddir)/mplapack/test/matgen -lmatgen_double -L$(top_builddir)/mplapack/reference -lmplapack_double
opt_libdepends = -L$(top_builddir)/mplapack/test/lin -llin_double_opt -llin_override_double_opt -L$(top_builddir)/mplapack/test/matgen -lmatgen_double_opt -L$(top_builddir)/mpblas/optimized/double -lmplapack_double_opt
endif

mplibs = -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -DMPLAPACK_BUILD_WITH_DOUBLE -I$(top_srcdir)/fable -DMPLAPACK_INTERNAL

xlintstR_double_SOURCES  = ../common/Rchkaa.cpp
xlintstR_double_opt_SOURCES  = ../common/Rchkaa.cpp
xlintstR_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_double_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_double_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_double_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_double_CXXFLAGS = $(cxxflags)
xlintstR_double_opt_CXXFLAGS = $(cxxflags)

xlintstC_double_SOURCES  = ../common/Cchkaa.cpp
xlintstC_double_opt_SOURCES  = ../common/Cchkaa.cpp
xlintstC_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_double_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_double_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_double_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_double_CXXFLAGS = $(cxxflags)
xlintstC_double_opt_CXXFLAGS = $(cxxflags)

xlintstrfR_double_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_double_opt_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_double_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_double_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_double_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_double_CXXFLAGS = $(cxxflags)
xlintstrfR_double_opt_CXXFLAGS = $(cxxflags)

xlintstrfC_double_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_double_opt_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_double_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_double_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_double_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_double_CXXFLAGS = $(cxxflags)
xlintstrfC_double_opt_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_double_testdir) $(prefix)
endif
# -----------------------------------------------------------------------------
# Output directory: CPU_COMPILER_OS  (or CPU_COMPILER_mingw for MinGW builds)
#
# cpu_os_detection.sh prints "CPU|OS" (e.g. "Ryzen_Threadripper_3970X|ubuntu22_04").
# MPLAPACK_CC_TAG is set by AX_MPLAPACK_CC_TAG in configure.ac (e.g. "gcc_14_2_0").
#
# Native layout : Ryzen_Threadripper_3970X_gcc_14_2_0_ubuntu22_04/<test>.out
# MinGW  layout : Ryzen_Threadripper_3970X_gcc_14_2_0_windows/<test>.out
# -----------------------------------------------------------------------------

_CPU_OS_RAW := $(shell bash $(abs_top_srcdir)/misc/cpu_os_detection.sh)
_CPU_TAG    := $(word 1,$(subst |, ,$(_CPU_OS_RAW)))
_OS_TAG     := $(word 2,$(subst |, ,$(_CPU_OS_RAW)))

if IS_MINGW
OUTDIR = $(_CPU_TAG)_windows_$(MPLAPACK_CC_TAG)
else
OUTDIR = $(_CPU_TAG)_$(_OS_TAG)_$(MPLAPACK_CC_TAG)
endif

OUTDIR_PATH = ../results/$(OUTDIR)/double



# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <OUTDIR>/<test>.out)
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

# All output logs: <OUTDIR>/<test>.out
DOUBLE_LIN_OUTS = $(foreach t,$(DOUBLE_LIN_TESTS),$(OUTDIR_PATH)/$(t).out)
CLEANFILES = $(DOUBLE_LIN_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(DOUBLE_LIN_OUTS): FORCE

check-local: double-check

.PRECIOUS: $(DOUBLE_LIN_OUTS)

double-check: $$(DOUBLE_LIN_OUTS)
	@echo "double lin tests done."
	@echo "Output directory: $(OUTDIR_PATH)"

clean-local:
	rm -f $(DOUBLE_LIN_OUTS) *~ \#*
	@if test -d "$(OUTDIR_PATH)"; then rmdir -p "$(OUTDIR_PATH)" 2>/dev/null || true; fi

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

$(OUTDIR_PATH)/Rtest.out: $(srcdir)/../Rtest.in xlintstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK linear equation routines"; \
	set -e; { \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xlintstR_double$(EXEEXT) < $(srcdir)/../Rtest.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ctest.out: $(srcdir)/../Ctest.in xlintstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK linear equation routines"; \
	set -e; { \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xlintstC_double$(EXEEXT) < $(srcdir)/../Ctest.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rtest_rfp.out: $(srcdir)/../Rtest_rfp.in xlintstrfR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfR_double$(EXEEXT) < $(srcdir)/../Rtest_rfp.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ctest_rfp.out: $(srcdir)/../Ctest_rfp.in xlintstrfC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfC_double$(EXEEXT) < $(srcdir)/../Ctest_rfp.in; \
	} > $@ 2>&1
