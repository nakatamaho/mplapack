check_PROGRAMS = $(mplapack_binary128_test_PROGRAMS)

mplapack_binary128_test_PROGRAMS = xlintstR_binary128 xlintstC_binary128 xlintstrfR_binary128 xlintstrfC_binary128
mplapack_binary128_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(check_PROGRAMS): ../liblin_binary128.la ../liblin_override_binary128.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_binary128.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_binary128.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_binary128.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_binary128.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/lin -llin_binary128 -llin_override_binary128 -L$(top_builddir)/mplapack/test/matgen -lmatgen_binary128 -L$(top_builddir)/mplapack/reference -lmplapack_binary128
endif

if BINARY128_USE_QUADMATH
mplibs = -lm -lquadmath
else
mplibs = -lm
endif

cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_BINARY128___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/fable

xlintstR_binary128_SOURCES  = ../common/Rchkaa.cpp
xlintstR_binary128_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_binary128_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_binary128_CXXFLAGS = $(cxxflags)

xlintstC_binary128_SOURCES  = ../common/Cchkaa.cpp
xlintstC_binary128_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_binary128_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_binary128_CXXFLAGS = $(cxxflags)

xlintstrfR_binary128_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_binary128_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_binary128_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_binary128_CXXFLAGS = $(cxxflags)

xlintstrfC_binary128_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_binary128_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_binary128_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_binary128_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_binary128_testdir) $(prefix)
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

OUTDIR_PATH = ../results/$(OUTDIR)/binary128



# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <OUTDIR>/<test>.out)
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

# All output logs: <OUTDIR>/<test>.out
F128_LIN_OUTS = $(foreach t,$(F128_LIN_TESTS),$(OUTDIR_PATH)/$(t).out)
CLEANFILES = $(F128_LIN_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(F128_LIN_OUTS): FORCE

check-local: f128-check

.PRECIOUS: $(F128_LIN_OUTS)

f128-check: $$(F128_LIN_OUTS)
	@echo "binary128 lin tests done."
	@echo "Output directory: $(OUTDIR_PATH)"

clean-local:
	rm -f $(F128_LIN_OUTS) *~ \#*
	@if test -d "$(OUTDIR_PATH)"; then rmdir -p "$(OUTDIR_PATH)" 2>/dev/null || true; fi

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

$(OUTDIR_PATH)/Rtest.out: $(srcdir)/../Rtest.in xlintstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK linear equation routines"; \
	set -e; { \
		$(F128_EXTRA_ENV) $$wrapper $(WINE) ./xlintstR_binary128$(EXEEXT) < $(srcdir)/../Rtest.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ctest.out: $(srcdir)/../Ctest.in xlintstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK linear equation routines"; \
	set -e; { \
		$(F128_EXTRA_ENV) $$wrapper $(WINE) ./xlintstC_binary128$(EXEEXT) < $(srcdir)/../Ctest.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rtest_rfp.out: $(srcdir)/../Rtest_rfp.in xlintstrfR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(F128_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfR_binary128$(EXEEXT) < $(srcdir)/../Rtest_rfp.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ctest_rfp.out: $(srcdir)/../Ctest_rfp.in xlintstrfC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(F128_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfC_binary128$(EXEEXT) < $(srcdir)/../Ctest_rfp.in; \
	} > $@ 2>&1
