check_PROGRAMS = $(mplapack_dd_test_PROGRAMS) $(mplapack_dd_opt_test_PROGRAMS)

mplapack_dd_test_PROGRAMS = xlintstR_dd xlintstC_dd xlintstrfR_dd xlintstrfC_dd
mplapack_dd_opt_test_PROGRAMS = xlintstR_dd_opt xlintstC_dd_opt xlintstrfR_dd_opt xlintstrfC_dd_opt
mplapack_dd_testdir = $(prefix)/lib/$(target)/mplapack/test/lin
mplapack_dd_opt_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(mplapack_dd_test_PROGRAMS): ../liblin_dd.la ../liblin_override_dd.la
$(mplapack_dd_opt_test_PROGRAMS): ../liblin_dd_opt.la ../liblin_override_dd_opt.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_dd.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_dd.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_dd.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_dd.a,--no-whole-archive
opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_dd_opt.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_dd_opt.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_dd_opt.a,$(top_builddir)/mpblas/optimized/dd/.libs/libmplapack_dd_opt.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/lin -llin_dd -llin_override_dd -L$(top_builddir)/mplapack/test/matgen -lmatgen_dd -L$(top_builddir)/mplapack/reference -lmplapack_dd
opt_libdepends = -L$(top_builddir)/mplapack/test/lin -llin_dd_opt -llin_override_dd_opt -L$(top_builddir)/mplapack/test/matgen -lmatgen_dd_opt -L$(top_builddir)/mpblas/optimized/dd -lmplapack_dd_opt
endif

mplibs = -L$(QD_LIBDIR) -lqd -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_DD @qd_extra_cxxflags@ -DMPLAPACK_INTERNAL -I$(top_srcdir)/fable

xlintstR_dd_SOURCES  = ../common/Rchkaa.cpp
xlintstR_dd_opt_SOURCES  = ../common/Rchkaa.cpp
xlintstR_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_dd_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_dd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_dd_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_dd_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xlintstR_dd_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

xlintstC_dd_SOURCES  = ../common/Cchkaa.cpp
xlintstC_dd_opt_SOURCES  = ../common/Cchkaa.cpp
xlintstC_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_dd_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_dd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_dd_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_dd_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xlintstC_dd_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

xlintstrfR_dd_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_dd_opt_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_dd_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_dd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_dd_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_dd_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xlintstrfR_dd_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

xlintstrfC_dd_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_dd_opt_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_dd_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_dd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_dd_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_dd_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xlintstrfC_dd_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_dd_testdir) $(prefix)
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

OUTDIR_PATH = ../results/$(OUTDIR)/dd



# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <OUTDIR>/<test>.out)
# -----------------------------------------------------------------------------

DD_LIN_TESTS = Rtest Ctest Rtest_rfp Ctest_rfp

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
# -----------------------------------------------------------------------------

VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended (last assignment wins)
DD_EXTRA_ENV ?=

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
.PHONY: dd-check FORCE
FORCE:

# All output logs: <OUTDIR>/<test>.out
DD_LIN_OUTS = $(foreach t,$(DD_LIN_TESTS),$(OUTDIR_PATH)/$(t).out)
CLEANFILES = $(DD_LIN_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(DD_LIN_OUTS): FORCE

check-local: dd-check

.PRECIOUS: $(DD_LIN_OUTS)

dd-check: $$(DD_LIN_OUTS)
	@echo "DD lin tests done."
	@echo "Output directory: $(OUTDIR_PATH)"

clean-local:
	rm -f $(DD_LIN_OUTS) *~ \#*
	@if test -d "$(OUTDIR_PATH)"; then rmdir -p "$(OUTDIR_PATH)" 2>/dev/null || true; fi

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

$(OUTDIR_PATH)/Rtest.out: $(srcdir)/../Rtest.in xlintstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK linear equation routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstR_dd$(EXEEXT) < $(srcdir)/../Rtest.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ctest.out: $(srcdir)/../Ctest.in xlintstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK linear equation routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstC_dd$(EXEEXT) < $(srcdir)/../Ctest.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rtest_rfp.out: $(srcdir)/../Rtest_rfp.in xlintstrfR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfR_dd$(EXEEXT) < $(srcdir)/../Rtest_rfp.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ctest_rfp.out: $(srcdir)/../Ctest_rfp.in xlintstrfC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfC_dd$(EXEEXT) < $(srcdir)/../Ctest_rfp.in; \
	} > $@ 2>&1
