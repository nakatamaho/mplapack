check_PROGRAMS = $(mplapack_td_test_PROGRAMS) $(mplapack_td_opt_test_PROGRAMS)

# Own the test routines and common arrays in the executable only.
# Do not also link an override library: its finalizers would run twice.
TEST_OVERRIDE_SOURCES = ../common/Mxerbla.cpp ../common/iMlaenv.cpp ../common/Mxlaenv.cpp

mplapack_td_test_PROGRAMS = xlintstR_td xlintstC_td xlintstrfR_td xlintstrfC_td
mplapack_td_opt_test_PROGRAMS = xlintstR_td_opt xlintstC_td_opt xlintstrfR_td_opt xlintstrfC_td_opt
mplapack_td_testdir = $(prefix)/lib/$(target)/mplapack/test/lin
mplapack_td_opt_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(mplapack_td_test_PROGRAMS): ../liblin_td.la
$(mplapack_td_opt_test_PROGRAMS): ../liblin_td_opt.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_td.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_td.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_td.a,--no-whole-archive
opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_td_opt.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_td_opt.a,$(top_builddir)/mpblas/optimized/td/.libs/libmplapack_td_opt.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/lin -llin_td -L$(top_builddir)/mplapack/test/matgen -lmatgen_td -L$(top_builddir)/mplapack/reference -lmplapack_td
opt_libdepends = -L$(top_builddir)/mplapack/test/lin -llin_td_opt -L$(top_builddir)/mplapack/test/matgen -lmatgen_td_opt -L$(top_builddir)/mpblas/optimized/td -lmplapack_td_opt
endif

mplibs = -L$(QD_LIBDIR) -lqd $(QD_RPATH_LDFLAGS) -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_TD @qd_extra_cxxflags@ -DMPLAPACK_INTERNAL -I$(top_srcdir)/fable

xlintstR_td_SOURCES  = ../common/Rchkaa.cpp $(TEST_OVERRIDE_SOURCES)
xlintstR_td_opt_SOURCES  = ../common/Rchkaa.cpp $(TEST_OVERRIDE_SOURCES)
xlintstR_td_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_td_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_td_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_td_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_td_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xlintstR_td_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

xlintstC_td_SOURCES  = ../common/Cchkaa.cpp $(TEST_OVERRIDE_SOURCES)
xlintstC_td_opt_SOURCES  = ../common/Cchkaa.cpp $(TEST_OVERRIDE_SOURCES)
xlintstC_td_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_td_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_td_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_td_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_td_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xlintstC_td_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

xlintstrfR_td_SOURCES  = ../common/Rchkrfp.cpp $(TEST_OVERRIDE_SOURCES)
xlintstrfR_td_opt_SOURCES  = ../common/Rchkrfp.cpp $(TEST_OVERRIDE_SOURCES)
xlintstrfR_td_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_td_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_td_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_td_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_td_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xlintstrfR_td_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

xlintstrfC_td_SOURCES  = ../common/Cchkrfp.cpp $(TEST_OVERRIDE_SOURCES)
xlintstrfC_td_opt_SOURCES  = ../common/Cchkrfp.cpp $(TEST_OVERRIDE_SOURCES)
xlintstrfC_td_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_td_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_td_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_td_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_td_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xlintstrfC_td_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_td_testdir) $(prefix)
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

OUTDIR_PATH = ../results/$(OUTDIR)/td



# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <OUTDIR>/<test>.out)
# -----------------------------------------------------------------------------

TD_LIN_TESTS = Rtest Ctest Rtest_rfp Ctest_rfp

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
# -----------------------------------------------------------------------------

VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended (last assignment wins)
TD_EXTRA_ENV ?=

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
.PHONY: td-check FORCE
FORCE:

# All output logs: <OUTDIR>/<test>.out
TD_LIN_OUTS = $(foreach t,$(TD_LIN_TESTS),$(OUTDIR_PATH)/$(t).out)
CLEANFILES = $(TD_LIN_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(TD_LIN_OUTS): FORCE

check-local: td-check

.PRECIOUS: $(TD_LIN_OUTS)

td-check: $$(TD_LIN_OUTS)
	@echo "TD lin tests done."
	@echo "Output directory: $(OUTDIR_PATH)"

clean-local:
	rm -f $(TD_LIN_OUTS) *~ \#*
	@if test -d "$(OUTDIR_PATH)"; then rmdir -p "$(OUTDIR_PATH)" 2>/dev/null || true; fi

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

$(OUTDIR_PATH)/Rtest.out: $(srcdir)/../Rtest.in xlintstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK linear equation routines"; \
	set -e; { \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstR_td$(EXEEXT) < $(srcdir)/../Rtest.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ctest.out: $(srcdir)/../Ctest.in xlintstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK linear equation routines"; \
	set -e; { \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstC_td$(EXEEXT) < $(srcdir)/../Ctest.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rtest_rfp.out: $(srcdir)/../Rtest_rfp.in xlintstrfR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfR_td$(EXEEXT) < $(srcdir)/../Rtest_rfp.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ctest_rfp.out: $(srcdir)/../Ctest_rfp.in xlintstrfC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfC_td$(EXEEXT) < $(srcdir)/../Ctest_rfp.in; \
	} > $@ 2>&1
