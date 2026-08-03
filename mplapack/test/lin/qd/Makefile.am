check_PROGRAMS = $(mplapack_qd_test_PROGRAMS) $(mplapack_qd_opt_test_PROGRAMS)

mplapack_qd_test_PROGRAMS = xlintstR_qd xlintstC_qd xlintstrfR_qd xlintstrfC_qd
mplapack_qd_opt_test_PROGRAMS = xlintstR_qd_opt xlintstC_qd_opt xlintstrfR_qd_opt xlintstrfC_qd_opt
mplapack_qd_testdir = $(prefix)/lib/$(target)/mplapack/test/lin
mplapack_qd_opt_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(mplapack_qd_test_PROGRAMS): ../liblin_qd.la ../liblin_override_qd.la
$(mplapack_qd_opt_test_PROGRAMS): ../liblin_qd_opt.la ../liblin_override_qd_opt.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_qd.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_qd.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_qd.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_qd.a,--no-whole-archive
opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_qd_opt.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_qd_opt.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_qd_opt.a,$(top_builddir)/mpblas/optimized/qd/.libs/libmplapack_qd_opt.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/lin -llin_qd -llin_override_qd -L$(top_builddir)/mplapack/test/matgen -lmatgen_qd -L$(top_builddir)/mplapack/reference -lmplapack_qd
opt_libdepends = -L$(top_builddir)/mplapack/test/lin -llin_qd_opt -llin_override_qd_opt -L$(top_builddir)/mplapack/test/matgen -lmatgen_qd_opt -L$(top_builddir)/mpblas/optimized/qd -lmplapack_qd_opt
endif

mplibs = -L$(QD_LIBDIR) -lqd $(QD_RPATH_LDFLAGS) -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_QD @qd_extra_cxxflags@ -DMPLAPACK_INTERNAL -I$(top_srcdir)/fable

xlintstR_qd_SOURCES  = ../common/Rchkaa.cpp
xlintstR_qd_opt_SOURCES  = ../common/Rchkaa.cpp
xlintstR_qd_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_qd_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_qd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_qd_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_qd_CXXFLAGS = $(cxxflags)
xlintstR_qd_opt_CXXFLAGS = $(cxxflags)

xlintstC_qd_SOURCES  = ../common/Cchkaa.cpp
xlintstC_qd_opt_SOURCES  = ../common/Cchkaa.cpp
xlintstC_qd_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_qd_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_qd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_qd_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_qd_CXXFLAGS = $(cxxflags)
xlintstC_qd_opt_CXXFLAGS = $(cxxflags)

xlintstrfR_qd_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_qd_opt_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_qd_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_qd_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_qd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_qd_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_qd_CXXFLAGS = $(cxxflags)
xlintstrfR_qd_opt_CXXFLAGS = $(cxxflags)

xlintstrfC_qd_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_qd_opt_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_qd_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_qd_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_qd_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_qd_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_qd_CXXFLAGS = $(cxxflags)
xlintstrfC_qd_opt_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_qd_testdir) $(prefix)
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

OUTDIR_PATH = ../results/$(OUTDIR)/qd



# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <OUTDIR>/<test>.out)
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

# All output logs: <OUTDIR>/<test>.out
QD_LIN_OUTS = $(foreach t,$(QD_LIN_TESTS),$(OUTDIR_PATH)/$(t).out)
CLEANFILES = $(QD_LIN_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(QD_LIN_OUTS): FORCE

check-local: qd-check

.PRECIOUS: $(QD_LIN_OUTS)

qd-check: $$(QD_LIN_OUTS)
	@echo "QD lin tests done."
	@echo "Output directory: $(OUTDIR_PATH)"

clean-local:
	rm -f $(QD_LIN_OUTS) *~ \#*
	@if test -d "$(OUTDIR_PATH)"; then rmdir -p "$(OUTDIR_PATH)" 2>/dev/null || true; fi

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

$(OUTDIR_PATH)/Rtest.out: $(srcdir)/../Rtest.in xlintstR_qd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK linear equation routines"; \
	set -e; { \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstR_qd$(EXEEXT) < $(srcdir)/../Rtest.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ctest.out: $(srcdir)/../Ctest.in xlintstC_qd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK linear equation routines"; \
	set -e; { \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstC_qd$(EXEEXT) < $(srcdir)/../Ctest.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rtest_rfp.out: $(srcdir)/../Rtest_rfp.in xlintstrfR_qd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfR_qd$(EXEEXT) < $(srcdir)/../Rtest_rfp.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ctest_rfp.out: $(srcdir)/../Ctest_rfp.in xlintstrfC_qd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfC_qd$(EXEEXT) < $(srcdir)/../Ctest_rfp.in; \
	} > $@ 2>&1
