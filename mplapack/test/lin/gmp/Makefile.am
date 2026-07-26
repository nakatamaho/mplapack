check_PROGRAMS = $(mplapack_gmp_test_PROGRAMS) $(mplapack_gmp_opt_test_PROGRAMS)

mplapack_gmp_test_PROGRAMS = xlintstR_gmp xlintstC_gmp xlintstrfR_gmp xlintstrfC_gmp
mplapack_gmp_opt_test_PROGRAMS = xlintstR_gmp_opt xlintstC_gmp_opt xlintstrfR_gmp_opt xlintstrfC_gmp_opt
mplapack_gmp_testdir = $(prefix)/lib/$(target)/mplapack/test/lin
mplapack_gmp_opt_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(mplapack_gmp_test_PROGRAMS): ../liblin_gmp.la ../liblin_override_gmp.la
$(mplapack_gmp_opt_test_PROGRAMS): ../liblin_gmp_opt.la ../liblin_override_gmp_opt.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_gmp.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_gmp.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_gmp.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_gmp.a,--no-whole-archive
opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_gmp_opt.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_gmp_opt.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_gmp_opt.a,$(top_builddir)/mpblas/optimized/gmp/.libs/libmplapack_gmp_opt.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/lin -llin_gmp -llin_override_gmp -L$(top_builddir)/mplapack/test/matgen -lmatgen_gmp -L$(top_builddir)/mplapack/reference -lmplapack_gmp
opt_libdepends = -L$(top_builddir)/mplapack/test/lin -llin_gmp_opt -llin_override_gmp_opt -L$(top_builddir)/mplapack/test/matgen -lmatgen_gmp_opt -L$(top_builddir)/mpblas/optimized/gmp -lmplapack_gmp_opt
endif

mplibs = -L$(GMP_LIBDIR) -lgmp -lm
cxxflags = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/fable

xlintstR_gmp_SOURCES  = ../common/Rchkaa.cpp
xlintstR_gmp_opt_SOURCES  = ../common/Rchkaa.cpp
xlintstR_gmp_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_gmp_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_gmp_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_gmp_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_gmp_CXXFLAGS = $(cxxflags)
xlintstR_gmp_opt_CXXFLAGS = $(cxxflags)

xlintstC_gmp_SOURCES  = ../common/Cchkaa.cpp
xlintstC_gmp_opt_SOURCES  = ../common/Cchkaa.cpp
xlintstC_gmp_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_gmp_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_gmp_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_gmp_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_gmp_CXXFLAGS = $(cxxflags)
xlintstC_gmp_opt_CXXFLAGS = $(cxxflags)

xlintstrfR_gmp_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_gmp_opt_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_gmp_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_gmp_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_gmp_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_gmp_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_gmp_CXXFLAGS = $(cxxflags)
xlintstrfR_gmp_opt_CXXFLAGS = $(cxxflags)

xlintstrfC_gmp_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_gmp_opt_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_gmp_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_gmp_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_gmp_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_gmp_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_gmp_CXXFLAGS = $(cxxflags)
xlintstrfC_gmp_opt_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_gmp_testdir) $(prefix)
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

OUTDIR_PATH = ../results/$(OUTDIR)/gmp



# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <OUTDIR>/<test>.out)
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

# All output logs: <OUTDIR>/<test>.out
GMP_LIN_OUTS = $(foreach t,$(GMP_LIN_TESTS),$(OUTDIR_PATH)/$(t).out)
CLEANFILES = $(GMP_LIN_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(GMP_LIN_OUTS): FORCE

check-local: gmp-check

.PRECIOUS: $(GMP_LIN_OUTS)

gmp-check: $$(GMP_LIN_OUTS)
	@echo "GMP lin tests done."
	@echo "Output directory: $(OUTDIR_PATH)"

clean-local:
	rm -f $(GMP_LIN_OUTS) *~ \#*
	@if test -d "$(OUTDIR_PATH)"; then rmdir -p "$(OUTDIR_PATH)" 2>/dev/null || true; fi

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

$(OUTDIR_PATH)/Rtest.out: $(srcdir)/../Rtest.in xlintstR_gmp$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK linear equation routines"; \
	set -e; { \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xlintstR_gmp$(EXEEXT) < $(srcdir)/../Rtest.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ctest.out: $(srcdir)/../Ctest.in xlintstC_gmp$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK linear equation routines"; \
	set -e; { \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xlintstC_gmp$(EXEEXT) < $(srcdir)/../Ctest.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rtest_rfp.out: $(srcdir)/../Rtest_rfp.in xlintstrfR_gmp$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfR_gmp$(EXEEXT) < $(srcdir)/../Rtest_rfp.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ctest_rfp.out: $(srcdir)/../Ctest_rfp.in xlintstrfC_gmp$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK RFP prototype linear equation routines"; \
	set -e; { \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfC_gmp$(EXEEXT) < $(srcdir)/../Ctest_rfp.in; \
	} > $@ 2>&1
