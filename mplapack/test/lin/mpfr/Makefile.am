check_PROGRAMS = $(mplapack_mpfr_test_PROGRAMS) $(mplapack_mpfr_opt_test_PROGRAMS)

mplapack_mpfr_test_PROGRAMS = xlintstR_mpfr xlintstC_mpfr xlintstrfR_mpfr xlintstrfC_mpfr
mplapack_mpfr_opt_test_PROGRAMS = xlintstR_mpfr_opt xlintstC_mpfr_opt xlintstrfR_mpfr_opt xlintstrfC_mpfr_opt
mplapack_mpfr_testdir = $(prefix)/lib/$(target)/mplapack/test/lin
mplapack_mpfr_opt_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(mplapack_mpfr_test_PROGRAMS): ../liblin_mpfr.la ../liblin_override_mpfr.la
$(mplapack_mpfr_opt_test_PROGRAMS): ../liblin_mpfr_opt.la ../liblin_override_mpfr_opt.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_mpfr.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_mpfr.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_mpfr.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_mpfr.a,--no-whole-archive
opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_mpfr_opt.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_mpfr_opt.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_mpfr_opt.a,$(top_builddir)/mpblas/optimized/mpfr/.libs/libmplapack_mpfr_opt.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/lin -llin_mpfr -llin_override_mpfr -L$(top_builddir)/mplapack/test/matgen -lmatgen_mpfr -L$(top_builddir)/mplapack/reference -lmplapack_mpfr
opt_libdepends = -L$(top_builddir)/mplapack/test/lin -llin_mpfr_opt -llin_override_mpfr_opt -L$(top_builddir)/mplapack/test/matgen -lmatgen_mpfr_opt -L$(top_builddir)/mpblas/optimized/mpfr -lmplapack_mpfr_opt
endif

mplibs = -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpc -lmpfr -lgmp -lm
cxxflags = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_MPFR___ -D___MPLAPACK_INTERNAL___

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_mpfr_testdir) $(prefix)
endif
# -----------------------------------------------------------------------------
# Output directory: CPU_COMPILER_OS  (or CPU_COMPILER_mingw for MinGW builds)
#
# cpu_os_detection.sh prints "CPU|OS" (e.g. "Ryzen_Threadripper_3970X|ubuntu22_04").
# MPLAPACK_CC_TAG is set by AX_MPLAPACK_CC_TAG in configure.ac (e.g. "gcc_14_2_0").
#
# Native layout : Ryzen_Threadripper_3970X_gcc_14_2_0_ubuntu22_04/<test>.<profile>.out
# MinGW  layout : Ryzen_Threadripper_3970X_gcc_14_2_0_windows/<test>.<profile>.out
# -----------------------------------------------------------------------------

_CPU_OS_RAW := $(shell bash $(abs_top_srcdir)/misc/cpu_os_detection.sh)
_CPU_TAG    := $(word 1,$(subst |, ,$(_CPU_OS_RAW)))
_OS_TAG     := $(word 2,$(subst |, ,$(_CPU_OS_RAW)))

if IS_MINGW
OUTDIR = $(_CPU_TAG)_windows_$(MPLAPACK_CC_TAG)
else
OUTDIR = $(_CPU_TAG)_$(_OS_TAG)_$(MPLAPACK_CC_TAG)
endif

OUTDIR_PATH = ../results/$(OUTDIR)/mpfr




xlintstR_mpfr_SOURCES  = ../common/Rchkaa.cpp
xlintstR_mpfr_opt_SOURCES  = ../common/Rchkaa.cpp
xlintstR_mpfr_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_mpfr_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_mpfr_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_mpfr_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_mpfr_CXXFLAGS = $(cxxflags)
xlintstR_mpfr_opt_CXXFLAGS = $(cxxflags)

xlintstC_mpfr_SOURCES  = ../common/Cchkaa.cpp
xlintstC_mpfr_opt_SOURCES  = ../common/Cchkaa.cpp
xlintstC_mpfr_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_mpfr_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_mpfr_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_mpfr_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_mpfr_CXXFLAGS = $(cxxflags)
xlintstC_mpfr_opt_CXXFLAGS = $(cxxflags)

xlintstrfR_mpfr_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_mpfr_opt_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_mpfr_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_mpfr_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_mpfr_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_mpfr_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_mpfr_CXXFLAGS = $(cxxflags)
xlintstrfR_mpfr_opt_CXXFLAGS = $(cxxflags)

xlintstrfC_mpfr_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_mpfr_opt_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_mpfr_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_mpfr_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_mpfr_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_mpfr_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_mpfr_CXXFLAGS = $(cxxflags)
xlintstrfC_mpfr_opt_CXXFLAGS = $(cxxflags)

# -----------------------------------------------------------------------------
# MPFR emulation profiles
# -----------------------------------------------------------------------------

MPFR_PROFILES = default binary64 binary128

MPFR_ENV_default   =
MPFR_ENV_binary64  = MPLAPACK_MPFR_EMAX=1024  MPLAPACK_MPFR_EMIN=-1021  MPLAPACK_MPFR_PRECISION=53
MPFR_ENV_binary128 = MPLAPACK_MPFR_EMAX=16384 MPLAPACK_MPFR_EMIN=-16381 MPLAPACK_MPFR_PRECISION=113

MPFR_PROFILE_LABEL_default   = default
MPFR_PROFILE_LABEL_binary64  = binary64 (prec=53 emin=-1021 emax=1024)
MPFR_PROFILE_LABEL_binary128 = binary128 (prec=113 emin=-16381 emax=16384)

# -----------------------------------------------------------------------------
# Optional valgrind (toggle from command line)
# -----------------------------------------------------------------------------
VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended after profile env (last assignment wins)
MPFR_EXTRA_ENV ?=

# -----------------------------------------------------------------------------
# Cross-compilation support (wine for MinGW)
# -----------------------------------------------------------------------------
if IS_MINGW
WINE = wine
else
WINE =
endif

# -----------------------------------------------------------------------------
# Hook into Automake's "make check"
# -----------------------------------------------------------------------------

# Default: run all profiles in "make check"
MPFR_CHECK_TARGETS ?= mpfr-check-all

# GNU make: allow deferred expansion of prerequisites (Automake may reorder vars/rules)
.SECONDEXPANSION:

# Always re-run tests on "make check"
.PHONY: mpfr-check mpfr-check-all FORCE
FORCE:

# Also remove test logs on "make clean"
clean-local:
	rm -f $(MPFR_OUTFILES) *~ \#*
	@if test -d "$(OUTDIR_PATH)"; then rmdir -p "$(OUTDIR_PATH)" 2>/dev/null || true; fi

# Clean test logs on "make clean"
MPFR_OUTFILES = \
	$(MPFR_PROFILES:%=$(OUTDIR_PATH)/Rtest.%.out) \
	$(MPFR_PROFILES:%=$(OUTDIR_PATH)/Ctest.%.out) \
	$(MPFR_PROFILES:%=$(OUTDIR_PATH)/Rtest_rfp.%.out) \
	$(MPFR_PROFILES:%=$(OUTDIR_PATH)/Ctest_rfp.%.out)
CLEANFILES = $(MPFR_OUTFILES)

.PRECIOUS: $(MPFR_OUTFILES)

check-local: mpfr-check

mpfr-check: $(MPFR_CHECK_TARGETS)
	@echo "MPFR lin tests done: $(MPFR_CHECK_TARGETS)"
	@echo "Output directory: $(OUTDIR_PATH)"

# Run all profiles listed in MPFR_PROFILES
mpfr-check-all: $$(addprefix mpfr-check-,$$(MPFR_PROFILES))
	@echo "MPFR lin tests done for profiles: $(MPFR_PROFILES)"

# One rule for any profile name (default/binary64/binary128/...)
mpfr-check-%: $(OUTDIR_PATH)/Rtest.%.out $(OUTDIR_PATH)/Ctest.%.out $(OUTDIR_PATH)/Rtest_rfp.%.out $(OUTDIR_PATH)/Ctest_rfp.%.out
	@echo "MPFR profile '$*' done."

# -----------------------------------------------------------------------------
# Test rules (profile-abstracted)
# -----------------------------------------------------------------------------

$(OUTDIR_PATH)/Rtest.%.out: $(srcdir)/../Rtest.in xlintstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK linear equation routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xlintstR_mpfr$(EXEEXT) < $(srcdir)/../Rtest.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ctest.%.out: $(srcdir)/../Ctest.in xlintstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK linear equation routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xlintstC_mpfr$(EXEEXT) < $(srcdir)/../Ctest.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rtest_rfp.%.out: $(srcdir)/../Rtest_rfp.in xlintstrfR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Real MPLAPACK RFP prototype linear equation routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfR_mpfr$(EXEEXT) < $(srcdir)/../Rtest_rfp.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ctest_rfp.%.out: $(srcdir)/../Ctest_rfp.in xlintstrfC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Testing Complex MPLAPACK RFP prototype linear equation routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xlintstrfC_mpfr$(EXEEXT) < $(srcdir)/../Ctest_rfp.in; \
	} > $@ 2>&1
