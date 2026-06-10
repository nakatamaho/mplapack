check_PROGRAMS = $(mplapack_mpfr_test_PROGRAMS)

mplapack_mpfr_test_PROGRAMS = xeigtstR_mpfr xeigtstC_mpfr xdmdeigtstR_mpfr xdmdeigtstC_mpfr
mplapack_mpfr_testdir = $(prefix)/lib/$(target)/mplapack/test/eig

# -----------------------------------------------------------------------------
# MPFR emulation profiles
#
# MPFR uses 0.5 <= m < 1 for the mantissa, so IEEE-like emulation uses:
#   binary64  : emin=-1021,  emax=1024,  prec=53
#   binary128 : emin=-16381, emax=16384, prec=113
#
# Add a new profile by:
#   1) appending its name to MPFR_PROFILES
#   2) defining MPFR_ENV_<n> and MPFR_PROFILE_LABEL_<n>
# -----------------------------------------------------------------------------

MPFR_PROFILES = default binary64 binary128

MPFR_ENV_default   =
MPFR_ENV_binary64  = MPLAPACK_MPFR_EMAX=1024  MPLAPACK_MPFR_EMIN=-1021  MPLAPACK_MPFR_PRECISION=53
MPFR_ENV_binary128 = MPLAPACK_MPFR_EMAX=16384 MPLAPACK_MPFR_EMIN=-16381 MPLAPACK_MPFR_PRECISION=113

MPFR_PROFILE_LABEL_default   = default
MPFR_PROFILE_LABEL_binary64  = binary64 (prec=53 emin=-1021 emax=1024)
MPFR_PROFILE_LABEL_binary128 = binary128 (prec=113 emin=-16381 emax=16384)

$(check_PROGRAMS): ../libeig_mpfr.la ../libeig_override_mpfr.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_mpfr.a,$(top_builddir)/mplapack/test/eig/.libs/libeig_override_mpfr.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_mpfr.a,$(top_builddir)/mpblas/reference/.libs/libmpblas_mpfr.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_mpfr.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/eig -leig_mpfr -leig_override_mpfr -L$(top_builddir)/mplapack/test/matgen -lmatgen_mpfr -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr
endif

mplibs = -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpc -lmpfr -lgmpxx -lgmp -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/fable

xeigtstR_mpfr_SOURCES  = ../common/Rchkee.cpp
xeigtstR_mpfr_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_mpfr_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_mpfr_CXXFLAGS = $(cxxflags)

xeigtstC_mpfr_SOURCES  = ../common/Cchkee.cpp
xeigtstC_mpfr_LDFLAGS  = $(libdepends) $(mplibs) # -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_mpfr_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_mpfr_CXXFLAGS = $(cxxflags)

xdmdeigtstR_mpfr_SOURCES  = ../common/Rchkdmd.cpp
xdmdeigtstR_mpfr_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstR_mpfr_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstR_mpfr_CXXFLAGS = $(cxxflags)

xdmdeigtstC_mpfr_SOURCES  = ../common/Cchkdmd.cpp
xdmdeigtstC_mpfr_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstC_mpfr_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstC_mpfr_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_mpfr_testdir) $(prefix)
endif

# -----------------------------------------------------------------------------
# Output directory: CPU_OS/COMPILER_TAG
#
# cpu_os_detection.sh prints "CPU|OS" (e.g. "Ryzen_Threadripper_3970X|ubuntu22_04").
# MPLAPACK_CC_TAG is set by AX_MPLAPACK_CC_TAG in configure.ac (e.g. "gcc_14_2_0").
#
# Resulting layout:
#   Ryzen_Threadripper_3970X_ubuntu22_04/gcc_14_2_0/<test>.<profile>.out
# -----------------------------------------------------------------------------

_CPU_OS_RAW  := $(shell bash $(abs_top_srcdir)/misc/cpu_os_detection.sh)
_CPU_TAG     := $(word 1,$(subst |, ,$(_CPU_OS_RAW)))
_OS_TAG      := $(word 2,$(subst |, ,$(_CPU_OS_RAW)))

if IS_MINGW
OUTDIR = $(_CPU_TAG)_windows_$(MPLAPACK_CC_TAG)
else
OUTDIR = $(_CPU_TAG)_$(_OS_TAG)_$(MPLAPACK_CC_TAG)
endif

OUTDIR_PATH = ../results/$(OUTDIR)/mpfr


# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <OUTDIR>/<test>.<profile>.out)
# -----------------------------------------------------------------------------

MPFR_REAL_TESTS = \
	Rnep Rsep Rse2 Rsvd Rec Red Rgg Rgd Rsb Rsg Rbal Rbak Rgbal Rgbak Rbb Rglm Rgqr Rgsv Rcsd Rlse Rdmd

MPFR_COMPLEX_TESTS = \
	Cnep Csep Cse2 Csvd Cec Ced Cgg Cgd Csb Csg Cbal Cbak Cgbal Cgbak Cbb Cglm Cgqr Cgsv Ccsd Clse Cdmd

MPFR_ALL_TESTS = $(MPFR_REAL_TESTS) $(MPFR_COMPLEX_TESTS)

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
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
# Hook into Automake's "make check" (profile-abstracted; parallel like LIN)
# -----------------------------------------------------------------------------

# Default: run all profiles in "make check"
MPFR_CHECK_TARGETS = mpfr-check-all

# GNU make: allow deferred expansion of prerequisites (Automake may reorder vars/rules)
.SECONDEXPANSION:

# Force test execution even if *.out already exists
.PHONY: mpfr-check mpfr-check-all FORCE
FORCE:

# All output logs: <OUTDIR>/<test>.<profile>.out
MPFR_EIG_OUTS = $(foreach t,$(MPFR_ALL_TESTS),$(foreach p,$(MPFR_PROFILES),$(OUTDIR_PATH)/$(t).$(p).out))
CLEANFILES = $(MPFR_EIG_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(MPFR_EIG_OUTS): FORCE

check-local: mpfr-check

.PRECIOUS: $(MPFR_EIG_OUTS)

mpfr-check: $(MPFR_CHECK_TARGETS)
	@echo "MPFR eig tests done: $(MPFR_CHECK_TARGETS)"
	@echo "Output directory: $(OUTDIR_PATH)"

# Run all profiles listed in MPFR_PROFILES
mpfr-check-all: $$(addprefix mpfr-check-,$$(MPFR_PROFILES))
	@echo "MPFR eig tests done for profiles: $(MPFR_PROFILES)"

# One rule for any profile name (default/binary64/binary128/...)
mpfr-check-%: $$(foreach t,$$(MPFR_ALL_TESTS),$(OUTDIR_PATH)/$$(t).$$*.out)
	@echo "MPFR profile '$*' done."

clean-local:
	rm -f $(MPFR_EIG_OUTS) *~ \#*
	@if test -d "$(OUTDIR_PATH)"; then rmdir -p "$(OUTDIR_PATH)" 2>/dev/null || true; fi

# -----------------------------------------------------------------------------
# Test rules (profile-abstracted; one output per test per profile)
# -----------------------------------------------------------------------------

$(OUTDIR_PATH)/Rnep.%.out: $(srcdir)/../nep.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../nep.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rsep.%.out: $(srcdir)/../sep.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../sep.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rse2.%.out: $(srcdir)/../se2.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../se2.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rsvd.%.out: $(srcdir)/../svd.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../svd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rec.%.out: $(srcdir)/../Rec.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../Rec.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Red.%.out: $(srcdir)/../Red.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../Red.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgg.%.out: $(srcdir)/../Rgg.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGS: Testing Generalized Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../Rgg.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgd.%.out: $(srcdir)/../Rgd.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routine [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../Rgd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rsb.%.out: $(srcdir)/../Rsb.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SBB: Testing banded Schur decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../Rsb.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rsg.%.out: $(srcdir)/../Rsg.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SGG: Testing Generalized Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../Rsg.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rbal.%.out: $(srcdir)/../Rbal.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "RGEBAL: Testing the balancing of a REAL general matrix [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../Rbal.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rbak.%.out: $(srcdir)/../Rbak.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "RGEBAK: Testing the back transformation of a REAL balanced matrix [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../Rbak.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgbal.%.out: $(srcdir)/../Rgbal.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "RGGBAL: Testing the balancing of a pair of REAL general matrices [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../Rgbal.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgbak.%.out: $(srcdir)/../Rgbak.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "RGGBAK: Testing the back transformation of a pair of REAL balanced matrices [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../Rgbak.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rbb.%.out: $(srcdir)/../Rbb.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "RBB: Testing banded Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../Rbb.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rglm.%.out: $(srcdir)/../glm.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../glm.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgqr.%.out: $(srcdir)/../gqr.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../gqr.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgsv.%.out: $(srcdir)/../gsv.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../gsv.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rcsd.%.out: $(srcdir)/../csd.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../csd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rlse.%.out: $(srcdir)/../lse.in xeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_mpfr$(EXEEXT) < $(srcdir)/../lse.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cnep.%.out: $(srcdir)/../nep.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../nep.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Csep.%.out: $(srcdir)/../sep.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../sep.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cse2.%.out: $(srcdir)/../se2.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../se2.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Csvd.%.out: $(srcdir)/../svd.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../svd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cec.%.out: $(srcdir)/../Cec.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../Cec.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ced.%.out: $(srcdir)/../Ced.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../Ced.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgg.%.out: $(srcdir)/../Cgg.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGS: Testing Generalized Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../Cgg.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgd.%.out: $(srcdir)/../Cgd.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routine [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../Cgd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Csb.%.out: $(srcdir)/../Csb.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SBB: Testing banded Schur decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../Csb.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Csg.%.out: $(srcdir)/../Csg.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SGG: Testing Generalized Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../Csg.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cbal.%.out: $(srcdir)/../Cbal.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../Cbal.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cbak.%.out: $(srcdir)/../Cbak.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../Cbak.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgbal.%.out: $(srcdir)/../Cgbal.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../Cgbal.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgbak.%.out: $(srcdir)/../Cgbak.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../Cgbak.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cbb.%.out: $(srcdir)/../Cbb.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../Cbb.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cglm.%.out: $(srcdir)/../glm.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../glm.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgqr.%.out: $(srcdir)/../gqr.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../gqr.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgsv.%.out: $(srcdir)/../gsv.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../gsv.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ccsd.%.out: $(srcdir)/../csd.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../csd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Clse.%.out: $(srcdir)/../lse.in xeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_mpfr$(EXEEXT) < $(srcdir)/../lse.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rdmd.%.out: $(srcdir)/../Rdmd.in xdmdeigtstR_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DMD: Testing Dynamic Mode Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xdmdeigtstR_mpfr$(EXEEXT) < $(srcdir)/../Rdmd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cdmd.%.out: $(srcdir)/../Cdmd.in xdmdeigtstC_mpfr$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DMD: Testing COMPLEX Dynamic Mode Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper $(WINE) ./xdmdeigtstC_mpfr$(EXEEXT) < $(srcdir)/../Cdmd.in; \
	} > $@ 2>&1
