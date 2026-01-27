check_PROGRAMS = $(mplapack_mpfr_test_PROGRAMS)

mplapack_mpfr_test_PROGRAMS = xeigtstR_mpfr xeigtstC_mpfr
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
#   2) defining MPFR_ENV_<name> and MPFR_PROFILE_LABEL_<name>
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
xeigtstR_mpfr_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_mpfr_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_mpfr_CXXFLAGS = $(cxxflags)

xeigtstC_mpfr_SOURCES  = ../common/Cchkee.cpp
xeigtstC_mpfr_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_mpfr_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_mpfr_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapack_mpfr_testdir) $(prefix)
endif


# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <test>.<profile>.out)
# -----------------------------------------------------------------------------

MPFR_REAL_TESTS = \
	Rnep Rsep Rse2 Rsvd Rec Red Rgg Rgd Rsb Rsg Rbal Rbak Rgbal Rgbak Rbb Rglm Rgqr Rgsv Rcsd Rlse

MPFR_COMPLEX_TESTS = \
	Cnep Csep Cse2 Csvd Cec Ced Cgg Cgd Csb Csg Cbal Cbak Cgbal Cgbak Cbb Cglm Cgqr Cgsv Ccsd Clse

MPFR_ALL_TESTS = $(MPFR_REAL_TESTS) $(MPFR_COMPLEX_TESTS)

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
# -----------------------------------------------------------------------------

VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended after profile env (last assignment wins)
MPFR_EXTRA_ENV ?=

# -----------------------------------------------------------------------------
# Hook into Automake's "make check" (profile-abstracted; parallel like LIN)
# -----------------------------------------------------------------------------

# Default: run all profiles in "make check"
MPFR_CHECK_TARGETS ?= mpfr-check-all

# GNU make: allow deferred expansion of prerequisites (Automake may reorder vars/rules)
.SECONDEXPANSION:

# Force test execution even if *.out already exists
.PHONY: mpfr-check mpfr-check-all FORCE
FORCE:

# All output logs: <test>.<profile>.out
MPFR_EIG_OUTS = $(foreach t,$(MPFR_ALL_TESTS),$(foreach p,$(MPFR_PROFILES),$(t).$(p).out))
CLEANFILES = $(MPFR_EIG_OUTS)

check-local: mpfr-check

.PRECIOUS: $(MPFR_EIG_OUTS)

mpfr-check: $(MPFR_CHECK_TARGETS)
	@echo "MPFR eig tests done: $(MPFR_CHECK_TARGETS)"

# Run all profiles listed in MPFR_PROFILES
mpfr-check-all: $$(addprefix mpfr-check-,$$(MPFR_PROFILES))
	@echo "MPFR eig tests done for profiles: $(MPFR_PROFILES)"

# One rule for any profile name (default/binary64/binary128/...)
mpfr-check-%: $$(foreach t,$$(MPFR_ALL_TESTS),$$(t).$$*.out)
	@echo "MPFR profile '$*' done."

clean-local:
	rm -f *.out *~ \#*

# -----------------------------------------------------------------------------
# Test rules (profile-abstracted; one output per test per profile)
# -----------------------------------------------------------------------------

Rnep.%.out: ../nep.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../nep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../nep.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../nep.in; \
	} > $@ 2>&1

Rsep.%.out: ../sep.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../sep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../sep.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../sep.in; \
	} > $@ 2>&1

Rse2.%.out: ../se2.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../se2.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../se2.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../se2.in; \
	} > $@ 2>&1

Rsvd.%.out: ../svd.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../svd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../svd.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../svd.in; \
	} > $@ 2>&1

Rec.%.out: ../Rec.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rec.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rec.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rec.in; \
	} > $@ 2>&1

Red.%.out: ../Red.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Red.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Red.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Red.in; \
	} > $@ 2>&1

Rgg.%.out: ../Rgg.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGS: Testing Generalized Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgg.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgg.in; \
	} > $@ 2>&1

Rgd.%.out: ../Rgd.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routine [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgd.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgd.in; \
	} > $@ 2>&1

Rsb.%.out: ../Rsb.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SBB: Testing banded Schur decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsb.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsb.in; \
	} > $@ 2>&1

Rsg.%.out: ../Rsg.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SGG: Testing Generalized Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsg.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsg.in; \
	} > $@ 2>&1

Rbal.%.out: ../Rbal.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGEBAL: Testing the balancing of a real general matrix [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbal.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbal.in; \
	} > $@ 2>&1

Rbak.%.out: ../Rbak.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGEBAK: Testing the back transformation of a real balanced matrix [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbak.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbak.in; \
	} > $@ 2>&1

Rgbal.%.out: ../Rgbal.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGGBAL: Testing the balancing of a pair of real general matrices [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbal.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbal.in; \
	} > $@ 2>&1

Rgbak.%.out: ../Rgbak.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGGBAK: Testing the back transformation of a pair of real balanced matrices [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbak.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbak.in; \
	} > $@ 2>&1

Rbb.%.out: ../Rbb.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DBB: Testing banded Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbb.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbb.in; \
	} > $@ 2>&1

Rglm.%.out: ../glm.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../glm.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../glm.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../glm.in; \
	} > $@ 2>&1

Rgqr.%.out: ../gqr.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gqr.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gqr.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gqr.in; \
	} > $@ 2>&1

Rgsv.%.out: ../gsv.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gsv.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gsv.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gsv.in; \
	} > $@ 2>&1

Rcsd.%.out: ../csd.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../csd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../csd.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../csd.in; \
	} > $@ 2>&1

Rlse.%.out: ../lse.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../lse.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../lse.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../lse.in; \
	} > $@ 2>&1

Cnep.%.out: ../nep.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../nep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../nep.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../nep.in; \
	} > $@ 2>&1

Csep.%.out: ../sep.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../sep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../sep.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../sep.in; \
	} > $@ 2>&1

Cse2.%.out: ../se2.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../se2.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../se2.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../se2.in; \
	} > $@ 2>&1

Csvd.%.out: ../svd.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../svd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../svd.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../svd.in; \
	} > $@ 2>&1

Cec.%.out: ../Cec.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cec.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cec.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cec.in; \
	} > $@ 2>&1

Ced.%.out: ../Ced.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Ced.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Ced.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Ced.in; \
	} > $@ 2>&1

Cgg.%.out: ../Cgg.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGS: Testing Generalized Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgg.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgg.in; \
	} > $@ 2>&1

Cgd.%.out: ../Cgd.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routine [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgd.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgd.in; \
	} > $@ 2>&1

Csb.%.out: ../Csb.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SBB: Testing banded Schur decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csb.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csb.in; \
	} > $@ 2>&1

Csg.%.out: ../Csg.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SGG: Testing Generalized Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csg.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csg.in; \
	} > $@ 2>&1

Cbal.%.out: ../Cbal.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbal.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbal.in; \
	} > $@ 2>&1

Cbak.%.out: ../Cbak.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbak.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbak.in; \
	} > $@ 2>&1

Cgbal.%.out: ../Cgbal.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbal.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbal.in; \
	} > $@ 2>&1

Cgbak.%.out: ../Cgbak.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbak.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbak.in; \
	} > $@ 2>&1

Cbb.%.out: ../Cbb.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbb.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbb.in; \
	} > $@ 2>&1

Cglm.%.out: ../glm.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../glm.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../glm.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../glm.in; \
	} > $@ 2>&1

Cgqr.%.out: ../gqr.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gqr.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gqr.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gqr.in; \
	} > $@ 2>&1

Cgsv.%.out: ../gsv.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gsv.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gsv.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gsv.in; \
	} > $@ 2>&1

Ccsd.%.out: ../csd.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../csd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../csd.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../csd.in; \
	} > $@ 2>&1

Clse.%.out: ../lse.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(MPFR_PROFILE_LABEL_$*)]"; \
	echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../lse.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		echo "CMD: $(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../lse.in"; \
		$(MPFR_ENV_$*) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../lse.in; \
	} > $@ 2>&1
