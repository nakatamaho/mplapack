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
# Hook into Automake's "make check" (select profiles like LIN)
# -----------------------------------------------------------------------------

MPFR_EIG_OUTS_default = \
	Rnep.default.out Rsep.default.out Rse2.default.out Rsvd.default.out Rec.default.out \
	Red.default.out Rgg.default.out Rgd.default.out Rsb.default.out Rsg.default.out \
	Rbal.default.out Rbak.default.out Rgbal.default.out Rgbak.default.out Rbb.default.out \
	Rglm.default.out Rgqr.default.out Rgsv.default.out Rcsd.default.out Rlse.default.out \
	Cnep.default.out Csep.default.out Cse2.default.out Csvd.default.out Cec.default.out \
	Ced.default.out Cgg.default.out Cgd.default.out Csb.default.out Csg.default.out \
	Cbal.default.out Cbak.default.out Cgbal.default.out Cgbak.default.out Cbb.default.out \
	Cglm.default.out Cgqr.default.out Cgsv.default.out Ccsd.default.out Clse.default.out

MPFR_EIG_OUTS_binary64 = \
	Rnep.binary64.out Rsep.binary64.out Rse2.binary64.out Rsvd.binary64.out Rec.binary64.out \
	Red.binary64.out Rgg.binary64.out Rgd.binary64.out Rsb.binary64.out Rsg.binary64.out \
	Rbal.binary64.out Rbak.binary64.out Rgbal.binary64.out Rgbak.binary64.out Rbb.binary64.out \
	Rglm.binary64.out Rgqr.binary64.out Rgsv.binary64.out Rcsd.binary64.out Rlse.binary64.out \
	Cnep.binary64.out Csep.binary64.out Cse2.binary64.out Csvd.binary64.out Cec.binary64.out \
	Ced.binary64.out Cgg.binary64.out Cgd.binary64.out Csb.binary64.out Csg.binary64.out \
	Cbal.binary64.out Cbak.binary64.out Cgbal.binary64.out Cgbak.binary64.out Cbb.binary64.out \
	Cglm.binary64.out Cgqr.binary64.out Cgsv.binary64.out Ccsd.binary64.out Clse.binary64.out

MPFR_EIG_OUTS_binary128 = \
	Rnep.binary128.out Rsep.binary128.out Rse2.binary128.out Rsvd.binary128.out Rec.binary128.out \
	Red.binary128.out Rgg.binary128.out Rgd.binary128.out Rsb.binary128.out Rsg.binary128.out \
	Rbal.binary128.out Rbak.binary128.out Rgbal.binary128.out Rgbak.binary128.out Rbb.binary128.out \
	Rglm.binary128.out Rgqr.binary128.out Rgsv.binary128.out Rcsd.binary128.out Rlse.binary128.out \
	Cnep.binary128.out Csep.binary128.out Cse2.binary128.out Csvd.binary128.out Cec.binary128.out \
	Ced.binary128.out Cgg.binary128.out Cgd.binary128.out Csb.binary128.out Csg.binary128.out \
	Cbal.binary128.out Cbak.binary128.out Cgbal.binary128.out Cgbak.binary128.out Cbb.binary128.out \
	Cglm.binary128.out Cgqr.binary128.out Cgsv.binary128.out Ccsd.binary128.out Clse.binary128.out

MPFR_EIG_OUTS_all = $(MPFR_EIG_OUTS_default) $(MPFR_EIG_OUTS_binary64) $(MPFR_EIG_OUTS_binary128)

.PHONY: mpfr-check mpfr-check-all mpfr-check-default mpfr-check-binary64 mpfr-check-binary128

# Default: run only the default profile during "make check"
MPFR_CHECK_TARGETS ?= mpfr-check-default

check-local: mpfr-check

mpfr-check: $(MPFR_CHECK_TARGETS)
	@echo "MPFR eig tests done: $(MPFR_CHECK_TARGETS)"

mpfr-check-default: $(MPFR_EIG_OUTS_default)
	@echo "MPFR profile 'default' done."

mpfr-check-binary64: $(MPFR_EIG_OUTS_binary64)
	@echo "MPFR profile 'binary64' done."

mpfr-check-binary128: $(MPFR_EIG_OUTS_binary128)
	@echo "MPFR profile 'binary128' done."

mpfr-check-all: $(MPFR_EIG_OUTS_all)
	@echo "MPFR eig tests done for profiles: $(MPFR_PROFILES)"

clean-local:
	rm -f *.out *~ \#*

# -----------------------------------------------------------------------------
# Test rules (double-style; one output per test per profile)
# -----------------------------------------------------------------------------

Rnep.default.out: ../nep.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../nep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../nep.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../nep.in; \
	} > $@ 2>&1

Rnep.binary64.out: ../nep.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../nep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../nep.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../nep.in; \
	} > $@ 2>&1

Rnep.binary128.out: ../nep.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../nep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../nep.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../nep.in; \
	} > $@ 2>&1

Rsep.default.out: ../sep.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../sep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../sep.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../sep.in; \
	} > $@ 2>&1

Rsep.binary64.out: ../sep.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../sep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../sep.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../sep.in; \
	} > $@ 2>&1

Rsep.binary128.out: ../sep.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../sep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../sep.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../sep.in; \
	} > $@ 2>&1

Rse2.default.out: ../se2.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../se2.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../se2.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../se2.in; \
	} > $@ 2>&1

Rse2.binary64.out: ../se2.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../se2.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../se2.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../se2.in; \
	} > $@ 2>&1

Rse2.binary128.out: ../se2.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../se2.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../se2.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../se2.in; \
	} > $@ 2>&1

Rsvd.default.out: ../svd.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../svd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../svd.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../svd.in; \
	} > $@ 2>&1

Rsvd.binary64.out: ../svd.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../svd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../svd.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../svd.in; \
	} > $@ 2>&1

Rsvd.binary128.out: ../svd.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../svd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../svd.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../svd.in; \
	} > $@ 2>&1

Rec.default.out: ../Rec.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rec.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rec.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rec.in; \
	} > $@ 2>&1

Rec.binary64.out: ../Rec.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rec.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rec.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rec.in; \
	} > $@ 2>&1

Rec.binary128.out: ../Rec.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rec.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rec.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rec.in; \
	} > $@ 2>&1

Red.default.out: ../Red.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Red.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Red.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Red.in; \
	} > $@ 2>&1

Red.binary64.out: ../Red.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Red.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Red.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Red.in; \
	} > $@ 2>&1

Red.binary128.out: ../Red.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Red.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Red.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Red.in; \
	} > $@ 2>&1

Rgg.default.out: ../Rgg.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGS: Testing Generalized Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgg.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgg.in; \
	} > $@ 2>&1

Rgg.binary64.out: ../Rgg.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGS: Testing Generalized Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgg.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgg.in; \
	} > $@ 2>&1

Rgg.binary128.out: ../Rgg.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGS: Testing Generalized Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgg.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgg.in; \
	} > $@ 2>&1

Rgd.default.out: ../Rgd.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routine [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgd.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgd.in; \
	} > $@ 2>&1

Rgd.binary64.out: ../Rgd.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routine [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgd.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgd.in; \
	} > $@ 2>&1

Rgd.binary128.out: ../Rgd.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routine [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgd.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgd.in; \
	} > $@ 2>&1

Rsb.default.out: ../Rsb.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SBB: Testing banded Schur decomposition routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsb.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsb.in; \
	} > $@ 2>&1

Rsb.binary64.out: ../Rsb.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SBB: Testing banded Schur decomposition routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsb.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsb.in; \
	} > $@ 2>&1

Rsb.binary128.out: ../Rsb.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SBB: Testing banded Schur decomposition routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsb.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsb.in; \
	} > $@ 2>&1

Rsg.default.out: ../Rsg.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SGG: Testing Generalized Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsg.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsg.in; \
	} > $@ 2>&1

Rsg.binary64.out: ../Rsg.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SGG: Testing Generalized Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsg.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsg.in; \
	} > $@ 2>&1

Rsg.binary128.out: ../Rsg.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SGG: Testing Generalized Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsg.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rsg.in; \
	} > $@ 2>&1

Rbal.default.out: ../Rbal.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGEBAL: Testing the balancing of a real general matrix [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbal.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbal.in; \
	} > $@ 2>&1

Rbal.binary64.out: ../Rbal.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGEBAL: Testing the balancing of a real general matrix [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbal.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbal.in; \
	} > $@ 2>&1

Rbal.binary128.out: ../Rbal.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGEBAL: Testing the balancing of a real general matrix [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbal.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbal.in; \
	} > $@ 2>&1

Rbak.default.out: ../Rbak.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGEBAK: Testing the back transformation of a real balanced matrix [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbak.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbak.in; \
	} > $@ 2>&1

Rbak.binary64.out: ../Rbak.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGEBAK: Testing the back transformation of a real balanced matrix [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbak.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbak.in; \
	} > $@ 2>&1

Rbak.binary128.out: ../Rbak.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGEBAK: Testing the back transformation of a real balanced matrix [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbak.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbak.in; \
	} > $@ 2>&1

Rgbal.default.out: ../Rgbal.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGGBAL: Testing the balancing of a pair of real general matrices [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbal.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbal.in; \
	} > $@ 2>&1

Rgbal.binary64.out: ../Rgbal.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGGBAL: Testing the balancing of a pair of real general matrices [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbal.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbal.in; \
	} > $@ 2>&1

Rgbal.binary128.out: ../Rgbal.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGGBAL: Testing the balancing of a pair of real general matrices [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbal.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbal.in; \
	} > $@ 2>&1

Rgbak.default.out: ../Rgbak.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGGBAK: Testing the back transformation of a pair of real balanced matrices [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbak.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbak.in; \
	} > $@ 2>&1

Rgbak.binary64.out: ../Rgbak.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGGBAK: Testing the back transformation of a pair of real balanced matrices [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbak.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbak.in; \
	} > $@ 2>&1

Rgbak.binary128.out: ../Rgbak.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGGBAK: Testing the back transformation of a pair of real balanced matrices [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbak.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rgbak.in; \
	} > $@ 2>&1

Rbb.default.out: ../Rbb.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DBB: Testing banded Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbb.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbb.in; \
	} > $@ 2>&1

Rbb.binary64.out: ../Rbb.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DBB: Testing banded Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbb.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbb.in; \
	} > $@ 2>&1

Rbb.binary128.out: ../Rbb.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DBB: Testing banded Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbb.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../Rbb.in; \
	} > $@ 2>&1

Rglm.default.out: ../glm.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../glm.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../glm.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../glm.in; \
	} > $@ 2>&1

Rglm.binary64.out: ../glm.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../glm.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../glm.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../glm.in; \
	} > $@ 2>&1

Rglm.binary128.out: ../glm.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../glm.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../glm.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../glm.in; \
	} > $@ 2>&1

Rgqr.default.out: ../gqr.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gqr.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gqr.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gqr.in; \
	} > $@ 2>&1

Rgqr.binary64.out: ../gqr.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gqr.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gqr.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gqr.in; \
	} > $@ 2>&1

Rgqr.binary128.out: ../gqr.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gqr.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gqr.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gqr.in; \
	} > $@ 2>&1

Rgsv.default.out: ../gsv.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gsv.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gsv.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gsv.in; \
	} > $@ 2>&1

Rgsv.binary64.out: ../gsv.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gsv.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gsv.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gsv.in; \
	} > $@ 2>&1

Rgsv.binary128.out: ../gsv.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gsv.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gsv.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../gsv.in; \
	} > $@ 2>&1

Rcsd.default.out: ../csd.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../csd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../csd.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../csd.in; \
	} > $@ 2>&1

Rcsd.binary64.out: ../csd.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../csd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../csd.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../csd.in; \
	} > $@ 2>&1

Rcsd.binary128.out: ../csd.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../csd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../csd.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../csd.in; \
	} > $@ 2>&1

Rlse.default.out: ../lse.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../lse.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../lse.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../lse.in; \
	} > $@ 2>&1

Rlse.binary64.out: ../lse.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../lse.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../lse.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../lse.in; \
	} > $@ 2>&1

Rlse.binary128.out: ../lse.in xeigtstR_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../lse.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../lse.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstR_mpfr < ../lse.in; \
	} > $@ 2>&1

Cnep.default.out: ../nep.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../nep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../nep.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../nep.in; \
	} > $@ 2>&1

Cnep.binary64.out: ../nep.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../nep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../nep.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../nep.in; \
	} > $@ 2>&1

Cnep.binary128.out: ../nep.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../nep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../nep.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../nep.in; \
	} > $@ 2>&1

Csep.default.out: ../sep.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../sep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../sep.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../sep.in; \
	} > $@ 2>&1

Csep.binary64.out: ../sep.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../sep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../sep.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../sep.in; \
	} > $@ 2>&1

Csep.binary128.out: ../sep.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../sep.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../sep.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../sep.in; \
	} > $@ 2>&1

Cse2.default.out: ../se2.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../se2.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../se2.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../se2.in; \
	} > $@ 2>&1

Cse2.binary64.out: ../se2.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../se2.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../se2.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../se2.in; \
	} > $@ 2>&1

Cse2.binary128.out: ../se2.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../se2.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../se2.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../se2.in; \
	} > $@ 2>&1

Csvd.default.out: ../svd.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../svd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../svd.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../svd.in; \
	} > $@ 2>&1

Csvd.binary64.out: ../svd.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../svd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../svd.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../svd.in; \
	} > $@ 2>&1

Csvd.binary128.out: ../svd.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../svd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../svd.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../svd.in; \
	} > $@ 2>&1

Cec.default.out: ../Cec.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cec.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cec.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cec.in; \
	} > $@ 2>&1

Cec.binary64.out: ../Cec.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cec.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cec.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cec.in; \
	} > $@ 2>&1

Cec.binary128.out: ../Cec.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cec.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cec.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cec.in; \
	} > $@ 2>&1

Ced.default.out: ../Ced.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Ced.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Ced.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Ced.in; \
	} > $@ 2>&1

Ced.binary64.out: ../Ced.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Ced.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Ced.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Ced.in; \
	} > $@ 2>&1

Ced.binary128.out: ../Ced.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Ced.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Ced.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Ced.in; \
	} > $@ 2>&1

Cgg.default.out: ../Cgg.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGS: Testing Generalized Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgg.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgg.in; \
	} > $@ 2>&1

Cgg.binary64.out: ../Cgg.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGS: Testing Generalized Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgg.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgg.in; \
	} > $@ 2>&1

Cgg.binary128.out: ../Cgg.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGS: Testing Generalized Nonsymmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgg.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgg.in; \
	} > $@ 2>&1

Cgd.default.out: ../Cgd.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routine [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgd.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgd.in; \
	} > $@ 2>&1

Cgd.binary64.out: ../Cgd.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routine [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgd.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgd.in; \
	} > $@ 2>&1

Cgd.binary128.out: ../Cgd.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routine [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgd.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgd.in; \
	} > $@ 2>&1

Csb.default.out: ../Csb.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SBB: Testing banded Schur decomposition routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csb.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csb.in; \
	} > $@ 2>&1

Csb.binary64.out: ../Csb.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SBB: Testing banded Schur decomposition routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csb.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csb.in; \
	} > $@ 2>&1

Csb.binary128.out: ../Csb.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SBB: Testing banded Schur decomposition routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csb.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csb.in; \
	} > $@ 2>&1

Csg.default.out: ../Csg.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SGG: Testing Generalized Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csg.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csg.in; \
	} > $@ 2>&1

Csg.binary64.out: ../Csg.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SGG: Testing Generalized Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csg.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csg.in; \
	} > $@ 2>&1

Csg.binary128.out: ../Csg.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SGG: Testing Generalized Symmetric Eigenvalue Problem routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csg.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csg.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Csg.in; \
	} > $@ 2>&1

Cbal.default.out: ../Cbal.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbal.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbal.in; \
	} > $@ 2>&1

Cbal.binary64.out: ../Cbal.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbal.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbal.in; \
	} > $@ 2>&1

Cbal.binary128.out: ../Cbal.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbal.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbal.in; \
	} > $@ 2>&1

Cbak.default.out: ../Cbak.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbak.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbak.in; \
	} > $@ 2>&1

Cbak.binary64.out: ../Cbak.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbak.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbak.in; \
	} > $@ 2>&1

Cbak.binary128.out: ../Cbak.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbak.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbak.in; \
	} > $@ 2>&1

Cgbal.default.out: ../Cgbal.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbal.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbal.in; \
	} > $@ 2>&1

Cgbal.binary64.out: ../Cgbal.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbal.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbal.in; \
	} > $@ 2>&1

Cgbal.binary128.out: ../Cgbal.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbal.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbal.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbal.in; \
	} > $@ 2>&1

Cgbak.default.out: ../Cgbak.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbak.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbak.in; \
	} > $@ 2>&1

Cgbak.binary64.out: ../Cgbak.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbak.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbak.in; \
	} > $@ 2>&1

Cgbak.binary128.out: ../Cgbak.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbak.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbak.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cgbak.in; \
	} > $@ 2>&1

Cbb.default.out: ../Cbb.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbb.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbb.in; \
	} > $@ 2>&1

Cbb.binary64.out: ../Cbb.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbb.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbb.in; \
	} > $@ 2>&1

Cbb.binary128.out: ../Cbb.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbb.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbb.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../Cbb.in; \
	} > $@ 2>&1

Cglm.default.out: ../glm.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../glm.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../glm.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../glm.in; \
	} > $@ 2>&1

Cglm.binary64.out: ../glm.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../glm.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../glm.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../glm.in; \
	} > $@ 2>&1

Cglm.binary128.out: ../glm.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../glm.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../glm.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../glm.in; \
	} > $@ 2>&1

Cgqr.default.out: ../gqr.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gqr.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gqr.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gqr.in; \
	} > $@ 2>&1

Cgqr.binary64.out: ../gqr.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gqr.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gqr.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gqr.in; \
	} > $@ 2>&1

Cgqr.binary128.out: ../gqr.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gqr.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gqr.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gqr.in; \
	} > $@ 2>&1

Cgsv.default.out: ../gsv.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gsv.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gsv.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gsv.in; \
	} > $@ 2>&1

Cgsv.binary64.out: ../gsv.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gsv.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gsv.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gsv.in; \
	} > $@ 2>&1

Cgsv.binary128.out: ../gsv.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gsv.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gsv.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../gsv.in; \
	} > $@ 2>&1

Ccsd.default.out: ../csd.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../csd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../csd.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../csd.in; \
	} > $@ 2>&1

Ccsd.binary64.out: ../csd.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../csd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../csd.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../csd.in; \
	} > $@ 2>&1

Ccsd.binary128.out: ../csd.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../csd.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../csd.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../csd.in; \
	} > $@ 2>&1

Clse.default.out: ../lse.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(MPFR_PROFILE_LABEL_default)]"; \
	echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../lse.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_default) ==="; \
		echo "CMD: $(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../lse.in"; \
		$(MPFR_ENV_default) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../lse.in; \
	} > $@ 2>&1

Clse.binary64.out: ../lse.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(MPFR_PROFILE_LABEL_binary64)]"; \
	echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../lse.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary64) ==="; \
		echo "CMD: $(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../lse.in"; \
		$(MPFR_ENV_binary64) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../lse.in; \
	} > $@ 2>&1

Clse.binary128.out: ../lse.in xeigtstC_mpfr
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(MPFR_PROFILE_LABEL_binary128)]"; \
	echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../lse.in > $@ 2>&1"; \
	set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_binary128) ==="; \
		echo "CMD: $(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../lse.in"; \
		$(MPFR_ENV_binary128) $(MPFR_EXTRA_ENV) $$wrapper ./xeigtstC_mpfr < ../lse.in; \
	} > $@ 2>&1

