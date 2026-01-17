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

MPFR_PROFILES = current binary64 binary128

MPFR_ENV_current  =
MPFR_ENV_binary64  = MPLAPACK_MPFR_EMAX=1024  MPLAPACK_MPFR_EMIN=-1021  MPLAPACK_MPFR_PRECISION=53
MPFR_ENV_binary128 = MPLAPACK_MPFR_EMAX=16384 MPLAPACK_MPFR_EMIN=-16381 MPLAPACK_MPFR_PRECISION=113

MPFR_PROFILE_LABEL_current  = current
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
# "make check": run per profile in parallel
# -----------------------------------------------------------------------------

.PHONY: mpfr-check mpfr-check-%

check: mpfr-check

mpfr-check: $(MPFR_PROFILES:%=mpfr-check-%)
	@echo "MPFR eig tests done for profiles: $(MPFR_PROFILES)"

# Pattern target: builds all <test>.<profile>.out for the given profile.
mpfr-check-%: $(addsuffix .%.out,$(MPFR_ALL_TESTS))
	@echo "MPFR profile '$*' done."

clean-local:
	rm -f *.out *~ \#*

# -----------------------------------------------------------------------------
# Rule generator (keeps recipes identical; avoids copy/paste per profile)
# -----------------------------------------------------------------------------

define MPFR_EIG_RULE
$1.$2.out: $3 $4
	@echo "$5: $6 [$(MPFR_PROFILE_LABEL_$2)]"
	@set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$2) ==="; \
		$$(MPFR_ENV_$2) ./$4 < $3; \
	} > $$@ 2>&1
endef

define GEN_MPFR_EIG_TEST
$$(foreach p,$$(MPFR_PROFILES),$$(eval $$(call MPFR_EIG_RULE,$1,$$(p),$2,$3,$4,$5)))
endef

# -----------------------------------------------------------------------------
# Real tests
# -----------------------------------------------------------------------------

$(eval $(call GEN_MPFR_EIG_TEST,Rnep,../nep.in,xeigtstR_mpfr,NEP,Testing Nonsymmetric Eigenvalue Problem routines))
$(eval $(call GEN_MPFR_EIG_TEST,Rsep,../sep.in,xeigtstR_mpfr,SEP,Testing Symmetric Eigenvalue Problem routines))
$(eval $(call GEN_MPFR_EIG_TEST,Rse2,../se2.in,xeigtstR_mpfr,SEP,Testing Symmetric Eigenvalue Problem routines))
$(eval $(call GEN_MPFR_EIG_TEST,Rsvd,../svd.in,xeigtstR_mpfr,SVD,Testing Singular Value Decomposition routines))
$(eval $(call GEN_MPFR_EIG_TEST,Rec,../Rec.in,xeigtstR_mpfr,DEC,Testing Eigen Condition Routines))
$(eval $(call GEN_MPFR_EIG_TEST,Red,../Red.in,xeigtstR_mpfr,DEV,Testing Nonsymmetric Eigenvalue Driver))
$(eval $(call GEN_MPFR_EIG_TEST,Rgg,../Rgg.in,xeigtstR_mpfr,GGS,Testing Generalized Nonsymmetric Eigenvalue Problem routines))
$(eval $(call GEN_MPFR_EIG_TEST,Rgd,../Rgd.in,xeigtstR_mpfr,DGD,Testing Nonsymmetric Generalized Eigenvalue Problem driver routine))
$(eval $(call GEN_MPFR_EIG_TEST,Rsb,../Rsb.in,xeigtstR_mpfr,SBB,Testing banded Schur decomposition routines))
$(eval $(call GEN_MPFR_EIG_TEST,Rsg,../Rsg.in,xeigtstR_mpfr,SGG,Testing Generalized Symmetric Eigenvalue Problem routines))
$(eval $(call GEN_MPFR_EIG_TEST,Rbal,../Rbal.in,xeigtstR_mpfr,DGEBAL,Testing the balancing of a real general matrix))
$(eval $(call GEN_MPFR_EIG_TEST,Rbak,../Rbak.in,xeigtstR_mpfr,DGEBAK,Testing the back transformation of a real balanced matrix))
$(eval $(call GEN_MPFR_EIG_TEST,Rgbal,../Rgbal.in,xeigtstR_mpfr,DGGBAL,Testing the balancing of a pair of real general matrices))
$(eval $(call GEN_MPFR_EIG_TEST,Rgbak,../Rgbak.in,xeigtstR_mpfr,DGGBAK,Testing the back transformation of a pair of real balanced matrices))
$(eval $(call GEN_MPFR_EIG_TEST,Rbb,../Rbb.in,xeigtstR_mpfr,DBB,Testing banded Singular Value Decomposition routines))
$(eval $(call GEN_MPFR_EIG_TEST,Rglm,../glm.in,xeigtstR_mpfr,GLM,Testing Generalized Linear Regression Model routines))
$(eval $(call GEN_MPFR_EIG_TEST,Rgqr,../gqr.in,xeigtstR_mpfr,GQR,Testing Generalized QR and RQ factorization routines))
$(eval $(call GEN_MPFR_EIG_TEST,Rgsv,../gsv.in,xeigtstR_mpfr,GSV,Testing Generalized Singular Value Decomposition routines))
$(eval $(call GEN_MPFR_EIG_TEST,Rcsd,../csd.in,xeigtstR_mpfr,CSD,Testing CS Decomposition routines))
$(eval $(call GEN_MPFR_EIG_TEST,Rlse,../lse.in,xeigtstR_mpfr,LSE,Testing Constrained Linear Least Squares routines))

# -----------------------------------------------------------------------------
# Complex tests
# -----------------------------------------------------------------------------

$(eval $(call GEN_MPFR_EIG_TEST,Cnep,../nep.in,xeigtstC_mpfr,NEP,Testing Nonsymmetric Eigenvalue Problem routines))
$(eval $(call GEN_MPFR_EIG_TEST,Csep,../sep.in,xeigtstC_mpfr,SEP,Testing Symmetric Eigenvalue Problem routines))
$(eval $(call GEN_MPFR_EIG_TEST,Cse2,../se2.in,xeigtstC_mpfr,SEP,Testing Symmetric Eigenvalue Problem routines))
$(eval $(call GEN_MPFR_EIG_TEST,Csvd,../svd.in,xeigtstC_mpfr,SVD,Testing Singular Value Decomposition routines))
$(eval $(call GEN_MPFR_EIG_TEST,Cec,../Cec.in,xeigtstC_mpfr,DEC,Testing Eigen Condition Routines))
$(eval $(call GEN_MPFR_EIG_TEST,Ced,../Ced.in,xeigtstC_mpfr,DEV,Testing Nonsymmetric Eigenvalue Driver))
$(eval $(call GEN_MPFR_EIG_TEST,Cgg,../Cgg.in,xeigtstC_mpfr,GGS,Testing Generalized Nonsymmetric Eigenvalue Problem routines))
$(eval $(call GEN_MPFR_EIG_TEST,Cgd,../Cgd.in,xeigtstC_mpfr,CGD,Testing Nonsymmetric Generalized Eigenvalue Problem driver routine))
$(eval $(call GEN_MPFR_EIG_TEST,Csb,../Csb.in,xeigtstC_mpfr,SBB,Testing banded Schur decomposition routines))
$(eval $(call GEN_MPFR_EIG_TEST,Csg,../Csg.in,xeigtstC_mpfr,SGG,Testing Generalized Symmetric Eigenvalue Problem routines))
$(eval $(call GEN_MPFR_EIG_TEST,Cbal,../Cbal.in,xeigtstC_mpfr,CGEBAL,Testing the balancing of a COMPLEX general matrix))
$(eval $(call GEN_MPFR_EIG_TEST,Cbak,../Cbak.in,xeigtstC_mpfr,CGEBAK,Testing the back transformation of a COMPLEX balanced matrix))
$(eval $(call GEN_MPFR_EIG_TEST,Cgbal,../Cgbal.in,xeigtstC_mpfr,CGGBAL,Testing the balancing of a pair of COMPLEX general matrices))
$(eval $(call GEN_MPFR_EIG_TEST,Cgbak,../Cgbak.in,xeigtstC_mpfr,CGGBAK,Testing the back transformation of a pair of COMPLEX balanced matrices))
$(eval $(call GEN_MPFR_EIG_TEST,Cbb,../Cbb.in,xeigtstC_mpfr,CBB,Testing banded Singular Value Decomposition routines))
$(eval $(call GEN_MPFR_EIG_TEST,Cglm,../glm.in,xeigtstC_mpfr,GLM,Testing Generalized Linear Regression Model routines))
$(eval $(call GEN_MPFR_EIG_TEST,Cgqr,../gqr.in,xeigtstC_mpfr,GQR,Testing Generalized QR and RQ factorization routines))
$(eval $(call GEN_MPFR_EIG_TEST,Cgsv,../gsv.in,xeigtstC_mpfr,GSV,Testing Generalized Singular Value Decomposition routines))
$(eval $(call GEN_MPFR_EIG_TEST,Ccsd,../csd.in,xeigtstC_mpfr,CSD,Testing CS Decomposition routines))
$(eval $(call GEN_MPFR_EIG_TEST,Clse,../lse.in,xeigtstC_mpfr,LSE,Testing Constrained Linear Least Squares routines))
