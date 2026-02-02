check_PROGRAMS = $(mplapack_dd_test_PROGRAMS)

mplapack_dd_test_PROGRAMS = xeigtstR_dd xeigtstC_dd
mplapack_dd_testdir = $(prefix)/lib/$(target)/mplapack/test/eig

$(check_PROGRAMS): ../libeig_dd.la ../libeig_override_dd.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_dd.a,$(top_builddir)/mplapack/test/eig/.libs/libeig_override_dd.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_dd.a,$(top_builddir)/mpblas/reference/.libs/libmpblas_dd.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_dd.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/eig -leig_dd -leig_override_dd -L$(top_builddir)/mplapack/test/matgen -lmatgen_dd -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(top_builddir)/mpblas/reference -lmpblas_dd
endif

mplibs = -L$(QD_LIBDIR) -lqd -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/fable

xeigtstR_dd_SOURCES  = ../common/Rchkee.cpp
xeigtstR_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_dd_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_dd_CXXFLAGS = $(cxxflags)

xeigtstC_dd_SOURCES  = ../common/Cchkee.cpp
xeigtstC_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_dd_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_dd_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapack_dd_testdir) $(prefix)
endif

# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <test>.out)
# -----------------------------------------------------------------------------

DD_REAL_TESTS = \
	Rnep Rsep Rse2 Rsvd Rec Red Rgg Rgd Rsb Rsg Rbal Rbak Rgbal Rgbak Rbb Rglm Rgqr Rgsv Rcsd Rlse

DD_COMPLEX_TESTS = \
	Cnep Csep Cse2 Csvd Cec Ced Cgg Cgd Csb Csg Cbal Cbak Cgbal Cgbak Cbb Cglm Cgqr Cgsv Ccsd Clse

DD_ALL_TESTS = $(DD_REAL_TESTS) $(DD_COMPLEX_TESTS)

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

# All output logs: <test>.out
DD_EIG_OUTS = $(foreach t,$(DD_ALL_TESTS),$(t).out)
CLEANFILES = $(DD_EIG_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(DD_EIG_OUTS): FORCE

check-local: dd-check

.PRECIOUS: $(DD_EIG_OUTS)

dd-check: $$(DD_EIG_OUTS)
	@echo "DD eig tests done."

clean-local:
	rm -f *.out *~ \#*

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

Rnep.out: ../nep.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../nep.in; \
	} > $@ 2>&1

Rsep.out: ../sep.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../sep.in; \
	} > $@ 2>&1

Rse2.out: ../se2.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../se2.in; \
	} > $@ 2>&1

Rsvd.out: ../svd.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../svd.in; \
	} > $@ 2>&1

Rec.out: ../Rec.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../Rec.in; \
	} > $@ 2>&1

Red.out: ../Red.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../Red.in; \
	} > $@ 2>&1

Rgg.out: ../Rgg.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGS: Testing Generalized Nonsymmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../Rgg.in; \
	} > $@ 2>&1

Rgd.out: ../Rgd.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routine"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../Rgd.in; \
	} > $@ 2>&1

Rsb.out: ../Rsb.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SBB: Testing banded Schur decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../Rsb.in; \
	} > $@ 2>&1

Rsg.out: ../Rsg.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SGG: Testing Generalized Symmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../Rsg.in; \
	} > $@ 2>&1

Rbal.out: ../Rbal.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGEBAL: Testing the balancing of a real general matrix"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../Rbal.in; \
	} > $@ 2>&1

Rbak.out: ../Rbak.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGEBAK: Testing the back transformation of a real balanced matrix"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../Rbak.in; \
	} > $@ 2>&1

Rgbal.out: ../Rgbal.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGGBAL: Testing the balancing of a pair of real general matrices"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../Rgbal.in; \
	} > $@ 2>&1

Rgbak.out: ../Rgbak.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGGBAK: Testing the back transformation of a pair of real balanced matrices"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../Rgbak.in; \
	} > $@ 2>&1

Rbb.out: ../Rbb.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DBB: Testing banded Singular Value Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../Rbb.in; \
	} > $@ 2>&1

Rglm.out: ../glm.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../glm.in; \
	} > $@ 2>&1

Rgqr.out: ../gqr.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../gqr.in; \
	} > $@ 2>&1

Rgsv.out: ../gsv.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../gsv.in; \
	} > $@ 2>&1

Rcsd.out: ../csd.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../csd.in; \
	} > $@ 2>&1

Rlse.out: ../lse.in xeigtstR_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < ../lse.in; \
	} > $@ 2>&1

Cnep.out: ../nep.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../nep.in; \
	} > $@ 2>&1

Csep.out: ../sep.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../sep.in; \
	} > $@ 2>&1

Cse2.out: ../se2.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../se2.in; \
	} > $@ 2>&1

Csvd.out: ../svd.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../svd.in; \
	} > $@ 2>&1

Cec.out: ../Cec.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../Cec.in; \
	} > $@ 2>&1

Ced.out: ../Ced.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../Ced.in; \
	} > $@ 2>&1

Cgg.out: ../Cgg.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGS: Testing Generalized Nonsymmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../Cgg.in; \
	} > $@ 2>&1

Cgd.out: ../Cgd.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routine"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../Cgd.in; \
	} > $@ 2>&1

Csb.out: ../Csb.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SBB: Testing banded Schur decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../Csb.in; \
	} > $@ 2>&1

Csg.out: ../Csg.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SGG: Testing Generalized Symmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../Csg.in; \
	} > $@ 2>&1

Cbal.out: ../Cbal.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../Cbal.in; \
	} > $@ 2>&1

Cbak.out: ../Cbak.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../Cbak.in; \
	} > $@ 2>&1

Cgbal.out: ../Cgbal.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../Cgbal.in; \
	} > $@ 2>&1

Cgbak.out: ../Cgbak.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../Cgbak.in; \
	} > $@ 2>&1

Cbb.out: ../Cbb.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../Cbb.in; \
	} > $@ 2>&1

Cglm.out: ../glm.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../glm.in; \
	} > $@ 2>&1

Cgqr.out: ../gqr.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../gqr.in; \
	} > $@ 2>&1

Cgsv.out: ../gsv.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../gsv.in; \
	} > $@ 2>&1

Ccsd.out: ../csd.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../csd.in; \
	} > $@ 2>&1

Clse.out: ../lse.in xeigtstC_dd$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < ../lse.in; \
	} > $@ 2>&1
