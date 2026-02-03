check_PROGRAMS = $(mplapack_double_test_PROGRAMS)

mplapack_double_test_PROGRAMS = xeigtstR_double xeigtstC_double
mplapack_double_testdir = $(prefix)/lib/$(target)/mplapack/test/eig

$(check_PROGRAMS): ../libeig_double.la ../libeig_override_double.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_double.a,$(top_builddir)/mplapack/test/eig/.libs/libeig_override_double.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_double.a,$(top_builddir)/mpblas/reference/.libs/libmpblas_double.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_double.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/eig -leig_double -leig_override_double -L$(top_builddir)/mplapack/test/matgen -lmatgen_double -L$(top_builddir)/mplapack/reference -lmplapack_double -L$(top_builddir)/mpblas/reference -lmpblas_double
endif
mplibs = -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DOUBLE___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/fable

xeigtstR_double_SOURCES  = ../common/Rchkee.cpp
xeigtstR_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_double_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_double_CXXFLAGS = $(cxxflags)

xeigtstC_double_SOURCES  = ../common/Cchkee.cpp
xeigtstC_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_double_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_double_CXXFLAGS = $(cxxflags)


install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapack_double_testdir) $(prefix)
endif

# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <test>.out)
# -----------------------------------------------------------------------------

DOUBLE_REAL_TESTS = \
	Rnep Rsep Rse2 Rsvd Rec Red Rgg Rgd Rsb Rsg Rbal Rbak Rgbal Rgbak Rbb Rglm Rgqr Rgsv Rcsd Rlse

DOUBLE_COMPLEX_TESTS = \
	Cnep Csep Cse2 Csvd Cec Ced Cgg Cgd Csb Csg Cbal Cbak Cgbal Cgbak Cbb Cglm Cgqr Cgsv Ccsd Clse

DOUBLE_ALL_TESTS = $(DOUBLE_REAL_TESTS) $(DOUBLE_COMPLEX_TESTS)

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
# -----------------------------------------------------------------------------

VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended before running the test (last assignment wins)
DOUBLE_EXTRA_ENV ?=
DOUBLE_LABEL ?= double

# -----------------------------------------------------------------------------
# Cross-compilation support (wine for MinGW)
# -----------------------------------------------------------------------------
if IS_MINGW
WINE = wine
else
WINE =
endif

# -----------------------------------------------------------------------------
# Hook into Automake's "make check" (parallel like LIN)
# -----------------------------------------------------------------------------

# GNU make: allow deferred expansion of prerequisites (Automake may reorder vars/rules)
.SECONDEXPANSION:

# Force test execution even if *.out already exists
.PHONY: double-check FORCE
FORCE:

DOUBLE_EIG_OUTS = $(foreach t,$(DOUBLE_ALL_TESTS),$(t).out)
CLEANFILES = $(DOUBLE_EIG_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(DOUBLE_EIG_OUTS): FORCE

check-local: double-check

.PRECIOUS: $(DOUBLE_EIG_OUTS)

double-check: $(DOUBLE_EIG_OUTS)
	@echo "$(DOUBLE_LABEL) eig tests done: $(DOUBLE_ALL_TESTS)"

clean-local:
	rm -f *.out *~ \#*

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

Rnep.out: ../nep.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../nep.in; \
	} > $@ 2>&1
Rsep.out: ../sep.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../sep.in; \
	} > $@ 2>&1
Rse2.out: ../se2.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../se2.in; \
	} > $@ 2>&1
Rsvd.out: ../svd.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../svd.in; \
	} > $@ 2>&1
Rec.out: ../Rec.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../Rec.in; \
	} > $@ 2>&1
Red.out: ../Red.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../Red.in; \
	} > $@ 2>&1
Rgg.out: ../Rgg.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../Rgg.in; \
	} > $@ 2>&1
Rgd.out: ../Rgd.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../Rgd.in; \
	} > $@ 2>&1
Rsb.out: ../Rsb.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSB: Testing Symmetric Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../Rsb.in; \
	} > $@ 2>&1
Rsg.out: ../Rsg.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../Rsg.in; \
	} > $@ 2>&1
Rbal.out: ../Rbal_double.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebal: Testing the balancing of a general matrix [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../Rbal_double.in; \
	} > $@ 2>&1
Rbak.out: ../Rbak.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebak: Testing the back transformation of balanced matrix [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../Rbak.in; \
	} > $@ 2>&1
Rgbal.out: ../Rgbal.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbal: Testing the balancing of a pair of general matrices [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../Rgbal.in; \
	} > $@ 2>&1
Rgbak.out: ../Rgbak.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbak: Testing the back transformation of a pair of balanced matrices [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../Rgbak.in; \
	} > $@ 2>&1
Rbb.out: ../Rbb.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DBB: Testing banded Singular Value Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../Rbb.in; \
	} > $@ 2>&1
Rglm.out: ../glm.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../glm.in; \
	} > $@ 2>&1
Rgqr.out: ../gqr.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../gqr.in; \
	} > $@ 2>&1
Rgsv.out: ../gsv.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../gsv.in; \
	} > $@ 2>&1
Rcsd.out: ../csd.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../csd.in; \
	} > $@ 2>&1
Rlse.out: ../lse.in xeigtstR_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < ../lse.in; \
	} > $@ 2>&1
Cnep.out: ../nep.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../nep.in; \
	} > $@ 2>&1
Csep.out: ../sep.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../sep.in; \
	} > $@ 2>&1
Cse2.out: ../se2.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../se2.in; \
	} > $@ 2>&1
Csvd.out: ../svd.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../svd.in; \
	} > $@ 2>&1
Cec.out: ../Cec.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CEC: Testing COMPLEX Eigen Condition Routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../Cec.in; \
	} > $@ 2>&1
Ced.out: ../Ced.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CES: Testing COMPLEX Nonsymmetric Schur Form Driver [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../Ced.in; \
	} > $@ 2>&1
Cgg.out: ../Cgg.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGG: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../Cgg.in; \
	} > $@ 2>&1
Cgd.out: ../Cgd.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem driver routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../Cgd.in; \
	} > $@ 2>&1
Csb.out: ../Csb.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CHB: Testing Hermitian Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../Csb.in; \
	} > $@ 2>&1
Csg.out: ../Csg.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../Csg.in; \
	} > $@ 2>&1
Cbal.out: ../Cbal_double.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../Cbal_double.in; \
	} > $@ 2>&1
Cbak.out: ../Cbak.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../Cbak.in; \
	} > $@ 2>&1
Cgbal.out: ../Cgbal.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../Cgbal.in; \
	} > $@ 2>&1
Cgbak.out: ../Cgbak.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../Cgbak.in; \
	} > $@ 2>&1
Cbb.out: ../Cbb.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../Cbb.in; \
	} > $@ 2>&1
Cglm.out: ../glm.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../glm.in; \
	} > $@ 2>&1
Cgqr.out: ../gqr.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../gqr.in; \
	} > $@ 2>&1
Cgsv.out: ../gsv.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../gsv.in; \
	} > $@ 2>&1
Ccsd.out: ../csd.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../csd.in; \
	} > $@ 2>&1
Clse.out: ../lse.in xeigtstC_double$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < ../lse.in; \
	} > $@ 2>&1
