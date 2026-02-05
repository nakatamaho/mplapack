check_PROGRAMS = $(mplapack__Float64x_test_PROGRAMS)

mplapack__Float64x_test_PROGRAMS = xeigtstR_binary80 xeigtstC_binary80
mplapack__Float64x_testdir = $(prefix)/lib/$(target)/mplapack/test/eig

$(check_PROGRAMS): ../libeig_binary80.la ../libeig_override_binary80.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_binary80.a,$(top_builddir)/mplapack/test/eig/.libs/libeig_override_binary80.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_binary80.a,$(top_builddir)/mpblas/reference/.libs/libmpblas_binary80.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_binary80.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/eig -leig_binary80 -leig_override_binary80 -L$(top_builddir)/mplapack/test/matgen -lmatgen_binary80 -L$(top_builddir)/mplapack/reference -lmplapack_binary80 -L$(top_builddir)/mpblas/reference -lmpblas_binary80
endif

mplibs = -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_BINARY80___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/fable

xeigtstR__Float64x_SOURCES  = ../common/Rchkee.cpp
xeigtstR__Float64x_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR__Float64x_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR__Float64x_CXXFLAGS = $(cxxflags)

xeigtstC__Float64x_SOURCES  = ../common/Cchkee.cpp
xeigtstC__Float64x_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC__Float64x_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC__Float64x_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapack__Float64x_testdir) $(prefix)
endif

# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <test>.out)
# -----------------------------------------------------------------------------

FLOAT64X_REAL_TESTS = \
	Rnep Rsep Rse2 Rsvd Rec Red Rgg Rgd Rsb Rsg Rbal Rbak Rgbal Rgbak Rbb Rglm Rgqr Rgsv Rcsd Rlse

FLOAT64X_COMPLEX_TESTS = \
	Cnep Csep Cse2 Csvd Cec Ced Cgg Cgd Csb Csg Cbal Cbak Cgbal Cgbak Cbb Cglm Cgqr Cgsv Ccsd Clse

FLOAT64X_ALL_TESTS = $(FLOAT64X_REAL_TESTS) $(FLOAT64X_COMPLEX_TESTS)

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
# -----------------------------------------------------------------------------

VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended before running the test (last assignment wins)
FLOAT64X_EXTRA_ENV ?=
FLOAT64X_LABEL ?= _Float64x

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
.PHONY: float64x-check FORCE
FORCE:

FLOAT64X_EIG_OUTS = $(foreach t,$(FLOAT64X_ALL_TESTS),$(t).out)
CLEANFILES = $(FLOAT64X_EIG_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(FLOAT64X_EIG_OUTS): FORCE

check-local: float64x-check

.PRECIOUS: $(FLOAT64X_EIG_OUTS)

float64x-check: $(FLOAT64X_EIG_OUTS)
	@echo "_Float64x eig tests done: $(FLOAT64X_ALL_TESTS)"

clean-local:
	rm -f *.out *~ \#*

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

Rnep.out: ../nep.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../nep.in; \
	} > $@ 2>&1

Rsep.out: ../sep.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../sep.in; \
	} > $@ 2>&1

Rse2.out: ../se2.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../se2.in; \
	} > $@ 2>&1

Rsvd.out: ../svd.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../svd.in; \
	} > $@ 2>&1

Rec.out: ../Rec.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../Rec.in; \
	} > $@ 2>&1

Red.out: ../Red.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../Red.in; \
	} > $@ 2>&1

Rgg.out: ../Rgg.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../Rgg.in; \
	} > $@ 2>&1

Rgd.out: ../Rgd.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../Rgd.in; \
	} > $@ 2>&1

Rsb.out: ../Rsb.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSB: Testing Symmetric Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../Rsb.in; \
	} > $@ 2>&1

Rsg.out: ../Rsg.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../Rsg.in; \
	} > $@ 2>&1

Rbal.out: ../Rbal.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebal: Testing the balancing of a general matrix [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../Rbal.in; \
	} > $@ 2>&1

Rbak.out: ../Rbak.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebak: Testing the back transformation of balanced matrix [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../Rbak.in; \
	} > $@ 2>&1

Rgbal.out: ../Rgbal.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbal: Testing the balancing of a pair of general matrices [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../Rgbal.in; \
	} > $@ 2>&1

Rgbak.out: ../Rgbak.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbak: Testing the back transformation of a pair of balanced matrices [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../Rgbak.in; \
	} > $@ 2>&1

Rbb.out: ../Rbb.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DBB: Testing banded Singular Value Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../Rbb.in; \
	} > $@ 2>&1

Rglm.out: ../glm.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../glm.in; \
	} > $@ 2>&1

Rgqr.out: ../gqr.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../gqr.in; \
	} > $@ 2>&1

Rgsv.out: ../gsv.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../gsv.in; \
	} > $@ 2>&1

Rcsd.out: ../csd.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../csd.in; \
	} > $@ 2>&1

Rlse.out: ../lse.in xeigtstR_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary80$(EXEEXT) < ../lse.in; \
	} > $@ 2>&1

Cnep.out: ../nep.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../nep.in; \
	} > $@ 2>&1

Csep.out: ../sep.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../sep.in; \
	} > $@ 2>&1

Cse2.out: ../se2.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../se2.in; \
	} > $@ 2>&1

Csvd.out: ../svd.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../svd.in; \
	} > $@ 2>&1

Cec.out: ../Cec.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CEC: Testing COMPLEX Eigen Condition Routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../Cec.in; \
	} > $@ 2>&1

Ced.out: ../Ced.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CES: Testing COMPLEX Nonsymmetric Schur Form Driver [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../Ced.in; \
	} > $@ 2>&1

Cgg.out: ../Cgg.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGG: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../Cgg.in; \
	} > $@ 2>&1

Cgd.out: ../Cgd.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem driver routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../Cgd.in; \
	} > $@ 2>&1

Csb.out: ../Csb.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CHB: Testing Hermitian Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../Csb.in; \
	} > $@ 2>&1

Csg.out: ../Csg.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../Csg.in; \
	} > $@ 2>&1

Cbal.out: ../Cbal.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../Cbal.in; \
	} > $@ 2>&1

Cbak.out: ../Cbak.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../Cbak.in; \
	} > $@ 2>&1

Cgbal.out: ../Cgbal.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../Cgbal.in; \
	} > $@ 2>&1

Cgbak.out: ../Cgbak.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../Cgbak.in; \
	} > $@ 2>&1

Cbb.out: ../Cbb.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../Cbb.in; \
	} > $@ 2>&1

Cglm.out: ../glm.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../glm.in; \
	} > $@ 2>&1

Cgqr.out: ../gqr.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../gqr.in; \
	} > $@ 2>&1

Cgsv.out: ../gsv.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../gsv.in; \
	} > $@ 2>&1
Ccsd.out: ../csd.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../csd.in; \
	} > $@ 2>&1

Clse.out: ../lse.in xeigtstC_binary80$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary80$(EXEEXT) < ../lse.in; \
	} > $@ 2>&1
