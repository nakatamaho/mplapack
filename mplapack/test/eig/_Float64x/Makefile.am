check_PROGRAMS = $(mplapack__Float64x_test_PROGRAMS)

mplapack__Float64x_test_PROGRAMS = xeigtstR__Float64x xeigtstC__Float64x
mplapack__Float64x_testdir = $(prefix)/lib/$(target)/mplapack/test/eig

$(check_PROGRAMS): ../libeig__Float64x.la ../libeig_override__Float64x.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig__Float64x.a,$(top_builddir)/mplapack/test/eig/.libs/libeig_override__Float64x.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen__Float64x.a,$(top_builddir)/mpblas/reference/.libs/libmpblas__Float64x.a,$(top_builddir)/mplapack/reference/.libs/libmplapack__Float64x.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/eig -leig__Float64x -leig_override__Float64x -L$(top_builddir)/mplapack/test/matgen -lmatgen__Float64x -L$(top_builddir)/mplapack/reference -lmplapack__Float64x -L$(top_builddir)/mpblas/reference -lmpblas__Float64x
endif

mplibs = -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH__FLOAT64X___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/fable

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

Rnep.out: ../nep.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../nep.in; \
	} > $@ 2>&1

Rsep.out: ../sep.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../sep.in; \
	} > $@ 2>&1

Rse2.out: ../se2.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../se2.in; \
	} > $@ 2>&1

Rsvd.out: ../svd.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../svd.in; \
	} > $@ 2>&1

Rec.out: ../Rec.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../Rec.in; \
	} > $@ 2>&1

Red.out: ../Red.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../Red.in; \
	} > $@ 2>&1

Rgg.out: ../Rgg.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../Rgg.in; \
	} > $@ 2>&1

Rgd.out: ../Rgd.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../Rgd.in; \
	} > $@ 2>&1

Rsb.out: ../Rsb.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSB: Testing Symmetric Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../Rsb.in; \
	} > $@ 2>&1

Rsg.out: ../Rsg.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../Rsg.in; \
	} > $@ 2>&1

Rbal.out: ../Rbal.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebal: Testing the balancing of a general matrix [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../Rbal.in; \
	} > $@ 2>&1

Rbak.out: ../Rbak.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebak: Testing the back transformation of balanced matrix [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../Rbak.in; \
	} > $@ 2>&1

Rgbal.out: ../Rgbal.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbal: Testing the balancing of a pair of general matrices [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../Rgbal.in; \
	} > $@ 2>&1

Rgbak.out: ../Rgbak.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbak: Testing the back transformation of a pair of balanced matrices [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../Rgbak.in; \
	} > $@ 2>&1

Rbb.out: ../Rbb.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DBB: Testing banded Singular Value Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../Rbb.in; \
	} > $@ 2>&1

Rglm.out: ../glm.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../glm.in; \
	} > $@ 2>&1

Rgqr.out: ../gqr.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../gqr.in; \
	} > $@ 2>&1

Rgsv.out: ../gsv.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../gsv.in; \
	} > $@ 2>&1

Rcsd.out: ../csd.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../csd.in; \
	} > $@ 2>&1

Rlse.out: ../lse.in xeigtstR__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR__Float64x$(EXEEXT) < ../lse.in; \
	} > $@ 2>&1

Cnep.out: ../nep.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../nep.in; \
	} > $@ 2>&1

Csep.out: ../sep.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../sep.in; \
	} > $@ 2>&1

Cse2.out: ../se2.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../se2.in; \
	} > $@ 2>&1

Csvd.out: ../svd.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../svd.in; \
	} > $@ 2>&1

Cec.out: ../Cec.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CEC: Testing COMPLEX Eigen Condition Routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../Cec.in; \
	} > $@ 2>&1

Ced.out: ../Ced.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CES: Testing COMPLEX Nonsymmetric Schur Form Driver [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../Ced.in; \
	} > $@ 2>&1

Cgg.out: ../Cgg.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGG: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../Cgg.in; \
	} > $@ 2>&1

Cgd.out: ../Cgd.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem driver routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../Cgd.in; \
	} > $@ 2>&1

Csb.out: ../Csb.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CHB: Testing Hermitian Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../Csb.in; \
	} > $@ 2>&1

Csg.out: ../Csg.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../Csg.in; \
	} > $@ 2>&1

Cbal.out: ../Cbal.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../Cbal.in; \
	} > $@ 2>&1

Cbak.out: ../Cbak.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../Cbak.in; \
	} > $@ 2>&1

Cgbal.out: ../Cgbal.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../Cgbal.in; \
	} > $@ 2>&1

Cgbak.out: ../Cgbak.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../Cgbak.in; \
	} > $@ 2>&1

Cbb.out: ../Cbb.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../Cbb.in; \
	} > $@ 2>&1

Cglm.out: ../glm.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../glm.in; \
	} > $@ 2>&1

Cgqr.out: ../gqr.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../gqr.in; \
	} > $@ 2>&1

Cgsv.out: ../gsv.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../gsv.in; \
	} > $@ 2>&1
Ccsd.out: ../csd.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../csd.in; \
	} > $@ 2>&1

Clse.out: ../lse.in xeigtstC__Float64x$(EXEEXT)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(FLOAT64X_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT64X_LABEL) ==="; \
		$(FLOAT64X_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC__Float64x$(EXEEXT) < ../lse.in; \
	} > $@ 2>&1
