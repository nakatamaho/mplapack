check_PROGRAMS = $(mplapack_qd_test_PROGRAMS)

mplapack_qd_test_PROGRAMS = xeigtstR_qd xeigtstC_qd
mplapack_qd_testdir = $(prefix)/lib/$(target)/mplapack/test/eig

$(check_PROGRAMS): ../libeig_qd.la ../libeig_override_qd.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_qd.a,$(top_builddir)/mplapack/test/eig/.libs/libeig_override_qd.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_qd.a,$(top_builddir)/mpblas/reference/.libs/libmpblas_qd.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_qd.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/eig -leig_qd -leig_override_qd -L$(top_builddir)/mplapack/test/matgen -lmatgen_qd -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(top_builddir)/mpblas/reference -lmpblas_qd
endif

mplibs = -L$(QD_LIBDIR) -lqd -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/fable

xeigtstR_qd_SOURCES  = ../common/Rchkee.cpp
xeigtstR_qd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_qd_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_qd_CXXFLAGS = $(cxxflags)

xeigtstC_qd_SOURCES  = ../common/Cchkee.cpp
xeigtstC_qd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_qd_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_qd_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapack_qd_testdir) $(prefix)
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

_CPU_OS_RAW := $(shell bash $(abs_top_builddir)/misc/cpu_os_detection.sh)
_CPU_TAG    := $(word 1,$(subst |, ,$(_CPU_OS_RAW)))
_OS_TAG     := $(word 2,$(subst |, ,$(_CPU_OS_RAW)))

if IS_MINGW
OUTDIR = $(_CPU_TAG)_$(MPLAPACK_CC_TAG)_windows
else
OUTDIR = $(_CPU_TAG)_$(MPLAPACK_CC_TAG)_$(_OS_TAG)
endif

$(OUTDIR):
	$(MKDIR_P) $@


# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <OUTDIR>/<test>.out)
# -----------------------------------------------------------------------------

QD_REAL_TESTS = \
	Rnep Rsep Rse2 Rsvd Rec Red Rgg Rgd Rsb Rsg Rbal Rbak Rgbal Rgbak Rbb Rglm Rgqr Rgsv Rcsd Rlse

QD_COMPLEX_TESTS = \
	Cnep Csep Cse2 Csvd Cec Ced Cgg Cgd Csb Csg Cbal Cbak Cgbal Cgbak Cbb Cglm Cgqr Cgsv Ccsd Clse

QD_ALL_TESTS = $(QD_REAL_TESTS) $(QD_COMPLEX_TESTS)

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
# -----------------------------------------------------------------------------

VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended before running the test (last assignment wins)
QD_EXTRA_ENV ?=
QD_LABEL ?= qd

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
.PHONY: qd-check FORCE
FORCE:

QD_EIG_OUTS = $(foreach t,$(QD_ALL_TESTS),$(OUTDIR)/$(t).out)
CLEANFILES = $(QD_EIG_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(QD_EIG_OUTS): FORCE

check-local: qd-check

.PRECIOUS: $(QD_EIG_OUTS)

qd-check: $(QD_EIG_OUTS)
	@echo "$(QD_LABEL) eig tests done: $(QD_ALL_TESTS)"
	@echo "Output directory: $(OUTDIR)"

clean-local:
	rm -f $(QD_EIG_OUTS) *~ \#*
	@if test -d "$(OUTDIR)"; then rmdir -p "$(OUTDIR)" 2>/dev/null || true; fi

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

$(OUTDIR)/Rnep.out: ../nep.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../nep.in; \
	} > $@ 2>&1
$(OUTDIR)/Rsep.out: ../sep.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../sep.in; \
	} > $@ 2>&1
$(OUTDIR)/Rse2.out: ../se2.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../se2.in; \
	} > $@ 2>&1
$(OUTDIR)/Rsvd.out: ../svd.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../svd.in; \
	} > $@ 2>&1
$(OUTDIR)/Rec.out: ../Rec.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../Rec.in; \
	} > $@ 2>&1
$(OUTDIR)/Red.out: ../Red.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../Red.in; \
	} > $@ 2>&1
$(OUTDIR)/Rgg.out: ../Rgg.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../Rgg.in; \
	} > $@ 2>&1
$(OUTDIR)/Rgd.out: ../Rgd.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../Rgd.in; \
	} > $@ 2>&1
$(OUTDIR)/Rsb.out: ../Rsb.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSB: Testing Symmetric Eigenvalue Problem routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../Rsb.in; \
	} > $@ 2>&1
$(OUTDIR)/Rsg.out: ../Rsg.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../Rsg.in; \
	} > $@ 2>&1
$(OUTDIR)/Rbal.out: ../Rbal.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebal: Testing the balancing of a general matrix [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../Rbal.in; \
	} > $@ 2>&1
$(OUTDIR)/Rbak.out: ../Rbak.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebak: Testing the back transformation of balanced matrix [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../Rbak.in; \
	} > $@ 2>&1
$(OUTDIR)/Rgbal.out: ../Rgbal.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbal: Testing the balancing of a pair of general matrices [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../Rgbal.in; \
	} > $@ 2>&1
$(OUTDIR)/Rgbak.out: ../Rgbak.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbak: Testing the back transformation of a pair of balanced matrices [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../Rgbak.in; \
	} > $@ 2>&1
$(OUTDIR)/Rbb.out: ../Rbb.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DBB: Testing banded Singular Value Decomposition routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../Rbb.in; \
	} > $@ 2>&1
$(OUTDIR)/Rglm.out: ../glm.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../glm.in; \
	} > $@ 2>&1
$(OUTDIR)/Rgqr.out: ../gqr.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../gqr.in; \
	} > $@ 2>&1
$(OUTDIR)/Rgsv.out: ../gsv.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../gsv.in; \
	} > $@ 2>&1
$(OUTDIR)/Rcsd.out: ../csd.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../csd.in; \
	} > $@ 2>&1
$(OUTDIR)/Rlse.out: ../lse.in xeigtstR_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_qd$(EXEEXT) < ../lse.in; \
	} > $@ 2>&1
$(OUTDIR)/Cnep.out: ../nep.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../nep.in; \
	} > $@ 2>&1
$(OUTDIR)/Csep.out: ../sep.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../sep.in; \
	} > $@ 2>&1
$(OUTDIR)/Cse2.out: ../se2.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../se2.in; \
	} > $@ 2>&1
$(OUTDIR)/Csvd.out: ../svd.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../svd.in; \
	} > $@ 2>&1
$(OUTDIR)/Cec.out: ../Cec.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CEC: Testing COMPLEX Eigen Condition Routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../Cec.in; \
	} > $@ 2>&1
$(OUTDIR)/Ced.out: ../Ced.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CES: Testing COMPLEX Nonsymmetric Schur Form Driver [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../Ced.in; \
	} > $@ 2>&1
$(OUTDIR)/Cgg.out: ../Cgg.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGG: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../Cgg.in; \
	} > $@ 2>&1
$(OUTDIR)/Cgd.out: ../Cgd.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem driver routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../Cgd.in; \
	} > $@ 2>&1
$(OUTDIR)/Csb.out: ../Csb.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CHB: Testing Hermitian Eigenvalue Problem routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../Csb.in; \
	} > $@ 2>&1
$(OUTDIR)/Csg.out: ../Csg.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../Csg.in; \
	} > $@ 2>&1
$(OUTDIR)/Cbal.out: ../Cbal.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../Cbal.in; \
	} > $@ 2>&1
$(OUTDIR)/Cbak.out: ../Cbak.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../Cbak.in; \
	} > $@ 2>&1
$(OUTDIR)/Cgbal.out: ../Cgbal.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../Cgbal.in; \
	} > $@ 2>&1
$(OUTDIR)/Cgbak.out: ../Cgbak.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../Cgbak.in; \
	} > $@ 2>&1
$(OUTDIR)/Cbb.out: ../Cbb.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../Cbb.in; \
	} > $@ 2>&1
$(OUTDIR)/Cglm.out: ../glm.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../glm.in; \
	} > $@ 2>&1
$(OUTDIR)/Cgqr.out: ../gqr.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../gqr.in; \
	} > $@ 2>&1
$(OUTDIR)/Cgsv.out: ../gsv.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../gsv.in; \
	} > $@ 2>&1
$(OUTDIR)/Ccsd.out: ../csd.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../csd.in; \
	} > $@ 2>&1
$(OUTDIR)/Clse.out: ../lse.in xeigtstC_qd$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(QD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(QD_LABEL) ==="; \
		$(QD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_qd$(EXEEXT) < ../lse.in; \
	} > $@ 2>&1
