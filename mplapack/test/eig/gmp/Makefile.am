check_PROGRAMS = $(mplapack_gmp_test_PROGRAMS)

mplapack_gmp_test_PROGRAMS = xeigtstR_gmp xeigtstC_gmp
mplapack_gmp_testdir = $(prefix)/lib/$(target)/mplapack/test/eig

$(check_PROGRAMS): ../libeig_gmp.la ../libeig_override_gmp.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_gmp.a,$(top_builddir)/mplapack/test/eig/.libs/libeig_override_gmp.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_gmp.a,$(top_builddir)/mpblas/reference/.libs/libmpblas_gmp.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_gmp.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/eig -leig_gmp -leig_override_gmp -L$(top_builddir)/mplapack/test/matgen -lmatgen_gmp -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp 
endif

mplibs = -L$(GMP_LIBDIR) -lgmpxx -lgmp -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/fable

xeigtstR_gmp_SOURCES  = ../common/Rchkee.cpp
xeigtstR_gmp_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_gmp_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_gmp_CXXFLAGS = $(cxxflags)

xeigtstC_gmp_SOURCES  = ../common/Cchkee.cpp
xeigtstC_gmp_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_gmp_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_gmp_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapack_gmp_testdir) $(prefix)
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

GMP_REAL_TESTS = \
	Rnep Rsep Rse2 Rsvd Rec Red Rgg Rgd Rsb Rsg Rbal Rbak Rgbal Rgbak Rbb Rglm Rgqr Rgsv Rlse

GMP_COMPLEX_TESTS = \
	Cnep Csep Cse2 Csvd Cec Ced Cgg Cgd Csb Csg Cbal Cbak Cgbal Cgbak Cbb Cglm Cgqr Cgsv Clse

GMP_ALL_TESTS = $(GMP_REAL_TESTS) $(GMP_COMPLEX_TESTS)

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
# -----------------------------------------------------------------------------

VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended before running the test (last assignment wins)
GMP_EXTRA_ENV ?=
GMP_LABEL ?= gmp

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
.PHONY: gmp-check FORCE
FORCE:

GMP_EIG_OUTS = $(foreach t,$(GMP_ALL_TESTS),$(OUTDIR)/$(t).out)
CLEANFILES = $(GMP_EIG_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(GMP_EIG_OUTS): FORCE

check-local: gmp-check

.PRECIOUS: $(GMP_EIG_OUTS)

gmp-check: $(GMP_EIG_OUTS)
	@echo "$(GMP_LABEL) eig tests done: $(GMP_ALL_TESTS)"
	@echo "Output directory: $(OUTDIR)"

clean-local:
	rm -f $(GMP_EIG_OUTS) *~ \#*
	@if test -d "$(OUTDIR)"; then rmdir -p "$(OUTDIR)" 2>/dev/null || true; fi

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

$(OUTDIR)/Rnep.out: ../nep.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../nep.in; \
	} > $@ 2>&1
$(OUTDIR)/Rsep.out: ../sep.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../sep.in; \
	} > $@ 2>&1
$(OUTDIR)/Rse2.out: ../se2.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../se2.in; \
	} > $@ 2>&1
$(OUTDIR)/Rsvd.out: ../svd.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../svd.in; \
	} > $@ 2>&1
$(OUTDIR)/Rec.out: ../Rec.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../Rec.in; \
	} > $@ 2>&1
$(OUTDIR)/Red.out: ../Red.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../Red.in; \
	} > $@ 2>&1
$(OUTDIR)/Rgg.out: ../Rgg.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../Rgg.in; \
	} > $@ 2>&1
$(OUTDIR)/Rgd.out: ../Rgd.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../Rgd.in; \
	} > $@ 2>&1
$(OUTDIR)/Rsb.out: ../Rsb.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSB: Testing Symmetric Eigenvalue Problem routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../Rsb.in; \
	} > $@ 2>&1
$(OUTDIR)/Rsg.out: ../Rsg.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../Rsg.in; \
	} > $@ 2>&1
$(OUTDIR)/Rbal.out: ../Rbal.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebal: Testing the balancing of a general matrix [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../Rbal.in; \
	} > $@ 2>&1
$(OUTDIR)/Rbak.out: ../Rbak.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebak: Testing the back transformation of balanced matrix [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../Rbak.in; \
	} > $@ 2>&1
$(OUTDIR)/Rgbal.out: ../Rgbal.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbal: Testing the balancing of a pair of general matrices [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../Rgbal.in; \
	} > $@ 2>&1
$(OUTDIR)/Rgbak.out: ../Rgbak.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbak: Testing the back transformation of a pair of balanced matrices [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../Rgbak.in; \
	} > $@ 2>&1
$(OUTDIR)/Rbb.out: ../Rbb.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DBB: Testing banded Singular Value Decomposition routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../Rbb.in; \
	} > $@ 2>&1
$(OUTDIR)/Rglm.out: ../glm.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../glm.in; \
	} > $@ 2>&1
$(OUTDIR)/Rgqr.out: ../gqr.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../gqr.in; \
	} > $@ 2>&1
$(OUTDIR)/Rgsv.out: ../gsv.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../gsv.in; \
	} > $@ 2>&1
$(OUTDIR)/Rlse.out: ../lse.in xeigtstR_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_gmp$(EXEEXT) < ../lse.in; \
	} > $@ 2>&1
$(OUTDIR)/Cnep.out: ../nep.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../nep.in; \
	} > $@ 2>&1
$(OUTDIR)/Csep.out: ../sep.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../sep.in; \
	} > $@ 2>&1
$(OUTDIR)/Cse2.out: ../se2.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../se2.in; \
	} > $@ 2>&1
$(OUTDIR)/Csvd.out: ../svd.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../svd.in; \
	} > $@ 2>&1
$(OUTDIR)/Cec.out: ../Cec.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CEC: Testing COMPLEX Eigen Condition Routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../Cec.in; \
	} > $@ 2>&1
$(OUTDIR)/Ced.out: ../Ced.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CES: Testing COMPLEX Nonsymmetric Schur Form Driver [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../Ced.in; \
	} > $@ 2>&1
$(OUTDIR)/Cgg.out: ../Cgg.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGG: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../Cgg.in; \
	} > $@ 2>&1
$(OUTDIR)/Cgd.out: ../Cgd.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem driver routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../Cgd.in; \
	} > $@ 2>&1
$(OUTDIR)/Csb.out: ../Csb.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CHB: Testing Hermitian Eigenvalue Problem routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../Csb.in; \
	} > $@ 2>&1
$(OUTDIR)/Csg.out: ../Csg.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../Csg.in; \
	} > $@ 2>&1
$(OUTDIR)/Cbal.out: ../Cbal.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../Cbal.in; \
	} > $@ 2>&1
$(OUTDIR)/Cbak.out: ../Cbak.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../Cbak.in; \
	} > $@ 2>&1
$(OUTDIR)/Cgbal.out: ../Cgbal.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../Cgbal.in; \
	} > $@ 2>&1
$(OUTDIR)/Cgbak.out: ../Cgbak.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../Cgbak.in; \
	} > $@ 2>&1
$(OUTDIR)/Cbb.out: ../Cbb.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../Cbb.in; \
	} > $@ 2>&1
$(OUTDIR)/Cglm.out: ../glm.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../glm.in; \
	} > $@ 2>&1
$(OUTDIR)/Cgqr.out: ../gqr.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../gqr.in; \
	} > $@ 2>&1
$(OUTDIR)/Cgsv.out: ../gsv.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../gsv.in; \
	} > $@ 2>&1
$(OUTDIR)/Clse.out: ../lse.in xeigtstC_gmp$(EXEEXT) | $(OUTDIR)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(GMP_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(GMP_LABEL) ==="; \
		$(GMP_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_gmp$(EXEEXT) < ../lse.in; \
	} > $@ 2>&1
