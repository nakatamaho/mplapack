check_PROGRAMS = $(mplapack_double_test_PROGRAMS) $(mplapack_double_opt_test_PROGRAMS)

mplapack_double_test_PROGRAMS = xeigtstR_double xeigtstC_double xdmdeigtstR_double xdmdeigtstC_double
mplapack_double_opt_test_PROGRAMS = xeigtstR_double_opt xeigtstC_double_opt xdmdeigtstR_double_opt xdmdeigtstC_double_opt
mplapack_double_testdir = $(prefix)/lib/$(target)/mplapack/test/eig
mplapack_double_opt_testdir = $(prefix)/lib/$(target)/mplapack/test/eig

$(mplapack_double_test_PROGRAMS): ../libeig_double.la ../libeig_override_double.la
$(mplapack_double_opt_test_PROGRAMS): ../libeig_double_opt.la ../libeig_override_double_opt.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_double.a,$(top_builddir)/mplapack/test/eig/.libs/libeig_override_double.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_double.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_double.a,--no-whole-archive
opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_double_opt.a,$(top_builddir)/mplapack/test/eig/.libs/libeig_override_double_opt.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_double_opt.a,$(top_builddir)/mpblas/optimized/double/.libs/libmplapack_double_opt.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/eig -leig_double -leig_override_double -L$(top_builddir)/mplapack/test/matgen -lmatgen_double -L$(top_builddir)/mplapack/reference -lmplapack_double
opt_libdepends = -L$(top_builddir)/mplapack/test/eig -leig_double_opt -leig_override_double_opt -L$(top_builddir)/mplapack/test/matgen -lmatgen_double_opt -L$(top_builddir)/mpblas/optimized/double -lmplapack_double_opt
endif
mplibs = -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_DOUBLE -DMPLAPACK_INTERNAL -I$(top_srcdir)/fable

xeigtstR_double_SOURCES  = ../common/Rchkee.cpp
xeigtstR_double_opt_SOURCES  = ../common/Rchkee.cpp
xeigtstR_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_double_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_double_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_double_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_double_CXXFLAGS = $(cxxflags)
xeigtstR_double_opt_CXXFLAGS = $(cxxflags)

xeigtstC_double_SOURCES  = ../common/Cchkee.cpp
xeigtstC_double_opt_SOURCES  = ../common/Cchkee.cpp
xeigtstC_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_double_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_double_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_double_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_double_CXXFLAGS = $(cxxflags)
xeigtstC_double_opt_CXXFLAGS = $(cxxflags)

xdmdeigtstR_double_SOURCES  = ../common/Rchkdmd.cpp
xdmdeigtstR_double_opt_SOURCES  = ../common/Rchkdmd.cpp
xdmdeigtstR_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstR_double_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstR_double_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstR_double_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstR_double_CXXFLAGS = $(cxxflags)
xdmdeigtstR_double_opt_CXXFLAGS = $(cxxflags)

xdmdeigtstC_double_SOURCES  = ../common/Cchkdmd.cpp
xdmdeigtstC_double_opt_SOURCES  = ../common/Cchkdmd.cpp
xdmdeigtstC_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstC_double_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstC_double_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstC_double_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstC_double_CXXFLAGS = $(cxxflags)
xdmdeigtstC_double_opt_CXXFLAGS = $(cxxflags)


install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_double_testdir) $(prefix)
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

_CPU_OS_RAW := $(shell bash $(abs_top_srcdir)/misc/cpu_os_detection.sh)
_CPU_TAG    := $(word 1,$(subst |, ,$(_CPU_OS_RAW)))
_OS_TAG     := $(word 2,$(subst |, ,$(_CPU_OS_RAW)))

if IS_MINGW
OUTDIR = $(_CPU_TAG)_windows_$(MPLAPACK_CC_TAG)
else
OUTDIR = $(_CPU_TAG)_$(_OS_TAG)_$(MPLAPACK_CC_TAG)
endif

OUTDIR_PATH = ../results/$(OUTDIR)/double



# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <OUTDIR>/<test>.out)
# -----------------------------------------------------------------------------

DOUBLE_REAL_TESTS = \
	Rnep Rsep Rse2 Rsvd Rec Red Rgg Rgd Rsb Rsg Rbal Rbak Rgbal Rgbak Rbb Rglm Rgqr Rgsv Rcsd Rlse Rdmd

DOUBLE_COMPLEX_TESTS = \
	Cnep Csep Cse2 Csvd Cec Ced Cgg Cgd Csb Csg Cbal Cbak Cgbal Cgbak Cbb Cglm Cgqr Cgsv Ccsd Clse Cdmd

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

DOUBLE_EIG_OUTS = $(foreach t,$(DOUBLE_ALL_TESTS),$(OUTDIR_PATH)/$(t).out)
CLEANFILES = $(DOUBLE_EIG_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(DOUBLE_EIG_OUTS): FORCE

check-local: double-check

.PRECIOUS: $(DOUBLE_EIG_OUTS)

double-check: $(DOUBLE_EIG_OUTS)
	@echo "$(DOUBLE_LABEL) eig tests done: $(DOUBLE_ALL_TESTS)"
	@echo "Output directory: $(OUTDIR_PATH)"

clean-local:
	rm -f $(DOUBLE_EIG_OUTS) *~ \#*
	@if test -d "$(OUTDIR_PATH)"; then rmdir -p "$(OUTDIR_PATH)" 2>/dev/null || true; fi

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

$(OUTDIR_PATH)/Rnep.out: $(srcdir)/../nep.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../nep.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rsep.out: $(srcdir)/../sep.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../sep.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rse2.out: $(srcdir)/../se2.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../se2.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rsvd.out: $(srcdir)/../svd.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../svd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rec.out: $(srcdir)/../Rec.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../Rec.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Red.out: $(srcdir)/../Red.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../Red.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rgg.out: $(srcdir)/../Rgg.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../Rgg.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rgd.out: $(srcdir)/../Rgd.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../Rgd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rsb.out: $(srcdir)/../Rsb.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSB: Testing Symmetric Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../Rsb.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rsg.out: $(srcdir)/../Rsg.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../Rsg.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rbal.out: $(srcdir)/../Rbal.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebal: Testing the balancing of a general matrix [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../Rbal.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rbak.out: $(srcdir)/../Rbak.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebak: Testing the back transformation of balanced matrix [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../Rbak.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rgbal.out: $(srcdir)/../Rgbal.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbal: Testing the balancing of a pair of general matrices [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../Rgbal.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rgbak.out: $(srcdir)/../Rgbak.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbak: Testing the back transformation of a pair of balanced matrices [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../Rgbak.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rbb.out: $(srcdir)/../Rbb.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DBB: Testing banded Singular Value Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../Rbb.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rglm.out: $(srcdir)/../glm.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../glm.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rgqr.out: $(srcdir)/../gqr.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../gqr.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rgsv.out: $(srcdir)/../gsv.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../gsv.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rcsd.out: $(srcdir)/../csd.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../csd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rlse.out: $(srcdir)/../lse.in xeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_double$(EXEEXT) < $(srcdir)/../lse.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cnep.out: $(srcdir)/../nep.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../nep.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Csep.out: $(srcdir)/../sep.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../sep.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cse2.out: $(srcdir)/../se2.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../se2.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Csvd.out: $(srcdir)/../svd.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../svd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cec.out: $(srcdir)/../Cec.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CEC: Testing COMPLEX Eigen Condition Routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../Cec.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Ced.out: $(srcdir)/../Ced.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CES: Testing COMPLEX Nonsymmetric Schur Form Driver [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../Ced.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cgg.out: $(srcdir)/../Cgg.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGG: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../Cgg.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cgd.out: $(srcdir)/../Cgd.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem driver routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../Cgd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Csb.out: $(srcdir)/../Csb.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CHB: Testing Hermitian Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../Csb.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Csg.out: $(srcdir)/../Csg.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../Csg.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cbal.out: $(srcdir)/../Cbal.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../Cbal.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cbak.out: $(srcdir)/../Cbak.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../Cbak.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cgbal.out: $(srcdir)/../Cgbal.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../Cgbal.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cgbak.out: $(srcdir)/../Cgbak.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../Cgbak.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cbb.out: $(srcdir)/../Cbb.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../Cbb.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cglm.out: $(srcdir)/../glm.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../glm.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cgqr.out: $(srcdir)/../gqr.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../gqr.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cgsv.out: $(srcdir)/../gsv.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../gsv.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Ccsd.out: $(srcdir)/../csd.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../csd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Clse.out: $(srcdir)/../lse.in xeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_double$(EXEEXT) < $(srcdir)/../lse.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rdmd.out: $(srcdir)/../Rdmd.in xdmdeigtstR_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DMD: Testing Dynamic Mode Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xdmdeigtstR_double$(EXEEXT) < $(srcdir)/../Rdmd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cdmd.out: $(srcdir)/../Cdmd.in xdmdeigtstC_double$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DMD: Testing COMPLEX Dynamic Mode Decomposition routines [$(DOUBLE_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(DOUBLE_LABEL) ==="; \
		$(DOUBLE_EXTRA_ENV) $$wrapper $(WINE) ./xdmdeigtstC_double$(EXEEXT) < $(srcdir)/../Cdmd.in; \
	} > $@ 2>&1
