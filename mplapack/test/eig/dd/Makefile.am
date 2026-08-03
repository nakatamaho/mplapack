check_PROGRAMS = $(mplapack_dd_test_PROGRAMS) $(mplapack_dd_opt_test_PROGRAMS)

mplapack_dd_test_PROGRAMS = xeigtstR_dd xeigtstC_dd xdmdeigtstR_dd xdmdeigtstC_dd
mplapack_dd_opt_test_PROGRAMS = xeigtstR_dd_opt xeigtstC_dd_opt xdmdeigtstR_dd_opt xdmdeigtstC_dd_opt
mplapack_dd_testdir = $(prefix)/lib/$(target)/mplapack/test/eig
mplapack_dd_opt_testdir = $(prefix)/lib/$(target)/mplapack/test/eig

$(mplapack_dd_test_PROGRAMS): ../libeig_dd.la ../libeig_override_dd.la
$(mplapack_dd_opt_test_PROGRAMS): ../libeig_dd_opt.la ../libeig_override_dd_opt.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_dd.a,$(top_builddir)/mplapack/test/eig/.libs/libeig_override_dd.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_dd.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_dd.a,--no-whole-archive
opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_dd_opt.a,$(top_builddir)/mplapack/test/eig/.libs/libeig_override_dd_opt.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_dd_opt.a,$(top_builddir)/mpblas/optimized/dd/.libs/libmplapack_dd_opt.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/eig -leig_dd -leig_override_dd -L$(top_builddir)/mplapack/test/matgen -lmatgen_dd -L$(top_builddir)/mplapack/reference -lmplapack_dd
opt_libdepends = -L$(top_builddir)/mplapack/test/eig -leig_dd_opt -leig_override_dd_opt -L$(top_builddir)/mplapack/test/matgen -lmatgen_dd_opt -L$(top_builddir)/mpblas/optimized/dd -lmplapack_dd_opt
endif

mplibs = -L$(QD_LIBDIR) -lqd $(QD_RPATH_LDFLAGS) -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -I$(QD_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_DD @qd_extra_cxxflags@ -DMPLAPACK_INTERNAL -I$(top_srcdir)/fable

xeigtstR_dd_SOURCES  = ../common/Rchkee.cpp
xeigtstR_dd_opt_SOURCES  = ../common/Rchkee.cpp
xeigtstR_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_dd_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_dd_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_dd_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_dd_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xeigtstR_dd_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

xeigtstC_dd_SOURCES  = ../common/Cchkee.cpp
xeigtstC_dd_opt_SOURCES  = ../common/Cchkee.cpp
xeigtstC_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_dd_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_dd_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_dd_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_dd_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xeigtstC_dd_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

xdmdeigtstR_dd_SOURCES  = ../common/Rchkdmd.cpp
xdmdeigtstR_dd_opt_SOURCES  = ../common/Rchkdmd.cpp
xdmdeigtstR_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstR_dd_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstR_dd_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstR_dd_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstR_dd_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xdmdeigtstR_dd_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

xdmdeigtstC_dd_SOURCES  = ../common/Cchkdmd.cpp
xdmdeigtstC_dd_opt_SOURCES  = ../common/Cchkdmd.cpp
xdmdeigtstC_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstC_dd_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstC_dd_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstC_dd_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstC_dd_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xdmdeigtstC_dd_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_dd_testdir) $(prefix)
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

OUTDIR_PATH = ../results/$(OUTDIR)/dd



# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <OUTDIR>/<test>.out)
# -----------------------------------------------------------------------------

DD_REAL_TESTS = \
	Rnep Rsep Rse2 Rsvd Rec Red Rgg Rgd Rsb Rsg Rbal Rbak Rgbal Rgbak Rbb Rglm Rgqr Rgsv Rcsd Rlse Rdmd

DD_COMPLEX_TESTS = \
	Cnep Csep Cse2 Csvd Cec Ced Cgg Cgd Csb Csg Cbal Cbak Cgbal Cgbak Cbb Cglm Cgqr Cgsv Ccsd Clse Cdmd

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

# All output logs: <OUTDIR>/<test>.out
DD_EIG_OUTS = $(foreach t,$(DD_ALL_TESTS),$(OUTDIR_PATH)/$(t).out)
CLEANFILES = $(DD_EIG_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(DD_EIG_OUTS): FORCE

check-local: dd-check

.PRECIOUS: $(DD_EIG_OUTS)

dd-check: $$(DD_EIG_OUTS)
	@echo "DD eig tests done."
	@echo "Output directory: $(OUTDIR_PATH)"

clean-local:
	rm -f $(DD_EIG_OUTS) *~ \#*
	@if test -d "$(OUTDIR_PATH)"; then rmdir -p "$(OUTDIR_PATH)" 2>/dev/null || true; fi

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

$(OUTDIR_PATH)/Rnep.out: $(srcdir)/../nep.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../nep.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rsep.out: $(srcdir)/../sep.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../sep.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rse2.out: $(srcdir)/../se2.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../se2.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rsvd.out: $(srcdir)/../svd.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../svd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rec.out: $(srcdir)/../Rec.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../Rec.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Red.out: $(srcdir)/../Red.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../Red.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgg.out: $(srcdir)/../Rgg.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGS: Testing Generalized Nonsymmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../Rgg.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgd.out: $(srcdir)/../Rgd.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routine"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../Rgd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rsb.out: $(srcdir)/../Rsb.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SBB: Testing banded Schur decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../Rsb.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rsg.out: $(srcdir)/../Rsg.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SGG: Testing Generalized Symmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../Rsg.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rbal.out: $(srcdir)/../Rbal.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGEBAL: Testing the balancing of a real general matrix"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../Rbal.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rbak.out: $(srcdir)/../Rbak.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGEBAK: Testing the back transformation of a real balanced matrix"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../Rbak.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgbal.out: $(srcdir)/../Rgbal.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGGBAL: Testing the balancing of a pair of real general matrices"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../Rgbal.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgbak.out: $(srcdir)/../Rgbak.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGGBAK: Testing the back transformation of a pair of real balanced matrices"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../Rgbak.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rbb.out: $(srcdir)/../Rbb.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DBB: Testing banded Singular Value Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../Rbb.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rglm.out: $(srcdir)/../glm.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../glm.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgqr.out: $(srcdir)/../gqr.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../gqr.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgsv.out: $(srcdir)/../gsv.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../gsv.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rcsd.out: $(srcdir)/../csd.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../csd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rlse.out: $(srcdir)/../lse.in xeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_dd$(EXEEXT) < $(srcdir)/../lse.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cnep.out: $(srcdir)/../nep.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../nep.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Csep.out: $(srcdir)/../sep.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../sep.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cse2.out: $(srcdir)/../se2.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../se2.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Csvd.out: $(srcdir)/../svd.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../svd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cec.out: $(srcdir)/../Cec.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../Cec.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ced.out: $(srcdir)/../Ced.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../Ced.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgg.out: $(srcdir)/../Cgg.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GGS: Testing Generalized Nonsymmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../Cgg.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgd.out: $(srcdir)/../Cgd.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routine"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../Cgd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Csb.out: $(srcdir)/../Csb.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SBB: Testing banded Schur decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../Csb.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Csg.out: $(srcdir)/../Csg.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SGG: Testing Generalized Symmetric Eigenvalue Problem routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../Csg.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cbal.out: $(srcdir)/../Cbal.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../Cbal.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cbak.out: $(srcdir)/../Cbak.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../Cbak.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgbal.out: $(srcdir)/../Cgbal.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../Cgbal.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgbak.out: $(srcdir)/../Cgbak.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../Cgbak.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cbb.out: $(srcdir)/../Cbb.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../Cbb.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cglm.out: $(srcdir)/../glm.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../glm.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgqr.out: $(srcdir)/../gqr.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../gqr.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgsv.out: $(srcdir)/../gsv.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../gsv.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ccsd.out: $(srcdir)/../csd.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../csd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Clse.out: $(srcdir)/../lse.in xeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_dd$(EXEEXT) < $(srcdir)/../lse.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rdmd.out: $(srcdir)/../Rdmd.in xdmdeigtstR_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DMD: Testing Dynamic Mode Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xdmdeigtstR_dd$(EXEEXT) < $(srcdir)/../Rdmd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cdmd.out: $(srcdir)/../Cdmd.in xdmdeigtstC_dd$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DMD: Testing COMPLEX Dynamic Mode Decomposition routines"; \
	set -e; { \
		$(DD_EXTRA_ENV) $$wrapper $(WINE) ./xdmdeigtstC_dd$(EXEEXT) < $(srcdir)/../Cdmd.in; \
	} > $@ 2>&1
