check_PROGRAMS = $(mplapack_td_test_PROGRAMS) $(mplapack_td_opt_test_PROGRAMS)

# Own the test routines and common arrays in the executable only.
# Do not also link an override library: its finalizers would run twice.
TEST_OVERRIDE_SOURCES = ../common/Mxerbla.cpp ../common/iMlaenv.cpp ../common/Mxlaenv.cpp

mplapack_td_test_PROGRAMS = xeigtstR_td xeigtstC_td xdmdeigtstR_td xdmdeigtstC_td
mplapack_td_opt_test_PROGRAMS = xeigtstR_td_opt xeigtstC_td_opt xdmdeigtstR_td_opt xdmdeigtstC_td_opt
mplapack_td_testdir = $(prefix)/lib/$(target)/mplapack/test/eig
mplapack_td_opt_testdir = $(prefix)/lib/$(target)/mplapack/test/eig

$(mplapack_td_test_PROGRAMS): ../libeig_td.la
$(mplapack_td_opt_test_PROGRAMS): ../libeig_td_opt.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_td.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_td.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_td.a,--no-whole-archive
opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_td_opt.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_td_opt.a,$(top_builddir)/mpblas/optimized/td/.libs/libmplapack_td_opt.a,--no-whole-archive
dmd_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_td.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_td.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_td.a,--no-whole-archive
dmd_opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_td_opt.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_td_opt.a,$(top_builddir)/mpblas/optimized/td/.libs/libmplapack_td_opt.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/eig -leig_td -L$(top_builddir)/mplapack/test/matgen -lmatgen_td -L$(top_builddir)/mplapack/reference -lmplapack_td
opt_libdepends = -L$(top_builddir)/mplapack/test/eig -leig_td_opt -L$(top_builddir)/mplapack/test/matgen -lmatgen_td_opt -L$(top_builddir)/mpblas/optimized/td -lmplapack_td_opt
dmd_libdepends = -L$(top_builddir)/mplapack/test/eig -leig_td -L$(top_builddir)/mplapack/test/matgen -lmatgen_td -L$(top_builddir)/mplapack/reference -lmplapack_td
dmd_opt_libdepends = -L$(top_builddir)/mplapack/test/eig -leig_td_opt -L$(top_builddir)/mplapack/test/matgen -lmatgen_td_opt -L$(top_builddir)/mpblas/optimized/td -lmplapack_td_opt
endif

mplibs = -L$(QD_LIBDIR) -lqd $(QD_RPATH_LDFLAGS) -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -I$(QD_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_TD @qd_extra_cxxflags@ -DMPLAPACK_INTERNAL -I$(top_srcdir)/fable

xeigtstR_td_SOURCES  = ../common/Rchkee.cpp $(TEST_OVERRIDE_SOURCES)
xeigtstR_td_opt_SOURCES  = ../common/Rchkee.cpp $(TEST_OVERRIDE_SOURCES)
xeigtstR_td_LDFLAGS  = $(libdepends) $(mplibs) #-L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_td_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) #-L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_td_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_td_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_td_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xeigtstR_td_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

xeigtstC_td_SOURCES  = ../common/Cchkee.cpp $(TEST_OVERRIDE_SOURCES)
xeigtstC_td_opt_SOURCES  = ../common/Cchkee.cpp $(TEST_OVERRIDE_SOURCES)
xeigtstC_td_LDFLAGS  = $(libdepends) $(mplibs) #-L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_td_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) #-L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_td_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_td_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_td_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xeigtstC_td_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

xdmdeigtstR_td_SOURCES  = ../common/Rchkdmd.cpp ../common/Mxerbla.cpp ../common/Mxlaenv.cpp
xdmdeigtstR_td_opt_SOURCES  = ../common/Rchkdmd.cpp ../common/Mxerbla.cpp ../common/Mxlaenv.cpp
xdmdeigtstR_td_LDFLAGS  = $(dmd_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstR_td_opt_LDFLAGS  = $(dmd_opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstR_td_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstR_td_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstR_td_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xdmdeigtstR_td_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

xdmdeigtstC_td_SOURCES  = ../common/Cchkdmd.cpp ../common/Mxerbla.cpp ../common/Mxlaenv.cpp
xdmdeigtstC_td_opt_SOURCES  = ../common/Cchkdmd.cpp ../common/Mxerbla.cpp ../common/Mxlaenv.cpp
xdmdeigtstC_td_LDFLAGS  = $(dmd_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstC_td_opt_LDFLAGS  = $(dmd_opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstC_td_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstC_td_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstC_td_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)
xdmdeigtstC_td_opt_CXXFLAGS = $(cxxflags) $(DD_CXXFLAGS)

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_td_testdir) $(prefix)
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

OUTDIR_PATH = ../results/$(OUTDIR)/td



# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <OUTDIR>/<test>.out)
# -----------------------------------------------------------------------------

TD_REAL_TESTS = \
	Rnep Rsep Rse2 Rsvd Rec Red Rgg Rgd Rsb Rsg Rbal Rbak Rgbal Rgbak Rbb Rglm Rgqr Rgsv Rcsd Rlse Rdmd

TD_COMPLEX_TESTS = \
	Cnep Csep Cse2 Csvd Cec Ced Cgg Cgd Csb Csg Cbal Cbak Cgbal Cgbak Cbb Cglm Cgqr Cgsv Ccsd Clse Cdmd

TD_ALL_TESTS = $(TD_REAL_TESTS) $(TD_COMPLEX_TESTS)

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
# -----------------------------------------------------------------------------

VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended before running the test (last assignment wins)
TD_EXTRA_ENV ?=
TD_LABEL ?= td

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
.PHONY: td-check FORCE
FORCE:

TD_EIG_OUTS = $(foreach t,$(TD_ALL_TESTS),$(OUTDIR_PATH)/$(t).out)
CLEANFILES = $(TD_EIG_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(TD_EIG_OUTS): FORCE

check-local: td-check

.PRECIOUS: $(TD_EIG_OUTS)

td-check: $(TD_EIG_OUTS)
	@echo "$(TD_LABEL) eig tests done: $(TD_ALL_TESTS)"
	@echo "Output directory: $(OUTDIR_PATH)"

clean-local:
	rm -f $(TD_EIG_OUTS) *~ \#*
	@if test -d "$(OUTDIR_PATH)"; then rmdir -p "$(OUTDIR_PATH)" 2>/dev/null || true; fi

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

$(OUTDIR_PATH)/Rnep.out: $(srcdir)/../nep.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../nep.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rsep.out: $(srcdir)/../sep.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../sep.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rse2.out: $(srcdir)/../se2.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../se2.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rsvd.out: $(srcdir)/../svd.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../svd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rec.out: $(srcdir)/../Rec.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../Rec.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Red.out: $(srcdir)/../Red.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../Red.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rgg.out: $(srcdir)/../Rgg.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../Rgg.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rgd.out: $(srcdir)/../Rgd.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../Rgd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rsb.out: $(srcdir)/../Rsb.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSB: Testing Symmetric Eigenvalue Problem routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../Rsb.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rsg.out: $(srcdir)/../Rsg.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../Rsg.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rbal.out: $(srcdir)/../Rbal.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebal: Testing the balancing of a general matrix [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../Rbal.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rbak.out: $(srcdir)/../Rbak.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebak: Testing the back transformation of balanced matrix [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../Rbak.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rgbal.out: $(srcdir)/../Rgbal.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbal: Testing the balancing of a pair of general matrices [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../Rgbal.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rgbak.out: $(srcdir)/../Rgbak.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbak: Testing the back transformation of a pair of balanced matrices [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../Rgbak.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rbb.out: $(srcdir)/../Rbb.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DBB: Testing banded Singular Value Decomposition routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../Rbb.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rglm.out: $(srcdir)/../glm.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../glm.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rgqr.out: $(srcdir)/../gqr.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../gqr.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rgsv.out: $(srcdir)/../gsv.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../gsv.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rcsd.out: $(srcdir)/../csd.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../csd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Rlse.out: $(srcdir)/../lse.in xeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_td$(EXEEXT) < $(srcdir)/../lse.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cnep.out: $(srcdir)/../nep.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../nep.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Csep.out: $(srcdir)/../sep.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../sep.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cse2.out: $(srcdir)/../se2.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../se2.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Csvd.out: $(srcdir)/../svd.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../svd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cec.out: $(srcdir)/../Cec.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CEC: Testing COMPLEX Eigen Condition Routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../Cec.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Ced.out: $(srcdir)/../Ced.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CES: Testing COMPLEX Nonsymmetric Schur Form Driver [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../Ced.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cgg.out: $(srcdir)/../Cgg.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGG: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../Cgg.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cgd.out: $(srcdir)/../Cgd.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem driver routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../Cgd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Csb.out: $(srcdir)/../Csb.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CHB: Testing Hermitian Eigenvalue Problem routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../Csb.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Csg.out: $(srcdir)/../Csg.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../Csg.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cbal.out: $(srcdir)/../Cbal.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../Cbal.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cbak.out: $(srcdir)/../Cbak.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../Cbak.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cgbal.out: $(srcdir)/../Cgbal.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../Cgbal.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cgbak.out: $(srcdir)/../Cgbak.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../Cgbak.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cbb.out: $(srcdir)/../Cbb.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../Cbb.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cglm.out: $(srcdir)/../glm.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../glm.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cgqr.out: $(srcdir)/../gqr.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../gqr.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cgsv.out: $(srcdir)/../gsv.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../gsv.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Ccsd.out: $(srcdir)/../csd.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../csd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Clse.out: $(srcdir)/../lse.in xeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_td$(EXEEXT) < $(srcdir)/../lse.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rdmd.out: $(srcdir)/../Rdmd.in xdmdeigtstR_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DMD: Testing Dynamic Mode Decomposition routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xdmdeigtstR_td$(EXEEXT) < $(srcdir)/../Rdmd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cdmd.out: $(srcdir)/../Cdmd.in xdmdeigtstC_td$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DMD: Testing COMPLEX Dynamic Mode Decomposition routines [$(TD_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(TD_LABEL) ==="; \
		$(TD_EXTRA_ENV) $$wrapper $(WINE) ./xdmdeigtstC_td$(EXEEXT) < $(srcdir)/../Cdmd.in; \
	} > $@ 2>&1
