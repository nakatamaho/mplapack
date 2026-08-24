check_PROGRAMS = $(mplapack_binary128_test_PROGRAMS) $(mplapack_binary128_opt_test_PROGRAMS)

mplapack_binary128_test_PROGRAMS = xeigtstR_binary128 xeigtstC_binary128 xdmdeigtstR_binary128 xdmdeigtstC_binary128
mplapack_binary128_opt_test_PROGRAMS = xeigtstR_binary128_opt xeigtstC_binary128_opt xdmdeigtstR_binary128_opt xdmdeigtstC_binary128_opt
mplapack_binary128_testdir = $(prefix)/lib/$(target)/mplapack/test/eig
mplapack_binary128_opt_testdir = $(prefix)/lib/$(target)/mplapack/test/eig

$(mplapack_binary128_test_PROGRAMS): ../libeig_binary128.la
$(mplapack_binary128_opt_test_PROGRAMS): ../libeig_binary128_opt.la
xeigtstR_binary128 xeigtstC_binary128: ../libeig_override_binary128.la
xeigtstR_binary128_opt xeigtstC_binary128_opt: ../libeig_override_binary128_opt.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_binary128.a,$(top_builddir)/mplapack/test/eig/.libs/libeig_override_binary128.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_binary128.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_binary128.a,--no-whole-archive
opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_binary128_opt.a,$(top_builddir)/mplapack/test/eig/.libs/libeig_override_binary128_opt.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_binary128_opt.a,$(top_builddir)/mpblas/optimized/binary128/.libs/libmplapack_binary128_opt.a,--no-whole-archive
dmd_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_binary128.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_binary128.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_binary128.a,--no-whole-archive
dmd_opt_libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/eig/.libs/libeig_binary128_opt.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_binary128_opt.a,$(top_builddir)/mpblas/optimized/binary128/.libs/libmplapack_binary128_opt.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/eig -leig_binary128 -leig_override_binary128 -L$(top_builddir)/mplapack/test/matgen -lmatgen_binary128 -L$(top_builddir)/mplapack/reference -lmplapack_binary128
opt_libdepends = -L$(top_builddir)/mplapack/test/eig -leig_binary128_opt -leig_override_binary128_opt -L$(top_builddir)/mplapack/test/matgen -lmatgen_binary128_opt -L$(top_builddir)/mpblas/optimized/binary128 -lmplapack_binary128_opt
dmd_libdepends = -L$(top_builddir)/mplapack/test/eig -leig_binary128 -L$(top_builddir)/mplapack/test/matgen -lmatgen_binary128 -L$(top_builddir)/mplapack/reference -lmplapack_binary128
dmd_opt_libdepends = -L$(top_builddir)/mplapack/test/eig -leig_binary128_opt -L$(top_builddir)/mplapack/test/matgen -lmatgen_binary128_opt -L$(top_builddir)/mpblas/optimized/binary128 -lmplapack_binary128_opt
endif

if BINARY128_USE_QUADMATH
mplibs = -lm -lquadmath
else
mplibs = -lm
endif
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_BINARY128 -DMPLAPACK_INTERNAL -I$(top_srcdir)/fable

xeigtstR_binary128_SOURCES  = ../common/Rchkee.cpp
xeigtstR_binary128_opt_SOURCES  = ../common/Rchkee.cpp
xeigtstR_binary128_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_binary128_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_binary128_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_binary128_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_binary128_CXXFLAGS = $(cxxflags)
xeigtstR_binary128_opt_CXXFLAGS = $(cxxflags)

xeigtstC_binary128_SOURCES  = ../common/Cchkee.cpp
xeigtstC_binary128_opt_SOURCES  = ../common/Cchkee.cpp
xeigtstC_binary128_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_binary128_opt_LDFLAGS  = $(opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_binary128_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_binary128_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_binary128_CXXFLAGS = $(cxxflags)
xeigtstC_binary128_opt_CXXFLAGS = $(cxxflags)

xdmdeigtstR_binary128_SOURCES  = ../common/Rchkdmd.cpp ../common/Mxerbla.cpp ../common/Mxlaenv.cpp
xdmdeigtstR_binary128_opt_SOURCES  = ../common/Rchkdmd.cpp ../common/Mxerbla.cpp ../common/Mxlaenv.cpp
xdmdeigtstR_binary128_LDFLAGS  = $(dmd_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstR_binary128_opt_LDFLAGS  = $(dmd_opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstR_binary128_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstR_binary128_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstR_binary128_CXXFLAGS = $(cxxflags)
xdmdeigtstR_binary128_opt_CXXFLAGS = $(cxxflags)

xdmdeigtstC_binary128_SOURCES  = ../common/Cchkdmd.cpp ../common/Mxerbla.cpp ../common/Mxlaenv.cpp
xdmdeigtstC_binary128_opt_SOURCES  = ../common/Cchkdmd.cpp ../common/Mxerbla.cpp ../common/Mxlaenv.cpp
xdmdeigtstC_binary128_LDFLAGS  = $(dmd_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstC_binary128_opt_LDFLAGS  = $(dmd_opt_libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xdmdeigtstC_binary128_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstC_binary128_opt_LDFLAGS += -Wl,-rpath,$(libdir)
xdmdeigtstC_binary128_CXXFLAGS = $(cxxflags)
xdmdeigtstC_binary128_opt_CXXFLAGS = $(cxxflags)

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapack_binary128_testdir) $(prefix)
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

OUTDIR_PATH = ../results/$(OUTDIR)/binary128



# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <OUTDIR>/<test>.out)
# -----------------------------------------------------------------------------

FLOAT128_REAL_TESTS = \
	Rnep Rsep Rse2 Rsvd Rec Red Rgg Rgd Rsb Rsg Rbal Rbak Rgbal Rgbak Rbb Rglm Rgqr Rgsv Rcsd Rlse Rdmd

FLOAT128_COMPLEX_TESTS = \
	Cnep Csep Cse2 Csvd Cec Ced Cgg Cgd Csb Csg Cbal Cbak Cgbal Cgbak Cbb Cglm Cgqr Cgsv Ccsd Clse Cdmd

FLOAT128_ALL_TESTS = $(FLOAT128_REAL_TESTS) $(FLOAT128_COMPLEX_TESTS)

# -----------------------------------------------------------------------------
# Optional valgrind wrapper (toggle from command line)
# -----------------------------------------------------------------------------

VALGRIND ?=
VALGRIND_FLAGS ?= --show-leak-kinds=all --track-origins=yes --leak-check=full --error-exitcode=99

# Extra env appended before running the test (last assignment wins)
FLOAT128_EXTRA_ENV ?=
FLOAT128_LABEL ?= binary128

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
.PHONY: float128-check FORCE
FORCE:

FLOAT128_EIG_OUTS = $(foreach t,$(FLOAT128_ALL_TESTS),$(OUTDIR_PATH)/$(t).out)
CLEANFILES = $(FLOAT128_EIG_OUTS)

# Re-run all test commands on every "make check" even if *.out already exists.
$(FLOAT128_EIG_OUTS): FORCE

check-local: float128-check

.PRECIOUS: $(FLOAT128_EIG_OUTS)

float128-check: $(FLOAT128_EIG_OUTS)
	@echo "binary128 eig tests done: $(FLOAT128_ALL_TESTS)"
	@echo "Output directory: $(OUTDIR_PATH)"

clean-local:
	rm -f $(FLOAT128_EIG_OUTS) *~ \#*
	@if test -d "$(OUTDIR_PATH)"; then rmdir -p "$(OUTDIR_PATH)" 2>/dev/null || true; fi

# -----------------------------------------------------------------------------
# Test rules (one output per test)
# -----------------------------------------------------------------------------

$(OUTDIR_PATH)/Rnep.out: $(srcdir)/../nep.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../nep.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rsep.out: $(srcdir)/../sep.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../sep.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rse2.out: $(srcdir)/../se2.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../se2.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rsvd.out: $(srcdir)/../svd.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../svd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rec.out: $(srcdir)/../Rec.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEC: Testing Eigen Condition Routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../Rec.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Red.out: $(srcdir)/../Red.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DEV: Testing Nonsymmetric Eigenvalue Driver [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../Red.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgg.out: $(srcdir)/../Rgg.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../Rgg.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgd.out: $(srcdir)/../Rgd.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../Rgd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rsb.out: $(srcdir)/../Rsb.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSB: Testing Symmetric Eigenvalue Problem routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../Rsb.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rsg.out: $(srcdir)/../Rsg.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../Rsg.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rbal.out: $(srcdir)/../Rbal.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebal: Testing the balancing of a general matrix [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../Rbal.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rbak.out: $(srcdir)/../Rbak.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rgebak: Testing the back transformation of balanced matrix [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../Rbak.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgbal.out: $(srcdir)/../Rgbal.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbal: Testing the balancing of a pair of general matrices [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../Rgbal.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgbak.out: $(srcdir)/../Rgbak.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "Rggbak: Testing the back transformation of a pair of balanced matrices [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../Rgbak.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rbb.out: $(srcdir)/../Rbb.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DBB: Testing banded Singular Value Decomposition routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../Rbb.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rglm.out: $(srcdir)/../glm.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../glm.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgqr.out: $(srcdir)/../gqr.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../gqr.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rgsv.out: $(srcdir)/../gsv.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../gsv.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rcsd.out: $(srcdir)/../csd.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../csd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rlse.out: $(srcdir)/../lse.in xeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstR_binary128$(EXEEXT) < $(srcdir)/../lse.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cnep.out: $(srcdir)/../nep.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "NEP: Testing Nonsymmetric Eigenvalue Problem routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../nep.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Csep.out: $(srcdir)/../sep.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../sep.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cse2.out: $(srcdir)/../se2.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SEP: Testing Symmetric Eigenvalue Problem routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../se2.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Csvd.out: $(srcdir)/../svd.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "SVD: Testing Singular Value Decomposition routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../svd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cec.out: $(srcdir)/../Cec.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CEC: Testing COMPLEX Eigen Condition Routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../Cec.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ced.out: $(srcdir)/../Ced.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CES: Testing COMPLEX Nonsymmetric Schur Form Driver [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../Ced.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgg.out: $(srcdir)/../Cgg.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGG: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../Cgg.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgd.out: $(srcdir)/../Cgd.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGD: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem driver routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../Cgd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Csb.out: $(srcdir)/../Csb.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CHB: Testing Hermitian Eigenvalue Problem routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../Csb.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Csg.out: $(srcdir)/../Csg.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSG: Testing Symmetric Generalized Eigenvalue Problem routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../Csg.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cbal.out: $(srcdir)/../Cbal.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAL: Testing the balancing of a COMPLEX general matrix [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../Cbal.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cbak.out: $(srcdir)/../Cbak.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGEBAK: Testing the back transformation of a COMPLEX balanced matrix [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../Cbak.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgbal.out: $(srcdir)/../Cgbal.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAL: Testing the balancing of a pair of COMPLEX general matrices [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../Cgbal.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgbak.out: $(srcdir)/../Cgbak.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../Cgbak.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cbb.out: $(srcdir)/../Cbb.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CBB: Testing banded Singular Value Decomposition routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../Cbb.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cglm.out: $(srcdir)/../glm.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GLM: Testing Generalized Linear Regression Model routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../glm.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgqr.out: $(srcdir)/../gqr.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GQR: Testing Generalized QR and RQ factorization routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../gqr.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Cgsv.out: $(srcdir)/../gsv.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "GSV: Testing Generalized Singular Value Decomposition routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../gsv.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Ccsd.out: $(srcdir)/../csd.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "CSD: Testing CS Decomposition routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../csd.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Clse.out: $(srcdir)/../lse.in xeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "LSE: Testing Constrained Linear Least Squares routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xeigtstC_binary128$(EXEEXT) < $(srcdir)/../lse.in; \
	} > $@ 2>&1

$(OUTDIR_PATH)/Rdmd.out: $(srcdir)/../Rdmd.in xdmdeigtstR_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DMD: Testing Dynamic Mode Decomposition routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xdmdeigtstR_binary128$(EXEEXT) < $(srcdir)/../Rdmd.in; \
	} > $@ 2>&1
$(OUTDIR_PATH)/Cdmd.out: $(srcdir)/../Cdmd.in xdmdeigtstC_binary128$(EXEEXT)
	@$(MKDIR_P) $(OUTDIR_PATH)
	@wrapper=""; \
	if test -n "$(VALGRIND)"; then \
		wrapper="$(LIBTOOL) --mode=execute valgrind $(VALGRIND_FLAGS)"; \
	fi; \
	echo "DMD: Testing COMPLEX Dynamic Mode Decomposition routines [$(FLOAT128_LABEL)]"; \
	set -e; { \
		echo "=== precision: $(FLOAT128_LABEL) ==="; \
		$(FLOAT128_EXTRA_ENV) $$wrapper $(WINE) ./xdmdeigtstC_binary128$(EXEEXT) < $(srcdir)/../Cdmd.in; \
	} > $@ 2>&1
