check_PROGRAMS = $(mplapack_mpfr_test_PROGRAMS)

mplapack_mpfr_test_PROGRAMS = xlintstR_mpfr xlintstC_mpfr xlintstrfR_mpfr xlintstrfC_mpfr
mplapack_mpfr_testdir = $(prefix)/lib/$(target)/mplapack/test/lin

$(check_PROGRAMS): ../liblin_mpfr.la ../liblin_override_mpfr.la

if IS_MINGW
libdepends = -Wl,--allow-multiple-definition -Wl,--whole-archive,$(top_builddir)/mplapack/test/lin/.libs/liblin_mpfr.a,$(top_builddir)/mplapack/test/lin/.libs/liblin_override_mpfr.a,$(top_builddir)/mplapack/test/matgen/.libs/libmatgen_mpfr.a,$(top_builddir)/mpblas/reference/.libs/libmpblas_mpfr.a,$(top_builddir)/mplapack/reference/.libs/libmplapack_mpfr.a,--no-whole-archive
else
libdepends = -L$(top_builddir)/mplapack/test/lin -llin_mpfr -llin_override_mpfr -L$(top_builddir)/mplapack/test/matgen -lmatgen_mpfr -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr
endif

mplibs = -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpc -lmpfr -lgmpxx -lgmp -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/fable -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -D___MPLAPACK_INTERNAL___

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapack_mpfr_testdir) $(prefix)
endif

xlintstR_mpfr_SOURCES  = ../common/Rchkaa.cpp
xlintstR_mpfr_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstR_mpfr_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstR_mpfr_CXXFLAGS = $(cxxflags) 

xlintstC_mpfr_SOURCES  = ../common/Cchkaa.cpp
xlintstC_mpfr_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstC_mpfr_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstC_mpfr_CXXFLAGS = $(cxxflags)

xlintstrfR_mpfr_SOURCES  = ../common/Rchkrfp.cpp
xlintstrfR_mpfr_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfR_mpfr_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfR_mpfr_CXXFLAGS = $(cxxflags)

xlintstrfC_mpfr_SOURCES  = ../common/Cchkrfp.cpp
xlintstrfC_mpfr_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xlintstrfC_mpfr_LDFLAGS += -Wl,-rpath,$(libdir)
xlintstrfC_mpfr_CXXFLAGS = $(cxxflags)

# -----------------------------------------------------------------------------
# MPFR emulation profiles
# -----------------------------------------------------------------------------

MPFR_PROFILES = current binary64 binary128

MPFR_ENV_current  =
MPFR_ENV_binary64  = MPLAPACK_MPFR_EMAX=1024  MPLAPACK_MPFR_EMIN=-1021  MPLAPACK_MPFR_PRECISION=53
MPFR_ENV_binary128 = MPLAPACK_MPFR_EMAX=16384 MPLAPACK_MPFR_EMIN=-16381 MPLAPACK_MPFR_PRECISION=113

MPFR_PROFILE_LABEL_current  = current
MPFR_PROFILE_LABEL_binary64  = binary64 (prec=53 emin=-1021 emax=1024)
MPFR_PROFILE_LABEL_binary128 = binary128 (prec=113 emin=-16381 emax=16384)

# -----------------------------------------------------------------------------
# Test sets (base names; outputs are <test>.<profile>.out)
# -----------------------------------------------------------------------------

MPFR_LIN_REAL_TESTS = Rtest Rtest_rfp
MPFR_LIN_COMPLEX_TESTS = Ctest Ctest_rfp
MPFR_LIN_ALL_TESTS = $(MPFR_LIN_REAL_TESTS) $(MPFR_LIN_COMPLEX_TESTS)

# -----------------------------------------------------------------------------
# Hook into Automake's "make check"
# -----------------------------------------------------------------------------

.PHONY: mpfr-check mpfr-check-%

# DO NOT define "check:" here. Use check-local to avoid target redefinition merges.
check-local: mpfr-check

mpfr-check: $(MPFR_PROFILES:%=mpfr-check-%)
	@echo "MPFR lin tests done for profiles: $(MPFR_PROFILES)"

# Pattern target: builds all <test>.<profile>.out for the given profile.
mpfr-check-%: $(addsuffix .%.out,$(MPFR_LIN_ALL_TESTS))
	@echo "MPFR profile '$*' done."

clean-local:
	rm -f *.out *~ \#*

# -----------------------------------------------------------------------------
# Rule generator (GNU make extensions)
# -----------------------------------------------------------------------------

define MPFR_LIN_RULE
$1.$2.out: $3 $4
	@echo "$5 [$(MPFR_PROFILE_LABEL_$2)]"
	@set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$2) ==="; \
		$$(MPFR_ENV_$2) ./$4 < $3; \
	} > $$@ 2>&1
endef

define GEN_MPFR_LIN_TEST
$$(foreach p,$$(MPFR_PROFILES),$$(eval $$(call MPFR_LIN_RULE,$1,$$(p),$2,$3,$4)))
endef

# -----------------------------------------------------------------------------
# Linear-equation tests (lin/)
# -----------------------------------------------------------------------------

$(eval $(call GEN_MPFR_LIN_TEST,Rtest,../Rtest.in,xlintstR_mpfr,Testing Real MPLAPACK linear equation routines))
$(eval $(call GEN_MPFR_LIN_TEST,Ctest,../Ctest.in,xlintstC_mpfr,Testing Complex MPLAPACK linear equation routines))
$(eval $(call GEN_MPFR_LIN_TEST,Rtest_rfp,../Rtest_rfp.in,xlintstrfR_mpfr,Testing Real MPLAPACK RFP prototype linear equation routines))
$(eval $(call GEN_MPFR_LIN_TEST,Ctest_rfp,../Ctest_rfp.in,xlintstrfC_mpfr,Testing Complex MPLAPACK RFP prototype linear equation routines))
