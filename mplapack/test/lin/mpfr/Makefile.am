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

MPFR_PROFILES = default binary64 binary128

MPFR_ENV_default   =
MPFR_ENV_binary64  = MPLAPACK_MPFR_EMAX=1024  MPLAPACK_MPFR_EMIN=-1021  MPLAPACK_MPFR_PRECISION=53
MPFR_ENV_binary128 = MPLAPACK_MPFR_EMAX=16384 MPLAPACK_MPFR_EMIN=-16381 MPLAPACK_MPFR_PRECISION=113

MPFR_PROFILE_LABEL_default   = default
MPFR_PROFILE_LABEL_binary64  = binary64 (prec=53 emin=-1021 emax=1024)
MPFR_PROFILE_LABEL_binary128 = binary128 (prec=113 emin=-16381 emax=16384)

# -----------------------------------------------------------------------------
# Outputs to build during "make check"
# -----------------------------------------------------------------------------

MPFR_LIN_OUTS = \
  $(addprefix Rtest.,$(addsuffix .out,$(MPFR_PROFILES))) \
  $(addprefix Ctest.,$(addsuffix .out,$(MPFR_PROFILES))) \
  $(addprefix Rtest_rfp.,$(addsuffix .out,$(MPFR_PROFILES))) \
  $(addprefix Ctest_rfp.,$(addsuffix .out,$(MPFR_PROFILES)))

.PHONY: mpfr-check mpfr-check-%

check-local: mpfr-check

mpfr-check: $(MPFR_LIN_OUTS)
	@echo "MPFR lin tests done for profiles: $(MPFR_PROFILES)"

mpfr-check-%: Rtest.%.out Ctest.%.out Rtest_rfp.%.out Ctest_rfp.%.out
	@echo "MPFR profile '$*' done."

# -----------------------------------------------------------------------------
# Test rules (double-style, profile is the stem)
# -----------------------------------------------------------------------------

Rtest.%.out: ../Rtest.in xlintstR_mpfr
	@echo "Testing Real MPLAPACK linear equation routines: $(MPFR_ENV_$*) ./xlintstR_mpfr < ../Rtest.in > $@ 2>&1"
	@set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) ./xlintstR_mpfr < ../Rtest.in; \
	} > $@ 2>&1

Ctest.%.out: ../Ctest.in xlintstC_mpfr
	@echo "Testing Complex MPLAPACK linear equation routines: $(MPFR_ENV_$*) ./xlintstC_mpfr < ../Ctest.in > $@ 2>&1"
	@set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) ./xlintstC_mpfr < ../Ctest.in; \
	} > $@ 2>&1

Rtest_rfp.%.out: ../Rtest_rfp.in xlintstrfR_mpfr
	@echo "Testing Real MPLAPACK RFP prototype linear equation routines: $(MPFR_ENV_$*) ./xlintstrfR_mpfr < ../Rtest_rfp.in > $@ 2>&1"
	@set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) ./xlintstrfR_mpfr < ../Rtest_rfp.in; \
	} > $@ 2>&1

Ctest_rfp.%.out: ../Ctest_rfp.in xlintstrfC_mpfr
	@echo "Testing Complex MPLAPACK RFP prototype linear equation routines: $(MPFR_ENV_$*) ./xlintstrfC_mpfr < ../Ctest_rfp.in > $@ 2>&1"
	@set -e; { \
		echo "=== MPFR profile: $(MPFR_PROFILE_LABEL_$*) ==="; \
		$(MPFR_ENV_$*) ./xlintstrfC_mpfr < ../Ctest_rfp.in; \
	} > $@ 2>&1

clean-local:
	rm -f *.out *~ \#*
