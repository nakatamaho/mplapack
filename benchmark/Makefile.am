bench_PROGRAMS =

# -----------------------------------------------------------------------------
# Static data: gnuplot templates (.plt.in) installed alongside binaries.
# go.*.sh are generated at build time from .sh.in; not listed here.
# -----------------------------------------------------------------------------
bench_DATA = \
    Raxpy1.plt.in  Raxpy2.plt.in  Raxpy3.plt.in  \
    Rcopy1.plt.in  Rcopy2.plt.in  Rcopy3.plt.in  \
    Rdot1.plt.in   Rdot2.plt.in   Rdot3.plt.in   \
    Rgemm1.plt.in  Rgemm2.plt.in  Rgemm3.plt.in  \
    Rgemm_cuda.plt.in              \
    Rgemv1.plt.in  Rgemv2.plt.in  Rgemv3.plt.in  \
    Rsyrk1.plt.in  Rsyrk2.plt.in  Rsyrk3.plt.in  \
    Rsyrk_cuda.plt.in              \
    Rgetrf1.plt.in Rgetrf2.plt.in Rgetrf3.plt.in \
    Rgetrf_cuda.plt.in             \
    Rpotrf1.plt.in Rpotrf2.plt.in Rpotrf3.plt.in \
    Rpotrf_cuda.plt.in

# -----------------------------------------------------------------------------
# Source templates for go.*.sh (included in "make dist").
# -----------------------------------------------------------------------------
EXTRA_DIST = \
    go.Raxpy.sh.in       \
    go.Rcopy.sh.in       \
    go.Rdot.sh.in        \
    go.Rgemm.sh.in       \
    go.Rgemm_cuda.sh.in  \
    go.Rgemv.sh.in       \
    go.Rsyrk.sh.in       \
    go.Rsyrk_cuda.sh.in  \
    go.Rgetrf.sh.in      \
    go.Rgetrf_cuda.sh.in \
    go.Rpotrf.sh.in      \
    go.Rpotrf_cuda.sh.in

# -----------------------------------------------------------------------------
# Generated scripts (built in builddir, cleaned on "make clean").
# %%CPU_OS_DETECTION_SH%% -> builddir path here; overwritten at install time.
# %%CC_TAG%%              -> baked in at build time (never changes per install).
# %%PREFIX%%              -> builddir abs path here; overwritten at install time.
# -----------------------------------------------------------------------------
BUILT_SOURCES = \
    go.Raxpy.sh       \
    go.Rcopy.sh       \
    go.Rdot.sh        \
    go.Rgemm.sh       \
    go.Rgemm_cuda.sh  \
    go.Rgemv.sh       \
    go.Rsyrk.sh       \
    go.Rsyrk_cuda.sh  \
    go.Rgetrf.sh      \
    go.Rgetrf_cuda.sh \
    go.Rpotrf.sh      \
    go.Rpotrf_cuda.sh

CLEANFILES = $(BUILT_SOURCES)

go.%.sh: $(srcdir)/go.%.sh.in
	$(SED) \
	  -e "s|%%CPU_OS_DETECTION_SH%%|$(abs_top_srcdir)/misc/cpu_os_detection.sh|g" \
	  -e "s|%%CC_TAG%%|$(MPLAPACK_CC_TAG)|g" \
	  -e "s|%%PREFIX%%|$(abs_top_builddir)|g" \
	  $< > $@
	chmod +x $@

# -----------------------------------------------------------------------------
# Install destination
# -----------------------------------------------------------------------------
benchdir = $(prefix)/lib/$(target)/mplapack/benchmark

# -----------------------------------------------------------------------------
# Backend-specific build rules.
# Each Makefile.{backend}.am contains both BLAS and LAPACK benchmark targets.
# -----------------------------------------------------------------------------
if ENABLE_MPFR
include Makefile.mpfr.am
endif
if ENABLE_GMP
include Makefile.gmp.am
endif
if ENABLE_QD
include Makefile.qd.am
endif
if ENABLE_DD
include Makefile.dd.am
endif
if ENABLE_CUDA
include Makefile.dd_cuda.am
endif
if ENABLE_DOUBLE
include Makefile.double.am
endif
if ENABLE_BINARY128
include Makefile.binary128.am
endif
if ENABLE_BINARY80
include Makefile.binary80.am
endif
if ENABLE_BENCHMARK
include Makefile.blas.am
endif

# -----------------------------------------------------------------------------
# Install go.*.sh: copy generated scripts, then rewrite builddir paths to
# installed paths.  %%CC_TAG%% is already baked in; only paths need updating.
# -----------------------------------------------------------------------------
install-exec-hook:
	$(MKDIR_P) $(DESTDIR)$(benchdir)
	for f in $(BUILT_SOURCES); do \
	  $(INSTALL_SCRIPT) $$f $(DESTDIR)$(benchdir)/$$f; \
	  $(SED) -i \
	    -e "s|%%CPU_OS_DETECTION_SH%%|$(prefix)/misc/cpu_os_detection.sh|g" \
	    -e "s|%%PREFIX%%|$(prefix)|g" \
	    $(DESTDIR)$(benchdir)/$$f; \
	done
	bash $(top_builddir)/misc/add_ldlibpath.sh $(prefix) $(DESTDIR)$(benchdir)/go.*.sh
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(DESTDIR)$(benchdir) $(prefix)
endif

.PHONY: benchmarks
benchmarks: $(bench_PROGRAMS) $(BUILT_SOURCES)
	@echo "=== Starting MPLAPACK benchmarks ==="
	@failed=""; \
	for sh in $(BUILT_SOURCES); do \
	  echo "--- $$sh ---"; \
	  bash $$sh || failed="$$failed $$sh"; \
	done; \
	if [ -n "$$failed" ]; then \
	  echo "WARNING: the following scripts reported errors:$$failed" >&2; \
	fi; \
	echo "=== Done ==="
