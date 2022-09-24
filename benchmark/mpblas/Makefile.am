bench_PROGRAMS =
bench_DATA = \
Raxpy1.plt.in Raxpy2.plt.in Raxpy3.plt.in \
Rcopy1.plt.in Rcopy2.plt.in Rcopy3.plt.in \
Rdot1.plt.in  Rdot2.plt.in  Rdot3.plt.in \
Rgemm1.plt.in Rgemm2.plt.in Rgemm3.plt.in \
Rgemm_cuda.plt.in \
Rgemv1.plt.in Rgemv2.plt.in Rgemv3.plt.in \
Rsyrk1.plt.in Rsyrk2.plt.in Rsyrk3.plt.in \
Rsyrk_cuda.plt.in \
go.Raxpy.sh go.Rcopy.sh go.Rdot.sh go.Rgemm.sh \
go.Rgemm_cuda.sh \
go.Rgemv.sh go.Rsyrk.sh \
go.Rsyrk_cuda.sh

benchdir=$(prefix)/lib/$(target)/mplapack/benchmark

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

if ENABLE_DOUBLE
include Makefile.double.am
endif

if ENABLE__FLOAT128
include Makefile._Float128.am
endif

if ENABLE__FLOAT64X
include Makefile._Float64x.am
endif

if ENABLE_BENCHMARK
include Makefile.blas.am
endif

install-data-hook:
	bash $(top_builddir)/misc/add_ldlibpath.sh $(prefix) $(benchdir)/go.*.sh
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(benchdir) $(prefix)
endif
