bench_PROGRAMS =
bench_DATA = Rgetrf1.plt.in Rgetrf2.plt.in Rgetrf3.plt.in Rgetrf_cuda.plt.in \
Rpotrf1.plt.in Rpotrf2.plt.in Rpotrf3.plt.in Rpotrf_cuda.plt.in \
go.Rgetrf.sh go.Rgetrf_cuda.sh go.Rpotrf.sh go.Rpotrf_cuda.sh

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

if ENABLE_CUDA
include Makefile.dd_cuda.am
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
include Makefile.lapack.am
endif

install-data-hook:
	bash $(top_builddir)/misc/add_ldlibpath.sh $(prefix) $(benchdir)/go.*.sh
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(benchdir) $(prefix)
endif
