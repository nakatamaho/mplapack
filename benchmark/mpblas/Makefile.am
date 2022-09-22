bench_PROGRAMS =
bench_DATA = \
Raxpy1.plt Raxpy2.plt Raxpy3.plt Rcopy1.plt Rcopy2.plt Rcopy3.plt Rdot1.plt Rdot2.plt Rdot3.plt \
Rgemm1.plt Rgemm2.plt Rgemm3.plt \
Rgemm_cuda.plt \
Rgemv1.plt Rgemv2.plt Rgemv3.plt \
Rsyrk1.plt Rsyrk2.plt Rsyrk3.plt \
Rsyrk_cuda.plt \
go.Raxpy.sh go.Rcopy.sh go.Rdot.sh go.Rgemm.sh \
go.Rgemm_cuda.sh \
go.Rgemv.sh go.Rsyrk.sh \
go.Rsyrk_cuda.sh

benchdir=$(prefix)/lib/$(target)/mplapack/benchmark

if !IS_MINGW

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

endif
