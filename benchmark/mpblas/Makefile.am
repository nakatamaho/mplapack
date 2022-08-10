bench_PROGRAMS =
bench_DATA = Raxpy1.plt Raxpy2.plt \
Rgemm1.plt Rgemm2.plt Rgemm_cuda.plt \
Rsyrk1.plt Rsyrk2.plt Rsyrk_cuda.plt \
go.Raxpy.sh go.Rgemm.sh go.Rgemm_cuda.sh go.Rsyrk.sh go.Rsyrk_cuda.sh

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

if DO_BLASBENCH
include Makefile.blas.am
endif

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(benchdir) $(prefix)
endif

endif
