bench_PROGRAMS =
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

endif


