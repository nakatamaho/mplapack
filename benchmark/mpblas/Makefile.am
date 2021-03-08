bench_PROGRAMS =
benchdir=$(prefix)/share/mplapack/bench

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

if ENABLE_BINARY128
include Makefile.binary128.am
endif

if DO_BLASBENCH
include Makefile.blas.am
endif
