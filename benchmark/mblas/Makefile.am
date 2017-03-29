bench_PROGRAMS =
benchdir=$(prefix)/share/mpack/bench

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

if ENABLE___FLOAT128
include Makefile.__float128.am
endif

if DO_BLASBENCH
include Makefile.blas.am
endif
