bench_PROGRAMS =
benchdir=$(prefix)/share/mplapack/bench

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

if DO_LAPACKBENCH
include Makefile.lapack.am
endif

endif
