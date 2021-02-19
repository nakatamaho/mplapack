SUBDIRS = external mpfrc++ mblas mlapack 
if !DISABLE_BENCHMARK
SUBDIRS += benchmark 
endif
if !DISABLE_EXAMPLES
SUBDIRS += examples
endif

ACLOCAL_AMFLAGS= -I m4

mpack_includedir = $(includedir)/mpack
mpack_include_HEADERS = \
    include/mpack_config.h \
    include/blas.h \
    include/lapack.h \
    include/blaswrapper.h

if ENABLE_DD
mpack_include_HEADERS += \
    include/dd_complex.h \
    include/mblas_dd.h \
    include/mlapack_dd.h \
    include/mutils_dd.h 
endif

if ENABLE_GMP
mpack_include_HEADERS += \
    include/mblas_gmp.h \
    include/mlapack_gmp.h \
    include/mutils_gmp.h \
    include/mpc_class.h
endif

if ENABLE_QD
mpack_include_HEADERS += \
    include/mblas_qd.h \
    include/mlapack_qd.h \
    include/mutils_qd.h \
    include/qd_complex.h
endif

if ENABLE_MPFR
mpack_include_HEADERS += \
    include/mblas_mpfr.h \
    include/mlapack_mpfr.h \
    include/mutils_mpfr.h \
    mpfrc++/mpcomplex.h \
    mpfrc++/mpreal.h
endif

if ENABLE_DOUBLE
mpack_include_HEADERS += \
    include/mblas_double.h \
    include/mlapack_double.h \
    include/mutils_double.h
endif

if ENABLE_LONGDOUBLE
mpack_include_HEADERS += \
     include/mblas_longdouble.h \
     include/mlapack_longdouble.h \
     include/mutils_longdouble.h
endif

if ENABLE___FLOAT128
mpack_include_HEADERS += \
    include/mblas___float128.h \
    include/mlapack___float128.h \
    include/mutils___float128.h
endif

noinst_HEADERS = include/mpack_debug.h include/mutils.h include/mblas.h include/mlapack.h include/mlapack_generic.h include/mpack_benchmark.h include/mpack_benchmark___float128.h include/mpack_benchmark_dd.h include/mpack_benchmark_double.h include/mpack_benchmark_gmp.h include/mpack_benchmark_mpfr.h include/mpack_benchmark_qd.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = mpack.pc

EXTRA_DIST = \
misc/cudalt.py \
misc/build_mpack_for_mic.sh \
misc/prepare_mpacklibs_for_mic.sh \
mblas/optimized/dd/cuda/dd_real_cuda.h \
external/gmp/download/*.tar.bz2 \
external/gmp/Makefile.am \
external/gmp/Makefile.in \
external/mpfr/download/*.tar.bz2 \
external/mpfr/Makefile.am \
external/mpfr/Makefile.in \
external/mpc/download/*.tar.gz \
external/mpc/Makefile.am \
external/mpc/Makefile.in \
external/qd/download/*.tar.gz \
external/qd/Makefile.am \
external/qd/Makefile.in \
external/qd/patches/patch-mingw \
external/qd/patches/patch-noopt \
external/lapack/Makefile.am \
external/lapack/Makefile.in \
external/lapack/download/*.tgz \
external/lapack/patches/patch-mingw \
external/lapack/patches/patch-mic \
mpfrc++/copying.lib mpfrc++/copying mpfrc++/makefile.org \
mpfrc++/org/copying mpfrc++/org/copying.lib mpfrc++/org/example.cpp \
mpfrc++/org/makefile  mpfrc++/org/mpfrc++.tar.gz \
mpfrc++/org/mpreal.cpp  mpfrc++/org/mpreal.h

if !DISABLE_EXAMPLES
EXTRA_DIST += \
examples/mblas/Makefile.mingw \
examples/mblas/Makefile.freebsd \
examples/mblas/Makefile.linux \
examples/mblas/Makefile.macosx \
examples/mlapack/Makefile.mingw \
examples/mlapack/Makefile.freebsd \
examples/mlapack/Makefile.linux \
examples/mlapack/Makefile.macosx
endif

