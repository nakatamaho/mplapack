SUBDIRS = external mpblas mplapack
if !DISABLE_BENCHMARK
SUBDIRS += benchmark 
endif
if !DISABLE_EXAMPLES
SUBDIRS += examples
endif

ACLOCAL_AMFLAGS= -I m4

mplapack_includedir = $(includedir)/mplapack
mplapack_include_HEADERS = \
    include/mplapack_config.h \
    include/blas.h \
    include/lapack.h \
    include/blaswrapper.h

if ENABLE_DD
mplapack_include_HEADERS += \
    include/dd_complex.h \
    include/mpblas_dd.h \
    include/mplapack_dd.h \
    include/mplapack_utils_dd.h 
endif

if ENABLE_GMP
mplapack_include_HEADERS += \
    include/mpblas_gmp.h \
    include/mplapack_gmp.h \
    include/mplapack_utils_gmp.h \
    include/mpc_class.h
endif

if ENABLE_QD
mplapack_include_HEADERS += \
    include/mpblas_qd.h \
    include/mplapack_qd.h \
    include/mplapack_utils_qd.h \
    include/qd_complex.h
endif

if ENABLE_MPFR
mplapack_include_HEADERS += \
    include/mpblas_mpfr.h \
    include/mplapack_mpfr.h \
    include/mplapack_utils_mpfr.h \
    mpfrc++/mpcomplex.h \
    mpfrc++/mpreal.h
endif

if ENABLE_DOUBLE
mplapack_include_HEADERS += \
    include/mpblas_double.h \
    include/mplapack_double.h \
    include/mplapack_utils_double.h
endif

if ENABLE__FLOAT64X
mplapack_include_HEADERS += \
     include/mpblas__Float64x.h \
     include/mplapack__Float64x.h \
     include/mplapack_utils__Float64x.h.in
endif

if ENABLE__FLOAT128
mplapack_include_HEADERS += \
    include/mpblas__Float128.h \
    include/mplapack__Float128.h \
    include/mplapack_utils__Float128.h.in
endif

noinst_HEADERS = include/mplapack_debug.h include/mplapack_utils.h include/mpblas.h include/mplapack.h include/mplapack_generic.h include/mplapack_benchmark.h include/mplapack_benchmark__Float128.h include/mplapack_benchmark_dd.h include/mplapack_benchmark_double.h include/mplapack_benchmark_gmp.h include/mplapack_benchmark_mpfr.h include/mplapack_benchmark_qd.h include/mplapack_compare_debug.h include/mplapack_print_double.h include/mplapack_print.h include/mplapack_matgen.h include/mplapack_matgen__Float128.h include/mplapack_matgen__Float64x.h include/mplapack_matgen_dd.h include/mplapack_matgen_double.h include/mplapack_matgen_gmp.h include/mplapack_matgen_mpfr.hinclude/mplapack_matgen_qd.h include/mplapack_lin.h include/mplapack_lin__Float128.h include/mplapack_lin__Float64x.h include/mplapack_lin_dd.h include/mplapack_lin_double.h include/mplapack_lin_gmp.h include/mplapack_lin_mpfr.h include/mplapack_lin_qd.h include/mplapack_eig.h include/mplapack_eig__Float128.h include/mplapack_eig__Float64x.h include/mplapack_eig_dd.h include/mplapack_eig_double.h include/mplapack_eig_gmp.h include/mplapack_eig_mpfr.h include/mplapack_eig_qd.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = mplapack.pc

EXTRA_DIST = \
misc/cudalt.py \
misc/build_mplapack_for_mic.sh \
misc/prepare_mplapacklibs_for_mic.sh \
mpblas/optimized/dd/cuda/dd_real_cuda.h \
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
external/qd/patches/patch-* \
external/qd/patches/config* \
external/lapack/Makefile.am \
external/lapack/Makefile.in \
external/lapack/download/*.tgz \
external/lapack/download/*.tar.gz \
external/lapack/patches/patch-* \
external/lapack/patches/print.c \
mpfrc++/copying.lib mpfrc++/copying mpfrc++/makefile.org \
mpfrc++/org/copying mpfrc++/org/copying.lib mpfrc++/org/example.cpp \
mpfrc++/org/makefile mpfrc++/org/mpfrc++.tar.gz \
mpfrc++/org/mpreal.cpp mpfrc++/org/mpreal.h
