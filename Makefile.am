SUBDIRS = external mpblas mplapack
if ENABLE_BENCHMARK
SUBDIRS += benchmark 
endif
if !DISABLE_EXAMPLES
SUBDIRS += examples
endif

ACLOCAL_AMFLAGS= -I m4

mplapack_includedir = $(includedir)/mplapack
mplapack_include_HEADERS = \
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
    include/mplapack_utils_double.h \
    include/mplapack_double.h
endif

if ENABLE__FLOAT64X
mplapack_include_HEADERS += \
    include/mpblas__Float64x.h \
    include/mplapack__Float64x.h
endif

if ENABLE__FLOAT128
mplapack_include_HEADERS += \
    include/mpblas__Float128.h \
    include/mplapack__Float128.h
endif

nodist_mplapack_include_HEADERS = include/mplapack_config.h
if ENABLE__FLOAT64X
nodist_mplapack_include_HEADERS += include/mplapack_utils__Float64x.h
endif
if ENABLE__FLOAT128
nodist_mplapack_include_HEADERS += include/mplapack_utils__Float128.h
endif

noinst_HEADERS = include/mplapack_debug.h include/mplapack_utils.h include/mpblas.h include/mplapack.h include/mplapack_generic.h include/mplapack_benchmark.h include/mplapack_benchmark__Float128.h include/mplapack_benchmark_dd.h include/mplapack_benchmark_double.h include/mplapack_benchmark_gmp.h include/mplapack_benchmark_mpfr.h include/mplapack_benchmark_qd.h include/mplapack_compare_debug.h include/mplapack_print_double.h include/mplapack_print.h include/mplapack_matgen.h include/mplapack_matgen__Float128.h include/mplapack_matgen__Float64x.h include/mplapack_matgen_dd.h include/mplapack_matgen_double.h include/mplapack_matgen_gmp.h include/mplapack_matgen_mpfr.h include/mplapack_matgen_qd.h include/mplapack_lin.h include/mplapack_lin__Float128.h include/mplapack_lin__Float64x.h include/mplapack_lin_dd.h include/mplapack_lin_double.h include/mplapack_lin_gmp.h include/mplapack_lin_mpfr.h include/mplapack_lin_qd.h include/mplapack_eig.h include/mplapack_eig__Float128.h include/mplapack_eig__Float64x.h include/mplapack_eig_dd.h include/mplapack_eig_double.h include/mplapack_eig_gmp.h include/mplapack_eig_mpfr.h include/mplapack_eig_qd.h include/mplapack_common_sslct.h include/mplapack_common_mn.h include/mplapack_benchmark.h include/mplapack_benchmark__Float128.h include/mplapack_benchmark__Float64x.h include/mplapack_benchmark_dd.h include/mplapack_benchmark_double.h include/mplapack_benchmark_gmp.h include/mplapack_benchmark_mpfr.h include/mplapack_benchmark_qd.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = mplapack.pc

EXTRA_DIST = \
include/mplapack_config.h.in \
benchmark/benchall.sh \
benchmark/mpblas/Raxpy1.plt.in \
benchmark/mpblas/Raxpy2.plt.in \
benchmark/mpblas/Raxpy3.plt.in \
benchmark/mpblas/Rcopy1.plt.in \
benchmark/mpblas/Rcopy2.plt.in \
benchmark/mpblas/Rcopy3.plt.in \
benchmark/mpblas/Rdot1.plt.in \
benchmark/mpblas/Rdot2.plt.in \
benchmark/mpblas/Rdot3.plt.in \
benchmark/mpblas/Rgemv1.plt.in \
benchmark/mpblas/Rgemv2.plt.in \
benchmark/mpblas/Rgemv3.plt.in \
benchmark/mpblas/Rgemm1.plt.in \
benchmark/mpblas/Rgemm2.plt.in \
benchmark/mpblas/Rgemm3.plt.in \
benchmark/mpblas/Rgemm_cuda.plt.in \
benchmark/mpblas/Rsyrk1.plt.in \
benchmark/mpblas/Rsyrk2.plt.in \
benchmark/mpblas/Rsyrk3.plt.in \
benchmark/mpblas/Rsyrk_cuda.plt.in \
benchmark/mpblas/go.Raxpy.sh \
benchmark/mpblas/go.Rcopy.sh \
benchmark/mpblas/go.Rdot.sh \
benchmark/mpblas/go.Rgemv.sh \
benchmark/mpblas/go.Rgemm.sh \
benchmark/mpblas/go.Rgemm_cuda.sh \
benchmark/mpblas/go.Rsyrk.sh \
benchmark/mpblas/go.Rsyrk_cuda.sh \
benchmark/mplapack/Rgetrf1.plt.in \
benchmark/mplapack/Rgetrf2.plt.in \
benchmark/mplapack/Rgetrf3.plt.in \
benchmark/mplapack/Rpotrf1.plt.in \
benchmark/mplapack/Rpotrf2.plt.in \
benchmark/mplapack/Rpotrf3.plt.in \
benchmark/mplapack/Rgetrf_cuda.plt.in \
benchmark/mplapack/Rpotrf_cuda.plt.in \
benchmark/mplapack/go.Rgetrf.sh \
benchmark/mplapack/go.Rgetrf_cuda.sh \
benchmark/mplapack/go.Rpotrf.sh \
benchmark/mplapack/go.Rpotrf_cuda.sh \
misc/cudalt.py \
misc/build_mplapack_for_mic.sh \
misc/prepare_mplapacklibs_for_mic.sh \
misc/fix_dylib_macOS.sh \
misc/add_ldlibpath.sh \
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
external/lapack/Makefile.am \
external/lapack/Makefile.in \
external/lapack/download/*.tar.gz \
external/lapack/patches/patch-* \
external/lapack/patches/print.c \
mpfrc++/copying.lib mpfrc++/copying mpfrc++/makefile.org \
mpfrc++/org/copying mpfrc++/org/copying.lib mpfrc++/org/example.cpp \
mpfrc++/org/makefile mpfrc++/org/mpfrc++.tar.gz \
mpfrc++/org/mpreal.cpp mpfrc++/org/mpreal.h \
examples/mpblas/Makefile.freebsd.in \
examples/mpblas/Makefile.linux.in \
examples/mpblas/Makefile.linux_cuda.in \
examples/mpblas/Makefile.macos.in \
examples/mpblas/Makefile.mingw.in \
examples/mplapack/*/Makefile.freebsd.in \
examples/mplapack/*/Makefile.linux.in \
examples/mplapack/*/Makefile.macos.in \
examples/mplapack/*/Makefile.mingw.in \
examples/mplapack/04_NonsymmetricEigenproblems/Matrix_*.txt \
examples/mplapack/05_SingularValueDecomposition/Matrix_*.txt \
mplapack/test/eig/Cbak.in \
mplapack/test/eig/Cbal.in \
mplapack/test/eig/Cbb.in \
mplapack/test/eig/Cec.in \
mplapack/test/eig/Ced.in \
mplapack/test/eig/Cgbak.in \
mplapack/test/eig/Cgbal.in \
mplapack/test/eig/Cgd.in \
mplapack/test/eig/Cgg.in \
mplapack/test/eig/Csb.in \
mplapack/test/eig/Csg.in \
mplapack/test/eig/Rbak.in \
mplapack/test/eig/Rbal.in \
mplapack/test/eig/Rbb.in \
mplapack/test/eig/Rec.in \
mplapack/test/eig/Red.in \
mplapack/test/eig/Rgbak.in \
mplapack/test/eig/Rgbal.in \
mplapack/test/eig/Rgd.in \
mplapack/test/eig/Rgg.in \
mplapack/test/eig/Rsb.in \
mplapack/test/eig/Rsg.in \
mplapack/test/eig/csd.in \
mplapack/test/eig/glm.in \
mplapack/test/eig/gqr.in \
mplapack/test/eig/gsv.in \
mplapack/test/eig/lse.in \
mplapack/test/eig/nep.in \
mplapack/test/eig/se2.in \
mplapack/test/eig/sep.in \
mplapack/test/eig/svd.in \
mplapack/test/eig/test_eig_all.sh \
mplapack/test/eig/test_eig_all_mingw.sh \
mplapack/test/lin/Ctest.in \
mplapack/test/lin/Ctest_rfp.in \
mplapack/test/lin/Rtest.in \
mplapack/test/lin/Rtest_rfp.in \
mplapack/test/lin/test_lin_all.sh \
mplapack/test/lin/test_lin_all_mingw.sh
