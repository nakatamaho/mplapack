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

nodist_mplapack_include_HEADERS =
if ENABLE__FLOAT64X
nodist_mplapack_include_HEADERS += include/mplapack_utils__Float64x.h
endif
if ENABLE__FLOAT128
nodist_mplapack_include_HEADERS += include/mplapack_utils__Float128.h
endif

noinst_HEADERS = include/mplapack_utils__Float128.h.in include/mplapack_utils__Float64x.h.in include/mplapack_debug.h include/mplapack_utils.h include/mpblas.h include/mplapack.h include/mplapack_generic.h include/mplapack_benchmark.h include/mplapack_benchmark__Float128.h include/mplapack_benchmark_dd.h include/mplapack_benchmark_double.h include/mplapack_benchmark_gmp.h include/mplapack_benchmark_mpfr.h include/mplapack_benchmark_qd.h include/mplapack_compare_debug.h include/mplapack_print_double.h include/mplapack_print.h include/mplapack_matgen.h include/mplapack_matgen__Float128.h include/mplapack_matgen__Float64x.h include/mplapack_matgen_dd.h include/mplapack_matgen_double.h include/mplapack_matgen_gmp.h include/mplapack_matgen_mpfr.h include/mplapack_matgen_qd.h include/mplapack_lin.h include/mplapack_lin__Float128.h include/mplapack_lin__Float64x.h include/mplapack_lin_dd.h include/mplapack_lin_double.h include/mplapack_lin_gmp.h include/mplapack_lin_mpfr.h include/mplapack_lin_qd.h include/mplapack_eig.h include/mplapack_eig__Float128.h include/mplapack_eig__Float64x.h include/mplapack_eig_dd.h include/mplapack_eig_double.h include/mplapack_eig_gmp.h include/mplapack_eig_mpfr.h include/mplapack_eig_qd.h mplapack/test/fable/fem/do.hpp mplapack/test/fable/fem/intrinsics.hpp mplapack/test/fable/fem/arr_size.hpp mplapack/test/fable/fem/main.hpp mplapack/test/fable/fem/arr_and_str_indices.hpp mplapack/test/fable/fem/open_chain.hpp mplapack/test/fable/fem/format.hpp mplapack/test/fable/fem/arr.hpp mplapack/test/fable/fem/major_types.hpp mplapack/test/fable/fem/data.hpp mplapack/test/fable/fem/utils mplapack/test/fable/fem/utils/string_to_double_fmt.hpp mplapack/test/fable/fem/utils/int_types.hpp mplapack/test/fable/fem/utils/char.hpp mplapack/test/fable/fem/utils/double_to_string.hpp mplapack/test/fable/fem/utils/string.hpp mplapack/test/fable/fem/utils/simple_streams.hpp mplapack/test/fable/fem/utils/string_to_double.hpp mplapack/test/fable/fem/utils/path.hpp mplapack/test/fable/fem/utils/misc.hpp mplapack/test/fable/fem/utils/real_as_string.hpp mplapack/test/fable/fem/utils/random.hpp mplapack/test/fable/fem/utils/token.hpp mplapack/test/fable/fem/utils/equivalence.hpp mplapack/test/fable/fem/str_arr_ref.hpp mplapack/test/fable/fem/size_t.hpp mplapack/test/fable/fem/read.hpp mplapack/test/fable/fem/error_utils.hpp mplapack/test/fable/fem/str_ref.hpp mplapack/test/fable/fem/file_positioning_chain.hpp mplapack/test/fable/fem/zero.hpp mplapack/test/fable/fem/io_exceptions.hpp mplapack/test/fable/fem/variant.hpp mplapack/test/fable/fem/intrinsics_extra.hpp mplapack/test/fable/fem/inquire_chain.hpp mplapack/test/fable/fem/common.hpp mplapack/test/fable/fem/arr_ref.hpp mplapack/test/fable/fem/data_of_type.hpp mplapack/test/fable/fem/str.hpp mplapack/test/fable/fem/data_type_star.hpp mplapack/test/fable/fem/star.hpp mplapack/test/fable/fem/dimension.hpp mplapack/test/fable/fem/close_chain.hpp mplapack/test/fable/fem/stop.hpp mplapack/test/fable/fem/io.hpp mplapack/test/fable/fem/write.hpp mplapack/test/fable/tbxx/time_accu.hpp mplapack/test/fable/tbxx/optional_copy.hpp mplapack/test/fable/tbxx/error_utils.hpp mplapack/test/fable/tbxx/pretty_type_name.hpp mplapack/test/fable/tbxx/libc_backtrace.hpp mplapack/test/fable/noexcept_false.hpp mplapack/test/fable/fem.hpp mplapack/test/fable/LICENSE.txt include/mplapack_common_sslct.h include/mplapack_common_mn.h include/mplapack_benchmark.h include/mplapack_benchmark__Float128.h include/mplapack_benchmark__Float64x.h include/mplapack_benchmark_dd.h include/mplapack_benchmark_double.h include/mplapack_benchmark_gmp.h include/mplapack_benchmark_mpfr.h include/mplapack_benchmark_qd.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = mplapack.pc

EXTRA_DIST = \
benchmark/benchall.sh \
benchmark/mpblas/Raxpy1.plt \
benchmark/mpblas/Raxpy2.plt \
benchmark/mpblas/Raxpy3.plt \
benchmark/mpblas/Rcopy1.plt \
benchmark/mpblas/Rcopy2.plt \
benchmark/mpblas/Rcopy3.plt \
benchmark/mpblas/Rdot1.plt \
benchmark/mpblas/Rdot2.plt \
benchmark/mpblas/Rdot3.plt \
benchmark/mpblas/Rgemv1.plt \
benchmark/mpblas/Rgemv2.plt \
benchmark/mpblas/Rgemv3.plt \
benchmark/mpblas/Rgemm1.plt \
benchmark/mpblas/Rgemm2.plt \
benchmark/mpblas/Rgemm3.plt \
benchmark/mpblas/Rgemm_cuda.plt \
benchmark/mpblas/Rsyrk1.plt \
benchmark/mpblas/Rsyrk2.plt \
benchmark/mpblas/Rsyrk3.plt \
benchmark/mpblas/Rsyrk_cuda.plt \
benchmark/mpblas/go.Raxpy.sh \
benchmark/mpblas/go.Rcopy.sh \
benchmark/mpblas/go.Rdot.sh \
benchmark/mpblas/go.Rgemv.sh \
benchmark/mpblas/go.Rgemm.sh \
benchmark/mpblas/go.Rgemm_cuda.sh \
benchmark/mpblas/go.Rsyrk.sh \
benchmark/mpblas/go.Rsyrk_cuda.sh \
benchmark/mplapack/Rgetrf1.plt \
benchmark/mplapack/Rgetrf2.plt \
benchmark/mplapack/Rgetrf3.plt \
benchmark/mplapack/Rpotrf1.plt \
benchmark/mplapack/Rpotrf2.plt \
benchmark/mplapack/Rpotrf3.plt \
benchmark/mplapack/Rgetrf_cuda.plt \
benchmark/mplapack/Rpotrf_cuda.plt \
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
mplapack/test/eig/Cbal_double.in \
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
mplapack/test/eig/Rbal_double.in \
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