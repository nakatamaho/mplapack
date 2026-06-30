SUBDIRS = external mpblas mplapack
if ENABLE_BENCHMARK
SUBDIRS += benchmark
endif
if !DISABLE_EXAMPLES
SUBDIRS += examples
endif

ACLOCAL_AMFLAGS = -I m4

miscdir = $(prefix)/misc
dist_misc_SCRIPTS = misc/cpu_os_detection.sh

mplapack_includedir = $(includedir)/mplapack
mplapack_include_HEADERS = \
    include/blas.h \
    include/lapack.h \
    include/blaswrapper.h

if ENABLE_DD
mplapack_include_HEADERS += \
    include/dd_complex.h \
    include/mpblas_dd.h \
    include/mplapack_arithmetic_params_dd.h \
    include/mplapack_dd.h \
    include/mplapack_utils_dd.h
endif

if ENABLE_GMP
mplapack_include_HEADERS += \
    include/mpblas_gmp.h \
    include/mplapack_arithmetic_params_gmp.h \
    include/mplapack_gmp.h \
    include/mplapack_gmp_transcendents.h \
    include/mplapack_utils_gmp.h \
    include/mpc_class.h
endif

if ENABLE_QD
mplapack_include_HEADERS += \
    include/mpblas_qd.h \
    include/mplapack_arithmetic_params_qd.h \
    include/mplapack_qd.h \
    include/mplapack_utils_qd.h \
    include/qd_complex.h
endif

if ENABLE_MPFR
mplapack_include_HEADERS += \
    include/mpblas_mpfr.h \
    include/mplapack_arithmetic_params_mpfr.h \
    include/mplapack_mpfr.h \
    include/mplapack_utils_mpfr.h \
    mpfrc++/mpcomplex.h \
    mpfrc++/mpreal.h
endif

if ENABLE_DOUBLE
mplapack_include_HEADERS += \
    include/mpblas_double.h \
    include/mplapack_arithmetic_params_double.h \
    include/mplapack_utils_double.h \
    include/mplapack_double.h
endif

if ENABLE_BINARY80
mplapack_include_HEADERS += \
    include/mpblas_binary80.h \
    include/mplapack_arithmetic_params_binary80.h \
    include/mplapack_binary80.h \
    include/mplapack_utils_binary80.h
endif

if ENABLE_BINARY128
mplapack_include_HEADERS += \
    include/mpblas_binary128.h \
    include/mplapack_arithmetic_params_binary128.h \
    include/mplapack_binary128.h \
    include/mplapack_utils_binary128.h
endif

nodist_mplapack_include_HEADERS = include/mplapack_config.h
if ENABLE_BINARY80
nodist_mplapack_include_HEADERS += include/mplapack_utils_binary80.h
endif
if ENABLE_BINARY128
nodist_mplapack_include_HEADERS += include/mplapack_utils_binary128.h
endif

# Internal headers used by tests/benchmarks; keep them in dist tarballs.
noinst_HEADERS = \
    include/mplapack_arithmetic_params.h \
    include/mplapack_debug.h \
    include/mplapack_utils.h \
    include/mpblas.h \
    include/mplapack.h \
    include/mplapack_generic.h \
    include/mplapack_compare_debug.h \
    include/mplapack_print_double.h \
    include/mplapack_hex_helpers.h \
    include/mplapack_print.h \
    include/mplapack_common_sslct.h \
    include/mplapack_common_mn.h \
    include/mplapack_benchmark.h \
    include/mplapack_benchmark_binary128.h \
    include/mplapack_benchmark_binary80.h \
    include/mplapack_benchmark_dd.h \
    include/mplapack_benchmark_double.h \
    include/mplapack_benchmark_gmp.h \
    include/mplapack_benchmark_mpfr.h \
    include/mplapack_benchmark_qd.h \
    include/mplapack_matgen.h \
    include/mplapack_matgen_binary128.h \
    include/mplapack_matgen_binary80.h \
    include/mplapack_matgen_dd.h \
    include/mplapack_matgen_double.h \
    include/mplapack_matgen_gmp.h \
    include/mplapack_matgen_mpfr.h \
    include/mplapack_matgen_qd.h \
    include/mplapack_lin.h \
    include/mplapack_lin_binary128.h \
    include/mplapack_lin_binary80.h \
    include/mplapack_lin_dd.h \
    include/mplapack_lin_double.h \
    include/mplapack_lin_gmp.h \
    include/mplapack_lin_mpfr.h \
    include/mplapack_lin_qd.h \
    include/mplapack_eig.h \
    include/mplapack_eig_binary128.h \
    include/mplapack_eig_binary80.h \
    include/mplapack_eig_dd.h \
    include/mplapack_eig_double.h \
    include/mplapack_eig_gmp.h \
    include/mplapack_eig_mpfr.h \
    include/mplapack_eig_qd.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = mplapack.pc

EXTRA_DIST = \
benchmark/Rgetrf1.plt.in \
benchmark/Rgetrf2.plt.in \
benchmark/Rgetrf3.plt.in \
benchmark/Rgetrf_cuda.plt.in \
benchmark/Rpotrf1.plt.in \
benchmark/Rpotrf2.plt.in \
benchmark/Rpotrf3.plt.in \
benchmark/Rpotrf_cuda.plt.in \
benchmark/benchall.sh \
examples/mplapack/04_NonsymmetricEigenproblems/Matrix_01.txt \
examples/mplapack/04_NonsymmetricEigenproblems/Matrix_02.txt \
examples/mplapack/04_NonsymmetricEigenproblems/Matrix_03.txt \
examples/mplapack/04_NonsymmetricEigenproblems/Matrix_04.txt \
examples/mplapack/04_NonsymmetricEigenproblems/Matrix_05.txt \
examples/mplapack/04_NonsymmetricEigenproblems/Matrix_06.txt \
examples/mplapack/04_NonsymmetricEigenproblems/Matrix_07.txt \
examples/mplapack/05_SingularValueDecomposition/Matrix_01.txt \
examples/mplapack/05_SingularValueDecomposition/Matrix_02.txt \
examples/mplapack/05_SingularValueDecomposition/Matrix_03.txt \
examples/mpblas/Makefile.linux.inteloneAPI.in \
examples/mpblas/generic/Makefile.linux.inteloneAPI.in \
examples/mplapack/generic/Makefile.linux.inteloneAPI.in \
examples/mplapack/00_LinearEquations/Makefile.linux.inteloneAPI.in \
examples/mplapack/03_SymmetricEigenproblems/Makefile.linux.inteloneAPI.in \
examples/mplapack/04_NonsymmetricEigenproblems/Makefile.linux.inteloneAPI.in \
examples/mplapack/05_SingularValueDecomposition/Makefile.linux.inteloneAPI.in \
external/dlfcn-win32/Makefile.in \
external/dlfcn-win32/Makefile.am \
external/dlfcn-win32/download/dlfcn-win32-1.4.2.tar.gz \
external/distfiles.sha256 \
external/gmp/Makefile.am \
external/gmp/Makefile.in \
external/gmp/download/gmp-6.3.0.tar.xz \
external/gmp/patches/gcc-15.patch \
external/gmp/patches/gmp-6.3.0_libtool-2.4.7.patch \
external/lapack/Makefile.am \
external/lapack/Makefile.in \
external/lapack/download/lapack-3.12.1.tar.gz \
external/lapack/download/lapack-3.9.1.tar.gz \
external/lapack/patches/3.12.1/patch-cabs1_unification \
external/lapack/patches/3.12.1/patch-cabsmax \
external/lapack/patches/3.12.1/patch-fix-testparameters \
external/lapack/patches/3.12.1/patch-iseed \
external/lapack/patches/3.12.1/patch-mingw \
external/lapack/patches/3.12.1/patch-print \
external/lapack/patches/3.12.1/patch-redundant-dble-real \
external/lapack/patches/3.12.1/patch-stebz-error-propagation \
external/lapack/patches/3.12.1/patch-trsyl3-ldswork-overwrite \
external/lapack/patches/3.12.1/patch-xERRTSQR.f \
external/lapack/patches/3.12.1/patch-xLAQP2RK.f \
external/lapack/patches/3.12.1/patch-xpgvx-info-to-iinfo \
external/lapack/patches/3.12.1/patch-zunbdb3-conjugate-x21 \
external/lapack/patches/3.12.1/patch-xxxGVD.f \
external/lapack/patches/3.9.1/patch-ddrvst.f \
external/lapack/patches/3.9.1/patch-ddrvst2stg.f \
external/lapack/patches/3.9.1/patch-dlatm4.f \
external/lapack/patches/3.9.1/patch-iseed \
external/lapack/patches/3.9.1/patch-mic \
external/lapack/patches/3.9.1/patch-mingw \
external/lapack/patches/3.9.1/patch-print \
external/lapack/patches/3.9.1/patch-zchkhe_aa_2stage.f \
external/lapack/patches/3.9.1/patch-zlqt04.f \
external/lapack/patches/3.9.1/patch-zlqt05.f \
external/lapack/patches/3.9.1/patch-ztsqr01.f \
external/lapack/patches/3.9.1/print.c \
external/lapack/patches/print.c \
external/mpc/Makefile.am \
external/mpc/Makefile.in \
external/mpc/download/mpc-1.4.1.tar.xz \
external/mpfr/Makefile.am \
external/mpfr/Makefile.in \
external/mpfr/download/allpatches \
external/mpfr/download/mpfr-4.2.2.tar.bz2 \
external/openblas/Makefile.am \
external/openblas/Makefile.in \
external/openblas/download/OpenBLAS-0.3.33.tar.gz \
external/openblas/patches/patch-mingw \
external/qd/Makefile.am \
external/qd/Makefile.in \
external/qd/download/qd-2.3.24.tar.gz \
external/qd/patches/patch-configure \
external/qd/patches/patch-mingw \
external/qd/patches/patch-no-doc \
external/qd/patches/patch-noopt \
external/qd/patches/patch-overflow \
fable/fem.hpp \
fable/noexcept_false.hpp \
fable/fem/do.hpp \
fable/fem/intrinsics.hpp \
fable/fem/arr_size.hpp \
fable/fem/main.hpp \
fable/fem/arr_and_str_indices.hpp \
fable/fem/open_chain.hpp \
fable/fem/format.hpp \
fable/fem/arr.hpp \
fable/fem/major_types.hpp \
fable/fem/data.hpp \
fable/fem/utils/string_to_double_fmt.hpp \
fable/fem/utils/int_types.hpp \
fable/fem/utils/char.hpp \
fable/fem/utils/double_to_string.hpp \
fable/fem/utils/string.hpp \
fable/fem/utils/simple_streams.hpp \
fable/fem/utils/string_to_double.hpp \
fable/fem/utils/path.hpp \
fable/fem/utils/misc.hpp \
fable/fem/utils/real_as_string.hpp \
fable/fem/utils/random.hpp \
fable/fem/utils/token.hpp \
fable/fem/utils/equivalence.hpp \
fable/fem/str_arr_ref.hpp \
fable/fem/size_t.hpp \
fable/fem/read.hpp \
fable/fem/error_utils.hpp \
fable/fem/str_ref.hpp \
fable/fem/file_positioning_chain.hpp \
fable/fem/zero.hpp \
fable/fem/io_exceptions.hpp \
fable/fem/variant.hpp \
fable/fem/intrinsics_extra.hpp \
fable/fem/inquire_chain.hpp \
fable/fem/common.hpp \
fable/fem/arr_ref.hpp \
fable/fem/data_of_type.hpp \
fable/fem/str.hpp \
fable/fem/data_type_star.hpp \
fable/fem/star.hpp \
fable/fem/dimension.hpp \
fable/fem/close_chain.hpp \
fable/fem/stop.hpp \
fable/fem/io.hpp \
fable/fem/write.hpp \
fable/tbxx/error_utils.hpp \
fable/tbxx/libc_backtrace.hpp \
fable/tbxx/optional_copy.hpp \
fable/tbxx/pretty_type_name.hpp \
fable/tbxx/time_accu.hpp \
include/mplapack_symbol_resolver.h \
misc/add_ldlibpath.sh \
misc/ccache_set_min_maxsize.sh \
misc/filter-example-programs.sh \
release/buildtest_tier1_macos_amd64.sh \
release/buildtest_tier1_macos_arm64.sh \
release/buildtest_tier1_mingw64_amd64.sh \
release/buildtest_tier1_ubuntu_amd64.sh \
release/buildtest_tier1_ubuntu_inteloneapi_amd64.sh \
release/buildtest_tier1_ubuntu_arm64.sh \
release/buildtest_tier1_debian_i386.sh \
release/buildtest_tier3_ubuntu2604_cxxstd_amd64.sh \
release/check-installed-benchmarks.sh \
release/check-installed-examples.sh \
release/docker/common/entrypoint.sh \
release/docker/common/tarball-smoke.sh \
release/docker/distcheck/Dockerfile.debian-i386 \
release/docker/distcheck/Dockerfile.intel \
release/docker/distcheck/Dockerfile.mingw \
release/docker/distcheck/Dockerfile.ubuntu \
release/docker/distcheck/Dockerfile.ubuntu-nvidia \
release/docker/distcheck/Dockerfile.ubuntu24.04.nvidia \
release/docker/distcheck/run-ubuntu-distcheck.sh \
release/docker/distcheck/run-ubuntu-intel-distcheck.sh \
release/docker/distcheck/run-ubuntu-mingw64-distcheck.sh \
release/docker/distcheck/run-ubuntu-nvidia-distcheck.sh \
release/docker/matrix/Dockerfile.alpine \
release/docker/matrix/Dockerfile.cuda \
release/docker/matrix/Dockerfile.debian \
release/docker/matrix/Dockerfile.intel \
release/docker/matrix/Dockerfile.mingw \
release/docker/matrix/Dockerfile.redhat \
release/docker/matrix/Dockerfile.redhat-el8 \
release/docker/matrix/Dockerfile.suse \
release/docker/matrix/run-ubuntu-intel-build.sh \
release/docker/matrix/run-ubuntu-mingw64-build.sh \
release/docker/tarball/Dockerfile.ubuntu \
release/docker/tier3/Dockerfile.ubuntu-cxxstd \
release/docker/tier3/run-ubuntu-cxxstd-build.sh \
release/make-source-snapshot.sh \
release/run-remote-linux-docker.sh \
release/run-remote-tarball-docker-worker.sh \
release/run-remote-macos.sh \
release/run-remote-targets-by-host.sh \
release/clean-remote-workdirs.sh \
misc/cudalt.py \
misc/fix_dylib_macOS.sh \
misc/num_diff.py \
misc/reconfig.bullseye.sh \
misc/reconfig.develop.sh \
misc/reconfig.macOS.sh \
misc/reconfig.rocky10.sh \
misc/reconfig.ubuntu22.04.intel.sh \
misc/reconfig.ubuntu22.04.mingw64.sh \
misc/reconfig.ubuntu22.04.nvidia.sh \
misc/reconfig.ubuntu22.04.sh \
misc/reconfig.ubuntu24.04.intel.sh \
misc/reconfig.ubuntu24.04.mingw64.sh \
misc/reconfig.ubuntu24.04.nvidia.sh \
misc/reconfig.ubuntu24.04.sh \
misc/reconfig.ubuntu26.04.c++17.sh \
misc/reconfig.ubuntu26.04.c++20.sh \
misc/reconfig.ubuntu26.04.c++23.sh \
misc/reconfig.ubuntu26.04.c++26.sh \
misc/reconfig.ubuntu26.04.gnuc++20.sh \
misc/reconfig.ubuntu26.04.gnuc++23.sh \
misc/reconfig.ubuntu26.04.gnuc++26.sh \
misc/reconfig.ubuntu26.04.intel.sh \
misc/reconfig.ubuntu26.04.mingw64.sh \
misc/reconfig.ubuntu26.04.nvidia.sh \
misc/reconfig.ubuntu26.04.sh \
misc/summarize_mplapack_tests.py \
mpblas/optimized/dd/cuda/dd_real_cuda.h \
mpfrc++/copying.lib mpfrc++/copying mpfrc++/makefile.org \
mpfrc++/org/copying mpfrc++/org/copying.lib mpfrc++/org/example.cpp \
mpfrc++/org/makefile mpfrc++/org/mpfrc++.tar.gz \
mpfrc++/org/mpreal.cpp mpfrc++/org/mpreal.h \
mplapack/test/eig/Cbak.in \
mplapack/test/eig/Cbal.in \
mplapack/test/eig/Cbb.in \
mplapack/test/eig/Cdmd.in \
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
mplapack/test/eig/Rdmd.in \
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
mplapack/test/lin/Ctest.in \
mplapack/test/lin/Ctest_rfp.in \
mplapack/test/lin/Rtest.in \
mplapack/test/lin/Rtest_rfp.in

uninstall-hook:
	-rm -f $(DESTDIR)$(infodir)/dir

distclean-local:
	rm -f $(top_builddir)/external/openblas/Makefile
