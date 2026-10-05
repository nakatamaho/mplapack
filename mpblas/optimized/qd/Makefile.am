include $(top_srcdir)/mplapack/reference/sources.am

SUBDIRS =
if ENABLE_QD_CUDA
SUBDIRS += cuda
endif

# Basename shadowing intentionally omits Raxpy.cpp, Rcopy.cpp, Rdot.cpp,
# and Rgemm.cpp from the reference BLAS list. The public dispatchers below
# replace them; *_ref.cpp and openmp/*_omp.cpp remain non-shadowing helpers.
MPLAPACK_QD_OPTIMIZED_SOURCES =

lib_LTLIBRARIES = libmplapack_qd_opt.la

libmplapack_qd_opt_la_CPPFLAGS = -I. -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) $(OPENMP_CXXFLAGS) -DMPLAPACK_BUILD_WITH_QD
# Rgemm_blocked_omp.cpp needs uncontracted quad-double arithmetic.
libmplapack_qd_opt_la_CXXFLAGS = $(DD_CXXFLAGS)
libmplapack_qd_opt_la_SOURCES = \
Rcopy.cpp Rcopy_ref.cpp openmp/Rcopy_omp.cpp \
Rdot.cpp  Rdot_ref.cpp  openmp/Rdot_omp.cpp \
Rgemm.cpp Rgemm_ref.cpp openmp/Rgemm_blocked_omp.cpp openmp/Rgemm_NN_omp.cpp openmp/Rgemm_NT_omp.cpp openmp/Rgemm_TN_omp.cpp openmp/Rgemm_TT_omp.cpp \
Raxpy.cpp Raxpy_ref.cpp openmp/Raxpy_omp.cpp \
../../reference/CRrot.cpp \
../../reference/CRscal.cpp \
../../reference/Caxpy.cpp \
../../reference/Ccopy.cpp \
../../reference/Cdotc.cpp \
../../reference/Cdotu.cpp \
../../reference/Cgbmv.cpp \
../../reference/Cgemm.cpp \
../../reference/Cgemmtr.cpp \
../../reference/Cgemv.cpp \
../../reference/Cgerc.cpp \
../../reference/Cgeru.cpp \
../../reference/Chbmv.cpp \
../../reference/Chemm.cpp \
../../reference/Chemv.cpp \
../../reference/Cher.cpp \
../../reference/Cher2.cpp \
../../reference/Cher2k.cpp \
../../reference/Cherk.cpp \
../../reference/Chpmv.cpp \
../../reference/Chpr.cpp \
../../reference/Chpr2.cpp \
../../reference/Crotg.cpp \
../../reference/Cscal.cpp \
../../reference/Cswap.cpp \
../../reference/Csymm.cpp \
../../reference/Csyr2k.cpp \
../../reference/Csyrk.cpp \
../../reference/Ctbmv.cpp \
../../reference/Ctbsv.cpp \
../../reference/Ctpmv.cpp \
../../reference/Ctpsv.cpp \
../../reference/Ctrmm.cpp \
../../reference/Ctrmv.cpp \
../../reference/Ctrsm.cpp \
../../reference/Ctrsv.cpp \
../../reference/Mlsame.cpp \
../../reference/Mxerbla.cpp \
../../reference/RCabs1.cpp \
../../reference/RCasum.cpp \
../../reference/RCnrm2.cpp \
../../reference/Rasum.cpp \
../../reference/Rgbmv.cpp \
../../reference/Rgemmtr.cpp \
../../reference/Rgemv.cpp \
../../reference/Rger.cpp \
../../reference/Rnrm2.cpp \
../../reference/Rrot.cpp \
../../reference/Rrotg.cpp \
../../reference/Rrotm.cpp \
../../reference/Rrotmg.cpp \
../../reference/Rsbmv.cpp \
../../reference/Rscal.cpp \
../../reference/Rspmv.cpp \
../../reference/Rspr.cpp \
../../reference/Rspr2.cpp \
../../reference/Rswap.cpp \
../../reference/Rsymm.cpp \
../../reference/Rsymv.cpp \
../../reference/Rsyr.cpp \
../../reference/Rsyr2.cpp \
../../reference/Rsyr2k.cpp \
../../reference/Rsyrk.cpp \
../../reference/Rtbmv.cpp \
../../reference/Rtbsv.cpp \
../../reference/Rtpmv.cpp \
../../reference/Rtpsv.cpp \
../../reference/Rtrmm.cpp \
../../reference/Rtrmv.cpp \
../../reference/Rtrsm.cpp \
../../reference/Rtrsv.cpp \
../../reference/iCamax.cpp \
../../reference/iRamax.cpp \
../../reference/mplapackinit.cpp \
$(MPLAPACK_OPT_REFERENCE_SOURCES) \
$(MPLAPACK_QD_OPTIMIZED_SOURCES)

libmplapack_qd_opt_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO) -Wl,-rpath,$(libdir)
libmplapack_qd_opt_la_LIBADD = $(QD_LIBADD) $(OPENMP_LIBS)
if IS_MACOS
libmplapack_qd_opt_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmplapack_qd_opt_la_LDFLAGS += -no-undefined
endif

# Optimized Rgemm (../common/Rgemm_*_common.h): the blocked engine is
# bit-identical to openmp/Rgemm_*_omp.cpp; Winograd and Ozaki are bounded.
check_PROGRAMS = Rgemm_opt_test Rgemm_gpu_host_test
Rgemm_opt_test_SOURCES = test/Rgemm_opt_test.cpp
Rgemm_opt_test_CPPFLAGS = $(libmplapack_qd_opt_la_CPPFLAGS) -DRGEMM_TEST_QD
Rgemm_opt_test_CXXFLAGS = $(DD_CXXFLAGS)
Rgemm_opt_test_LDFLAGS = -static $(QD_RPATH_LDFLAGS)
Rgemm_opt_test_LDADD = libmplapack_qd_opt.la $(libmplapack_qd_opt_la_LIBADD)
# The CUDA Rgemm code run on the host (no CUDA needed).
Rgemm_gpu_host_test_SOURCES = test/Rgemm_gpu_host_test.cpp cuda/Rgemm_bridge_cuda.cpp
Rgemm_gpu_host_test_CPPFLAGS = $(libmplapack_qd_opt_la_CPPFLAGS) -DRGEMM_TEST_QD -DMPLAPACK_GEMM_CUDA_HOST_EMULATION
Rgemm_gpu_host_test_CXXFLAGS = $(DD_CXXFLAGS)
Rgemm_gpu_host_test_LDFLAGS = -static $(QD_RPATH_LDFLAGS)
Rgemm_gpu_host_test_LDADD = libmplapack_qd_opt.la $(libmplapack_qd_opt_la_LIBADD)
TESTS = Rgemm_opt_test Rgemm_gpu_host_test
