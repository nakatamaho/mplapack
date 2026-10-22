include $(top_srcdir)/mplapack/reference/accelerator-sources.am

# opencl/Rgemm.cpp defines public Rgemm and shadows ../Rgemm.cpp.
# Rgemm_ref.cpp remains as the CPU fallback helper.
MPLAPACK_BINARY128_OPENCL_LAPACK_OPTIMIZED_SOURCES =

MPLAPACK_BINARY128_OPENCL_CPPFLAGS = -I. -I.. -I$(top_srcdir)/include $(OPENMP_CXXFLAGS) -DMPLAPACK_BUILD_WITH_BINARY128
AM_CPPFLAGS = $(MPLAPACK_BINARY128_OPENCL_CPPFLAGS)

lib_LTLIBRARIES = libmplapack_binary128_opt_opencl.la
libmplapack_binary128_opt_opencl_la_CPPFLAGS = $(MPLAPACK_BINARY128_OPENCL_CPPFLAGS)
libmplapack_binary128_opt_opencl_la_CXXFLAGS = $(OPENCL_CFLAGS)
libmplapack_binary128_opt_opencl_la_SOURCES = \
../Rcopy.cpp ../Rcopy_ref.cpp ../openmp/Rcopy_omp.cpp \
../Rdot.cpp  ../Rdot_ref.cpp  ../openmp/Rdot_omp.cpp \
../Rgemm_ref.cpp ../openmp/Rgemm_NN_omp.cpp ../openmp/Rgemm_NT_omp.cpp ../openmp/Rgemm_TN_omp.cpp ../openmp/Rgemm_TT_omp.cpp \
../Raxpy.cpp ../Raxpy_ref.cpp ../openmp/Raxpy_omp.cpp \
../../common/blas/Rsyrk.cpp ../../common/blas/Rsyrk_omp.cpp \
../../common/blas/Rsyr2k.cpp ../../common/blas/Rsyr2k_omp.cpp \
../../common/blas/Rsymm.cpp ../../common/blas/Rsymm_omp.cpp \
../../common/blas/Rtrmm.cpp ../../common/blas/Rtrmm_omp.cpp \
../../common/blas/Rtrsm.cpp ../../common/blas/Rtrsm_omp.cpp \
../../common/blas/Rgemv.cpp ../../common/blas/Rgemv_omp.cpp \
../../common/blas/Rsymv.cpp ../../common/blas/Rsymv_omp.cpp \
../../common/blas/Rtrmv.cpp ../../common/blas/Rtrmv_omp.cpp \
../../common/blas/Rger.cpp ../../common/blas/Rger_omp.cpp \
../../common/blas/Rsyr.cpp ../../common/blas/Rsyr_omp.cpp \
../../common/blas/Rsyr2.cpp ../../common/blas/Rsyr2_omp.cpp \
../../common/blas/Rscal.cpp ../../common/blas/Rscal_omp.cpp \
../../common/blas/Rswap.cpp ../../common/blas/Rswap_omp.cpp \
../../common/blas/Rrot.cpp ../../common/blas/Rrot_omp.cpp \
../../common/blas/Rgemmtr.cpp ../../common/blas/Rgemmtr_omp.cpp \
../../common/blas/Cgemm.cpp ../../common/blas/Cgemm_omp.cpp \
../../common/blas/Csymm.cpp ../../common/blas/Csymm_omp.cpp \
../../common/blas/Chemm.cpp ../../common/blas/Chemm_omp.cpp \
../../common/blas/Ctrmm.cpp ../../common/blas/Ctrmm_omp.cpp \
../../common/blas/Ctrsm.cpp ../../common/blas/Ctrsm_omp.cpp \
../../common/blas/Cgemv.cpp ../../common/blas/Cgemv_omp.cpp \
../../common/blas/Cgeru.cpp ../../common/blas/Cgeru_omp.cpp \
../../common/blas/Cgerc.cpp ../../common/blas/Cgerc_omp.cpp \
../../common/blas/Caxpy.cpp ../../common/blas/Caxpy_omp.cpp \
../../common/blas/Ccopy.cpp ../../common/blas/Ccopy_omp.cpp \
../../common/blas/Cscal.cpp ../../common/blas/Cscal_omp.cpp \
../../common/blas/CRscal.cpp ../../common/blas/CRscal_omp.cpp \
../../common/blas/Cswap.cpp ../../common/blas/Cswap_omp.cpp \
../../common/blas/CRrot.cpp ../../common/blas/CRrot_omp.cpp \
../../common/blas/Csyrk.cpp ../../common/blas/Csyrk_omp.cpp \
../../common/blas/Cherk.cpp ../../common/blas/Cherk_omp.cpp \
../../common/blas/Csyr2k.cpp ../../common/blas/Csyr2k_omp.cpp \
../../common/blas/Cher2k.cpp ../../common/blas/Cher2k_omp.cpp \
../../common/blas/Cgemmtr.cpp ../../common/blas/Cgemmtr_omp.cpp \
../../common/blas/Rtrsv.cpp ../../common/blas/Rtrsv_omp.cpp \
../../common/blas/Ctrsv.cpp ../../common/blas/Ctrsv_omp.cpp \
../../common/blas/Cher.cpp ../../common/blas/Cher_omp.cpp \
../../common/blas/Cher2.cpp ../../common/blas/Cher2_omp.cpp \
../../../reference/Cdotc.cpp \
../../../reference/Cdotu.cpp \
../../common/blas/Cgbmv.cpp ../../common/blas/Cgbmv_omp.cpp \
../../common/blas/Chbmv.cpp ../../common/blas/Chbmv_omp.cpp \
../../common/blas/Chemv.cpp ../../common/blas/Chemv_omp.cpp \
../../common/blas/Chpmv.cpp ../../common/blas/Chpmv_omp.cpp \
../../common/blas/Chpr.cpp ../../common/blas/Chpr_omp.cpp \
../../common/blas/Chpr2.cpp ../../common/blas/Chpr2_omp.cpp \
../../../reference/Crotg.cpp \
../../common/blas/Ctbmv.cpp ../../common/blas/Ctbmv_omp.cpp \
../../common/blas/Ctbsv.cpp ../../common/blas/Ctbsv_omp.cpp \
../../common/blas/Ctpmv.cpp ../../common/blas/Ctpmv_omp.cpp \
../../common/blas/Ctpsv.cpp ../../common/blas/Ctpsv_omp.cpp \
../../common/blas/Ctrmv.cpp ../../common/blas/Ctrmv_omp.cpp \
../../../reference/Mlsame.cpp \
../../../reference/Mxerbla.cpp \
../../../reference/RCabs1.cpp \
../../../reference/RCasum.cpp \
../../../reference/RCnrm2.cpp \
../../../reference/Rasum.cpp \
../../common/blas/Rgbmv.cpp ../../common/blas/Rgbmv_omp.cpp \
../../../reference/Rnrm2.cpp \
../../../reference/Rrotg.cpp \
../../../reference/Rrotm.cpp \
../../../reference/Rrotmg.cpp \
../../common/blas/Rsbmv.cpp ../../common/blas/Rsbmv_omp.cpp \
../../common/blas/Rspmv.cpp ../../common/blas/Rspmv_omp.cpp \
../../common/blas/Rspr.cpp ../../common/blas/Rspr_omp.cpp \
../../common/blas/Rspr2.cpp ../../common/blas/Rspr2_omp.cpp \
../../common/blas/Rtbmv.cpp ../../common/blas/Rtbmv_omp.cpp \
../../common/blas/Rtbsv.cpp ../../common/blas/Rtbsv_omp.cpp \
../../common/blas/Rtpmv.cpp ../../common/blas/Rtpmv_omp.cpp \
../../common/blas/Rtpsv.cpp ../../common/blas/Rtpsv_omp.cpp \
../../../reference/iCamax.cpp \
../../../reference/iRamax.cpp \
../../../reference/mplapackinit.cpp \
$(MPLAPACK_ACCEL_REFERENCE_SOURCES) \
../../../../mplapack/optimized/common/Rlaswp.cpp \
../../../../mplapack/optimized/common/Claswp.cpp \
$(MPLAPACK_BINARY128_OPENCL_LAPACK_OPTIMIZED_SOURCES) \
Rgemm.cpp OclHelper.hpp cpp_kernel.h

libmplapack_binary128_opt_opencl_la_LDFLAGS = $(AM_LDFLAGS) $(OPENCL_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO) -Wl,-rpath,$(libdir)
libmplapack_binary128_opt_opencl_la_LIBADD = $(OPENCL_LIBS) $(OPENMP_LIBS)
if BINARY128_USE_QUADMATH
libmplapack_binary128_opt_opencl_la_LIBADD += -lquadmath
endif
if IS_MACOS
libmplapack_binary128_opt_opencl_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmplapack_binary128_opt_opencl_la_LDFLAGS += -no-undefined
endif
