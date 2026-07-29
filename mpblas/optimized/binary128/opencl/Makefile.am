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
../../../reference/CRrot.cpp \
../../../reference/CRscal.cpp \
../../../reference/Caxpy.cpp \
../../../reference/Ccopy.cpp \
../../../reference/Cdotc.cpp \
../../../reference/Cdotu.cpp \
../../../reference/Cgbmv.cpp \
../../../reference/Cgemm.cpp \
../../../reference/Cgemmtr.cpp \
../../../reference/Cgemv.cpp \
../../../reference/Cgerc.cpp \
../../../reference/Cgeru.cpp \
../../../reference/Chbmv.cpp \
../../../reference/Chemm.cpp \
../../../reference/Chemv.cpp \
../../../reference/Cher.cpp \
../../../reference/Cher2.cpp \
../../../reference/Cher2k.cpp \
../../../reference/Cherk.cpp \
../../../reference/Chpmv.cpp \
../../../reference/Chpr.cpp \
../../../reference/Chpr2.cpp \
../../../reference/Crotg.cpp \
../../../reference/Cscal.cpp \
../../../reference/Cswap.cpp \
../../../reference/Csymm.cpp \
../../../reference/Csyr2k.cpp \
../../../reference/Csyrk.cpp \
../../../reference/Ctbmv.cpp \
../../../reference/Ctbsv.cpp \
../../../reference/Ctpmv.cpp \
../../../reference/Ctpsv.cpp \
../../../reference/Ctrmm.cpp \
../../../reference/Ctrmv.cpp \
../../../reference/Ctrsm.cpp \
../../../reference/Ctrsv.cpp \
../../../reference/Mlsame.cpp \
../../../reference/Mxerbla.cpp \
../../../reference/RCabs1.cpp \
../../../reference/RCasum.cpp \
../../../reference/RCnrm2.cpp \
../../../reference/Rasum.cpp \
../../../reference/Rgbmv.cpp \
../../../reference/Rgemmtr.cpp \
../../../reference/Rgemv.cpp \
../../../reference/Rger.cpp \
../../../reference/Rnrm2.cpp \
../../../reference/Rrot.cpp \
../../../reference/Rrotg.cpp \
../../../reference/Rrotm.cpp \
../../../reference/Rrotmg.cpp \
../../../reference/Rsbmv.cpp \
../../../reference/Rscal.cpp \
../../../reference/Rspmv.cpp \
../../../reference/Rspr.cpp \
../../../reference/Rspr2.cpp \
../../../reference/Rswap.cpp \
../../../reference/Rsymm.cpp \
../../../reference/Rsymv.cpp \
../../../reference/Rsyr.cpp \
../../../reference/Rsyr2.cpp \
../../../reference/Rsyr2k.cpp \
../../../reference/Rsyrk.cpp \
../../../reference/Rtbmv.cpp \
../../../reference/Rtbsv.cpp \
../../../reference/Rtpmv.cpp \
../../../reference/Rtpsv.cpp \
../../../reference/Rtrmm.cpp \
../../../reference/Rtrmv.cpp \
../../../reference/Rtrsm.cpp \
../../../reference/Rtrsv.cpp \
../../../reference/iCamax.cpp \
../../../reference/iRamax.cpp \
../../../reference/mplapackinit.cpp \
$(MPLAPACK_ACCEL_REFERENCE_SOURCES) \
$(MPLAPACK_BINARY128_OPENCL_LAPACK_OPTIMIZED_SOURCES) \
Rgemm.cpp OclHelper.hpp cpp_kernel.h

libmplapack_binary128_opt_opencl_la_LDFLAGS = $(AM_LDFLAGS) $(OPENCL_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO) -Wl,-rpath,$(libdir)
libmplapack_binary128_opt_opencl_la_LIBADD = $(OPENCL_LIBS)
if BINARY128_USE_QUADMATH
libmplapack_binary128_opt_opencl_la_LIBADD += -lquadmath
endif
if IS_MACOS
libmplapack_binary128_opt_opencl_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmplapack_binary128_opt_opencl_la_LDFLAGS += -fopenmp -no-undefined
endif
