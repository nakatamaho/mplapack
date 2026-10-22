include $(top_srcdir)/mplapack/reference/accelerator-sources.am

# Layered basename shadowing:
# - cuda/Rgemm.cpp replaces ../Rgemm.cpp; it computes eligible calls on
#   the GPU and otherwise runs the same OpenMP code as ../Rgemm.cpp.
# - cuda/R{syrk,syr2k,symm,trmm,trsm,gemv,symv,trmv,ger,syr,syr2,scal,swap,
#   rot}.cpp replace ../../common/blas/R*.cpp (the GPU: blas3_bridge_cuda.cpp;
#   Level 1 and 2 only with MPLAPACK_MPFR_CUDA_L12=1); the CPU fallback is
#   ../../common/blas/R*_omp.cpp (lane-parallel OpenMP and the reference code
#   under the name R*_reference_cpu).
# - *_cuda.* files are helpers and never shadow.
MPLAPACK_MPFR_CUDA_LAPACK_OPTIMIZED_SOURCES =

# The CUDA sources are relocatable device code (Rgemm_rt_cuda.cu calls the
# cu_mpfr device functions of libmpc_cuda.a); Rgemm_dlink_cuda.lo is their
# device link, so that C++ programs can link the library normally.
NVCC_MPFR_CUDA_FLAGS = $(NVCC_GLIBC_FEATURE_FLAGS) -std=c++17 -rdc=true -fmad=false -I$(MPC_CUDA_INCLUDEDIR) -I$(srcdir) -I$(top_builddir)/include -I$(top_srcdir)/include $(MPFR_CUDA_GENCODE)

.cu.o:
	$(NVCC) $(NVCC_MPFR_CUDA_FLAGS) -o $@ -c $< $(CUDA_CFLAGS)

.cu.lo:
	$(top_srcdir)/misc/cudalt.py $@ $(NVCC) $(NVCC_MPFR_CUDA_FLAGS) --compiler-options=\"$(CFLAGS) $(DEFAULT_INCLUDES) $(CPPFLAGS) \" -c $<

Rgemm_dlink_cuda.lo: Rgemm_device_cuda.lo Rgemm_rt_cuda.lo blas3_device_cuda.lo Cblas_mpfr_device_cuda.lo
	$(top_srcdir)/misc/cudalt.py $@ $(NVCC) $(MPFR_CUDA_GENCODE) -dlink .libs/Rgemm_device_cuda.o .libs/Rgemm_rt_cuda.o .libs/blas3_device_cuda.o .libs/Cblas_mpfr_device_cuda.o $(MPC_CUDA_LIBRARY)

# Must match the architecture libmpc_cuda.a was built for (--with-mpfr-cuda-arch).
MPFR_CUDA_GENCODE = -gencode arch=compute_$(MPFR_CUDA_ARCH),code=sm_$(MPFR_CUDA_ARCH) -m64 -O3
NVCC_GLIBC_FEATURE_FLAGS = -U_GNU_SOURCE -D_XOPEN_SOURCE=700

MPLAPACK_MPFR_CUDA_CPPFLAGS = $(GMPFRXX_MKII_CPPFLAGS) -I. -I$(srcdir) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(MPC_CUDA_INCLUDEDIR) $(OPENMP_CXXFLAGS) -DMPLAPACK_BUILD_WITH_MPFR
AM_CPPFLAGS = $(MPLAPACK_MPFR_CUDA_CPPFLAGS)

lib_LTLIBRARIES = libmplapack_mpfr_opt_cuda.la
libmplapack_mpfr_opt_cuda_la_CPPFLAGS = $(MPLAPACK_MPFR_CUDA_CPPFLAGS)
libmplapack_mpfr_opt_cuda_la_SOURCES = \
../Rcopy.cpp ../Rcopy_ref.cpp ../openmp/Rcopy_omp.cpp \
../Rdot.cpp  ../Rdot_ref.cpp  ../openmp/Rdot_omp.cpp \
../Rgemm_ref.cpp ../openmp/Rgemm_NN_omp.cpp ../openmp/Rgemm_NT_omp.cpp ../openmp/Rgemm_TN_omp.cpp ../openmp/Rgemm_TT_omp.cpp \
../Raxpy.cpp ../Raxpy_ref.cpp ../openmp/Raxpy_omp.cpp \
../../common/blas/Rsyrk_omp.cpp \
../../common/blas/Rsyr2k_omp.cpp \
../../common/blas/Rsymm_omp.cpp \
../../common/blas/Rtrmm_omp.cpp \
../../common/blas/Rtrsm_omp.cpp \
../../common/blas/Rgemv_omp.cpp Rgemv.cpp \
../../common/blas/Rsymv_omp.cpp Rsymv.cpp \
../../common/blas/Rtrmv_omp.cpp Rtrmv.cpp \
../../common/blas/Rger_omp.cpp Rger.cpp \
../../common/blas/Rsyr_omp.cpp Rsyr.cpp \
../../common/blas/Rsyr2_omp.cpp Rsyr2.cpp \
../../common/blas/Rscal_omp.cpp Rscal.cpp \
../../common/blas/Rswap_omp.cpp Rswap.cpp \
../../common/blas/Rrot_omp.cpp Rrot.cpp \
../../common/blas_cuda/Rgemmtr.cpp ../../common/blas/Rgemmtr_omp.cpp \
../../common/blas_cuda/Cgemm.cpp ../../common/blas/Cgemm_omp.cpp \
../../common/blas_cuda/Csymm.cpp ../../common/blas/Csymm_omp.cpp \
../../common/blas_cuda/Chemm.cpp ../../common/blas/Chemm_omp.cpp \
../../common/blas_cuda/Ctrmm.cpp ../../common/blas/Ctrmm_omp.cpp \
../../common/blas_cuda/Ctrsm.cpp ../../common/blas/Ctrsm_omp.cpp \
../../common/blas_cuda/Cgemv.cpp ../../common/blas/Cgemv_omp.cpp \
../../common/blas_cuda/Cgeru.cpp ../../common/blas/Cgeru_omp.cpp \
../../common/blas_cuda/Cgerc.cpp ../../common/blas/Cgerc_omp.cpp \
../../common/blas_cuda/Caxpy.cpp ../../common/blas/Caxpy_omp.cpp \
../../common/blas_cuda/Ccopy.cpp ../../common/blas/Ccopy_omp.cpp \
../../common/blas_cuda/Cscal.cpp ../../common/blas/Cscal_omp.cpp \
../../common/blas_cuda/CRscal.cpp ../../common/blas/CRscal_omp.cpp \
../../common/blas_cuda/Cswap.cpp ../../common/blas/Cswap_omp.cpp \
../../common/blas_cuda/CRrot.cpp ../../common/blas/CRrot_omp.cpp \
../../common/blas_cuda/Csyrk.cpp ../../common/blas/Csyrk_omp.cpp \
../../common/blas_cuda/Cherk.cpp ../../common/blas/Cherk_omp.cpp \
../../common/blas_cuda/Csyr2k.cpp ../../common/blas/Csyr2k_omp.cpp \
../../common/blas_cuda/Cher2k.cpp ../../common/blas/Cher2k_omp.cpp \
../../common/blas_cuda/Cgemmtr.cpp ../../common/blas/Cgemmtr_omp.cpp \
../../common/blas/Rtrsv.cpp ../../common/blas/Rtrsv_omp.cpp \
../../common/blas/Ctrsv.cpp ../../common/blas/Ctrsv_omp.cpp \
../../common/blas_cuda/Cher.cpp ../../common/blas/Cher_omp.cpp \
../../common/blas_cuda/Cher2.cpp ../../common/blas/Cher2_omp.cpp \
../../../reference/Cdotc.cpp \
../../../reference/Cdotu.cpp \
../../../reference/Cgbmv.cpp \
../../../reference/Chbmv.cpp \
../../../reference/Chemv.cpp \
../../../reference/Chpmv.cpp \
../../../reference/Chpr.cpp \
../../../reference/Chpr2.cpp \
../../../reference/Crotg.cpp \
../../../reference/Ctbmv.cpp \
../../../reference/Ctbsv.cpp \
../../../reference/Ctpmv.cpp \
../../../reference/Ctpsv.cpp \
../../../reference/Ctrmv.cpp \
../../../reference/Mlsame.cpp \
../../../reference/Mxerbla.cpp \
../../../reference/RCabs1.cpp \
../../../reference/RCasum.cpp \
../../../reference/RCnrm2.cpp \
../../../reference/Rasum.cpp \
../../../reference/Rgbmv.cpp \
../../../reference/Rnrm2.cpp \
../../../reference/Rrotg.cpp \
../../../reference/Rrotm.cpp \
../../../reference/Rrotmg.cpp \
../../../reference/Rsbmv.cpp \
../../../reference/Rspmv.cpp \
../../../reference/Rspr.cpp \
../../../reference/Rspr2.cpp \
../../../reference/Rtbmv.cpp \
../../../reference/Rtbsv.cpp \
../../../reference/Rtpmv.cpp \
../../../reference/Rtpsv.cpp \
../../../reference/iCamax.cpp \
../../../reference/iRamax.cpp \
../../../reference/mplapackinit.cpp \
$(MPLAPACK_ACCEL_REFERENCE_SOURCES) \
../../../../mplapack/optimized/common/Rlaswp.cpp \
../../../../mplapack/optimized/common/Claswp.cpp \
$(MPLAPACK_MPFR_CUDA_LAPACK_OPTIMIZED_SOURCES) \
Rgemm.cpp Rgemm_bridge_cuda.cpp Rgemm_device_cuda.cu Rgemm_kernel_cuda.h \
Rgemm_rt_cuda.cu Rgemm_rt_cuda.h Rgemm_winograd_cuda.h mpfr_pack_cuda.h \
Rsyrk.cpp Rsyr2k.cpp Rsymm.cpp Rtrmm.cpp Rtrsm.cpp \
blas3_bridge_cuda.cpp blas3_device_cuda.cu blas3_kernel_cuda.h \
Cblas_mpfr_bridge_cuda.cpp Cblas_mpfr_device_cuda.cu Cblas_mpfr_dev_types.h

libmplapack_mpfr_opt_cuda_la_LDFLAGS = $(AM_LDFLAGS) $(CUDA_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO) -Wl,-rpath,$(libdir)
libmplapack_mpfr_opt_cuda_la_LIBADD = Rgemm_dlink_cuda.lo $(MPC_CUDA_LIBRARY) $(CUDA_LIBS) -lcudart $(MPLAPACK_MPFR_LIBS) $(OPENMP_LIBS)
EXTRA_libmplapack_mpfr_opt_cuda_la_DEPENDENCIES = Rgemm_dlink_cuda.lo
CLEANFILES = Rgemm_dlink_cuda.lo .libs/Rgemm_dlink_cuda.o Rgemm_dlink_cuda.o
if IS_MACOS
libmplapack_mpfr_opt_cuda_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmplapack_mpfr_opt_cuda_la_LDFLAGS += -fopenmp -no-undefined
endif

EXTRA_DIST = README.md Rgemm_host_cuda.cpp test/Rgemm_cuda_test.cpp test/Rgemm_winograd_cuda_test.cpp test/blas3_cuda_test.cpp
