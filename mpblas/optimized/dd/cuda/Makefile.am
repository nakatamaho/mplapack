include $(top_srcdir)/mplapack/reference/accelerator-sources.am

# Layered basename shadowing:
# - cuda/Rgemm.cpp replaces ../Rgemm.cpp; it computes eligible calls on the
#   GPU (Rgemm_bridge_cuda.cpp, Rgemm_gpu_cuda.cu) and otherwise runs the CPU
#   code of ../Rgemm.cpp.
# - ../../common/blas_cuda/R*.cpp replace ../../common/blas/R*.cpp: Level
#   1/2/3 BLAS on the GPU (Rblas_bridge_cuda.cpp), else the CPU code of
#   ../../common/blas/R*_omp.cpp.
# - Rsyrk_legacy_cuda.cu holds the 2010-2011 Rsyrk kernels
#   (MPLAPACK_DD_CUDA_LEGACY=1).
# - *_cuda.* files and split kernels are helpers and never shadow.
MPLAPACK_DD_CUDA_LAPACK_OPTIMIZED_SOURCES =

# -fmad=false: the double-double device arithmetic must not be contracted.
NVCC_DD_CUDA_FLAGS = -std=c++17 -fmad=false -I$(srcdir) -I$(QD_INCLUDEDIR)

.cu.o:
	$(NVCC) $(NVCC_GLIBC_FEATURE_FLAGS) $(NVCC_DD_CUDA_FLAGS) $(HOPPER_FLAGS) -o $@ -c $< $(CUDA_CFLAGS)

.cu.lo:
	$(top_srcdir)/misc/cudalt.py $@ $(NVCC) $(NVCC_GLIBC_FEATURE_FLAGS) $(NVCC_DD_CUDA_FLAGS) --compiler-options=\"$(CFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \" $(AMPERE_FLAGS) $(HOPPER_FLAGS) -c $<

AMPERE_FLAGS = -gencode arch=compute_80,code=sm_80 -m64 -O3 --ptxas-options=-v
HOPPER_FLAGS = -gencode arch=compute_90,code=sm_90 -m64 -O3 --ptxas-options=-v
NVCC_GLIBC_FEATURE_FLAGS = -U_GNU_SOURCE -D_XOPEN_SOURCE=700

MPLAPACK_DD_CUDA_CPPFLAGS = -I. -I.. -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) $(OPENMP_CXXFLAGS) -DMPLAPACK_BUILD_WITH_DD
AM_CPPFLAGS = $(MPLAPACK_DD_CUDA_CPPFLAGS)

lib_LTLIBRARIES = libmplapack_dd_opt_cuda.la
libmplapack_dd_opt_cuda_la_CPPFLAGS = $(MPLAPACK_DD_CUDA_CPPFLAGS)
libmplapack_dd_opt_cuda_la_CXXFLAGS = $(DD_CXXFLAGS)
libmplapack_dd_opt_cuda_la_SOURCES = \
../Rcopy.cpp ../Rcopy_ref.cpp ../openmp/Rcopy_omp.cpp \
../Rdot.cpp  ../Rdot_ref.cpp  ../openmp/Rdot_omp.cpp \
../Rgemm_ref.cpp ../openmp/Rgemm_blocked_omp.cpp ../openmp/Rgemm_NN_omp.cpp ../openmp/Rgemm_NT_omp.cpp ../openmp/Rgemm_TN_omp.cpp ../openmp/Rgemm_TT_omp.cpp \
../Raxpy.cpp ../Raxpy_ref.cpp ../openmp/Raxpy_omp.cpp \
../../common/blas/Rsyrk_omp.cpp ../../common/blas_cuda/Rsyrk.cpp \
../../common/blas/Rsyr2k_omp.cpp ../../common/blas_cuda/Rsyr2k.cpp \
../../common/blas/Rsymm_omp.cpp ../../common/blas_cuda/Rsymm.cpp \
../../common/blas/Rtrmm_omp.cpp ../../common/blas_cuda/Rtrmm.cpp \
../../common/blas/Rtrsm_omp.cpp ../../common/blas_cuda/Rtrsm.cpp \
../../common/blas/Rgemv_omp.cpp ../../common/blas_cuda/Rgemv.cpp \
../../common/blas/Rsymv_omp.cpp ../../common/blas_cuda/Rsymv.cpp \
../../common/blas/Rtrmv_omp.cpp ../../common/blas_cuda/Rtrmv.cpp \
../../common/blas/Rger_omp.cpp ../../common/blas_cuda/Rger.cpp \
../../common/blas/Rsyr_omp.cpp ../../common/blas_cuda/Rsyr.cpp \
../../common/blas/Rsyr2_omp.cpp ../../common/blas_cuda/Rsyr2.cpp \
../../common/blas/Rscal_omp.cpp ../../common/blas_cuda/Rscal.cpp \
../../common/blas/Rswap_omp.cpp ../../common/blas_cuda/Rswap.cpp \
../../common/blas/Rrot_omp.cpp ../../common/blas_cuda/Rrot.cpp \
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
../../../reference/Cdotc.cpp \
../../../reference/Cdotu.cpp \
../../../reference/Cgbmv.cpp \
../../../reference/Chbmv.cpp \
../../../reference/Chemv.cpp \
../../../reference/Cher.cpp \
../../../reference/Cher2.cpp \
../../../reference/Chpmv.cpp \
../../../reference/Chpr.cpp \
../../../reference/Chpr2.cpp \
../../../reference/Crotg.cpp \
../../../reference/Ctbmv.cpp \
../../../reference/Ctbsv.cpp \
../../../reference/Ctpmv.cpp \
../../../reference/Ctpsv.cpp \
../../../reference/Ctrmv.cpp \
../../../reference/Ctrsv.cpp \
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
../../../reference/Rtrsv.cpp \
../../../reference/iCamax.cpp \
../../../reference/iRamax.cpp \
../../../reference/mplapackinit.cpp \
$(MPLAPACK_ACCEL_REFERENCE_SOURCES) \
../../../../mplapack/optimized/common/Rlaswp.cpp \
../../../../mplapack/optimized/common/Claswp.cpp \
$(MPLAPACK_DD_CUDA_LAPACK_OPTIMIZED_SOURCES) \
Rgemm.cpp Rgemm_bridge_cuda.cpp Rgemm_gpu_cuda.cu ../../common/blas_cuda/Rblas_bridge_cuda.cpp \
Rgemm_cuda.cu Rgemm_fermi.cu Rgemm_tesla.cu \
Rsyrk_legacy_cuda.cu Rsyrk_cuda.cu dd_real_cuda.h

libmplapack_dd_opt_cuda_la_LDFLAGS = $(AM_LDFLAGS) $(CUDA_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO) -Wl,-rpath,$(libdir)
libmplapack_dd_opt_cuda_la_LIBADD = $(CUDA_LIBS) -lcublas -lcudart $(QD_LIBADD) $(OPENMP_LIBS)
if IS_MACOS
libmplapack_dd_opt_cuda_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmplapack_dd_opt_cuda_la_LDFLAGS += -fopenmp -no-undefined
endif

EXTRA_DIST = \
Rgemm_fermi_NN_0.cu \
Rgemm_fermi_NN_p.cu \
Rgemm_fermi_NT_0.cu \
Rgemm_fermi_NT_p.cu \
Rgemm_fermi_TN_0.cu \
Rgemm_fermi_TN_p.cu \
Rgemm_fermi_TT_0.cu \
Rgemm_fermi_TT_p.cu \
Rgemm_tesla_NN_0.cu \
Rgemm_tesla_NN_p.cu \
Rgemm_tesla_NT_0.cu \
Rgemm_tesla_NT_p.cu \
Rgemm_tesla_TN_0.cu \
Rgemm_tesla_TN_p.cu \
Rgemm_tesla_TT_0.cu \
Rgemm_tesla_TT_p.cu \
Rsyrk_NL_0.cu \
Rsyrk_NL_p.cu \
Rsyrk_NU_0.cu \
Rsyrk_NU_p.cu \
Rsyrk_TL_0.cu \
Rsyrk_TL_p.cu \
Rsyrk_TU_0.cu \
Rsyrk_TU_p.cu
