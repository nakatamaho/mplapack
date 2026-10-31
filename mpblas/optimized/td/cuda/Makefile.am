include $(top_srcdir)/mplapack/reference/accelerator-sources.am

# Layered basename shadowing:
# - cuda/Rgemm.cpp replaces ../Rgemm.cpp; it computes eligible calls on the
#   GPU (Rgemm_bridge_cuda.cpp, Rgemm_gpu_cuda.cu) and otherwise runs the CPU
#   code of ../Rgemm.cpp.
# - ../../common/blas_cuda/R*.cpp replace ../../common/blas/R*.cpp: Level
#   1/2/3 BLAS on the GPU (Rblas_bridge_cuda.cpp), else the CPU code of
#   ../../common/blas/R*_omp.cpp.
# - *_cuda.* files are helpers and never shadow.
MPLAPACK_TD_CUDA_LAPACK_OPTIMIZED_SOURCES =

# -fmad=false: the triple-double device arithmetic must not be contracted.
NVCC_TD_CUDA_FLAGS = -std=c++17 -fmad=false -I$(srcdir) -I$(top_builddir)/include -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)

.cu.o:
	$(NVCC) $(NVCC_GLIBC_FEATURE_FLAGS) $(NVCC_TD_CUDA_FLAGS) $(AMPERE_FLAGS) $(HOPPER_FLAGS) -o $@ -c $< $(CUDA_CFLAGS)

.cu.lo:
	$(top_srcdir)/misc/cudalt.py $@ $(NVCC) $(NVCC_GLIBC_FEATURE_FLAGS) $(NVCC_TD_CUDA_FLAGS) --compiler-options=\"$(CFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CPPFLAGS) \" $(AMPERE_FLAGS) $(HOPPER_FLAGS) -c $<

Rpotrf_gpu_cuda.lo Rpotrf_gpu_cuda.o Rgetrf_gpu_cuda.lo Rgetrf_gpu_cuda.o Rsolver_gpu_cuda.lo Rsolver_gpu_cuda.o: NVCC_TD_CUDA_FLAGS += -Xcompiler=-ffp-contract=off

# Header dependencies of the .cu objects, written by nvcc through
# misc/cudalt.py (automake does not track the .cu.lo suffix rule).
-include $(wildcard *.lo.d)

AMPERE_FLAGS = -gencode arch=compute_80,code=sm_80 -m64 -O3 --ptxas-options=-v
HOPPER_FLAGS = -gencode arch=compute_90,code=sm_90 -m64 -O3 --ptxas-options=-v
NVCC_GLIBC_FEATURE_FLAGS = -U_GNU_SOURCE -D_XOPEN_SOURCE=700

MPLAPACK_TD_CUDA_CPPFLAGS = -I. -I.. -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) $(OPENMP_CXXFLAGS) -DMPLAPACK_BUILD_WITH_TD -DMPLAPACK_ACCELERATOR_LIBRARY=1
AM_CPPFLAGS = $(MPLAPACK_TD_CUDA_CPPFLAGS)

lib_LTLIBRARIES = libmplapack_td_opt_cuda.la
libmplapack_td_opt_cuda_la_CPPFLAGS = $(MPLAPACK_TD_CUDA_CPPFLAGS)
libmplapack_td_opt_cuda_la_CXXFLAGS = $(DD_CXXFLAGS)
libmplapack_td_opt_cuda_la_SOURCES = \
../Rcopy.cpp ../Rcopy_ref.cpp ../openmp/Rcopy_omp.cpp \
../Rdot.cpp  ../Rdot_ref.cpp  ../openmp/Rdot_omp.cpp \
../../common/blas/Rdot_par.cpp \
../Rgemm_ref.cpp ../../common/blas/Rgemm_blocked_omp.cpp ../openmp/Rgemm_NN_omp.cpp ../openmp/Rgemm_NT_omp.cpp ../openmp/Rgemm_TN_omp.cpp ../openmp/Rgemm_TT_omp.cpp \
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
../../common/blas/Rgemmtr_omp.cpp ../../common/blas_cuda/Rgemmtr.cpp \
../../common/blas/Cgemm_omp.cpp ../../common/blas_cuda/Cgemm.cpp \
../../common/blas/Csymm_omp.cpp ../../common/blas_cuda/Csymm.cpp \
../../common/blas/Chemm_omp.cpp ../../common/blas_cuda/Chemm.cpp \
../../common/blas/Ctrmm_omp.cpp ../../common/blas_cuda/Ctrmm.cpp \
../../common/blas/Ctrsm_omp.cpp ../../common/blas_cuda/Ctrsm.cpp \
../../common/blas/Cgemv_omp.cpp ../../common/blas_cuda/Cgemv.cpp \
../../common/blas/Cgeru_omp.cpp ../../common/blas_cuda/Cgeru.cpp \
../../common/blas/Cgerc_omp.cpp ../../common/blas_cuda/Cgerc.cpp \
../../common/blas/Caxpy_omp.cpp ../../common/blas_cuda/Caxpy.cpp \
../../common/blas/Ccopy_omp.cpp ../../common/blas_cuda/Ccopy.cpp \
../../common/blas/Cscal_omp.cpp ../../common/blas_cuda/Cscal.cpp \
../../common/blas/CRscal_omp.cpp ../../common/blas_cuda/CRscal.cpp \
../../common/blas/Cswap_omp.cpp ../../common/blas_cuda/Cswap.cpp \
../../common/blas/CRrot_omp.cpp ../../common/blas_cuda/CRrot.cpp \
../../common/blas/Csyrk_omp.cpp ../../common/blas_cuda/Csyrk.cpp \
../../common/blas/Cherk_omp.cpp ../../common/blas_cuda/Cherk.cpp \
../../common/blas/Csyr2k_omp.cpp ../../common/blas_cuda/Csyr2k.cpp \
../../common/blas/Cher2k_omp.cpp ../../common/blas_cuda/Cher2k.cpp \
../../common/blas/Cgemmtr_omp.cpp ../../common/blas_cuda/Cgemmtr.cpp \
../../common/blas/Rtrsv.cpp ../../common/blas/Rtrsv_omp.cpp \
../../common/blas/Ctrsv.cpp ../../common/blas/Ctrsv_omp.cpp \
../../common/blas/Cher_omp.cpp ../../common/blas_cuda/Cher.cpp \
../../common/blas/Cher2_omp.cpp ../../common/blas_cuda/Cher2.cpp \
../../common/blas/Cdotc.cpp ../../common/blas/Cdotc_omp.cpp \
../../common/blas/Cdotu.cpp ../../common/blas/Cdotu_omp.cpp \
../../common/blas_cuda/Cgbmv.cpp ../../common/blas/Cgbmv_omp.cpp \
../../common/blas_cuda/Chbmv.cpp ../../common/blas/Chbmv_omp.cpp \
../../common/blas_cuda/Chemv.cpp ../../common/blas/Chemv_omp.cpp \
../../common/blas_cuda/Chpmv.cpp ../../common/blas/Chpmv_omp.cpp \
../../common/blas_cuda/Chpr.cpp ../../common/blas/Chpr_omp.cpp \
../../common/blas_cuda/Chpr2.cpp ../../common/blas/Chpr2_omp.cpp \
../../../reference/Crotg.cpp \
../../common/blas_cuda/Ctbmv.cpp ../../common/blas/Ctbmv_omp.cpp \
../../common/blas/Ctbsv.cpp ../../common/blas/Ctbsv_omp.cpp \
../../common/blas_cuda/Ctpmv.cpp ../../common/blas/Ctpmv_omp.cpp \
../../common/blas/Ctpsv.cpp ../../common/blas/Ctpsv_omp.cpp \
../../common/blas_cuda/Ctrmv.cpp ../../common/blas/Ctrmv_omp.cpp \
../../../reference/Mlsame.cpp \
../../../reference/Mxerbla.cpp \
../../../reference/RCabs1.cpp \
../../common/blas/RCasum.cpp ../../common/blas/RCasum_omp.cpp \
../../common/blas/RCnrm2.cpp ../../common/blas/RCnrm2_omp.cpp \
../../common/blas/Rasum.cpp ../../common/blas/Rasum_omp.cpp \
../../common/blas_cuda/Rgbmv.cpp ../../common/blas/Rgbmv_omp.cpp \
../../common/blas/Rnrm2.cpp ../../common/blas/Rnrm2_omp.cpp \
../../../reference/Rrotg.cpp \
../../../reference/Rrotm.cpp \
../../../reference/Rrotmg.cpp \
../../common/blas_cuda/Rsbmv.cpp ../../common/blas/Rsbmv_omp.cpp \
../../common/blas_cuda/Rspmv.cpp ../../common/blas/Rspmv_omp.cpp \
../../common/blas_cuda/Rspr.cpp ../../common/blas/Rspr_omp.cpp \
../../common/blas_cuda/Rspr2.cpp ../../common/blas/Rspr2_omp.cpp \
../../common/blas_cuda/Rtbmv.cpp ../../common/blas/Rtbmv_omp.cpp \
../../common/blas/Rtbsv.cpp ../../common/blas/Rtbsv_omp.cpp \
../../common/blas_cuda/Rtpmv.cpp ../../common/blas/Rtpmv_omp.cpp \
../../common/blas/Rtpsv.cpp ../../common/blas/Rtpsv_omp.cpp \
../../common/blas/iCamax.cpp ../../common/blas/iCamax_omp.cpp \
../../common/blas/iRamax.cpp ../../common/blas/iRamax_omp.cpp \
../../../reference/mplapackinit.cpp \
$(MPLAPACK_ACCEL_REFERENCE_SOURCES) \
../../../../mplapack/optimized/common/Rlaswp.cpp \
../../../../mplapack/optimized/common/Claswp.cpp \
Rgetrs.cpp Rpotrs.cpp Rgesv.cpp Rposv.cpp Cgetrs.cpp Cpotrs.cpp Cgesv.cpp Cposv.cpp Rsolver_gpu_cuda.cu \
Rgetrf.cpp Cgetrf.cpp Rgetrf_gpu_cuda.cu \
../../../../mplapack/optimized/common/Rgetrf2.cpp \
../../../../mplapack/optimized/common/Cgetrf2.cpp \
Rpotrf.cpp Cpotrf.cpp Rpotrf_gpu_cuda.cu \
$(MPLAPACK_TD_CUDA_LAPACK_OPTIMIZED_SOURCES) \
Rgemm.cpp Rgemm_bridge_cuda.cpp Rgemm_gpu_cuda.cu ../../common/blas_cuda/Rblas_bridge_cuda.cpp ../../common/blas_cuda/Cblas_bridge_cuda.cpp

libmplapack_td_opt_cuda_la_LDFLAGS = $(AM_LDFLAGS) $(CUDA_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO) -Wl,-rpath,$(libdir)
libmplapack_td_opt_cuda_la_LIBADD = $(CUDA_LIBS) -lcublas -lcudart $(QD_LIBADD) $(OPENMP_LIBS)
if IS_MACOS
libmplapack_td_opt_cuda_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmplapack_td_opt_cuda_la_LDFLAGS += -fopenmp -no-undefined
endif

EXTRA_DIST = README.md

check_PROGRAMS = solver_gpu_host_test_td
TESTS = $(check_PROGRAMS)
solver_gpu_host_test_td_SOURCES = ../../../../mplapack/optimized/common/test/Rsolver_gpu_test.cpp
solver_gpu_host_test_td_CPPFLAGS = $(MPLAPACK_TD_CUDA_CPPFLAGS) -DMPLAPACK_GEMM_CUDA_HOST_EMULATION
solver_gpu_host_test_td_CXXFLAGS = $(DD_CXXFLAGS)
solver_gpu_host_test_td_LDADD = ../libmplapack_td_opt.la $(QD_LIBADD) $(OPENMP_LIBS)
solver_gpu_host_test_td_LDFLAGS = -static -Wl,--disable-new-dtags -Wl,-rpath,$(QD_LIBDIR)
TESTS_ENVIRONMENT = OMP_NUM_THREADS=1

clean-local:
	rm -f *.lo.d
