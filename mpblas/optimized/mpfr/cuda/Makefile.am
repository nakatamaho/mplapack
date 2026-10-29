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

# Header dependencies of the .cu objects, written by nvcc through
# misc/cudalt.py (automake does not track the .cu.lo suffix rule).
-include $(wildcard *.lo.d)

MPFR_CBLAS_DEVICE_LO = Cblas_mpfr_device_cuda.lo Cblas_mpfr_device_rt_cuda.lo Cblas_mpfr_device_rt64_cuda.lo Cblas_mpfr_device_rt128_cuda.lo Cblas_mpfr_device_rt256_cuda.lo
Rgemm_dlink_cuda.lo: Rgemm_device_cuda.lo Rgemm_rt_cuda.lo blas3_device_cuda.lo $(MPFR_CBLAS_DEVICE_LO)
	$(top_srcdir)/misc/cudalt.py $@ $(NVCC) $(MPFR_CUDA_GENCODE) -dlink .libs/Rgemm_device_cuda.o .libs/Rgemm_rt_cuda.o .libs/blas3_device_cuda.o $(MPFR_CBLAS_DEVICE_LO:%.lo=.libs/%.o) $(MPC_CUDA_LIBRARY)

# Must match the architectures libmpc_cuda.a was built for (--with-mpfr-cuda-arch,
# e.g. 80,90 for A100 and H100).
MPFR_CUDA_GENCODE = $(MPFR_CUDA_GENCODES) -m64 -O3
NVCC_GLIBC_FEATURE_FLAGS = -U_GNU_SOURCE -D_XOPEN_SOURCE=700

MPLAPACK_MPFR_CUDA_CPPFLAGS = $(GMPFRXX_MKII_CPPFLAGS) -I. -I$(srcdir) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(MPC_CUDA_INCLUDEDIR) $(OPENMP_CXXFLAGS) -DMPLAPACK_BUILD_WITH_MPFR -DMPLAPACK_ACCELERATOR_LIBRARY=1
AM_CPPFLAGS = $(MPLAPACK_MPFR_CUDA_CPPFLAGS)

lib_LTLIBRARIES = libmplapack_mpfr_opt_cuda.la
libmplapack_mpfr_opt_cuda_la_CPPFLAGS = $(MPLAPACK_MPFR_CUDA_CPPFLAGS)
libmplapack_mpfr_opt_cuda_la_SOURCES = \
../Rcopy.cpp ../Rcopy_ref.cpp ../openmp/Rcopy_omp.cpp \
../Rdot.cpp  ../Rdot_ref.cpp  ../openmp/Rdot_omp.cpp \
../../common/blas/Rdot_par.cpp \
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
../../../../mplapack/optimized/common/Rgetrf.cpp \
../../../../mplapack/optimized/common/Cgetrf.cpp \
../../../../mplapack/optimized/common/Rgetrf2.cpp \
../../../../mplapack/optimized/common/Cgetrf2.cpp \
../../../../mplapack/optimized/common/Rpotrf.cpp \
../../../../mplapack/optimized/common/Cpotrf.cpp \
$(MPLAPACK_MPFR_CUDA_LAPACK_OPTIMIZED_SOURCES) \
Rgemm.cpp Rgemm_bridge_cuda.cpp Rgemm_device_cuda.cu Rgemm_kernel_cuda.h \
Rgemm_rt_cuda.cu Rgemm_rt_cuda.h Rgemm_winograd_cuda.h mpfr_pack_cuda.h \
Rsyrk.cpp Rsyr2k.cpp Rsymm.cpp Rtrmm.cpp Rtrsm.cpp \
blas3_bridge_cuda.cpp blas3_device_cuda.cu blas3_kernel_cuda.h \
Cblas_mpfr_bridge_cuda.cpp Cblas_mpfr_device_cuda.cu Cblas_mpfr_dev_types.h Cblas_mpfr_device_cuda.h \
Cblas_mpfr_device_rt_cuda.cu Cblas_mpfr_device_rt64_cuda.cu Cblas_mpfr_device_rt128_cuda.cu Cblas_mpfr_device_rt256_cuda.cu

libmplapack_mpfr_opt_cuda_la_LDFLAGS = $(AM_LDFLAGS) $(CUDA_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO) -Wl,-rpath,$(libdir)
libmplapack_mpfr_opt_cuda_la_LIBADD = Rgemm_dlink_cuda.lo $(MPC_CUDA_LIBRARY) $(CUDA_LIBS) -lcudart $(MPLAPACK_MPFR_LIBS) $(OPENMP_LIBS)
EXTRA_libmplapack_mpfr_opt_cuda_la_DEPENDENCIES = Rgemm_dlink_cuda.lo
CLEANFILES = Rgemm_dlink_cuda.lo .libs/Rgemm_dlink_cuda.o Rgemm_dlink_cuda.o

clean-local:
	rm -f *.lo.d
if IS_MACOS
libmplapack_mpfr_opt_cuda_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmplapack_mpfr_opt_cuda_la_LDFLAGS += -fopenmp -no-undefined
endif

EXTRA_DIST = README.md Rgemm_host_cuda.cpp test/Rgemm_cuda_test.cpp test/Rgemm_winograd_cuda_test.cpp test/blas3_cuda_test.cpp \
    test/cu_fround_test.cpp test/cu_fround_device_test.cu test/cu_fround_bench.cu test/Rtrsm_round_bench.cpp test/run_round_bench.py
