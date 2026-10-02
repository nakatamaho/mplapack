include $(top_srcdir)/mplapack/reference/accelerator-sources.am

# Layered basename shadowing:
# - cuda/Rgemm.cpp replaces ../Rgemm.cpp; it computes 512/1024-bit calls on
#   the GPU and otherwise runs the same OpenMP code as ../Rgemm.cpp.
# - *_cuda.* files are helpers and never shadow.
MPLAPACK_MPFR_CUDA_LAPACK_OPTIMIZED_SOURCES =

# The CUDA sources are relocatable device code (Rgemm_rt_cuda.cu calls the
# cu_mpfr device functions of libmpc_cuda.a); Rgemm_dlink_cuda.lo is their
# device link, so that C++ programs can link the library normally.
NVCC_MPFR_CUDA_FLAGS = $(NVCC_GLIBC_FEATURE_FLAGS) -std=c++17 -rdc=true -fmad=false -I$(MPC_CUDA_INCLUDEDIR) -I$(srcdir) $(MPFR_CUDA_GENCODE)

.cu.o:
	$(NVCC) $(NVCC_MPFR_CUDA_FLAGS) -o $@ -c $< $(CUDA_CFLAGS)

.cu.lo:
	$(top_srcdir)/misc/cudalt.py $@ $(NVCC) $(NVCC_MPFR_CUDA_FLAGS) --compiler-options=\"$(CFLAGS) $(DEFAULT_INCLUDES) $(CPPFLAGS) \" -c $<

Rgemm_dlink_cuda.lo: Rgemm_device_cuda.lo Rgemm_rt_cuda.lo
	$(top_srcdir)/misc/cudalt.py $@ $(NVCC) $(MPFR_CUDA_GENCODE) -dlink .libs/Rgemm_device_cuda.o .libs/Rgemm_rt_cuda.o $(MPC_CUDA_LIBRARY)

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
$(MPLAPACK_MPFR_CUDA_LAPACK_OPTIMIZED_SOURCES) \
Rgemm.cpp Rgemm_bridge_cuda.cpp Rgemm_device_cuda.cu Rgemm_kernel_cuda.h \
Rgemm_rt_cuda.cu Rgemm_rt_cuda.h

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

EXTRA_DIST = README.md Rgemm_host_cuda.cpp test/Rgemm_cuda_test.cpp
