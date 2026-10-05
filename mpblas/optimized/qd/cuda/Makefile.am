include $(top_srcdir)/mplapack/reference/accelerator-sources.am

# Layered basename shadowing:
# - cuda/Rgemm.cpp replaces ../Rgemm.cpp; it computes eligible calls on the
#   GPU (Rgemm_bridge_cuda.cpp, Rgemm_gpu_cuda.cu) and otherwise runs the CPU
#   code of ../Rgemm.cpp.
# - *_cuda.* files are helpers and never shadow.
MPLAPACK_QD_CUDA_LAPACK_OPTIMIZED_SOURCES =

# -fmad=false: the quad-double device arithmetic must not be contracted.
NVCC_QD_CUDA_FLAGS = -std=c++17 -fmad=false -I$(srcdir) -I$(QD_INCLUDEDIR)

.cu.o:
	$(NVCC) $(NVCC_GLIBC_FEATURE_FLAGS) $(NVCC_QD_CUDA_FLAGS) $(HOPPER_FLAGS) -o $@ -c $< $(CUDA_CFLAGS)

.cu.lo:
	$(top_srcdir)/misc/cudalt.py $@ $(NVCC) $(NVCC_GLIBC_FEATURE_FLAGS) $(NVCC_QD_CUDA_FLAGS) --compiler-options=\"$(CFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CPPFLAGS) \" $(AMPERE_FLAGS) $(HOPPER_FLAGS) -c $<

AMPERE_FLAGS = -gencode arch=compute_80,code=sm_80 -m64 -O3 --ptxas-options=-v
HOPPER_FLAGS = -gencode arch=compute_90,code=sm_90 -m64 -O3 --ptxas-options=-v
NVCC_GLIBC_FEATURE_FLAGS = -U_GNU_SOURCE -D_XOPEN_SOURCE=700

MPLAPACK_QD_CUDA_CPPFLAGS = -I. -I.. -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) $(OPENMP_CXXFLAGS) -DMPLAPACK_BUILD_WITH_QD
AM_CPPFLAGS = $(MPLAPACK_QD_CUDA_CPPFLAGS)

lib_LTLIBRARIES = libmplapack_qd_opt_cuda.la
libmplapack_qd_opt_cuda_la_CPPFLAGS = $(MPLAPACK_QD_CUDA_CPPFLAGS)
libmplapack_qd_opt_cuda_la_CXXFLAGS = $(DD_CXXFLAGS)
libmplapack_qd_opt_cuda_la_SOURCES = \
../Rcopy.cpp ../Rcopy_ref.cpp ../openmp/Rcopy_omp.cpp \
../Rdot.cpp  ../Rdot_ref.cpp  ../openmp/Rdot_omp.cpp \
../Rgemm_ref.cpp ../openmp/Rgemm_blocked_omp.cpp ../openmp/Rgemm_NN_omp.cpp ../openmp/Rgemm_NT_omp.cpp ../openmp/Rgemm_TN_omp.cpp ../openmp/Rgemm_TT_omp.cpp \
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
$(MPLAPACK_QD_CUDA_LAPACK_OPTIMIZED_SOURCES) \
Rgemm.cpp Rgemm_bridge_cuda.cpp Rgemm_gpu_cuda.cu

libmplapack_qd_opt_cuda_la_LDFLAGS = $(AM_LDFLAGS) $(CUDA_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO) -Wl,-rpath,$(libdir)
libmplapack_qd_opt_cuda_la_LIBADD = $(CUDA_LIBS) -lcublas -lcudart $(QD_LIBADD) $(OPENMP_LIBS)
if IS_MACOS
libmplapack_qd_opt_cuda_la_LDFLAGS += -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
if IS_MINGW
libmplapack_qd_opt_cuda_la_LDFLAGS += -fopenmp -no-undefined
endif
