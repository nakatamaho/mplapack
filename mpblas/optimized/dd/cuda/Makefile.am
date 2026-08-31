include $(top_srcdir)/mplapack/reference/accelerator-sources.am

# Layered basename shadowing:
# - cuda/Rgemm.cu replaces ../Rgemm.cpp.
# - cuda/Rsyrk.cu replaces ../../../reference/Rsyrk.cpp.
# - *_cuda.cu and split kernels are helpers and never shadow.
MPLAPACK_DD_CUDA_LAPACK_OPTIMIZED_SOURCES =

.cu.o:
	$(NVCC) $(NVCC_GLIBC_FEATURE_FLAGS) $(HOPPER_FLAGS) -o $@ -c $< $(CUDA_CFLAGS)

.cu.lo:
	$(top_srcdir)/misc/cudalt.py $@ $(NVCC) $(NVCC_GLIBC_FEATURE_FLAGS) --compiler-options=\"$(CFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \" $(AMPERE_FLAGS) $(HOPPER_FLAGS) -c $<

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
$(MPLAPACK_DD_CUDA_LAPACK_OPTIMIZED_SOURCES) \
Rgemm.cu Rgemm_cuda.cu Rgemm_fermi.cu Rgemm_tesla.cu \
Rsyrk.cu Rsyrk_cuda.cu dd_real_cuda.h

libmplapack_dd_opt_cuda_la_LDFLAGS = $(AM_LDFLAGS) $(CUDA_LDFLAGS) -version-info $(MPLAPACK_LT_VERSION_INFO) -Wl,-rpath,$(libdir)
libmplapack_dd_opt_cuda_la_LIBADD = $(CUDA_LIBS) -lcudart $(QD_LIBADD)
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
