# Telsa code has also been added, but it won't run. We need to tune up this Makefile. Source code itself
# is fine, though.
# TODO:do not hardcode to run fermi routine in Rgemm.cu and Rgemm_cuda.cu
#     :cleanup Makefile so that we pass correct compilation flags. 

.cu.o:
	$(NVCC) $(VOLTA_FLAGS) -o $@ -c $< $(CUDA_CFLAGS)

.cu.lo:
	$(top_srcdir)/misc/cudalt.py $@ $(NVCC) --compiler-options=\"$(CFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \" $(VOLTA_FLAGS) -c $<

lib_LTLIBRARIES = libmpblas_dd_cuda.la

LINK = $(LIBTOOL) --mode=link $(CC) -o $@ $(CUDA_LDFLAGS) $(CUDA_LIBS)

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

TESLA_FLAGS = -gencode arch=compute_13,code=sm_13 -m64 -O3 --ptxas-options=-v 
FERMI_FLAGS = -gencode arch=compute_20,code=sm_21 -m64 -O3 --ptxas-options=-v
VOLTA_FLAGS = -gencode arch=compute_70,code=sm_70 -m64 -O3 --ptxas-options=-v

libmpblas_dd_cuda_la_SOURCES = Rgemm.cu Rgemm_cuda.cu Rgemm_fermi.cu Rgemm_tesla.cu \
				Rsyrk.cu Rsyrk_cuda.cu

