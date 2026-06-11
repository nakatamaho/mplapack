# Telsa code has also been added, but it won't run. We need to tune up this Makefile. Source code itself
# is fine, though.
# TODO:do not hardcode to run fermi routine in Rgemm.cu and Rgemm_cuda.cu
#     :cleanup Makefile so that we pass correct compilation flags. 

lib_LTLIBRARIES = libmpblas__Float128_opt_opencl.la

#.cpp.o:
#	$(CC) -o $@ -c $< $(OPENCL_CFLAGS)

AM_CXXFLAGS = $(OPENCL_CFLAGS)

LINK = $(LIBTOOL) --mode=link $(CC) -o $@ $(OPENCL_LDFLAGS) $(OPECL_LIBS)

EXTRA_DIST = Rgemm_opencl.cpp

libmpblas__Float128_opt_opencl_la_SOURCES = Rgemm_opencl.cpp OclHelper.hpp cpp_kernel.h fp_code4_def.hpp
