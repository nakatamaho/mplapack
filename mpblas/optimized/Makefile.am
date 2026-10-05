SUBDIRS = 

if ENABLE_GMP
SUBDIRS += gmp
endif

if ENABLE_MPFR
SUBDIRS += mpfr
endif

if ENABLE_QD
SUBDIRS += qd
endif

if ENABLE_DD
SUBDIRS += dd
endif

if ENABLE_DOUBLE
SUBDIRS += double
endif

if ENABLE_BINARY128
SUBDIRS += binary128
endif

if ENABLE_BINARY80
SUBDIRS += binary80
endif

# Header-only Rgemm code shared by the dd, qd and binary128 libraries.
EXTRA_DIST = common/README.md \
common/Rgemm_blocked_common.h common/Rgemm_cpu_common.h \
common/Rgemm_cuda_bridge_common.h common/Rgemm_cuda_common.cuh \
common/Rgemm_cuda_ops_common.h common/Rgemm_ozaki_common.h \
common/Rgemm_qd_ops_common.h common/Rgemm_simd_common.h \
common/Rgemm_simd_kernels_common.h common/Rgemm_winograd_common.h \
common/test/Rgemm_opt_test.cpp common/test/Rgemm_gpu_test.cpp
