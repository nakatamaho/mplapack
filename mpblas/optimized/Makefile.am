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
if ENABLE_TD
SUBDIRS += td
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

# Header-only Rgemm code shared by the dd, td, qd and binary128 libraries,
# and the lane-parallel BLAS headers shared by every library (common/blas).
EXTRA_DIST = common/README.md \
common/Rblas_lanes_common.h common/Rblas_calls_common.h common/Rblas_cpu_common.h \
common/Rblas_cuda_common.cuh common/Rblas_cuda_bridge_common.h common/Rblas_simd_common.h \
common/Rblas_split_common.h common/Cblas_simd_common.h \
common/Rblas_dev_types.h common/Cblas_cuda_common.cuh common/Cblas_cuda_jobs_common.h \
common/Rblas_dev_base.h common/Rblas_gpu_operand_common.h common/blas_cuda/Cblas_bridge_body.h \
common/blas_dev/l2_rows_dev.h \
common/blas_dev/CRrot_dev.h \
common/blas_dev/CRscal_dev.h \
common/blas_dev/Caxpy_dev.h \
common/blas_dev/Ccopy_dev.h \
common/blas_dev/Cgemm_dev.h \
common/blas_dev/Cgemmtr_columns_dev.h \
common/blas_dev/Cgemv_dev.h \
common/blas_dev/Cgerc_dev.h \
common/blas_dev/Cgeru_dev.h \
common/blas_dev/Chemm_dev.h \
common/blas_dev/Cher2_columns_dev.h \
common/blas_dev/Cher2k_columns_dev.h \
common/blas_dev/Cher_columns_dev.h \
common/blas_dev/Cherk_columns_dev.h \
common/blas_dev/Cscal_dev.h \
common/blas_dev/Cswap_dev.h \
common/blas_dev/Csymm_dev.h \
common/blas_dev/Csyr2k_columns_dev.h \
common/blas_dev/Csyrk_columns_dev.h \
common/blas_dev/Ctrmm_dev.h \
common/blas_dev/Ctrsm_dev.h \
common/blas_dev/Rgemmtr_columns_dev.h \
common/blas_dev/blas_dev.h \
common/test/Cblas_gpu_test.cpp \
common/blas/Rgemmtr_columns.h common/blas/Csyrk_columns.h common/blas/Cherk_columns.h common/blas/Csyr2k_columns.h common/blas/Cher2k_columns.h common/blas/Cgemmtr_columns.h common/blas/Cher_columns.h common/blas/Cher2_columns.h \
common/Rgemm_blocked_common.h common/Rgemm_cpu_common.h \
common/Rgemm_cuda_bridge_common.h common/Rgemm_cuda_common.cuh \
common/Rgemm_cuda_ops_common.h common/Rgemm_ozaki_common.h \
common/Rgemm_qd_ops_common.h common/Rgemm_simd_common.h \
common/Rgemm_simd_kernels_common.h common/Rgemm_winograd_common.h \
common/test/Rgemm_opt_test.cpp common/test/Rgemm_gpu_test.cpp
