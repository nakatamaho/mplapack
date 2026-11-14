# Each family also enters with restore enabled. The probes themselves switch
# both policies in-process and restore this environment for the numerical sweep.
if(MPLAPACK_BUILD_TESTS AND TARGET mplapack_mpfr_opt_cuda)
  add_executable(mpfr_cuda_stack_guard_host_test
    "${PROJECT_SOURCE_DIR}/mpblas/optimized/mpfr/cuda/test/stack_guard_test.cpp")
  target_include_directories(mpfr_cuda_stack_guard_host_test BEFORE PRIVATE
    "${PROJECT_SOURCE_DIR}/mpblas/optimized/mpfr/cuda/test/stack_guard_stub")
  target_compile_features(mpfr_cuda_stack_guard_host_test PRIVATE cxx_std_17)
  add_test(NAME mpfr_cuda_stack_guard_host COMMAND mpfr_cuda_stack_guard_host_test)
  function(mplapack_mpfr_restore_test name executable)
    add_test(NAME ${name} COMMAND ${executable} ${ARGN})
    set_tests_properties(${name} PROPERTIES SKIP_RETURN_CODE 77
      ENVIRONMENT "OMP_NUM_THREADS=1;MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=0;MPLAPACK_MPFR_CUDA_STACK_RESTORE=1;MPLAPACK_MPFR_CUDA_MIN_MNK=0;MPFRXX_DEFAULT_PRECISION_BITS=512")
  endfunction()
  mplapack_mpfr_restore_test(mpfr_cuda_solver_device_256_restore mpfr_cuda_solver_device_test 256)
  mplapack_mpfr_restore_test(mpfr_cuda_Rfactor_device_restore mpfr_cuda_Rfactor_device_test 256)
  mplapack_mpfr_restore_test(mpfr_cuda_Cfactor_device_restore mpfr_cuda_Cfactor_device_test 256)
  mplapack_mpfr_restore_test(mpfr_cuda_inverse_device_512_restore mpfr_cuda_inverse_device_test)
  mplapack_mpfr_restore_test(mpfr_cuda_Rgemm_device_runtime_restore mpfr_cuda_Rgemm_device_test 200)
  mplapack_mpfr_restore_test(mpfr_cuda_Cblas_device_fixed_restore mpfr_cuda_Cblas_device_test)
  mplapack_mpfr_restore_test(mpfr_cuda_blas3_device_restore mpfr_cuda_blas3_device_test --fixed-only)
endif()
