# MPLAPACK "compare" tests, ported from the autotools mplapack/test/compare
# harness.  Each common/<Routine>.test.cpp exercises one MPLAPACK routine and
# compares it against the reference Fortran LAPACK/BLAS (dgesv_ etc., declared in
# lapack.h/blas.h via F77_FUNC); it exit(1)s on a mismatch.
#
# Like mpblas/test, the support code always uses MPFR, so this directory needs
# the MPFR backend and reference BLAS + LAPACK; it is skipped otherwise.

if(NOT MPLAPACK_ENABLE_MPFR)
  message(STATUS "mplapack/test/compare: skipped (needs the MPFR backend)")
  return()
endif()
if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)
  message(STATUS "mplapack/test/compare: skipped (no reference BLAS/LAPACK found)")
  return()
endif()

set(_common  "${CMAKE_CURRENT_SOURCE_DIR}/common")
set(_mptest_common "${CMAKE_SOURCE_DIR}/mpblas/test/common")

# Shared support (same files as mpblas/test): xerbla_ for the reference library,
# the MPFR-based random/compare helpers, and the override Mxerbla which wins from
# the executable's own objects via --allow-multiple-definition.
set(_support
  "${_mptest_common}/xerbla.cpp"
  "${_mptest_common}/mplapack.test.cpp"
  "${_mptest_common}/Mxerbla.override.cpp")

file(GLOB _compare_tests CONFIGURE_DEPENDS "${_common}/*.test.cpp")

set(_backends gmp mpfr qd dd double binary80 binary128)
foreach(_b ${_backends})
  string(TOUPPER ${_b} _B)
  if(NOT MPLAPACK_ENABLE_${_B})
    continue()
  endif()

  foreach(_src ${_compare_tests})
    get_filename_component(_stem "${_src}" NAME_WE)
    set(_exe "${_stem}_compare_${_b}")

    set(_extra "")
    if(_stem STREQUAL "Mutils")
      set(_extra "${_common}/ilaenv.cpp")   # Mutils.test needs the test ILAENV
    endif()

    add_executable(${_exe} "${_src}" ${_support} ${_extra})
    # Exactly one backend macro (set explicitly; the mplapack/mpblas archives are
    # linked by file path below so their target macros do not leak in).
    target_compile_definitions(${_exe} PRIVATE
      ___MPLAPACK_BUILD_WITH_${_B}___ ___MPLAPACK_INTERNAL___)
    target_include_directories(${_exe} PRIVATE
      "${CMAKE_BINARY_DIR}/include"
      "${CMAKE_SOURCE_DIR}/include"
      "${CMAKE_SOURCE_DIR}/mpfrc++")
    target_link_libraries(${_exe} PRIVATE
      "$<TARGET_FILE:mplapack_${_b}>"
      "$<TARGET_FILE:mplapack_mpfr>"
      "$<TARGET_FILE:mpblas_${_b}>"
      "$<TARGET_FILE:mpblas_mpfr>"
      MPC::MPC MPFR::MPFR GMP::GMPXX
      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})   # reference Fortran LAPACK + BLAS
    add_dependencies(${_exe}
      mplapack_${_b} mplapack_mpfr mpblas_${_b} mpblas_mpfr)
    if(MPLAPACK_GFORTRAN_LIBRARY)
      target_link_libraries(${_exe} PRIVATE ${MPLAPACK_GFORTRAN_LIBRARY})
    endif()
    if(TARGET OpenMP::OpenMP_CXX)
      target_link_libraries(${_exe} PRIVATE OpenMP::OpenMP_CXX)
    endif()
    if(MPLAPACK_LD_ALLOW_MULTIDEF)
      target_link_options(${_exe} PRIVATE "-Wl,--allow-multiple-definition")
    endif()
    set_target_properties(${_exe} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${_b}")

    add_test(NAME compare_${_b}_${_stem} COMMAND ${_exe})
  endforeach()
endforeach()
