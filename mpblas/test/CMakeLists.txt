# MPBLAS correctness tests, ported from the autotools mpblas/test harness.
#
# Each common/<Routine>.test.cpp is a self-contained program that exercises one
# MPBLAS routine and compares it against the reference Fortran BLAS (dgemm_ etc.,
# declared in blas.h via F77_FUNC); it exit(1)s on a mismatch.  The test support
# code (mplapack.test.cpp) always uses MPFR for random-number generation, so
# these tests require the MPFR backend and a reference BLAS.
#
# Requirements (otherwise this directory is skipped):
#   * the MPFR backend (MPLAPACK_ENABLE_MPFR)
#   * a reference Fortran BLAS (find_package(BLAS)) + its Fortran runtime

if(NOT MPLAPACK_ENABLE_MPFR)
  message(STATUS "mpblas/test: skipped (needs the MPFR backend for its test support)")
  return()
endif()
if(NOT BLAS_FOUND)
  message(STATUS "mpblas/test: skipped (no reference BLAS found)")
  return()
endif()

set(_common "${CMAKE_CURRENT_SOURCE_DIR}/common")

# Support sources compiled into every test program (mirrors test_srcdepends):
#   xerbla.cpp           - the xerbla_ the reference BLAS calls on bad args
#   mplapack.test.cpp    - random generators / comparison helpers (uses MPFR)
#   Mxerbla.override.cpp - test version of Mxerbla; in the executable's own
#                          objects it wins over the backend's via
#                          --allow-multiple-definition
set(_support
  "${_common}/xerbla.cpp"
  "${_common}/mplapack.test.cpp"
  "${_common}/Mxerbla.override.cpp")

# The Fortran runtime that the reference BLAS needs (FCLIBS in autotools).
find_library(MPLAPACK_GFORTRAN_LIBRARY NAMES gfortran)

# Every routine in common/ gets a program per enabled backend.
file(GLOB _mpblas_tests CONFIGURE_DEPENDS "${_common}/*.test.cpp")

set(_backends gmp mpfr qd dd double binary80 binary128)
foreach(_b ${_backends})
  string(TOUPPER ${_b} _B)
  if(NOT MPLAPACK_ENABLE_${_B})
    continue()
  endif()

  foreach(_src ${_mpblas_tests})
    get_filename_component(_stem "${_src}" NAME_WE)   # e.g. Rgemm (from Rgemm.test.cpp)
    set(_exe "${_stem}_test_${_b}")

    add_executable(${_exe} "${_src}" ${_support})
    # Set exactly one backend macro explicitly.  We must NOT obtain it by linking
    # the mpblas_<b> target, because the test also needs the mpblas_mpfr backend
    # and linking both targets would push BOTH ___MPLAPACK_BUILD_WITH_*___ macros
    # onto the compile line (conflicting REAL/COMPLEX typedefs).
    target_compile_definitions(${_exe} PRIVATE
      ___MPLAPACK_BUILD_WITH_${_B}___ ___MPLAPACK_INTERNAL___)
    target_include_directories(${_exe} PRIVATE
      "${CMAKE_BINARY_DIR}/include"          # generated mplapack_config.h
      "${CMAKE_SOURCE_DIR}/include"
      "${CMAKE_SOURCE_DIR}/mpfrc++")         # mpreal/mpcomplex (always used)
    # Link the mpblas archives by file path so they contribute object code only,
    # without their target usage requirements (i.e. without the backend macro).
    target_link_libraries(${_exe} PRIVATE
      "$<TARGET_FILE:mpblas_${_b}>"
      "$<TARGET_FILE:mpblas_mpfr>"
      MPC::MPC MPFR::MPFR GMP::GMPXX   # mpfr/mpc/gmp headers + libs (no macro)
      ${BLAS_LIBRARIES})              # reference Fortran BLAS to compare against
    add_dependencies(${_exe} mpblas_${_b} mpblas_mpfr)
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

    add_test(NAME mpblas_${_b}_${_stem} COMMAND ${_exe})
  endforeach()
endforeach()
